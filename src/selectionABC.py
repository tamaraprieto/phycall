import pickle
import argparse
import os
import sys

import numpy as np
import pandas as pd
import pymc as pm
from sklearn.metrics import auc


def wf_selection_auc(fitness, aoo, population_size=100, max_generations=10000):
    while True:
        # Initialize population
        pop = [np.zeros(population_size, dtype=np.int8)]
        pop[0][0] = 1

        # Pre-calculate fitness probabilities
        fit_probs = np.array([1 - fitness, fitness], dtype=np.float64)

        for i in range(max_generations):
            # Use fitness values directly
            probs = fit_probs[pop[-1]]

            # Normalize probabilities
            probs /= probs.sum()

            # Update population
            new_pop = np.random.choice(pop[-1], population_size, p=probs)

            # Append new population
            pop.append(new_pop)

            # Break conditions
            total_ones = new_pop.sum()
            if total_ones == 0 or total_ones == population_size:
                break
        # Recursion to enforce outgrowth
        if total_ones != 0:
            break

    pop = np.array(pop)
    traj = pop.mean(axis=1)
    wfs_auc = auc(np.arange(len(traj)), traj) * aoo
    return wfs_auc


def sim(rng, fitness, aoo, size=None):
    if any([fitness.item() < 0.5, fitness.item() > 1.0, aoo.item() < 0]):
        return 0
    return wf_selection_auc(fitness, aoo)


def parse_ltt_file(path):
    ltt = pd.read_csv(path)
    if "Unnamed: 2" in ltt.columns:
        ltt.drop(columns=["Unnamed: 2"], inplace=True)
    ltt["lineages_norm"] = ltt["lineages"] / ltt["lineages"].max()
    ltt["time_yr"] = ltt["time"]  # / 365
    return auc(ltt["time_yr"], ltt["lineages_norm"]) * abs(ltt["time_yr"].min())


def run_abc(
    aoo_prior: int,
    observed_ltt: float,
    chains: int = 6,
    draws: int = 10000,
    epsilon: float = 1,
):
    with pm.Model() as model:
        fitness = pm.Uniform(
            "fitness",
            lower=0.5,
            upper=1.0,
        )
        aoo = pm.Uniform(
            "aoo",
            lower=0,
            upper=aoo_prior,
        )
        s = pm.Simulator(
            "s",
            sim,
            params=(
                fitness,
                aoo,
            ),
            epsilon=epsilon,
            observed=observed_ltt,
        )

        idata = pm.sample_smc(progressbar=True, chains=chains, draws=draws)
        idata.extend(pm.sample_posterior_predictive(idata))
    return idata


def main():
    parser = argparse.ArgumentParser(
        description="A tool for running Approximate Bayesian Computation (ABC) simulations of Wright-Fisher with selection clone dynamics.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input_ltt",
        required=True,
        type=str,
        help="Path to the observed Lineage Through Time (LTT) file.",
    )
    parser.add_argument(
        "-o",
        "--output_path",
        required=True,
        type=str,
        help="Path to the output file where the result of the simulation will be stored.",
    )
    parser.add_argument(
        "-a",
        "--age_of_patient",
        required=True,
        type=float,
        help="Age of the patient in years. This will be used as the upper limit of the prior distribution for the age of onset of the clone.",
    )
    parser.add_argument(
        "-e",
        "--epsilon",
        required=False,
        type=float,
        default=1,
        help="ABC epsilon value. This parameter controls the tolerance of the ABC rejection algorithm.",
    )
    parser.add_argument(
        "-c",
        "--chains",
        required=False,
        type=int,
        default=6,
        help="Number of Sequential Monte Carlo (SMC) chains to be used in the simulation.",
    )
    parser.add_argument(
        "-d",
        "--draws",
        required=False,
        type=int,
        default=10000,
        help="Number of draws per SMC chain.",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.input_ltt):
        print(f"Input file {args.input_ltt} does not exist.")
        sys.exit(1)

    if (
        not os.path.isdir(os.path.dirname(args.output_path))
        and os.path.dirname(args.output_path) != ""
    ):
        print(f"Output directory {os.path.dirname(args.output_path)} does not exist.")
        sys.exit(1)

    if (
        args.age_of_patient <= 0
        or args.epsilon <= 0
        or args.chains <= 0
        or args.draws <= 0
    ):
        print("Age of patient, epsilon, chains and draws should be positive numbers.")
        sys.exit(1)

    observed = parse_ltt_file(args.input_ltt)
    idata = run_abc(
        aoo_prior=args.age_of_patient,
        observed_ltt=observed,
        epsilon=args.epsilon,
        chains=args.chains,
        draws=args.draws,
    )

    with open(args.output_path + ".pkl", "wb") as f:
        pickle.dump(idata, f)


if __name__ == "__main__":
    main()
