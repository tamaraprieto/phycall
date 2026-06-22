#!/usr/bin/python

def usage():
    print("""
    This script is intended for modifying files from vcftools that contain 
    biallelic information and convert them to fasta format. After the vcf file has
    been exported using the vcf-to-tab program from VCFTools, and transposed in R,
    or Excel, this script will change biallelic information at each site to only 
    one nucleotide using UIPAC conventions when the alleles are different 
    (heterozygous) or to the available nucleotide if both alleles are the same. 
    If one alle is present and the other one is missing, the script will change 
    the site to the available allele. All of these changes will be saved to a 
    new file in fasta format.
    
    
    written by Simon Uribe-Convers - www.simonuribe.com
    October 23rd, 2017
    
    To use this script, type: python3.6 VCF-to-Tab_to_Fasta_IUPAC_Converter.py VCF-to-Tab_file Output_file
    """)


import sys
import __future__

if __name__ == "__main__":
    if len(sys.argv) != 3:
        usage()
        print("~~~~Error~~~~\n\nCorrect usage: python3.6 "+sys.argv[0]+" VCF-to-Tab file + Output file")
        sys.exit("Missing either the VCF-to-Tab and/or the output files!")

filename = open(sys.argv[1], "r")

outfile = open(sys.argv[2] + ".fasta", "w")


# def IUPAC_converter(filename, outfile):
# NEXUS FORMAT HAS ? INSTEAD OF N
IUPAC_Codes = { "G/G" : "GG", "C/C" : "CC", "T/T" : "TT", "A/A" : "AA", 
"G/T" : "GT", "T/G" : "TG", "A/C" : "AC", "C/A" : "CA", "C/G" : "CG",
"G/C" : "GC", "A/G" : "AG", "G/A" : "GA", "A/T" : "AT", "T/A" : "TA",
"C/T" : "CT", "T/C" : "TC", "./." : "??", "N/." : "??", "./N" : "??"}

for line in filename:
    species_name = line.strip().split(" ")[0]
    data = line.strip().split(" ")[1:]
    new_data = [IUPAC_Codes[i] for i in data]
    # print(new_data)
    new_data2 = "".join(new_data)
    outfile.write(">" + species_name + "\n" + new_data2 + "\n")
    
print("\n\n~~~~\n\nResults can be found in %s.fasta\n\n~~~~" %sys.argv[2])

sys.exit(0)
