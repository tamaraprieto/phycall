#!/usr/bin/env python
# coding: utf-8

# In[13]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import collections
import random
import pickle
import seaborn as sns
import scipy
import os
from matplotlib import cm
import matplotlib.colors as mcolors
import csv



# # 1. Input

# Input patient name and folders paths

# In[ ]:


# Change patient name here
# patient = '4295_hg38'
# dataset = '4295_hg38'

# or 
# patient = '417'
# dataset = '417_hg38'

# patient = '445'
# dataset = '445_hg38'


patient = 'Invitro'
dataset = 'Invitro_hg38'

# Whether or not to use Julien's optimization algorithm for allelic imbalance calculation
use_algo = False


# In[3]:

# Valid for all but invitro
filter_SNPS=10
filter_counts=4*filter_SNPS
coverage_threshold = 10


# In[ ]:


format_type = 1 # folder formating changed from RESULTS/Sample/BAF_XXk to Sample/RESULTS/BAF/BAF_XX000


if dataset=='4295_hg38':
    #4. Get the cell order
    treeidx=["D9",  "A11", "C5",  "G8",  "D5",  "E7",  "D4",  "F4",  "D3",  "E9",  "A2",  "D8",  "A3",  "F9",  "F10", "H6",  "A8",  "B10", "D1",  "H2",  "C1",  "D2",  "G7",  "G1",  "H11", "B2",  "A7",  "F6",
    "A4",  "C2",  "E11", "D6",  "C4",  "F7",  "E8",  "F11", "G6",  "B3",  "D11", "H5", "H8",  "C9",  "D10", "C8", "C11", "B7",  "E6",  "G4",  "F5",  "C7",  "H10", "G10", "B6",  "F2",  "E4",  "E3", 
    "A9",  "D7",  "E10", "B5",  "F8",  "E1",  "C10", "F1",  "G5",  "G11", "B8",  "C6",  "E5",  "B4",  "G2",  "H7",  "A5",  "A6",  "E2",  "B1",  "G9",  "G3",  "A10", "H1",  "F3",  "H9",  "H4",  "A1", "C3"]

    #remove the noisy cells
    treeidx_corrected = list(treeidx)
    noisy_cells = ['F6','H9','B7','D7'] #D7 because of Ginkgo ploidy...
    for cell in noisy_cells:
        treeidx_corrected.remove(cell)

    truth_regions = {'1': {'chr': 6, 'bounds': (101e6,129e6), 'cells': ['D5','E7','D4','F4','D3','E9','A2','D8','A3','F9','F10','H6','A8','B10','D1','H2','C1','D2','G7','H11','B2','A7','A4','C2','E11','D6','C4','F7','E8','F11','G6','B3','D11','H5','H8','C9','D10','C8','C11','E6','G4','F5']},
                 '2': {'chr': 11, 'bounds': (40e6,50e6), 'cells':['A5','A6','E2','B1','G9','G3','A10','H1','F3','H4']}}


    binsize_dict={'2_5M': 2500000,
                '500k': 500000,
                '250k': 250000,
                '100k': 100000,
                '50k': 50000,
                '25k': 25000,
                '10k': 10000}

    binsize_list = ['10k', '25k', '50k', '100k', '250k', '500k', '2_5M']



    if use_algo:
        dir_reads='/gpfs/commons/home/jlaval/PreprocessingResults/'
        results_dir = '/gpfs/commons/home/jlaval/RESULTS/'

        BAF_loc = '{}{}/BAF_{}/'#.format(results_dir,dataset,binsize)
    else:
        dir_reads = '/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/WGS/Invitro/RESULTS/'

        results_dir = '/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/WGS/Invitro/RESULTS/BAF/'
        BAF_loc = '{}/BAF_{}/'   # no _no_slide_no_algo suffix — folders don't have it
        format_type = 2

    healthy = '{}_UndeterminedBarcode'.format(dataset) # 22X depth, SNPS with less than 10X coverage are removed??

elif dataset=='445_hg38':
    #445_hg38
    treeidx = ['B8', 'G6', 'E8', 'C11', 'C5', 'C6', 'H11', 'F7', 'G10', 'B6',
                'H5', 'D3', 'E5', 'A11', 'E3', 'D10', 'D4', 'H10', 'E2', 'E10', 
                'H3', 'A8', 'H7', 'D7', 'B2', 'D6', 'H6', 'F4', 'H4', 'D8', 'G9', 
                'A3', 'E12']
    #remove the noisy cells
    treeidx_corrected = list(treeidx)
    # noisy_cells = ['A3', 'G9', 'C5', 'E2']
    noisy_cells = []
    for cell in noisy_cells:
        treeidx_corrected.remove(cell)


    truth_regions = {'1': {'chr': 7, 'bounds': (25000,30150000), 'cells': ['B8', 'G6', 'E8', 'C11', 'C5', 'C6', 'H11', 'F7', 'G10', 'B6',
                                                                        'H5', 'D3', 'E5', 'A11', 'E3', 'D10', 'D4', 'H10', 'E2', 'E10', 
                                                                        'H3', 'A8', 'H7', 'D7', 'B2', 'D6', 'H6']}
                 }


    binsize_dict={'2_5M': 2500000,
                '500k': 500000,
                '250k': 250000,
                '100k': 100000,
                '50k': 50000,
                '25k': 25000,
                '10k': 10000}
                
    binsize_list = ['10k', '25k', '50k', '100k', '250k', '500k', '2_5M']

    dir_reads= '/gpfs/commons/home/jquentin/tamara_project/test_tamara/PreprocessingResults/'
    results_dir = '/gpfs/commons/home/jquentin/tamara_project/test_tamara/RESULTS/'

    if use_algo:
        BAF_loc = '{}{}/BAF_{}/'#.format(results_dir,dataset,binsize)
    else:
        BAF_loc = '{}{}/BAF_{}_no_slide_no_algo/'#.format(RESULTS_dir,dataset,binsize)

    #healthy bulk
    healthy = '445_445T' # to check

elif dataset=='417_hg38':
    #445_hg38
    treeidx = [
                #'368B',
                 'C6', 'D8', 'H8', 'G5', 'A5', 'A12', 'E5', 'E7', 'G2', 'G6', 'D7', 'E8', 'C10', 'E6', 'F3',
                'H3', 'H7', 'B9', 'G9', 'D9', 'A2', 'C5', 'B5', 'B10']

    #remove the noisy cells
    treeidx_corrected = list(treeidx)
    # noisy_cells = ['B5', 'F3', 'B9', 'G2', 'C10']
    noisy_cells = []
    for cell in noisy_cells:
        treeidx_corrected.remove(cell)


    truth_regions = {'1': {'chr': 7, 'bounds': (29150000,51500000), 'cells': [
               # '368B',
                'C6', 'D8', 'H8', 'G5', 'A5', 'A12', 'E5', 'E7', 'G2', 'G6', 'D7', 'E8', 'C10', 'E6', 'F3',
                'H3', 'H7', 'B9', 'G9']}
                 }


    binsize_dict={'2_5M': 2500000,
                '500k': 500000,
                '250k': 250000,
                '100k': 100000,
               '50k': 50000,
                '25k': 25000,
                '10k': 10000}
                
    binsize_list = ['10k', '25k', 
    '50k',
     '100k', '250k', '500k', '2_5M']

    dir_reads= '/gpfs/commons/home/jquentin/tamara_project/test_tamara/PreprocessingResults/'
    results_dir = '/gpfs/commons/home/jquentin/tamara_project/test_tamara/RESULTS/'

    if use_algo:
        BAF_loc = '{}{}/BAF_{}/'#.format(results_dir,dataset,binsize)
    else:
        BAF_loc = '{}{}/BAF_{}_no_slide_no_algo/'#.format(RESULTS_dir,dataset,binsize)

    #healthy bulk
    healthy = '417_417T' # to check

elif dataset=='4084_hg38':
    print("Error, no big deletions were found for this event")
    assert False
    
    #445_hg38
    treeidx = [XX]
    #remove the noisy cells
    treeidx_corrected = list(treeidx)
    # noisy_cells = ['A3', 'G9', 'C5', 'E2']
    noisy_cells = []
    for cell in noisy_cells:
        treeidx_corrected.remove(cell)


    truth_regions = {'1': {'chr': XX, 'bounds': (XX,XX), 'cells': []}
                 }


    binsize_dict={'25k': 25000}

    binsize_list = ['25k']

    dir_reads= '/gpfs/commons/home/jquentin/tamara_project/test_tamara/PreprocessingResults/'
    results_dir = '/gpfs/commons/home/jquentin/tamara_project/test_tamara/RESULTS/'

    if use_algo:
        BAF_loc = '{}{}/BAF_{}/'#.format(results_dir,dataset,binsize)
    else:
        BAF_loc = '{}{}/BAF_{}_no_slide_no_algo/'#.format(RESULTS_dir,dataset,binsize)


    #healthy bulk
    healthy = XX # to check



elif dataset=='Invitro_hg38':
    
    #Invitro_hg38
    treeidx = ['A1', 'A2', 'B1', 'B2', 'Bulk', 'C1', 'C2', 'D1', 'D2', 'E1', 'E2', 'F1', 'F2', 'G1', 'G2', 'H1', 'H2']
    treeidx = ['A1', 'A2', 'B1', 'B2', 'C1', 'C2', 'D1', 'D2', 'E1', 'E2', 'F1', 'F2', 'G1', 'G2', 'H1', 'H2']

    #remove the noisy cells
    treeidx_corrected = list(treeidx)
    # noisy_cells = ['A3', 'G9', 'C5', 'E2']
    noisy_cells = []
    for cell in noisy_cells:
        treeidx_corrected.remove(cell)


    truth_regions = {'1': {'chr': 13, 'bounds': (41875000,42400000), 'cells': ['D1', 'D2', 'E1', 'E2', 'H1', 'H2']}
                 }

    binsize_dict={'100k': 100000,
                '50k': 50000,
                '25k': 25000,
                '10k': 10000,
                '5k': 5000}

    binsize_list = ['5k', '10k', '25k', '50k', '100k']



    dir_reads='/gpfs/commons/groups/landau_lab/tprieto/jlaval/PreprocessingResults_backup/' #unused

    results_dir = '/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/WGS/Invitro/RESULTS/BAF/'

    if use_algo:
        # not done
        raise ValueError
    else:
        BAF_loc = '{}/BAF_{}/' #change of format
        format_type = 2

    #healthy bulk
    healthy = 'Invitro_Bulk'

    filter_SNPS=2
    filter_counts=10
    coverage_threshold = 10



# In[5]:


cwd = os.getcwd()
save_dir = '{}/{}/'.format(cwd,dataset)

if not os.path.exists(save_dir):
   # Create a new directory because it does not exist
   os.makedirs(save_dir)


# In[6]:


save_dir


# In[14]:


with open('{}/cell_indexes.csv'.format(dataset), 'w') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(treeidx)


# # 2. Run AUC processing

# In[7]:


def extract_BAF(open_dir: str, chr: str, start: int, end: int, bin_size: int, patient: str, treeidx_corrected):
    """
    Extracts BAF outputs from pkl files
    """
    print(open_dir+'dict_chr{}.pkl'.format(chr))
    with open(open_dir+'dict_chr{}.pkl'.format(chr), 'rb') as f:
        dict_BAF_chr = pickle.load(f, encoding='latin1')


    # BAF = np.vstack([dict_BAF_chr[patient+'_'+cell]['BAF'] for cell in treeidx_corrected])
    # SNPS = np.vstack([dict_BAF_chr[patient+'_'+cell]['SNP'] for cell in treeidx_corrected])
    # COUNTS = np.vstack([dict_BAF_chr[patient+'_'+cell]['COUNTS'] for cell in treeidx_corrected])

    # output is different whether we use slides or not
    if len(dict_BAF_chr[patient+'_'+treeidx_corrected[0]]['BAF'][0])==1:
        BAF = np.swapaxes(np.hstack(list(np.array(dict_BAF_chr[patient+'_'+cell]['BAF']) for cell in treeidx_corrected)), 0, 1)
        SNPS = np.swapaxes(np.hstack(list(np.array(dict_BAF_chr[patient+'_'+cell]['SNP']) for cell in treeidx_corrected)), 0, 1)
        COUNTS = np.swapaxes(np.hstack(list(np.array(dict_BAF_chr[patient+'_'+cell]['COUNTS']) for cell in treeidx_corrected)), 0, 1)
        #keep track of the chromosome and position
        POSITION=bin_size*np.arange(BAF.shape[1])
        CHROMOSOME=np.array([chr for k in range(BAF.shape[1])])

    else:
        BAF = np.vstack([dict_BAF_chr[patient+'_'+cell]['BAF'] for cell in treeidx_corrected])
        SNPS = np.vstack([dict_BAF_chr[patient+'_'+cell]['SNP'] for cell in treeidx_corrected])
        COUNTS = np.vstack([dict_BAF_chr[patient+'_'+cell]['COUNTS'] for cell in treeidx_corrected])
        #keep track of the chromosome and position
        POSITION=bin_size*np.arange(BAF.shape[1])
        CHROMOSOME=np.array([chr for k in range(BAF.shape[1])])

    #subsample the region between start and end
    idx_start = int(np.ceil(start / bin_size))
    idx_end = int(end // bin_size)

    BAF = BAF[:,idx_start:idx_end]
    SNPS = SNPS[:,idx_start:idx_end]
    COUNTS = COUNTS = COUNTS[:,idx_start:idx_end]

    return dict_BAF_chr, idx_start, idx_end, BAF, SNPS, COUNTS, POSITION, CHROMOSOME


# In[8]:


def extract_BAF_mat(BAF: np.array, SNPS: np.array, COUNTS: np.array, filter_SNPS: int, filter_counts: int):
    """
    Create BAF matrix from array input by filtering for SNP count, read count and taking maximum if sliding is used
    """
    #apply filters and select a unique BAF value per bin
    BAF_mat = np.zeros(BAF.shape)
    mask_mat_snps = np.ones(BAF.shape)
    mask_mat_counts = np.ones(BAF.shape)
    
    for i in range(BAF.shape[0]):
        for j in range(BAF.shape[1]):

            baf_slides = BAF[i,j]
            snps_slides = SNPS[i,j]
            counts_slides = COUNTS[i,j]
            
            if type(snps_slides)==list:
                # #filter by number of SNPS and counts
                for k,(snp,count) in enumerate(zip(snps_slides,counts_slides)):
                    if snp < filter_SNPS or count < filter_counts:
                        # baf_slides[k]=0.5 ##ISSUE: THE BAF ARRAY ITSELF IS MODIFIED, FIND A COPY METHOD IN PYTHON 2
                        baf_slides[k]=np.nan

                #take the maximum imbalance
                idx=np.argmax(np.abs(np.array(baf_slides)-0.5))
                baf = baf_slides[idx]

            else:
                if snps_slides < filter_SNPS:
                    # Filter these away
                    baf_slides=np.nan
                else:
                    mask_mat_snps[i,j]=0

                if counts_slides < filter_counts:
                    # Filter these away
                    baf_slides=np.nan
                else:
                    mask_mat_counts[i,j]=0

                baf = baf_slides

            BAF_mat[i,j]=baf

    #counts masked only if not already masked by snps
    mask_mat_counts = np.logical_and(mask_mat_counts, 1-mask_mat_snps)

    return BAF_mat, mask_mat_counts, mask_mat_snps



# In[9]:


def returnROC_Curve(BAF_mat, positive_cells_idx, negative_cells_idx, mask_mat_snps, mask_mat_counts):
    """
    Returns ROC coordinates
    """
    #separate the positive and negative cells
    BAF_positive_list = []
    BAF_negative_list = []

    for k in range(BAF_mat.shape[0]):
        if k in positive_cells_idx:
            BAF_positive_list.append(BAF_mat[k])
        if k in negative_cells_idx:
            BAF_negative_list.append(BAF_mat[k])

    BAF_positive = np.vstack(BAF_positive_list)
    BAF_negative = np.vstack(BAF_negative_list)

    #same in the mask distribution
    mask_mat = np.logical_or(mask_mat_snps, mask_mat_counts)
    
    mask_positive_list = []
    mask_negative_list = []

    for k in range(mask_mat.shape[0]):
        if k in positive_cells_idx:
            mask_positive_list.append(mask_mat[k])
        if k in negative_cells_idx:
            mask_negative_list.append(mask_mat[k])

    mask_positive = np.vstack(mask_positive_list)
    mask_negative = np.vstack(mask_negative_list)

    list_TPR = []
    list_FPR = []

    for threshold in range(-1,51):

        thresh = threshold / 100.

        #threshold the BAF
        CNVs_positive = np.abs(BAF_positive-0.5) > thresh
        CNVs_negative = np.abs(BAF_negative-0.5) > thresh

        TP=np.sum(np.logical_and(CNVs_positive,1-mask_positive))
        FN=np.sum(np.logical_and(1-CNVs_positive,1-mask_positive))
        TN=np.sum(np.logical_and(1-CNVs_negative,1-mask_negative))
        FP=np.sum(np.logical_and(CNVs_negative,1-mask_negative))

        TPR = float(TP)/(TP+FN)
        FPR = float(FP)/(FP+TN)

        list_TPR.append(TPR)
        list_FPR.append(FPR)

    return list_TPR, list_FPR
    


# In[10]:


#### Compute perf

dict_perf = {region: {binsize : {'TPR': [], 'FPR': [], 'filtered': {}} for binsize in binsize_list} for region in truth_regions.keys()}

for region in truth_regions.keys():
    
    print(region)
    chr=truth_regions[region]['chr']
    start,end=truth_regions[region]['bounds']
    positive_cells=truth_regions[region]['cells']
    positive_cells_idx = [treeidx_corrected.index(cell) for cell in positive_cells] # we keep only the cells where we know there is an issue as the positive cells
    
    for binsize in binsize_dict.keys():
        
        print(binsize)
        bin_size=binsize_dict[binsize]

        #first, stack the BAF in a single matrix
        # open_dir = results_dir+'{}/BAF_{}/'.format(dataset,binsize)
        if format_type==1:
            open_dir = BAF_loc.format(results_dir, dataset, binsize)
        else:
            open_dir = BAF_loc.format(results_dir, bin_size)
        

        dict_BAF_chr, idx_start, idx_end, BAF, SNPS, COUNTS, POSITION, CHROMOSOME = extract_BAF(open_dir, chr, start, end, bin_size, patient, treeidx_corrected)

        BAF_mat, mask_mat_counts, mask_mat_snps = extract_BAF_mat(BAF, SNPS, COUNTS, filter_SNPS, filter_counts)

        # Save BAF matrix to csv
        imbalance_BAF_mat = np.where(BAF_mat<0.5, 1-BAF_mat, BAF_mat)
        np.savetxt("{}/region_{}_chr{}_{}_imbalance.csv".format(dataset, region, chr, binsize), imbalance_BAF_mat, delimiter="\t")

        negative_cells_idx = [k for k in range(mask_mat_counts.shape[0]) if k not in positive_cells_idx]


        #Sampling for boostrapping estimation of confidence interval
        N = 100
        positive_cells_k_fold = [np.random.choice(np.array(positive_cells_idx), replace=True, size=len(positive_cells_idx)) for k in range(N)]
        negative_cells_k_fold = [np.random.choice(np.array(negative_cells_idx), replace=True, size=len(negative_cells_idx)) for k in range(N)]

        all_tprs = []
        all_fprs = []
        for fold in range(N):
            list_TPR, list_FPR = returnROC_Curve(BAF_mat, list(positive_cells_k_fold[fold]), list(negative_cells_k_fold[fold]), mask_mat_snps, mask_mat_counts)
            all_tprs.append(list_TPR) 
            all_fprs.append(list_FPR)
        
        
        boostrap_var = np.var(np.stack( all_tprs, axis=0), axis=0)*N/(N-1)

        dict_perf[region][binsize]['TPR']=np.mean(np.stack( all_tprs, axis=0), axis=0)
        dict_perf[region][binsize]['TPR_lower']=np.mean(np.stack( all_tprs, axis=0), axis=0)-2*np.sqrt(boostrap_var)
        dict_perf[region][binsize]['TPR_higher']=np.mean(np.stack( all_tprs, axis=0), axis=0)+2*np.sqrt(boostrap_var)
        
        dict_perf[region][binsize]['FPR']=list_FPR
        # dict_perf[region][binsize]['filtered']=filtered


# In[11]:


### plot AUC and filtered bins on a two axes plot
### do one plot for each region

c1=[0.7,0.2,0.2,1]
c2=[0.2,0.2,0.7,1]
c3=[0.1,0.1,0.1,1]
cmap=[c1,c2,c3] 
width=0.3
offset=12

color_list = list(mcolors.TABLEAU_COLORS.keys())
plt.close()

sub_bin_size_list = binsize_list#['100k']
# sub_bin_size_list = ['10k', '25k', '50k']

for region, perf in dict_perf.items():

    fig,axs = plt.subplots(ncols=1,figsize=(5,5), dpi=200)

    for i, binsize in enumerate(sub_bin_size_list):
        FPR,TPR = perf[binsize]['FPR'],perf[binsize]['TPR']
        AUC = sum([0.5*(TPR[k+1]+TPR[k])*(FPR[k]-FPR[k+1]) for k in range(len(FPR)-1)]) #the lists are in decreasing order
        color_i = color_list[i]
        axs.plot(FPR,TPR,label='{}: {}'.format(binsize,round(AUC, 2)), color=color_i)
        lower, upper = perf[binsize]['TPR_lower'],perf[binsize]['TPR_higher']
        axs.fill_between(FPR, lower, upper, color=color_i, alpha=0.1)
        # save curves

        output_plot = os.path.join(dataset, "ROC_curve_data/")
        if not os.path.exists(output_plot):
            os.makedirs(output_plot)

        # Curve
        pd.DataFrame(data={"X":FPR, "Y":TPR}).to_csv(os.path.join(output_plot, "chr{}_event_{}_curve.csv".format(truth_regions[region]['chr'], binsize)), index=None)
        # Upper curve
        pd.DataFrame(data={"X":FPR, "Y":upper}).to_csv(os.path.join(output_plot, "chr{}_event_{}_upper.csv".format(truth_regions[region]['chr'], binsize)), index=None)
        # Lower curve
        pd.DataFrame(data={"X":FPR, "Y":lower}).to_csv(os.path.join(output_plot, "chr{}_event_{}_lower.csv".format(truth_regions[region]['chr'], binsize)), index=None)


    #plot the f(x)=x line in '- - -'
    axs.plot([0,1],[0,1],'--',color='k')
    axs.set_ylim([0,1.01])
    axs.set_xlim([0,1])  
    axs.legend(title='AUC:')
    axs.set_xlabel('FPR')
    axs.set_ylabel('TPR')
    axs.set_title('{}, Region: chr{}, {:.2E} - {:.2E}'.format(dataset,truth_regions[region]['chr'], truth_regions[region]['bounds'][0], truth_regions[region]['bounds'][1]), fontsize=10)
    axs.grid('off')
    
    fig.tight_layout()

    plt.savefig(save_dir+'{}_AUC_perf_chr{}_filter_SNPS_{}.pdf'.format(dataset,truth_regions[region]['chr'],filter_SNPS))
    plt.show()


# In[ ]:




