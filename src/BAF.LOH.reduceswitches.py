import numpy as np
import pandas as pd
import pickle
import os
import sys

# Make sure I am running python2
print(sys.version)

# Load variables
PATIENT=sys.argv[1]
ref=sys.argv[2]
LOCATION=sys.argv[3]
PREPROCESSING=sys.argv[4]
bin_size=sys.argv[5]
sliding=sys.argv[6]
switch=sys.argv[7]
print(PATIENT)
print(ref)
print(LOCATION)
print(PREPROCESSING)
print(bin_size)
dataset="{}_{}".format(PATIENT,ref)
MAINDIR=LOCATION


# 
chr_list = range(1,23)
filter_SNPS = [1, 3, 10]
# filter_counts = 4*filter_SNPS #if one copy is deleted the expected coverage on each SNP is 10 / 2 : allow 40...
thresh = 0.37 # Threshold for calling imbalance

dir_cell_names = "{}cell_names_indexes_{}.tsv".format(PREPROCESSING,dataset)
print(dir_cell_names)
df_names = pd.read_csv(dir_cell_names, sep='\t', header=None).iloc[:,1:3]
df_names.columns=['INDEX', 'NAME']
#names are sorted alphabetically in df, sort with respect to name to create the correspondance
df_names = df_names.sort_values(by='NAME',axis=0)
cell_names = df_names['NAME'].values.tolist()
print(cell_names)

# If using /gpfs/commons/groups/landau_lab/tprieto/jlaval/PreprocessingResults/
#if PATIENT=='4295':
#    cell_names = [cell for cell in cell_names if cell not in ['4295_4272T1', '4295_NonB']]
#    cell_names = ['4295_hg38_' + cell.split('_')[-1] for cell in cell_names if cell not in ['4295_4272T1']]
#    cell_names = cell_names + ['4295_hg38_UndeterminedBarcode']


open_dir = '{}BAF_{}_sliding{}_switch{}/'.format(MAINDIR, bin_size, sliding, switch)
save_dir = open_dir

#first, stack the BAF in a single matrix
#print(dataset, binsize)
BAF_list = []
SNPS_list = []
COUNTS_list = []
CHROMSOME_list = []
POSITION_list = []
CHROMOSOME_list = []
for chr in chr_list:
    
    print(chr), #in Python2, the coma is used to print different iterations on the same line
    with open(open_dir+'dict_chr{}.pkl'.format(chr), 'rb') as f:
        # print(f)
        #dict_BAF_chr = pickle.load(f, encoding='latin1')	
        dict_BAF_chr = pickle.load(f)

    if len(dict_BAF_chr[cell_names[0]]['BAF'][0])==1:
        BAF = np.swapaxes(np.hstack(np.array(dict_BAF_chr[cell]['BAF']) for cell in cell_names), 0, 1)
        SNPS = np.swapaxes(np.hstack(np.array(dict_BAF_chr[cell]['SNP']) for cell in cell_names), 0, 1)
        COUNTS = np.swapaxes(np.hstack(np.array(dict_BAF_chr[cell]['COUNTS']) for cell in cell_names), 0, 1)
        #keep track of the chromosome and position
        POSITION=float(bin_size)*np.arange(BAF.shape[1])
        CHROMOSOME=np.array([chr for k in range(BAF.shape[1])])

    else:
        BAF = np.vstack([dict_BAF_chr[cell]['BAF'] for cell in cell_names])
        SNPS = np.vstack([dict_BAF_chr[cell]['SNP'] for cell in cell_names])
        COUNTS = np.vstack([dict_BAF_chr[cell]['COUNTS'] for cell in cell_names])
        #keep track of the chromosome and position
        POSITION=float(bin_size)*np.arange(BAF.shape[1])
        CHROMOSOME=np.array([chr for k in range(BAF.shape[1])])
    
    BAF_list.append(BAF)
    SNPS_list.append(SNPS)
    COUNTS_list.append(COUNTS)
    POSITION_list.append(POSITION)
    CHROMOSOME_list.append(CHROMOSOME)

BAF = np.hstack(BAF_list)
SNPS = np.hstack(SNPS_list)
COUNTS = np.hstack(COUNTS_list)
POSITION=np.hstack(POSITION_list)
CHROMOSOME=np.hstack(CHROMOSOME_list)




for filter_i in filter_SNPS:

    filter_counts = 4*filter_i #if one copy is deleted the expected coverage on each SNP is 10 / 2 : allow 40...
    
    #apply filters and select a unique BAF value per bin
    BAF_mat = np.zeros(BAF.shape)
    for i in range(BAF.shape[0]):
        for j in range(BAF.shape[1]):
            
            baf_slides = BAF[i,j]
            snps_slides = SNPS[i,j]
            counts_slides = COUNTS[i,j]
            
            if type(snps_slides)==list:
                # filter by number of SNPS and counts
                for k,(snp,count) in enumerate(zip(snps_slides,counts_slides)):
                    if snp < filter_i or count < filter_counts:
                        # baf_slides[k]=0.5 ##ISSUE: THE BAF ARRAY ITSELF IS MODIFIED, FIND A COPY METHOD IN PYTHON 2
                        baf_slides[k]=np.nan

                #take the maximum imbalance
                idx=np.argmax(np.abs(np.array(baf_slides)-0.5))
                baf = baf_slides[idx]

            else:
                if snps_slides < filter_i or counts_slides < filter_counts:
                    baf_slides=np.nan
                baf = baf_slides
            
            BAF_mat[i,j]=baf


    #Save the BAF
    dict_BAF = {cell: BAF_mat[k] for k,cell in enumerate(cell_names)}
    dict_BAF['CHR']=CHROMOSOME
    dict_BAF['POS']=POSITION

    print(save_dir)
    with open(save_dir+'BAF_numSNPS_' + str(filter_i) + '.pkl', 'wb') as f:
        pickle.dump(dict_BAF, f)


BAF_mat = np.zeros(BAF.shape)
for i in range(BAF.shape[0]):
    for j in range(BAF.shape[1]):
        
        baf_slides = BAF[i,j]
        snps_slides = SNPS[i,j]
        counts_slides = COUNTS[i,j]
        
        if type(snps_slides)==list:
            #take the maximum imbalance
            idx=np.argmax(np.abs(np.array(baf_slides)-0.5))
            baf = baf_slides[idx]

        else:
            baf = baf_slides
        
        BAF_mat[i,j]=baf

#Save the BAF
dict_BAF = {cell: BAF_mat[k] for k,cell in enumerate(cell_names)}
dict_BAF['CHR']=CHROMOSOME
dict_BAF['POS']=POSITION

print(save_dir)
with open(save_dir+'BAF.pkl', 'wb') as f:
    pickle.dump(dict_BAF, f)


#open the BAF
with open(open_dir+'BAF.pkl', 'rb') as f:
    dict_BAF = pickle.load(f)

BAF_mat = np.vstack([dict_BAF[cell] for cell in cell_names])

#apply threshold for CNVs 
CNVs = np.abs(BAF_mat-0.5) > thresh
#and save it
dict_CNVs = {cell: CNVs[k] for k,cell in enumerate(cell_names)}
dict_CNVs['CHR']=CHROMOSOME
dict_CNVs['POS']=POSITION
with open(save_dir+'CNVs.pkl', 'wb') as f:
    pickle.dump(dict_CNVs, f)
#also save it as csv
df_CNVs = pd.DataFrame.from_dict(dict_CNVs, orient='columns')
df_CNVs.to_csv(save_dir+'CNVs.tsv',sep='\t',index=None)



# We threshold the BAF values to get allele-specific calls
BAF_mat[(0.5-thresh >= BAF_mat) & (BAF_mat >= 0.0)] = 0.0
BAF_mat[(0.5+thresh >= BAF_mat) & (BAF_mat >= 0.5-thresh)] = 0.5
BAF_mat[(1.0 >= BAF_mat) & (BAF_mat >= 0.5+thresh)] = 1.0

#Save the BAF calls
dict_BAF = {cell: BAF_mat[k] for k,cell in enumerate(cell_names)}
dict_BAF['CHR']=CHROMOSOME
dict_BAF['POS']=POSITION

with open(save_dir+'CNVs_allele_specific.pkl', 'wb') as f:
    pickle.dump(dict_BAF, f)
#also save it as csv
df_BAF = pd.DataFrame.from_dict(dict_BAF, orient='columns')
df_BAF.to_csv(save_dir+'CNVs_allele_specific.tsv',sep='\t',index=None)
save_dir
