#!/usr/bin/env Rscript
library(ape)
library(tidyr)
library(dplyr)
library(tibble)
library(data.table)
library(ggtree)

# get the arguments
args = commandArgs(trailingOnly=TRUE)
#patient <- "4295"
patient <- args[1]

main_dir <- "/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/2021-07-19_BioSkryb/"

GetExomeVariants <- function(patient,
                             mindepth4vaf=8,
                             distance2indel=10,
                             average_vaf_threshold=0.3,
                             variance_vaf_threshold=0.1,
                             perc_missing_threshold=0.4) {

  if (patient=="Invitro"){
      mydir <- paste0("/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/",
                   patient)
  } else {
      mydir <- paste0("/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/WGS/",patient,"/")
  }
  
  var_table <- data.table::fread(file = paste0(main_dir,"SomaticCalls/","Mutect",patient,".annotations.bed"))

# do not calculate average including bulk because VAF is not expected to be 0.5 in bulk (subclonal)
  bulk_samples <- c("445-380T","445-445T", 
                    "4292-4295T1","4292-4295T2","4292-4295NonB", 
                    "4292-4272T1","4292-4272T2","4295-4272B",
                    "417-368B","417-368T","417-417T",
                    "4084-4072T","4084-4072NonB","4084-4084CD19","4084-4084T")


  vaf_matrix <- var_table %>%
    # Remove positions which do not PASS or are clustered_events
    dplyr::filter(FILTER %in% c("PASS","clustered_events")) %>%
    # there are several VAFs for the same positions in some cells, we will keep just the first one (after the comma)
    dplyr::mutate(across(all_of(matches(":VAF")), ~ as.numeric(gsub(",.*","",.)))) %>%
    tidyr::pivot_longer(cols = matches(":DP|:VAF"),
                          names_to="sample",values_to="value") %>%
    tidyr::separate(col=sample, into=c("sample","metric"),sep=":") %>%
    dplyr::rename(dp_total=DP) %>%
    tidyr::pivot_wider(names_from = metric, values_from=value,values_fill = NA) %>%
    dplyr::mutate(type=ifelse(nchar(REF)==1 & nchar(ALT)==1,"snv","indel")) %>%
    # Rename samples
    dplyr::mutate(sample=gsub("MRD-BALL-PTA-NEXTERA-Exome-","4084-",gsub("_S.*","",sample))) %>%
    # If the coverage is lower than 8 then do not trust the VAF value
    dplyr::mutate(VAF=ifelse(DP<=mindepth4vaf,NA,VAF)) %>% 
    dplyr::group_by(`CHROM:POS`) %>%
    dplyr::mutate(vaf_pres=ifelse(VAF>0 & !(sample %in% bulk_samples),VAF,NA)) %>%
    # calculate the mean vaf of cells carrying a variant
    dplyr::mutate(average_vaf=mean(as.numeric(vaf_pres), na.rm = T)) %>% 
    # calculate the variance of the vaf in cells carrying the variant
    dplyr::mutate(variance_vaf=var(as.numeric(vaf_pres), na.rm = T)) #%>% 
 

  # Sort columns for the heatmap (cells)
  tree_distance <- readRDS(
    file = paste0(mydir,"RESULTS/isec_somatic-germline/FILTERED/MutationMapping/TreeMut.genome.rds"))  
  tree_distance@phylo$tip.label
  non_classified <- setdiff(gsub("-","_",gsub("_.*","",gsub(":VAF","",vaf_matrix$sample))),tree_distance@phylo$tip.label)
  ordered_tree_labels <- gsub("_","-",ggtree::get_taxa_name(ggtree(tree_distance@phylo)))
  all_cells <- c(rev(ordered_tree_labels),gsub("_","-",non_classified))


 output <-  vaf_matrix %>%
    dplyr::mutate(sequencing=ifelse(sample %in% gsub("_","-",non_classified), "not in tree","in tree")) %>%
    dplyr::mutate(pos=as.numeric(gsub(".*:","",`CHROM:POS`))) %>%
    dplyr::mutate(chrom=gsub(":.*","",`CHROM:POS`)) %>%
    dplyr::group_by(`CHROM:POS`) %>%
    # Remove variants close to an indel
    dplyr::mutate(dist_close_lag=ifelse(chrom==lag(chrom),pos-lag(pos),NA)) %>%
    dplyr::mutate(dist_close_lead=ifelse(chrom==lead(chrom),pos-lead(pos),NA)) %>%
   # when there are NAs the filtering does not work well
       dplyr::mutate(dist_close_lag=ifelse(is.nan(dist_close_lag) | is.na(dist_close_lag), as.numeric(distance2indel+1), dist_close_lag)) %>%
         dplyr::mutate(dist_close_lead=ifelse(is.nan(dist_close_lead) | is.na(dist_close_lead), as.numeric(distance2indel+1), dist_close_lead)) %>%
    # Remove clustered events close to an indel
     dplyr::filter(!(dist_close_lag<=distance2indel & lag(type)=="indel")) %>%
     dplyr::filter(!(dist_close_lead<=distance2indel & lead(type)=="indel")) %>%
    # Remove indels
      dplyr::filter(type=="snv") %>%
      # Remove variants with a mean VAF lower than 25% (could be damage). 30% to make sure
      dplyr::filter(average_vaf>average_vaf_threshold) %>%
      dplyr::filter(variance_vaf<variance_vaf_threshold) %>% # it is exome, the variance in vaf should not be high
    # Remove mutations which show a 1 vaf (somatic should be HET)
        dplyr::filter(average_vaf<1) %>%
    # remove variants with more than 50% missing data
    dplyr::mutate(missing=ifelse(is.na(VAF),1,0)) %>%
    dplyr::mutate(cell=1) %>%
   dplyr::mutate(perc_missing=sum(missing)/sum(cell)) %>%
    dplyr::filter(perc_missing<perc_missing_threshold) %>%
      dplyr::mutate(sample=factor(as.factor(sample),levels=all_cells)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(id=paste0(Gene.refGene,"_",gsub(":NM_.*","",gsub(".*c.","c.",AAChange.refGene)))) 
 
 # Sort rows in the heatmap 
 output$id <- reorder(output$id, -output$average_vaf)
 
  saveRDS(output, paste0(main_dir,"SomaticCalls/","Mutect",patient,".filtered.rds"))

  output %>%
    dplyr::select(`CHROM:POS`,REF,ALT,id,type) %>%
    unique() %>%
    tidyr::separate(`CHROM:POS`,into=c("chr","pos"),sep=":") %>%
    dplyr::mutate(start=as.numeric(pos)-1) %>%
    dplyr::mutate(info=gsub("_\\.$","_noncoding",gsub(":p.","_p.",paste0(chr,":",pos,REF,">",ALT,"_",id)))) %>%
    dplyr::select(chr,start,pos,info,type) %>%
    dplyr::arrange(chr,start) %>%
  write.table(paste0(main_dir,"SomaticCalls/","Mutect",patient,".filtered.bed"),
              sep = "\t", quote = F, row.names = F,col.names = F)
    
}


GetExomeVariants(patient)


