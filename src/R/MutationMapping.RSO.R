#!/usr/bin/env Rscript

# This script requires a phylogeny
# AD counts for the sites
# Cell metadata


library(treemut)
library(ape)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggtree)
library(data.table)

# treemut was designed to enable users to map clonal somatic mutations to a specified tree.  The
# method requires:
#   - phylogeny
#   - the mutant read counts and depth for each sample/loci in the form of matrices
#   - error rate per single cell


# my tests

patient <- "eso29"
# load the tree (old chr1)
mydir <- paste0("/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/",patient,"/results//MutationMapping/")
treename <- "/../CellPhy-01/CellPhy.GT10+FO+E.raxml.bestTree"
# load the tree (deepvar tree)
mydir <- "/gpfs/commons/groups/landau_lab/jzinno/ultima/eso-ug/annovar/cellphy-eso/"
treename <- paste0(patient,".annofilt_dp.CellPhy.raxml.support")
# File containing AD counts 
 
####################################################
# Obtain the patient name and phylogeny directory #
####################################################

args = commandArgs(trailingOnly=TRUE)
patient <- args[1]
mydir <- args[2]
treename <- paste0(patient,".annofilt_dp.CellPhy.raxml.support")

# Rooting messes up the mapping. Right now is working, be careful if modifying
MutationMappingTreeMut <- function(tree, mtr, depth, myoutgroup="t17",
                                   cross_sample_contamination=0,
                                   pval_thres=0.05, maxIt_pval=5,
                                   outputdir=".", filename="TreeMut.pdf"){


      #tree <- ape::drop.tip(tree,"zeros")
      #tree <- ape::root(tree, myoutgroup)
      tree=ape::bind.tree(tree,ape::read.tree(text="(zeros:0);"))
      tree<-multi2di(tree)
      tree <- ape::root(tree, "zeros")
      # Alternative below (it does not reroot the tree, do not reroot it on zeros)
      #tree <- phytools::bind.tip(tree, tip.label="zeros", edge.length=0, where=castor::find_root(tree), position=0)

      # Create node names (very important for downstream processing)
      tree <- ape::makeNodeLabel(tree, method = "number", prefix = "N")

      # FUNCTION TO MAP THE MUTATIONS TO THE TREE
      res=treemut::assign_to_tree(tree,mtr,depth,error_rate=cross_sample_contamination, maxits = maxIt_pval)

      # Code to map the mutations to the tree and visualize it
      mydataRaw <- cbind(res$tree$edge, original.branch.length=tree$edge.length,
                      res$df$df) %>%
        dplyr::full_join(cbind(res$summary,mutid=rownames(mtr)) %>%
                           dplyr::mutate(all=ifelse(pval>pval_thres,1,1)) %>%
        # pval: a heuristic pvalue assessing the hypothesis that the mutation is consistent with the provided tree topology
          # null hypothesis: pval>0.05 -> consistent
          # reject hypothesis: pval<0.05 -> inconsistent
                           dplyr::mutate(signif=ifelse(as.numeric(pval)>=pval_thres | as.numeric(p_else_where)==0,1,0)) %>%
                           dplyr::rename(id=edge_ml), by="id") %>%
        # if pval close to 1 then create another list of mutations
        dplyr::mutate(mutid2=ifelse(signif==1,mutid,NA)) %>%
        dplyr::rename(label2=label) %>%
        dplyr::rename(parent=`1`) %>%
        dplyr::rename(node=`2`) %>%
        dplyr::left_join(res$tree %>% as_tibble(), by=c("node","parent"))

      vaf=mtr/depth
      
      
      # Select a few sites to annotate the tree in case there are too many
      # Try to select those with intermediate allele frequencies
      if (length(unique(mydataRaw$mutid2))>500){
        # select sites based on frequency
        sites_to_plot <- unique(mydataRaw$mutid2[!is.na(mydataRaw$mutid2) &
                                                   # Remove singletons from the plot
                                                   as.numeric(ifelse(is.na(mydataRaw$mutid2),0,(gsub(".*_","",mydataRaw$mutid2))))>1 &
                                                   # Remove potential germline SNPs
                                                   #as.numeric(ifelse(is.na(mydataRaw$mutid2),0,(gsub(".*_","",mydataRaw$mutid2))))<as.numeric(length(res$tree$tip.label)-10) &
                                                   #as.vector(rowMeans(vaf, na.rm=T)>0.005)  &
                                                   mydataRaw$mutid2!="UNKNOWN" #&
                                                   #as.vector(rowMeans(vaf, na.rm=T)<0.6)
                                                 ])
        length(sites_to_plot)
        sites_to_plot_pre <- sites_to_plot
        # Remove number of cells from the name
        rownames(vaf)=gsub("_.*","",rownames(vaf))
        sites_to_plot <- gsub("_.*","",sites_to_plot)
        # select 20 random sites
        #sites_to_plot <- mydataRaw$mutid2[sample(x = which(!is.na(mydataRaw$mutid2) & mydataRaw$mut_count>3),size = 20)]
      } else {
        sites_to_plot <- mydataRaw$mutid2[which(!is.na(mydataRaw$mutid2))]
        sites_to_plot_pre <- sites_to_plot
      }
      
      
      mydata <- mydataRaw %>%
        # Create a third list with a few mutations
        dplyr::mutate(mutid3=ifelse(mutid2 %in% sites_to_plot_pre,mutid2,NA)) %>%
        # Create a list of mutations to map to each node
        dplyr::group_by(parent,node,label) %>%
        dplyr::mutate(mutList = gsub(",$","",gsub("NA","",gsub("NA,","",paste(mutid2, collapse = ","))))) %>%
        dplyr::mutate(mutList=ifelse(mutList=="NA",NA,mutList)) %>%
        dplyr::mutate(mutList=ifelse(mutList=="",NA,mutList)) %>%
        dplyr::mutate(mutList2 = gsub(",$","",gsub("NA","",gsub("NA,","",paste(mutid3, collapse = ","))))) %>%
        dplyr::mutate(mutList2=ifelse(mutList2=="NA",NA,mutList2)) %>%
        dplyr::mutate(mutList2=ifelse(mutList2=="",NA,mutList2))

      mydatasummarizedperbranch <- mydata %>%
        # remove non-unique columns to obtain branches alone
        dplyr::select(-pval,-p_else_where,-mutid,-mutid2,-mutid3,-signif,-all) %>% unique()

      newtree <- res$tree %>%
        tibble::as_tibble() %>%
        dplyr::left_join(mydatasummarizedperbranch, by=c("parent","node","branch.length","label")) %>%
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # OPTIONAL!!!Replace new branch length by original branch length
        #dplyr::mutate(branch.length=original.branch.length) %>%
        treeio::as.treedata()
      tree_withzeros <- newtree

      # Remove zeros and related data
      #newtree <- treeio::drop.tip(newtree, tip = "zeros")

      
      newtree@phylo$tip.label<- gsub("zeros","Germline",newtree@phylo$tip.label)
      tree_plot_new <- ggtree(newtree) +
        geom_tiplab(size=3) +
        geom_rootedge(1000) +
        geom_label(aes(x = branch, label=mutList2), size=1.5, hjust=0)# +
        #theme_tree2()
      
      vaf = cbind(vaf, Germline = 0)
      vaf_matrix_toplot <- t(vaf)[,sites_to_plot]

      mygheatmap <- gheatmap(tree_plot_new, vaf_matrix_toplot,
                             offset=0.09, # distance from tips to the matrix
                             #offset=20000, # distance from tips to the matrix
                             width=10, # default 1
                             colnames_angle=90,
                             hjust = 1, # centering colnames
                             colnames_offset_y = 0.4,
                             #colnames_position = "bottom",
                             font.size = 1, # matrix colnames
                             color="black",
                             #legend_title="",
                             colnames=T) +
        scale_fill_gradient2(low = "grey",mid="cornflowerblue", high = "darkred", na.value = "lightgrey")
      ggsave(plot = mygheatmap + vexpand(.15, direction = -1), 
             filename = paste0(outputdir,sep="/",filename), dpi=300, width = 15, height = 8,limitsize = FALSE)


      return(list(tree_withzeros=tree_withzeros,tree=newtree, mutations=mydata %>% dplyr::select(mutid,pval,p_else_where,mutList,mutList2) %>% dplyr::filter(!is.na(mutid))))
}


MapMutations <- function(patient, suffix2=""){

      # Load the tree
      tree_final <- read.tree(paste(mydir,treename,sep=""))
  
      suffix=".genome"
      # Load the count matrices for the genome
      ad_depth_table <- data.table::fread(paste0(mydir,"AD.bed"), header = T)  %>%
        tidyr::pivot_longer(cols = matches(paste0(patient,'_'))) %>%
        tidyr::separate(col = value,into=c("ref","alt"), convert=T) %>%
        dplyr::mutate(cov=as.numeric(ref+alt)) %>%
        dplyr::mutate(ref = ifelse(ref==0 & alt==0, NA,ref)) %>% #avoid division of 0 by 0
        dplyr::mutate(vaf=as.numeric(alt/(ref+alt)),vaf) %>%
        # dplyr::mutate(cov=as.numeric(ref+alt,)) %>%
        # number of cells in which the allele frequency is higher than 0
        dplyr::mutate(presabs=ifelse(vaf>0,1,0)) %>%
        # number of cells in which the allele frequency is not NA
        dplyr::mutate(cov_sc=ifelse(is.na(vaf),NA,1)) %>%
        dplyr::group_by(`CHROM:POS`) %>%
        dplyr::mutate(NumCells=sum(presabs, na.rm=T))  %>%
        dplyr::mutate(NumCellsWithCov=sum(cov_sc, na.rm=T))  %>%
        dplyr::select(-presabs,-cov_sc) %>%
        ungroup()
       
       # Bed of the genomic mutations
	     tag="VariantSitesForPhylogeny.informative.all"

       bed <- data.table::fread(paste(mydir,"../",tag,".bed",sep=""), col.names = c("CHROM","START","END","Type","Gene","Ogaini","num_cell_alt"), sep=" ", na.strings = "NA") %>%
        # add a counter to each gene (same gene with several mutations)
        group_by(Gene,Type) %>% mutate(counter = row_number(Gene)) %>%
        dplyr::ungroup() %>%   
        dplyr::mutate(Gene=ifelse(!is.na(Ogaini),paste0(Gene,"-ogacor",Ogaini),Gene)) %>%
         dplyr::mutate(Gene=ifelse(!is.na(Gene),paste0(Gene,"_",num_cell_alt),Gene)) %>%
        dplyr::mutate(pos=paste0(CHROM,":",END)) %>%
         dplyr::mutate(Gene=gsub(".*:","",Gene)) %>%
        dplyr::select(pos,Type,Gene)

        
      mtr <- ad_depth_table %>%
        #dplyr::filter(NumCells>1) %>%
        dplyr::select(`CHROM:POS`,name,alt) %>%
        tidyr::pivot_wider(names_from = name, values_from=alt) %>%
        tibble::column_to_rownames("CHROM:POS")  #%>% head(100)
      depth <- ad_depth_table %>%
        #dplyr::filter(NumCells>1) %>%
        dplyr::select(`CHROM:POS`,name,cov) %>%
        tidyr::pivot_wider(names_from = name, values_from=cov) %>%
        tibble::column_to_rownames("CHROM:POS") #%>% head(100)
      # Transform into matrices (avoid issues later on)
      mtr <- data.matrix(mtr)
      depth <- data.matrix(depth)
      # change mutation names
      # The name of the mutations annotated depends on the matrix mtr rownames
      rownames(mtr)<- bed$Gene  # They have to be unique, non NA values allowed
      rownames(depth)<- bed$Gene        

      # Obtain cross-sample contamination percentage
      # I think it has to be in the order of the counts and alt alleles
      # What it does the percentage is reducing the allele frequency at which the alternative allele might be found
      crossample_cont <- 
          cbind(Sample=colnames(mtr),in_tree=rep(1,length(colnames(mtr)))) %>% as.data.frame() %>%
          dplyr::full_join(data.table::fread("/gpfs/commons/groups/landau_lab/ResolveOME/StartDir/eso29/metadata/metadata.txt", header = T) %>%
                             dplyr::select(`sample name`, `contamination prop`) %>%
                             dplyr::mutate(Sample=`sample name`)) %>%
          dplyr::filter(in_tree==1) %>%
          rowwise() %>%
          dplyr::mutate(crosscont=if_else(is.na(`contamination prop`),0,`contamination prop`)) %>%
          dplyr::mutate(crossample_cont=`contamination prop`) %>%
        dplyr::mutate(crossample_cont=round(crossample_cont, digits = 2)) %>%
        dplyr::select(crossample_cont) %>% unlist()
      
      # MUTATIONS ARE ONLY ASSIGNED ONCE (hardassigned after pval calculations), get NEW LENGTHS
     results <- MutationMappingTreeMut(tree = tree_final, 
                            mtr = mtr, depth = depth,
                            myoutgroup = "zeros", 
                            outputdir = mydir,
                            filename = paste0("TreeMut",suffix,".pdf"),
                            # Error rate. This is intended as a coarse error measure and should not take more than 4 distinct values. Useful if a sample is mildly contaminated error->0.1
                            cross_sample_contamination = crossample_cont,
                            #cross_sample_contamination = 0,
                            pval_thres = 0.05,  # I am also selecting all sites for which the hypothesis of being consistent is rejected (pval<0.05) but with 0 probability anywhere else in the phylogeny (p_else_where=0)
                            maxIt_pval=100
                            )
    results$tree
    saveRDS(object = results$tree_withzeros, file = paste0(mydir,
            "TreeMutWithZeros.genome.rds"))    
    saveRDS(object = results$tree, file = paste0(mydir,
            "TreeMut.genome.rds"))
    results$mutations
}


# I need to fix this (probably mapping again)
CreateBedFileMutationInPhylogeny <- function (patient, suffix=".genome"){
  

  print ("Writing bed file...")
  tree_distance <- readRDS(file = paste0(mydir,
            "TreeMut.genome.rds"))
  unlist(strsplit(tree_distance@data$mutList,split = ",")) %>%
    as.data.frame() %>%
    cbind(branch_node=rep(x=tree_distance@data$label2,  times=unlist(lapply(strsplit(tree_distance@data$mutList,split = ","), length)))) %>%
    dplyr::rename("chrpos"=".") %>%
    tidyr::separate(col = "chrpos",into=c("CHROM","END"), convert=T) %>%
    dplyr::mutate(POS=END-1) %>%
    dplyr::select(CHROM,POS,END,branch_node) %>%
    dplyr::filter(!is.na(CHROM)) %>%
    write.table(
                  file=paste0(mydir,"VariantSitesCellPhyAfterMutationMapping",suffix,".bed"),
                       quote = F, row.names = F,
                       col.names = F, sep = "\t")
}

patient <- "eso29"
print("> Analyzing genome")
MapMutations(patient = patient)
print("> Create bed file genome")

# Any branch with a length of 0 after SNV fitting should be collapsed into a polytomy...??? implemented in Coorens et al. 
CreateBedFileMutationInPhylogeny(patient=patient)
#print("Analyzing exome")
#MapMutations(patient = patient,genome = F)
#print("> Create bed file exome")
#CreateBedFileMutationInPhylogeny(patient=patient,suffix=".exome")


