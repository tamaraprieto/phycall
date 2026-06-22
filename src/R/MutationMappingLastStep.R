#!/usr/bin/env Rscript
library(treemut)
library(ape)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggtree)
library(data.table)

patient <- "4295"
# get the arguments
#args = commandArgs(trailingOnly=TRUE)
#patient <- args[1]

#############################################
# Function to map the mutations to the tree #
# called below                              #
#############################################

MutationMappingTreeMut <- function(tree, mtr, depth, myoutgroup="t17",
                                   cross_sample_contamination=0,
                                   pval_thres=0.05, pval_elsewhere=-0.05, maxIt_pval=5){


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
      mydataRaw <- cbind(res$tree$edge, original.branch.length=tree$edge.length, new.branch.length=res$tree$edge.length,
                      res$df$df) %>%
        dplyr::full_join(cbind(res$summary,mutid=rownames(mtr)) %>%
                           dplyr::mutate(all=ifelse(pval>pval_thres,1,1)) %>%
        # pval: a heuristic pvalue assessing the hypothesis that the mutation is consistent with the provided tree topology
          # null hypothesis: pval>0.05 -> consistent
          # reject hypothesis: pval<0.05 -> inconsistent
                           dplyr::mutate(signif=ifelse(as.numeric(pval)>=pval_thres | as.numeric(p_else_where)>=pval_elsewhere,1,0)) %>%
                           dplyr::rename(id=edge_ml), by="id") %>%
        # if pval close to 1 then create another list of mutations
        dplyr::mutate(mutid2=ifelse(signif==1,mutid,NA)) %>%
        dplyr::rename(label2=label) %>%
        #dplyr::rename(parent=`1`) %>%
        #dplyr::rename(node=`2`) %>%
        dplyr::left_join(res$tree %>% as_tibble(), by=c("node","parent"))

      # Select a few sites to annotate the tree in case there are too many
      #if (length(unique(mydataRaw$mutid2))>100){
      #  sites_to_plot <- mydataRaw$mutid2[sample(x = which(!is.na(mydataRaw$mutid2) & mydataRaw$mut_count>3),size = 20)]
      #} else {
        sites_to_plot <- mydataRaw$mutid2[which(!is.na(mydataRaw$mutid2))]
      #}

      mydata <- mydataRaw %>%
        # Create a third list with a few mutations
        dplyr::mutate(mutid3=ifelse(mutid2 %in% sites_to_plot,mutid2,NA)) %>%
        # Create a list of mutations to map to each node
        dplyr::group_by(parent,node,label) %>%
        dplyr::mutate(mutList = gsub(",$","",gsub("NA","",gsub("NA,","",paste(mutid2, collapse = "\n"))))) %>%
        dplyr::mutate(mutList=ifelse(mutList=="NA",NA,mutList)) %>%
        dplyr::mutate(mutList=ifelse(mutList=="",NA,mutList)) %>%
        dplyr::mutate(mutList2 = gsub(",$","",gsub("NA","",gsub("NA,","",paste(mutid3, collapse = "\n"))))) %>%
        dplyr::mutate(mutList2=ifelse(mutList2=="NA",NA,mutList2)) %>%
        dplyr::mutate(mutList2=ifelse(mutList2=="",NA,mutList2))

      mydatasummarizedperbranch <- mydata %>%
        # remove non-unique columns to obtain branches alone
        dplyr::select(-pval,-p_else_where,-mutid,-mutid2,-mutid3,-signif,-all) %>% unique()

      newtree <- res$tree %>%
        tibble::as_tibble() %>%
        dplyr::left_join(mydatasummarizedperbranch, by=c("parent","node","branch.length","label")) %>%
        dplyr::mutate(branch.length=original.branch.length) %>%
        treeio::as.treedata()
      tree_withzeros <- newtree

      # Remove zeros and related data
      #newtree <- treeio::drop.tip(newtree, tip = "zeros")

      vaf <- mtr/depth
      vaf = cbind(vaf, zeros = 0)

      return(list(tree_withzeros=tree_withzeros,
                  tree=newtree,
                  vafmatrix=t(vaf)[,sites_to_plot], 
                  mutations=mydata %>% 
                    dplyr::select(mutid,pval,p_else_where,mutList,mutList2) %>% 
                    dplyr::filter(!is.na(mutid))))
}


######################
# Load tree and data #
######################


if (patient=="Invitro"){
  mydir=paste0("/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/",patient)
} else {
  mydir=paste0("/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/WGS/",patient)
}
  
# Load the tree with accurate branch lengths (after mutation mapping)
 tree_final <- readRDS(file = paste0(mydir,"/RESULTS/isec_somatic-germline/FILTERED/MutationMapping/",
                    "TreeMut.genome.rds"))
 tree_final <- tree_final@phylo

 

      ###################
      # Set up outgroup #
      ################### 
     
      #if (patient=="4295"){
      #  outgroup="4295_A11"
      #} else if (patient=="417"){
      #  outgroup <- "417_B10"
      #} else if (patient=="4084"){
      #  outgroup="4084_F1"
      #} else if (patient=="445") {
      #  outgroup="445_A3"
      #} else if (patient=="Invitro"){
      #  outgroup="Invitro_A1"
      #} else {
      #  outgroup=""
      #}

      ######## CNVs #########
      interval <- "50000" #10kb
      # maybe map with parsimony the 0s and 1s?? give it a try. I might need to correct dropout in some way
      #cnvs <- fread(paste0(mydir,"/RESULTS/BAF/BAF_",interval,"/CNVs.collapsedbycytogeneticbandtp.tsv"))
      # transform to 0.5 if not...
      #######################
      
      
      # the counts of the reference and alternative reads I think
      # is better to get them from the bam file than from the vcf for the exome
     #!!!! I have to recalculate them with the new calling
      alleliccounts <- data.table::fread(paste0("/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/2021-07-19_BioSkryb/CollectAllelicCountsExomeCalls/AllelicCounts.",patient,".txt"), col.names = c("chrom","pos","ref","alt","REF","ALT","name"))
      #bed2 <- data.table::fread(paste0("/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/2021-07-19_BioSkryb/SomaticCalls/",patient,"_SNP_filtered.bed"), col.names = c("chrom","start","pos","gene","vartype")) %>%
      bed2 <- data.table::fread(paste0("/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/2021-07-19_BioSkryb/SomaticCalls/Mutect",patient,".filtered.bed"), col.names = c("chrom","start","pos","gene","vartype")) %>%
        dplyr::select(-start)
      ad_depth_table <- alleliccounts %>%
        dplyr::left_join(bed2) %>%
        dplyr::mutate("CHROM:POS"=paste(chrom,pos,sep=":")) %>%
        dplyr::mutate(cov=as.numeric(ref+alt)) %>%
        dplyr::mutate(vaf=as.numeric(alt/(ref+alt))) %>%
        # number of cells in which the allele frequency is higher than 0
        dplyr::mutate(presabs=ifelse(vaf>0,1,0)) %>%
        # number of cells in which the allele frequency is not NA
        dplyr::mutate(cov_sc=ifelse(is.na(vaf),NA,1)) %>%
        dplyr::group_by(`CHROM:POS`) %>%
        dplyr::mutate(NumCells=sum(presabs, na.rm=T))  %>%
        dplyr::mutate(NumCellsWithCov=sum(cov_sc, na.rm=T))  %>%
        dplyr::select(-presabs,-cov_sc) %>%
        ungroup() %>%
        dplyr::select(`CHROM:POS`,name,ref,alt,cov,vaf,NumCells,NumCellsWithCov,REF,ALT,gene) %>%
        dplyr::filter(name %in% tree_final$tip.label)
      bed <- bed2 %>%
        dplyr::mutate(pos=paste0(chrom,":",pos)) %>%
        dplyr::select(-chrom) %>%
        dplyr::rename(Type=vartype) %>%
        dplyr::rename(Gene=gene) %>%
        group_by(Gene) %>% mutate(counter = row_number(Gene)) %>%
        dplyr::ungroup() %>%   dplyr::mutate(Gene=ifelse(!is.na(Gene),paste(Gene,counter,sep="."),NA)) %>%
        dplyr::mutate(Gene=gsub(".1$","",Gene)) %>%
        dplyr::select(pos,Type,Gene)
      ad_depth_table <- ad_depth_table %>%
        dplyr::select(-REF,-ALT,-`CHROM:POS`) %>%
        # simplify gene name
        dplyr::mutate(gene=gsub("_noncoding","",gsub("_c.*_p."," ",gsub(".*>(A|C|T|G)_","",gene))))
        
      mtr <- ad_depth_table %>%
        #dplyr::filter(NumCells>1) %>%
        dplyr::select(gene,name,alt) %>%
        tidyr::pivot_wider(names_from = name, values_from=alt) %>%
        tibble::column_to_rownames("gene") #%>% head(100)
      depth <- ad_depth_table %>%
        #dplyr::filter(NumCells>1) %>%
        dplyr::select(gene,name,cov) %>%
        tidyr::pivot_wider(names_from = name, values_from=cov) %>%
        tibble::column_to_rownames("gene") #%>% head(100)      
      # Transform into matrices (avoid issues later on)
      mtr <- data.matrix(mtr)
      depth <- data.matrix(depth)
        

      # Obtain cross-sample contamination percentage
      # I think it has to be in the order of the counts and alt alleles
      # What it does the percentage is reducing the allele frequency at which the alternative allele might be found
      # !!!!!!CAREFUL!!!I think this assumes all mutations I am trying to map are heterozygous!!!
      crossample_cont <- 
          cbind(Sample=colnames(mtr),in_tree=rep(1,length(colnames(mtr)))) %>% as.data.frame() %>%
          dplyr::full_join(data.table::fread("https://docs.google.com/spreadsheets/d/e/2PACX-1vR91fISOmmJc4WDAsnx9LOXd80YirilNxPSoFBktYrMTWh8fjQpczb3tFH3QXrY30dMeSVD9_mlP4YB/pub?gid=0&single=true&output=tsv") %>%
                             dplyr::filter(Patient==patient) %>%
                             dplyr::select(Sample_tag, `Contamination % estimated from WGS`,`Contamination % estimated from WES or lowpass WGS`) %>%
                             dplyr::mutate(Sample=paste0(patient,"_",Sample_tag))) %>%
          dplyr::filter(in_tree==1) %>%
          rowwise() %>%
          dplyr::mutate(crosscont=max(`Contamination % estimated from WES or lowpass WGS`,`Contamination % estimated from WES or lowpass WGS`, na.rm = T)) %>%
          dplyr::mutate(crosscont=if_else(is.na(crosscont),0,crosscont)) %>%
          dplyr::mutate(crossample_cont=as.numeric(crosscont)/100) %>%
        dplyr::mutate(crossample_cont=round(crossample_cont, digits = 2)) %>%
          dplyr::select(crossample_cont) %>% unlist()
      

      # MUTATIONS ARE ONLY ASSIGNED ONCE (hardassigned after pval calculations)
      # I would actually use the mapping from different branches to map more than once
     results <- MutationMappingTreeMut(tree = tree_final, 
                            mtr = mtr, depth = depth,
                            #myoutgroup = "zeros", 
                            cross_sample_contamination = crossample_cont,
                            pval_thres = 0.05,  # I am also selecting all sites for which the hypothesis of being consistent is rejected (pval<0.05) but with 0 probability anywhere else in the phylogeny (p_else_where=0)
                            maxIt_pval=100
                            )
    results$tree
    saveRDS(object = results$tree_withzeros, file = paste0(mydir,"/RESULTS/isec_somatic-germline/FILTERED/MutationMapping/",
            "TreeMutWithZeros.genome.exomemapped.rds"))    
    saveRDS(object = results$tree, file = paste0(mydir,"/RESULTS/isec_somatic-germline/FILTERED/MutationMapping/",
            "TreeMut.genome.exomemapped.rds"))
 
    
    ###################################
    # Plot tree with mapped mutations #
    ###################################
    
    filename = paste0("TreeMut.genome.withexomemapped.",patient,".pdf")
    
    if (patient=="Invitro"){
      mydir=paste0("/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/",patient)
    } else {
      mydir=paste0("/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/WGS/",patient)
    }
    outputdir = paste0(
      mydir,
      "/RESULTS/isec_somatic-germline/FILTERED/MutationMapping/")
    
    
    results$tree@phylo <- ape::root(results$tree@phylo, paste0("zeros"))
    tree_plot_new <- ggtree(results$tree) +
      geom_tiplab(size=1.5) +
      geom_label(aes(x=branch, label=mutList2, col="darkred",fontface = 'bold', alpha = 0.2), size=1, angle=0, hjust=0.5,show_guide=F) +
      #theme_tree2() +
      labs(fill="VAF")
    
    mygheatmap <- gheatmap(tree_plot_new, results$vafmatrix,
                           #offset=2.8, # distance from tips to the matrix
                           offset=130, # distance from tips to the matrix
                           width=0.3, # default 1
                           colnames_angle=90,
                           hjust = 1, # centering colnames
                           colnames_offset_y = 0.4,
                           #colnames_position = "bottom",
                           font.size = 1.5, # matrix colnames
                           color="black",
                           #legend_title="",
                           colnames=T) +
      geom_rootedge(100) +
      scale_fill_gradient2(low = "grey",mid="skyblue2", high = "brown3", na.value = "lightgrey") +
      vexpand(.1, direction = -1) +
      theme(legend.position = 'top', legend.title.align = 0, 
            legend.title = element_text(size=6), 
            legend.text = element_text(size=4),
            legend.key.size = unit(0.3, 'cm')) + 
      scale_alpha(guide = 'none') 
    mygheatmap
    ggsave(plot = mygheatmap, filename = paste0(outputdir,sep="/",filename), dpi=300, height = 6, width = 8)
    
    
    
