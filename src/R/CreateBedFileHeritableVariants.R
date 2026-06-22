#!/usr/bin/env Rscript
library(data.table)
library(ape)
library(dplyr)
library(tidyr)
library(PATH)

# get the arguments
args = commandArgs(trailingOnly=TRUE)
patient <- args[1]

mydir=paste0("/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/WGS/",patient)


ad <- fread(paste0(mydir,"/RESULTS//isec_somatic-germline/FILTERED/MutationMapping/Mutect.AD.bed"))
#head(ad)
#dim(ad)
#mytreemapped <- readRDS("/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/Invitro/Invitro_UG/RESULTS/isec_somatic-germline/FILTERED/MutationMapping/TreeMut.genome.rds")
#colnames(ad) <- gsub("-","",gsub("-L1-.*","",gsub("DY-","",colnames(ad))))

# Load the tree
  tree_raxml <- read.tree(paste(mydir,"/RESULTS/CellPhy-03/CellPhy.GT10+FO+E.nobulks.noCNVs.raxml.bestTree",sep=""))


Winv <- PATH:::inv.tree.dist(tree_raxml)
vafs <- ad[,-1] %>% as.matrix()
vafs_mat <- matrix(as.numeric(vafs),    # Convert to numeric matrix
                  ncol = ncol(vafs))
colnames(vafs_mat) <- colnames(vafs)
vafs_mat <- vafs_mat %>%
                       as.data.frame() %>%
                       dplyr::select(tree_raxml$tip.label)


# xcor is too slow with more than 10000 variable
# for that reason I am only going to use the other function below

parM <- function(d,w, break.length=100) {
    
  
    y <- c(seq(1,ncol(d), break.length), (ncol(d)+1) )
    y2 <- y[-1] - 1
    y1 <- y[-length(y)]
    
    mfunc <- function(b, e, d1=d, w1=w) {
        M <- xcor(d1[,b:e ], w1)
        x1 <- diag(M$Z)
        x2 <- diag(M$Morans.I)
        out <- data.frame("Z"=x1, "I"=x2)
        return(out)
    }
    
    out <- do.call(rbind,
                   parallel::mcmapply(function(b1,e1) mfunc(b=b1,e=e1, d, w), b1=y1, e1=y2,
                            mc.cores = parallel::detectCores(), SIMPLIFY = F))
    return(out)
}

morans <- parM(d = vafs_mat %>%
                       as.matrix() %>% t(),
     w = Winv, break.length = 100)


#ggplot(morans) +
#  geom_density(aes(x=Z)) +
#  labs(x="Autocorrelation (Z-score)") +
#  theme_classic()

torescue <- which(morans$Z>=2)
toremove <- which(morans$Z<2)
paste0("Number of mutations to rescue: ", length(torescue),". Number of mutations to ignore: ",length(toremove) ,sep="")

# Create a bed file to intersect the vcf with
ad$`CHROM:POS`[which(morans$Z>=2)] %>%
  as.data.frame() %>%
  dplyr::rename(chrompos=".") %>%
  tidyr::separate(chrompos, into=c("chrom","end"), sep=":") %>%
  dplyr::mutate(start=as.numeric(end)-1) %>%
  dplyr::select(chrom, start,end) %>%
  write.table(file = paste0(mydir,"/RESULTS/isec_somatic-germline/FILTERED/MutationMapping/FilteredToRescue.bed",sep=""), quote = F, sep = "\t",row.names = F, col.names = F)

