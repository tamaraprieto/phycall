source("/gpfs/commons/home/srajagopalan/PTA/treemut_PTA_package/code.R")
source("/gpfs/commons/home/srajagopalan/PTA/treemut_PTA_package/cna.R")
# library("kableExtra")
# library("ggtree")
# library(ggtree)
# library(data.table)
# library(ape)
# library(spam)
# library(INLA)
# library(phylodyn)
# library(phangorn)
# library(ComplexHeatmap)
# library(phytools)
# #library(treeio)
# library(ggplot2)
# library(tidyverse)
# library(ggnewscale)
# library(lsa)
# library(seqinr)
# library(binom)
# library(castor)
# library(corrplot)
# library(ggpubr)

args = commandArgs(trailingOnly=TRUE)

input=args[1]
output=args[2]
print(input)
print(output)

### First, we need to read in the treemut object, and add on things to make it look like the actual PD object for wraptreefit ####

inital_treemut_RDS_object <- readRDS(input)

# removing empty row/node associated with zeros 
new_data_table <- inital_treemut_RDS_object@data
# tama: added this for my data
new_data_table <-  new_data_table %>%
  dplyr::select(-mutList, -mutList2)
empty_row <- which(is.na(new_data_table$original.branch.length))
new_data_table <- new_data_table[-c(empty_row),] # remove row with NA values
new_data_table$node[which(new_data_table$node >= empty_row)] <- new_data_table$node[which(new_data_table$node > empty_row)] - 1 # correcting node id's after removing zeros
inital_treemut_RDS_object@data <- new_data_table


#First thing is making an agedf file - this file 

#Actual agedf
#PD$pdx$agedf
PD=readRDS("/gpfs/commons/home/srajagopalan/PTA/treemut_PTA_package/PD6629.RDS")


PD_test <- list()
agedf <- as.data.frame.matrix(matrix(NA, nrow = length(inital_treemut_RDS_object@phylo$tip.label), ncol = 8))
colnames(agedf) <- colnames(PD$pdx$agedf)
agedf$tip.label <- inital_treemut_RDS_object@phylo$tip.label
agedf$age_at_sample_pcy <- rnorm(length(agedf$tip.label), mean = 16, sd = 0.001)  # patient 4295: age at sampling = 16
agedf$telo_mean_length <- rnorm(length(agedf$tip.label), mean = 1416, sd = 100)
agedf$age_at_sample_pcy[86] <- 0.000001 
agedf$telo_mean_length[86] <- NA
agedf$patient <- '4295'
agedf$per.sample.sensitivity.hybrid <- 1
agedf$per.sample.sensitivity.reg <- 1
agedf$driver <- 'JAK2'
agedf$driver3 <- 'JAK2'

# #making tree_ml object
# tree_ml <- list()
# tree_ml$edge <- inital_treemut_RDS_object$tree$edge
# tree_ml$edge.length <- inital_treemut_RDS_object$tree$edge.length
# tree_ml$Nnode <- inital_treemut_RDS_object$tree$Nnode
# tree_ml$tip.label <- inital_treemut_RDS_object$tree$tip.label
# tree_ml$label <- inital_treemut_RDS_object$df$df$profile
# tree_ml$el.snv.local.filtered <- inital_treemut_RDS_object$df$df$mut_count
# tree_ml$per.branch.sensitivity.hybrid.multi <- rep(1, length(tree_ml$edge.length))

tree_ml <- inital_treemut_RDS_object@phylo
tree_ml$label <- inital_treemut_RDS_object@data$profile.x
tree_ml$el.snv.local.filtered <- inital_treemut_RDS_object@data$edge_length
tree_ml$per.branch.sensitivity.hybrid.multi <- rep(1, length(tree_ml$edge.length))

#making nodes object - let's start with one small node 
nodes <- as.data.frame.matrix(matrix(NA, nrow = 1, ncol = length(PD$nodes)))
colnames(nodes) <- colnames(PD$nodes)
nodes$node <- 136
nodes$driver <- 'SULF2, CDHR1'
nodes$status <- 1
nodes$driver2 <- 'SULF2'
nodes$driver3 <- 'SULF2:CDHR1'
nodes$child_count <- node.depth(inital_treemut_RDS_object@phylo)[136]

#now we need to create the object itself
PD_test$pdx <- list()
PD_test$pdx$agedf <- agedf
PD_test$pdx$tree_ml <- tree_ml
PD_test$nodes <- nodes
PD_test$localx.correction2 <- 1.01476

PD_test$patient <- '4295'

#so now we have our object ready
NITER=1000
#The following populates a list PD$fit$<treemodel>$<altmodel|nullmodel>.  The function wraptreefit wraps a call to the rtreefit package.
treemodel = "poisson_tree"

PD_test$pdx$tree_ml <- PD_test$pdx$tree_ml %>% purrr::list_modify("node.label" = NULL)

PD_test=wraptreefit(PD_test,niter=NITER,b.fit.null = TRUE, method=treemodel)

saveRDS(PD_test, output)
