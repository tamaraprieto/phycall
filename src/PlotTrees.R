library(RColorBrewer)
library(data.table)
library(tidyr)
library(ggplot2)
#library(ggpubr)
#library(ggrepel)
library(ggtree)
library(dplyr)
#library(circlize)
#library(ComplexHeatmap)
library(phytools)


# for patient in 4295 445 417 4084
# do
#  scp ny:/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/WGS/${patient}/RESULTS/CellPhy-03/CellPhy.GT10+FO+E.nobulks.noCNVs.raxml.support Downloads/CellPhy.GT10+FO+E.${patient}.nobulks.noCNVs.raxml.support
#  scp ny:/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/WGS/${patient}/RESULTS/isec_somatic-germline/FILTERED/MutationMapping/TreeMutWithZeros.genome.rds  Downloads/TreeMutWithZeros.${patient}.genome.rds
# done

# rename genome files 

#################
# Mutation tree #
#################

patient <- "4295"

# No D9 with ADO
# mytree <- readRDS("Downloads/TreeMut4295_Time_Tree_Null.RDS")


# Branch lengths in number of mutations

# raxml tree
mytree<- read.tree(paste0("Downloads/CellPhy.GT10+FO+E.",patient,".nobulks.noCNVs.raxml.support",sep=""))


# mutation mapping refined tree
mytree <- readRDS(paste0("Downloads/TreeMutWithZeros.",patient,".genome.16thMay23.rds"))
# mytree <- readRDS(paste0("Downloads/TreeMutWithZeros.",patient,".March28th23.genome.rds"))
mytree <- mytree@phylo
mytree$tip.label <- gsub(".*_","",mytree$tip.label)
mytree <- ape::drop.tip(mytree, tip = c("zeros"))


outgroup<- "A11"
#outgroup<- "A1"
mytree <- ape::root(mytree, outgroup)


# Validation of exome using healthy cells
groupInfo <- fread("https://docs.google.com/spreadsheets/d/e/2PACX-1vR91fISOmmJc4WDAsnx9LOXd80YirilNxPSoFBktYrMTWh8fjQpczb3tFH3QXrY30dMeSVD9_mlP4YB/pub?gid=0&single=true&output=csv") %>%
  dplyr::select(Patient,Sample_tag,`Other info`) %>%
  dplyr::filter(`Other info`=="Healthy" | `Other info`=="Almost healthy") %>%
  dplyr::filter(Patient==patient) %>%
  dplyr::rename(Healthy=Sample_tag) %>%
  dplyr::select(Healthy) %>%
  as.list()


# Some cell names are not in the plot
unassigned <- list(setdiff(mytree$tip.label,c(unlist(groupInfo))))
groupInfo[length(groupInfo)+1] <- unassigned
names(groupInfo) <- c(names(groupInfo)[1:(length(groupInfo)-1)],"unassigned")
finaltree <- groupOTU(mytree, groupInfo)


# Colors for groupInfo for clusters
mycols <- brewer.pal(length(groupInfo), "Dark2")
#mycols_mod <- mycols[c(1,2,3,4,5,7,8,9,10)]
#mycols_mod <- c(mycols_mod,"lightgrey")

# Colors for healthy / tumor (1 and 3)
# Second color works for the bootstrap annotations
mycols_mod <- c("lightskyblue","maroon4","forestgreen")
#   mycols_mod <- c("lightskyblue","maroon4")


treefigure <- ggtree(finaltree) + 
  layout_dendrogram() +  # bottom to top
  geom_rootedge(0.1) +
  geom_tiplab(size=3, angle=90, hjust=1) +
  geom_tippoint(size=2) +
  #geom_tippoint(size=2, aes(color=group)) +
  #geom_nodelab(col="forestgreen",angle=0,size=3, hjust = 0.5, vjust=0, show.legend=FALSE) +
  scale_color_manual(values = mycols_mod) +
  #theme_dendrogram() +
  #theme_tree2() +
  theme(legend.position = "none")
treefigure
ggsave(plot = treefigure, filename = paste("Downloads/4295.Numbermutations.png",sep=""), height = 6, width = 12)


treefigure <- ggtree(finaltree) + 
  layout_dendrogram() +  # bottom to top
  geom_rootedge(0.01) +
  geom_tiplab(size=3, angle=90, hjust=1) +
  geom_tippoint(size=2, aes(color=group)) +
  #geom_nodelab(col="forestgreen",size=3, hjust = 0.5, vjust=1.3, show.legend=FALSE) +
  scale_color_manual(values = mycols_mod) +
  #theme_dendrogram() +
  theme(legend.position = "none")

treefigure
ggsave(plot = treefigure, filename = paste("Downloads/4295.Numbermutations.HealthyTumorColor.png",sep=""), height = 6, width = 12)


#############
# Time tree #
#############

#finaltree
#finaltreeultra <- phytools::force.ultrametric(finaltree)
#finaltreeultra <- readRDS("Downloads/")


patient <- "4084"
patient <- "417"
patient <- "445"
patient <- "4295"
# Load the tree
if (patient=="4295"){
  outgroup="A11"
} else if (patient=="417"){
  outgroup <- "B10"
} else if (patient=="4084"){
  outgroup="F1"
} else if (patient=="445") {
  outgroup="A3"
} else {
  outgroup=""
}

# No D9 with ADO
# mytree <- readRDS("Downloads/TreeMut4295_Time_Tree_Null.RDS")



if (patient=="4295") {
  # There is something wrong with Sheng tree, we will have to take a look
  
  #mytree <- readRDS(paste0("Downloads/TreeMut",patient,"_Time_Tree_Null_WithouthADO.RDS",sep=""))
  mytree <- readRDS("/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/timetrees/TreeMut4295_Time_Tree_Null.RDS")
  mytree <- mytree$fit$poisson_tree$nullmodel$ultratree
  #mytree <- mytree$fit$poisson_tree$nullmodel$altmodel$ultratree
} else{
  mytree <- readRDS(paste0("Downloads/TreeMut",patient,"_Time_Tree_Null.RDS",sep=""))
}


mytree$tip.label <- gsub(".*_","",mytree$tip.label)
mytree <- ape::drop.tip(mytree, tip = c("zeros"))


mytree <- ape::root(mytree, outgroup)


# Validation of exome using healthy cells
groupInfo <- data.table::fread("https://docs.google.com/spreadsheets/d/e/2PACX-1vR91fISOmmJc4WDAsnx9LOXd80YirilNxPSoFBktYrMTWh8fjQpczb3tFH3QXrY30dMeSVD9_mlP4YB/pub?gid=0&single=true&output=csv") %>%
  dplyr::select(Patient,Sample_tag,`Other info`) %>%
  dplyr::filter(`Other info`=="Healthy" | `Other info`=="Almost healthy") %>%
  dplyr::filter(Patient==patient) %>%
  dplyr::rename(Healthy=Sample_tag) %>%
  dplyr::select(Healthy) %>%
  as.list()


# Some cell names are not in the plot
unassigned <- list(setdiff(mytree$tip.label,c(unlist(groupInfo))))
groupInfo[length(groupInfo)+1] <- unassigned
names(groupInfo) <- c(names(groupInfo)[1:(length(groupInfo)-1)],"unassigned")
finaltree <- groupOTU(mytree, groupInfo)


# Colors for groupInfo for clusters
mycols <- brewer.pal(length(groupInfo), "Dark2")
#mycols_mod <- mycols[c(1,2,3,4,5,7,8,9,10)]
#mycols_mod <- c(mycols_mod,"lightgrey")

# Colors for healthy / tumor (1 and 3)
# Second color works for the bootstrap annotations
mycols_mod <- c("lightskyblue","maroon4","forestgreen")
#   mycols_mod <- c("lightskyblue","maroon4")






treefigure <- ggtree(finaltree) + 
  #layout_dendrogram() +  # bottom to top
  geom_rootedge(0.01) +
  #geom_tiplab(size=3, angle=90, hjust=1) +
  geom_tippoint(size=2, aes(color=group)) +
  #geom_nodelab(col="forestgreen",size=3, hjust = 0.5, vjust=1.3, show.legend=FALSE) +
  scale_color_manual(values = mycols_mod) +
  #theme_dendrogram() +
  theme(legend.position = "none") +
  theme_tree2(legend.position = "none")
treefigure
ggsave(plot = treefigure, filename = paste("Downloads/",patient,".TimeTree.png",sep=""), height = 2, width = 10)


# Phylogenetic independent contrasts



#######################
# Demographic history #
#######################


# Population size analysis
library(phylodyn) 

###########
# Invitro #
###########

    patient <- "Invitro"
    # old, no A1 and A2
    #mytree <- readRDS(paste0("Downloads/TreeMutWithZeros_Time_Tree.RDS",sep=""))
    mytree <- readRDS(paste0("Downloads/TreeMutWithZerosInvitro_Time_TreewithA1A2.RDS",sep=""))
    mytree$tip.label <- gsub(paste0(patient,"_"),"",mytree$tip.label)
    # Remove one of the single cell expansions, only one subtree with a single-cell origin
    #mytree <- ape::drop.tip(phy = mytree, tip = c("zeros","A1","A2","D1","D2","H1","H2","E1","E2"))
    #mytree <- ape::drop.tip(phy = mytree, tip = c("A1","A2","zeros"))
    mytree <- ape::drop.tip(phy = mytree, tip = c("zeros"))
    
    
    
    tree_plot<- ggtree(mytree) +
      geom_tiplab(size=2) +
      scale_color_manual(values = c("lightblue")) +
      theme_tree2()
    tree_plot
    
    phy_ultra.p <- mytree
    phylodyn = BNPR(data = phy_ultra.p, lengthout = 100) # without sampling model
    myplot_bnpr <- phylodyn$summary %>%
      dplyr::mutate(newtime=-time*365) %>%
    ggplot(aes(x=newtime,y=mean)) +
    #ggplot(aes(x=newtime,y=mean)) +
      #geom_point(data = data.frame(coal_times=-phylodyn$samp_times*365 %>% c()) %>% dplyr::mutate(ypos=1) , aes(coal_times, ypos)) +
      geom_line(col="lightblue", size=1, linetype=2)  +
      geom_ribbon(aes(ymin=quant0.025,ymax=quant0.975),alpha=0.2, fill="grey") +
      geom_point(data = data.frame(coal_times=-phylodyn$samp_times*365%>% c(),
                                   freq=as.factor(phylodyn$n_sampled%>% c())) %>% dplyr::mutate(ypos=0.00001), 
                 aes(coal_times, ypos,col=freq), size=2) +
      scale_y_continuous(trans='log2',breaks = c(100,1e4,1e5,5e5,1e6,1e8,10e9,10e12), labels = c("100","10,000","50,000","100,000","1,000,000","100,000,000","1,000,000,000","1,000,000,000,000")) +
      #scale_x_continuous(limits = c(-58,0), breaks = c(-123,-91,-57,-24,0), labels = c("0","32\nsingle-cell\nseeding","66","99","123")) + # labels = c("123","91","57","24","0")
      scale_x_continuous(breaks = c(-123,-91,-57,-24,0), labels = c("0","32\nsingle-cell\nseeding","66","99","123")) + # labels = c("123","91","57","24","0") +
      theme_bw() +
      scale_color_brewer(palette = "Set1") +
      #scale_color_distiller(palette = 4,direction = 1)
      theme(axis.text.x = element_text(angle=45, hjust = 1)) +
      labs(y = "Effective cell population size\n(in log scale)",x="Cell culture time\n(days)",col="Cells sampled")# +
      #coord_cartesian(xlim = c(-58, 0),ylim=c(1,10e7))


ggpubr::ggarrange(tree_plot + theme(plot.margin = unit(c(0,34,0,20), "mm")), 
                  myplot_bnpr,
                  ncol = 1, heights = c(0.2,1))

plot_BNPR(phylodyn, 
          main=paste0(patient),
          col = rgb(0.829, 0.680, 0.306))


#################################
# Simulate sequences under the coalescent with sampling #
#################################



# We can get the decay FACTOR from: 
doubling_time=20
growth_rate <- log(x=2,base = exp(1))/doubling_time
# growth rate or exponential decay FACTOR
growth_rate
# We can get the exponential decay RATE from the FACTOR
# decay_factor = 1 - decay_rate
decay_rate <- 1 - growth_rate
decay_rate

new_times <- as.numeric(names(table(phylodyn$samp_times)))
new_times2 <- rev(-(new_times - max(new_times)))
descendants_f <- coalsim(traj = exp_traj,
                samp_times = new_times2[4],
                n_sampled = 2,
                rate = decay_rate)

descendants_g <- coalsim(traj = exp_traj,
                     samp_times = new_times2[4],
                     n_sampled = 2,
                     rate = decay_rate)


ancestors_fg <- coalsim(traj = exp_traj,
                     samp_times = new_times2[3],
                     n_sampled = 2,
                     rate = decay_rate)


ancestors_fg_tree <- generate_newick(ancestors_fg)$newick
descendants_g_tree <- generate_newick(descendants_g)$newick
descendants_f_tree <- generate_newick(descendants_f)$newick



ggtree(descendants_f_tree) +
  theme_tree2()
ggtree(descendants_g_tree) +
  theme_tree2()


ggtree(ancestors_fg_tree) +
  theme_tree2()



#################
# Patients BALL #
#################

patient <- "4295"
patient <- "417"
patient <- "445"
patient <- "4084"


RunPhylodyn <- function(patient){
  

        mytree <- readRDS(paste0("Downloads/TreeMut",patient,"_Time_Tree_Null.RDS",sep=""))
        length(mytree$tip.label)
        # remove healthy cells
        healthy_cells <- fread("https://docs.google.com/spreadsheets/d/e/2PACX-1vR91fISOmmJc4WDAsnx9LOXd80YirilNxPSoFBktYrMTWh8fjQpczb3tFH3QXrY30dMeSVD9_mlP4YB/pub?gid=0&single=true&output=csv") %>%
          dplyr::select(Patient,Sample_tag,`Other info`) %>%
          dplyr::filter(`Other info`=="Healthy" | `Other info`=="Almost healthy") %>%
          dplyr::filter(Patient==patient) %>%
          dplyr::rename(Healthy=Sample_tag) %>%
          dplyr::mutate(Healthy=paste0(patient,"_",Healthy)) %>%
          dplyr::select(Healthy) %>%
          unlist()  %>% as.character()
        mytree <- ape::drop.tip(phy = mytree, tip = c(healthy_cells,"zeros"))
        length(mytree$tip.label)
        
        phy_ultra.p <- mytree


        ##################
        # mcmc (slower)  #
        ##################
        # mcmc_sampling()
        
        
        ############################
        # INLA approximation: BNPR #
        ############################
        
        #!!!!!!!!!!!!!!!!!!!
        # What should be the lengthout? #
        #################################
        # lengthout default is 100
        # In van egeren they used the time of the tumor. 28 out of 34 years old
        # if I use a small interval then the confidence intervals increase a lot
        # phylodyn = BNPR(data = phy_ultra.p, lengthout = max(mytree$edge.length)) # without sampling model
        phylodyn = BNPR(data = phy_ultra.p, lengthout = 100) # without sampling model
        
        
        
        # In cases where the frequency of sampling is related to 
        # effective population size, including a sampling time model
        # provides additional accuracy and precision
        # In our case, we are not sampling more when the
        # population is bigger so it is not useful
        # phylodyn_sampl = BNPR_PS(data = phy_ultra.p, lengthout = 100) # with sampling model. 
        
        
        phylodyn_selected <- phylodyn
        #phylodyn_selected <- phylodyn_sampl
        #plot_BNPR(phylodyn_selected, 
        #          main=paste0(patient),
        #          col = rgb(0.829, 0.680, 0.306))
        
        
        results <- phylodyn_selected$summary %>%
          as.data.frame() %>%
          dplyr::mutate(Patient=patient)
        
    return(results)    
}


mytable <- RunPhylodyn(patient = "4295")
ggplot(mytable,aes(x=-time,y=mean)) +
  geom_line(col="plum3", size=1)  +
  geom_ribbon(aes(ymin=quant0.025,ymax=quant0.975),alpha=0.2, fill="plum3") +
  scale_y_continuous(trans='log2', breaks = c(1,100,1e4,1e6,5e7,1e9,1e12,40e12), labels = c("1","100","10,000","1 million","50 million","1 billion","1 trillion","40 trillion")) +
  scale_x_continuous(limits = c(-16,0),breaks = c(-16,-12,-10,-8,-6,-4,-2,0), labels = c("2004\n(birth)","2008","2010","2012","2014","2016","2018","2020\n(diagnosis)")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  labs(y = "Effective tumor cell\npopulation size",x="Time in years", title = paste0("Patient A"))
ggsave(filename = paste("Documents/Figs.BALL.PTA/Phylodynamics.4295.pdf",sep=""), dpi = 300, height = 3, width = 6)


mytable_all <- do.call("rbind.data.frame",lapply(X = c("417","445"),RunPhylodyn))
ggplot(mytable_all %>% dplyr::mutate(Patient=factor(as.factor(Patient),levels=c("417","445"))),
       aes(x=-time,y=mean,col=Patient, fill=Patient)) +
  geom_line(size=2)  +
  geom_ribbon(aes(ymin=quant0.025,ymax=quant0.975),alpha=0.2) +
  scale_y_continuous(trans='log2', breaks = c(1,100,1e4,1e6,1e9,1e12,40e12), labels = c("1","100","10,000","1 million","1 billion","1 trillion","40 trillion")) +
  scale_x_continuous(breaks = c(-16,-14,-12,-10,-8,-6,-4,-2,0), labels = c("2004","2006","2008","2010","2012","2014","2016","2018","2020\n(diagnosis)")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust = 1), 
        axis.title.y  = element_text(hjust = 0), legend.position = "right") +
  labs(y = "Effective tumor cell population size",x="Time in years") +
  scale_fill_manual(values=c("firebrick3","mediumaquamarine","olivedrab4")) +
  scale_color_manual(values=c("firebrick3","mediumaquamarine","olivedrab4"))
ggsave(filename = paste("Documents/Figs.BALL.PTA/Phylodynamics.all.pdf",sep=""), dpi = 300, height = 5, width = 6)


mytable_all

ggplot(mytable,aes(x=-time,y=log(mean))) +
  geom_line()  +
  geom_ribbon(aes(ymin=log(quant0.025),ymax=log(quant0.975)),alpha=0.2) +
  labs(y = "Effective population size\n(log scaled)",x="Time\n(in years from present to past)")



# ape::skyline implements the generalized 
# skyline method for isochronous genealogies
test <- ape::skyline(phy_ultra.p)
plot(test, show.years=TRUE, subst.rate=0.63, present.year = 2022, col=c(grey(.8),1))


#######################
# Patients esophageal #
#######################


tree_branchinnummutations <- readRDS(paste0("Downloads/phylogeny_annotated.RDS"))
tree_branchinnummutations <- treeio::drop.tip("germline")
#tree_branchinnummutations <- treeio::drop.tip("zeros")

# update info about exome mutations
#treemut_RDS_exome <- tree_branchinnummutations %>% as_tibble() %>%
  #dplyr::mutate(ogaini=ifelse(grepl("ogacor",mutList2),"previously\n described\ndriver","not previously\ndescribed\ndriver")) %>%
  #dplyr::mutate(mutListexome=gsub("-ogacor([A-Z]|[a-z][0-1])+","",gsub("_[0-9]","",gsub(",","\n",gsub(";","",gsub(";chr([0-9]|X|Y|:)+",";",mutList2)))))) %>%
  #dplyr::mutate(mutListexome2=ifelse(grepl("not",ogaini),NA,mutListexome)) %>%
  #treeio::as.treedata()
#tree_branchinnummutations<- treemut_RDS_exome

mycols_mod <- c("maroon","darkcyan","lightgrey")

ggtree(tree_branchinnummutations, col="darkgrey") +
  geom_tippoint(aes(color=celltype), size=2) +
  geom_nodelab(aes(label=node)) +
  #geom_tiplab() +
  scale_color_manual(values = mycols_mod) 


# SUBSAMPLE A SUBTREE WITHOUT IMMUNE
# function of treeio to a subtree 

# obtain descendant nodes (do not get tips)
# this is good because I don't want the terminal node mutations
immune_nodes <- phytools::getDescendants(tree_branchinnummutations@phylo, 103)
otherepi_nodes <- setdiff(tree_branchinnummutations@data$node,immune_nodes)
immune_cells <- tree_branchinnummutations %>% as_tibble() %>% dplyr::select(label,celltype) %>% dplyr::filter(celltype=="CD45+") %>% dplyr::select(label) %>% unlist()
# get subclones and draw the subtrees
immune_clone <- treeio::tree_subset(tree = tree_branchinnummutations,
                                           node=103)  # it doesn't matter which node as long as there are no immune to mess up, keep increasing
eso.confidence.tree <- treeio::drop.tip(tree_branchinnummutations, tip = immune_cells)


mytree <- eso.confidence.tree@phylo
length(mytree$tip.label)
phy_ultra.p <- mytree


#!!!!!!!!!!!!!!!!!!!
# What should be the lengthout? #
#################################
# lengthout default is 100
# In van egeren they used the time of the tumor. 28 out of 34 years old
# if I use a small interval then the confidence intervals increase a lot
# phylodyn = BNPR(data = phy_ultra.p, lengthout = max(mytree$edge.length)) # without sampling model
phylodyn = phylodyn::BNPR(data = phy_ultra.p, lengthout = 100) # without sampling model

# In cases where the frequency of sampling is related to 
# effective population size, including a sampling time model
# provides additional accuracy and precision
# In our case, we are not sampling more when the
# population is bigger so it is not useful
# phylodyn_sampl = BNPR_PS(data = phy_ultra.p, lengthout = 100) # with sampling model. 


phylodyn_selected <- phylodyn
#phylodyn_selected <- phylodyn_sampl
plot_BNPR(phylodyn_selected, 
          main=paste0("cu02"),
          col = rgb(0.829, 0.680, 0.306))


mytable <- phylodyn_selected$summary %>%
  as.data.frame() %>%
  dplyr::mutate(Patient="cu02")

ggplot(mytable,aes(x=-time,y=mean)) +
  geom_line(col="plum3", size=1)  +
  geom_ribbon(aes(ymin=quant0.025,ymax=quant0.975),alpha=0.2, fill="plum3") +
  #scale_y_continuous(trans='log2', breaks = c(1,100,1e4,1e6,5e7,1e9,1e12,40e12), labels = c("1","100","10,000","1 million","50 million","1 billion","1 trillion","40 trillion")) +
  #scale_x_continuous(limits = c(-16,0),breaks = c(-16,-12,-10,-8,-6,-4,-2,0), labels = c("2004\n(birth)","2008","2010","2012","2014","2016","2018","2020\n(diagnosis)")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  labs(y = "Effective tumor cell\npopulation size",x="Time in years", title = paste0("Patient A"))
ggsave(filename = paste("Documents/Figs.BALL.PTA/Phylodynamics.4295.pdf",sep=""), dpi = 300, height = 3, width = 6)


mytable_all <- do.call("rbind.data.frame",lapply(X = c("417","445"),RunPhylodyn))
ggplot(mytable_all %>% dplyr::mutate(Patient=factor(as.factor(Patient),levels=c("417","445"))),
       aes(x=-time,y=mean,col=Patient, fill=Patient)) +
  geom_line(size=2)  +
  geom_ribbon(aes(ymin=quant0.025,ymax=quant0.975),alpha=0.2) +
  scale_y_continuous(trans='log2', breaks = c(1,100,1e4,1e6,1e9,1e12,40e12), labels = c("1","100","10,000","1 million","1 billion","1 trillion","40 trillion")) +
  scale_x_continuous(breaks = c(-16,-14,-12,-10,-8,-6,-4,-2,0), labels = c("2004","2006","2008","2010","2012","2014","2016","2018","2020\n(diagnosis)")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust = 1), 
        axis.title.y  = element_text(hjust = 0), legend.position = "right") +
  labs(y = "Effective tumor cell population size",x="Time in years") +
  scale_fill_manual(values=c("firebrick3","mediumaquamarine","olivedrab4")) +
  scale_color_manual(values=c("firebrick3","mediumaquamarine","olivedrab4"))
ggsave(filename = paste("Documents/Figs.eso.RSO/Phylodynamics.cu02.pdf",sep=""), dpi = 300, height = 5, width = 6)




###################################
# my plots to understand dynamics #
###################################


y <- c(2,4,8,16,32,64,128)
y <- c(2,4,8,16,32,32,32)
x <- c(1,2,3,4,5,6,7) 
ggplot(data = cbind(x,y) %>% as.data.frame() ,aes(x,y)) +
  geom_line()
ggplot(data = cbind(x,y) %>% as.data.frame() ,aes(x,log(y))) +
  geom_line()




