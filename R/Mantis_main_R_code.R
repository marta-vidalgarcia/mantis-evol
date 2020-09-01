# MANTIS-EVOL 
rm(list=ls())

#### 0. Libraries and data loading ####
library(ape)
library(phytools)
library(geiger)
library(geomorph)
library(phangorn)
library(phylobase)
library(picante)
library(caper)


phyloTime <- read.nexus("BeastTree")
phyloTimeLadderized <- (ladderize(phyloTime))
plot(phyloTimeLadderized, cex=0.5, no.margin = T)

# write.tree(phyloTimeLadderized, "phyloTimeLadderized.nwk")
# phyloTimeLadderized <- read.tree("phyloTimeLadderized.nwk")

tree <-phyloTimeLadderized

tree_pre <- tree

data_no_order <- read.csv("Final_characters_CSV.csv", header=T, row.names = 1)

head(data_no_order)

row.names(data_no_order)

tree$tip.label

tree_pre<-ladderize(tree_pre)
Overlap <- treedata(tree_pre, data_no_order)

data_pre_order <- as.data.frame(Overlap$data)


data <- data_pre_order[match(Overlap$phy$tip.label,rownames(data_pre_order)),] 
tree <- Overlap$phy


#### 1. Evolutionary distinctiveness and testing bias ####
evol_distinct <- evol.distinct(tree, type = "equal.splits", scale = FALSE, use.branch.lengths = TRUE)
evol_distinct_reduced <- evol.distinct(tree_reduced, type = "equal.splits", scale = FALSE, use.branch.lengths = TRUE)
evol_d <- evol_distinct$w
names(evol_d) <- as.character(evol_distinct$Species)
evol_d_r <- evol_distinct_reduced$w
names(evol_d_r) <- evol_distinct_reduced$Species

# Authors (are they biased in whether they preferentially work with species with deimatic behaviour or not?)
fitPagel(tree, data$Presence_display, data$Edmunds)
fitPagel(tree, data$primary_defense_special_resemblance_masquerade, data$Edmunds)
fitPagel(tree, data$primary_defense_general_resemblance_crypsis, data$Edmunds) 
fitPagel(tree, data$World, data$Edmunds) 


# Authors (are they biased in whether they report species that are evolutionarily distinct?)
authors <- data$behav_author
names(authors) <- row.names(data)
phylANOVA(tree, authors, evol_d, nsim=1000, posthoc=TRUE, p.adj="holm")

# And only Edmunds?
edmunds <- data$Edmunds
names(edmunds) <- row.names(data)
phylANOVA(tree, edmunds, evol_d, nsim=1000, posthoc=TRUE, p.adj="holm")
dim(subset(data, Edmunds==1))[1] # And he described 31 of the behaviours
dim(subset(data, Edmunds==1))[1]/length(edmunds)*100 # 53.5% of all papers used in this study

data_ed <- subset(data, Edmunds==1)
tip_ed <- row.names(subset(data, Edmunds==0))
tree_ed <- drop.tip(tree, tip_ed)
plot(tree_ed)
length(tree_ed$tip.label) == dim(subset(data, Edmunds==1))[1] 

tip <- row.names(subset(data, Edmunds==1))
pos <- vector("numeric", length = length(tip))
for (i in 1:length(tip)){
  pos[i] <- which(names(evol_d) == tip[i])
}
evol_d_ed <- evol_d[pos]

# Display complexity
display_complexity <- data$display_complexity
names(display_complexity) <- row.names(data)
phylANOVA(tree, display_complexity, evol_d, nsim=1000, posthoc=TRUE, p.adj="holm")

display_complexity_ed <- data_ed$display_complexity
names(display_complexity_ed) <- row.names(data_ed)
phylANOVA(tree_ed, display_complexity_ed, evol_d_ed, nsim=1000, posthoc=TRUE, p.adj="holm")

# Presence display
Presence_display <- data$Presence_display
names(Presence_display) <- row.names(data)
phylANOVA(tree, Presence_display, evol_d, nsim=1000, posthoc=TRUE, p.adj="holm")

Presence_display_ed <- data_ed$Presence_display
names(Presence_display_ed) <- row.names(data_ed)
phylANOVA(tree_ed, Presence_display_ed, evol_d_ed, nsim=1000, posthoc=TRUE, p.adj="holm")


#### 2. ANCESTRAL STATE RECONSTRUCTION ####

display_complexity <- data$display_complexity
names (display_complexity) <-row.names(data)

fitER<-ace(display_complexity,tree,model="ER",type="discrete")
fitER$rates
fitER$loglik

fitSYM<-ace(display_complexity,tree,model="SYM",type="discrete")
fitSYM$loglik


fitARD<-ace(display_complexity,tree,model="ARD",type="discrete")
fitARD$loglik

1-pchisq(2*abs(fitARD$loglik - fitER$loglik), 1)



fitER <- rerootingMethod(tree, display_complexity, model = "ER")
print(fitER)

plotTree(tree,fsize=0.8,ftype="i")

cols<-setNames(palette()[1:length(unique(display_complexity))],sort(unique(display_complexity)))
tiplabels(pie=to.matrix(display_complexity,sort(unique(display_complexity))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(tree)),fsize=0.8)


plotTree(tree,fsize=0.8,ftype="i")
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
           piecol = c("azure", "aquamarine4", "aquamarine3", "darkcyan", "darkblue"), cex = 0.6)
tiplabels(pie = to.matrix(display_complexity, sort(unique(display_complexity))), piecol = c("azure", "aquamarine4", "aquamarine3", "darkcyan", "darkblue"), 
          cex = 0.3)
add.scale.bar()



Presence_display <- data$Presence_display
names (Presence_display) <-row.names(data)

fitER<-ace(Presence_display,tree,model="ER",type="discrete")
fitER

plotTree(tree,fsize=0.8,ftype="i")

cols<-setNames(palette()[1:length(unique(Presence_display))],sort(unique(Presence_display)))
tiplabels(pie=to.matrix(Presence_display,sort(unique(Presence_display))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(tree)),fsize=0.8)

nodelabels(node=1:tree$Nnode+Ntip(tree),
           pie=fitER$lik.anc,piecol=cols,cex=0.5)



#### 3. Statistical analyses ####

#### 3.1. fitPagel ####
# Function to test for correlated evolution of binary traits

fitPagel(tree, data$primary_defense_special_resemblance_masquerade, data$Presence_display)
fitPagel(tree, data$primary_defense_special_resemblance_masquerade, data$secondary_defence_present)
fitPagel(tree, data$secondary_defence_present, data$primary_defense_general_resemblance_crypsis)
fitPagel(tree, data$primary_defense_special_resemblance_masquerade, data$Presence_display)
fitPagel(tree, data$Presence_display, data$World)
fitPagel(tree, data$primary_defense_special_resemblance_masquerade, data$World)


# Binary traits reduced dataset: wings_display	arms_display	wings_colours	arms_colours	abdomen_colours	wings_eyespots	arms_eyespots	
# display_includes_sound_abd_or_wings	mouth_display


data_no_order_reduced <- subset(data, Presence_display == 1)


Overlap_reduced <- treedata(tree_pre, data_no_order_reduced)

data_pre_order_reduced <- as.data.frame(Overlap$data)
data_reduced <- data_pre_order_reduced[match(Overlap_reduced$phy$tip.label,rownames(data_pre_order_reduced)),] 
tree_reduced <- Overlap_reduced$phy

fitPagel(tree_reduced, data_reduced$primary_defense_special_resemblance_masquerade, data_reduced$wings_display)
fitPagel(tree_reduced, data_reduced$primary_defense_special_resemblance_masquerade, data_reduced$arms_display)
fitPagel(tree_reduced, data_reduced$primary_defense_special_resemblance_masquerade, data_reduced$wings_colours)
fitPagel(tree_reduced, data_reduced$primary_defense_special_resemblance_masquerade, data_reduced$arms_colours)
fitPagel(tree_reduced, data_reduced$primary_defense_special_resemblance_masquerade, data_reduced$abdomen_colours)
fitPagel(tree_reduced, data_reduced$primary_defense_special_resemblance_masquerade, data_reduced$wings_eyespots)
fitPagel(tree_reduced, data_reduced$primary_defense_special_resemblance_masquerade, data_reduced$arms_eyespots)
fitPagel(tree_reduced, data_reduced$primary_defense_special_resemblance_masquerade, data_reduced$display_includes_sound_abd_or_wings)
fitPagel(tree_reduced, data_reduced$primary_defense_special_resemblance_masquerade, data_reduced$mouth_display)

fitPagel(tree_reduced, data_reduced$Presence_display, data_reduced$wings_display)
fitPagel(tree_reduced, data_reduced$Presence_display, data_reduced$arms_display)
fitPagel(tree_reduced, data_reduced$Presence_display, data_reduced$wings_colours)
fitPagel(tree_reduced, data_reduced$Presence_display, data_reduced$arms_colours)
fitPagel(tree_reduced, data_reduced$Presence_display, data_reduced$abdomen_colours)
fitPagel(tree_reduced, data_reduced$Presence_display, data_reduced$wings_eyespots)
fitPagel(tree_reduced, data_reduced$Presence_display, data_reduced$arms_eyespots)
fitPagel(tree_reduced, data_reduced$Presence_display, data_reduced$display_includes_sound_abd_or_wings)
fitPagel(tree_reduced, data_reduced$Presence_display, data_reduced$mouth_display)

fitPagel(tree_reduced, data_reduced$World, data_reduced$wings_display)
fitPagel(tree_reduced, data_reduced$World, data_reduced$arms_display)
fitPagel(tree_reduced, data_reduced$World, data_reduced$wings_colours)
fitPagel(tree_reduced, data_reduced$World, data_reduced$arms_colours)
fitPagel(tree_reduced, data_reduced$World, data_reduced$abdomen_colours)
fitPagel(tree_reduced, data_reduced$World, data_reduced$wings_eyespots)
fitPagel(tree_reduced, data_reduced$World, data_reduced$arms_eyespots)
fitPagel(tree_reduced, data_reduced$World, data_reduced$display_includes_sound_abd_or_wings)
fitPagel(tree_reduced, data_reduced$World, data_reduced$mouth_display)

fitPagel(tree_reduced, data_reduced$wings_display, data_reduced$arms_display)
fitPagel(tree_reduced, data_reduced$wings_display, data_reduced$wings_colours)
fitPagel(tree_reduced, data_reduced$wings_display, data_reduced$arms_colours)
fitPagel(tree_reduced, data_reduced$wings_display, data_reduced$abdomen_colours)
fitPagel(tree_reduced, data_reduced$wings_display, data_reduced$wings_eyespots)
fitPagel(tree_reduced, data_reduced$wings_display, data_reduced$arms_eyespots)
fitPagel(tree_reduced, data_reduced$wings_display, data_reduced$display_includes_sound_abd_or_wings)
fitPagel(tree_reduced, data_reduced$wings_display, data_reduced$mouth_display)

fitPagel(tree_reduced, data_reduced$arms_display, data_reduced$wings_colours)
fitPagel(tree_reduced, data_reduced$arms_display, data_reduced$arms_colours)
fitPagel(tree_reduced, data_reduced$arms_display, data_reduced$abdomen_colours)
fitPagel(tree_reduced, data_reduced$arms_display, data_reduced$wings_eyespots)
fitPagel(tree_reduced, data_reduced$arms_display, data_reduced$arms_eyespots)
fitPagel(tree_reduced, data_reduced$arms_display, data_reduced$display_includes_sound_abd_or_wings)
fitPagel(tree_reduced, data_reduced$arms_display, data_reduced$mouth_display)

fitPagel(tree_reduced, data_reduced$wings_colours, data_reduced$arms_colours)
fitPagel(tree_reduced, data_reduced$wings_colours, data_reduced$abdomen_colours)
fitPagel(tree_reduced, data_reduced$wings_colours, data_reduced$wings_eyespots)
fitPagel(tree_reduced, data_reduced$wings_colours, data_reduced$arms_eyespots)
fitPagel(tree_reduced, data_reduced$wings_colours, data_reduced$display_includes_sound_abd_or_wings)
fitPagel(tree_reduced, data_reduced$wings_colours, data_reduced$mouth_display)

fitPagel(tree_reduced, data_reduced$arms_colours, data_reduced$abdomen_colours)
fitPagel(tree_reduced, data_reduced$arms_colours, data_reduced$wings_eyespots)
fitPagel(tree_reduced, data_reduced$arms_colours, data_reduced$arms_eyespots)
fitPagel(tree_reduced, data_reduced$arms_colours, data_reduced$display_includes_sound_abd_or_wings)
fitPagel(tree_reduced, data_reduced$arms_colours, data_reduced$mouth_display)

fitPagel(tree_reduced, data_reduced$abdomen_colours, data_reduced$wings_eyespots)
fitPagel(tree_reduced, data_reduced$abdomen_colours, data_reduced$arms_eyespots)
fitPagel(tree_reduced, data_reduced$abdomen_colours, data_reduced$display_includes_sound_abd_or_wings)
fitPagel(tree_reduced, data_reduced$abdomen_colours, data_reduced$mouth_display)

fitPagel(tree_reduced, data_reduced$wings_eyespots, data_reduced$arms_eyespots)
fitPagel(tree_reduced, data_reduced$wings_eyespots, data_reduced$display_includes_sound_abd_or_wings)
fitPagel(tree_reduced, data_reduced$wings_eyespots, data_reduced$mouth_display)

fitPagel(tree_reduced, data_reduced$arms_eyespots, data_reduced$display_includes_sound_abd_or_wings)
fitPagel(tree_reduced, data_reduced$arms_eyespots, data_reduced$mouth_display)

fitPagel(tree_reduced, data_reduced$display_includes_sound_abd_or_wings, data_reduced$mouth_display)


#Evolutionary distinctiveness
evol_distinct <- evol.distinct(tree, type = "equal.splits", scale = FALSE, use.branch.lengths = TRUE)
evol_distinct_reduced <- evol.distinct(tree_reduced, type = "equal.splits", scale = FALSE, use.branch.lengths = TRUE)
evol_d <- evol_distinct$w
names(evol_d) <- as.character(evol_distinct$Species)
evol_d_r <- evol_distinct_reduced$w
names(evol_d_r) <- evol_distinct_reduced$Species

# PhyloANOVAs with evol_d
length(data_reduced$abdomen_colours)
phylANOVA(tree, data$behaviour_species_family, evol_d, nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree, data$Region, evol_d, nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree, data$World, evol_d, nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree, data$primary_defense_special_resemblance_masquerade, evol_d, nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree, data$masquerade_model, evol_d, nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree, data$primary_defense_general_resemblance_crypsis, evol_d, nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree, data$secondary_defence_present, evol_d, nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree, data$display_complexity, evol_d, nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree, data$display_complexity_only_1_4, evol_d, nsim=1000, posthoc=TRUE, p.adj="holm")

phylANOVA(tree_reduced, data_reduced$wings_display, evol_d_r, nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree_reduced, data_reduced$arms_display, evol_d_r, nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree_reduced, data_reduced$wings_colours, evol_d_r, nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree_reduced, data_reduced$arms_colours, evol_d_r, nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree_reduced, data_reduced$abdomen_colours, evol_d_r, nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree_reduced, data_reduced$wings_eyespots, evol_d_r, nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree_reduced, data_reduced$arms_eyespots, evol_d_r, nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree_reduced, data_reduced$display_includes_sound_abd_or_wings, evol_d_r, nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree_reduced, data_reduced$mouth_display, evol_d_r, nsim=1000, posthoc=TRUE, p.adj="holm")


#PhyloANOVAs with dimorphism body size

length(data$DIMORPHISM_bodylength)
length(data_reduced$DIMORPHISM_bodylength)

phylANOVA(tree, data$behaviour_species_family, as.numeric(as.character(data$DIMORPHISM_bodylength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree, data$Region, as.numeric(as.character(data$DIMORPHISM_bodylength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree, data$World, as.numeric(as.character(data$DIMORPHISM_bodylength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree, data$primary_defense_special_resemblance_masquerade, as.numeric(as.character(data$DIMORPHISM_bodylength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree, data$masquerade_model, as.numeric(as.character(data$DIMORPHISM_bodylength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree, data$primary_defense_general_resemblance_crypsis, as.numeric(as.character(data$DIMORPHISM_bodylength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree, data$secondary_defence_present, as.numeric(as.character(data$DIMORPHISM_bodylength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree, data$display_complexity, as.numeric(as.character(data$DIMORPHISM_bodylength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree, data$display_complexity_only_1_4, as.numeric(as.character(data$DIMORPHISM_bodylength)), nsim=1000, posthoc=TRUE, p.adj="holm")

phylANOVA(tree_reduced, data_reduced$wings_display, as.numeric(as.character(data_reduced$DIMORPHISM_bodylength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree_reduced, data_reduced$arms_display, as.numeric(as.character(data_reduced$DIMORPHISM_bodylength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree_reduced, data_reduced$wings_colours, as.numeric(as.character(data_reduced$DIMORPHISM_bodylength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree_reduced, data_reduced$arms_colours, as.numeric(as.character(data_reduced$DIMORPHISM_bodylength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree_reduced, data_reduced$abdomen_colours, as.numeric(as.character(data_reduced$DIMORPHISM_bodylength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree_reduced, data_reduced$wings_eyespots, as.numeric(as.character(data_reduced$DIMORPHISM_bodylength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree_reduced, data_reduced$arms_eyespots, as.numeric(as.character(data_reduced$DIMORPHISM_bodylength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree_reduced, data_reduced$display_includes_sound_abd_or_wings, as.numeric(as.character(data_reduced$DIMORPHISM_bodylength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree_reduced, data_reduced$mouth_display, as.numeric(as.character(data_reduced$DIMORPHISM_bodylength)), nsim=1000, posthoc=TRUE, p.adj="holm")

#PhyloANOVAs with dimorphism pronotum length
phylANOVA(tree, data$behaviour_species_family, as.numeric(as.character(data$DIMORPHISM_pronlength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree, data$Region, as.numeric(as.character(data$DIMORPHISM_pronlength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree, data$World, as.numeric(as.character(data$DIMORPHISM_pronlength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree, data$primary_defense_special_resemblance_masquerade, as.numeric(as.character(data$DIMORPHISM_pronlength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree, data$masquerade_model, as.numeric(as.character(data$DIMORPHISM_pronlength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree, data$primary_defense_general_resemblance_crypsis, as.numeric(as.character(data$DIMORPHISM_pronlength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree, data$secondary_defence_present, as.numeric(as.character(data$DIMORPHISM_pronlength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree, data$display_complexity, as.numeric(as.character(data$DIMORPHISM_pronlength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree, data$display_complexity_only_1_4, as.numeric(as.character(data$DIMORPHISM_pronlength)), nsim=1000, posthoc=TRUE, p.adj="holm")

phylANOVA(tree_reduced, data_reduced$wings_display, as.numeric(as.character(data_reduced$DIMORPHISM_pronlength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree_reduced, data_reduced$arms_display, as.numeric(as.character(data_reduced$DIMORPHISM_pronlength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree_reduced, data_reduced$wings_colours, as.numeric(as.character(data_reduced$DIMORPHISM_pronlength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree_reduced, data_reduced$arms_colours, as.numeric(as.character(data_reduced$DIMORPHISM_pronlength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree_reduced, data_reduced$abdomen_colours, as.numeric(as.character(data_reduced$DIMORPHISM_pronlength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree_reduced, data_reduced$wings_eyespots, as.numeric(as.character(data_reduced$DIMORPHISM_pronlength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree_reduced, data_reduced$arms_eyespots, as.numeric(as.character(data_reduced$DIMORPHISM_pronlength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree_reduced, data_reduced$display_includes_sound_abd_or_wings, as.numeric(as.character(data_reduced$DIMORPHISM_pronlength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(tree_reduced, data_reduced$mouth_display, as.numeric(as.character(data_reduced$DIMORPHISM_pronlength)), nsim=1000, posthoc=TRUE, p.adj="holm")

#PhyloANOVAs with dimorphism forewing length
dim(data)
head(data[,c(9:19,38:42)])
dat_fwng <- as.data.frame(na.omit(data[, c(1:19, 19,41)]))
Overlap_dimorphism <-treedata(tree, dat_fwng)
Overlap_dimorphism$phy
dat_fwng_r <- as.data.frame(na.omit(data_reduced[, c(1:19,41)]))
Overlap_dimorphism_r <-treedata(tree, dat_fwng_r)
Overlap_dimorphism_r$phy
length(dat_fwng$behaviour_species_family)
length(dat_fwng_r$wings_display)

phylANOVA(Overlap_dimorphism$phy, dat_fwng$behaviour_species_family, as.numeric(as.character(dat_fwng$DIMORPHISM_forewinglength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$Region, as.numeric(as.character(dat_fwng$DIMORPHISM_forewinglength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$World, as.numeric(as.character(dat_fwng$DIMORPHISM_forewinglength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_special_resemblance_masquerade, as.numeric(as.character(dat_fwng$DIMORPHISM_forewinglength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$masquerade_model, as.numeric(as.character(dat_fwng$DIMORPHISM_forewinglength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_general_resemblance_crypsis, as.numeric(as.character(dat_fwng$DIMORPHISM_forewinglength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$secondary_defence_present, as.numeric(as.character(dat_fwng$DIMORPHISM_forewinglength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$display_complexity, as.numeric(as.character(dat_fwng$DIMORPHISM_forewinglength)), nsim=1000, posthoc=TRUE, p.adj="holm")

phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_display, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewinglength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_display, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewinglength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_colours, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewinglength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_colours, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewinglength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$abdomen_colours, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewinglength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_eyespots, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewinglength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_eyespots, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewinglength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$display_includes_sound_abd_or_wings, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewinglength)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$mouth_display, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewinglength)), nsim=1000, posthoc=TRUE, p.adj="holm")


#PhyloANOVAs with dimorphism forewing length divided by body length (RATIO)
dim(data)
head(data[,40:44])
dat_fwng <- as.data.frame(na.omit(data[, c(1:19,38)]))
Overlap_dimorphism <-treedata(tree, dat_fwng)
Overlap_dimorphism$phy
dat_fwng_r <- as.data.frame(na.omit(data_reduced[, c(1:19,38)]))
Overlap_dimorphism_r <-treedata(tree, dat_fwng_r)
Overlap_dimorphism_r$phy
length(dat_fwng$behaviour_species_family)
length(dat_fwng_r$wings_display)

phylANOVA(Overlap_dimorphism$phy, dat_fwng$behaviour_species_family, as.numeric(as.character(dat_fwng$DIMORPHISM_forewing_BD_ratio)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$Region, as.numeric(as.character(dat_fwng$DIMORPHISM_forewing_BD_ratio)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$World, as.numeric(as.character(dat_fwng$DIMORPHISM_forewing_BD_ratio)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_special_resemblance_masquerade, as.numeric(as.character(dat_fwng$DIMORPHISM_forewing_BD_ratio)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$masquerade_model, as.numeric(as.character(dat_fwng$DIMORPHISM_forewing_BD_ratio)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_general_resemblance_crypsis, as.numeric(as.character(dat_fwng$DIMORPHISM_forewing_BD_ratio)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$secondary_defence_present, as.numeric(as.character(dat_fwng$DIMORPHISM_forewing_BD_ratio)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$display_complexity, as.numeric(as.character(dat_fwng$DIMORPHISM_forewing_BD_ratio)), nsim=1000, posthoc=TRUE, p.adj="holm")

phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_display, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewing_BD_ratio)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_display, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewing_BD_ratio)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_colours, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewing_BD_ratio)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_colours, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewing_BD_ratio)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$abdomen_colours, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewing_BD_ratio)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_eyespots, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewing_BD_ratio)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_eyespots, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewing_BD_ratio)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$display_includes_sound_abd_or_wings, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewing_BD_ratio)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$mouth_display, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewing_BD_ratio)), nsim=1000, posthoc=TRUE, p.adj="holm")



#PhyloANOVAs with dimorphism forewing length divided by pronotum length (RATIO)
dim(data)
head(data[,40:44])
dat_fwng <- as.data.frame(na.omit(data[, c(1:9, 19,44)]))
Overlap_dimorphism <-treedata(tree, dat_fwng)
Overlap_dimorphism$phy
dat_fwng_r <- as.data.frame(na.omit(data_reduced[, c(10:18,44)]))
Overlap_dimorphism_r <-treedata(tree, dat_fwng_r)
Overlap_dimorphism_r$phy
length(dat_fwng$behaviour_species_family)
length(dat_fwng_r$wings_display)

phylANOVA(Overlap_dimorphism$phy, dat_fwng$behaviour_species_family, as.numeric(as.character(dat_fwng$DIMORPHISM_forewinglength_pronotum)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$Region, as.numeric(as.character(dat_fwng$DIMORPHISM_forewinglength_pronotum)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$World, as.numeric(as.character(dat_fwng$DIMORPHISM_forewinglength_pronotum)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_special_resemblance_masquerade, as.numeric(as.character(dat_fwng$DIMORPHISM_forewinglength_pronotum)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$masquerade_model, as.numeric(as.character(dat_fwng$DIMORPHISM_forewinglength_pronotum)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_general_resemblance_crypsis, as.numeric(as.character(dat_fwng$DIMORPHISM_forewinglength_pronotum)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$secondary_defence_present, as.numeric(as.character(dat_fwng$DIMORPHISM_forewinglength_pronotum)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$display_complexity, as.numeric(as.character(dat_fwng$DIMORPHISM_forewinglength_pronotum)), nsim=1000, posthoc=TRUE, p.adj="holm")

phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_display, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewinglength_pronotum)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_display, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewinglength_pronotum)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_colours, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewinglength_pronotum)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_colours, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewinglength_pronotum)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$abdomen_colours, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewinglength_pronotum)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_eyespots, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewinglength_pronotum)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_eyespots, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewinglength_pronotum)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$display_includes_sound_abd_or_wings, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewinglength_pronotum)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$mouth_display, as.numeric(as.character(dat_fwng_r$DIMORPHISM_forewinglength_pronotum)), nsim=1000, posthoc=TRUE, p.adj="holm")


#PhyloANOVAs with bodylength male
dim(data)
head(data[,22:39])
dat_fwng <- as.data.frame(na.omit(data[, c(1:19, 19,22)]))
Overlap_dimorphism <-treedata(tree, dat_fwng)
Overlap_dimorphism$phy
dat_fwng_r <- as.data.frame(na.omit(data_reduced[, c(1:19,22)]))
Overlap_dimorphism_r <-treedata(tree, dat_fwng_r)
Overlap_dimorphism_r$phy
length(dat_fwng$behaviour_species_family)
length(dat_fwng_r$wings_display)

phylANOVA(Overlap_dimorphism$phy, dat_fwng$behaviour_species_family, as.numeric(as.character(dat_fwng$bodylength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$Region, as.numeric(as.character(dat_fwng$bodylength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$World, as.numeric(as.character(dat_fwng$bodylength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_special_resemblance_masquerade, as.numeric(as.character(dat_fwng$bodylength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$masquerade_model, as.numeric(as.character(dat_fwng$bodylength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_general_resemblance_crypsis, as.numeric(as.character(dat_fwng$bodylength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$secondary_defence_present, as.numeric(as.character(dat_fwng$bodylength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$display_complexity, as.numeric(as.character(dat_fwng$bodylength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")

phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_display, as.numeric(as.character(dat_fwng_r$bodylength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_display, as.numeric(as.character(dat_fwng_r$bodylength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_colours, as.numeric(as.character(dat_fwng_r$bodylength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_colours, as.numeric(as.character(dat_fwng_r$bodylength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$abdomen_colours, as.numeric(as.character(dat_fwng_r$bodylength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_eyespots, as.numeric(as.character(dat_fwng_r$bodylength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_eyespots, as.numeric(as.character(dat_fwng_r$bodylength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$display_includes_sound_abd_or_wings, as.numeric(as.character(dat_fwng_r$bodylength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$mouth_display, as.numeric(as.character(dat_fwng_r$bodylength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")


#PhyloANOVAs with bodylength female
dim(data)
head(data[,22:39])
dat_fwng <- as.data.frame(na.omit(data[, c(1:19, 19,20)]))
Overlap_dimorphism <-treedata(tree, dat_fwng)
Overlap_dimorphism$phy
dat_fwng_r <- as.data.frame(na.omit(data_reduced[, c(1:19,20)]))
Overlap_dimorphism_r <-treedata(tree, dat_fwng_r)
Overlap_dimorphism_r$phy
length(dat_fwng$behaviour_species_family)
length(dat_fwng_r$wings_display)

phylANOVA(Overlap_dimorphism$phy, dat_fwng$behaviour_species_family, as.numeric(as.character(dat_fwng$bodylength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$Region, as.numeric(as.character(dat_fwng$bodylength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$World, as.numeric(as.character(dat_fwng$bodylength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_special_resemblance_masquerade, as.numeric(as.character(dat_fwng$bodylength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$masquerade_model, as.numeric(as.character(dat_fwng$bodylength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_general_resemblance_crypsis, as.numeric(as.character(dat_fwng$bodylength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$secondary_defence_present, as.numeric(as.character(dat_fwng$bodylength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$display_complexity, as.numeric(as.character(dat_fwng$bodylength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")

phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_display, as.numeric(as.character(dat_fwng_r$bodylength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_display, as.numeric(as.character(dat_fwng_r$bodylength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_colours, as.numeric(as.character(dat_fwng_r$bodylength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_colours, as.numeric(as.character(dat_fwng_r$bodylength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$abdomen_colours, as.numeric(as.character(dat_fwng_r$bodylength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_eyespots, as.numeric(as.character(dat_fwng_r$bodylength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_eyespots, as.numeric(as.character(dat_fwng_r$bodylength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$display_includes_sound_abd_or_wings, as.numeric(as.character(dat_fwng_r$bodylength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$mouth_display, as.numeric(as.character(dat_fwng_r$bodylength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")


#PhyloANOVAs with pronotumlength female
dim(data)
head(data[,22:39])
dat_fwng <- as.data.frame(na.omit(data[, c(1:19, 19,24)]))
Overlap_dimorphism <-treedata(tree, dat_fwng)
Overlap_dimorphism$phy
dat_fwng_r <- as.data.frame(na.omit(data_reduced[, c(1:18,24)]))
Overlap_dimorphism_r <-treedata(tree, dat_fwng_r)
Overlap_dimorphism_r$phy
length(dat_fwng$behaviour_species_family)
length(dat_fwng_r$wings_display)

phylANOVA(Overlap_dimorphism$phy, dat_fwng$behaviour_species_family, as.numeric(as.character(dat_fwng$pronlength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$Region, as.numeric(as.character(dat_fwng$pronlength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$World, as.numeric(as.character(dat_fwng$pronlength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_special_resemblance_masquerade, as.numeric(as.character(dat_fwng$pronlength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$masquerade_model, as.numeric(as.character(dat_fwng$pronlength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_general_resemblance_crypsis, as.numeric(as.character(dat_fwng$pronlength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$secondary_defence_present, as.numeric(as.character(dat_fwng$pronlength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$display_complexity, as.numeric(as.character(dat_fwng$pronlength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")

phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_display, as.numeric(as.character(dat_fwng_r$pronlength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_display, as.numeric(as.character(dat_fwng_r$pronlength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_colours, as.numeric(as.character(dat_fwng_r$pronlength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_colours, as.numeric(as.character(dat_fwng_r$pronlength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$abdomen_colours, as.numeric(as.character(dat_fwng_r$pronlength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_eyespots, as.numeric(as.character(dat_fwng_r$pronlength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_eyespots, as.numeric(as.character(dat_fwng_r$pronlength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$display_includes_sound_abd_or_wings, as.numeric(as.character(dat_fwng_r$pronlength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$mouth_display, as.numeric(as.character(dat_fwng_r$pronlength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")




#PhyloANOVAs with pronotumlength male
dim(data)
head(data[,22:26])
dat_fwng <- as.data.frame(na.omit(data[, c(1:19, 19,26)]))
Overlap_dimorphism <-treedata(tree, dat_fwng)
Overlap_dimorphism$phy
dat_fwng_r <- as.data.frame(na.omit(data_reduced[, c(1:18,26)]))
Overlap_dimorphism_r <-treedata(tree, dat_fwng_r)
Overlap_dimorphism_r$phy
length(dat_fwng$behaviour_species_family)
length(dat_fwng_r$wings_display)

phylANOVA(Overlap_dimorphism$phy, dat_fwng$behaviour_species_family, as.numeric(as.character(dat_fwng$pronlength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$Region, as.numeric(as.character(dat_fwng$pronlength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$World, as.numeric(as.character(dat_fwng$pronlength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_special_resemblance_masquerade, as.numeric(as.character(dat_fwng$pronlength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$masquerade_model, as.numeric(as.character(dat_fwng$pronlength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_general_resemblance_crypsis, as.numeric(as.character(dat_fwng$pronlength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$secondary_defence_present, as.numeric(as.character(dat_fwng$pronlength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$display_complexity, as.numeric(as.character(dat_fwng$pronlength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")

phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_display, as.numeric(as.character(dat_fwng_r$pronlength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_display, as.numeric(as.character(dat_fwng_r$pronlength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_colours, as.numeric(as.character(dat_fwng_r$pronlength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_colours, as.numeric(as.character(dat_fwng_r$pronlength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$abdomen_colours, as.numeric(as.character(dat_fwng_r$pronlength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_eyespots, as.numeric(as.character(dat_fwng_r$pronlength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_eyespots, as.numeric(as.character(dat_fwng_r$pronlength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$display_includes_sound_abd_or_wings, as.numeric(as.character(dat_fwng_r$pronlength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$mouth_display, as.numeric(as.character(dat_fwng_r$pronlength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")


#PhyloANOVAs with forewinglength female
dim(data)
head(data[,22:39])
dat_fwng <- as.data.frame(na.omit(data[, c(1:19, 19,28)]))
Overlap_dimorphism <-treedata(tree, dat_fwng)
Overlap_dimorphism$phy
dat_fwng_r <- as.data.frame(na.omit(data_reduced[, c(1:18,28)]))
Overlap_dimorphism_r <-treedata(tree, dat_fwng_r)
Overlap_dimorphism_r$phy
length(dat_fwng$behaviour_species_family)
length(dat_fwng_r$wings_display)

phylANOVA(Overlap_dimorphism$phy, dat_fwng$behaviour_species_family, as.numeric(as.character(dat_fwng$forewinglength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$Region, as.numeric(as.character(dat_fwng$forewinglength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$World, as.numeric(as.character(dat_fwng$forewinglength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_special_resemblance_masquerade, as.numeric(as.character(dat_fwng$forewinglength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$masquerade_model, as.numeric(as.character(dat_fwng$forewinglength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_general_resemblance_crypsis, as.numeric(as.character(dat_fwng$forewinglength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$secondary_defence_present, as.numeric(as.character(dat_fwng$forewinglength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$display_complexity, as.numeric(as.character(dat_fwng$forewinglength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")

phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_display, as.numeric(as.character(dat_fwng_r$forewinglength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_display, as.numeric(as.character(dat_fwng_r$forewinglength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_colours, as.numeric(as.character(dat_fwng_r$forewinglength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_colours, as.numeric(as.character(dat_fwng_r$forewinglength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$abdomen_colours, as.numeric(as.character(dat_fwng_r$forewinglength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_eyespots, as.numeric(as.character(dat_fwng_r$forewinglength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_eyespots, as.numeric(as.character(dat_fwng_r$forewinglength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$display_includes_sound_abd_or_wings, as.numeric(as.character(dat_fwng_r$forewinglength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$mouth_display, as.numeric(as.character(dat_fwng_r$forewinglength_f_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")



#PhyloANOVAs with forewinglength male
dim(data)
head(data[,22:39])
dat_fwng <- as.data.frame(na.omit(data[, c(1:19, 19,30)]))
Overlap_dimorphism <-treedata(tree, dat_fwng)
Overlap_dimorphism$phy
dat_fwng_r <- as.data.frame(na.omit(data_reduced[, c(1:18,30)]))
Overlap_dimorphism_r <-treedata(tree, dat_fwng_r)
Overlap_dimorphism_r$phy
length(dat_fwng$behaviour_species_family)
length(dat_fwng_r$wings_display)

phylANOVA(Overlap_dimorphism$phy, dat_fwng$behaviour_species_family, as.numeric(as.character(dat_fwng$forewinglength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$Region, as.numeric(as.character(dat_fwng$forewinglength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$World, as.numeric(as.character(dat_fwng$forewinglength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_special_resemblance_masquerade, as.numeric(as.character(dat_fwng$forewinglength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$masquerade_model, as.numeric(as.character(dat_fwng$forewinglength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_general_resemblance_crypsis, as.numeric(as.character(dat_fwng$forewinglength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$secondary_defence_present, as.numeric(as.character(dat_fwng$forewinglength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$display_complexity, as.numeric(as.character(dat_fwng$forewinglength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")

phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_display, as.numeric(as.character(dat_fwng_r$forewinglength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_display, as.numeric(as.character(dat_fwng_r$forewinglength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_colours, as.numeric(as.character(dat_fwng_r$forewinglength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_colours, as.numeric(as.character(dat_fwng_r$forewinglength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$abdomen_colours, as.numeric(as.character(dat_fwng_r$forewinglength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_eyespots, as.numeric(as.character(dat_fwng_r$forewinglength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_eyespots, as.numeric(as.character(dat_fwng_r$forewinglength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$display_includes_sound_abd_or_wings, as.numeric(as.character(dat_fwng_r$forewinglength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$mouth_display, as.numeric(as.character(dat_fwng_r$forewinglength_m_ave)), nsim=1000, posthoc=TRUE, p.adj="holm")



#PhyloANOVAs with forewinglength_pronotum_F female
dim(data)
head(data[,22:39])
dat_fwng <- as.data.frame(na.omit(data[, c(1:9, 19,34)]))
Overlap_dimorphism <-treedata(tree, dat_fwng)
Overlap_dimorphism$phy
dat_fwng_r <- as.data.frame(na.omit(data_reduced[, c(10:18,34)]))
Overlap_dimorphism_r <-treedata(tree, dat_fwng_r)
Overlap_dimorphism_r$phy
length(dat_fwng$behaviour_species_family)
length(dat_fwng_r$wings_display)

phylANOVA(Overlap_dimorphism$phy, dat_fwng$behaviour_species_family, as.numeric(as.character(dat_fwng$forewinglength_pronotum_F)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$Region, as.numeric(as.character(dat_fwng$forewinglength_pronotum_F)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$World, as.numeric(as.character(dat_fwng$forewinglength_pronotum_F)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_special_resemblance_masquerade, as.numeric(as.character(dat_fwng$forewinglength_pronotum_F)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$masquerade_model, as.numeric(as.character(dat_fwng$forewinglength_pronotum_F)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_general_resemblance_crypsis, as.numeric(as.character(dat_fwng$forewinglength_pronotum_F)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$secondary_defence_present, as.numeric(as.character(dat_fwng$forewinglength_pronotum_F)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$display_complexity, as.numeric(as.character(dat_fwng$forewinglength_pronotum_F)), nsim=1000, posthoc=TRUE, p.adj="holm")

phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_display, as.numeric(as.character(dat_fwng_r$forewinglength_pronotum_F)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_display, as.numeric(as.character(dat_fwng_r$forewinglength_pronotum_F)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_colours, as.numeric(as.character(dat_fwng_r$forewinglength_pronotum_F)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_colours, as.numeric(as.character(dat_fwng_r$forewinglength_pronotum_F)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$abdomen_colours, as.numeric(as.character(dat_fwng_r$forewinglength_pronotum_F)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_eyespots, as.numeric(as.character(dat_fwng_r$forewinglength_pronotum_F)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_eyespots, as.numeric(as.character(dat_fwng_r$forewinglength_pronotum_F)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$display_includes_sound_abd_or_wings, as.numeric(as.character(dat_fwng_r$forewinglength_pronotum_F)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$mouth_display, as.numeric(as.character(dat_fwng_r$forewinglength_pronotum_F)), nsim=1000, posthoc=TRUE, p.adj="holm")


#PhyloANOVAs with forewinglength_pronotum_M male
dim(data)
head(data[,22:39])
dat_fwng <- as.data.frame(na.omit(data[, c(1:19, 19,37)]))
Overlap_dimorphism <-treedata(tree, dat_fwng)
Overlap_dimorphism$phy
dat_fwng_r <- as.data.frame(na.omit(data_reduced[, c(1:18,37)]))
Overlap_dimorphism_r <-treedata(tree, dat_fwng_r)
Overlap_dimorphism_r$phy
length(dat_fwng$behaviour_species_family)
length(dat_fwng_r$wings_display)

phylANOVA(Overlap_dimorphism$phy, dat_fwng$behaviour_species_family, as.numeric(as.character(dat_fwng$forewinglength_pronotum_M)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$Region, as.numeric(as.character(dat_fwng$forewinglength_pronotum_M)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$World, as.numeric(as.character(dat_fwng$forewinglength_pronotum_M)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_special_resemblance_masquerade, as.numeric(as.character(dat_fwng$forewinglength_pronotum_M)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$masquerade_model, as.numeric(as.character(dat_fwng$forewinglength_pronotum_M)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_general_resemblance_crypsis, as.numeric(as.character(dat_fwng$forewinglength_pronotum_M)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$secondary_defence_present, as.numeric(as.character(dat_fwng$forewinglength_pronotum_M)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$display_complexity, as.numeric(as.character(dat_fwng$forewinglength_pronotum_M)), nsim=1000, posthoc=TRUE, p.adj="holm")

phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_display, as.numeric(as.character(dat_fwng_r$forewinglength_pronotum_M)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_display, as.numeric(as.character(dat_fwng_r$forewinglength_pronotum_M)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_colours, as.numeric(as.character(dat_fwng_r$forewinglength_pronotum_M)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_colours, as.numeric(as.character(dat_fwng_r$forewinglength_pronotum_M)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$abdomen_colours, as.numeric(as.character(dat_fwng_r$forewinglength_pronotum_M)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_eyespots, as.numeric(as.character(dat_fwng_r$forewinglength_pronotum_M)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_eyespots, as.numeric(as.character(dat_fwng_r$forewinglength_pronotum_M)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$display_includes_sound_abd_or_wings, as.numeric(as.character(dat_fwng_r$forewinglength_pronotum_M)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$mouth_display, as.numeric(as.character(dat_fwng_r$forewinglength_pronotum_M)), nsim=1000, posthoc=TRUE, p.adj="holm")


#PhyloANOVAs with forewinglength_bodylength_F female
dim(data)
head(data[,22:39])
dat_fwng <- as.data.frame(na.omit(data[, c(1:19, 19,36)]))
Overlap_dimorphism <-treedata(tree, dat_fwng)
Overlap_dimorphism$phy
dat_fwng_r <- as.data.frame(na.omit(data_reduced[, c(1:18,36)]))
Overlap_dimorphism_r <-treedata(tree, dat_fwng_r)
Overlap_dimorphism_r$phy
length(dat_fwng$behaviour_species_family)
length(dat_fwng_r$wings_display)

phylANOVA(Overlap_dimorphism$phy, dat_fwng$behaviour_species_family, as.numeric(as.character(dat_fwng$forewinglength_bodylength_f)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$Region, as.numeric(as.character(dat_fwng$forewinglength_bodylength_f)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$World, as.numeric(as.character(dat_fwng$forewinglength_bodylength_f)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_special_resemblance_masquerade, as.numeric(as.character(dat_fwng$forewinglength_bodylength_f)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$masquerade_model, as.numeric(as.character(dat_fwng$forewinglength_bodylength_f)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_general_resemblance_crypsis, as.numeric(as.character(dat_fwng$forewinglength_bodylength_f)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$secondary_defence_present, as.numeric(as.character(dat_fwng$forewinglength_bodylength_f)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$display_complexity, as.numeric(as.character(dat_fwng$forewinglength_bodylength_f)), nsim=1000, posthoc=TRUE, p.adj="holm")

phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_display, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_f)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_display, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_f)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_colours, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_f)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_colours, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_f)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$abdomen_colours, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_f)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_eyespots, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_f)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_eyespots, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_f)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$display_includes_sound_abd_or_wings, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_f)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$mouth_display, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_f)), nsim=1000, posthoc=TRUE, p.adj="holm")



#PhyloANOVAs with forewinglength_bodylength_M male
dim(data)
head(data[,22:39])
dat_fwng <- as.data.frame(na.omit(data[, c(1:19, 19,37)]))
Overlap_dimorphism <-treedata(tree, dat_fwng)
Overlap_dimorphism$phy
dat_fwng_r <- as.data.frame(na.omit(data_reduced[, c(1:18,37)]))
Overlap_dimorphism_r <-treedata(tree, dat_fwng_r)
Overlap_dimorphism_r$phy
length(dat_fwng$behaviour_species_family)
length(dat_fwng_r$wings_display)

phylANOVA(Overlap_dimorphism$phy, dat_fwng$behaviour_species_family, as.numeric(as.character(dat_fwng$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$Region, as.numeric(as.character(dat_fwng$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$World, as.numeric(as.character(dat_fwng$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_special_resemblance_masquerade, as.numeric(as.character(dat_fwng$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$masquerade_model, as.numeric(as.character(dat_fwng$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_general_resemblance_crypsis, as.numeric(as.character(dat_fwng$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$secondary_defence_present, as.numeric(as.character(dat_fwng$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$display_complexity, as.numeric(as.character(dat_fwng$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")

phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_display, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_display, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_colours, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_colours, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$abdomen_colours, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_eyespots, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_eyespots, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$display_includes_sound_abd_or_wings, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$mouth_display, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")


#PhyloANOVAs with forewinglength_pronotum_m male
dim(data)
head(data[,22:39])
dat_fwng <- as.data.frame(na.omit(data[, c(1:9, 19,39)]))
Overlap_dimorphism <-treedata(tree, dat_fwng)
Overlap_dimorphism$phy
dat_fwng_r <- as.data.frame(na.omit(data_reduced[, c(10:18,39)]))
Overlap_dimorphism_r <-treedata(tree, dat_fwng_r)
Overlap_dimorphism_r$phy
length(dat_fwng$behaviour_species_family)
length(dat_fwng_r$wings_display)

phylANOVA(Overlap_dimorphism$phy, dat_fwng$behaviour_species_family, as.numeric(as.character(dat_fwng$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$Region, as.numeric(as.character(dat_fwng$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$World, as.numeric(as.character(dat_fwng$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_special_resemblance_masquerade, as.numeric(as.character(dat_fwng$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$masquerade_model, as.numeric(as.character(dat_fwng$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$primary_defense_general_resemblance_crypsis, as.numeric(as.character(dat_fwng$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$secondary_defence_present, as.numeric(as.character(dat_fwng$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism$phy, dat_fwng$display_complexity, as.numeric(as.character(dat_fwng$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")

phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_display, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_display, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_colours, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_colours, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$abdomen_colours, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$wings_eyespots, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$arms_eyespots, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$display_includes_sound_abd_or_wings, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")
phylANOVA(Overlap_dimorphism_r$phy, dat_fwng_r$mouth_display, as.numeric(as.character(dat_fwng_r$forewinglength_bodylength_m)), nsim=1000, posthoc=TRUE, p.adj="holm")


#### 3.2. Phylogenetic signal ####

# Is there strong phylo singal in presence of display?

presence_display <- data$Presence_display
names(presence_display) <- row.names(data)

# Using fitDiscrete. Pagel’s lambda is a multiplier of the off-diagonal elements of a variance–covariance matrix, which best fits the distribution of 
  # data at the tips of the phylogeny. Values vary between 0 (phylogenetic independence) and 1 (traits evolve according to a Brownian process)
Presence_display_lambda <- fitDiscrete(tree, presence_display, transform="lambda")
Presence_display_lambda
Presence_display_lambda$opt$lambda

  # Randomizations of a character on a tree  and calculate a p value from https://github.com/juliema/publications/blob/master/BrueeliaMS/Maddison.Slatkin.R
  # Will tell you whether there are LESS or MORE evolutionary transitions than expected by randomization process. gives you a p-value to go with
  # that and an histogram with an arrow showing where your observed number of transitions fall agains the randomised ones.
phylo.signal.disc <- function(trait,phy,rep = 999,cost=NULL){
    lev <- attributes(factor(trait))$levels
    if (length(lev) == length(trait))
      stop("Are you sure this variable is categorical?")
    if(is.null(cost)){
      cost1 <- 1-diag(length(lev))
    }
    else {
      if (length(lev) != dim(cost)[1])
        stop("Dimensions of the character state transition matrix do not agree with the number of levels")
      cost1<- t(cost)
    }
    dimnames(cost1) <- list(lev,lev)
    trait <- as.numeric(trait)
    attributes(trait)$names <- phy$tip
    NULL.MODEL <- matrix(NA,rep,1)
    obs <- t(data.frame(trait))
    obs <- phyDat(t(obs),type="USER",levels=attributes(factor(obs))$levels)
    OBS <- parsimony(phy,obs,method="sankoff",cost=cost1)
    for (i in 1:rep){
      null <- sample(as.numeric(trait))
      attributes(null)$names <- attributes(trait)$names
      null <- t(data.frame(null))
      null <- phyDat(t(null),type="USER",levels=attributes(factor(null))$levels)
      NULL.MODEL[i,]<-parsimony(phy,null,method="sankoff",cost=cost1)
      P.value <- sum(OBS >= NULL.MODEL)/(rep + 1)
    }
    par(mfrow=c(1,2))
    hist(NULL.MODEL,xlab="Transitions.in.Randomizations",xlim=c(min(c(min(NULL.MODEL,OBS-1))),max(NULL.MODEL)+1))
    arrows(OBS,rep/10,OBS,0,angle=20,col="red",lwd=4)
    phy$tip.label <- rep(".",length(trait))
    plot(phy,tip.col=trait+10,cex=250/length(trait),font=1)
    title("Character states")
    par(mfrow=c(1,1))
    
    OUTPUT1 <- t(data.frame(Number.of.Levels = length(attributes(factor(trait))$levels), Evolutionary.Transitions.Observed=OBS,Evolutionary.Transitions.Randomization.Median=median(NULL.MODEL),Evolutionary.Transitions.Randomization.Min=min(NULL.MODEL),Evolutionary.Transitions.Randomization.Max=max(NULL.MODEL),P.value))
    
    if(is.null(cost)){
      list(.Randomization.Results=OUTPUT1,.Levels= lev,.Costs.of.character.state.transition.UNORDERED.PARSIMONY = t(cost1))
    }
    else {
      list(.Randomization.Results=OUTPUT1,.Levels= lev,.Costs.of.character.state.transition.FROM.ROW.TO.COL = t(cost1))        }
  }


# I like this function because as well as telling you whether the evolutionary trasnitions and distribution of the character across the tree is random
# (which you would know with lambda anyway), it tells you the number of evolutionary transitions, and the number of transitions that would be 
# randomly be expected (Plus you get a nice plot with it)
Presence_display <- phylo.signal.disc(presence_display, tree)
Presence_display$.Randomization.Results 


# 3. Is there strong phylogenetic signal in the World?
World <- data$World
names(World) <- row.names(data)

World_lambda <- fitDiscrete(tree, World, transform="lambda")
World_lambda


World_d <- phylo.signal.disc(World, tree)
World_d 


# 3. Is there strong phylogenetic signal in the secondary_defence_present?
secondary_defence_present <- data$secondary_defence_present
names(secondary_defence_present) <- row.names(data)

secondary_defence_present_lambda <- fitDiscrete(tree, secondary_defence_present, transform="lambda")
secondary_defence_present_lambda


secondary_defence_present_d <- phylo.signal.disc(secondary_defence_present, tree)
secondary_defence_present_d 


2# 3. Is there strong phylogenetic signal in the complexity of the display (all)?
complexity_display <- data$display_complexity
names(complexity_display) <- row.names(data)

Complexity_display_lambda <- fitDiscrete(tree, complexity_display, transform="lambda")
Complexity_display_lambda
Complexity_display_lambda$opt$lambda 

Complexity_display_d <- phylo.signal.disc(complexity_display, tree)
Complexity_display_d 


# 3. Is there strong phylogenetic signal in the complexity of the display (between 1 and 4)? 
Model_3_data <- as.data.frame(na.omit(cbind(row.names(data), as.numeric(as.character(data$display_complexity)))))
row.names(Model_3_data)<-Model_3_data$V1
Overlap_3 <-treedata(tree, Model_3_data)


complexity_display <- Model_3_data$V2
names(complexity_display) <- row.names(Model_3_data)

Complexity_display_lambda <- fitDiscrete(Overlap_3$phy, complexity_display, transform="lambda")
Complexity_display_lambda
Complexity_display_lambda$opt$lambda 

Complexity_display <- phylo.signal.disc(complexity_display, Overlap_3$phy)
Complexity_display 


# Is there strong phylogenetic signal in primary_defense_special_resemblance? 
primary_defense_special_resemblance <- data$primary_defense_special_resemblance_masquerade

names(primary_defense_special_resemblance) <- row.names(data)

primary_defense_special_resemblance_lambda <- fitDiscrete(tree, primary_defense_special_resemblance, transform="lambda")
primary_defense_special_resemblance_lambda
primary_defense_special_resemblance_lambda$opt$lambda 

Primary_defense_special_resemblance <- phylo.signal.disc(primary_defense_special_resemblance, tree)
Primary_defense_special_resemblance$.Randomization.Results 

# Is there strong phylogenetic signal in primary_defense_general_resemblance? 
primary_defense_general_resemblance <- data$primary_defense_general_resemblance_crypsis
names(primary_defense_general_resemblance) <- row.names(data)

primary_defense_general_resemblance_lambda <- fitDiscrete(tree, primary_defense_general_resemblance, transform="lambda")
primary_defense_general_resemblance_lambda 

Primary_defense_general_resemblance <- phylo.signal.disc(primary_defense_general_resemblance, tree)
Primary_defense_general_resemblance$.Randomization.Results 
  
head(data)

head(data_reduced)

# REDUCED DATASET
# data_reduced$wings_display
wings_display <- data_reduced$wings_display
names(wings_display) <- row.names(data_reduced)

wings_display_lambda <- fitDiscrete(tree_reduced, wings_display, transform="lambda")
wings_display_lambda


wings_display_d <- phylo.signal.disc(wings_display, tree_reduced)
wings_display_d 

# data_reduced$arms_display
arms_display <- data_reduced$arms_display
names(arms_display) <- row.names(data_reduced)

arms_display_lambda <- fitDiscrete(tree_reduced, arms_display, transform="lambda")
arms_display_lambda


arms_display_d <- phylo.signal.disc(arms_display, tree_reduced)
arms_display_d 


#2 data_reduced$wings_colours
wings_colours <- data_reduced$wings_colours
names(wings_colours) <- row.names(data_reduced)

wings_colours_lambda <- fitDiscrete(tree_reduced, wings_colours, transform="lambda")
wings_colours_lambda


wings_colours_d <- phylo.signal.disc(wings_colours, tree_reduced)
wings_colours_d 

# data_reduced$arms_colours
arms_colours <- data_reduced$arms_colours
names(arms_colours) <- row.names(data_reduced)

arms_colours_lambda <- fitDiscrete(tree_reduced, arms_colours, transform="lambda")
arms_colours_lambda


arms_colours_d <- phylo.signal.disc(arms_colours, tree_reduced)
arms_colours_d 

# data_reduced$abdomen_colours

abdomen_colours <- data_reduced$abdomen_colours
names(abdomen_colours) <- row.names(data_reduced)

abdomen_colours_lambda <- fitDiscrete(tree_reduced, abdomen_colours, transform="lambda")
abdomen_colours_lambda


abdomen_colours_d <- phylo.signal.disc(abdomen_colours, tree_reduced)
abdomen_colours_d 


# data_reduced$wings_eyespots
wings_eyespots <- data_reduced$wings_eyespots
names(wings_eyespots) <- row.names(data_reduced)

wings_eyespots_lambda <- fitDiscrete(tree_reduced, wings_eyespots, transform="lambda")
wings_eyespots_lambda


wings_eyespots_d <- phylo.signal.disc(wings_eyespots, tree_reduced)
wings_eyespots_d 

# data_reduced$arms_eyespots
arms_eyespots <- data_reduced$arms_eyespots
names(arms_eyespots) <- row.names(data_reduced)

arms_eyespots_lambda <- fitDiscrete(tree_reduced, arms_eyespots, transform="lambda")
arms_eyespots_lambda


arms_eyespots_d <- phylo.signal.disc(arms_eyespots, tree_reduced)
arms_eyespots_d 

# data_reduced$display_includes_sound_abd_or_wings
display_includes_sound_abd_or_wings <- data_reduced$display_includes_sound_abd_or_wings
names(display_includes_sound_abd_or_wings) <- row.names(data_reduced)

display_includes_sound_abd_or_wings_lambda <- fitDiscrete(tree_reduced, display_includes_sound_abd_or_wings, transform="lambda")
display_includes_sound_abd_or_wings_lambda


display_includes_sound_abd_or_wings_d <- phylo.signal.disc(display_includes_sound_abd_or_wings, tree_reduced)
display_includes_sound_abd_or_wings_d 

# data_reduced$mouth_display
mouth_display <- data_reduced$mouth_display
names(mouth_display) <- row.names(data_reduced)

mouth_display_lambda <- fitDiscrete(tree_reduced, mouth_display, transform="lambda")
mouth_display_lambda


mouth_display_d <- phylo.signal.disc(mouth_display, tree_reduced)
mouth_display_d 


# fitContinuous
head(data[,22:39])
head(data[,40:44])

fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 22)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 24)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 26)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 28)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 30)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 32)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 34)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 35)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 36)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 37)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 38)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 39)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 40)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 41)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 42)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 43)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 44)]))

Overlap_fC <-treedata(tree, fC_bodylength_f_ave)
Overlap_fC$phy
bodylength_f_ave <- Overlap_fC$data
fC <- as.numeric(as.character(Overlap_fC$data[,2]))
names(fC) <-Overlap_fC$data[,1]

fitContinuous(phy = Overlap_fC$phy, dat = fC, model = "lambda", bounds=c(exp(-500), 1))

?phylosig

fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 22,23)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 24,25)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 26,27)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 28,29)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 30,31)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 32,33)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 34)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 35)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 36)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 37)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 38)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 39)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 40)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 41)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 42)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 43)]))
fC_bodylength_f_ave <- as.data.frame(na.omit(data[, c(2, 44)]))


Overlap_fC <-treedata(tree, fC_bodylength_f_ave)
Overlap_fC$phy
bodylength_f_ave <- Overlap_fC$data
fC <- as.numeric(as.character(Overlap_fC$data[,2]))
names(fC) <-Overlap_fC$data[,1]

phylosig(Overlap_fC$phy, fC, method="K", test=TRUE, nsim=1000, se=NULL)


SE_vector <- as.numeric(as.character(Overlap_fC$data[,3]))
names(SE_vector) <-Overlap_fC$data[,1]

phylosig(Overlap_fC$phy, fC, method="K", test=TRUE, nsim=1000, se=SE_vector)


#### 3.3 Phylogenetic ANOVAs ####

head(data)

dim(data)
head(data[,21:43])

Species <- row.names(data)
data_with_evol_d <- as.data.frame(cbind(data, evol_d, Species))

str(data_with_evol_d)

PGLS_data <- as.data.frame(na.omit(data_with_evol_d[, c(2, 21, 23, 25, 27, 29, 31, 33:45)]))
dim(PGLS_data) # 50 x 18

Overlap_PGLS <-treedata(tree, PGLS_data)
Overlap_PGLS$phy
head(Overlap_PGLS$data)
PGLS_tree <- Overlap_PGLS$phy

PGLS_data_no_order <- as.data.frame(Overlap_PGLS$data)
PGLS_data_order <- PGLS_data_no_order[match(PGLS_tree$tip.label,rownames(PGLS_data_no_order)),] 

PGLS_bodylength_f_ave <- as.numeric(as.character(PGLS_data_order[,2]))
names(PGLS_bodylength_f_ave) <- row.names(PGLS_data_order)

PGLS_bodylength_m_ave <- as.numeric(as.character(PGLS_data_order[,3]))
names(PGLS_bodylength_m_ave) <- row.names(PGLS_data_order)

PGLS_pronlength_f_ave <- as.numeric(as.character(PGLS_data_order[,4]))
names(PGLS_pronlength_f_ave) <- row.names(PGLS_data_order)

PGLS_pronlength_m_ave <- as.numeric(as.character(PGLS_data_order[,5]))
names(PGLS_pronlength_m_ave) <- row.names(PGLS_data_order)

PGLS_forewinglength_f_ave <- as.numeric(as.character(PGLS_data_order[,6]))
names(PGLS_forewinglength_f_ave) <- row.names(PGLS_data_order)

PGLS_forewinglength_m_ave <- as.numeric(as.character(PGLS_data_order[,7]))
names(PGLS_forewinglength_m_ave) <- row.names(PGLS_data_order)

PGLS_forewinglength_pronotum_F <- as.numeric(as.character(PGLS_data_order[,8]))
names(PGLS_forewinglength_pronotum_F) <- row.names(PGLS_data_order)

PGLS_forewinglength_pronotum_M <- as.numeric(as.character(PGLS_data_order[,9]))
names(PGLS_forewinglength_pronotum_M) <- row.names(PGLS_data_order)

PGLS_bodylength_pronotum_F <- as.numeric(as.character(PGLS_data_order[,10]))
names(PGLS_bodylength_pronotum_F) <- row.names(PGLS_data_order)

PGLS_bodylength_pronotum_M <- as.numeric(as.character(PGLS_data_order[,11]))
names(PGLS_bodylength_pronotum_M) <- row.names(PGLS_data_order)

PGLS_forewinglength_bodylength_f <- as.numeric(as.character(PGLS_data_order[,12]))
names(PGLS_forewinglength_bodylength_f) <- row.names(PGLS_data_order)

PGLS_forewinglength_bodylength_m <- as.numeric(as.character(PGLS_data_order[,13]))
names(PGLS_forewinglength_bodylength_m) <- row.names(PGLS_data_order)

PGLS_DIMORPHISM_forewing_BD_ratio <- as.numeric(as.character(PGLS_data_order[,14]))
names(PGLS_DIMORPHISM_forewing_BD_ratio) <- row.names(PGLS_data_order)

PGLS_DIMORPHISM_bodylength <- as.numeric(as.character(PGLS_data_order[,15]))
names(PGLS_DIMORPHISM_bodylength) <- row.names(PGLS_data_order)

PGLS_DIMORPHISM_pronlength <- as.numeric(as.character(PGLS_data_order[,16]))
names(PGLS_DIMORPHISM_pronlength) <- row.names(PGLS_data_order)


PGLS_DIMORPHISM_forewinglength <-as.numeric(as.character(PGLS_data_order[,17]))
names(PGLS_DIMORPHISM_forewinglength) <- row.names(PGLS_data_order)


PGLS_DIMORPHISM_forewinglength_pronotum <- as.numeric(as.character(PGLS_data_order[,18]))
names(PGLS_DIMORPHISM_forewinglength_pronotum) <- row.names(PGLS_data_order)


EVOL_DIS <- as.numeric(as.character(PGLS_data_order[,19]))
names(EVOL_DIS) <- row.names(PGLS_data_order)


?comparative.data

PGLS_comparative_dataset <- comparative.data(Overlap_PGLS$phy, PGLS_data_no_order, Species, vcv=TRUE, vcv.dim=3)

model_1 <- pgls(log(PGLS_bodylength_f_ave) ~ log(EVOL_DIS),  PGLS_comparative_dataset, lambda='ML')
summary(model_1)
coef(model_1)

model_2 <- pgls(log(PGLS_bodylength_m_ave) ~ log(EVOL_DIS),  PGLS_comparative_dataset, lambda='ML')
summary(model_2)

model_3 <- pgls(log(PGLS_pronlength_f_ave) ~ log(EVOL_DIS),  PGLS_comparative_dataset, lambda='ML')
summary(model_3)

model_4 <- pgls(log(PGLS_pronlength_m_ave) ~ log(EVOL_DIS),  PGLS_comparative_dataset, lambda='ML')
summary(model_4)

model_5 <- pgls(log(PGLS_forewinglength_f_ave) ~ log(EVOL_DIS),  PGLS_comparative_dataset, lambda='ML')
summary(model_5)

model_6 <- pgls(log(PGLS_forewinglength_m_ave) ~ log(EVOL_DIS),  PGLS_comparative_dataset, lambda='ML')
summary(model_6)

model_7 <- pgls(log(PGLS_forewinglength_pronotum_F) ~ log(EVOL_DIS),  PGLS_comparative_dataset, lambda='ML')
summary(model_7)

model_8 <- pgls(log(PGLS_forewinglength_pronotum_M) ~ log(EVOL_DIS),  PGLS_comparative_dataset, lambda='ML')
summary(model_8)

model_9 <- pgls(log(PGLS_bodylength_pronotum_F) ~ log(EVOL_DIS),  PGLS_comparative_dataset, lambda='ML')
summary(model_9)

model_10 <- pgls(log(PGLS_bodylength_pronotum_M) ~ log(EVOL_DIS),  PGLS_comparative_dataset, lambda='ML')
summary(model_10)

model_11 <- pgls(log(PGLS_forewinglength_bodylength_f) ~ log(EVOL_DIS),  PGLS_comparative_dataset, lambda='ML')
summary(model_11)

model_12 <- pgls(log(PGLS_forewinglength_bodylength_m) ~ log(EVOL_DIS),  PGLS_comparative_dataset, lambda='ML')
summary(model_12)

model_13 <- pgls(PGLS_DIMORPHISM_forewing_BD_ratio ~ log(EVOL_DIS),  PGLS_comparative_dataset, lambda='ML')
summary(model_13)

model_14 <- pgls(PGLS_DIMORPHISM_bodylength ~ log(EVOL_DIS),  PGLS_comparative_dataset, lambda='ML')
summary(model_14)

model_15 <- pgls(PGLS_DIMORPHISM_pronlength ~ log(EVOL_DIS),  PGLS_comparative_dataset, lambda='ML')
summary(model_15)

model_16 <- pgls(PGLS_DIMORPHISM_forewinglength ~ log(EVOL_DIS),  PGLS_comparative_dataset, lambda='ML')
summary(model_16)

model_17 <- pgls(PGLS_DIMORPHISM_forewinglength_pronotum ~ log(EVOL_DIS),  PGLS_comparative_dataset, lambda='ML')
summary(model_17)


# Not included in the paper
# procD.pgls(log(PGLS_bodylength_f_ave) ~ log(EVOL_DIS),correlation=bm.PGLS, PGLS_tree, iter = 999)
# procD.pgls(log(PGLS_bodylength_m_ave) ~ log(EVOL_DIS),correlation=bm.PGLS, PGLS_tree, iter = 999)
# procD.pgls(log(PGLS_pronlength_f_ave) ~ log(EVOL_DIS),correlation=bm.PGLS, PGLS_tree, iter = 999)
# procD.pgls(log(PGLS_pronlength_m_ave) ~ log(EVOL_DIS),correlation=bm.PGLS, PGLS_tree, iter = 999)
# procD.pgls(log(PGLS_forewinglength_f_ave) ~ log(EVOL_DIS),correlation=bm.PGLS, PGLS_tree, iter = 999)
# procD.pgls(log(PGLS_forewinglength_m_ave) ~ log(EVOL_DIS),correlation=bm.PGLS, PGLS_tree, iter = 999)
# procD.pgls(log(PGLS_forewinglength_pronotum_F) ~ log(EVOL_DIS),correlation=bm.PGLS, PGLS_tree, iter = 999)
# procD.pgls(log(PGLS_forewinglength_pronotum_M) ~ log(EVOL_DIS),correlation=bm.PGLS, PGLS_tree, iter = 999)
# procD.pgls(log(PGLS_bodylength_pronotum_F) ~ log(EVOL_DIS),correlation=bm.PGLS, PGLS_tree, iter = 999)
# procD.pgls(log(PGLS_bodylength_pronotum_M) ~ log(EVOL_DIS),correlation=bm.PGLS, PGLS_tree, iter = 999)
# procD.pgls(log(PGLS_forewinglength_bodylength_f) ~ log(EVOL_DIS),correlation=bm.PGLS, PGLS_tree, iter = 999)
# procD.pgls(log(PGLS_forewinglength_bodylength_m) ~ log(EVOL_DIS),correlation=bm.PGLS, PGLS_tree, iter = 999)
# procD.pgls(PGLS_DIMORPHISM_forewing_BD_ratio ~ log(EVOL_DIS),correlation=bm.PGLS, PGLS_tree, iter = 999)
# procD.pgls(PGLS_DIMORPHISM_bodylength ~ log(EVOL_DIS),correlation=bm.PGLS, PGLS_tree, iter = 999)
# procD.pgls(PGLS_DIMORPHISM_pronlength ~ log(EVOL_DIS),correlation=bm.PGLS, PGLS_tree, iter = 999)
# procD.pgls(PGLS_DIMORPHISM_forewinglength ~ log(EVOL_DIS),correlation=bm.PGLS, PGLS_tree, iter = 999)
# procD.pgls(PGLS_bodylength_f_ave ~ log(EVOL_DIS),correlation=bm.PGLS, PGLS_tree, iter = 999)




#### 3.4 Phylogenetic Generalized Least Squares - BROWNIAN MOTION ####

bm.PGLS<-corBrownian(phy=PGLS_tree)

bm.gls_1<-gls(log(PGLS_bodylength_f_ave) ~ log(EVOL_DIS),correlation=bm.PGLS)
summary(bm.gls_1)
plot(log(PGLS_bodylength_f_ave) ~ log(EVOL_DIS))
abline(a = coef(ou.gls_1)[1], b = coef(ou.gls_1)[2], col = "green")
abline(a = coef(bm.gls_1)[1], b = coef(bm.gls_1)[2], col = "blue")

bm.gls_zero_1 <- gls(log(PGLS_bodylength_f_ave)~1, correlation=bm.PGLS)
1-(as.numeric((logLik(bm.gls_1)-1)/logLik(bm.gls_zero_1)))


bm.gls_1<-gls(log(PGLS_bodylength_m_ave) ~ log(EVOL_DIS),correlation=bm.PGLS)
summary(bm.gls_1)
plot(log(PGLS_bodylength_m_ave) ~ log(EVOL_DIS))
abline(a = coef(bm.gls_1)[1], b = coef(bm.gls_1)[2])
bm.gls_zero_1 <- gls(log(PGLS_bodylength_m_ave)~1, correlation=bm.PGLS)
1-(as.numeric((logLik(bm.gls_1)-1)/logLik(bm.gls_zero_1)))


bm.gls_2<-gls(log(PGLS_pronlength_f_ave)~ log(EVOL_DIS),correlation=bm.PGLS)
summary(bm.gls_2)
plot(log(PGLS_pronlength_f_ave) ~ log(EVOL_DIS))
abline(a = coef(bm.gls_2)[1], b = coef(bm.gls_2)[2])
bm.gls_zero_2 <- gls(log(PGLS_pronlength_f_ave)~1, correlation=bm.PGLS)
1-(as.numeric((logLik(bm.gls_2)-1)/logLik(bm.gls_zero_2)))

bm.gls_2<-gls(log(PGLS_pronlength_m_ave) ~ log(EVOL_DIS),correlation=bm.PGLS)
summary(bm.gls_2)
plot(log(PGLS_pronlength_m_ave) ~ log(EVOL_DIS))
abline(a = coef(bm.gls_2)[1], b = coef(bm.gls_2)[2])
bm.gls_zero_2 <- gls(log(PGLS_pronlength_m_ave)~1, correlation=bm.PGLS)
1-(as.numeric((logLik(bm.gls_2)-1)/logLik(bm.gls_zero_2)))

bm.gls_3<-gls(log(PGLS_forewinglength_f_ave)~ log(EVOL_DIS),correlation=bm.PGLS)
summary(bm.gls_3)
plot(log(PGLS_forewinglength_f_ave) ~ log(EVOL_DIS))
abline(a = coef(bm.gls_3)[1], b = coef(bm.gls_3)[2])
bm.gls_zero_3 <- gls(log(PGLS_forewinglength_f_ave)~1, correlation=bm.PGLS)
1-(as.numeric((logLik(bm.gls_2)-1)/logLik(bm.gls_zero_3)))


bm.gls_3<-gls(log(PGLS_forewinglength_m_ave)~ log(EVOL_DIS),correlation=bm.PGLS)
summary(bm.gls_3)
plot(log(PGLS_forewinglength_m_ave) ~ log(EVOL_DIS))
abline(a = coef(bm.gls_3)[1], b = coef(bm.gls_3)[2])
bm.gls_zero_3 <- gls(log(PGLS_forewinglength_m_ave)~1, correlation=bm.PGLS)
1-(as.numeric((logLik(bm.gls_3)-1)/logLik(bm.gls_zero_3)))


bm.gls_4<-gls(log(PGLS_forewinglength_pronotum_F) ~ log(EVOL_DIS),correlation=bm.PGLS)
summary(bm.gls_4)
plot(log(PGLS_forewinglength_pronotum_F) ~ log(EVOL_DIS))
abline(a = coef(bm.gls_4)[1], b = coef(bm.gls_4)[2])
bm.gls_zero_4 <- gls(log(PGLS_forewinglength_pronotum_F)~1, correlation=bm.PGLS)
1-(as.numeric((logLik(bm.gls_4)-1)/logLik(bm.gls_zero_4)))


bm.gls_4<-gls(log(PGLS_forewinglength_pronotum_M)~ log(EVOL_DIS),correlation=bm.PGLS)
summary(bm.gls_4)
plot(log(PGLS_forewinglength_pronotum_M) ~ log(EVOL_DIS))
abline(a = coef(bm.gls_4)[1], b = coef(bm.gls_4)[2])
bm.gls_zero_4 <- gls(log(PGLS_forewinglength_pronotum_M)~1, correlation=bm.PGLS)
1-(as.numeric((logLik(bm.gls_4)-1)/logLik(bm.gls_zero_4)))

bm.gls_5<-gls(log(PGLS_bodylength_pronotum_F) ~ log(EVOL_DIS),correlation=bm.PGLS)
summary(bm.gls_5)
plot(log(PGLS_bodylength_pronotum_F) ~ log(EVOL_DIS))
abline(a = coef(bm.gls_5)[1], b = coef(bm.gls_5)[2])
bm.gls_zero_5 <- gls(log(PGLS_bodylength_pronotum_F)~1, correlation=bm.PGLS)
1-(as.numeric((logLik(bm.gls_5)-1)/logLik(bm.gls_zero_5)))


bm.gls_5<-gls(log(PGLS_bodylength_pronotum_M)~ log(EVOL_DIS),correlation=bm.PGLS)
summary(bm.gls_5)
plot(log(PGLS_bodylength_pronotum_M) ~ log(EVOL_DIS))
abline(a = coef(bm.gls_5)[1], b = coef(bm.gls_5)[2])
bm.gls_zero_5 <- gls(log(PGLS_bodylength_pronotum_M)~1, correlation=bm.PGLS)
1-(as.numeric((logLik(bm.gls_5)-1)/logLik(bm.gls_zero_5)))


bm.gls_6<-gls(log(PGLS_forewinglength_bodylength_f) ~ log(EVOL_DIS),correlation=bm.PGLS)
summary(bm.gls_6)
plot(log(PGLS_forewinglength_bodylength_f) ~ log(EVOL_DIS))
abline(a = coef(bm.gls_6)[1], b = coef(bm.gls_6)[2])
bm.gls_zero_6 <- gls(log(PGLS_forewinglength_bodylength_f)~1, correlation=bm.PGLS)
1-(as.numeric((logLik(bm.gls_6)-1)/logLik(bm.gls_zero_6)))


bm.gls_6<-gls(log(PGLS_forewinglength_bodylength_m)~ log(EVOL_DIS),correlation=bm.PGLS)
summary(bm.gls_6)
plot(log(PGLS_forewinglength_bodylength_m) ~ log(EVOL_DIS))
abline(a = coef(bm.gls_6)[1], b = coef(bm.gls_6)[2])
bm.gls_zero_6 <- gls(log(PGLS_forewinglength_bodylength_m)~1, correlation=bm.PGLS)
1-(as.numeric((logLik(bm.gls_6)-1)/logLik(bm.gls_zero_6)))


bm.gls_7<-gls(PGLS_DIMORPHISM_forewing_BD_ratio ~ log(EVOL_DIS),correlation=bm.PGLS)
summary(bm.gls_7)
plot(PGLS_DIMORPHISM_forewing_BD_ratio ~ log(EVOL_DIS))
abline(a = coef(bm.gls_7)[1], b = coef(bm.gls_7)[2])
bm.gls_zero_7 <- gls(PGLS_DIMORPHISM_forewing_BD_ratio~1, correlation=bm.PGLS)
1-(as.numeric((logLik(bm.gls_7)-1)/logLik(bm.gls_zero_7)))

bm.gls_8<-gls(PGLS_DIMORPHISM_bodylength~ log(EVOL_DIS),correlation=bm.PGLS)
summary(bm.gls_8)
plot(PGLS_DIMORPHISM_bodylength ~ log(EVOL_DIS))
abline(a = coef(bm.gls_8)[1], b = coef(bm.gls_8)[2])
bm.gls_zero_8 <- gls(PGLS_DIMORPHISM_bodylength~1, correlation=bm.PGLS)
1-(as.numeric((logLik(bm.gls_8)-1)/logLik(bm.gls_zero_8)))


bm.gls_9<-gls(PGLS_DIMORPHISM_pronlength ~ log(EVOL_DIS),correlation=bm.PGLS)
summary(bm.gls_9)
plot(PGLS_DIMORPHISM_pronlength ~ log(EVOL_DIS))
abline(a = coef(bm.gls_9)[1], b = coef(bm.gls_9)[2])
bm.gls_zero_9 <- gls(PGLS_DIMORPHISM_pronlength~1, correlation=bm.PGLS)
1-(as.numeric((logLik(bm.gls_9)-1)/logLik(bm.gls_zero_9)))


bm.gls_10<-gls(PGLS_DIMORPHISM_forewinglength ~ log(EVOL_DIS),correlation=bm.PGLS)
summary(bm.gls_10)
plot(PGLS_DIMORPHISM_forewinglength ~ log(EVOL_DIS))
abline(a = coef(bm.gls_10)[1], b = coef(bm.gls_10)[2])
bm.gls_zero_10 <- gls(PGLS_DIMORPHISM_forewinglength~1, correlation=bm.PGLS)
1-(as.numeric((logLik(bm.gls_10)-1)/logLik(bm.gls_zero_10)))

bm.gls_11<-gls(PGLS_DIMORPHISM_forewinglength_pronotum ~ log(EVOL_DIS),correlation=bm.PGLS)
summary(bm.gls_11)
plot(PGLS_DIMORPHISM_forewinglength_pronotum ~ log(EVOL_DIS))
abline(a = coef(bm.gls_11)[1], b = coef(bm.gls_11)[2])
bm.gls_zero_11 <- gls(PGLS_DIMORPHISM_forewinglength_pronotum~1, correlation=bm.PGLS)
1-(as.numeric((logLik(bm.gls_11)-1)/logLik(bm.gls_zero_11)))



#### 3.5 Phylogenetic Generalized Least Squares - Ornstein–Uhlenbeck ####
ou.PGLS<-corMartins(7.1, phy=PGLS_tree, fixed=TRUE)

ou.gls_1<-gls(log(PGLS_bodylength_f_ave) ~ log(EVOL_DIS),correlation=ou.PGLS)
summary(ou.gls_1)
plot(log(PGLS_bodylength_f_ave) ~ log(EVOL_DIS))
abline(a = coef(ou.gls_1)[1], b = coef(ou.gls_1)[2], col = "green")
abline(a = coef(bm.gls_1)[1], b = coef(bm.gls_1)[2], col = "blue")

ou.gls_zero_1 <- gls(log(PGLS_bodylength_f_ave)~1, correlation=ou.PGLS)
1-(as.numeric((logLik(ou.gls_1)-1)/logLik(ou.gls_zero_1)))


ou.gls_1<-gls(log(PGLS_bodylength_m_ave) ~ log(EVOL_DIS),correlation=ou.PGLS)
summary(ou.gls_1)
plot(log(PGLS_bodylength_m_ave) ~ log(EVOL_DIS))
abline(a = coef(ou.gls_1)[1], b = coef(ou.gls_1)[2], col = "green")
abline(a = coef(bm.gls_1)[1], b = coef(bm.gls_1)[2], col = "blue")
ou.gls_zero_1 <- gls(log(PGLS_bodylength_m_ave)~1, correlation=ou.PGLS)
1-(as.numeric((logLik(ou.gls_1)-1)/logLik(ou.gls_zero_1)))


ou.gls_2<-gls(log(PGLS_pronlength_f_ave)~ log(EVOL_DIS),correlation=ou.PGLS)
summary(ou.gls_2)
plot(log(PGLS_pronlength_f_ave) ~ log(EVOL_DIS))
abline(a = coef(ou.gls_2)[1], b = coef(ou.gls_2)[2], col = "green")
abline(a = coef(bm.gls_2)[1], b = coef(bm.gls_2)[2], col = "blue")
ou.gls_zero_2 <- gls(log(PGLS_pronlength_f_ave)~1, correlation=ou.PGLS)
1-(as.numeric((logLik(ou.gls_2)-1)/logLik(ou.gls_zero_2)))

ou.gls_2<-gls(log(PGLS_pronlength_m_ave) ~ log(EVOL_DIS),correlation=ou.PGLS)
summary(ou.gls_2)
plot(log(PGLS_pronlength_m_ave) ~ log(EVOL_DIS))
abline(a = coef(ou.gls_2)[1], b = coef(ou.gls_2)[2], col = "green")
abline(a = coef(bm.gls_2)[1], b = coef(bm.gls_2)[2], col = "blue")
ou.gls_zero_2 <- gls(log(PGLS_pronlength_m_ave)~1, correlation=ou.PGLS)
1-(as.numeric((logLik(ou.gls_2)-1)/logLik(ou.gls_zero_2)))

ou.gls_3<-gls(log(PGLS_forewinglength_f_ave)~ log(EVOL_DIS),correlation=ou.PGLS)
summary(ou.gls_3)
plot(log(PGLS_forewinglength_f_ave) ~ log(EVOL_DIS))
abline(a = coef(ou.gls_3)[1], b = coef(ou.gls_3)[2], col = "green")
abline(a = coef(bm.gls_3)[1], b = coef(bm.gls_3)[2], col = "blue")
ou.gls_zero_3 <- gls(log(PGLS_forewinglength_f_ave)~1, correlation=ou.PGLS)
1-(as.numeric((logLik(ou.gls_2)-1)/logLik(ou.gls_zero_3)))


ou.gls_3<-gls(log(PGLS_forewinglength_m_ave)~ log(EVOL_DIS),correlation=ou.PGLS)
summary(ou.gls_3)
plot(log(PGLS_forewinglength_m_ave) ~ log(EVOL_DIS))
abline(a = coef(ou.gls_3)[1], b = coef(ou.gls_3)[2], col = "green")
abline(a = coef(bm.gls_3)[1], b = coef(bm.gls_3)[2], col = "blue")
ou.gls_zero_3 <- gls(log(PGLS_forewinglength_m_ave)~1, correlation=ou.PGLS)
1-(as.numeric((logLik(ou.gls_3)-1)/logLik(ou.gls_zero_3)))


ou.gls_4<-gls(log(PGLS_forewinglength_pronotum_F) ~ log(EVOL_DIS),correlation=ou.PGLS)
summary(ou.gls_4)
plot(log(PGLS_forewinglength_pronotum_F) ~ log(EVOL_DIS))
abline(a = coef(ou.gls_4)[1], b = coef(ou.gls_4)[2], col = "green")
abline(a = coef(bm.gls_4)[1], b = coef(bm.gls_4)[2], col = "blue")
ou.gls_zero_4 <- gls(log(PGLS_forewinglength_pronotum_F)~1, correlation=ou.PGLS)
1-(as.numeric((logLik(ou.gls_4)-1)/logLik(ou.gls_zero_4)))


ou.gls_4<-gls(log(PGLS_forewinglength_pronotum_M)~ log(EVOL_DIS),correlation=ou.PGLS)
summary(ou.gls_4)
plot(log(PGLS_forewinglength_pronotum_M) ~ log(EVOL_DIS))
abline(a = coef(ou.gls_4)[1], b = coef(ou.gls_4)[2], col = "green")
abline(a = coef(bm.gls_4)[1], b = coef(bm.gls_4)[2], col = "blue")
ou.gls_zero_4 <- gls(log(PGLS_forewinglength_pronotum_M)~1, correlation=ou.PGLS)
1-(as.numeric((logLik(ou.gls_4)-1)/logLik(ou.gls_zero_4)))

ou.gls_5<-gls(log(PGLS_bodylength_pronotum_F) ~ log(EVOL_DIS),correlation=ou.PGLS)
summary(ou.gls_5)
plot(log(PGLS_bodylength_pronotum_F) ~ log(EVOL_DIS))
abline(a = coef(ou.gls_5)[1], b = coef(ou.gls_5)[2], col = "green")
abline(a = coef(bm.gls_5)[1], b = coef(bm.gls_5)[2], col = "blue")
ou.gls_zero_5 <- gls(log(PGLS_bodylength_pronotum_F)~1, correlation=ou.PGLS)
1-(as.numeric((logLik(ou.gls_5)-1)/logLik(ou.gls_zero_5)))


ou.gls_5<-gls(log(PGLS_bodylength_pronotum_M)~ log(EVOL_DIS),correlation=ou.PGLS)
summary(ou.gls_5)
plot(log(PGLS_bodylength_pronotum_M) ~ log(EVOL_DIS))
abline(a = coef(ou.gls_5)[1], b = coef(ou.gls_5)[2], col = "green")
abline(a = coef(bm.gls_5)[1], b = coef(bm.gls_5)[2], col = "blue")
ou.gls_zero_5 <- gls(log(PGLS_bodylength_pronotum_M)~1, correlation=ou.PGLS)
1-(as.numeric((logLik(ou.gls_5)-1)/logLik(ou.gls_zero_5)))


ou.gls_6<-gls(log(PGLS_forewinglength_bodylength_f) ~ log(EVOL_DIS),correlation=ou.PGLS)
summary(ou.gls_6)
plot(log(PGLS_forewinglength_bodylength_f) ~ log(EVOL_DIS))
abline(a = coef(ou.gls_6)[1], b = coef(ou.gls_6)[2], col = "green")
abline(a = coef(bm.gls_6)[1], b = coef(bm.gls_6)[2], col = "blue")
ou.gls_zero_6 <- gls(log(PGLS_forewinglength_bodylength_f)~1, correlation=ou.PGLS)
1-(as.numeric((logLik(ou.gls_6)-1)/logLik(ou.gls_zero_6)))


ou.gls_6<-gls(log(PGLS_forewinglength_bodylength_m)~ log(EVOL_DIS),correlation=ou.PGLS)
summary(ou.gls_6)
plot(log(PGLS_forewinglength_bodylength_m) ~ log(EVOL_DIS))
abline(a = coef(ou.gls_6)[1], b = coef(ou.gls_6)[2], col = "green")
abline(a = coef(bm.gls_6)[1], b = coef(bm.gls_6)[2], col = "blue")
ou.gls_zero_6 <- gls(log(PGLS_forewinglength_bodylength_m)~1, correlation=ou.PGLS)
1-(as.numeric((logLik(ou.gls_6)-1)/logLik(ou.gls_zero_6)))


ou.gls_7<-gls(PGLS_DIMORPHISM_forewing_BD_ratio ~ log(EVOL_DIS),correlation=ou.PGLS)
summary(ou.gls_7)
plot(PGLS_DIMORPHISM_forewing_BD_ratio ~ log(EVOL_DIS))
abline(a = coef(ou.gls_7)[1], b = coef(ou.gls_7)[2], col = "green")
abline(a = coef(bm.gls_7)[1], b = coef(bm.gls_7)[2], col = "blue")
ou.gls_zero_7 <- gls(PGLS_DIMORPHISM_forewing_BD_ratio~1, correlation=ou.PGLS)
1-(as.numeric((logLik(ou.gls_7)-1)/logLik(ou.gls_zero_7)))

ou.gls_8<-gls(PGLS_DIMORPHISM_bodylength~ log(EVOL_DIS),correlation=ou.PGLS)
summary(ou.gls_8)
plot(PGLS_DIMORPHISM_bodylength ~ log(EVOL_DIS))
abline(a = coef(ou.gls_8)[1], b = coef(ou.gls_8)[2], col = "green")
abline(a = coef(bm.gls_8)[1], b = coef(bm.gls_8)[2], col = "blue")
ou.gls_zero_8 <- gls(PGLS_DIMORPHISM_bodylength~1, correlation=ou.PGLS)
1-(as.numeric((logLik(ou.gls_8)-1)/logLik(ou.gls_zero_8)))


ou.gls_9<-gls(PGLS_DIMORPHISM_pronlength ~ log(EVOL_DIS),correlation=ou.PGLS)
summary(ou.gls_9)
plot(PGLS_DIMORPHISM_pronlength ~ log(EVOL_DIS))
abline(a = coef(ou.gls_9)[1], b = coef(ou.gls_9)[2], col = "green")
abline(a = coef(bm.gls_9)[1], b = coef(bm.gls_9)[2], col = "blue")
ou.gls_zero_9 <- gls(PGLS_DIMORPHISM_pronlength~1, correlation=ou.PGLS)
1-(as.numeric((logLik(ou.gls_9)-1)/logLik(ou.gls_zero_9)))


ou.gls_10<-gls(PGLS_DIMORPHISM_forewinglength ~ log(EVOL_DIS),correlation=ou.PGLS)
summary(ou.gls_10)
plot(PGLS_DIMORPHISM_forewinglength ~ log(EVOL_DIS))
abline(a = coef(ou.gls_10)[1], b = coef(ou.gls_10)[2], col = "green")
abline(a = coef(bm.gls_10)[1], b = coef(bm.gls_10)[2], col = "blue")
ou.gls_zero_10 <- gls(PGLS_DIMORPHISM_forewinglength~1, correlation=ou.PGLS)
1-(as.numeric((logLik(ou.gls_10)-1)/logLik(ou.gls_zero_10)))

ou.gls_11<-gls(PGLS_DIMORPHISM_forewinglength_pronotum ~ log(EVOL_DIS),correlation=ou.PGLS)
summary(ou.gls_11)
plot(PGLS_DIMORPHISM_forewinglength_pronotum ~ log(EVOL_DIS))
abline(a = coef(ou.gls_11)[1], b = coef(ou.gls_11)[2], col = "green")
abline(a = coef(bm.gls_11)[1], b = coef(bm.gls_11)[2], col = "blue")
ou.gls_zero_11 <- gls(PGLS_DIMORPHISM_forewinglength_pronotum~1, correlation=ou.PGLS)
1-(as.numeric((logLik(ou.gls_11)-1)/logLik(ou.gls_zero_11)))


#### 4. Plots ####

library(RColorBrewer)
tree<-reorder(tree,"cladewise")
X<-data.frame(data$wings_display, data$arms_display, data$wings_colours, data$arms_colours)

X<-X[tree$tip.label,]
plotTree(tree,plot=FALSE)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
plotTree(tree,lwd=1,ylim=c(0,obj$y.lim[2]*1.05),xlim=c(0,obj$x.lim[2]*1.3),
         ftype="off")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
h<-max(obj$xx)
fsize<-0.6
for(i in 1:Ntip(tree)){
  lines(c(obj$xx[i],h),rep(obj$yy[i],2),lty="dotted")
  text(h,obj$yy[i],tree$tip.label[i],cex=fsize,pos=4,font=3,offset=0.1)
}
s<-max(fsize*strwidth(tree$tip.label))
start.x<-1.05*h+s

palettes<-c("PuOr","PuOr", "PuOr", "PuOr")
cols<-list()
for(i in 1:ncol(X)){
  text(start.x,max(obj$yy)+1,paste("trait",colnames(X)[i]),pos=4,srt=60,
       cex=0.8,offset=0)
  cols[[i]]<-setNames(sample(brewer.pal(max(3,length(levels(X[[i]]))),
                                        palettes[i]),length(levels(X[[i]]))),levels(X[[i]]))
  for(j in 1:nrow(X)){
    xy<-c(start.x,obj$yy[j])
    y<-c(xy[2]-0.5,xy[2]+0.5,xy[2]+0.5,xy[2]-0.5)
    asp<-(par()$usr[2]-par()$usr[1])/(par()$usr[4]-par()$usr[3])*
      par()$pin[2]/par()$pin[1]
    x<-c(xy[1]-0.5*asp,xy[1]-0.5*asp,xy[1]+0.5*asp,xy[1]+0.5*asp)
    polygon(x,y,col=cols[[i]][as.character(X[[i]][j])])
  }
  start.x<-start.x+2*asp
}
start.y<-max(obj$yy)
for(i in 1:ncol(X)){
  text(start.x,start.y,paste("trait",colnames(X)[i]),pos=4,cex=0.9,
       offset=0)
  add.simmap.legend(colors=cols[[i]],shape="square",prompt=FALSE,
                    x=start.x,y=start.y-2*strheight("W")*0.9,fsize=0.9)
  start.y<-start.y-1.5*0.9*strheight("W")*(length(cols[[i]])-1)-6
}
 title("Fig 3 A")

 tree<-reorder(tree,"cladewise")
 X<-data.frame(data$abdomen_colours, data$wings_eyespots, data$arms_eyespots, 
               data$display_includes_sound_abd_or_wings, data$mouth_display)
 
 X<-X[tree$tip.label,]
 plotTree(tree,plot=FALSE)
 obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
 plotTree(tree,lwd=1,ylim=c(0,obj$y.lim[2]*1.05),xlim=c(0,obj$x.lim[2]*1.3),
          ftype="off")
 obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
 h<-max(obj$xx)
 fsize<-0.6
 for(i in 1:Ntip(tree)){
   lines(c(obj$xx[i],h),rep(obj$yy[i],2),lty="dotted")
   text(h,obj$yy[i],tree$tip.label[i],cex=fsize,pos=4,font=3,offset=0.1)
 }
 s<-max(fsize*strwidth(tree$tip.label))
 start.x<-1.05*h+s
 palettes<-c("PuOr","PuOr", "PuOr", "PuOr", "PuOr")
 cols<-list()
 for(i in 1:ncol(X)){
   text(start.x,max(obj$yy)+1,paste("trait",colnames(X)[i]),pos=4,srt=60,
        cex=0.8,offset=0)
   cols[[i]]<-setNames(sample(brewer.pal(max(3,length(levels(X[[i]]))),
                                         palettes[i]),length(levels(X[[i]]))),levels(X[[i]]))
   for(j in 1:nrow(X)){
     xy<-c(start.x,obj$yy[j])
     y<-c(xy[2]-0.5,xy[2]+0.5,xy[2]+0.5,xy[2]-0.5)
     asp<-(par()$usr[2]-par()$usr[1])/(par()$usr[4]-par()$usr[3])*
       par()$pin[2]/par()$pin[1]
     x<-c(xy[1]-0.5*asp,xy[1]-0.5*asp,xy[1]+0.5*asp,xy[1]+0.5*asp)
     polygon(x,y,col=cols[[i]][as.character(X[[i]][j])])
   }
   start.x<-start.x+2*asp
 }
 start.y<-max(obj$yy)
 for(i in 1:ncol(X)){
   text(start.x,start.y,paste("trait",colnames(X)[i]),pos=4,cex=0.9,
        offset=0)
   add.simmap.legend(colors=cols[[i]],shape="square",prompt=FALSE,
                     x=start.x,y=start.y-2*strheight("W")*0.9,fsize=0.9)
   start.y<-start.y-1.5*0.9*strheight("W")*(length(cols[[i]])-1)-6
 }
 title("Fig 3B")
 
 
boxplot(evol_d~data$display_complexity)
title("Display complexity vs evolutionary distinctiveness")

boxplot(evol_d~data$Presence_display)
title("Display presence vs evolutionary distinctiveness")

boxplot(evol_d~data$World)
title("Geographic distribution vs evolutionary distinctiveness")

boxplot(evol_d~data$Region)
title("Geographic distribution vs evolutionary distinctiveness")

boxplot(evol_d~data$Resemblance)
title("Ressemblance vs evolutionary distinctiveness")




tree$tip.label
data_no_order <- read.csv("body_measurements.csv", row.names = 1)

Overlap <- treedata(tree, data_no_order)

tree<- Overlap$phy
data_pre_order <- as.data.frame(Overlap$data)

morpho <- data_pre_order[match(Overlap$phy$tip.label,rownames(data_pre_order)),] 
morpho <- data_no_order

morpho$bodylength_f_ave
X_data <- morpho[,c(1,3)]
Y_data <- morpho[,c(5,7)]
Z_data <- morpho[,c(9,11)]

X_data <- morpho[,1]
names(X_data) <- row.names(morpho)


pdf("Body length_amended.pdf")
plotTree.barplot(tree,X_data,cex=0.6, tip.labels = TRUE, method="plotTree", args.barplot=list(beside=TRUE,space=c(0,1.2)))
title("Body length_amended")
dev.off()

pdf("Pronotum length_amended.pdf")
plotTree.barplot(tree,Y_data, cex=0.6, args.barplot=list(beside=TRUE,space=c(0,1.2)))
title("Pronotum length_amended")
dev.off()

pdf("Forewing length_amended.pdf")
plotTree.barplot(tree,Z_data,cex=0.6, args.barplot=list(beside=TRUE,space=c(0,1.2)))
title("Forewing length_amended")
dev.off()

tree.test <-reorder(tree)

# Not included
# boxplot(evol_d~data$wings_display)
# boxplot(evol_d~data$arms_display)
# boxplot(evol_d~data$wings_colours)
# boxplot(evol_d~data$arms_colours)
# boxplot(evol_d~data$abdomen_colours)
# boxplot(evol_d~data$wings_eyespots)
# boxplot(evol_d~data$arms_eyespots)
# pronlength_m_ave<-as.numeric(as.character(data$pronlength_m_ave))
# names(pronlength_m_ave)<-row.names(data)
# boxplot(pronlength_m_ave~data$Presence_display)
# bodylength_m_ave<-as.numeric(as.character(data$bodylength_m_ave))
# names(bodylength_m_ave)<-row.names(data)
# boxplot(bodylengthh_m_ave~data$Presence_display)
# 
# wings_eyespots <- data$wings_eyespots
# names(wings_eyespots) <- row.names(data)
# procD.pgls(evol_d~wings_eyespots, tree, iter = 999)


Overlap_cicle <- treedata(tree, data)

primary_defence_Y <- Overlap_cicle$data[,3]
display_complexity_Y <- Overlap_cicle$data[,18]

obj <- contMap(tree,evol_d,  plot=FALSE)
obj<-setMap(obj,invert=TRUE)
plot(obj,type= 'fan', fsize=c(0.6,1),outline=FALSE,lwd=c(3,7),leg.txt="evol_distinct", ftype="i", offset= 2, open.angle=0, show.tip.label = FALSE)
plot.phylo(tree, type='i', show.tip.label = FALSE)
plot(tree, "fan", show.tip.label = FALSE, open.angle = 0)
cols<-setNames(palette()[1:length(unique(secondary_defence))],sort(unique(secondary_defence)))
ring(primary_defence_Y, tree, "ring", offset = 1)
ring(primary_defence_Y, tree, "ring", offset = 2)
ring(primary_defence_Y, tree, "ring", offset = 3)
ring(primary_defence_Y, tree, "ring", offset = 4)
ring(display_complexity_Y, tree, style = "ring", offset = 5)
ring(display_complexity_Y, tree, style = "ring", offset = 6)
ring(display_complexity_Y, tree, style = "ring", offset = 7)
ring(display_complexity_Y, tree, style = "ring", offset = 8)


# MODELS

# 1.1 males 
logreg1.1 <- glm(secondary_defence_present ~ bodylength_m_ave, data=data, family=binomial(link=logit))
summary(logreg1.1)
plot(allEffects(logreg1.1), type="response")
range(data$bodylength_m_ave)
xbodylength_m <- seq(14, 118, 0.01)
ybodylength_m <- predict(logreg1.1, list(bodylength_m_ave = xbodylength_m), type="response")
plot(data$bodylength_m_ave, data$secondary_defence_present, xlab = "body length males (average)", ylab = "defensive display") # plot data
lines(xbodylength_m, ybodylength_m) # add prediction line

pdf("bodysizem_display.pdf", onefile=TRUE)
plot(data$bodylength_m_ave, data$secondary_defence_present, xlab = "body length males (average)", ylab = "defensive display",xlim=c(0,140)) # plot data
lines(xbodylength_m, ybodylength_m) # add prediction line
dev.off()

# 1.2 females
logreg1.2 <- glm(secondary_defence_present ~ bodylength_f_ave, data=data, family=binomial(link=logit))
summary(logreg1.2)
plot(allEffects(logreg1.2), type="response")
range(data$bodylength_f_ave, na.rm=TRUE)
xbodylength_f <- seq(15, 127, 0.01)
ybodylength_f <- predict(logreg1.2, list(bodylength_f_ave = xbodylength_f), type="response")
plot(data$bodylength_f_ave, data$secondary_defence_present, xlab = "body length females (average)", ylab = "defensive display") # plot data
lines(xbodylength_f, ybodylength_f) # add prediction line

pdf("bodysizef_display.pdf", onefile=TRUE)
plot(data$bodylength_f_ave, data$secondary_defence_present, xlab = "body length females (average)", ylab = "defensive display", xlim=c(0,140)) # plot data
lines(xbodylength_f, ybodylength_f) # add prediction line
dev.off()

