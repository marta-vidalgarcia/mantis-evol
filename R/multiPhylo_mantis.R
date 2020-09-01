
#### 0. Load libraries and data ####
library(ape)
library(phytools)
library(geiger)
library(ggplot2)

phyloTime <- read.nexus("BeastTree") # Load a ultrametric tree
phyloTimeLadderized <- (ladderize(phyloTime))  # Ladderization
plot(phyloTimeLadderized, cex=0.5, no.margin = T) #Plot ladderized and rescaled tree
add.scale.bar() # Add a simple scale bar indicating the scale for the branches in your tree

write.tree(phyloTimeLadderized, "phyloTimeLadderized.nwk")
phyloTimeLadderized <- read.tree("phyloTimeLadderized.nwk")

tree<-phyloTimeLadderized

ntrees <- read.nexus("BeastTrees")

# species<-tmp[,1]
# keep.tip<-function(tree,tip) drop.tip(tree,setdiff(tree$tip.label,tip))
# ntrees<-lapply(ntrees,keep.tip,tip=species)

"multiPhylo"-> class(ntrees)

random.trees <- sample(ntrees,size=1000) # randomly re-sample 1,000 trees

ntrees <- random.trees
plotTree(ntrees[[2]])

data <- read.csv("Final_characters_CSV_amended.csv", header=T, row.names = 1)

Overlap <- vector("list", length(ntrees))

#### 1. Calculate lambda across trees and plot distribution ####

lambda <- matrix(data = NA, nrow = length(ntrees), ncol = 2)
for (i in 1:length(ntrees)){
  Overlap[[i]] <- treedata(ntrees[[i]], data[,19:20])
  lambda[i,2] <- fitDiscrete(Overlap[[i]]$phy, Overlap[[i]]$data[,2], transform="lambda")$opt$lambda # Presence display
  lambda[i,1] <- fitDiscrete(Overlap[[i]]$phy, Overlap[[i]]$data[,1], transform="lambda")$opt$lambda # Display complexity
}

lambda <- as.data.frame(lambda)
colnames(lambda) <- c("Display_complexity", "Presence_display")
write.csv(x = lambda, "lambda.csv")

# These figures have been included as a panel in the ESM_figures
# Change y axis to count instead of density
pdf("Phylogenetic_signal_1000trees_Presence_display.pdf", width=15, height=15, useDingbats = FALSE)
ggplot(lambda, aes(x = Presence_display)) + geom_density(aes(y = ..count..), fill = "lightgray") +
  geom_vline(aes(xintercept = mean(Presence_display)), 
             linetype = "dashed", size = 0.6,
             color = "blue")
dev.off()

# Change y axis to count instead of density
pdf("Phylogenetic_signal_1000trees_Display_complexity.pdf", width=15, height=15, useDingbats = FALSE)
ggplot(lambda, aes(x = Display_complexity)) + geom_density(aes(y = ..count..), fill = "lightgray") +
  geom_vline(aes(xintercept = mean(Display_complexity)), 
             linetype = "dashed", size = 0.6,
             color = "blue")
dev.off()

