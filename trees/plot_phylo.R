### Script to plot out some trees from the Ameiva data

# On MedBow:
#  module load gcc/14.2.0 r/4.4.0 gdal/3.7.3 rstudio/2024.04.2


# load up some packages
library(ape)
library(phytools)
library(LEA)
library(vcfR)
library(plotrix)



# Load up the paths we need
iq_file <- "/project/inbreh/ameiva/iqtree_out/ameiva_dn_c92_nolowcov.tre"
coords <- read.csv("/project/inbreh/ameiva/metadata/coords_spec_info.csv", header=TRUE, row.names=NULL) # read in geographic coordinates for the samples
outdir <- "/project/inbreh/ameiva/iqtree_out"
# set the path to the snmf project for exsul only
snmf_proj_path <- "/project/inbreh/ameiva/ipyrad_out/ameiva_dn_c92_exsulONLY_outfiles/ameiva_dn_c92_exsulONLYLINKED.snmfProject"
path_vcf <- "/project/inbreh/ameiva/ipyrad_out/ameiva_dn_c92_exsulONLY_outfiles/ameiva_dn_c92_exsulONLY.vcf"


#plot IQTree
iq_tree <- read.tree(iq_file)
# root the tree on the outgroup
rooted_tree <- root(iq_tree, outgroup = "38_Cuba", resolve.root = TRUE)

plot(rooted_tree, show.tip.label = FALSE)

edgelabels() # add labels on the edges (branches) to get their indexes for the next step.



rooted_tree$edge.length[1]

# take part of the length of edge 144 and assign it to the length of edge 1, which is zero right now
# Previously did this by halving it, comments will still reflect that

# Get the current length of edge 144
edge_144_length <- rooted_tree$edge.length[144]

# Halve the length of edge 144
new_edge_144_length <- edge_144_length * 0.8

# Get the current length of edge 1
edge_1_length <- rooted_tree$edge.length[1]

# Assign the new halved length to edge 1
rooted_tree$edge.length[144] <- new_edge_144_length

# Add the remaining length to edge 144
rooted_tree$edge.length[1] <- edge_1_length + (edge_144_length - new_edge_144_length)


plot(rooted_tree, show.tip.label = FALSE)
# this looks good


pdf(file = paste0(outdir, "/IQtree_support.pdf"), width = 14, height = 20)
plot(rooted_tree)
nodelabels(text = rooted_tree$node.label, adj = c(1.05, 1.5), cex = 0.7, frame = "none")
add.scale.bar()
dev.off()



# Now that we made that look nice, let's actually just drop that tip out entirely to make it plot on the map better:
ingroup_tree <- drop.tip(rooted_tree, tip = "38_Cuba")



# Plot tree to map:

# some set up 

# Get the coordinates for just the individuals we want to plot and put in the same order as the tip labels - we read in the coordinates up top when we read in the tree file
ind_names <- ingroup_tree$tip.label # get individual names from the tree
ind_names[which(!ind_names %in% coords[,"gen_dat_num"])] # check if there are any individuals not in the coordinates


# match up the coordinates to the order of the individuals in the tree
match_coords <- match(ind_names, coords[,"gen_dat_num"]) # get index of the coords that are in ind_names
tree_coords <- coords[match_coords,] # get just the coordinates for individuals in the tree in the order they show up in the tree
rownames(tree_coords) <- tree_coords[, "gen_dat_num"] # for phylo.to.map, we need the individual names from the tree to be the row names of the coordinates, not a column
tree_coords <- tree_coords[, c("Latitude", "Longitude")] # make the tree_coords object only 2 columns, Latitude and Longitude



# Plot the tree to the map
pdf(file = paste0(outdir, "/IQtree_map.pdf"), width = 14, height = 12)
phylo.to.map(ingroup_tree, tree_coords, fsize = 0.15)
dev.off()


# Let's prune this to only exsul for easier viewing
plot(ingroup_tree, cex = 0.5)
nodelabels()
# Get the MRCA of taxa to be dropped
mrca_node <- getMRCA(ingroup_tree, c("PRBVI2H09_Guania", "PRBVI2H06_Protestant"))
# Get all descendant tips of that node:
to_remove <- extract.clade(ingroup_tree, mrca_node)$tip.label
# Drop those taxa from the tree
exsul_tree <- drop.tip(ingroup_tree, to_remove)
plot(exsul_tree)


# plot just this to map:
# Get the coordinates for just the individuals we want to plot and put in the same order as the tip labels - we read in the coordinates up top when we read in the tree file
ind_names <- exsul_tree$tip.label # get individual names from the tree
ind_names[which(!ind_names %in% coords[,"gen_dat_num"])] # check if there are any individuals not in the coordinates


# match up the coordinates to the order of the individuals in the tree
match_coords <- match(ind_names, coords[,"gen_dat_num"]) # get index of the coords that are in ind_names
tree_coords <- coords[match_coords,] # get just the coordinates for individuals in the tree in the order they show up in the tree
rownames(tree_coords) <- tree_coords[, "gen_dat_num"] # for phylo.to.map, we need the individual names from the tree to be the row names of the coordinates, not a column
tree_coords <- tree_coords[, c("Latitude", "Longitude")] # make the tree_coords object only 2 columns, Latitude and Longitude



pdf(file = paste0(outdir, "/IQtree_Exsul_map.pdf"), width = 14, height = 12)
phylo.to.map(exsul_tree, tree_coords, fsize = 0.5, lty = "solid", colors = "black")
dev.off()




#####################################################################################################################################################################
#####################################################################################################################################################################
############        Jittered coorindates map with pies
#####################################################################################################################################################################
#####################################################################################################################################################################




## Add SNMF pies onto the exsul:
# laod the snmf project
obj.at <- load.snmfProject(snmf_proj_path)

# make a list of colors:
colors_6<-c("V1" = "red", "V2" = "blue", "V3" = "white", "V4" = "purple", "V5" = "pink", "V6" = "yellow")



# read in the vcf file to get the order of the individual names from ipyrad to make sure we can match these up to the SNMF qmatrix correctly
gendata_all <- read.vcfR(path_vcf) # read in all of the genetic data from the vcf file
gendata <- vcfR2genind(gendata_all) # convert to genind format
ind_names <- rownames(gendata@tab)





# Create a small function to nudge overlapping coordinates
nudge_coords <- function(data, nudge_by = 0.01) {
  duplicated_coords <- duplicated(data[, c("Longitude", "Latitude")]) |
    duplicated(data[, c("Longitude", "Latitude")], fromLast = TRUE)
  data[duplicated_coords, "Longitude"] <- data[duplicated_coords, "Longitude"] + 
    runif(sum(duplicated_coords), -nudge_by, nudge_by)
  data[duplicated_coords, "Latitude"] <- data[duplicated_coords, "Latitude"] + 
    runif(sum(duplicated_coords), -nudge_by, nudge_by)
  return(data)
}

jittered_coords <- nudge_coords(coords, nudge_by = 0.03) # jitter the coords a bit


# Match jittered coordinates to tree tips
match_coords <- match(exsul_tree$tip.label, jittered_coords[,"gen_dat_num"])
tree_coords_jittered <- jittered_coords[match_coords, c("Latitude", "Longitude")]
rownames(tree_coords_jittered) <- jittered_coords[match_coords, "gen_dat_num"]

# Match jittered coordinates for SNMF pies
match_coords_snmf <- match(ind_names, jittered_coords[,"gen_dat_num"])
snmf_coords_jittered <- jittered_coords[match_coords_snmf,]





# Set up values of k to plot
k_vals <- 2:6
  
pdf(file = paste0(outdir, "/IQ_piemap_exsul_jittered.pdf"), width = 14, height = 12)


# make pdf of cross-entropy plot
plot(obj.at, col = "lightblue", cex = 1.2, pch = 19)


#Plot these in a loop within a pdf
for(k in k_vals){

  # set a value of k, the number of clusters, we want to use 2 here

  
  # confirm cross entropy values for K are consist. across runs
  ce <- cross.entropy(obj.at, K = k) 
  best.run <- which.min(ce) # find the run with the lowest cross validation error
  
  ## Get the snmf Q matrix from the best run at the best k
  qmatrix <- Q(obj.at, K = k, run = best.run)
  admix <- as.data.frame(qmatrix)
 
  # Create a combined data frame with jittered coordinates and admixture proportions
  for_pies_jittered <- cbind(admix, snmf_coords_jittered)
  

  # Get the right number of colors
  pie_colors <- colors_6[1:ncol(admix)]
  
  # Plot tree on map with jittered coordinates
  phylo.to.map(exsul_tree, tree_coords_jittered, lty = "solid", colors = "black", fsize = 0.5)
  
  # Add jittered pies
  for (i in 1:nrow(for_pies_jittered)) {
    x <- for_pies_jittered$Longitude[i]
    y <- for_pies_jittered$Latitude[i]
    proportions <- as.numeric(for_pies_jittered[i, 1:k])
    # pie_colors <- c("red", "blue")  
    
    # Plot floating pie at jittered position
    floating.pie(xpos = x, ypos = y, radius = 0.01, x = proportions, col = pie_colors)
  }
  title(paste0("exsul_k_", k), line = -1.2)

  
}  
dev.off()

