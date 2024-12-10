#### Script to do population assignment on RADseq data for Carribean Ameiva


# On MedicineBow, need to load gdal module - required for mapping packages
# module load gcc/14.2.0 r/4.4.0 gdal/3.7.3 rstudio/2024.04.2



### load up relevant packages
library(adegenet)
library(LEA)
library(plotrix)
library(mapdata)
library(rworldmap)
library(ggplot2)
library(scatterpie)
library(dplyr)

## Set up an object to contain the path to the main directory with the data and then set that as the working directory
main_dir<-"/project/inbreh/ameiva/ipyrad_out/ameiva_dn_c92_no_outgroup_outfiles"
setwd(main_dir)

 # path to coordinates file
coords_file <- "/project/inbreh/ameiva/metadata/coords_spec_info.csv"

## Set up an output directory
sNMF_out_dir<-"/project/inbreh/ameiva/popstr_out"  # specify a full path to the directory
if(!dir.exists(sNMF_out_dir)){ # check if the directory  exists and then only create it if it does not
  dir.create(sNMF_out_dir)
}


## Specify all of the assemblies that we want to run sNMF on - 
##    looped from older versions of this script for other assemblies that 
##    used multiple assemblies
all_assemblies<-c(
  # "ameiva_dn_c92_nolowcov"
  "ameiva_dn_c92_no_outgroup"
)


#### Some overall setup for mapping and plotting

## Getout map data for area of interest
world <- map_data("worldHires") # Mexico data

# Filter for Caribbean region (Cuba, Virgin Islands, etc.)
to_map <- world %>%
  filter(long > -68 & long < -64, lat > 17 & lat < 20)  # Approximate coordinates for the Caribbean


# make a list of colors:
colors_6 <- c("V1" = "red", "V2" = "blue", "V3" = "white", "V4" = "purple", "V5" = "pink", "V6" = "yellow")

## Read in coordinates for plotting farther down
setwd(main_dir)
coords <- read.csv(coords_file, header=TRUE, row.names=NULL) # coordinates of everything I sequenced and many I didn't



######################################################################################################################
## Loop to run sNMF over all assemblies
######################################################################################################################

for(species in all_assemblies){   ### if we want to loop over all assemblies, this line and line starting "all_assemblies<-c" should be uncommented, as well as final "}" if doing a single assembly, comment these lines instead
  ###########################################################
  ## Set up paths to input files
  ###########################################################
  setwd(main_dir)
  path_ugeno<-paste0(main_dir,"/", species,".ugeno")
  path_ustr<-paste0(main_dir,"/", species,".ustr")
  
  # snmf requires the geno file to have the extension .geno - the geno file of unlinked snps has ugeno
  #   as above, copy the geno and make one with the extension .u.geno
  path_geno<-gsub(".ugeno", ".u.geno", path_ugeno)  # Use a regular expression substitution to generate the new file name
  file.copy(path_ugeno, path_geno) # do the copying with the new name
  
  
  # Run sNMF using 1 to 10 ancestral populations and evaluate the fit of different k values to the data using cross entropy criterion
  # before running snmf, check if it's already been run
  if(dir.exists(gsub("geno", "snmf", basename(path_geno)))){
    obj.at<-load.snmfProject(gsub("geno", "snmfProject", basename(path_geno))) # if it has, just load up the results
  }else{ # otherwise, run sNMF
    obj.at <- snmf(input.file = path_geno,  # input file is the .geno format file. We set up the path to this above
                   K = 1:10, # we will test for k=1 through 10
                   ploidy = 2, 
                   entropy = T, # use the cross entropy criterion for assessing the best k value
                   repetitions = 10, # Run 10 independent replicate analyses
                   CPU = 1, 
                   project = "new", tolerance = 0.00001, iterations = 500)
  }
  
  setwd(sNMF_out_dir)
  
  # make pdf of cross-entropy plot
  plot(obj.at, col = "lightblue", cex = 1.2, pch = 19)
  
  # look at outstats
  outstats <- summary(obj.at)
  outstats # take a look
  
  # Plot k=2 through k=6 for all
  k_plot<-2:9
  
  
  ## This code block reads in the ustr file to get individual names in the order they show up
  ##    in data files, since geno files don't have ind names in them - this is not a great way
  ##    to do it, and is a holdover from when I used the ustr for other stuff that I've removed from this script,
  ##    but it works, so it stays - if I was building this ground-up again, I'd do this differently
  ##
  geno_txt<-readLines(path_ugeno)
  nums_snps<-length(geno_txt)
  num_ind<-length(strsplit(geno_txt[[1]], "")[[1]])
  ## quirk of read.structure function is that it requires the strucure file to have the file extension “.stru” - do some copying to make a new file with this extension
  path_stru<-gsub(".ustr", ".stru", path_ustr)  # Use a regular expression substitution to generate the new file name
  file.copy(path_ustr, path_stru) # make a copy of the file with the new name
  # Now we can read in this file
  ustr<-read.structure(path_stru, n.ind=num_ind, n.loc=nums_snps, onerowperind = FALSE, col.lab=1, col.pop=0, NA.char="-9", pop=NULL, ask=FALSE, quiet=FALSE)
  ind_names<-rownames(ustr@tab) ## get the individual names in the order that they show up in the various files - this is important farther down for getting coordinates into the right order for plotting
  
  
  
  ### For  down below, get the geographic coordinates sorted out
  ## make sure there aren't any individuals that don't have coordinates
  ind_names[which(!ind_names %in% coords[,"gen_dat_num"])]
  # match up the coordinates to the order of the individuals from snmf
  match_coords<-match(ind_names, coords[,"gen_dat_num"])
  snmf_coords<-coords[match_coords,]
  
  
  pdf(file=paste0(species,"_SNMF_plots", ".pdf"), width=6, height=5)
  # put in the cross-entropy plot at the start
  plot(obj.at, col = "lightblue", cex = 1.2, pch = 19)
  #### use a loop to plot various different k values 
  for(i in k_plot){
    # confirm cross entropy values for K are consist. across runs
    ce <- cross.entropy(obj.at, K = i) 
    ce # pretty similar
    best.run <- which.min(ce) # find the run with the lowest cross validation error
    
    ## Get the snmf Q matrix from the best run at the best k
    qmatrix <- Q(obj.at, K = i, run = best.run)
    admix<-as.data.frame(qmatrix)
    
    # get the coordinate and admix data into a single dataframe
    for_pies <- cbind(snmf_coords, admix)
    
    
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
    
    for_pies <- nudge_coords(for_pies, nudge_by = 0.03) # jitter the coords a bit
    
    
    
    # Get the right number of colors
    colors <- colors_6[1:ncol(admix)]
    
    
    # plot it out
    snmf_plot <- ggplot(to_map, aes(long, lat, group = group)) + # map out the US & Mexico
      geom_polygon(data = to_map, fill = "grey90", color = "black", size = 0.2) + # make them polygons
      geom_scatterpie(data = for_pies, aes(x=Longitude, y=Latitude, group = gen_dat_num, r = 0.008), cols = grep("^V", colnames(for_pies), value = TRUE), size = 0.1) + # plot the pies - use grep to get the column names that start with V, these are the admix proportions
      scale_fill_manual(values = colors) +
      guides(fill="none") + # get rid of the legend for admixture
      theme_minimal() +
      labs(title=paste0(species,"_SNMF_K",i), x ="Longitude", y = "Longitude") +
      coord_map("moll") # Mollweide projection

    
    # Save current margins to restore later
    old_par <- par()
    
    # Increase the bottom margin
    par(mar = c(8, 4, 4, 2))  # c(bottom, left, top, right)
    
    
    # Seems might have to plot barchart here:
    bp <- barchart(obj.at, K = i, run = best.run,
             border = "black", space = 0,
             col = colors,
             xlab = NULL,
             ylab = "Ancestry proportions",
             main = paste0(species," SNMF K ",i))
    axis(1, at = 1:length(bp$order),
         labels = ind_names[bp$order], las=2,
         cex.axis = 0.3)
    
    par(old_par)
    
    
    print(snmf_plot)
    
  }
  dev.off()
}
