#### Script to plot output from Admixture for Carribean Ameiva


# On MedicineBow, need to load gdal module - required for mapping packages
# module load gcc/14.2.0 r/4.4.0 gdal/3.7.3 rstudio/2024.04.2



### load up relevant packages
library(plotrix)
library(mapdata)
library(rworldmap)
library(ggplot2)
library(scatterpie)
library(dplyr)
library(vcfR)
library(stringr)
library(tidyr)


## Set up an object to contain the path to the main directory with the data and then set that as the working directory
main_dir<-"/project/inbreh/ameiva/admix_out/"
setwd(main_dir)

# path to coordinates file
coords_file <- "/project/inbreh/ameiva/metadata/coords_spec_info.csv"

# path to ipyrad output to read in genetic data for individual names
ipyrad_out_dir <- "/project/inbreh/ameiva/ipyrad_out/"

## Specify all of the assemblies that we want to run sNMF on - 
##    looped from older versions of this script for other assemblies that 
##    used multiple assemblies
all_assemblies<-c(
  # "ameiva_dn_c92_no_outgroup",
  # "ameiva_dn_c92_exsul",
  "ameiva_dn_c92_exsulONLY"
)



#### Some overall setup for mapping and plotting

## Getout map data for area of interest
world <- map_data("worldHires") # Mexico data

# Filter for Caribbean region (Cuba, Virgin Islands, etc.)
to_map <- world %>%
  filter(long > -68 & long < -64, lat > 17 & lat < 20)  # Approximate coordinates for the Caribbean


# make a list of colors:
colors_6 <- c("V1" = "red", "V2" = "blue", "V3" = "white", "V4" = "purple", "V5" = "pink", 
              "V6" = "yellow", "V7" = "black", "V8" = "green", "V9" = "orange", V10 = "cyan")


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
  setwd(paste0(main_dir, species)) # move into directory for assembly
  
  # read in genetic data from ipyrad out to get individual names in order
  path_vcf<-paste0(ipyrad_out_dir, species, "_outfiles/", species,".vcf")
  gendata_all<-read.vcfR(path_vcf) # read in all of the genetic data
  gendata<-vcfR2genlight(gendata_all) # make the genetic data a biallelic matrix of alleles in genlight format
  ind_names<-gendata@ind.names ## get the individual names in the order that they show up in the various files - this is important farther down for getting coordinates into the right order for plotting

  
  # read in the cross validation summary
  rawcv <- readLines("CV_summ.txt")
  
  # Extract K values and corresponding error values
  k_values <- as.numeric(str_extract(rawcv, "(?<=K=)\\d+"))
  errors <- as.numeric(str_extract(rawcv, "(?<=: )\\d+\\.\\d+"))
  
  cv <- data.frame(K = k_values, Error = errors)

  # plot out cross validation
  plot(cv, col = "lightblue", cex = 1.2, pch = 19)
  
  
  # list out the Q files
  qfiles <- normalizePath(list.files(pattern = ".Q$", full.names = TRUE))
  
  

  # Plot k 2 through 10
  k_plot<-2:10
  

  
  ### For  down below, get the geographic coordinates sorted out
  ## make sure there aren't any individuals that don't have coordinates
  ind_names[which(!ind_names %in% coords[,"gen_dat_num"])]
  # match up the coordinates to the order of the individuals from admixture
  match_coords<-match(ind_names, coords[,"gen_dat_num"])
  admix_coords<-coords[match_coords,]
  
  setwd(main_dir)
  pdf(file=paste0(species,"_Admixture_plots", ".pdf"), width=6, height=5)
  # put in the cross-entropy plot at the start
  plot(cv, col = "lightblue", cex = 1.2, pch = 19)
  #### use a loop to plot various different k values 
  for(i in k_plot){
    
    ## Get the Q matrix for this value of k
    qfile <- grep(paste0("\\.", i, "\\.Q$"), qfiles, value = TRUE)
    qmatrix <- as.data.frame(read.table(qfile))
        
    # get the coordinate and admix data into a single dataframe
    for_pies <- cbind(admix_coords, qmatrix)
    
    
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
    colors <- colors_6[1:ncol(qmatrix)]
    
    
    # plot it out
    admix_plot <- ggplot(to_map, aes(long, lat, group = group)) + # map out the US & Mexico
      geom_polygon(data = to_map, fill = "grey90", color = "black", size = 0.2) + # make them polygons
      geom_scatterpie(data = for_pies, aes(x=Longitude, y=Latitude, group = gen_dat_num, r = 0.008), cols = grep("^V", colnames(for_pies), value = TRUE), size = 0.1) + # plot the pies - use grep to get the column names that start with V, these are the admix proportions
      scale_fill_manual(values = colors) +
      guides(fill="none") + # get rid of the legend for admixture
      theme_minimal() +
      labs(title=paste0(species,"_Admixture_K",i), x ="Longitude", y = "Longitude") +
      coord_map("moll") # Mollweide projection

    
    # Save current margins to restore later
    old_par <- par()
    
    # Increase the bottom margin
    par(mar = c(8, 4, 4, 2))  # c(bottom, left, top, right)
    

    # Seems might have to plot barchart here:
    
    # Convert data to long format for ggplot
    for_pies$individual <- ind_names
    
    long_data <- for_pies %>%
      pivot_longer(cols = starts_with("V"), 
                   names_to = "Cluster", 
                   values_to = "Proportion")
    
    # Plot the barchart
    barchart_plot <- ggplot(long_data, aes(x = individual, 
                                           y = Proportion, 
                                           fill = Cluster)) +
      geom_bar(stat = "identity", position = "stack", color = "black") +
      scale_fill_manual(values = colors) +
      labs(title = paste0(species, " Admixture K ", i),
           x = "Individuals", 
           y = "Ancestry proportions") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
            legend.position = "none") # Remove legend
    
    
    # Display the plot
    print(barchart_plot)
    
    
    par(old_par)
    
    print(admix_plot)
    
  }
  dev.off()
}
