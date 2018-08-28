# Read coordinates
# Load 'readxl' package to read .xlsx file
if(!is.element("readxl",rownames(installed.packages()))){
  install.packages("readxl")  
}
library(readxl)

# ROI coordinates can be found here:
# https://www.jonathanpower.net/2011-neuron-bigbrain.html
xls <- read_xlsx("Neuron_consensus_264.xlsx",skip=1)

# Round MNI coordinates
xls[,7:9] <- round(xls[,7:9])

# Load 'fslr' package
if(!is.element("fslr",rownames(installed.packages()))){
  install.packages("fslr")  
}
library(fslr)

# Read T1w ICBM-152 1mm
t1 <- file.path(fsldir(),"data/standard/MNI152_T1_1mm_brain.nii.gz")
t1_nii <- readnii(file.path(fsldir(),"data/standard/MNI152_T1_1mm_brain.nii.gz"))
# Get XYZ volume coordinates from MNI coordinates
xyz_coor <- cbind((xls[,7]-t1_nii@srow_x[4])*t1_nii@srow_x[1],
                   (xls[,8]-t1_nii@srow_y[4])*t1_nii@srow_y[2],
                   (xls[,9]-t1_nii@srow_z[4])*t1_nii@srow_z[3])

# Create spheres with 'fslmaths' with 5mm radius
roi_n <- nrow(xyz_coor)
for(ii in 1:roi_n){
  cat(paste("\nCreating sphere",ii,"\n"))
  # Create and empty NIfTI object with T1 standard dimensions
  coor_nii <- fslmaths(t1,opts="-mul 0")
  # Note that we have to add one voxel to the coordinates
  # cause 'oro.nifti' axis starts in one instead of zero (like FSL axis)
  xyz <- unlist(xyz_coor[ii,])+1
  coor_nii[xyz[1],xyz[2],xyz[3]] <- 1
  # Create sphere (5mm radious)
  fslmaths(coor_nii,outfile = paste0("ROI_",ii),opts="-kernel sphere 5 -fmean -bin")
}

# Merge all volumes
rois <- list.files(pattern = "ROI_")
# Lenght of this objects should be 264
# Merge list
fslmerge(rois,"t","power264_4D")

# Verify if there is any overlap between spheres
sum_nii <- fslmaths("power264_4D",opts = paste("-Tmean -mul",roi_n))
table(sum_nii)
# As we can see there is no overlap between ROIs
# (at 5mm radious, even consididering spheres were grided)

# So.. let's create the labelled atlas
# Clear NIfTI objects from workspace
rm(avg_nii,t1_nii)
# First, create empty volume
atlas_nii <- fslmaths(t1,opts="-mul 0")
# Add every sphere
for(ii in 1:roi_n){
  cat(paste("\nAdding sphere",ii,"\n"))
  # Create and empty NIfTI object with T1 standard dimensions
  roi_nii <- readnii(paste0("ROI_",ii))
  # Add labeling
  roi_nii <- roi_nii*ii
  # Add sphere
  atlas_nii <- atlas_nii+roi_nii
}
# Write results into a NIfTI
writeNIfTI(atlas_nii,"power264_1mm")

# Remove intermediate files
file.remove(rois)

# Now, the atlas can be resize to any other dimension with 'flirt_apply'
# For example to 2mm
t1_2mm <- file.path(fsldir(),"data/standard/MNI152_T1_2mm_brain.nii.gz")
# Flirt Apply
flirt_apply(infile = "power264_1mm", reffile = t1_2mm,
            initmat = paste0(fsldir(),"/etc/flirtsch/ident.mat"),
            outfile = "power264_2mm", opts = "-interp nearestneighbour")
