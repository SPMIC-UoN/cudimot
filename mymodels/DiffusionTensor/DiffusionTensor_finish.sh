#!/bin/sh
#
#   Moises Hernandez-Fernandez - FMRIB Image Analysis Group
#
#   Copyright (C) 2004 University of Oxford
#
#   SHCOPYRIGHT

# Script to generate FA, MD, V1, L1, etc...
Usage() {
    echo ""
    echo "Usage: DiffusionTensor_finish.sh <subject_directory>"
    echo ""
    exit 1
}

[ "$1" = "" ] && Usage

directory=$1
cd ${directory}

mv $directory/Param_0_samples.nii.gz $directory/S0_samples.nii.gz
$FSLDIR/bin/fslmaths $directory/S0_samples.nii.gz -Tmean $directory/mean_S0


$FSLDIR/bin/fslmaths $directory/Param_1_samples.nii.gz -Tmean $directory/Param_1_samples.nii.gz
$FSLDIR/bin/fslmaths $directory/Param_2_samples.nii.gz -Tmean $directory/Param_2_samples.nii.gz
$FSLDIR/bin/fslmaths $directory/Param_3_samples.nii.gz -Tmean $directory/Param_3_samples.nii.gz
$FSLDIR/bin/fslmaths $directory/Param_4_samples.nii.gz -Tmean $directory/Param_4_samples.nii.gz
$FSLDIR/bin/fslmaths $directory/Param_5_samples.nii.gz -Tmean $directory/Param_5_samples.nii.gz
$FSLDIR/bin/fslmaths $directory/Param_6_samples.nii.gz -Tmean $directory/Param_6_samples.nii.gz

$FSLDIR/bin/fslmerge -t $directory/Tensor $directory/Param_1_samples.nii.gz $directory/Param_2_samples.nii.gz $directory/Param_3_samples.nii.gz $directory/Param_4_samples.nii.gz $directory/Param_5_samples.nii.gz $directory/Param_6_samples.nii.gz
$FSLDIR/bin/fslmaths $directory/Tensor -tensor_decomp $directory/DTI
