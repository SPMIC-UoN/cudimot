#!/bin/sh
#
#   Moises Hernandez-Fernandez - FMRIB Image Analysis Group
#
#   Copyright (C) 2004 University of Oxford
#
#   SHCOPYRIGHT

Usage() {
    echo ""
    echo "Usage: Ball_1_Sticks_Multiexponential_finish.sh <subject_directory>"
    echo ""
    echo "expects to find all the estimatedParameters and nodif_brain_mask in subject directory"
    echo ""
    exit 1
}

[ "$1" = "" ] && Usage

directory=$1
cd ${directory}

mv $directory/Param_0_samples.nii.gz $directory/S0_samples.nii.gz
mv $directory/Param_1_samples.nii.gz $directory/d_samples.nii.gz
mv $directory/Param_2_samples.nii.gz $directory/d_std_samples.nii.gz
mv $directory/Param_3_samples.nii.gz $directory/f1_samples.nii.gz
mv $directory/Param_4_samples.nii.gz $directory/th1_samples.nii.gz
mv $directory/Param_5_samples.nii.gz $directory/ph1_samples.nii.gz

$FSLDIR/bin/fslmaths $directory/S0_samples.nii.gz -Tmean $directory/mean_S0
$FSLDIR/bin/fslmaths $directory/d_samples.nii.gz -Tmean $directory/mean_d
$FSLDIR/bin/fslmaths $directory/f1_samples.nii.gz -Tmean $directory/mean_f1

$FSLDIR/bin/make_dyadic_vectors $directory/th1_samples $directory/ph1_samples $directory/nodif_brain_mask.nii.gz dyads1



