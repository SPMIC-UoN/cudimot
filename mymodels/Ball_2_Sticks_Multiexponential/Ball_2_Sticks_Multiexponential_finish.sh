#!/bin/sh
#
#   Moises Hernandez-Fernandez - FMRIB Image Analysis Group
#
#   Copyright (C) 2004 University of Oxford
#
#   SHCOPYRIGHT

Usage() {
    echo ""
    echo "Usage: Ball_2_Sticks_Multiexponential_finish.sh <subject_directory>"
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
mv $directory/Param_6_samples.nii.gz $directory/f2_samples.nii.gz
mv $directory/Param_7_samples.nii.gz $directory/th2_samples.nii.gz
mv $directory/Param_8_samples.nii.gz $directory/ph2_samples.nii.gz

$FSLDIR/bin/fslmaths $directory/S0_samples.nii.gz -Tmean $directory/mean_S0
$FSLDIR/bin/fslmaths $directory/d_samples.nii.gz -Tmean $directory/mean_d
$FSLDIR/bin/fslmaths $directory/d_std_samples.nii.gz -Tmean $directory/mean_d_std
$FSLDIR/bin/fslmaths $directory/f1_samples.nii.gz -Tmean $directory/mean_f1
$FSLDIR/bin/fslmaths $directory/f2_samples.nii.gz -Tmean $directory/mean_f2

$FSLDIR/bin/make_dyadic_vectors $directory/th1_samples $directory/ph1_samples $directory/nodif_brain_mask.nii.gz dyads1
$FSLDIR/bin/make_dyadic_vectors $directory/th2_samples $directory/ph2_samples $directory/nodif_brain_mask.nii.gz dyads2
${FSLDIR}/bin/maskdyads $directory/dyads2 $directory/mean_f2

${FSLDIR}/bin/fslmaths $directory/mean_f2 -div $directory/mean_f1 $directory/mean_f2_f1
${FSLDIR}/bin/fslmaths $directory/dyads2_thr0.05 -mul $directory/mean_f2_f1 $directory/dyads2_thr0.05_modf2
${FSLDIR}/bin/imrm $directory/mean_f2_f1


${FSLDIR}/bin/fslmaths $directory/mean_f1 -mul 0 $directory/mean_fsum
${FSLDIR}/bin/fslmaths $directory/mean_fsum -add $directory/mean_f1 $directory/mean_fsum
${FSLDIR}/bin/fslmaths $directory/mean_fsum -add $directory/mean_f2 $directory/mean_fsum



