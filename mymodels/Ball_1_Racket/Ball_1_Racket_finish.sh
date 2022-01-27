#!/bin/sh
#
#   Moises Hernandez-Fernandez - FMRIB Image Analysis Group
#
#   Copyright (C) 2004 University of Oxford
#
#   SHCOPYRIGHT

bindir=${CUDIMOT}/bin

Usage() {
    echo ""
    echo "Usage: Ball_1_Racket_finish.sh <subject_directory>"
    echo ""
    echo "expects to find all the estimatedParameters and nodif_brain_mask in subject directory"
    echo ""
    exit 1
}

[ "$1" = "" ] && Usage

directory=$1

mv $directory/Param_0_samples.nii.gz $directory/d_samples.nii.gz
mv $directory/Param_1_samples.nii.gz $directory/faniso_samples.nii.gz
mv $directory/Param_2_samples.nii.gz $directory/k1_samples.nii.gz
mv $directory/Param_3_samples.nii.gz $directory/k2_samples.nii.gz
mv $directory/Param_4_samples.nii.gz $directory/th_samples.nii.gz
mv $directory/Param_5_samples.nii.gz $directory/ph_samples.nii.gz
mv $directory/Param_6_samples.nii.gz $directory/psi_samples.nii.gz


Two_div_pi=0.636619772367581

$FSLDIR/bin/fslmaths $directory/d_samples.nii.gz -Tmean $directory/mean_d

$FSLDIR/bin/fslmaths $directory/faniso_samples.nii.gz -Tmean $directory/mean_faniso

$FSLDIR/bin/fslmaths $directory/k1_samples.nii.gz -Tmean $directory/mean_k1
$FSLDIR/bin/fslmaths $directory/k2_samples.nii.gz -Tmean $directory/mean_k2

$FSLDIR/bin/make_dyadic_vectors $directory/th_samples $directory/ph_samples $directory/nodif_brain_mask.nii.gz  $directory/dyads1

${bindir}/getFanningOrientation $directory/th_samples $directory/ph_samples $directory/psi_samples $directory/nodif_brain_mask.nii.gz  $directory/fanningDir

${FSLDIR}/bin/fslmaths $directory/mean_k2 -recip -atan -mul $Two_div_pi $directory/ODIp
${FSLDIR}/bin/fslmaths $directory/mean_k1 -recip -atan -mul $Two_div_pi $directory/ODIs
${FSLDIR}/bin/fslmaths $directory/mean_k1 -recip $directory/ODItot_temp1
${FSLDIR}/bin/fslmaths $directory/mean_k2 -recip $directory/ODItot_temp2
${FSLDIR}/bin/fslmaths $directory/ODItot_temp1 -mul $directory/ODItot_temp2 $directory/ODItot_temp1
${FSLDIR}/bin/fslmaths $directory/ODItot_temp1 -sqrt -abs -atan -mul $Two_div_pi $directory/ODItot
rm $directory/ODItot_temp1.nii.gz
rm $directory/ODItot_temp2.nii.gz

#Dispersion Anisotropy Index
${FSLDIR}/bin/fslmaths $directory/mean_k1 -sub $directory/mean_k2 -div $directory/mean_k2 -atan -mul $Two_div_pi $directory/DA
