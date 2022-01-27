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
    echo "Usage: NODDI_Bingham_finish.sh <output_directory>"
    echo ""
    echo "expects to find all the estimatedParameters and nodif_brain_mask in subject directory"
    echo ""
    exit 1
}

[ "$1" = "" ] && Usage

directory=$1
cd ${directory}

mv $directory/Param_0_samples.nii.gz $directory/fiso_samples.nii.gz
mv $directory/Param_1_samples.nii.gz $directory/fintra_samples.nii.gz
mv $directory/Param_2_samples.nii.gz $directory/kappa_samples.nii.gz
mv $directory/Param_3_samples.nii.gz $directory/beta_samples.nii.gz
mv $directory/Param_4_samples.nii.gz $directory/th_samples.nii.gz
mv $directory/Param_5_samples.nii.gz $directory/ph_samples.nii.gz
mv $directory/Param_6_samples.nii.gz $directory/psi_samples.nii.gz

${FSLDIR}/bin/fslmaths $directory/kappa_samples.nii.gz -sub $directory/beta_samples.nii.gz $directory/k2_samples.nii.gz 

Two_div_pi=0.636619772367581

$FSLDIR/bin/fslmaths $directory/fiso_samples.nii.gz -Tmean $directory/mean_fiso
$FSLDIR/bin/fslmaths $directory/fintra_samples.nii.gz -Tmean $directory/mean_fintra
$FSLDIR/bin/fslmaths $directory/kappa_samples.nii.gz -Tmean $directory/mean_kappa
$FSLDIR/bin/fslmaths $directory/beta_samples.nii.gz -Tmean $directory/mean_beta
$FSLDIR/bin/fslmaths $directory/k2_samples.nii.gz -Tmean $directory/mean_k2

$FSLDIR/bin/make_dyadic_vectors $directory/th_samples $directory/ph_samples $directory/nodif_brain_mask.nii.gz dyads1

${bindir}/getFanningOrientation $directory/th_samples $directory/ph_samples $directory/psi_samples $directory/nodif_brain_mask.nii.gz  $directory/fanningDir

${FSLDIR}/bin/fslmaths $directory/mean_k2 -recip -atan -mul $Two_div_pi $directory/ODIp
${FSLDIR}/bin/fslmaths $directory/mean_kappa -recip -atan -mul $Two_div_pi $directory/ODIs
${FSLDIR}/bin/fslmaths $directory/mean_kappa -recip $directory/ODItot_temp1
${FSLDIR}/bin/fslmaths $directory/mean_k2 -recip $directory/ODItot_temp2
${FSLDIR}/bin/fslmaths $directory/ODItot_temp1 -mul $directory/ODItot_temp2 $directory/ODItot_temp1
${FSLDIR}/bin/fslmaths $directory/ODItot_temp1 -sqrt -abs -atan -mul $Two_div_pi $directory/ODItot
rm $directory/ODItot_temp1.nii.gz
rm $directory/ODItot_temp2.nii.gz

#Dispersion Anisotropy Index
${FSLDIR}/bin/fslmaths $directory/mean_beta -div $directory/mean_k2 -atan -mul $Two_div_pi $directory/DA

