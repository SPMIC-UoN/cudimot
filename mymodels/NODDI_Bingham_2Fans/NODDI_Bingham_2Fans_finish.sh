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
    echo "Usage: NODDI_Bingham_2Fans_finish.sh <output_directory>"
    echo ""
    echo "expects to find all the estimatedParameters and nodif_brain_mask in subject directory"
    echo ""
    exit 1
}

[ "$1" = "" ] && Usage

directory=$1
cd ${directory}

mv $directory/Param_0_samples.nii.gz $directory/fiso_samples.nii.gz
mv $directory/Param_1_samples.nii.gz $directory/fFan2_samples.nii.gz
mv $directory/Param_2_samples.nii.gz $directory/fintra_f1_samples.nii.gz
mv $directory/Param_3_samples.nii.gz $directory/kappa_f1_samples.nii.gz
mv $directory/Param_4_samples.nii.gz $directory/beta_f1_samples.nii.gz
mv $directory/Param_5_samples.nii.gz $directory/th_f1_samples.nii.gz
mv $directory/Param_6_samples.nii.gz $directory/ph_f1_samples.nii.gz
mv $directory/Param_7_samples.nii.gz $directory/psi_f1_samples.nii.gz
mv $directory/Param_8_samples.nii.gz $directory/fintra_f2_samples.nii.gz
mv $directory/Param_9_samples.nii.gz $directory/kappa_f2_samples.nii.gz
mv $directory/Param_10_samples.nii.gz $directory/beta_f2_samples.nii.gz
mv $directory/Param_11_samples.nii.gz $directory/th_f2_samples.nii.gz
mv $directory/Param_12_samples.nii.gz $directory/ph_f2_samples.nii.gz
mv $directory/Param_13_samples.nii.gz $directory/psi_f2_samples.nii.gz

${FSLDIR}/bin/fslmaths $directory/kappa_f1_samples.nii.gz -sub $directory/beta_f1_samples.nii.gz $directory/k2_f1_samples.nii.gz 
${FSLDIR}/bin/fslmaths $directory/kappa_f2_samples.nii.gz -sub $directory/beta_f2_samples.nii.gz $directory/k2_f2_samples.nii.gz 

Two_div_pi=0.636619772367581

$FSLDIR/bin/fslmaths $directory/fiso_samples.nii.gz -Tmean $directory/mean_fiso
$FSLDIR/bin/fslmaths $directory/fFan2_samples.nii.gz -Tmean $directory/mean_fFan2
$FSLDIR/bin/fslmaths $directory/fintra_f1_samples.nii.gz -Tmean $directory/mean_fintra_f1
$FSLDIR/bin/fslmaths $directory/kappa_f1_samples.nii.gz -Tmean $directory/mean_kappa_f1
$FSLDIR/bin/fslmaths $directory/beta_f1_samples.nii.gz -Tmean $directory/mean_beta_f1
$FSLDIR/bin/fslmaths $directory/k2_f1_samples.nii.gz -Tmean $directory/mean_k2_f1
$FSLDIR/bin/fslmaths $directory/fintra_f2_samples.nii.gz -Tmean $directory/mean_fintra_f2
$FSLDIR/bin/fslmaths $directory/kappa_f2_samples.nii.gz -Tmean $directory/mean_kappa_f2
$FSLDIR/bin/fslmaths $directory/beta_f2_samples.nii.gz -Tmean $directory/mean_beta_f2
$FSLDIR/bin/fslmaths $directory/k2_f2_samples.nii.gz -Tmean $directory/mean_k2_f2

$FSLDIR/bin/make_dyadic_vectors $directory/th_f1_samples $directory/ph_f1_samples $directory/nodif_brain_mask.nii.gz dyads1
$FSLDIR/bin/make_dyadic_vectors $directory/th_f2_samples $directory/ph_f2_samples $directory/nodif_brain_mask.nii.gz dyads2

${bindir}/getFanningOrientation $directory/th_f1_samples $directory/ph_f1_samples $directory/psi_f1_samples $directory/nodif_brain_mask.nii.gz  $directory/fanningDir_f1
${bindir}/getFanningOrientation $directory/th_f2_samples $directory/ph_f2_samples $directory/psi_f2_samples $directory/nodif_brain_mask.nii.gz  $directory/fanningDir_f2

${FSLDIR}/bin/fslmaths $directory/mean_k2_f1 -recip -atan -mul $Two_div_pi $directory/ODIp_f1
${FSLDIR}/bin/fslmaths $directory/mean_kappa_f1 -recip -atan -mul $Two_div_pi $directory/ODIs_f1
${FSLDIR}/bin/fslmaths $directory/mean_kappa_f1 -recip $directory/ODItot_temp1
${FSLDIR}/bin/fslmaths $directory/mean_k2_f1 -recip $directory/ODItot_temp2
${FSLDIR}/bin/fslmaths $directory/ODItot_temp1 -mul $directory/ODItot_temp2 $directory/ODItot_temp1
${FSLDIR}/bin/fslmaths $directory/ODItot_temp1 -sqrt -abs -atan -mul $Two_div_pi $directory/ODItot_f1
rm $directory/ODItot_temp1.nii.gz
rm $directory/ODItot_temp2.nii.gz
${FSLDIR}/bin/fslmaths $directory/mean_k2_f2 -recip -atan -mul $Two_div_pi $directory/ODIp_f2
${FSLDIR}/bin/fslmaths $directory/mean_kappa_f2 -recip -atan -mul $Two_div_pi $directory/ODIs_f2
${FSLDIR}/bin/fslmaths $directory/mean_kappa_f2 -recip $directory/ODItot_temp1
${FSLDIR}/bin/fslmaths $directory/mean_k2_f2 -recip $directory/ODItot_temp2
${FSLDIR}/bin/fslmaths $directory/ODItot_temp1 -mul $directory/ODItot_temp2 $directory/ODItot_temp1
${FSLDIR}/bin/fslmaths $directory/ODItot_temp1 -sqrt -abs -atan -mul $Two_div_pi $directory/ODItot_f2
rm $directory/ODItot_temp1.nii.gz
rm $directory/ODItot_temp2.nii.gz


#Dispersion Anisotropy Index
${FSLDIR}/bin/fslmaths $directory/mean_beta_f1 -div $directory/mean_k2_f1 -atan -mul $Two_div_pi $directory/DA_f1
${FSLDIR}/bin/fslmaths $directory/mean_beta_f2 -div $directory/mean_k2_f2 -atan -mul $Two_div_pi $directory/DA_f2

