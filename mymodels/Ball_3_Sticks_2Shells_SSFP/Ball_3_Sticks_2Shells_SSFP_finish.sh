#!/bin/sh
# modelname=Ball_3_Sticks_2Shells_SSFP
# Author: MKS mohamed.selim@nottingham.ac.uk"


if [ ! $# -eq 1 ]; then
  echo 1>&2 "$0: requires DW-SSFP bedpostX results directory"
  exit 2
fi

outdir=$1
cd $outdir

mkdir -p cudimot_Param
mv Param* cudimot_Param 2> /dev/null

fslmaths cudimot_Param/Param_0_samples -Tmean mean_d1samples
fslmaths cudimot_Param/Param_1_samples -Tmean mean_d2samples

fslmaths cudimot_Param/Param_2_samples -Tmean mean_f1samples
fslmaths cudimot_Param/Param_5_samples -Tmean mean_f2samples
fslmaths cudimot_Param/Param_8_samples -Tmean mean_f3samples # added 3rd fibre

make_dyadic_vectors cudimot_Param/Param_3_samples cudimot_Param/Param_4_samples nodif_brain_mask dyads1
make_dyadic_vectors cudimot_Param/Param_6_samples cudimot_Param/Param_7_samples nodif_brain_mask dyads2
make_dyadic_vectors cudimot_Param/Param_9_samples cudimot_Param/Param_10_samples nodif_brain_mask dyads3 # added 3rd fibre

fslmaths cudimot_Param/Param_11_samples -Tmean mean_S0samples

fslmaths mean_f2samples -thr 0.05 -bin grot
fslmaths dyads2 -mas grot dyads2_thr0.05
rm grot.nii.gz

# added 3rd fibre
fslmaths mean_f3samples -thr 0.05 -bin grot
fslmaths dyads3 -mas grot dyads3_thr0.05
rm grot.nii.gz

fslmaths mean_f1samples -add  mean_f2samples -add  mean_f2samples mean_fsumsamples
fslmaths mean_f2samples -div mean_f1samples mean_f2_f1
fslmaths mean_f2samples -div mean_f1samples mean_f3_f1

mv cudimot_Param/Param_0_samples.nii.gz ./merged_d1samples.nii.gz
mv cudimot_Param/Param_1_samples.nii.gz ./merged_d2samples.nii.gz
mv cudimot_Param/Param_2_samples.nii.gz ./merged_f1samples.nii.gz
mv cudimot_Param/Param_3_samples.nii.gz ./merged_th1samples.nii.gz
mv cudimot_Param/Param_4_samples.nii.gz ./merged_ph1samples.nii.gz
mv cudimot_Param/Param_5_samples.nii.gz ./merged_f2samples.nii.gz
mv cudimot_Param/Param_6_samples.nii.gz ./merged_th2samples.nii.gz
mv cudimot_Param/Param_7_samples.nii.gz ./merged_ph2samples.nii.gz
mv cudimot_Param/Param_8_samples.nii.gz ./merged_f3samples.nii.gz # added 3rd fibre
mv cudimot_Param/Param_9_samples.nii.gz ./merged_th3samples.nii.gz # added 3rd fibre
mv cudimot_Param/Param_10_samples.nii.gz ./merged_ph3samples.nii.gz # added 3rd fibre
mv cudimot_Param/Param_11_samples.nii.gz ./merged_S0samples.nii.gz

#rm -r ./cudimot_Param

echo "Done"
