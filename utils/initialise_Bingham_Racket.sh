#!/bin/sh
#
#   Moises Hernandez-Fernandez - FMRIB Image Analysis Group
#
#   Copyright (C) 2004 University of Oxford
#
#   Part of FSL - FMRIB's Software Library
#   http://www.fmrib.ox.ac.uk/fsl
#   fsl@fmrib.ox.ac.uk
#
#   Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
#   Imaging of the Brain), Department of Clinical Neurology, Oxford
#   University, Oxford, UK
#

export LC_ALL=C
#
# Script for initialise Psi angle and k2tok1 in Ball&Racket

Usage() {
    echo ""
    echo "Usage: initialise_Bingham <bindir> <DTI_output_dir> <mask_NIfTI_file> <output_dir>"
    echo ""
    exit 1
}

[ "$4" = "" ] && Usage

bindir=$1
PathDTI=$2
mask=$3
outputdir=$4

${bindir}/initialise_Psi ${PathDTI}/dtifit_V1.nii.gz ${PathDTI}/dtifit_V2.nii.gz ${mask} ${outputdir}/initialPsi

#k2tok1 = (eigs(1)/eigs(2))^2;
#if L2 is too low-> k2tok1 = 1
#${FSLDIR}/bin/fslmaths ${PathDTI}/dtifit_L2.nii.gz -uthr 0.0001 -bin ${outputdir}/temp1
#${FSLDIR}/bin/fslmaths ${PathDTI}/dtifit_L3.nii.gz -mul ${outputdir}/temp1 ${outputdir}/temp1
#${FSLDIR}/bin/fslmaths ${PathDTI}/dtifit_L2.nii.gz -thr 0.0001 -bin ${outputdir}/temp2
#${FSLDIR}/bin/fslmaths ${PathDTI}/dtifit_L2.nii.gz -mul ${outputdir}/temp2 -mul 1.1 -add ${outputdir}/temp1 ${outputdir}/temp1

#${FSLDIR}/bin/fslmaths ${PathDTI}/dtifit_L3.nii.gz -div ${outputdir}/temp1 -sqr ${outputdir}/k2tok1

#rm ${outputdir}/temp1.nii.gz
#rm ${outputdir}/temp2.nii.gz
