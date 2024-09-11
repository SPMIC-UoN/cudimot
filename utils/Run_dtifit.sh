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
# Run Dtifit and calculate the spherical coordinates from the tensors

bindir=${CUDIMOT}/bin

Usage() {
    echo ""
    echo "Usage: Run_dtifit <subject_directory> <output>"
    echo ""
    echo "expects to find data, nodif_brain_mask, bvals and bvecs in subject directory"
    echo ""
    exit 1
}

[ "$3" = "" ] && Usage

subjdir=$1
output=$2

${FSLDIR}/bin/dtifit -k ${subjdir}/data -m ${subjdir}/nodif_brain_mask -r ${subjdir}/bvecs -b ${subjdir}/bvals -o ${output}/Dtifit/dtifit --save_tensor

PathDTI=${output}/Dtifit

# calculate Spherical Coordinates: th1 and ph1
${bindir}/cart2spherical ${PathDTI}/dtifit_V1 ${PathDTI}/dtifit_V1
${bindir}/cart2spherical ${PathDTI}/dtifit_V2 ${PathDTI}/dtifit_V2

# for f2
${FSLDIR}/bin/fslmaths ${PathDTI}/dtifit_FA -div 2 ${PathDTI}/dtifit_FA_div2
${FSLDIR}/bin/fslmaths ${PathDTI}/dtifit_FA -div 4 ${PathDTI}/dtifit_FA_div4
