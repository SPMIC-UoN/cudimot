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

export LC_ALL=C
#
# Script for initialise k2 from k1 and ration k1/k2 in Ball & Rackets

Usage() {
    echo ""
    echo "Usage: initialise_dzepplin_perp.sh <dpar_file> <frac_file> <output_dper>"
    echo ""
    exit 1
}

[ "$3" = "" ] && Usage

dpar=$1
f=$2
output=$3

#dperp = dpar * (1-f)
${FSLDIR}/bin/fslmaths $dpar -mul $f $output
${FSLDIR}/bin/fslmaths $dpar -sub $output $output
 
