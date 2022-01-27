#!/bin/bash
. $FSLDIR/etc/fslconf/fsl.sh
export FSLCONFDIR=$FSLDIR/config
export FSLMACHTYPE=`$FSLDIR/etc/fslconf/fslmachtype.sh`
export CUDA=/opt/cuda-7.5
unset SGE_ROOT
