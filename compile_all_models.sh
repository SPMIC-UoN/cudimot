#!/bin/bash
if [ "x$FSLDIR" == "x" ]; then
  echo ""
  echo "Please set environment variable FSLDIR before compiling"
  echo ""
	exit 1
fi
if [ "x$FSLDEVDIR" == "x" ]; then
  echo ""
  echo "Please set environment variable FSLDEVDIR before compiling"
  echo ""
	exit 1
fi
if [ "x$CUDA" == "x" ]; then
  echo ""
  echo "Please set enviroment variable CUDA with the path to the version of CUDA to use"
  echo "For instance: export CUDA=/usr/local/cuda-8.0" 
  echo ""
	exit 1
fi

. $FSLDIR/etc/fslconf/fsl.sh
export FSLCONFDIR=$FSLDIR/config
export FSLMACHTYPE=`$FSLDIR/etc/fslconf/fslmachtype.sh`

set -e
models=`ls mymodels -I utils`
for m in $models; do
    echo
    echo ------------------------------
    echo -------- Compiling $m --------
    echo ------------------------------
    echo
    export modelname=$m
    make cleanall
    make install
done
