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

if [ ! $(type -P nvcc) ]; then
  echo ""
  echo "Please ensure that the nvcc compiler command is on your PATH"
  echo "For instance: export PATH=/usr/local/cuda-8.0/bin:$PATH"
  echo ""
  exit 1
fi

. $FSLDIR/etc/fslconf/fsl-devel.sh


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
