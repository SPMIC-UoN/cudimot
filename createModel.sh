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

if [ -z $modelname ]; then
    echo "Please set modelname"
    exit 1
fi

if [ -d mymodels/$modelname ]; then
    echo "This model has already been created. Please, select a different modelname or remove the directory of the model: mymodels/$modelname"
    exit 1
fi

cp -r mymodels/template mymodels/$modelname

echo ""
echo "Model created sucesfully in mymodels/$modelname."
echo "Please edit the files inside this directory to specify your model parameters and functions."
echo ""
