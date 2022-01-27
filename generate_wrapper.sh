#!/bin/bash

if [ "$FSLDEVDIR" = "" ]; then
    echo "Error: The variable FSLDEVDIR must be set"
    exit 1
fi
if [ ! -d "$FSLDEVDIR" ]; then
    echo "Error: The directory $FSLDEVDIR specified in FSLDEVDIR does not exist"
    exit 1
fi

# Copy priors file into the binary directory
cp mymodels/$modelname/modelpriors $FSLDEVDIR/bin/${modelname}_priors

# Generate Wrapper
scriptFile=$FSLDEVDIR/bin/cudimot_${modelname}.sh
echo "#!/bin/sh" > $scriptFile
echo "#" >> $scriptFile
echo "# Created with CUDIMOT: Copyright (C) 2004 University of Oxford" >> $scriptFile
echo "# Moises Hernandez-Fernandez - FMRIB Image Analysis Group" >> $scriptFile
echo "#" >> $scriptFile
echo "" >> $scriptFile

echo "modelname=${modelname}" >> $scriptFile
cat wrapper_template.sh >> $scriptFile
chmod 755 $scriptFile

# No text file with info about the model
#outDIR=`cat Makefile | grep DIR_objs= | cut -f 2 -d "="`
infoFile=$FSLDEVDIR/bin/$modelname".info"

echo > $infoFile
echo "---------------------------------------------" >> $infoFile
echo "---------------- NPARAMETERS  ---------------" >> $infoFile
echo "---------------------------------------------" >> $infoFile
cat mymodels/$modelname/modelparameters.h | grep define | grep -v INCLUDED >> $infoFile
echo >> $infoFile
cat mymodels/$modelname/modelparameters.cc | grep int >> $infoFile
echo "---------------------------------------------" >> $infoFile
echo >> $infoFile
cat mymodels/$modelname/modelparameters.h | grep type >> $infoFile
echo >> $infoFile
echo >> $infoFile

echo "-------------------------------------------------------------" >> $infoFile
echo "--- Predicted Signal, Constraints & Derivatives Functions ---" >> $infoFile
echo "--- FILE: "$functionsFile" ---" >> $infoFile
echo "-------------------------------------------------------------" >> $infoFile
cat mymodels/$modelname/modelfunctions.h >> $infoFile
