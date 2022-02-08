#!/bin/bash

# No text file with info about the model
#outDIR=`cat Makefile | grep DIR_objs= | cut -f 2 -d "="`
infoFile=$modelname".info"

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
