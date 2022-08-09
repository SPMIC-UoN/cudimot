#!/bin/bash

MODELNAME=$1
SRCDIR=`dirname "$0"`

# Generate Wrapper
scriptFile=cudimot_${MODELNAME}.sh
echo "#!/bin/sh" > $scriptFile
echo "#" >> $scriptFile
echo "# Created with CUDIMOT: Copyright (C) 2004 University of Oxford" >> $scriptFile
echo "# Moises Hernandez-Fernandez - FMRIB Image Analysis Group" >> $scriptFile
echo "#" >> $scriptFile
echo "" >> $scriptFile

echo "modelname=${MODELNAME}" >> $scriptFile
cat ${SRCDIR}/wrapper_template.sh >> $scriptFile
chmod 755 $scriptFile
