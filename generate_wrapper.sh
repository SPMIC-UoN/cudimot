#!/bin/bash

# Generate Wrapper
scriptFile=cudimot_${modelname}.sh
echo "#!/bin/sh" > $scriptFile
echo "#" >> $scriptFile
echo "# Created with CUDIMOT: Copyright (C) 2004 University of Oxford" >> $scriptFile
echo "# Moises Hernandez-Fernandez - FMRIB Image Analysis Group" >> $scriptFile
echo "#" >> $scriptFile
echo "" >> $scriptFile

echo "modelname=${modelname}" >> $scriptFile
cat wrapper_template.sh >> $scriptFile
chmod 755 $scriptFile
