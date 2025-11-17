#!/bin/sh
#
# Created with CUDIMOT: Copyright (C) 2004 University of Oxford
# Moises Hernandez-Fernandez - FMRIB Image Analysis Group
#

modelname=Ball_3_Sticks_1Shell_SSFP

bindir=${FSLDEVDIR}/bin
if [ ! -f $bindir/$modelname ]; then
    echo "Please set a correct modelname. Not binaries found for model: $modelname"
    exit 1
fi

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${FSLDIR}/lib

Usage() {
    echo ""
    echo "Usage: cuditmot <subject_directory> [options]"
    echo ""
    echo "expects to find data and nodif_brain_mask in subject directory"
    echo ""
    echo "<options>:"
    echo "-waitfor (job_ID)"
    echo "-Q (name of the GPU(s) queue, default cuda.q (defined in environment variable: FSLGECUDAQ)"
    echo "-NJOBS (number of jobs to queue, the data is divided in NJOBS parts, usefull for a GPU cluster, default 4)"
    echo "--no_LevMar (Do not run Levenberg-Marquardt)"
    echo "--runMCMC (Run MCMC)"
    echo "--CFP=filePath (Specify path of the file with the list of common fixed parameters ascii files)"
    echo "--FixP=filePath (Specify path of the file with the list of fixed parameters NIfTI files)"
    echo "-b (burnin period, default 5000)"
    echo "-j (number of jumps, default 1250)"
    echo "-s (sample every, default 25)"
    echo ""
    exit 1
}

make_absolute(){
    dir=$1;
    if [ -d ${dir} ]; then
	OLDWD=`pwd`
	cd ${dir}
	dir_all=`pwd`
	cd $OLDWD
    else
	dir_all=${dir}
    fi
    echo ${dir_all}
}

[ "$1" = "" ] && Usage
if [ "$modelname" = "" ]; then
    echo "Error: The variable modelname must be set"
    exit 1
fi

subjdir=`make_absolute $1`
subjdir=`echo $subjdir | sed 's/\/$/$/g'`

echo "---------------------------------------------------------------------------------"
echo "------------------------------------ CUDIMOT ------------------------------------"
echo "----------------------------- MODEL: $modelname -----------------------------"
echo "---------------------------------------------------------------------------------"
echo subjectdir is $subjdir

#parse option arguments
qsys=0
njobs=4
burnin=1000
njumps=1250
sampleevery=25
other=""
queue=""
wait=""

shift
while [ ! -z "$1" ]
do
  case "$1" in
  	  -waitfor) wait="-j $2";shift;;
      -Q) queue="-q $2";shift;;
      -NJOBS) njobs=$2;shift;;
      -b) burnin=$2;shift;;
      -j) njumps=$2;shift;;
      -s) sampleevery=$2;shift;;
      *) other=$other" "$1;;
  esac
  shift
done

#Set options
opts="--bi=$burnin --nj=$njumps --se=$sampleevery"
opts="$opts $other"

if [ $qsys -eq 0 ] && [ "x$SGE_ROOT" != "x" ]; then
	queue="-q $FSLGECUDAQ"
fi


#check that all required files exist

if [ ! -d $subjdir ]; then
	echo "subject directory $1 not found"
	exit 1
fi

if [ `${FSLDIR}/bin/imtest ${subjdir}/data` -eq 0 ]; then
	echo "${subjdir}/data not found"
	exit 1
fi

if [ `${FSLDIR}/bin/imtest ${subjdir}/nodif_brain_mask` -eq 0 ]; then
	echo "${subjdir}/nodif_brain_mask not found"
	exit 1
fi

if [ -e ${subjdir}.${modelname}/xfms/eye.mat ]; then
	echo "${subjdir} has already been processed: ${subjdir}.${modelname}." 
	echo "Delete or rename ${subjdir}.${modelname} before repeating the process."
	exit 1
fi

echo Making output directory structure

mkdir -p ${subjdir}.${modelname}/
mkdir -p ${subjdir}.${modelname}/diff_parts
mkdir -p ${subjdir}.${modelname}/logs
part=0

#mkdir -p ${subjdir}.${modelname}/logs/logs_gpu
partsdir=${subjdir}.${modelname}/diff_parts

echo Copying files to output directory

${FSLDIR}/bin/imcp ${subjdir}/nodif_brain_mask ${subjdir}.${modelname}
if [ `${FSLDIR}/bin/imtest ${subjdir}/nodif` = 1 ] ; then
    ${FSLDIR}/bin/fslmaths ${subjdir}/nodif -mas ${subjdir}/nodif_brain_mask ${subjdir}.${modelname}/nodif_brain
fi

############################################
########### following is for initialisation - added by MKS
############################################

########### use dtifit for diffusivity estimates
# echo "Queue Dtifit"
# PathDTI=${subjdir}.${modelname}/Dtifit
# dtifit_command="${bindir}/Run_dtifit.sh ${subjdir} ${subjdir}.${modelname} ${bindir}"
# #SGE
# dtifitProcess=`${FSLDIR}/bin/fsl_sub $queue -l $PathDTI/logs -N dtifit $dtifit_command`
echo "Running dtifit..."
PathDTI=${subjdir}.${modelname}/Dtifit
${bindir}/Run_dtifit.sh ${subjdir} ${subjdir}.${modelname} ${bindir}

########### check bedpostX files exist
if [ ! -d ${subjdir}.pvm ]; then
	echo "subject bedpostX directory $1 not found"
	exit 1
fi

# using the following as initialisation parameters:
# scaled MD from DTIFTI for d1,d2,d3
# dyads: dyads1 dyads2 dyads3 for theta and phi
# mean_f1samples mean_f2sample mean_f3sample for volumes
# mean_S0samples
for img in dyads1 dyads2 dyads3 f1 f2 f3 S0; do
  if [ `${FSLDIR}/bin/imtest ${subjdir}.pvm/${img}` -eq 0 ]; then
  	echo "${subjdir}.pvm/${img} not found"
  	exit 1
  fi
done


########### Create file to specify initialisation parameters
InitializationFile=$subjdir.${modelname}/InitializationParameters
InitializationDir=$subjdir.${modelname}/InitializationDir

mkdir -p $InitializationDir

# here, do some scaling of MD from dtifit to get approx estimates of d1, d2
fslmaths ${subjdir}.${modelname}/Dtifit/dtifit_MD -mul 2 ${InitializationDir}/d.nii.gz

# from spin-echo bedpostX
for img in dyads1 dyads2 dyads3 f1 f2 f3 S0; do
  ${FSLDIR}/bin/imcp ${subjdir}.pvm/${img} $InitializationDir
done

# convert dyads spherical to cartesian and use th/ph as params below
for img in dyads1 dyads2 dyads3; do
    echo "converting ${img} to spherical"
    fslpython /home/mszms12/dyad2sph.py -dyad ${InitializationDir}/${img}.nii.gz -mask ${subjdir}.${modelname}/nodif_brain_mask.nii.gz
done

echo ${InitializationDir}/d >> $InitializationFile # d P[0]
echo ${InitializationDir}/f1 >> $InitializationFile # f1 P[1]
echo ${InitializationDir}/th1 >> $InitializationFile # th1 P[2]
echo ${InitializationDir}/ph1 >> $InitializationFile # ph1 P[3]
echo ${InitializationDir}/f2 >> $InitializationFile # f2 P[4]
echo ${InitializationDir}/th2 >> $InitializationFile # th2 P[5]
echo ${InitializationDir}/ph2 >> $InitializationFile # ph2 P[6]
echo ${InitializationDir}/f3 >> $InitializationFile # f3 P[7]
echo ${InitializationDir}/th3 >> $InitializationFile # th3 P[8]
echo ${InitializationDir}/ph3 >> $InitializationFile # ph3 P[9]
echo ${InitializationDir}/S0 >> $InitializationFile # S0 P[10]

########### resume normal script
############################################
############################################
#Set more default options
opts=$opts" --data=${subjdir}/data --maskfile=$subjdir.${modelname}/nodif_brain_mask --partsdir=$partsdir --outputdir=$subjdir.${modelname} --forcedir"

# Split the dataset in parts
echo Pre-processing stage
	PreprocOpts=$opts" --idPart=0 --nParts=$njobs --logdir=$subjdir.${modelname}/logs/preProcess"
	preproc_command="$bindir/split_parts_${modelname} $PreprocOpts"

	#SGE
	preProcess=`${FSLDIR}/bin/fsl_sub $wait $queue -l ${subjdir}.${modelname}/logs -N ${modelname}_preproc $preproc_command`

echo Queuing Fitting model processing stage

	[ -f ${subjdir}.${modelname}/commands.txt ] && rm ${subjdir}.${modelname}/commands.txt

	part=0
	while [ $part -lt $njobs ]
	do
	    	partzp=`$FSLDIR/bin/zeropad $part 4`
	    
		Fitopts=$opts

		#${FSLDIR}/bin/
		echo "$bindir/${modelname} --idPart=$part --nParts=$njobs --logdir=$subjdir.${modelname}/logs/${modelname}_$partzp $Fitopts" >> ${subjdir}.${modelname}/commands.txt
	    
	    	part=$(($part + 1))
	done

	#SGE
	FitProcess=`${FSLDIR}/bin/fsl_sub $queue -N ${modelname} -j $preProcess -t ${subjdir}.${modelname}/commands.txt -l ${subjdir}.${modelname}/logs`

echo Queuing Post-processing stage
# Needs the parent directory where all the output parts are stored $subjdir.${modelname}
PostprocOpts=$opts" --idPart=0 --nParts=$njobs --logdir=$subjdir.${modelname}/logs/postProcess"

#${FSLDIR}/bin/
postproc_command="$bindir/merge_parts_${modelname} $PostprocOpts"

if [ $FitProcess -eq $FitProcess 2>/dev/null ]; then
    #if not error in the prevoius fsl_sub, Fitprocess is a number
    #SGE
    postProcess=`${FSLDIR}/bin/fsl_sub $queue -j $FitProcess -N ${modelname}_postproc_gpu -l ${subjdir}.${modelname}/logs $postproc_command`
fi

