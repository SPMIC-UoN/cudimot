#!/bin/sh
#
#   Moises Hernandez-Fernandez - FMRIB Image Analysis Group
#
#   Copyright (C) 2004 University of Oxford
#
#   SHCOPYRIGHT
#
# Pipeline for fitting Ball-2Sticks Multiexponential as Bedpostx

if [ "x$CUDIMOT" == "x" ]; then
	echo ""
	echo "Please, set enviroment variable CUDIMOT with the path where cuDIMOT is installed"
	echo "The path must contain a bin directory with binary files, i.e. \$CUDIMOT/bin"
	echo "For instance:   export CUDIMOT=/home/moises/CUDIMOT"
	echo ""
  exit 1
fi

bindir=${CUDIMOT}/bin

modelname=Ball_2_Sticks_Multiexponential
step1=Ball_2_Sticks

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
Usage() {
    echo ""
    echo "Usage: Pipeline_Ball_2_Sticks_Multiexponential.sh <subject_directory> [options]"
    echo ""
    echo "expects to find data and nodif_brain_mask in subject directory"
    echo ""
    echo "<options>:"
    echo "-Q (name of the GPU(s) queue, default cuda.q (defined in environment variable: FSLGECUDAQ)"
    echo "-NJOBS (number of jobs to queue, the data is divided in NJOBS parts, usefull for a GPU cluster, default 4)"
    echo "--runMCMC (if you want to run MCMC)"
    echo "-b (burnin period, default 5000)"
    echo "-j (number of jumps, default 1250)"
    echo "-s (sample every, default 25)"
    echo "--BIC_AIC (calculate BIC & AIC)"
    echo ""
    exit 1
}

[ "$1" = "" ] && Usage

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${FSLDIR}/lib

subjdir=`make_absolute $1`
subjdir=`echo $subjdir | sed 's/\/$/$/g'`

echo "---------------------------------------------------------------------------------"
echo "------------------------------------ CUDIMOT ------------------------------------"
echo "----------------------------- MODEL: $modelname -----------------------------"
echo "---------------------------------------------------------------------------------"
echo subjectdir is $subjdir

start=`date +%s`

#parse option arguments
njobs=4
burnin=1000
njumps=1250
sampleevery=25
other=""
queue=""
lastStepModelOpts=""

shift
while [ ! -z "$1" ]
do
  case "$1" in
      -Q) queue="-q $2";shift;;
      -NJOBS) njobs=$2;shift;;
      -b) burnin=$2;shift;;
      -j) njumps=$2;shift;;
      -s) sampleevery=$2;shift;;
      --runMCMC) lastStepModelOpts=$lastStepModelOpts" --runMCMC";;
      --BIC_AIC) lastStepModelOpts=$lastStepModelOpts" --BIC_AIC";;
      *) other=$other" "$1;;
  esac
  shift
done

#Set options
opts="--bi=$burnin --nj=$njumps --se=$sampleevery"
opts="$opts $other"

if [ "x$SGE_ROOT" != "x" ]; then
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
mkdir -p ${subjdir}.${modelname}/Dtifit
mkdir -p ${subjdir}.${modelname}/${step1}
mkdir -p ${subjdir}.${modelname}/${step1}/diff_parts
mkdir -p ${subjdir}.${modelname}/${step1}/logs

echo Copying files to output directory

${FSLDIR}/bin/imcp ${subjdir}/nodif_brain_mask ${subjdir}.${modelname}
if [ `${FSLDIR}/bin/imtest ${subjdir}/nodif` = 1 ] ; then
    ${FSLDIR}/bin/fslmaths ${subjdir}/nodif -mas ${subjdir}/nodif_brain_mask ${subjdir}.${modelname}/nodif_brain
fi

# Specify Common Fixed Parameters
CFP_file=$subjdir.${modelname}/CFP
cp ${subjdir}/bvecs $subjdir.${modelname}
cp ${subjdir}/bvals $subjdir.${modelname}
echo  $subjdir.${modelname}/bvecs > $CFP_file
echo  $subjdir.${modelname}/bvals >> $CFP_file

#Set more options
opts=$opts" --data=${subjdir}/data --maskfile=$subjdir.${modelname}/nodif_brain_mask --forcedir --CFP=$CFP_file"

##############################################################################
################################ First Dtifit  ###############################
##############################################################################
echo "Queue Dtifit"
PathDTI=${subjdir}.${modelname}/Dtifit
dtifit_command="${bindir}/Run_dtifit.sh ${subjdir} ${subjdir}.${modelname} ${bindir}"
#SGE
dtifitProcess=`${FSLDIR}/bin/fsl_sub $queue -l $PathDTI/logs -N dtifit $dtifit_command`

#############################################################
################ Ball & 2 Sticks (simple)  ##################
#############################################################
echo "${step1} fitting process"
PathStep1=$subjdir.${modelname}/${step1}

# Create file to specify initialisation parameters
InitializationFile=$PathStep1/InitializationParameters
echo ${PathDTI}/dtifit_S0.nii.gz  > $InitializationFile # S0
echo ${PathDTI}/dtifit_L1.nii.gz >> $InitializationFile # d  // empirically found that d~L1
echo ${PathDTI}/dtifit_FA.nii.gz >> $InitializationFile # f1
echo ${PathDTI}/dtifit_V1_th.nii.gz >> $InitializationFile # th1
echo ${PathDTI}/dtifit_V1_ph.nii.gz >> $InitializationFile # ph1
echo ${PathDTI}/dtifit_FA_div2.nii.gz >> $InitializationFile # f2
echo ${PathDTI}/dtifit_V2_th.nii.gz >> $InitializationFile # th2
echo ${PathDTI}/dtifit_V2_ph.nii.gz >> $InitializationFile # ph2

partsdir=$PathStep1/diff_parts
outputdir=$PathStep1
Step1Opts=$opts" --outputdir=$outputdir --partsdir=$partsdir --init_params=$InitializationFile"

postproc=`${bindir}/jobs_wrapper.sh $PathStep1 $dtifitProcess $step1 Fit $njobs $Step1Opts`
###############################################################
############ Ball & 2 Sticks - Multiexponential ###############
###############################################################
echo "$modelname fitting process"

# Create file to specify initialization parameters
InitializationFile=$subjdir.${modelname}/InitializationParameters
echo $PathStep1/Param_0_samples > $InitializationFile # S0
echo $PathStep1/Param_1_samples >> $InitializationFile # d
echo $PathStep1/Param_1_samples >> $InitializationFile # d_std = d
echo $PathStep1/Param_2_samples >> $InitializationFile # f1
echo $PathStep1/Param_3_samples >> $InitializationFile # th1
echo $PathStep1/Param_4_samples >> $InitializationFile # ph1
echo $PathStep1/Param_5_samples >> $InitializationFile # f2
echo $PathStep1/Param_6_samples >> $InitializationFile # th2
echo $PathStep1/Param_7_samples >> $InitializationFile # ph2

partsdir=${subjdir}.${modelname}/diff_parts
outputdir=${subjdir}.${modelname}
ModelOpts=$opts" --outputdir=$outputdir --partsdir=$partsdir --init_params=$InitializationFile $lastStepModelOpts"

postproc=`${bindir}/jobs_wrapper.sh ${subjdir}.${modelname} $postproc $modelname FitProcess $njobs $ModelOpts`

#########################################
############## PostProc #################
#########################################
finish_command="${bindir}/${modelname}_finish.sh ${subjdir}.${modelname}"
#SGE
finishProcess=`${FSLDIR}/bin/fsl_sub $queue -l ${subjdir}.${modelname}/logs -N ${modelname}_finish -j $postproc $finish_command`

endt=`date +%s`
runtime=$((endt-start))
#echo Runtime $runtime
echo Everything Queued
