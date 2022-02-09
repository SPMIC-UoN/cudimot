# This is the Makefile for the CUDIMOT project
#
# The CUDIMOT project provides ...
#
# A CUDA compiler and toolkit must be installed

include $(FSLCONFDIR)/default.mk

PROJNAME = cudimot
MODELNAME = NODDI_Bingham
MODELDIR = mymodels/${MODELNAME}

LIBS     = -lfsl-warpfns -lfsl-basisfield -lfsl-meshclass \
           -lfsl-newimage -lfsl-miscmaths -lfsl-NewNifti \
           -lfsl-utils -lfsl-newran -lfsl-znz -lfsl-cprob

CUDIMOT_CUDA_OBJS = \
    ${MODELDIR}/modelparameters.o \
    init_gpu.o \
    dMRI_Data.o \
    Model.o \
    Parameters.o \
    GridSearch.o \
    Levenberg_Marquardt.o \
    MCMC.o \
    BIC_AIC.o \
    getPredictedSignal.o

CUDIMOT_OBJS = \
    cudimot.o \
    cudimotoptions.o

USRINCFLAGS = -I${MODELDIR} -I$(FSLDIR)/include/armawrap
SCRIPTS  = \
    utils/Run_dtifit.sh \
    utils/jobs_wrapper.sh \
    utils/initialise_Bingham.sh \
    ${MODELDIR}/${MODELNAME}_finish.sh \
    ${MODELDIR}/Pipeline_${MODELNAME}.sh
XFILES   = \
    cart2spherical \
    getFanningOrientation \
    initialise_Psi \
    ${MODELNAME} \
    cudimot_${MODELNAME}.sh \
    merge_parts_${MODELNAME} \
    split_parts_${MODELNAME} \
    testFunctions_${MODELNAME}
DATAFILES = \
    ${MODELNAME}_priors \
    ${MODELNAME}.info

all: ${XFILES} ${DATAFILES}
clean: cleandata

cleandata:
	rm -rf ${DATAFILES}

%.o: %.cu
	${NVCC} -c ${NVCCFLAGS} $^

cart2spherical: utils/cart2spherical.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}

getFanningOrientation: utils/getFanningOrientation.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}

initialise_Psi: utils/initialise_Psi.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}

${MODELNAME}: cudimotoptions.o cudimot.cc $(CUDIMOT_CUDA_OBJS)
	${NVCC} ${NVCCFLAGS} ${NVCCLDFLAGS} -o $@ $^

merge_parts_${MODELNAME}: merge_parts.cc cudimotoptions.o $(CUDIMOT_CUDA_OBJS)
	${NVCC} ${NVCCFLAGS} ${NVCCLDFLAGS} -o $@ $^ -lboost_filesystem -lboost_system

split_parts_${MODELNAME}: split_parts.cc cudimotoptions.o $(CUDIMOT_CUDA_OBJS)
	${NVCC} ${NVCCFLAGS} $(NVCCLDFLAGS) -o $@ $^ -lboost_filesystem -lboost_system

testFunctions_${MODELNAME}: ${MODELDIR}/modelparameters.o testFunctions.cu
	${NVCC} ${NVCCFLAGS} ${NVCCLDFLAGS} -o $@ $^

cudimot_${MODELNAME}.sh : ${MODELNAME}
	./generate_wrapper.sh ${MODELNAME}

${MODELNAME}_priors:
	cp ${MODELDIR}/modelpriors $@

${MODELNAME}.info : cudimot_${MODELNAME}.sh
	./generate_info.sh ${MODELNAME}
