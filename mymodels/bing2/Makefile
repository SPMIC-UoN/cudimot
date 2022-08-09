# This is the Makefile for the CUDIMOT project
#
# The CUDIMOT project provides ...
#
# A CUDA compiler and toolkit must be installed

include $(FSLCONFDIR)/default.mk

PROJNAME = cudimot
MODELNAME = NODDI_Bingham
MODELDIR = /home/bbzmsc/fslconda/cudimot/mymodels/bing2
CUDIMOT_SRCDIR = $(MODELDIR)/../..

LIBS     = -lfsl-warpfns -lfsl-basisfield -lfsl-meshclass  \
           -lfsl-newimage -lfsl-miscmaths -lfsl-NewNifti   \
           -lfsl-utils -lfsl-newran -lfsl-znz -lfsl-cprob

USRINCFLAGS = -I$(MODELDIR) -I$(FSLDIR)/include/armawrap

CUDIMOT_CUDA_OBJS = \
    $(MODELDIR)/modelparameters.o \
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

SCRIPTS  = \
    $(CUDIMOT_SRCDIR)/utils/Run_dtifit.sh \
    $(CUDIMOT_SRCDIR)/utils/jobs_wrapper.sh \
    $(CUDIMOT_SRCDIR)/utils/initialise_Bingham.sh \
    $(MODELDIR)/$(MODELNAME)_finish.sh \
    $(MODELDIR)/Pipeline_$(MODELNAME).sh

XFILES   = \
    cart2spherical \
    getFanningOrientation \
    initialise_Psi \
    $(MODELNAME) \
    cudimot_$(MODELNAME).sh \
    merge_parts_$(MODELNAME) \
    split_parts_$(MODELNAME) \
    testFunctions_$(MODELNAME)

DATAFILES = \
    $(MODELNAME)_priors \
    $(MODELNAME).info

all: $(XFILES) $(DATAFILES)
clean: cleandata

cleandata:
	rm -rf $(DATAFILES)

%.o: $(CUDIMOT_SRCDIR)/%.cu
	$(NVCC) -c $(NVCCFLAGS) $^

%.o: $(CUDIMOT_SRCDIR)/%.cc
	$(CXX) -c $(CXXFLAGS) $^

cart2spherical: $(CUDIMOT_SRCDIR)/utils/cart2spherical.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

getFanningOrientation: $(CUDIMOT_SRCDIR)/utils/getFanningOrientation.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

initialise_Psi: $(CUDIMOT_SRCDIR)/utils/initialise_Psi.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(MODELNAME): $(CUDIMOT_SRCDIR)/cudimot.cc cudimotoptions.o $(CUDIMOT_CUDA_OBJS)
	$(NVCC) $(NVCCFLAGS) $(NVCCLDFLAGS) -o $@ $^

merge_parts_$(MODELNAME): $(CUDIMOT_SRCDIR)/merge_parts.cc cudimotoptions.o $(CUDIMOT_CUDA_OBJS)
	$(NVCC) $(NVCCFLAGS) $(NVCCLDFLAGS) -o $@ $^ -lboost_filesystem -lboost_system

split_parts_$(MODELNAME): $(CUDIMOT_SRCDIR)/split_parts.cc cudimotoptions.o $(CUDIMOT_CUDA_OBJS)
	$(NVCC) $(NVCCFLAGS) $(NVCCLDFLAGS) -o $@ $^ -lboost_filesystem -lboost_system

testFunctions_$(MODELNAME): $(CUDIMOT_SRCDIR)/testFunctions.cu $(MODELDIR)/modelparameters.o
	$(NVCC) $(NVCCFLAGS) $(NVCCLDFLAGS) -o $@ $^

cudimot_$(MODELNAME).sh : $(MODELNAME)
	$(CUDIMOT_SRCDIR)/generate_wrapper.sh $(MODELNAME)

$(MODELNAME).info : cudimot_$(MODELNAME).sh
	$(CUDIMOT_SRCDIR)/generate_info.sh $(MODELNAME)

$(MODELNAME)_priors:
	cp $(MODELDIR)/modelpriors $@
