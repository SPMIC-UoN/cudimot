include $(FSLCONFDIR)/default.mk

ifndef modelname
$(error Error: modelname has not been set)
endif
ifndef FSLDEVDIR
$(error Error: FSLDEVDIR has not been set)
endif

MAX_REGISTERS= -maxrregcount 64
TYPE := $(shell cat mymodels/${modelname}/modelparameters.h | grep MyType | cut -f 2 -d ' ')
ifeq ($(TYPE),float)
  MAX_REGISTERS= -maxrregcount 64
else
  MAX_REGISTERS= -maxrregcount 128
endif

MODELPATH= mymodels/$(modelname)

PROJNAME = CUDIMOT

# Name of main executable
CUDIMOT  = $(DIR_objs)/${modelname}


DIR_objs     = objs
USRINCFLAGS  = -I$(MODELPATH)
USRNVCCFLAGS =  $(MAX_REGISTERS)

LIBS     = -lfsl-warpfns -lfsl-basisfield -lfsl-meshclass -lfsl-newimage \
           -lfsl-utils -lfsl-miscmaths -lfsl-newran -lfsl-NewNifti \
           -lfsl-znz -lfsl-cprob -lboost_filesystem -lboost_system
CUDALIBS = -lcurand

OBJS    := modelparameters.o init_gpu.o dMRI_Data.o Model.o Parameters.o \
           GridSearch.o Levenberg_Marquardt.o MCMC.o BIC_AIC.o \
           getPredictedSignal.o cudimotoptions.o

SCRIPTS  = ${modelname}@info ${MODELPATH}/Pipeline_${modelname}.sh \
           ${MODELPATH}/${modelname}_finish.sh utils/Run_dtifit.sh \
           utils/jobs_wrapper.sh utils/initialise_Bingham.sh
XFILES  := cart2spherical getFanningOrientation initialise_Psi \
           split_parts_${modelname} ${modelname} \
           merge_parts_${modelname} testFunctions_${modelname} \
           cudimot_${modelname}.sh ${modelname}_priors $(modelname).info

XFILES := $(addprefix $(DIR_objs)/, $(XFILES))
OBJS   := $(addprefix $(DIR_objs)/, $(OBJS))

cleanall:
	rm -f $(DIR_objs)/*.o
	rm -f $(DIR_objs)/testFunctions_${modelname}
cleanbin:
	rm $(DIR_objs)/*
debugging:
	rm -f $(DIR_objs)/testFunctions_${modelname}
	rm -f $(DIR_objs)/Levenberg_Marquardt.o
	rm -f $(DIR_objs)/MCMC.o
	rm -f $(DIR_objs)/GridSearch.o
	make install
makedir:
	mkdir -p $(FSLDEVDIR)/bin
	mkdir -p $(DIR_objs)

all: makedir ${XFILES} ${DATAFILES}

$(DIR_objs)/cart2spherical:
	${CXX} -o $@ utils/cart2spherical.cc ${CXXFLAGS} ${LDFLAGS}
$(DIR_objs)/getFanningOrientation:
	${CXX} -o $@ utils/getFanningOrientation.cc ${CXXFLAGS} ${LDFLAGS}
$(DIR_objs)/initialise_Psi:
	${CXX} -o $@ utils/initialise_Psi.cc ${CXXFLAGS} ${LDFLAGS}

$(DIR_objs)/cudimotoptions.o:
	${CXX} ${CXXFLAGS} -c -o $@ cudimotoptions.cc

$(DIR_objs)/split_parts_${modelname}: ${OBJS}
	${NVCC} ${NVCCFLAGS} -o $@ split_parts.cc $^ ${NVCCLDFLAGS}

$(DIR_objs)/merge_parts_${modelname}: ${OBJS}
	${NVCC} ${NVCCFLAGS} -o $@ merge_parts.cc $^ ${NVCCLDFLAGS}

$(DIR_objs)/%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c -o $@ $<

$(DIR_objs)/%.o: ${MODELPATH}/%.cc
	$(NVCC) $(NVCCFLAGS) -c -o $@ $<

${CUDIMOT}:	$(DIR_objs)/cudimot.o ${OBJS}
	${NVCC} ${NVCCFLAGS} -o $@ $^ ${NVCCLDFLAGS}

$(DIR_objs)/testFunctions_${modelname}:
	$(NVCC) ${NVCCFLAGS} -o $@ $(MODELPATH)/modelparameters.cc testFunctions.cu ${NVCCLDFLAGS}

$(DIR_objs)/cudimot_${modelname}.sh: $(DIR_objs)/${modelname}
	./generate_wrapper.sh ${modelname}
	mv cudimot_${modelname}.sh $(DIR_objs)

$(DIR_objs)/${modelname}.info : $(DIR_objs)/cudimot_${modelname}.sh
	./generate_info.sh ${modelname} ${MODELPATH}
	mv ${modelname}.info $(DIR_objs)

$(DIR_objs)/${modelname}_priors:
	cp ${MODELPATH}/modelpriors $@
