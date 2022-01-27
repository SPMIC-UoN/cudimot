include $(FSLCONFDIR)/default.mk

ifndef modelname
$(error Error: modelname has not been set)
endif
ifndef FSLDEVDIR
$(error Error: FSLDEVDIR has not been set)
endif

NVCC = ${CUDA}/bin/nvcc
DIR_objs=./objs

MAX_REGISTERS= -maxrregcount 64
TYPE := $(shell cat mymodels/${modelname}/modelparameters.h | grep MyType | cut -f 2 -d ' ')
ifeq ($(TYPE),float)
  MAX_REGISTERS= -maxrregcount 64
else
  MAX_REGISTERS= -maxrregcount 128
endif

MODELPATH= mymodels/$(modelname)

NVCC_FLAGS = -I$(MODELPATH) -O3 -dc $(MAX_REGISTERS)
#-Xptxas -v 
#-G -lineinfo

CUDA_INC=-I${CUDA}/lib -I${CUDA}/lib64

CUDA_INC = -I. -I${FSLDIR}/extras/include/newmat -I${FSLDIR}/include
SM_20 = -gencode arch=compute_20,code=sm_20
SM_21 = -gencode arch=compute_20,code=sm_21
SM_30 = -gencode arch=compute_30,code=sm_30
SM_35 = -gencode arch=compute_35,code=sm_35
SM_37 = -gencode arch=compute_37,code=sm_37
SM_50 = -gencode arch=compute_50,code=sm_50
SM_52 = -gencode arch=compute_52,code=sm_52
SM_60 = -gencode arch=compute_60,code=sm_60
SM_61 = -gencode arch=compute_61,code=sm_61
SM_70 = -gencode arch=compute_70,code=sm_70

#for Realease
GPU_CARDs = $(SM_30) $(SM_35) $(SM_37) $(SM_50) $(SM_52) 
#$(SM_60) $(SM_61)

#for FMRIB
#GPU_CARDs = $(SM_37) $(SM_35)

PROJNAME = CUDIMOT

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_NEWRAN} -I${INC_CPROB} -I${INC_PROB} -I${INC_BOOST} -I${INC_ZLIB} -I$(MODELPATH) 
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_NEWRAN} -L${LIB_CPROB} -L${LIB_PROB} -L${LIB_ZLIB}

DLIBS = -lwarpfns -lbasisfield -lmeshclass -lnewimage -lutils -lmiscmaths -lnewmat -lnewran -lNewNifti -lznz -lcprob -lprob -lm -lz

CUDIMOT=$(DIR_objs)/${modelname}

CUDIMOT_CUDA_OBJS=$(DIR_objs)/modelparameters.o $(DIR_objs)/init_gpu.o $(DIR_objs)/dMRI_Data.o $(DIR_objs)/Model.o $(DIR_objs)/Parameters.o $(DIR_objs)/GridSearch.o $(DIR_objs)/Levenberg_Marquardt.o $(DIR_objs)/MCMC.o $(DIR_objs)/BIC_AIC.o $(DIR_objs)/getPredictedSignal.o

CUDIMOT_OBJS=$(DIR_objs)/link_cudimot_gpu.o $(DIR_objs)/cudimot.o $(DIR_objs)/cudimotoptions.o

SGEBEDPOST = bedpost
SGEBEDPOSTX = bedpostx bedpostx_postproc.sh bedpostx_preproc.sh bedpostx_single_slice.sh bedpostx_datacheck

SCRIPTS = ${modelname}@info utils/Run_dtifit.sh utils/jobs_wrapper.sh utils/initialise_Bingham.sh
FILES = cart2spherical getFanningOrientation initialise_Psi split_parts_${modelname} ${modelname} merge_parts_${modelname} testFunctions_${modelname}
XFILES=$(addprefix $(DIR_objs)/, $(FILES))

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

all: 	cleanall makedir ${XFILES}

$(DIR_objs)/cart2spherical: 
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ utils/cart2spherical.cc ${DLIBS} 
$(DIR_objs)/getFanningOrientation: 
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ utils/getFanningOrientation.cc ${DLIBS} 
$(DIR_objs)/initialise_Psi: 
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ utils/initialise_Psi.cc ${DLIBS} 

$(DIR_objs)/cudimotoptions.o:
	${CXX} ${CXXFLAGS} ${LDFLAGS} -c -o $@ cudimotoptions.cc ${DLIBS} 

$(DIR_objs)/split_parts_${modelname}: $(DIR_objs)/cudimotoptions.o $(DIR_objs)/link_cudimot_gpu.o
	${CXX} ${CXXFLAGS} $(USRINCFLAGS) ${LDFLAGS} -o $@ $(DIR_objs)/cudimotoptions.o $(DIR_objs)/link_cudimot_gpu.o split_parts.cc ${DLIBS} $(CUDIMOT_CUDA_OBJS) -lcudart -lboost_filesystem -lboost_system -L${CUDA}/lib64 -L${CUDA}/lib

$(DIR_objs)/merge_parts_${modelname}: $(DIR_objs)/cudimotoptions.o $(DIR_objs)/link_cudimot_gpu.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ $(DIR_objs)/cudimotoptions.o $(DIR_objs)/link_cudimot_gpu.o merge_parts.cc ${DLIBS} $(CUDIMOT_CUDA_OBJS) -lcudart -lboost_filesystem -lboost_system -L${CUDA}/lib64 -L${CUDA}/lib

$(DIR_objs)/init_gpu.o: 
		$(NVCC) $(GPU_CARDs) $(NVCC_FLAGS) -o $@ init_gpu.cu $(CUDA_INC)

$(DIR_objs)/dMRI_Data.o: 
		$(NVCC) $(GPU_CARDs) $(USRINCFLAGS) $(NVCC_FLAGS) -o $@ dMRI_Data.cu $(CUDA_INC)

$(DIR_objs)/Parameters.o: 
		$(NVCC) $(GPU_CARDs) $(USRINCFLAGS) $(NVCC_FLAGS) -o $@ Parameters.cu $(CUDA_INC)

$(DIR_objs)/Model.o: 
		$(NVCC) $(GPU_CARDs) $(USRINCFLAGS) $(NVCC_FLAGS) -o $@ Model.cu $(CUDA_INC)

$(DIR_objs)/modelparameters.o: 
		$(NVCC) $(GPU_CARDs) $(NVCC_FLAGS) $(MODELPATH)/modelparameters.cc -o $(DIR_objs)/modelparameters.o $(CUDA_INC)

$(DIR_objs)/GridSearch.o: 	
		$(NVCC) $(GPU_CARDs) $(USRINCFLAGS) $(NVCC_FLAGS) -o $@ GridSearch.cu $(CUDA_INC)

$(DIR_objs)/Levenberg_Marquardt.o: 	
		$(NVCC) $(GPU_CARDs) $(USRINCFLAGS) $(NVCC_FLAGS) -o $@ Levenberg_Marquardt.cu $(CUDA_INC)

$(DIR_objs)/MCMC.o: 	
		$(NVCC) $(GPU_CARDs) $(USRINCFLAGS) $(NVCC_FLAGS) -o $@ MCMC.cu $(CUDA_INC)

$(DIR_objs)/BIC_AIC.o: 	
		$(NVCC) $(GPU_CARDs) $(USRINCFLAGS) $(NVCC_FLAGS) -o $@ BIC_AIC.cu $(CUDA_INC)

$(DIR_objs)/getPredictedSignal.o: 	
		$(NVCC) $(GPU_CARDs) $(NVCC_FLAGS) -o $@ getPredictedSignal.cu $(CUDA_INC)

$(DIR_objs)/link_cudimot_gpu.o:	$(CUDIMOT_CUDA_OBJS)
		$(NVCC) $(GPU_CARDs) -dlink $(CUDIMOT_CUDA_OBJS) -o $@ -L${CUDA}/lib64 -L${CUDA}/lib

$(DIR_objs)/cudimot.o:
		$(NVCC) $(GPU_CARDs) $(USRINCFLAGS) $(NVCC_FLAGS) -o $@ cudimot.cc $(CUDA_INC)

${CUDIMOT}:	${CUDIMOT_OBJS}
		${CXX} ${CXXFLAGS} ${LDFLAGS} -o $(DIR_objs)/${modelname} ${CUDIMOT_OBJS} $(CUDIMOT_CUDA_OBJS) ${DLIBS} -lcudart -L${CUDA}/lib64 -L${CUDA}/lib
		./generate_wrapper.sh

$(DIR_objs)/testFunctions_${modelname}: 
	$(NVCC) $(GPU_CARDs) -I$(MODELPATH) -O3 $(MAX_REGISTERS) $(MODELPATH)/modelparameters.cc testFunctions.cu -o $(DIR_objs)/testFunctions_${modelname} $(CUDA_INC)

