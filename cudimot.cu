/*  cudimot.cc

    Moises Hernandez-Fernandez - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

#include <sys/time.h>
//#include <boost/filesystem/operations.hpp>
//#include <boost/filesystem/path.hpp>
#include "cudimotoptions.h"
#include "init_gpu.h"
#include "dMRI_Data.h"
#include "Model.h"
#include "Parameters.h"
#include "GridSearch.h"
#include "Levenberg_Marquardt.h"
#include "MCMC.h"

using namespace Cudimot;

double timeval_diff(struct timeval *a, struct timeval *b){
  return (double)(a->tv_sec +(double)a->tv_usec/1000000) - (double)(b->tv_sec +(double)b->tv_usec/1000000);
}

int main(int argc, char *argv[]){
  struct timeval t1,t2;
  double time;
  gettimeofday(&t1,NULL); 
  
  // Setup logging:
  Log& logger = LogSingleton::getInstance();
  cudimotOptions& opts = cudimotOptions::getInstance();
  opts.parse_command_line(argc,argv,logger);
  srand(opts.seed.value());  //randoms seed
  
  init_gpu();

  // Check if GridSearch, MCMC or LevMar flags
  if(opts.gridSearch.value()=="" && opts.no_LevMar.value() && !opts.runMCMC.value()){
    cerr << "CUDIMOT Error: You must select at least one method to fit the model: GridSearch, Levenberg_Marquardt or MCMC" << endl;
    exit(-1);
  }

  // Encapsulate dMRI data
  dMRI_Data<MyType> data;

  // get path of this binary to get the priors file
  char buf[1024];
  ssize_t count = readlink("/proc/self/exe",buf,sizeof(buf)-1);
  string bin_path(buf,(count > 0) ? count : 0 );
  string default_priors_file(bin_path+"_priors");
  
  Model<MyType> model(default_priors_file);
  
  Parameters<MyType> params(model,data);

  GridSearch<MyType> methodGridSearch(model.getNGridParams(),
  				      model.getGridParams(),
  				      model.getGridCombs(),
  				      model.getGrid(),
				      model.getBound_types(),
				      model.getBounds_min(),
				      model.getBounds_max());

  Levenberg_Marquardt<MyType> methodLM(model.getBound_types(),
				       model.getBounds_min(),
				       model.getBounds_max(),
				       model.getFixed());

  MCMC<MyType> methodMCMC(data.getNvoxFit_part(),
			  model.getBound_types(),
			  model.getBounds_min(),
			  model.getBounds_max(),
			  model.getPrior_types(),
			  model.getPriors_a(),
			  model.getPriors_b(),
			  model.getFixed());

  // Data and parameters are divided into parts => process each part
  for(int part=0;part<data.getNparts();part++){
    
    int part_size=0;    
    MyType* meas=data.getMeasPart(part,part_size);
    MyType* parameters_part = params.getParametersPart(part);
    
    if(opts.gridSearch.value()!=""){
      methodGridSearch.run(part_size,data.getNmeas(),
			   params.getTsize_CFP(),
			   params.getTsize_FixP(),
			   meas,parameters_part,
			   params.getCFP(),
			   params.getFixP_part(part));
    }

    if(!opts.no_LevMar.value()){
      methodLM.run(part_size,data.getNmeas(),
		   params.getTsize_CFP(),
		   params.getTsize_FixP(),
		   meas,parameters_part,
		   params.getCFP(),
		   params.getFixP_part(part));
    }
    
    if(!opts.runMCMC.value()){
      params.copyParamsPartGPU2Host(part);
      params.calculate_predictedSignal_BIC_AIC(0,part,meas);
    }else{
      methodMCMC.run(part_size,data.getNmeas(),
		     params.getTsize_CFP(),
		     params.getTsize_FixP(),
		     meas,parameters_part,
		     params.getCFP(),
		     params.getFixP_part(part),
		     params.getSamples(),
		     params.getTauSamples());
      
      params.copyParamsPartGPU2Host(part);
      params.copySamplesPartGPU2Host(part);
      params.calculate_predictedSignal_BIC_AIC(1,part,meas);
    }
  }
  if(!opts.runMCMC.value()){
    // only 1 sample, the value of parameters
    params.copyParams2Samples();
  }
  params.writeSamples();
  
  gettimeofday(&t2,NULL);
  time=timeval_diff(&t2,&t1);
  cout << endl << "Part processed in: " << time << " seconds" << endl;	
  
}

