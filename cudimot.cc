/*  cudimot.cc

    Moises Hernandez-Fernandez - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk

    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK


    LICENCE

    FMRIB Software Library, Release 6.0 (c) 2018, The University of
    Oxford (the "Software")

    The Software remains the property of the Oxford University Innovation
    ("the University").

    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.

    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.

    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.

    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. If you are
    interested in using the Software commercially, please contact Oxford
    University Innovation ("OUI"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    fsl@innovation.ox.ac.uk quoting Reference Project 9564, FSL.*/


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

