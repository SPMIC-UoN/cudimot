/*  cudimotOptions.h
    
    Moises Hernandez-Fernandez FMRIB Image Analysis Group
    
    Copyright (C) 1999-2010 University of Oxford  */

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

#if !defined(cudimotOptions_h)
#define cudimotOptions_h

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"
#include "utils/tracer_plus.h"

using namespace Utilities;

namespace Cudimot {
  
  class cudimotOptions {
  public:
    static cudimotOptions& getInstance();
    ~cudimotOptions() { delete gopt; }
    
    Option<bool> help;
    Option<std::string> logdir;
    Option<bool> forcedir;
    Option<std::string> datafile;
    Option<std::string> maskfile;
    Option<std::string> partsdir;
    Option<std::string> outputdir;
    Option<int> idPart;
    Option<int> nParts;
    Option<int> njumps;
    Option<int> nburn;
    Option<int> sampleevery;
    Option<int> updateproposalevery;
    Option<bool> no_updateproposal;
    Option<int> seed;
    Option<bool> no_LevMar;
    Option<bool> no_Marquardt;
    Option<int> iterLevMar;
    Option<std::string> gridSearch;
    Option<bool> runMCMC;
    Option<bool> rician;
    Option<bool> keepTmp;
    Option<bool> getPredictedSignal;
    Option<std::string> CFP;
    Option<std::string> FixP;
    Option<std::string> fixed;
    Option<std::string> init_params;
    Option<std::string> debug;
    Option<bool> BIC_AIC;
    FmribOption<std::string> priorsfile;
    
    void parse_command_line(int argc, char** argv,  Log& logger);
  
  private:
    cudimotOptions();  
    const cudimotOptions& operator=(cudimotOptions&);
    cudimotOptions(cudimotOptions&);
    
    OptionParser options; 
    
    static cudimotOptions* gopt;
  };
  
  inline cudimotOptions& cudimotOptions::getInstance(){
    if(gopt == NULL)
      gopt = new cudimotOptions();
    
    return *gopt;
  }
  
  inline cudimotOptions::cudimotOptions():
	help(std::string("-h,--help"), false,
		std::string("\tDisplay this message"),
		false, no_argument),
	logdir(std::string("--ld,--logdir"), std::string("logdir"),
		std::string("\tLog directory (default is logdir)"),
		false, requires_argument),
	forcedir(std::string("--forcedir"),false,
		std::string("\tUse the actual directory name given - i.e. don't add + to make a new directory"),
		false,no_argument),
	datafile(std::string("--data"), std::string("data"),
		std::string("\t\tData file"),
		true, requires_argument),  
	maskfile(std::string("--maskfile"), std::string("nodif_brain_mask"),
		std::string("\tMask file"),
		true, requires_argument),
	partsdir(std::string("--partsdir"), std::string(""),
		std::string("\tDirectory where different parts of the data/results will be stored"),
		true, requires_argument),
	outputdir(std::string("--outputdir"), std::string(""),
		std::string("\tDirectory where to write the output files"),
		true, requires_argument),
	idPart(std::string("--idPart"),0,
		std::string("\tNumber of the part of the dataset to process [0..N-1]"),
		true, requires_argument),
	nParts(std::string("--nParts"),0,
		std::string("\tTotal number of parts of the dataset [1..N]"),
		true, requires_argument),
	njumps(std::string("--nj,--njumps"),1250,
		std::string("\tNum of jumps to be made by MCMC (default is 1250)"),
		false,requires_argument),
	nburn(std::string("--bi,--burnin"),1000,
		std::string("\tTotal num of jumps at start of MCMC to be discarded (default is 1000)"),
		false,requires_argument),
	sampleevery(std::string("--se,--sampleevery"),25,
		std::string("\tNum of jumps for each sample (MCMC) (default is 25)"),
		false,requires_argument),
	updateproposalevery(std::string("--upe,--updateproposalevery"),40,
		std::string("Num of jumps for each update to the proposal density std (MCMC) (default is 40)"),
		false,requires_argument),
	no_updateproposal(std::string("--no_updateproposal"),false,
		std::string("Do not update the proposal density std during the recording step of MCMC"),
		false,no_argument),
	seed(std::string("--seed"),8219,
		std::string("\t\tSeed for pseudo random number generator"),
		false,requires_argument),
	no_LevMar(std::string("--no_LevMar"),false,
		std::string("\tDo not run Levenberg-Marquardt algorithm"),
		false, no_argument),
	no_Marquardt(std::string("--no_Marquardt"),false,
		std::string("\tDo not use Marquardt contribution in Levenberg algorithm"),
		false, no_argument),
	iterLevMar(std::string("--iterLevMar"),200,
		std::string("\tMaximum number of iterations in Levenberg(-Marquardt) algorithm (default is 200)"),
		false,requires_argument),
	gridSearch(std::string("--gridSearch"),std::string(""),
		std::string("\tRun gridSearch algorithm using the values specified in a file"),
		false,requires_argument),
   	runMCMC(std::string("--runMCMC"),false,
		std::string("\tRun MCMC algorithm and get distribution of estimates"),
		false,no_argument),
	rician(std::string("--rician"),false,
		std::string("\tUse Rician noise modelling in MCMC"),
		false,no_argument),
        keepTmp(std::string("--keepTmp"),false,
		std::string("\tDo not remove the temporal directory created for storing the data/results parts"),
		false,no_argument),
	getPredictedSignal(std::string("--getPredictedSignal"),false,
		std::string("Save the predicted signal by the model at the end"),
		false,no_argument),
	CFP(std::string("--CFP"), std::string(""),
		std::string("\t\tFile with a list of ASCCI files for specifying the common (to all voxels) fixed parameters of the model"),
		false, requires_argument),  
	FixP(std::string("--FixP"), std::string(""),
		std::string("\t\tFile with a list of NIfTI files for specifying the fixed parameters of the model"),
		false, requires_argument), 
	fixed(std::string("--fixed"), std::string(""),
		std::string("\t\tList of the Ids of the fixed parameters separated by commas. For example: --fixed=2,4"),
		false, requires_argument), 
	init_params(std::string("--init_params"), std::string(""),
		std::string("\tFile with a list of NIfTI files for the initialization of the model parameters"),
		false, requires_argument),
        debug(std::string("--debug"), std::string(""),
		std::string("\t\tSpecify a voxel for debugging. Some variables at certain steps of the algorithms will be printed (use few iterations)"),
		false, requires_argument),
	BIC_AIC(std::string("--BIC_AIC"), false,
	        std::string("\tCalculate Bayesian and Akaike Information Criteria at the end"),
		false, no_argument),
	priorsfile(std::string("--priors"), std::string(""),
		std::string("\tFile with parameters information (initialization, bounds and priors)"),
		false, requires_argument),

   options("CUDIMOT","YourModelName --help (for list of options)\n")
     {
       try {
	options.add(help);
	options.add(logdir);
	options.add(forcedir);
	options.add(datafile);
	options.add(maskfile);
	options.add(partsdir);
	options.add(outputdir);
	options.add(idPart);
	options.add(nParts);
	options.add(njumps);
	options.add(nburn);
	options.add(sampleevery);
	options.add(updateproposalevery);
	options.add(no_updateproposal);
	options.add(seed);
	options.add(gridSearch); 
	options.add(runMCMC); 
	options.add(no_LevMar);
	options.add(no_Marquardt);
	options.add(iterLevMar);
	options.add(rician);
	options.add(keepTmp);
	options.add(getPredictedSignal);
	options.add(CFP);
	options.add(FixP);
	options.add(fixed);
	options.add(init_params);
	options.add(debug);
	options.add(BIC_AIC);
	options.add(priorsfile);
     }
     catch(X_OptionError& e) {
       options.usage();
       std::cerr << std::endl << e.what() << std::endl;
     } 
     catch(std::exception &e) {
       std::cerr << e.what() << std::endl;
     }    
     
   }
}

#endif





