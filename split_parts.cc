/*  split_parts.cc

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

#include <sys/stat.h>
#include "boost/filesystem.hpp"
#include "newimage/newimageall.h"
#include "cudimotoptions.h"
#include "Model.h"

using namespace Cudimot;
using namespace boost::filesystem;

void save_part(Matrix data, string path, string name, int idpart){
  int nvox = data.Ncols();
  int nmeas = data.Nrows();
 
  string file_name;
  file_name = path+num2str(idpart)+"/"+name;
  
  ofstream out;
  out.open(file_name.data(), ios::out | ios::binary);
  out.write((char*)&nvox,4); // number of voxels
  out.write((char*)&nmeas,4); // number of measurements
  long size=nvox*nmeas*sizeof(Real);
  out.write((char*)&size,sizeof(long)); // number of bytes
  out.write((char*)&data(1,1),size); // information
  out.close(); 
}


int main(int argc, char *argv[]){
  Log& logger = LogSingleton::getInstance();
  cudimotOptions& opts = cudimotOptions::getInstance();
  opts.parse_command_line(argc,argv,logger);
  
   // Check if GridSearch, MCMC or LevMar flags
  if(opts.gridSearch.value()=="" && opts.no_LevMar.value() && !opts.runMCMC.value()){
    cerr << "CUDIMOT Error: You must select at least one method to fit the model: GridSearch, Levenberg_Marquardt or MCMC" << endl;
    exit(-1);
  }

  NEWIMAGE::volume4D<MyType> data;
  NEWIMAGE::volume<MyType> mask;
  read_volume4D(data,opts.datafile.value());
  read_volume(mask,opts.maskfile.value());
 
  Matrix dataM;
  dataM=data.matrix(mask);
  
  int nmeas = dataM.Nrows();
  if(nmeas<=0){
    cerr << "CUDIMOT Error: The number of diffusion-weighted measurements must be greater than 0 in the input volume\n" << endl;
    exit (EXIT_FAILURE);
  }
  
  int nvoxels=dataM.Ncols();
  if(nvoxels<=0){
    cerr << "CUDIMOT Error: The number of voxels must be greater than 0" << endl;
    exit (EXIT_FAILURE);
  }

  if(nvoxels<opts.nParts.value()){
    cerr << "CUDIMOT Error: The number of parts/jobs must be lower than number of voxels" << endl;
    exit (EXIT_FAILURE);
  }

  // Create directories for the different parts
  for(int i=0;i<(opts.nParts.value());i++){
    string dirpath=opts.partsdir.value()+"/part_"+num2str(i);

    create_directory(dirpath);
  }
  
  int size_part=nvoxels/opts.nParts.value();
  
  Matrix data_part;
  string out_path;
  string out_name("data");
  out_path.append(opts.partsdir.value());
  out_path.append("/part_");
  for(int i=0;i<(opts.nParts.value()-1);i++){
    data_part = dataM.SubMatrix(1,nmeas,i*size_part+1,(i+1)*size_part);
    save_part(data_part,out_path,out_name,i);
  }

  // last part
  data_part = dataM.SubMatrix(1,nmeas,(opts.nParts.value()-1)*size_part+1,nvoxels);
  save_part(data_part,out_path,out_name,(opts.nParts.value()-1));
  
  //////////////////////////////////////////////////////
  /// Initialization of parameters
  /// The user can provide nifti files for some parameters
  /// Divide into different parts
  //////////////////////////////////////////////////////
  char buf[1024]; // get path of this binary to get the priors file
  ssize_t count = readlink("/proc/self/exe",buf,sizeof(buf)-1);
  string bin_path(buf,(count > 0) ? count : 0 );
  string default_priors_file(bin_path.substr(0,bin_path.find_last_of("\\/")));
  string pattern("split_parts_");
  default_priors_file+=("/"+bin_path.substr(bin_path.find(pattern)+pattern.size())+"_priors");
  
  Model<MyType> model(default_priors_file);

  int nparams=model.getNparams();
  if (opts.init_params.set()){
    string filename(opts.init_params.value());
    ifstream file(filename.data());
    if (file.is_open()){
      string line;
      int id_param=-1;
      while(getline(file,line)){
	id_param++;
	if(id_param>=(nparams)){
	  cerr << "CUDIMOT Error: Too many lines for specifying the initialization of the parameters in the file: " << filename << ". The number of lines must match the number of parameters in this model: " << nparams << endl;
	  exit (EXIT_FAILURE);
	}

	if (!line.empty()){
	  // Read volume with values fot this parameter
	  string name_file(line);
	  NEWIMAGE::volume4D<MyType> param_vals;
	  read_volume4D(param_vals,name_file);

	  if(mask.xsize()!=param_vals.xsize() || mask.ysize()!=param_vals.ysize() || mask.zsize()!=param_vals.zsize()){
	    cerr << "CUDIMOT Error: The size of the mask and the volume used for initilizing the parameters: " << name_file << " does not match\n" << endl;
	    exit (EXIT_FAILURE);
	  }

	  Matrix paramM;
	  paramM=param_vals.matrix(mask);
	  
	  int size4dim = paramM.Nrows();
	  if(size4dim!=1){
	    cerr << "CUDIMOT Error: The volume used for initilize the parameters: " << name_file << " must be a 3D volume \n" << endl;
	    exit (EXIT_FAILURE);
	  }
	  
	  if(nvoxels!=paramM.Ncols()){
	    cerr << "CUDIMOT Error: The number of voxels in the data and the volume used for initilizing the parameters: " << name_file << " does not match\n" << endl;
	    exit (EXIT_FAILURE);
	  }
	  
	  Matrix param_part;
	  string out_name_p;
	  out_name_p.append("ParamInit_");
	  out_name_p.append(num2str(id_param));
	  
	  for(int i=0;i<(opts.nParts.value()-1);i++){
	    param_part = paramM.SubMatrix(1,1,i*size_part+1,(i+1)*size_part);
	    save_part(param_part,out_path,out_name_p,i);
	  }
	  // last part
	  param_part = paramM.SubMatrix(1,1,(opts.nParts.value()-1)*size_part+1,nvoxels);
	  save_part(param_part,out_path,out_name_p,(opts.nParts.value()-1));
	  
	}else{
	  // Empty line, initialise this parameter with default value if provided or zeros otherwise
	  Matrix param_part;
	  string out_name_p;
	  out_name_p.append("ParamInit_");
	  out_name_p.append(num2str(id_param));
	  
	  param_part.ReSize(1,size_part);
	  for(int i=0;i<(opts.nParts.value()-1);i++){
	    if(model.initProvided()){
	      param_part = model.getParam_init(id_param);
	    }else{
	      param_part = 0;
	    }
	    save_part(param_part,out_path,out_name_p,i);
	  }
	  // last part
	  int size_last_part=nvoxels-((opts.nParts.value()-1)*size_part);
	  param_part.ReSize(1,size_last_part); 
	  if(model.initProvided()){
	    param_part = model.getParam_init(id_param);
	  }else{
	    param_part = 0;
	  }
	  save_part(param_part,out_path,out_name_p,(opts.nParts.value()-1));
	}
	
      } //end lines
      
      if(id_param!=(nparams-1)){
	cerr << "CUDIMOT Error: The number of volumes (lines) provided for initializing the parameters in: " << filename.data() << " does not match the number of parameters of this model: " << nparams << ". If a parameter does not need initialization, its line can be empty.\n" << endl;
	exit(-1);
      }      
    }else{
      cerr << "CUDIMOT Error: Unable to open Initialization Parameter file: " << filename.data() << endl; 
      exit(-1);
    }
  }else{
    // Not initialization file provided. Initialise with default parameter values if provided or zeros otherwise. But then not file division is needed. So not need to do anything here.
  }

  //////////////////////////////////////////////////////
  /// Fixed parameters - FP
  /// The user can provide nifti files
  /// Divide into different parts
  //////////////////////////////////////////////////////
  int nFixP = model.getNFixP();
  if (nFixP>0 && !opts.FixP.set()){
    cerr << "CUDIMOT Error: Expected a list with " << nFixP << " NIfTI files for specifying the Fixed Parameters of the model: Use option --FixP"<< endl; 
    exit(-1);
  }
  
  if (nFixP){
    string filename(opts.FixP.value());
    ifstream file(filename.data());
    if (file.is_open()){
      string line;
      int id_FP=-1;
      while(getline(file,line) && id_FP<(nFixP-1)){
	id_FP++;
	if (!line.empty()){
	  // Read volume with values for this parameter
	  string name_file(line);
	  NEWIMAGE::volume4D<MyType> fixedParam;
	  read_volume4D(fixedParam,name_file);
	  
	  if(mask.xsize()!=fixedParam.xsize() || mask.ysize()!=fixedParam.ysize() || mask.zsize()!=fixedParam.zsize()){
	    cerr << "CUDIMOT Error: The size of the mask and the volume used for specifying the Fixed Parameters: " << name_file << " does not match" << endl;
	    exit (EXIT_FAILURE);
	  }

	  Matrix fixedParamM;
	  fixedParamM=fixedParam.matrix(mask);
	  
	  int size4dim = fixedParamM.Nrows();
	  if(size4dim!=model.getNFixP_size(id_FP)){
	    cerr << "CUDIMOT Error: The NIfTI file used for specifying the Fixed Parameters: " << name_file << " must have " << model.getNFixP_size(id_FP) << " volumes" << endl;
	    exit (EXIT_FAILURE);
	  }
	  
	  if(nvoxels!=fixedParamM.Ncols()){
	    cerr << "CUDIMOT Error: The number of voxels in the data and the volume used for specifying the Fixed Parameters: " << name_file << " does not match\n" << endl;
	    exit (EXIT_FAILURE);
	  }
	  
	  Matrix FP_part;
	  string out_name_p;
	  out_name_p.append("FixParam_");
	  out_name_p.append(num2str(id_FP));
	  
	  for(int i=0;i<(opts.nParts.value()-1);i++){
	    FP_part = fixedParamM.SubMatrix(1,size4dim,i*size_part+1,(i+1)*size_part);
	    save_part(FP_part,out_path,out_name_p,i);
	  }
	  // last part
	  FP_part = fixedParamM.SubMatrix(1,size4dim,(opts.nParts.value()-1)*size_part+1,nvoxels);
	  save_part(FP_part,out_path,out_name_p,(opts.nParts.value()-1));
	  
	}else{
	  // Empty line, initialise this parameter with default value if provided or zeros otherwise
	  cerr << "CUDIMOT Error: Please remove the empty lines in the file used for specifying the Fixed Parameters: " << filename.data() << endl; 
	  exit(-1);
	}
      } //end lines
      if(id_FP!=(nFixP-1)){
	cerr << "CUDIMOT Error: The number of NIfTI files (lines) provided in the file for specifying the Fixed Parameters: " << filename.data() << " does not match the number of parameters of this model: " << nFixP << endl; 
	exit(-1);
      }      
    }else{
      cerr << "CUDIMOT Error: Unable to open Fixed Parameter file: " << filename.data() << endl; 
      exit(-1);
    }
  }else{
    if(opts.FixP.set()){
      cerr << "CUDIMOT Error: This model does not need any Fixed Parameter, however the user has provided a Fixed Parameter file."  << endl; 
      exit(-1);
    }
  }
    
}
