/*  merge_parts.cc

    Moises Hernandez  - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

#include <sys/stat.h>
#include "boost/filesystem.hpp"
#include "armawrap/newmat.h"
#include "newimage/newimageall.h"
#include "cudimotoptions.h"
#include "dMRI_Data.h"
#include "Model.h"
    
using namespace Cudimot;
using namespace boost::filesystem;
using namespace std;
using namespace NEWMAT;
using MISCMATHS::num2str;

void join_Parts(NEWIMAGE::volume<MyType> mask, string directory_in, string name_in, string name_out, int nsamples, int nParts, float max, float min){
    
  Matrix result(nsamples,0);
  Matrix part;

  for(int i=0;i<nParts;i++){
    
    std::string file_name;
    file_name.assign(directory_in);
    file_name += num2str(i);
    file_name += "/"; 
    file_name += name_in; 
    
    std::ifstream in;
    long nbytes_file;
    int nvox_file,nsamples_file;
    in.open(file_name.data(), ios::in | ios::binary);
    in.read((char*)&nvox_file, 4);
    in.read((char*)&nsamples_file, 4);
    if(nsamples==-1){
      // Do not know id advance the number od data measurements
      nsamples=nsamples_file;
      result.ReSize(nsamples,0);
    }
    in.read((char*)&nbytes_file, sizeof(long));
    if(nvox_file<=0 || nsamples_file<=0 || nsamples_file!=nsamples ){
      cerr << "CUDIMOT Error: The amount of data in the intermediate output file: " << file_name.data() << " is not correct." << endl;
      exit(-1);
    }
    part.ReSize(nsamples_file,nvox_file);
    in.read((char*)&part(1,1), nbytes_file);
    in.close();
    result |= part;
  }
  NEWIMAGE::volume4D<MyType> tmp;
  tmp.setmatrix(result,mask);
  if(max==-10) max=tmp.max();
  if(min==-10) min=tmp.min(); 
  tmp.setDisplayMaximumMinimum(max,min);
  save_volume4D(tmp,name_out);
}

//////////////////////////////////////////////////////////
//       MERGE THE OUTPUTS FILES OF CUDIMOT
//////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  // Setup logging:
  Log& logger = LogSingleton::getInstance();
  cudimotOptions& opts = cudimotOptions::getInstance();
  opts.parse_command_line(argc,argv,logger);

  NEWIMAGE::volume<MyType> mask;
  read_volume(mask,opts.maskfile.value());
  
  int nsamples=0;
  if(!opts.runMCMC.value()){
    nsamples=1; // LM
  }else{
    nsamples=opts.njumps.value()/opts.sampleevery.value();
  }
  if(nsamples<=0){
    cerr << "CUDIMOT Error: The number of samples must be greater than 0" << endl;
    exit (EXIT_FAILURE);
  }

  // get path of this binary to get the priors file
  char buf[1024];
  ssize_t count = readlink("/proc/self/exe",buf,sizeof(buf)-1);
  string bin_path(buf,(count > 0) ? count : 0 );
  string default_priors_file(bin_path.substr(0,bin_path.find_last_of("\\/")));
  string pattern("merge_parts_");
  default_priors_file+=("/"+bin_path.substr(bin_path.find(pattern)+pattern.size())+"_priors");
  
  Model<MyType> model(default_priors_file);
  
  int nparams = model.getNparams();

  string path_in;
  path_in.append(opts.partsdir.value());
  path_in.append("/part_");

  string path_out;
  path_out.append(opts.outputdir.value());

  for(int par=0;par<nparams;par++){
    string file_name = "Param_" + num2str(par) + "_samples";
    std::string output_file=path_out+"/"+file_name;
            
    join_Parts(mask,path_in,file_name,output_file,nsamples,opts.nParts.value(),-10,-10);
  }

  // If getPredictedSignal, join the different parts
  if(opts.getPredictedSignal.value()){
    string file_name = "PredictedSignal";
    std::string output_file=path_out+"/"+file_name;
    
    // it does not know the number of measurements. Set to -1 and it will get the number from the first part
    join_Parts(mask,path_in,file_name,output_file,-1,opts.nParts.value(),-10,-10);
  }

  // If BIC/AIC, join the different parts
  if(opts.BIC_AIC.value()){
    string file_name = "BIC";
    std::string output_file=path_out+"/"+file_name;
    join_Parts(mask,path_in,file_name,output_file,1,opts.nParts.value(),-10,-10);
    file_name = "AIC";
    output_file=path_out+"/"+file_name;
    join_Parts(mask,path_in,file_name,output_file,1,opts.nParts.value(),-10,-10);
    
  }

  // If Rician Noise, join tau samples
  if(opts.rician.value()&&opts.runMCMC.value()){
    string file_name = "Tau_samples";
    std::string output_file=path_out+"/"+file_name;
    
    join_Parts(mask,path_in,file_name,output_file,nsamples,opts.nParts.value(),-10,-10);
  }

  // Delete the temporal files
  if(!opts.keepTmp.value()){
    string path_remove;
    remove_all(opts.partsdir.value());
  }


  return 0;
}

