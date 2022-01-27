/* Parameters.cu
   
   Moises Hernandez-Fernandez - FMRIB Image Analysis Group
   
   Copyright (C) 2005 University of Oxford */

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

#include "Parameters.h"

namespace Cudimot{
  
  template <typename T>
  Parameters<T>::Parameters(Model<T> model,dMRI_Data<T> dMRI_data):
    nparams(model.nparams),
    nFixP(model.nFixP), FixP_Tsize(model.FixP_Tsize),
    nCFP(model.nCFP),CFP_Tsize(model.CFP_Tsize),
    nvox(dMRI_data.nvox), nmeas(dMRI_data.nmeas),nparts(dMRI_data.nparts),
    size_part(dMRI_data.size_part),size_last_part(dMRI_data.size_last_part),
    nvoxFit_part(dMRI_data.nvoxFit_part),bic_aic(dMRI_data.nvoxFit_part)
  {
    
    Log& logger = LogSingleton::getInstance();
    cudimotOptions& opts = cudimotOptions::getInstance();
    
    //////////////////////////////////////////////////////
    /// Initialise parameters
    /// The user can provide nifti files for some parameters
    //////////////////////////////////////////////////////    
    params_host=new T[nvox*nparams];
    if(opts.getPredictedSignal.value()){
      predSignal_host=new T[nvox*nmeas];
    }
    if(opts.BIC_AIC.value()){
      BIC_host=new T[nvox];
      AIC_host=new T[nvox];
    }
    
    if (opts.init_params.set()){
      // If volumes provided for initialize parameters
      for(int idParam=0;idParam<nparams;idParam++){
	// Read binary file with values for this part, logfile
	string name_file;
	name_file.append(opts.partsdir.value());
	name_file.append("/part_");
	name_file.append(num2str(opts.idPart.value()));
	name_file.append("/ParamInit_");
	name_file.append(num2str(idParam));
	
	ifstream in;
	long nbytes;
	int nvox_file,nmeas_file;
	in.open(name_file.data(), ios::in | ios::binary);
	in.read((char*)&nvox_file, 4);
	in.read((char*)&nmeas_file, 4);
	in.read((char*)&nbytes, sizeof(long));
	
	if(nvox!=nvox_file || nmeas_file!=1){
	  cerr << "CUDIMOT Error: The amount of data in the input file " <<  name_file << " for initializing the parameters is not correct" << endl;
	  exit(-1);
	}
	
	Matrix Parameters_init;
	Parameters_init.ReSize(1,nvox);
	in.read((char*)&Parameters_init(1,1),nbytes);
	in.close();
	
	for(int i=0;i<nvox;i++){
	  params_host[i*nparams+idParam]=Parameters_init(1,i+1);
	}
      }
      
    }else{
      // If NOT volumes provided for initialise parameters, try default values
      for(int idParam=0;idParam<nparams;idParam++){
	for(int i=0;i<nvox;i++){
	  if(model.initProvided()){
	    params_host[i*nparams+idParam]=model.getParam_init(idParam);
	  }else{
	    params_host[i*nparams+idParam]=0;
	  }
	}
      }
    }
    
    //////////////////////////////////////////////////////
    /// Read Common Fixed Parameters (kxM, M:measurements)
    /// Provided by the user
    //////////////////////////////////////////////////////
    CFP_host= new T[nmeas*CFP_Tsize];
    
    if (opts.CFP.set()){
      int nCFP_set = 0; //to check that all the values are set
      int cumulative_size=0;
      string filename(opts.CFP.value());
      ifstream file(filename.data());
      if (file.is_open()){
	string line;
	while(getline(file,line)){
       	  if (!line.empty()){ // if is empty, ignore the line
	    // Read file with values
	    Matrix values;
	    values=read_ascii_matrix(line);
	    int param_size = model.CFP_sizes[nCFP_set];
	    if(values.Nrows()!=param_size){
	      cerr << "CUDIMOT Error: Common Fixed Parameters number " << nCFP_set << " and file " << line << " do not match the dimensions specified for this model" << endl;
	      exit(-1);
	    }
	    if(values.Ncols()!=nmeas){
	      cerr << "CUDIMOT Error:  Common Fixed Parameters number " << nCFP_set << " and file " << line << " do not match the number of measurements: " << nmeas << endl;
	      exit(-1);
	    }
	    
	    for(int i=0;i<param_size;i++){
	      for (int j=0;j<nmeas;j++){
		CFP_host[cumulative_size+i+j*CFP_Tsize]=values(i+1,j+1);
	      }
	    }
	    nCFP_set++;
	    cumulative_size+=param_size;
	    
	  }
	} //end lines
	
	if(nCFP_set!=nCFP){
	  cerr << "CUDIMOT Error: The number of Common Fixed Parameters provided in file: " << filename.data() << " is not correct for this model. The number of Common Fixed Parameters must be " << nCFP << endl;
	  exit(-1);
	}
      }else{
	cerr << "CUDIMOT Error: Unable to open file with Common Fixed Parameters: " << filename.data() << endl; 
	exit(-1);
      }
    }else{
      // No CFP provided 
      if(CFP_Tsize){
	cerr << "CUDIMOT Error: The user must provide a file with a list of common fixed parameters for this model: Use option --CFP" << endl;
	exit(-1);
      }
    }
    //////////////////////////////////////////////////////
    


    //////////////////////////////////////////////////////
    /// Read Fixed Parameters (Nvoxels x M, M:measurements)
    /// Provided by the user
    //////////////////////////////////////////////////////
    FixP_host= new T[nvox*FixP_Tsize];
    
    if (opts.FixP.set()){
      int nFP_set = 0; //to check that all the values are set
      int cumulativeFP=0;

      for(int FP=0;FP<nFixP;FP++){
	// Read binary file with values for this part, logfile
	string name_file;
	name_file.append(opts.partsdir.value());
	name_file.append("/part_");
	name_file.append(num2str(opts.idPart.value()));
	name_file.append("/FixParam_");
	name_file.append(num2str(FP));

	ifstream in;
	long nbytes;
	int nvox_file,nmeas_file;
	in.open(name_file.data(), ios::in | ios::binary);
	in.read((char*)&nvox_file, 4);
	in.read((char*)&nmeas_file, 4);
	in.read((char*)&nbytes, sizeof(long));
	    
	if(nvox!=nvox_file || nmeas_file!=model.getNFixP_size(nFP_set)){
	  cerr << "CUDIMOT Error: The amount of data in the intermediate file " <<  name_file << " with Fixed Parameters is not correct" << endl;
	  exit(-1);
	}
	    
	Matrix FixPars;
	FixPars.ReSize(nmeas_file,nvox);
	in.read((char*)&FixPars(1,1),nbytes);
	in.close();
	    
	for (int v=0;v<nvox;v++){
	  for(int m=0;m<nmeas_file;m++){
	    FixP_host[v*FixP_Tsize+cumulativeFP+m]=FixPars(m+1,v+1);
	  }
	}
	    
	nFP_set++;
	cumulativeFP+=nmeas_file;
	    
      }	
    }else{
      // No FP provided 
      if(FixP_Tsize){
	cerr << "CUDIMOT Error: The user must provide a file with a list of Fixed Parameters for this model: Use option --FixP" << endl;
	exit(-1);
      }
    }
    //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////
    
    /// Allocate GPU memory
    cudaMalloc((void**)&params_gpu,nvoxFit_part*nparams*sizeof(T));
    cudaMalloc((void**)&CFP_gpu,CFP_Tsize*nmeas*sizeof(T));
    cudaMalloc((void**)&FixP_gpu,nvoxFit_part*FixP_Tsize*sizeof(T));
    if(opts.getPredictedSignal.value()){
      cudaMalloc((void**)&predSignal_gpu,nvoxFit_part*nmeas*sizeof(T));
    }
    if(opts.BIC_AIC.value()){
      cudaMalloc((void**)&BIC_gpu,nvoxFit_part*sizeof(T));
      cudaMalloc((void**)&AIC_gpu,nvoxFit_part*sizeof(T));
    }
    sync_check("Allocating Model Parameters on GPU\n");
    
    // Copy Fixed Common Parameters from host to GPU
    cudaMemcpy(CFP_gpu,CFP_host,CFP_Tsize*nmeas*sizeof(T),cudaMemcpyHostToDevice);
    sync_check("Copying Common-Fixed Model Parameters to GPU\n");
        
    /// If MCMC: allocate memory in host and GPU for samples;
    if(opts.runMCMC.value()){
      nsamples=(opts.njumps.value()/opts.sampleevery.value());   
    }else{
      nsamples=1;
    }
    samples_host = new T[nsamples*nparams*nvox];
    cudaMalloc((void**)&samples_gpu,nvoxFit_part*nparams*nsamples*sizeof(T));
    if(opts.rician.value()){
      tau_samples_host=new T[nsamples*nvox];
      cudaMalloc((void**)&tau_samples_gpu,nvoxFit_part*nsamples*sizeof(T));
    }
    sync_check("Allocating Samples on GPU\n");
  }
  
  template <typename T>
  Parameters<T>::~Parameters(){}
  
  template <typename T>
  T* Parameters<T>::getParametersPart(int part){
    if(part>=nparts){
      cerr << "CUDIMOT Error: Trying to get an incorrect part of the Parameters: " << part << ". There are only " << nparts << " parts and index starts at 0." << endl;
      exit(-1);
    }
    
    int initial_pos=part*size_part*nparams;
    int size=size_part;
    if(part==(nparts-1)) 
      size=size_last_part; // last part
    // Copy from host to GPU
    cudaMemcpy(params_gpu,&params_host[initial_pos],size*nparams*sizeof(T),cudaMemcpyHostToDevice);
    sync_check("Copying Model Parameters to GPU\n");

    return params_gpu;
  }
 
  template <typename T>
  T* Parameters<T>::getSamples(){
    return samples_gpu;
  }

  template <typename T>
  int Parameters<T>::getTsize_CFP(){
    return CFP_Tsize;
  }
 
  template <typename T>
  int Parameters<T>::getTsize_FixP(){
    return FixP_Tsize;
  }
  
  template <typename T>
  T* Parameters<T>::getCFP(){
    return CFP_gpu;
  }

  template <typename T>
  T* Parameters<T>::getFixP_part(int part){

    if(part>=nparts){
      cerr << "CUDIMOT Error: Trying to get an incorrect part of the Fixed Parameters: " << part << ". There are only " << nparts << " parts and index starts at 0." << endl;
      exit(-1);
    }
    int initial_pos=part*size_part*FixP_Tsize;
    int size=size_part;
    if(part==(nparts-1)) 
      size=size_last_part; // last part

    // Copy Fixed Parameters from host to GPU
    cudaMemcpy(FixP_gpu,&FixP_host[initial_pos],size*FixP_Tsize*sizeof(T),cudaMemcpyHostToDevice);
    sync_check("Copying Fixed Model Parameters to GPU\n");
    return FixP_gpu;
  }
  
  template <typename T>
  void Parameters<T>::copyParamsPartGPU2Host(int part){
    cudimotOptions& opts = cudimotOptions::getInstance();

    if(part>=nparts){
      cerr << "CUDIMOT Error: Trying to store an incorrect part of Parameters: " << part << ". There are only " << nparts << " parts and index starts at 0." << endl;
      exit(-1);
    }
    
    int size=size_part; // this ignores the extra voxels added
    if(part==(nparts-1)){
      size=size_last_part; // this ignores the extra voxels added
    }
    int initial_pos=part*size_part*nparams;
    
    cudaMemcpy(&params_host[initial_pos],params_gpu,size*nparams*sizeof(T),cudaMemcpyDeviceToHost);
    sync_check("Copying Model Parameters from GPU\n");
  }
  
  template <typename T>
  void Parameters<T>::copyParams2Samples(){
    samples_host=params_host;
  }

  template <typename T>
  void Parameters<T>:: calculate_predictedSignal_BIC_AIC(int mode, int part, T* meas){
    cudimotOptions& opts = cudimotOptions::getInstance();
    int initial_pos;
    int size=size_part; // this ignores the extra voxels added
    if(part==(nparts-1)){
      size=size_last_part; // this ignores the extra voxels added
    }

    if(opts.getPredictedSignal.value()){
      if(mode==0){
	// Calculate Predicted Signal of this part from parameters
	PredictedSignal.run(nvoxFit_part,nmeas,1,CFP_Tsize,FixP_Tsize,params_gpu,CFP_gpu,getFixP_part(part),predSignal_gpu);
      }else{
	// Calculate Predicted Signal of this part from samples (MCMC)
	PredictedSignal.run(nvoxFit_part,nmeas,nsamples,CFP_Tsize,FixP_Tsize,samples_gpu,CFP_gpu,getFixP_part(part),predSignal_gpu);

      }
      // Copy Predicted Signal to host
      initial_pos=part*size_part*nmeas;
      cudaMemcpy(&predSignal_host[initial_pos],predSignal_gpu,size*nmeas*sizeof(T),cudaMemcpyDeviceToHost);
      sync_check("Copying Predicted Signal from GPU\n");
    }

    if(opts.BIC_AIC.value()){
      if(mode==0){
	// Calculate BIC/AIC of this part from parameters
	bic_aic.run(nvoxFit_part,nmeas,1,CFP_Tsize,FixP_Tsize,meas,params_gpu,CFP_gpu,getFixP_part(part),BIC_gpu,AIC_gpu,tau_samples_gpu);
      }else{
	// Calculate BIC/AIC of this part from samples (MCMC)
	bic_aic.run(nvoxFit_part,nmeas,nsamples,CFP_Tsize,FixP_Tsize,meas,samples_gpu,CFP_gpu,getFixP_part(part),BIC_gpu,AIC_gpu,tau_samples_gpu);
	
      }
      // Copy Predicted Signal to host
      initial_pos=part*size_part;
      cudaMemcpy(&BIC_host[initial_pos],BIC_gpu,size*sizeof(T),cudaMemcpyDeviceToHost);
      cudaMemcpy(&AIC_host[initial_pos],AIC_gpu,size*sizeof(T),cudaMemcpyDeviceToHost);
      sync_check("Copying BIC/AIC from GPU\n");
    }
  }
 
  template <typename T>
  void Parameters<T>::copySamplesPartGPU2Host(int part){
    cudimotOptions& opts = cudimotOptions::getInstance();
    if(part>=nparts){
      cerr << "CUDIMOT Error: Trying to store an incorrect part of Samples: " << part << ". There are only " << nparts << " parts and index starts at 0." << endl;
      exit(-1);
    }
    
    int size=size_part; // this ignores the extra voxels added
    if(part==(nparts-1)){
      size=size_last_part; // this ignores the extra voxels added
    }
    int initial_pos=part*size_part*nparams*nsamples;
    
    cudaMemcpy(&samples_host[initial_pos],samples_gpu,size*nparams*nsamples*sizeof(T),cudaMemcpyDeviceToHost);
    sync_check("Copying Samples from GPU\n");

    if(opts.rician.value()){
      initial_pos=part*size_part*nsamples;
      cudaMemcpy(&tau_samples_host[initial_pos],tau_samples_gpu,size*nsamples*sizeof(T),cudaMemcpyDeviceToHost);
      sync_check("Copying Tau Samples from GPU\n");
    }
  }
  
  template <typename T>
  void Parameters<T>::writeSamples(){
    Log& logger = LogSingleton::getInstance();
    cudimotOptions& opts = cudimotOptions::getInstance();

    // Create a Matrix [nsamples X nvoxels] for each parameter
    vector<Matrix> samples;
    samples.resize(nparams);
    
    for(int par=0;par<nparams;par++){
      samples[par].ReSize(nsamples,nvox);
      samples[par]=0;
    }
    // Copy samples to each Matrix
    for(int vox=0;vox<nvox;vox++){  
      for(int par=0;par<nparams;par++){
	for(int sam=0;sam<nsamples;sam++){
	  samples[par](sam+1,vox+1)=samples_host[vox*nparams*nsamples+par*nsamples+sam];
	}
      }
    }
    
    // Write to file
    for(int par=0;par<nparams;par++){
      string file_name;
      file_name.append(opts.partsdir.value());
      file_name.append("/part_");
      file_name.append(num2str(opts.idPart.value()));
      file_name.append("/Param_"+num2str(par)+"_samples");
      ofstream out;
      out.open(file_name.data(), ios::out | ios::binary);
      out.write((char*)&nvox,4); // number of voxels
      out.write((char*)&nsamples,4); // number of measurements
      long size=nvox*nsamples*sizeof(Real); //need Real here (NEWMAT Object!)
      out.write((char*)&size,sizeof(long)); // number of bytes
      out.write((char*)&samples[par](1,1),size);
      out.close();
    }

    // Rician Noise: Write tau samples
    if(opts.rician.value()){
      Matrix tau_samples;
      tau_samples.ReSize(nsamples,nvox);
      tau_samples=0;
      for(int vox=0;vox<nvox;vox++){  
	for(int sam=0;sam<nsamples;sam++){
	  tau_samples(sam+1,vox+1)=tau_samples_host[vox*nsamples+sam];
	}
      }
      string file_name;
      file_name.append(opts.partsdir.value());
      file_name.append("/part_");
      file_name.append(num2str(opts.idPart.value()));
      file_name.append("/Tau_samples");
      ofstream out;
      out.open(file_name.data(), ios::out | ios::binary);
      out.write((char*)&nvox,4); // number of voxels
      out.write((char*)&nsamples,4); // number of measurements
      long size=nvox*nsamples*sizeof(Real); //need Real here (NEWMAT Object!)
      out.write((char*)&size,sizeof(long)); // number of bytes
      out.write((char*)&tau_samples(1,1),size);
      out.close();
    }

    if(opts.getPredictedSignal.value()){
      Matrix PredSignalM;
      PredSignalM.ReSize(nmeas,nvox);
      PredSignalM=0;
      for(int vox=0;vox<nvox;vox++){  
	for(int mea=0;mea<nmeas;mea++){
	  PredSignalM(mea+1,vox+1)=predSignal_host[vox*nmeas+mea];
	}
      }
      string file_name;
      file_name.append(opts.partsdir.value());
      file_name.append("/part_");
      file_name.append(num2str(opts.idPart.value()));
      file_name.append("/PredictedSignal");
      ofstream out;
      out.open(file_name.data(), ios::out | ios::binary);
      out.write((char*)&nvox,4); // number of voxels
      out.write((char*)&nmeas,4); // number of measurements
      long size=nvox*nmeas*sizeof(Real); //need Real here (NEWMAT Object!)
      out.write((char*)&size,sizeof(long)); // number of bytes
      out.write((char*)&PredSignalM(1,1),size);
      out.close();
    }

    if(opts.BIC_AIC.value()){
      Matrix BICM;
      BICM.ReSize(1,nvox);
      BICM=0;
      for(int vox=0;vox<nvox;vox++){  
	BICM(1,vox+1)=BIC_host[vox];
      }
      string file_name;
      file_name.append(opts.partsdir.value());
      file_name.append("/part_");
      file_name.append(num2str(opts.idPart.value()));
      file_name.append("/BIC");
      ofstream out;
      out.open(file_name.data(), ios::out | ios::binary);
      out.write((char*)&nvox,4); // number of voxels
      int nm=1;
      out.write((char*)&nm,4); // number of measurements
      long size=nvox*1*sizeof(Real); //need Real here (NEWMAT Object!)
      out.write((char*)&size,sizeof(long)); // number of bytes
      out.write((char*)&BICM(1,1),size);
      out.close();

      Matrix AICM;
      AICM.ReSize(1,nvox);
      AICM=0;
      for(int vox=0;vox<nvox;vox++){  
	AICM(1,vox+1)=AIC_host[vox];
      }
      file_name.clear();
      file_name.append(opts.partsdir.value());
      file_name.append("/part_");
      file_name.append(num2str(opts.idPart.value()));
      file_name.append("/AIC");
      out.open(file_name.data(), ios::out | ios::binary);
      out.write((char*)&nvox,4); // number of voxels
      out.write((char*)&nm,4); // number of measurements
      out.write((char*)&size,sizeof(long)); // number of bytes
      out.write((char*)&AICM(1,1),size);
      out.close();
    }
  }

  template <typename T>
  T* Parameters<T>::getTauSamples(){
    return tau_samples_gpu;
  }
  
  // Explicit Instantiations of the template
  template class Parameters<float>;
  template class Parameters<double>;
}
