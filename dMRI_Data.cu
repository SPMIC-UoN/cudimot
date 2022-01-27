/* dMRI_Data.cu

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

#include "dMRI_Data.h"

namespace Cudimot{
  
  template <typename T>
  void dMRI_Data<T>::remove_NonPositive_entries(NEWMAT::ColumnVector& Voxdata){ 
    //Zero, Negative Entries can be obtained from spline interpolation 
    int pos; 
    float MinS=Voxdata.Minimum1(pos); 
    float MaxS=Voxdata.Maximum();
    if (MinS<=0 && MaxS>0){  
      //when there are some non-positive entries, but not all are zero
      vector<int> minpositions;
      while (MinS<=0){
	minpositions.push_back(pos);
	Voxdata(pos)=MaxS;    //temporarilly make the non-positive values Max
	MinS=Voxdata.Minimum1(pos);
      }
      MinS=Voxdata.Minimum(); //Now find the Minimum of positive entries
      for (unsigned int i=0; i<minpositions.size(); i++)
	Voxdata(minpositions[i])=MinS; //Replace non-positive entries with that minimum
    }
  }
  
  template <typename T>
  dMRI_Data<T>::dMRI_Data(){
    
    cudimotOptions& opts = cudimotOptions::getInstance();
    
    // Read the binary file with data (genereted previously in split_parts)
    ifstream in;
    long nbytes;
    string file_input;
    file_input.append(opts.partsdir.value());
    file_input.append("/part_");
    file_input.append(num2str(opts.idPart.value()));
    file_input.append("/data");
    
    in.open(file_input.data(), ios::in | ios::binary);
    in.read((char*)&nvox, 4);
    in.read((char*)&nmeas, 4);
    in.read((char*)&nbytes, sizeof(long));
    
    if(nvox<=0 || nmeas<=0){
      cerr << "CUDIMOT Error: The number of voxels and diffusion-weighted measurements in the input file must be greater than 0" << endl;
      exit (EXIT_FAILURE);
    }
    
    cout << "Number of Voxels to compute: " << nvox << endl;  
    cout << "Number of Measurements: " << nmeas << endl;  

    // Read diffusion-weighted measurements
    dataM.ReSize(nmeas,nvox);
    in.read((char*)&dataM(1,1),nbytes);
    in.close();
    
    // Data is divided into parts
    nparts=nvox/SIZE_PART;
    size_part=SIZE_PART;
    if(nvox%SIZE_PART) nparts++;
    size_last_part = nvox - ((nparts-1)*SIZE_PART);
    if(size_last_part<(SIZE_PART*0.5)){ 
      // if last part is too small, we distribute its voxels between the others parts
      if(nparts-1){ // More than 1 part
	size_part = size_part + size_last_part/(nparts-1);
	nparts--;
      }else{
	size_part = 0;
      }
      size_last_part = nvox - ((nparts-1)*size_part);
    }
    
    // Allocate memory on host and GPU for measurements
    int max_nvox =  max(size_part,size_last_part);
    // number of voxels can be a non-multiple of voxels per block, so somethreads could access to non-allocated memory. We use the closest upper multiple. The added voxels will be ignored.
    nvoxFit_part=int(max_nvox/MAX_VOXELS_BLOCK)*MAX_VOXELS_BLOCK;
    if(max_nvox%MAX_VOXELS_BLOCK) nvoxFit_part=nvoxFit_part+MAX_VOXELS_BLOCK;
    meas_host=new T[nvoxFit_part*nmeas];
    cudaMalloc((void**)&meas_gpu,nvoxFit_part*nmeas*sizeof(T));
    sync_check("Allocating dMRI_Data on the GPU");
  }
  
  template <typename T>
  dMRI_Data<T>::~dMRI_Data(){
    //cudaFree(meas_gpu);
    //sync_check("Deallocating dMRI_Data from GPU");
  }
  
  template <typename T>
  int dMRI_Data<T>::getNvoxFit_part() const{
    return nvoxFit_part;
  }
  
  template <typename T>
  int dMRI_Data<T>::getNmeas() const{
    return nmeas;
  }
  
  template <typename T>
  int dMRI_Data<T>::getNparts() const{
    return nparts;
  }
  
  // Returns size of part in the second parameter
  template <typename T>
  T* dMRI_Data<T>::getMeasPart(int part, int &sp){
    
    cudimotOptions& opts = cudimotOptions::getInstance();
    
    if(part>=nparts){
      cerr << "CUDIMOT Error: Trying to get an incorrect part of the data: " << part << ". There are only " << nparts << " parts and index starts at 0." << endl;
      exit(-1);
    }
    
    int size=size_part;
    int initial_vox=part*size_part;
    if(part==(nparts-1)){
      size=size_last_part;
    }
    
    cout << endl << endl << endl << "Part " << part+1 << " of " << nparts << ": processing " << size << " voxels" << endl;

    int vox=0;
    for(vox=0;vox<size;vox++){
      ColumnVector voxmeas;
      voxmeas=dataM.Column(initial_vox+vox+1);
      if(opts.rician.value()) remove_NonPositive_entries(voxmeas); //So that log(data) does not give infinity in the likelihood
      for(int m=0;m<nmeas;m++){
	meas_host[vox*nmeas+m]=voxmeas(m+1);
      }
    }
    // Fill with 0 the rest of the vector
    for(;vox<nvoxFit_part;vox++){
      for(int m=0;m<nmeas;m++){
	meas_host[vox*nmeas+m]=0;
      }
    }
    
    // Copy from host to GPU
    cudaMemcpy(meas_gpu,meas_host,nvoxFit_part*nmeas*sizeof(T),cudaMemcpyHostToDevice);
    sync_check("Copying dMRI_Data to GPU");
    sp=nvoxFit_part; 
    return meas_gpu;
  }
  
  template class dMRI_Data<float>;
  template class dMRI_Data<double>;
}
