/* Model.cu
   
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

// This class contains the information of the Model

#include "Model.h"

namespace Cudimot{
  
  // This method reads the list of parameters from a file. File must be provided at execution time.
  template <typename T>
  void Model<T>::Modelparser(string default_priors_file){
    cudimotOptions& opts = cudimotOptions::getInstance();

    string input_file;
    if(opts.priorsfile.value()!=""){
      input_file.assign(opts.priorsfile.value());
    }else{
      input_file.assign(default_priors_file);
    }
    ifstream file_parameters(input_file.data());
    //string mark_start="NEW_MODEL";
    string mark_comment="//";
    // marks for defining initialization values
    string Pinit_mark="P_init[";
    string bound_mark="bounds[";
    string prior_mark="prior[";
    string CFP_size_mark="CFP_size[";
    string FixP_size_mark="FixP_size[";
    string mark_end="]";
    string delimiter=",";
    string line;
    
    nparams=NPARAMS;
    provided_params_init=false;
    nCFP=0;
    CFP_Tsize=0;
    nFixP=0;
    FixP_Tsize=0;

    bound_types.resize(nparams);
    bounds_min.resize(nparams);
    bounds_max.resize(nparams);
    prior_types.resize(nparams);
    priors_a.resize(nparams);
    priors_b.resize(nparams);
    for(int p=0;p<nparams;p++){
      bound_types[p]=NOBOUNDS;
      bounds_min[p]=0;
      bounds_max[p]=0;
      prior_types[p]=NOPRIOR;
      priors_a[p]=0;
      priors_b[p]=0;
    }
    
    // File reading steps 0:Ignore, 1:Reading
    int reader_state=1;
    if (file_parameters.is_open()){
      while(getline(file_parameters,line)){
	
       	//if(!line.compare(mark_start)) reader_state=1; // Starting MARK
				if(reader_state){ 
	  			// Reading useful information
	  			int pos0, pos1;
	  
	 				if((pos0=line.find(mark_comment))>=0){
	    			// If line commented out using "//", ignore the line
	    			continue;
	  			}
	  
					//if((pos0=line.find(mark_NPARAMS))>=0){ 
					// Read Number of Parameters
					//pos0=pos0+mark_NPARAMS.length();
					//line.erase(0,pos0);
					//nparams=atoi(line.data());
					//}
	  
					if((pos0=line.find(Pinit_mark))>=0){
						// Read Parameter Initialization
						provided_params_init=true;
						pos0=pos0+Pinit_mark.length();
						pos1=line.rfind(mark_end);
						int len = pos1-pos0;
						string vals=line.substr(pos0,len);
						string token;
						while ((pos0 = vals.find(delimiter))>=0){
							params_init.push_back(strtod(vals.substr(0,pos0).data(),NULL));
							vals.erase(0,pos0+delimiter.length());
						}
						params_init.push_back(strtod(vals.data(),NULL));
						if(params_init.size()!=NPARAMS){
							cerr << "CUDIMOT Error: The number of values for initializating the parameters " << params_init.size() <<" do not match the number of parameters of this model " << NPARAMS << " in the file: "<< input_file << endl; 
							exit(-1);
						}
					}
	
					if((pos0=line.find(bound_mark))>=0){ 
						// Read bounds for a parameter
						pos0=pos0+bound_mark.length();
						pos1=line.rfind(mark_end);
						int len = pos1-pos0;
						int idBound=atoi(line.substr(pos0,len).data());

						if(idBound>=NPARAMS || idBound<0){
							cerr << "CUDIMOT Error: Wrong number of parameter [" << idBound << "] when specifying the bounds in file" << input_file << endl;
							exit(-1);
						}
						
						bound_types[idBound]=BMINMAX; // assume(min,max)
						// look for the parameters
						pos0=line.rfind("=")+1;
						pos1=line.rfind("(");
						len = pos1-pos0;
						string type =  line.substr(pos0,len);
						
						// First parameter
						pos0=pos1+1;
						pos1=line.rfind(",");
						len= pos1-pos0;
						if(len==0){ 
							// No first parameter
							bound_types[idBound]=BMAX;
						}else{
							bounds_min[idBound] = strtod(line.substr(pos0,len).data(),NULL);
						}

						// Second parameter
						pos0=pos1+1;
						pos1=line.rfind(")");
						len= pos1-pos0;
						if(len==0 && bound_types[idBound]==BMINMAX){ 
							// First parameter and No second parameter
							bound_types[idBound]=BMIN;
						}else if(len==0){
							cerr << "CUDIMOT Error: For bounded priors, at least a min or a max value must be specified: bounds[x]=(min,), bounds[x]=(,max) or bounds[x]=(min,max)" << type << endl; 
							exit(-1);
						}else{
							bounds_max[idBound] = strtod(line.substr(pos0,len).data(),NULL);
						}
					}
	  
					if((pos0=line.find(prior_mark))>=0){
						// Read a Prior
						pos0=pos0+prior_mark.length();
						pos1=line.rfind(mark_end);
						int len = pos1-pos0;
						int idPrior=atoi(line.substr(pos0,len).data());
						if(idPrior>=NPARAMS || idPrior<0){
							cerr << "CUDIMOT Error: Wrong number of parameter [" << idPrior << "] when specifying the priors in file" << input_file << endl;
							exit(-1);
						}
						
						pos0=line.rfind("=")+1;
						pos1=line.rfind("(");
						len = pos1-pos0;
						string type =  line.substr(pos0,len);
						
						if(type.compare("Gaussian")==0){
							prior_types[idPrior]=GAUSSPRIOR;
						}else if(type.compare("Gamma")==0){
							prior_types[idPrior]=GAMMAPRIOR;
						}else if(type.compare("ARD")==0){
							prior_types[idPrior]=ARDPRIOR;
						}else if(type.compare("sin")==0){
							prior_types[idPrior]=SINPRIOR;
						}else if(type.compare("custom")==0){
							prior_types[idPrior]=CUSTOM;
						}else{
							cerr << "CUDIMOT Error: Prior " << type << " is a not a recognised prior type" << endl; 
							exit(-1);
						}
						if(prior_types[idPrior]<=3){ 
							// look for the parameters
							// First parameter
							pos0=pos1+1;
							pos1=line.rfind(",");
							len= pos1-pos0;
							if(len==0){
								cerr << "CUDIMOT Error: A prior of type " << type << " needs a first argument" << endl; 
								exit(-1);
							}else{
								priors_a[idPrior] = strtod(line.substr(pos0,len).data(),NULL);
							}
							
							if(prior_types[idPrior]<=2){ 
								// Second parameter
								pos0=pos1+1;
								pos1=line.rfind(")");
								len= pos1-pos0;
								if(len==0){
									cerr << "CUDIMOT Error: A prior of type " << type << " needs a second argument" << endl; 
									exit(-1);
								}else{
									priors_b[idPrior] = strtod(line.substr(pos0,len).data(),NULL);
								}
							}
						}
					}
	  
					/*if((pos0=line.find(CFP_size_mark))>=0){ 
					// Read common fixed param Sizes
					pos0=pos0+CFP_size_mark.length();
					pos1=line.rfind(mark_end);
					int len = pos1-pos0;
					string vals=line.substr(pos0,len);
					string token;
					while ((pos0 = vals.find(delimiter))>=0){
						int size=atoi(vals.substr(0,pos0).data());
						CFP_sizes.push_back(size);
						CFP_Tsize+=size;
						nCFP++;
						vals.erase(0,pos0+delimiter.length());
					}
					int size=atoi(vals.data());
					CFP_sizes.push_back(size);
					CFP_Tsize+=size;
					nCFP++;
				}

				if((pos0=line.find(FixP_size_mark))>=0){ 
					// Read fixed param Sizes
					pos0=pos0+FixP_size_mark.length();
					pos1=line.rfind(mark_end);
					int len = pos1-pos0;
					string vals=line.substr(pos0,len);
					string token;
					while ((pos0 = vals.find(delimiter))>=0){
						int size=atoi(vals.substr(0,pos0).data());
						FixP_sizes.push_back(size);
						FixP_Tsize+=size;
						nFixP++;
						vals.erase(0,pos0+delimiter.length());
					}
					int size=atoi(vals.data());
					FixP_sizes.push_back(size);
					FixP_Tsize+=size;
					nFixP++;
					}*/
				}
			}
			file_parameters.close();

      nCFP=NCFP;
      for(int i=0;i<nCFP;i++){
				int size=MODEL::CFP_size[i];
				CFP_sizes.push_back(size);
				CFP_Tsize+=size;
      }

      nFixP=NFIXP;
      for(int i=0;i<nFixP;i++){
				int size=MODEL::FixP_size[i];
				FixP_sizes.push_back(size);
				FixP_Tsize+=size;
      }
            
      // CHECK
      if(!nparams){
				string input_file(opts.priorsfile.value());
				cerr << "CUDIMOT Error: Number of parameters must be greater than 0. Please check your parameter specification: " << input_file.data() << endl;
				exit(-1);
      }

      //if(nFixP =! NFIXED_PARAMS){
      //cerr << "CUDIMOT Error: Number of Fixed Parametersmust be greater than 0. Please check your parameter specification: " << input_file.data() << endl;
      //	exit(-1);
      //}

    }else{
      cerr << "CUDIMOT Error: Unable to open file with Priors information: " << input_file.data() << endl; 
      exit(-1);
    }

    fixed.resize(nparams);
    for(int p=0;p<nparams;p++) 
      fixed[p]=0;
    if(opts.fixed.value()!=""){
      stringstream ss(opts.fixed.value());
      vector<int> vect;
      int i;
      while (ss >> i){
				if(i>=nparams){
	   			cerr << "CUDIMOT Error: Wrong number of parameter " <<  i << " when specifying the fixed parameters" << endl;
	   			exit(-1);
				}
        vect.push_back(i);
        if (ss.peek() == ',')
	  		ss.ignore();
      }
      for (i=0; i<vect.size(); i++){
				fixed[vect[i]]=1;
      }
    }
  }

  // If gridSearch, this method reads from a file the values to search and set the grid with all the combinations. File can be provided at execution time
  template <typename T>
  void Model<T>::Parser_gridSearch(){
    cudimotOptions& opts = cudimotOptions::getInstance();
    string input_file(opts.gridSearch.value());
    ifstream file_parameters(input_file.data());
    string mark_comment="//";
    string search_mark="search[";
    string parid_end="]";
    
    string delimiter=",";
    string line;
    
    vector< vector <T> > gridTmp;
    gridTmp.resize(NPARAMS);

    if (file_parameters.is_open()){
      while(getline(file_parameters,line)){
				int pos0, pos1;
				if((pos0=line.find(mark_comment))>=0){
					// If line commented out using "//", ignore the line
					continue;
				}
			  
				if((pos0=line.find(search_mark))>=0){
					// Read Search values for one parameter
					pos0=pos0+search_mark.length();
					pos1=line.rfind(parid_end);
					int len = pos1-pos0;
					int idSearch=atoi(line.substr(pos0,len).data());
					if(idSearch>=NPARAMS || idSearch<0){
						cerr << "CUDIMOT Error: Wrong number of parameter [" << idSearch << "] when specifying the GridSearch values in file" << input_file << endl;
						exit(-1);
					}
					pos1=line.rfind("("); // ignore "]=("
					pos0=pos1+1; 
					pos1=line.rfind(")");
					len = pos1-pos0;
					string vals=line.substr(pos0,len);

					//string token;
					while ((pos0 = vals.find(delimiter))>=0){
						gridTmp[idSearch].push_back(strtod(vals.substr(0,pos0).data(),NULL));
						vals.erase(0,pos0+delimiter.length());
					}
					gridTmp[idSearch].push_back(strtod(vals.data(),NULL));
				}
      }
    }else{
      cerr << "CUDIMOT Error: Unable to open GridSearch file configuration: " << input_file.data() << endl; 
      exit(-1);
    }
    
    nGridParams=0;
    gridCombs=1;
    for(int i=0;i<NPARAMS;i++){
      if(gridTmp[i].size()){
				gridParams.push_back(i);
				nGridParams++;
				gridCombs*=gridTmp[i].size();
      }
    }
    if(nGridParams==0) gridCombs=0;
    grid = new T[gridCombs*nGridParams];
   
    int ncomb=0;
    T* comb = new T[nGridParams];
    set_grid(0,nGridParams,ncomb,comb,gridTmp,gridParams,grid);
  }

  template <typename T>
  void Model<T>::set_grid(int level,
			  int nGridParams,
			  int &ncomb,
			  T* comb,
			  vector <vector <T> > &gridTmp,
			  vector<int> gridParams,
			  T* grid)
  {
    for(int i=0;i<gridTmp[gridParams[level]].size();i++){
      if(level==(nGridParams-1)){
				comb[level] = gridTmp[gridParams[level]][i];
				for(int j=0;j<nGridParams;j++){
					grid[ncomb*nGridParams+j]=comb[j];
				}
				ncomb++;
      }else{
				comb[level] = gridTmp[gridParams[level]][i];
				set_grid(level+1,nGridParams,ncomb,comb,gridTmp,gridParams,grid);
      }
    }
      
  }
  
  template <typename T>
  Model<T>::Model(string default_priors_file){
    cudimotOptions& opts = cudimotOptions::getInstance();
    /// Read text file with parameters information (initialization, priors)
    Modelparser(default_priors_file);
    if(opts.gridSearch.value()!=""){
      Parser_gridSearch(); // if grid_search, set the grid
    }
  }
  
  template <typename T>
  Model<T>::~Model(){}

  template <typename T>
  int Model<T>:: getNparams(){
    return nparams;
  }

  template <typename T>
  int Model<T>:: getNFixP(){
    return nFixP;
  }

  template <typename T>
  int Model<T>:: getNFixP_size(int id_FP){
    return FixP_sizes[id_FP];
  }
  
  template <typename T>
  bool Model<T>:: initProvided(){
    return provided_params_init;
  }
  
  template <typename T>
  T Model<T>:: getParam_init(int id_param){
    return params_init[id_param];
  }
  
  template <typename T>
  vector<int> Model<T>::getBound_types(){
    return bound_types;
  }

  template <typename T>
  vector<T> Model<T>::getBounds_min(){
    return bounds_min;
  }
  
  template <typename T>
  vector<T> Model<T>::getBounds_max(){
    return bounds_max;
  }

  template <typename T>
  vector<int> Model<T>::getPrior_types(){
    return prior_types;
  }

  template <typename T>
  vector<T> Model<T>::getPriors_a(){
    return priors_a;
  }
  
  template <typename T>
  vector<T> Model<T>::getPriors_b(){
    return priors_b;
  }

  template <typename T>
  vector<int> Model<T>::getFixed(){
    return fixed;
  }

  template <typename T>
  int Model<T>::getNGridParams(){
    return nGridParams;
  }

  template <typename T>
  int Model<T>::getGridCombs(){
    return gridCombs;
  }

  template <typename T>
  vector<int> Model<T>::getGridParams(){
    return gridParams;
  }

  template <typename T>
  T* Model<T>::getGrid(){
    return grid;
  }

  // Explicit Instantiations of the template
  template class Model<float>;
  template class Model<double>;
}
