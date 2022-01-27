#ifndef CUDIMOT_MODEL_H_INCLUDED
#define CUDIMOT_MODEL_H_INCLUDED
/**
 *
 * \class Model
 *
 * \brief A class for managing the Model specified by the Designer
 *
 * This class contains the information (parameters) of the diffusion MRI model specified by a designer.
 *
 * \author Moises Hernandez-Fernandez - FMRIB Image Analysis Group
 *
 * \date March 2017
 *
 *
 *  Copyright (C) 2005 University of Oxford 
 */

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

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "newmat.h"
#include "newimage/newimageall.h"
#include "checkcudacalls.h"
#include "cudimotoptions.h"
#include "cudimot.h"
#include "modelparameters.h"

namespace Cudimot{
  
  template <typename T>
  class Parameters;
  
  template <typename T>
  class Model{
    friend class Parameters<T>;
    
  private:
    /**
     * Number of parameters in the model (to estimate)
     */
    int nparams;

    /**
     * Vector with a default values for initializing the model parameters
     */
    vector<T> params_init;

    /**
     * If activated, the Designer has provided default values for initializing the model parameters
     */
    bool provided_params_init;
    
    ////////////////////////////////////////////////////////////////////////
    /// FIXED PARAMETERS different for each voxel, nifti volumes (ex.T1) ///
    ////////////////////////////////////////////////////////////////////////
 
    /**
     * Number of Fixed Parameters in the model (different for each voxel). NIfTI files must be used for specifying the fixed parameters
     */
    int nFixP; 
    
    /**
     * Total size of volumes of Fixed Parameter
     */
    int FixP_Tsize; 
    
    /**
     * Vector with the size in volumes of each Fixed Parameter
     */
    vector<int> FixP_sizes;
    /////////////////////////////////////////////////////////////////

    
    /////////////////////////////////////////////////////////////////
    /// COMMON (to all voxels) FIXED PARAMETERS (bvals,bvecs,...) ///
    /////////////////////////////////////////////////////////////////
    
    /**
     * Number of Common Fixed Parameters in the model (common to all voxels)
     */
    int nCFP;
    
    /**
     * Total size of Common Fixed Parameters (without counting measurements)
     */
    int CFP_Tsize; 

    /**
     * Vector with the size of each Common Fixed Parameter (without counting measurements)
     */
    vector<int> CFP_sizes;
    /////////////////////////////////////////////////////////////////
    
    ///////////////////////
    /// Priors for MCMC ///
    ///////////////////////

    /**
     * Type of the bounds used in each parameter during MCMC
     * bound{NONE,BMIN,BMAX,BMINMAX};
     * 0: No bounds
     * 1: Bounded(min,)
     * 2: Bounded(,max)
     * 3: Bounded(min,max)
     */
    vector<int> bound_types;

    /**
     * Type of the prior used in each parameter during MCMC
     * prior{NOPRIOR,GAUSSPRIOR,GAMMAPRIOR,ARDPRIOR,SINPRIOR};
     * 0: No prior
     * 1: Gaussian(mean,sd)
     * 2: Gamma(alpha,beta)
     * 3: ARD(fudge_factor)
     * 4: sin()
     * 5: custom()
     */
    vector<int> prior_types;
    
    /**
     * Minimum bound
     */
    vector<T> bounds_min;

    /**
     * Maximum bound
     */
    vector<T> bounds_max;
    

    /**
     * First parameter of the priors (mean, alpha or fudge_factor)
     */
    vector<T> priors_a;
    
    /**
     * Second parameter of the priors (standard_deviation or beta)
     */
    vector<T> priors_b;

    /**
     * Vector to specify if a parameter is fixed. 0 non-fixed, 1 fixed
     */
    vector<int> fixed;
    
     /**
     * Number of parameter-value combination to try in GridSearch
     */
    int gridCombs;
    
    /**
     * Ids of parameters in the grid
     */
    vector<int> gridParams;

    /**
     * Number of parameters in each combination of the grid
     */
    int nGridParams;

    /**
     * Values to try in GridSearch
     */
    T* grid;
    
    /**
     * This method reads the model configuration file and sets the class attributes. File name must be provided at execution time
     */
    void Modelparser(string default_priors_file);

    /**
     * This method reads the values of each parameter to try in GridSearch (if used). File name must be provided at execution time
     */
    void Parser_gridSearch();

     /**
     * Recursive Method for setting the grid id gridSearch used
     */
    void set_grid(int level,int nGridParams, int &gridComb, T* comb, vector <vector <T> > &gridTmp, vector<int> gridParams,T* grid);
        
  public:
    /**
     * Constructor
     */
    Model(string default_priors_file);

    /**
     * Destructor
     */
    ~Model();

    /** 
     * @return The number of parameters in the model (to estimate)
     */
    int getNparams();

    /** 
     * @return The number of Fixed parameters in the model (NIfTI files must be used for specifying the fixed parameters)
     */
    int getNFixP();   

    /** 
     * @param Fixed parameter Identifier
     * @return The number of Fixed parameters in the model (NIfTI files must be used for specifying the fixed parameters)
     */
    int getNFixP_size(int id_FP);   

    /**
     * @return Indicates if the Designer has provided default values for initializing the model parameters
     */
    bool initProvided();

    /**
     * @param Parameter identifier 
     * @return The default value for initializing one of the model parameters
     */
    T getParam_init(int id_param);

    /**
     * @return Vector with the the type of each parameter bounds
     */
    vector<int> getBound_types();

    /**
     * @return Vector with the minimum bound of each parameter
     */
    vector<T> getBounds_min();
    
    /**
     * @return Vector with the maximum bound of each parameter
     */
    vector<T> getBounds_max();

    /**
     * @return Vector with the type of each parameter prior
     */
    vector<int> getPrior_types();

    /**
     * @return Vector with the first argument of each prior
     */
    vector<T> getPriors_a();
    
    /**
     * @return Vector with the second argument of each prior
     */
    vector<T> getPriors_b();

    /**
     * @return Vector for specifying if a parameter is fixed or not
     */
    vector<int> getFixed();

    /*
     * @return Number of parameters in each combination of the grid
     */
    int getNGridParams();

    /*
     * @return Array of paramer
     */
    vector<int> getGridParams();
    
    /*
     * @return Number of combinations in the grid
     */
    int getGridCombs();

    /*
     * @return The grid with all the combinations values for Grid Search
     */
    T* getGrid();

  };
}

#endif
