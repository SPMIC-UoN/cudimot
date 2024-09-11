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
    std::vector<T> params_init;

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
    std::vector<int> FixP_sizes;
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
    std::vector<int> CFP_sizes;
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
    std::vector<int> bound_types;

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
    std::vector<int> prior_types;
    
    /**
     * Minimum bound
     */
    std::vector<T> bounds_min;

    /**
     * Maximum bound
     */
    std::vector<T> bounds_max;
    

    /**
     * First parameter of the priors (mean, alpha or fudge_factor)
     */
    std::vector<T> priors_a;
    
    /**
     * Second parameter of the priors (standard_deviation or beta)
     */
    std::vector<T> priors_b;

    /**
     * Vector to specify if a parameter is fixed. 0 non-fixed, 1 fixed
     */
    std::vector<int> fixed;
    
     /**
     * Number of parameter-value combination to try in GridSearch
     */
    int gridCombs;
    
    /**
     * Ids of parameters in the grid
     */
    std::vector<int> gridParams;

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
    void Modelparser(std::string default_priors_file);

    /**
     * This method reads the values of each parameter to try in GridSearch (if used). File name must be provided at execution time
     */
    void Parser_gridSearch();

     /**
     * Recursive Method for setting the grid id gridSearch used
     */
    void set_grid(int level,int nGridParams, int &gridComb, T* comb, std::vector <std::vector <T> > &gridTmp, std::vector<int> gridParams,T* grid);
        
  public:
    /**
     * Constructor
     */
    Model(std::string default_priors_file);

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
    std::vector<int> getBound_types();

    /**
     * @return Vector with the minimum bound of each parameter
     */
    std::vector<T> getBounds_min();
    
    /**
     * @return Vector with the maximum bound of each parameter
     */
    std::vector<T> getBounds_max();

    /**
     * @return Vector with the type of each parameter prior
     */
    std::vector<int> getPrior_types();

    /**
     * @return Vector with the first argument of each prior
     */
    std::vector<T> getPriors_a();
    
    /**
     * @return Vector with the second argument of each prior
     */
    std::vector<T> getPriors_b();

    /**
     * @return Vector for specifying if a parameter is fixed or not
     */
    std::vector<int> getFixed();

    /*
     * @return Number of parameters in each combination of the grid
     */
    int getNGridParams();

    /*
     * @return Array of paramer
     */
    std::vector<int> getGridParams();
    
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
