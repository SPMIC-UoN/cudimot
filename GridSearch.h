#ifndef CUDIMOT_GRID_SEARCH_H_INCLUDED
#define CUDIMOT_GRID_SEARCH_H_INCLUDED

/**
 *
 * \class GridSearch
 *
 * \brief A class for initialising the parameters, using Grid Search on a GPU, before the fitting process of a model.
 *
 * \author Moises Hernandez-Fernandez - FMRIB Image Analysis Group
 *
 * \date April 2017
 *
 *
 *  Copyright (C) 2005 University of Oxford 
 */

#include <vector>
#include "gridOptions.h"
#include "cudimot.h"
#include "checkcudacalls.h"
#include "cudimotoptions.h"

using std::vector;

namespace Cudimot{
  template <typename T>
  class GridSearch{

  private:

    /**
     * All the combinations of parameter values to try on the GPU
     */
    T* grid_gpu;

    /**
     * Number of combinations in the grid (size is this * NPARAMS)
     */
    int gridCombs;
    
    /*
     * The ID of the parameters in the grid
     */
    int* gridParams_host;
    
    /**
     * Number of parameters in each combination of the grid
     */
    int nGridParams;

    /**
     * Type of each parameter bounds
     */
    int* bound_types_host;
    
    /**
     * Minimum bound of each parameter
     */
    float* bounds_min_host;
    
    /**
     * Maximum bound of each parameter
     */
    float* bounds_max_host;
    
    /**
     * Activate debugging messages for a voxel. It prints the evaluation of each set of parameters during GridSearch
     */
    bool DEBUG;
    
    /**
     * Number of the voxel (starts at 0) to debug if debugging is activated
     */
    int debugVOX;
    
  public:
    
    /**
     * Constructor
     * @param number of parameters in each combination of the grid
     * @param Ids of parameters in the grid
     * @param Number of parameter-value combination to try in GridSearch
     * @param the grid with all the combinations values for Grid Search
     * @param bound_types Vector with the type of each bound type
     * @param bounds_min Vector with the lower bound of each parameter
     * @param bounds_max Vector with the upper bound of each parameter
     */
    GridSearch(int nGridParams, vector<int> gridParams, int gridCombs, T* grid_host, vector<int> bou_types, vector<T> bou_min, vector<T> bou_max);

    /**
     * Run GridSearch on the GPU
     * @param nvox Number of voxels in the dataset
     * @param nmeas Number of measurements in the dataset
     * @param CFP_size Size (x M measurements) of the common fixed parameters of the model
     * @param FixP_size Size (x N voxels) of the fixed parameters of the model
     * @param meas Measurements of all the voxels of the dataset (on GPU)
     * @param params Value of the parameters to estimate of the model for all the voxels of the dataset (on GPU)
     * @param CFP Common (to all the voxels) fixed parameters of the model (on GPU). CFP_size*nmeas
     * @param FixP Fixed parameters of the model (on GPU). FixP_size*nvoxels
     */
    void run( int nvox, int nmeas,
	      int CFP_size, int FixP_size,
	      T* meas, T* params,
	      T* CFP, T* FixP);
  };
}

#endif
