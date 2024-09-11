#ifndef CUDIMOT_LEVENBERG_MARQUARDT_H_INCLUDED
#define CUDIMOT_LEVENBERG_MARQUARDT_H_INCLUDED

/**
 *
 * \class Levenberg_Marquardt
 *
 * \brief A class for fitting a model to some dataset on a GPU using Levenberg-Marquardt algorithm
 *
 * \author Moises Hernandez-Fernandez - FMRIB Image Analysis Group
 *
 * \date March 2017
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
  class Levenberg_Marquardt{

  private:

    /**
     * Maximum number of iterations
     */
    int max_iterations;
    
    /**
     * If activeted, use Marquardt contribution
     */
    bool Marquardt;

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
     * 0: Parameter non fixed, 1: Parameter fixed
     */
    int* fixed_host;

    /**
     * Activate debugging messages for a voxel. It prints the value of some variables at certain steps of Levenberg_Marquardt (Parameters, PredictedSignal, Derivatives, Gradient, Proposed Parameters)
     */
    bool DEBUG;

    /**
     * Number of the voxel (starts at 0) to debug if debugging is activated
     */
    int debugVOX;
 
  public:
    
    /**
     * Constructor
     * @param bound_types Vector with the type of each bound type
     * @param bounds_min Vector with the lower bound of each parameter
     * @param bounds_max Vector with the upper bound of each parameter
     * @param fixed Vector with information to know if parameters are fixed
     */
    Levenberg_Marquardt(vector<int> bound_types, vector<T> bounds_min, vector<T> bounds_max,vector<int> fixed);

    /**
     * Run Levenberg-Marquard on the GPU
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
