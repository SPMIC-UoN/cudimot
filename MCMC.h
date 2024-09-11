#ifndef CUDIMOT_MCMC_H_INCLUDED
#define CUDIMOT_MCMC_H_INCLUDED

/**
 *
 * \class MCMC
 *
 * \brief A class for fitting a model to some dataset on a GPU using MCMC algorithm
 *
 * \author Moises Hernandez-Fernandez - FMRIB Image Analysis Group
 *
 * \date March 2017
 *
 *
 *  Copyright (C) 2005 University of Oxford 
 */

#include <vector>
#include <curand_kernel.h>
#include <curand.h>
#include <curand_kernel.h>
#include "gridOptions.h"
#include "cudimot.h"
#include "checkcudacalls.h"
#include "cudimotoptions.h"

using std::vector;

namespace Cudimot{

  template <typename T>
  class MCMC{

  private:

    /** 
     * The number of voxels to fit in a part
     */
    int nvoxFit_part;

    /** 
     * Total number of jumps in MCMC to be discarded at the beginning (default is 5000)
     */
    int nburnin; 
    
    /** 
     * Number of jumps in MCMC after burnin period (default is 1250)
     */
    int njumps;
    
    /** 
     * Number of samples to record per parameter and voxel
     */
    int nsamples;

    /** 
     * Number of jumps between sample recording
     */
    int sampleevery;

    /** 
     * Number of jumps between updates of the proposal standard deviations
     */
    int updateproposalevery;

    /** 
     * if false, do not update the proposal standard deviations durinf the recording step of MCMC.
     */
    bool updateproposal;

    /** 
     * If activated, use Rician noise modelling
     */
    bool RicianNoise;

    /**
     * State of several random number generators on the GPU
     */
    curandState* randStates;

    /**
     * Standard Deviation of Proposal Distributions on the GPU
     */
     T* propSD;

    /**
     * Standard Deviation of Proposal Distributions for Tau parameter (rician noise) on the GPU
     */
    T* tau_propSD;
    
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
     * Type of each parameter prior
     */
    int* prior_types_host;
    
    /**
     * First argument of each parameter prior
     */
    float* priors_a_host;
    
    /**
     * Second argument of each parameter prior
     */
    float* priors_b_host;

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
     * @param nvoxFitpart Number of voxel of the data to fit
     * @param bound_types Vector with the type of each bound type
     * @param bounds_min Vector with the lower bound of each parameter
     * @param bounds_max Vector with the upper bound of each parameter
     * @param prior_types Vector with the type of each parameter type
     * @param prior_a Vector with the first argument of each prior
     * @param prior_b Vector with the second argument of each prior
     * @param fixed Vector with information to know if parameters are fixed
     */
    MCMC(int nvoxFitpart, 
	 vector<int> bound_types, vector<T> bounds_min, vector<T> bounds_max,
	 vector<int> prior_types, vector<T> priors_a, vector<T> prior_b,
	 vector<int> fixed);

    /**
     * Run MCMC algorithm on the GPU
     * @param nvox Number of voxels in the dataset
     * @param nmeas Number of measurements in the dataset
     * @param CFP_size Size (x M measurements) of the common fixed parameters of the model
     * @param FixP_size Size (x N voxels) of the fixed parameters of the model
     * @param meas Measurements of all the voxels of the dataset (on GPU)
     * @param params Value of the parameters to estimate of the model for all the voxels of the dataset (on GPU)
     * @param CFP Common (to all the voxels) fixed parameters of the model (on GPU). CFP_size*nmeas
     * @param FixP Fixed parameters of the model (on GPU). FixP_size*nvoxels
     * @param samples Samples of the parameters estimation will be stored here (on the GPU)
     */
    void run( int nvox, int nmeas, 
	      int CFP_size, int FixP_size,
	      T* meas, T* params, 
	      T* CFP, T* FixP,
	      T* samples, T* tau);
  };
}

#endif
