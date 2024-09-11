#ifndef CUDIMOT_BIC_AIC_H_INCLUDED
#define CUDIMOT_BIC_AIC_H_INCLUDED

/**
 *
 * \class BIC_AIC
 *
 * \brief A class for calculating the Bayesian & Akaike Information Criteria on the GPU
 *
 * \author Moises Hernandez-Fernandez - FMRIB Image Analysis Group
 *
 * \date May 2017
 *
 *
 *  Copyright (C) 2005 University of Oxford 
 */

#include "gridOptions.h"
#include "cudimot.h"
#include "checkcudacalls.h"
#include "cudimotoptions.h"

namespace Cudimot{
  
  template <typename T>
  class BIC_AIC{
    
  private:
    
    /** 
     * The number of voxels to fit in a part
     */
    int nvoxFit_part;
    
    /** 
     * If activated, use Rician noise modelling
     */
    bool RicianNoise;
            
  public:

    /**
     * Constructor
     * @param nvoxFitpart Number of voxel of the data to fit
     */
    BIC_AIC(int nvoxFitpart);

    /**
     * Calculate Bayesian Information Criterion
     * @param nvox Number of voxels in the dataset
     * @param nmeas Number of measurements in the dataset
     * @param nsamples Number of samples (can be several if MCMC is used). In this case, the mean of the samples is used.
     * @param CFP_size Size (x M measurements) of the common fixed parameters of the model
     * @param FixP_size Size (x N voxels) of the fixed parameters of the model
     * @param meas Measurements of all the voxels of the dataset (on GPU)
     * @param samples Value of the estimated parameters
     * @param CFP Common (to all the voxels) fixed parameters of the model (on GPU). CFP_size*nmeas
     * @param FixP Fixed parameters of the model (on GPU). FixP_size*nvoxels
     * @param BIC calculated BIC will be stored here (on the GPU)
     * @param AIC calculated AIC will be stored here (on the GPU)
     * @param tau if Rician, it contains the calue of tau for all the voxels
     */
    void run( int nvox, int nmeas, int nsamples,
	      int CFP_size, int FixP_size,
	      T* meas, T* samples, 
	      T* CFP, T* FixP,
	      T* BIC, T* AIC, T* tau);
  };
}

#endif
