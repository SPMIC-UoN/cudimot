#ifndef CUDIMOT_GET_PREDICTED_SIGNAL_H_INCLUDED
#define CUDIMOT_GET_PREDICTED_SIGNAL_H_INCLUDED

/**
 *
 * \class getPredictedSignal
 *
 * \brief A class for calculating the signal predicted by a model after been fitted on the GPU.
 *
 * \author Moises Hernandez-Fernandez - FMRIB Image Analysis Group
 *
 * \date March 2017
 *
 *
 *  Copyright (C) 2005 University of Oxford 
 */

#include <iostream>
#include "gridOptions.h"
#include "checkcudacalls.h"

using namespace std;

namespace Cudimot{
  template <typename T>
  class getPredictedSignal{

  private:
     
  public:
    
    /**
     * Constructor
     */
    getPredictedSignal();

    /**
     * Run the method to get the model predicted signal on the GPU
     * @param nvox Number of voxels in the dataset
     * @param nmeas Number of measurements in the dataset
     * @param nsamples Number of samples (can be several if MCMC is used). In this case, the mean of the samples is used.
     * @param CFP_size Size (without measurements) of the common fixed parameters of the model
     * @param FixP_size Size (x N voxels) of the fixed parameters of the model
     * @param samples Value of the estimated parameters
     * @param CFP Common (to all the voxels) fixed parameters of the model (on GPU). CFP_size*nmeas
     * @param FixP Fixed parameters of the model (on GPU). FixP_size*nvoxels
     * @param PredictedSignal The result values of the model predicted signal (on GPU)
     */
    void run( int nvox, int nmeas, int nsamples, int CFP_size, int FixP_size,
	      T* samples,T* CFP, T* FixP, T* PredictedSignal);
  };
}

#endif
