#ifndef CUDIMOT_PARAMETERS_H_INCLUDED
#define CUDIMOT_PARAMETERS_H_INCLUDED

/**
 *
 * \class Parameters
 *
 * \brief A class for managing the Parameters of a Model
 *
 * This class contains the value of the estimated parameters of a diffusion MRI model. It also contains the value of fixed parameters of the model.
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
#include "Model.h"
#include "dMRI_Data.h"
#include "getPredictedSignal.h"
#include "BIC_AIC.h"
#include "modelparameters.h"

namespace Cudimot{

  template <typename T>
  class Parameters{
    
  private:

    /**
     * Number of parameters in the model (to estimate)
     */
    int nparams;

    ////////////////////////////////////////////////////////////////////////
    /// FIXED PARAMETERS different for each voxel, nifti volumes (ex.T1) ///
    ////////////////////////////////////////////////////////////////////////

    /**
     * Number of Fixed Parameters in the model (different for each voxel)
     */
    int nFixP;

    /**
     * Total size of volumes of Fixed Parameter
     */
    int FixP_Tsize;

    /**
     * Value of the fixed parameters of the voxels of a single part
     */
    T* FixP_host;

    /**
     * Value of the fixed parameters of the voxels of a single part allocated on the GPU
     */
    T* FixP_gpu;
    
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
     * Value of the common (to all voxels) fixed parameters
     * CFP1_meas1, CFP2_meas1, CFP3_meas1 ... CFP1_meas2, CFP2_meas2 ...
     */
    T* CFP_host;

    /**
     * Value of the common (to all voxels) fixed parameters allocated on the GPU
     */
    T* CFP_gpu;

    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    
    /**
     * Number of Voxels included in the data
     */
    int nvox;

    /**
     * Number of diffusion-weighted measurements in the data
     */
    int nmeas;

    /**
     * The data is divided into parts before beeing processd on the GPU
     */
    int nparts;

    /**
     * Number of voxels of each part of the data
     */
    int size_part;
 
    /**
     * Number of voxels of the last part of the data
     */
    int size_last_part;

    /**
     * The number of voxels in a part can be a non-multiple of voxels per block, so some threads could access to non-allocated memory. We use the closest upper multiple. The added voxels will be ignored.
     */
    int nvoxFit_part;
    
    /**
     * Parameter (to estimate) values of all the voxels.
     * Voxel 0 [p0,p1,...], Voxel1 [p0,p1,...], etc...
     */
    T* params_host;
    
    /**
     * Parameter (to estimate) values of the voxels in a single part allocated on the GPU
     */
    T* params_gpu;

    // For MCMC step
    
    /**
     * Value of the samples recorded during MCMC in a single part
     */
    T* samples_host;

    /**
     * Value of the samples recorded during MCMC in a single part and allocated on the GPU
     */
    T* samples_gpu;

     /**
     * Number of samples recorded during MCMC
     */
    int nsamples;
    
    /**
     * Predicted Signal by the model (on the host)
     */
    T* predSignal_host;

    /**
     * Predicted Signal by the model (on the gpu)
     */
    T* predSignal_gpu;

    /**
     * Class with the method to get the Predicted Signal by the model (on the gpu)
     */
    getPredictedSignal<T> PredictedSignal;

    /**
     * Bayesian information criterion (on the host)
     */
    T* BIC_host;

    /**
     * Bayesian information criterion (on the gpu)
     */
    T* BIC_gpu;

    /**
     * Class with the method to get the Bayesian & Akaike information criteria (on the gpu)
     */
    BIC_AIC<T> bic_aic;

    /**
     * Bayesian information criterion (on the host)
     */
    T* AIC_host;

    /**
     * Bayesian information criterion (on the gpu)
     */
    T* AIC_gpu;

   
    /**
     * Tau. If rician noise, tau is is 1/sigma with sigma the scale parameter. Values for each voxel/sample on the host.
     */
    T* tau_samples_host;
    
    /**
     * Tau. If rician noise, tau is is 1/sigma with sigma the scale parameter. Values for each voxel/sample on the GPU.
     */
    T* tau_samples_gpu;

    
  public:

    /**
     * Constructor
     * @param model An instance of the class Model with information about the dMRI model
     * @param dMRI_data An instance of the class dMRI_data with information about the dMRI data to be analyzed
     */
    Parameters(Model<T> model,dMRI_Data<T> dMRI_data);    

    /**
     * Destructor
     */
    ~Parameters();

    /**
     * @param part A number to identify a part of the data
     * @return A pointer to the estimated parameter values of a part of the data (on the GPU)
     */
    T* getParametersPart(int part);

    /**
     * @return A pointer to the samples recorded during MCMC (on the GPU)
     */
    T* getSamples();

    /**
     * @return Total size of Common Fixed Parameters (without counting measurements)
     */
    int getTsize_CFP();

    /**
     * @return Total size of Fixed Parameters (without counting voxels)
     */
    int getTsize_FixP();

    /**
     * @return A pointer to the Common Fixed Parameters (on the GPU)
     */
    T* getCFP();

    /**
     * @return A pointer to the Fixed Parameters for a part of the data (on the GPU)
     */
    T* getFixP_part(int part);
    
    /**
     * Copies the value of the estimated parameters of a part from GPU to the host array with all the parameter values (at its correct position). 
     * @param part A number to identify a part of the data
     */
    void copyParamsPartGPU2Host(int part);

    /**
     * Makes the pointer to the samples to point to the value of the estimated parameters (used id MCMC is not run, i.e. only one sample per parameter)
     */
    void copyParams2Samples();

    /**
     * Copies the samples of the parameters of a part from GPU to the host array with all the samples values (at its correct position)
     * @param part A number to identify a part of the data
     */
    void copySamplesPartGPU2Host(int part);

    /**
     * Writes to a binary file the samples of the parameters, Including tau (rician noise), the predicted signal, the BIC and AIC if requested.
     */
    void writeSamples();

    /**
     * @return A pointer to the tau parameter samples of each voxel (on the GPU)
     */ 
    T* getTauSamples();

    /**
     * If the user ask for the predicted signal or BIC/AIC, calculates the predicted signal or/and the BIC/AIC for this part.
     * @param mode 0: from GridSearch or LevMar (1 sample), from MCMC (several samples, needs to calculate the mean)
     * @param part A number to identify a part of the data
     * @param meas Pointer to Data measurements on the GPU
     */
    void calculate_predictedSignal_BIC_AIC(int mode, int part,T* meas);
  };
}

#endif
