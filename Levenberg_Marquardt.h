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
