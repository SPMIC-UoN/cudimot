#ifndef CUDIMOT_dMRI_DATA_H_INCLUDED
#define CUDIMOT_dMRI_DATA_H_INCLUDED

/**
 *
 * \class dMRI_Data
 *
 * \brief A class for managing the dMRI Data to analyse
 *
 * This class contains the diffusion MRI data of a group of voxels for beeing analysed. A model will be fitted to this data on the GPU. A large number of voxels (>10000) should be used for exploiting the GPU.
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
#include "newmat.h"
#include "newimage/newimageall.h"
#include "checkcudacalls.h"
#include "gridOptions.h"
#include "cudimotoptions.h"

namespace Cudimot{

  template <typename T>
  class Parameters;

  template <typename T>
  class dMRI_Data{
    friend class Parameters<T>;

  private:
    /**
     * Number of diffusion-weighted measurements in the data
     */
    int nmeas;

    /**
     * Number of Voxels included in the data
     */
    int nvox;

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
     * Measurements for all the voxels of the data
     */
    Matrix dataM;

    /**
     * The number of voxels in a part can be a non-multiple of voxels per block, so some threads could access to non-allocated memory. We use the closest upper multiple. The added voxels will be ignored.
     */
    int nvoxFit_part;

    /**
     * Measurements of the voxels in a single part.
     * Voxel 0 [meas0, meas 1, ...], Voxel1 [meas0, meas 1, ...], etc...
     */
    T* meas_host;
    
    /**
     * Measurements of the voxels in a single part allocated on the GPU
     */
    T* meas_gpu;
    
    /**
     * Method to remove the negative measurements of the data of one voxel
     * @param Voxdata A vector with the measurements of one voxel
     */
    void remove_NonPositive_entries(NEWMAT::ColumnVector& Voxdata);
    
  public:
    /**
     * Constructor
     */
    dMRI_Data();
    
    /**
     * Destructor
     */
    ~dMRI_Data();

    /**
     * @return The number of voxels to fit in each part of the data
     */
    int getNvoxFit_part() const;
    
    /**
     * @return The number of measurements of the data (all the voxels have the same number of measurements)
     */
    int getNmeas() const;
    
    /**
     * @return The number of parts of the data
     */
    int getNparts() const;
    
    /**
     * @param part A number to identify a part of the data
     * @param part_size The method returns here the size of the part
     * @return The measurements of a part of the data (on the GPU)
    */
    T* getMeasPart(int part, int &part_size);
  };
}

#endif
