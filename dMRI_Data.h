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

#include <vector>
#include "armawrap/newmat.h"
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
    NEWMAT::Matrix dataM;

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
