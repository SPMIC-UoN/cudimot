#include "modelparameters.h"

///// Edit the final part {} /////
// You need to provide the size of each Common Fixed Parameter: [Data_Measurements x K]. Only K must be provided. In this example there are 2 CFP, bvecs with size [Data_Measurements x 3] and bvals with size [Data_Measurements x 1].
int MODEL::CFP_size[] = {3,1};

// You also need to provide the size of each Fixed Parameter: [Number_Voxels x K] and only K (size of 4th dimension of the volume) must be provided.
int MODEL::FixP_size[] = {};
//////////////////////////////////

