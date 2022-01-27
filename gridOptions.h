#define SIZE_PART 12800 	// Number of voxels to compute in each subpart of the volume. 
				// It uses about 500MB allocated in global memory
				// Could increase it, but almost the same performance
				// Enough to exploit parallelism and low memory requirements for old/laptop GPUs

#define MAX_VOXELS_BLOCK 8      // Maximum number of voxels per Block in all the fitting routines. If any rotine does not use this number, the it should use a lower and multiple of this one. The reason is because the tool allocate the voxels data for all the routines, and the allocated memory must be a multiple of VOXELS_BLOCK, so some empty voxels may be added to avoid bad memory accesses.

// If Voxels per block are 8, in total there will be 12800/8 = 1600 blocks.




