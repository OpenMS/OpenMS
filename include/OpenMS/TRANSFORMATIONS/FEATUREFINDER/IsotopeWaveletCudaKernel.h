// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Rene Hussong$
// --------------------------------------------------------------------------

//**************************
//Uses the sorting code provided by Alan Kaatz
//**************************


#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELETCUDAKERNEL_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELETCUDAKERNEL_H

#include <string>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cutil.h>


#ifndef LAMBDA_Q_0
//Quadratic Fit (maybe a overkill, since the dependency looks quite linear, at least in moderate mass ranges)
#define LAMBDA_Q_0 -0.137152573151174711f
#define LAMBDA_Q_1 0.851289601785403817e-3f
#define LAMBDA_Q_2 -0.2834469691e-7f
#define BLOCK_SIZE_MAX 512
#define REGISTERS_IN_USE 13
#define ONEOLOG2E 0.6931471806f
#endif

#ifndef NEUTRON_MASS
#define NEUTRON_MASS 1.00866491578f 
#define HALF_NEUTRON_MASS 0.5043325f
#define QUARTER_NEUTRON_MASS 0.252166228f
#define WAVELET_PERIODICITY 6.229209734f
#endif

#ifndef CUDA_SORT_PARAMS
#define CUDA_SORT_PARAMS
#define MAX_BLOCKS_PER_GRID 65535

// Shared memory sort kernel
#define ELEMENTS_SORT   512
#define THREADS_SORT    (ELEMENTS_SORT >> 2)
#define SORT_NUM        0xFFFFFE00

// Shared memory merge kernel
#define ELEMENTS_MERGE  1024
#define MERGE_NUM       0xFFFFFC00
#define THREADS_MERGE   (ELEMENTS_MERGE >> 2)

// Global memory merge kernel
#define ELEMENTS_GL     2
#define THREADS_GL      256

// Minimum size sortOnDevice() can handle
#define MIN_SORT_SIZE   ELEMENTS_SORT
#endif

using namespace std;

namespace OpenMS
{

	int checkCUDAError(const char *msg);

	void getExternalCudaTransforms (dim3 dimGrid, dim3 dimBlock, float* positions_dev, float* intensities_dev, int from_max_to_left, int from_max_to_right, float* result_dev, 
		const int charge, const int block_size, const float peak_cutoff_intercept, const float peak_cutoff_slope);
	
	int sortOnDevice(float *array, int* pos_indices, int numElements, int padding);

	void scoreOnDevice (int* sorted_positions_indices, float* trans_intensities,  float* pos, float* scores, 
		const int c, const int num_of_scores, const int overall_size, const float peak_cutoff_intercept, const float peak_cutoff_slope, const unsigned int max_peak_cutoff);

}

#endif
