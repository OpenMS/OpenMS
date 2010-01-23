// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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

namespace OpenMS
{
	int checkCUDAError(const char *msg);

	void getExternalCudaTransforms (dim3 dimGrid, dim3 dimBlock, float* positions_dev, float* intensities_dev, int from_max_to_left, int from_max_to_right, float* result_dev, 
		const int charge, const int to_load, const int to_compute, const int size, float* fwd2);
	
	int sortOnDevice(float *array, int* pos_indices, int numElements, int padding);

	void scoreOnDevice (int* sorted_positions_indices, float* trans_intensities,  float* pos, float* scores, 
		const int c, const int num_of_scores, const int overall_size, const unsigned int max_peak_cutoff, const float ampl_cutoff);

	void deriveOnDevice (float* spec, float* spec_pos, float* fwd, const int size, float* intensities_dev);
}

#endif
