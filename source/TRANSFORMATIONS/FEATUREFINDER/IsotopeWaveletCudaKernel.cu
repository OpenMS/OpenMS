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
// $Maintainer: Rene Hussong $
// $Authors: Rene Hussong $
// --------------------------------------------------------------------------


//**************************
//Uses the sorting code provided by Alan Kaatz
//http://courses.ece.uiuc.edu/ece498/al1/mps/MP5-TopWinners/kaatz/MP5-parallel_sort.zip
//**************************

using namespace std;

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletCudaKernel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletConstants.h>

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <iomanip>

#include <assert.h>
#include <cuda.h>


texture<float,1> trans_intensities_tex, pos_tex, int_tex;
texture<int, 1> sorted_positions_indices_tex;

namespace OpenMS
{
	
	int checkCUDAError(const char *msg)
	{
			cudaError_t err = cudaGetLastError();
			if( cudaSuccess != err) 
			{
					fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
					return (-1);
			};     
			return (0);       
	}

	
	__device__ float isotope_wavelet (float tz1, float mz)
	{
		//float fac (-(Constants::LAMBDA_Q_0 + Constants::LAMBDA_Q_1*mz + Constants::LAMBDA_Q_2*mz*mz));
		float fac (-(Constants::LAMBDA_L_0 + Constants::LAMBDA_L_1*mz));

		fac += (tz1-1)*__log2f(-fac)*Constants::ONEOLOG2E - lgammaf(tz1);
			
		return (__sinf((tz1-1)*Constants::WAVELET_PERIODICITY) * __expf(fac));
	}

	__device__ unsigned int getMzPeakCutOffAtMonoPos (float mass, unsigned int  z)
	{
		float m (mass*z);
		return ( m<Constants::BORDER_MZ_FIT99 ? 
			ceil((Constants::CUTOFF_FIT99_POLY_0+Constants::CUTOFF_FIT99_POLY_1*m+Constants::CUTOFF_FIT99_POLY_2*m*m))
				: ceil((Constants::CUTOFF_FIT99_POLY_3+Constants::CUTOFF_FIT99_POLY_4*m+Constants::CUTOFF_FIT99_POLY_5*m*m)));
	}
	
	__device__ unsigned int getNumPeakCutOff (float mass, unsigned int  z)
	{
		float m (mass*z);
		return ( m<Constants::BORDER_MZ_FIT99 ? 
			ceil((Constants::CUTOFF_FIT99_POLY_0+Constants::CUTOFF_FIT99_POLY_1*m+Constants::CUTOFF_FIT99_POLY_2*m*m-Constants::IW_QUARTER_NEUTRON_MASS))
				: ceil((Constants::CUTOFF_FIT99_POLY_3+Constants::CUTOFF_FIT99_POLY_4*m+Constants::CUTOFF_FIT99_POLY_5*m*m-Constants::IW_QUARTER_NEUTRON_MASS)));
	}

	__device__ unsigned int getNumPeakCutOff (float m)
	{
		return ( m<Constants::BORDER_MZ_FIT99 ? 
			ceil((Constants::CUTOFF_FIT99_POLY_0+Constants::CUTOFF_FIT99_POLY_1*m+Constants::CUTOFF_FIT99_POLY_2*m*m-Constants::IW_QUARTER_NEUTRON_MASS))
				: ceil((Constants::CUTOFF_FIT99_POLY_3+Constants::CUTOFF_FIT99_POLY_4*m+Constants::CUTOFF_FIT99_POLY_5*m*m-Constants::IW_QUARTER_NEUTRON_MASS)));
	}


	__global__ void ConvolutionIsotopeWaveletKernel(float* signal_pos, float* signal_int, const int from_max_to_left, const int from_max_to_right, float* result, 
		const unsigned int charge, const int to_load, const int to_compute, const int size)
	{
		// the device-shared memory storing one data block
		// This is currently hard-coded to 256 points, since we require two 4B floats for each
		// data point, leading to 2kB per block.
		__shared__ float signal_pos_block[Constants::CUDA_EXTENDED_BLOCK_SIZE_MAX];//[BLOCK_SIZE_MAX];
		__shared__ float signal_int_block[Constants::CUDA_EXTENDED_BLOCK_SIZE_MAX];//[BLOCK_SIZE_MAX];

		// load the data from device memory to shared memory. 
		// to distribute the loads as evenly as possible over the threads, each thread loads
		// the data point it will later compute in the output. the first wavelet_length threads
		// will also load the padding to the left of the signal, the last wavelet_length ones will 
		// load the padding to the right
		
		// we will silently ignore the first wavelet_length points in the output; these have to be
		// zero-padded by the calling function. our data organization is as follows: each block computes
		// a part of the output that is block_size-2*wavelet_length points long. For the computation, we
		// require wavelet_length points on the left and on the right so we can put the wavelet on all
		// points even at the boundary.
		//                 left padding,                                                      position of thread
		//                 ignored in output    the points computed by the previous blocks    in block
		int my_data_pos  = from_max_to_left    +  blockIdx.x*to_compute + threadIdx.x;

		int my_local_pos = threadIdx.x + from_max_to_left;

		//every thread with an ID smaller than the number of from_max_to_left loads the additional boundary points
		//at the left end
		if (threadIdx.x < from_max_to_left)
		{
			signal_pos_block[threadIdx.x] = signal_pos[my_data_pos-from_max_to_left];
			signal_int_block[threadIdx.x] = signal_int[my_data_pos-from_max_to_left];
			//printf ("PreLoading: %i\t%f\n", threadIdx.x, signal_pos_block[threadIdx.x]);
		}
			
		int additional_right_end_loads=0;
		while (my_local_pos + (additional_right_end_loads)*Constants::CUDA_BLOCK_SIZE_MAX < to_load)
		{
			signal_pos_block[my_local_pos+additional_right_end_loads*Constants::CUDA_BLOCK_SIZE_MAX] = signal_pos[my_data_pos+additional_right_end_loads*(Constants::CUDA_BLOCK_SIZE_MAX)];
			signal_int_block[my_local_pos+additional_right_end_loads*Constants::CUDA_BLOCK_SIZE_MAX] = signal_int[my_data_pos+additional_right_end_loads*(Constants::CUDA_BLOCK_SIZE_MAX)];
			//printf ("Loading: %i\t%f\n", my_local_pos+additional_right_end_loads*BLOCK_SIZE_MAX, signal_pos_block[my_local_pos+additional_right_end_loads*BLOCK_SIZE_MAX]);
			++additional_right_end_loads;
		};
		
		//wait until the shared data is loaded completely
		__syncthreads(); 

		if (threadIdx.x >= to_compute || 	my_data_pos - from_max_to_left >= size)
			return;

		float value = 0, boundary = getMzPeakCutOffAtMonoPos(signal_pos_block[my_local_pos], charge)/charge;
		float old=0, c_diff, current, old_pos = (my_local_pos-from_max_to_left-1) > 0 ? signal_pos_block[my_local_pos-from_max_to_left-1] 
			: signal_pos_block[my_local_pos-from_max_to_left]-(signal_pos[size-1]-signal_pos[size-2]); //i.e. min_spacing

		//printf ("*********************************************** %f\n",  signal_pos_block[my_local_pos]);

		for (int current_conv_pos = my_local_pos-from_max_to_left; 
						current_conv_pos < my_local_pos+from_max_to_right; 
							++current_conv_pos)
		{
			c_diff = signal_pos_block[current_conv_pos]-signal_pos_block[my_local_pos]+Constants::IW_QUARTER_NEUTRON_MASS/(float)charge;

			//Attention! The +1. has nothing to do with the charge, it is caused by the wavelet's formula (tz1).
			current = c_diff > 0 && c_diff <= boundary ? isotope_wavelet(c_diff*charge+1., signal_pos_block[current_conv_pos]*charge)*signal_int_block[current_conv_pos] : 0;
			
			value += 0.5*(current + old)*(signal_pos_block[current_conv_pos]-old_pos);
			
			#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
				
				//if (trunc(signal_pos_block[my_local_pos]) == 556 && charge == 2)
				{
					//printf ("%i \t %f \t %f  \t %f \t %f \t %f \t %f \t %f\n", current_conv_pos, signal_pos_block[current_conv_pos], c_diff, current, old, c_diff, signal_int_block[current_conv_pos], boundary);
					printf ("%f  \t %f \t %f \t %f \t %f \t %f\n", signal_pos_block[current_conv_pos], c_diff > 0 && c_diff <= boundary ? isotope_wavelet(c_diff*charge+1., signal_pos_block[current_conv_pos]*charge): 0, value, signal_int_block[current_conv_pos], c_diff, boundary);

				};
			#endif

			old = current;
			old_pos = signal_pos_block[current_conv_pos];
		};
		
		//printf ("################################################\n");

		result[my_data_pos] = value;
	}


	__global__ void ConvolutionIsotopeWaveletKernelTexture(const int from_max_to_left, const int from_max_to_right, float* result, 
		const unsigned int charge, const int size)
	{
		int my_data_pos  = from_max_to_left    +  blockIdx.x*blockDim.x + threadIdx.x;

		if (my_data_pos - from_max_to_left >= size)
			return;

		float value = 0, boundary = getMzPeakCutOffAtMonoPos(tex1Dfetch(pos_tex, my_data_pos), charge)/charge;
		float old=0, c_diff, current, old_pos = (my_data_pos-from_max_to_left-1) > 0 ? tex1Dfetch(pos_tex, my_data_pos-from_max_to_left-1) 
				: tex1Dfetch(pos_tex,my_data_pos-from_max_to_left)-(tex1Dfetch(pos_tex, size-1)-tex1Dfetch(pos_tex, size-2)); //i.e. min_spacing

		for (int current_conv_pos = my_data_pos-from_max_to_left; 
						current_conv_pos < my_data_pos+from_max_to_right; 
							++current_conv_pos)
		{
			c_diff =  tex1Dfetch(pos_tex, current_conv_pos)- tex1Dfetch(pos_tex, my_data_pos)+Constants::IW_QUARTER_NEUTRON_MASS/(float)charge;

			//Attention! The +1. has nothing to do with the charge, it is caused by the wavelet's formula (tz1).
			current = c_diff > 0 && c_diff <= boundary ? isotope_wavelet(c_diff*charge+1., tex1Dfetch(pos_tex, current_conv_pos)*charge)*tex1Dfetch(int_tex, current_conv_pos) : 0;
			value += 0.5*(current + old)*(tex1Dfetch(pos_tex, current_conv_pos)-old_pos);
			
			old = current;
			old_pos = tex1Dfetch(pos_tex, current_conv_pos);
		};

		result[my_data_pos] = value;
	}


	__global__ void getDerivatives (float* spec, float* spec_pos, float* fwd2, const int size, float* intensities_dev)
	{
		int i = threadIdx.x + blockIdx.x * blockDim.x;
	
		if ((i+2>=size && i<size) || i==0)
		{
			fwd2[i] = 0;
			return;
		};
	
		float share = spec[i+1], share_pos = spec_pos[i+1];
		float bwd = (share-spec[i])/(share_pos-spec_pos[i]);
		float fwd = (spec[i+2]-share)/(spec_pos[i+2]-share_pos);

		if (bwd>=0 && fwd<=0 && share <= intensities_dev[i+1])
		{
			fwd2[i+1] = spec[i+1];
		};
	}


	void deriveOnDevice (float* spec, float* spec_pos, float* fwd, const int size, float* intensities_dev)
	{
		dim3 blockDim (Constants::CUDA_BLOCK_SIZE_MAX);
		dim3 gridDim ((int)(ceil)(size/(float)Constants::CUDA_BLOCK_SIZE_MAX));
		getDerivatives<<<gridDim, blockDim>>> (spec, spec_pos, fwd, size, intensities_dev);
		cudaThreadSynchronize();
		checkCUDAError("deriveOnDevice");
	}


	void getExternalCudaTransforms (dim3 dimGrid, dim3 dimBlock, float* positions_dev, float* intensities_dev, int from_max_to_left, int from_max_to_right, float* result_dev, 
		const int charge, const int to_load, const int to_compute, const int size, float* fwd2) 
	{
		if (to_load < Constants::CUDA_EXTENDED_BLOCK_SIZE_MAX)
		{	
			ConvolutionIsotopeWaveletKernel<<<dimGrid,dimBlock>>> (positions_dev, intensities_dev, from_max_to_left, from_max_to_right, result_dev, charge, to_load, to_compute, size);
			cudaThreadSynchronize();
			checkCUDAError("ConvolutionIsotopeWaveletKernel");
			deriveOnDevice (result_dev, positions_dev, fwd2, size, intensities_dev);
		}
		else
		{
			dimBlock = dim3(Constants::CUDA_TEXTURE_THREAD_LIMIT);
			dimGrid = dim3((int)ceil(size/(float)dimBlock.x));
			cudaBindTexture(0, int_tex, intensities_dev, (size+from_max_to_left+from_max_to_right)*sizeof(float));
			cudaBindTexture(0, pos_tex, positions_dev, (size+from_max_to_left+from_max_to_right)*sizeof(float));

			ConvolutionIsotopeWaveletKernelTexture<<<dimGrid,dimBlock>>> (from_max_to_left, from_max_to_right, result_dev, charge, size);
			cudaThreadSynchronize();
			checkCUDAError("ConvolutionIsotopeWaveletKernelTexture");
			deriveOnDevice (result_dev, positions_dev, fwd2, size, intensities_dev);
			
			cudaUnbindTexture(int_tex);
			cudaUnbindTexture(pos_tex);
		};
	}


	__device__ inline void swap(float &a, float &b, int &c, int &d) 
	{
			float tmp (a);
			a = b;
			b = tmp;
				
			int tmp2 (c);
			c = d;
			d = tmp2;
	}

	__global__ void sharedMemMerge(float *array, int *pos, int k) {

			__shared__ float shmem[Constants::CUDA_ELEMENTS_MERGE];
			__shared__ int posshmem[Constants::CUDA_ELEMENTS_MERGE];

			int tmp = blockIdx.x * Constants::CUDA_ELEMENTS_MERGE + threadIdx.x;

			float data = array[tmp];
			float data2 = array[tmp + (Constants::CUDA_ELEMENTS_MERGE / 2)];

			float data3 = array[tmp + Constants::CUDA_THREADS_MERGE];
			float data4 = array[tmp + Constants::CUDA_THREADS_MERGE + (Constants::CUDA_ELEMENTS_MERGE / 2)];
			
			int posdata = pos[tmp];
			int posdata2 = pos[tmp + (Constants::CUDA_ELEMENTS_MERGE / 2)];

			int posdata3 = pos[tmp + Constants::CUDA_THREADS_MERGE];
			int posdata4 = pos[tmp + Constants::CUDA_THREADS_MERGE + (Constants::CUDA_ELEMENTS_MERGE / 2)];

			int dir = k & (blockIdx.x * (Constants::CUDA_ELEMENTS_MERGE));


			if (dir == 0) {
					if (data > data2) {  // ascending
							shmem[threadIdx.x] = data2;
							shmem[threadIdx.x + (Constants::CUDA_ELEMENTS_MERGE / 2)] = data;
							posshmem[threadIdx.x] = posdata2;
							posshmem[threadIdx.x + (Constants::CUDA_ELEMENTS_MERGE / 2)] = posdata;
					} else {
							shmem[threadIdx.x] = data;
							shmem[threadIdx.x + (Constants::CUDA_ELEMENTS_MERGE / 2)] = data2;
							posshmem[threadIdx.x] = posdata;
							posshmem[threadIdx.x + (Constants::CUDA_ELEMENTS_MERGE / 2)] = posdata2;
					}

					if (data3 > data4) {  // ascending
							shmem[threadIdx.x + Constants::CUDA_THREADS_MERGE] = data4;
							shmem[threadIdx.x + (Constants::CUDA_ELEMENTS_MERGE / 2) + Constants::CUDA_THREADS_MERGE] = data3;
							posshmem[threadIdx.x + Constants::CUDA_THREADS_MERGE] = posdata4;
							posshmem[threadIdx.x + (Constants::CUDA_ELEMENTS_MERGE / 2) + Constants::CUDA_THREADS_MERGE] = posdata3;
					} else {
							shmem[threadIdx.x + Constants::CUDA_THREADS_MERGE] = data3;
							shmem[threadIdx.x + (Constants::CUDA_ELEMENTS_MERGE / 2) + Constants::CUDA_THREADS_MERGE] = data4;
							posshmem[threadIdx.x + Constants::CUDA_THREADS_MERGE] = posdata3;
							posshmem[threadIdx.x + (Constants::CUDA_ELEMENTS_MERGE / 2) + Constants::CUDA_THREADS_MERGE] = posdata4;
					}
			} else {
					if (data < data2) {  // descending
							shmem[threadIdx.x] = data2;
							shmem[threadIdx.x + (Constants::CUDA_ELEMENTS_MERGE / 2)] = data;
							posshmem[threadIdx.x] = posdata2;
							posshmem[threadIdx.x + (Constants::CUDA_ELEMENTS_MERGE / 2)] = posdata;
					} else {
							shmem[threadIdx.x] = data;
							shmem[threadIdx.x + (Constants::CUDA_ELEMENTS_MERGE / 2)] = data2;							
							posshmem[threadIdx.x] = posdata;
							posshmem[threadIdx.x + (Constants::CUDA_ELEMENTS_MERGE / 2)] = posdata2;
					}

					if (data3 < data4) {  // descending
							shmem[threadIdx.x + Constants::CUDA_THREADS_MERGE] = data4;
							shmem[threadIdx.x + (Constants::CUDA_ELEMENTS_MERGE / 2) + Constants::CUDA_THREADS_MERGE] = data3;
							posshmem[threadIdx.x + Constants::CUDA_THREADS_MERGE] = posdata4;
							posshmem[threadIdx.x + (Constants::CUDA_ELEMENTS_MERGE / 2) + Constants::CUDA_THREADS_MERGE] = posdata3;
					} else {
							shmem[threadIdx.x + Constants::CUDA_THREADS_MERGE] = data3;
							shmem[threadIdx.x + (Constants::CUDA_ELEMENTS_MERGE / 2) + Constants::CUDA_THREADS_MERGE] = data4;							
							posshmem[threadIdx.x + Constants::CUDA_THREADS_MERGE] = posdata3;
							posshmem[threadIdx.x + (Constants::CUDA_ELEMENTS_MERGE / 2) + Constants::CUDA_THREADS_MERGE] = posdata4;
					}
			}



			int j = 256, s = Constants::CUDA_MERGE_NUM >> 2; 


			int x = threadIdx.x + (s & threadIdx.x);
			int y = x + j;
			__syncthreads();

			if (dir == 0) {
					if (shmem[x] > shmem[y]) {  // ascending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			} else {
					if (shmem[x] < shmem[y]) {  // descending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);							
					}
			}

			x = (threadIdx.x + Constants::CUDA_THREADS_MERGE) + ((threadIdx.x + Constants::CUDA_THREADS_MERGE) & s);
			y = x + j;

			if (dir == 0) {
					if (shmem[x] > shmem[y]) {  // ascending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);							
					}
			} else {
					if (shmem[x] < shmem[y]) {  // descending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);							
					}
			}

			j >>= 1;
			s >>= 1;


			x = threadIdx.x + (s & threadIdx.x);
			y = x + j;
			__syncthreads();

			if (dir == 0) {
					if (shmem[x] > shmem[y]) {  // ascending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);							
					}
			} else {
					if (shmem[x] < shmem[y]) {  // descending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			}

			x = (threadIdx.x + Constants::CUDA_THREADS_MERGE) + ((threadIdx.x + Constants::CUDA_THREADS_MERGE) & s);
			y = x + j;

			if (dir == 0) {
					if (shmem[x] > shmem[y]) {  // ascending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			} else {
					if (shmem[x] < shmem[y]) {  // descending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			}

			j >>= 1;
			s >>= 1;


			x = threadIdx.x + (s & threadIdx.x);
			y = x + j;
			__syncthreads();

			if (dir == 0) {
					if (shmem[x] > shmem[y]) {  // ascending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			} else {
					if (shmem[x] < shmem[y]) {  // descending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			}

			x = (threadIdx.x + Constants::CUDA_THREADS_MERGE) + ((threadIdx.x + Constants::CUDA_THREADS_MERGE) & s);
			y = x + j;

			if (dir == 0) {
					if (shmem[x] > shmem[y]) {  // ascending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			} else {
					if (shmem[x] < shmem[y]) {  // descending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			}
			
			j >>= 1;
			s >>= 1;


			x = threadIdx.x + (s & threadIdx.x);
			y = x + j;
			__syncthreads();

			if (dir == 0) {
					if (shmem[x] > shmem[y]) {  // ascending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			} else {
					if (shmem[x] < shmem[y]) {  // descending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			}

			x = (threadIdx.x + Constants::CUDA_THREADS_MERGE) + ((threadIdx.x + Constants::CUDA_THREADS_MERGE) & s);
			y = x + j;

			if (dir == 0) {
					if (shmem[x] > shmem[y]) {  // ascending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			} else {
					if (shmem[x] < shmem[y]) {  // descending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			}

			j >>= 1;
			s >>= 1;


			x = threadIdx.x + (s & threadIdx.x);
			y = x + j;
			__syncthreads();

			if (dir == 0) {
					if (shmem[x] > shmem[y]) {  // ascending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			} else {
					if (shmem[x] < shmem[y]) {  // descending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			}

			x = (threadIdx.x + Constants::CUDA_THREADS_MERGE) + ((threadIdx.x + Constants::CUDA_THREADS_MERGE) & s);
			y = x + j;

			if (dir == 0) {
					if (shmem[x] > shmem[y]) {  // ascending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			} else {
					if (shmem[x] < shmem[y]) {  // descending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			}

			j >>= 1;
			s >>= 1;


			x = threadIdx.x + (s & threadIdx.x);
			y = x + j;
			__syncthreads();

			if (dir == 0) {
					if (shmem[x] > shmem[y]) {  // ascending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			} else {
					if (shmem[x] < shmem[y]) {  // descending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			}

			x = (threadIdx.x + Constants::CUDA_THREADS_MERGE) + ((threadIdx.x + Constants::CUDA_THREADS_MERGE) & s);
			y = x + j;

			if (dir == 0) {
					if (shmem[x] > shmem[y]) {  // ascending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			} else {
					if (shmem[x] < shmem[y]) {  // descending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			}

			j >>= 1;
			s >>= 1;


			x = threadIdx.x + (s & threadIdx.x);
			y = x + j;
			__syncthreads();

			if (dir == 0) {
					if (shmem[x] > shmem[y]) {  // ascending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			} else {
					if (shmem[x] < shmem[y]) {  // descending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			}

			x = (threadIdx.x + Constants::CUDA_THREADS_MERGE) + ((threadIdx.x + Constants::CUDA_THREADS_MERGE) & s);
			y = x + j;

			if (dir == 0) {
					if (shmem[x] > shmem[y]) {  // ascending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			} else {
					if (shmem[x] < shmem[y]) {  // descending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			}

			j >>= 1;
			s >>= 1;


			x = threadIdx.x + (s & threadIdx.x);
			y = x + j;
			__syncthreads();

			if (dir == 0) {
					if (shmem[x] > shmem[y]) {  // ascending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			} else {
					if (shmem[x] < shmem[y]) {  // descending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);	
					}
			}

			x = (threadIdx.x + Constants::CUDA_THREADS_MERGE) + ((threadIdx.x + Constants::CUDA_THREADS_MERGE) & s);
			y = x + j;

			if (dir == 0) {
					if (shmem[x] > shmem[y]) {  // ascending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			} else {
					if (shmem[x] < shmem[y]) {  // descending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			}

			j >>= 1;
			s >>= 1;


			x = threadIdx.x + (s & threadIdx.x);
			y = x + j;
			__syncthreads();

			if (dir == 0) {
					if (shmem[x] > shmem[y]) {  // ascending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			} else {
					if (shmem[x] < shmem[y]) {  // descending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			}

			x = (threadIdx.x + Constants::CUDA_THREADS_MERGE) + ((threadIdx.x + Constants::CUDA_THREADS_MERGE) & s);
			y = x + j;

			if (dir == 0) {
					if (shmem[x] > shmem[y]) {  // ascending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			} else {
					if (shmem[x] < shmem[y]) {  // descending
							swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
					}
			}


			__syncthreads();

			int i = blockIdx.x * Constants::CUDA_ELEMENTS_MERGE + threadIdx.x;
			array[i] = shmem[threadIdx.x];
			array[i + (Constants::CUDA_ELEMENTS_MERGE / 4)] = shmem[(Constants::CUDA_ELEMENTS_MERGE / 4) + threadIdx.x];
			array[i + (Constants::CUDA_ELEMENTS_MERGE / 2)] = shmem[(Constants::CUDA_ELEMENTS_MERGE / 2) + threadIdx.x];
			array[i + (3 * Constants::CUDA_ELEMENTS_MERGE / 4)] = shmem[(3 * Constants::CUDA_ELEMENTS_MERGE / 4) + threadIdx.x];

			pos[i] = posshmem[threadIdx.x];
			pos[i + (Constants::CUDA_ELEMENTS_MERGE / 4)] = posshmem[(Constants::CUDA_ELEMENTS_MERGE / 4) + threadIdx.x];
			pos[i + (Constants::CUDA_ELEMENTS_MERGE / 2)] = posshmem[(Constants::CUDA_ELEMENTS_MERGE / 2) + threadIdx.x];
			pos[i + (3 * Constants::CUDA_ELEMENTS_MERGE / 4)] = posshmem[(3 * Constants::CUDA_ELEMENTS_MERGE / 4) + threadIdx.x];
	}




	__global__ void mergeArray(float *array, int* pos, int j, int k, int s) {
			int tmp = (blockIdx.x * Constants::CUDA_THREADS_GL);
			int x = tmp +  threadIdx.x + (tmp & s);
			j += x;

			float data1 = array[x];
			float data2 = array[j];

			if ((x & k) == 0) {    // ascending
					if (data1 > data2) {
							swap(array[x], array[j], pos[x], pos[j]);
					}
			} else {                // descending
					if (data1 < data2) {
							swap(array[x], array[j], pos[x], pos[j]);
					}
			}
	}




	__global__ void sharedMemSort(float2 *array, int2 *pos) 
	{
			__shared__ float shmem[Constants::CUDA_ELEMENTS_SORT];
			__shared__ int posshmem[Constants::CUDA_ELEMENTS_SORT];

			float2 data = array[blockIdx.x * (Constants::CUDA_ELEMENTS_SORT / 2) + threadIdx.x];
			int2 posdata = pos[blockIdx.x * (Constants::CUDA_ELEMENTS_SORT / 2) + threadIdx.x];


			if ( (threadIdx.x & 1) == 0) {
					if (data.x > data.y) {  // ascending
							shmem[2 * threadIdx.x] = data.y;
							shmem[2 * threadIdx.x + 1] = data.x;
							posshmem[2 * threadIdx.x] = posdata.y;
							posshmem[2 * threadIdx.x + 1] = posdata.x;
					} else {
							shmem[2 * threadIdx.x] = data.x;
							shmem[2 * threadIdx.x + 1] = data.y;									
							posshmem[2 * threadIdx.x] = posdata.x;
							posshmem[2 * threadIdx.x + 1] = posdata.y;							
					}
			} else {
					if (data.x < data.y) {  // descending
							shmem[2 * threadIdx.x] = data.y;
							shmem[2 * threadIdx.x + 1] = data.x;							
							posshmem[2 * threadIdx.x] = posdata.y;
							posshmem[2 * threadIdx.x + 1] = posdata.x;

					} else {
							shmem[2 * threadIdx.x] = data.x;
							shmem[2 * threadIdx.x + 1] = data.y;							
							posshmem[2 * threadIdx.x] = posdata.x;
							posshmem[2 * threadIdx.x + 1] = posdata.y;
					}
			}


			for (int k = 4, r = 0xFFFFFFFC; k <= (Constants::CUDA_ELEMENTS_SORT / 2); k *= 2, r <<= 1) {

					for (int j = k >> 1, s = r >> 1; j > 0; j >>= 1, s >>= 1) {

							int x = threadIdx.x + (threadIdx.x & s);
							int y = x + j;
							
							__syncthreads();
							
							if ((x & k) == 0) {

									if (shmem[x] > shmem[y]) {  // ascending
											swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
									}
							} else {
									if (shmem[x] < shmem[y]) {  // descending
											swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
									}
							}
					}
			}


			data = array[blockIdx.x * (Constants::CUDA_ELEMENTS_SORT / 2) + Constants::CUDA_THREADS_SORT + threadIdx.x];
			posdata = pos[blockIdx.x * (Constants::CUDA_ELEMENTS_SORT / 2) + Constants::CUDA_THREADS_SORT + threadIdx.x];
			float* shmem2 = &shmem[Constants::CUDA_ELEMENTS_SORT / 2];
			int* posshmem2 = &posshmem[Constants::CUDA_ELEMENTS_SORT / 2];
			
			__syncthreads();

			if ( (threadIdx.x & 1) == 0) {
					if (data.x > data.y) {  // ascending
							shmem2[2 * threadIdx.x] = data.y;
							shmem2[2 * threadIdx.x + 1] = data.x;
							posshmem2[2 * threadIdx.x] = posdata.y;
							posshmem2[2 * threadIdx.x + 1] = posdata.x;
					} else {
							shmem2[2 * threadIdx.x] = data.x;
							shmem2[2 * threadIdx.x + 1] = data.y;							
							posshmem2[2 * threadIdx.x] = posdata.x;
							posshmem2[2 * threadIdx.x + 1] = posdata.y;
					}
			} else {
					if (data.x < data.y) {  // descending
							shmem2[2 * threadIdx.x] = data.y;
							shmem2[2 * threadIdx.x + 1] = data.x;
							posshmem2[2 * threadIdx.x] = posdata.y;
							posshmem2[2 * threadIdx.x + 1] = posdata.x;
					} else {
							shmem2[2 * threadIdx.x] = data.x;
							shmem2[2 * threadIdx.x + 1] = data.y;							
							posshmem2[2 * threadIdx.x] = posdata.x;
							posshmem2[2 * threadIdx.x + 1] = posdata.y;
					}
			}


			for (int k = 4, r = 0xFFFFFFFC; k <= (Constants::CUDA_ELEMENTS_SORT / 2); k *= 2, r <<= 1) {

					for (int j = k >> 1, s = r >> 1; j > 0; j >>= 1, s >>= 1) {

							int x = threadIdx.x + (threadIdx.x & s);
							int y = x + j;
							__syncthreads();

							if ((x & k) == 0) {
									if (shmem2[x] < shmem2[y]) {  // descending
											swap(shmem2[x], shmem2[y], posshmem2[x], posshmem2[y]);	
									}
							} else {
									if (shmem2[x] > shmem2[y]) {  // ascending
											swap(shmem2[x], shmem2[y], posshmem2[x], posshmem2[y]);
									}
							}
					}
			}


			if ((blockIdx.x & 1) == 0) {

					for (int j = Constants::CUDA_ELEMENTS_SORT / 2, s = Constants::CUDA_SORT_NUM >> 1; j > 0; j >>= 1, s >>= 1) {

							int x = threadIdx.x + (threadIdx.x & s);
							int y = x + j;
							__syncthreads();

							if (shmem[x] > shmem[y]) {  // ascending
									swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
							}

							x = (threadIdx.x + Constants::CUDA_THREADS_SORT) + ((threadIdx.x + Constants::CUDA_THREADS_SORT) & s);
							y = x + j;
			
							if (shmem[x] > shmem[y]) {  // ascending
									swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
							}
					}
			
			} else {
			
					for (int j = Constants::CUDA_ELEMENTS_SORT / 2, s = Constants::CUDA_SORT_NUM >> 1; j > 0; j >>= 1, s >>= 1) {

							int x = threadIdx.x + (threadIdx.x & s);
							int y = x + j;
							__syncthreads();

							if (shmem[x] < shmem[y]) {  // descending
									swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
							}

							x = (threadIdx.x + Constants::CUDA_THREADS_SORT) + ((threadIdx.x + Constants::CUDA_THREADS_SORT) & s);
							y = x + j;
			
							if (shmem[x] < shmem[y]) {  // descending
									swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
							}
					}
			
			}


			__syncthreads();

			int i = blockIdx.x * Constants::CUDA_ELEMENTS_SORT + threadIdx.x;
			((float*)array)[i] = shmem[threadIdx.x];
			((float*)array)[i + (Constants::CUDA_ELEMENTS_SORT / 4)] = shmem[(Constants::CUDA_ELEMENTS_SORT / 4) + threadIdx.x];
			((float*)array)[i + (Constants::CUDA_ELEMENTS_SORT / 2)] = shmem[(Constants::CUDA_ELEMENTS_SORT / 2) + threadIdx.x];
			((float*)array)[i + (3 * Constants::CUDA_ELEMENTS_SORT / 4)] = shmem[(3 * Constants::CUDA_ELEMENTS_SORT / 4) + threadIdx.x];

			((int*)pos)[i] = posshmem[threadIdx.x];
			((int*)pos)[i + (Constants::CUDA_ELEMENTS_SORT / 4)] = posshmem[(Constants::CUDA_ELEMENTS_SORT / 4) + threadIdx.x];
			((int*)pos)[i + (Constants::CUDA_ELEMENTS_SORT / 2)] = posshmem[(Constants::CUDA_ELEMENTS_SORT / 2) + threadIdx.x];
			((int*)pos)[i + (3 * Constants::CUDA_ELEMENTS_SORT / 4)] = posshmem[(3 * Constants::CUDA_ELEMENTS_SORT / 4) + threadIdx.x];
	}


	__global__ void findCutOffIndex (float* array, int* cut_off_index)
	{
		int my_index = blockIdx.x*blockDim.x + threadIdx.x;

		if (my_index+1 >= gridDim.x*blockDim.x)
		{
			//printf ("returning\n");
			return;
		};
		float first = array[my_index], second = array[my_index+1];
		//printf ("my_index: %i\t%f\t%f\n", my_index, first, second);

		if (first <=0 && second > 0)
		{
			//printf ("writing: %i",  my_index+1);
			*cut_off_index = my_index+1;
		};
	};


	int sortOnDevice(float *array, int* pos_indices, int numElements, int padding)
	{
		dim3 dimGridSharedMemSort((numElements / Constants::CUDA_ELEMENTS_SORT) - (padding / Constants::CUDA_ELEMENTS_SORT), 1, 1);
    dim3 dimBlockSharedMemSort(Constants::CUDA_THREADS_SORT, 1, 1);

    dim3 dimGridMergeArray(numElements / (Constants::CUDA_THREADS_GL * Constants::CUDA_ELEMENTS_GL), 1, 1);
    dim3 dimBlockMergeArray(Constants::CUDA_THREADS_GL, 1, 1);

    dim3 dimGridSharedMemMerge(numElements / Constants::CUDA_ELEMENTS_MERGE, 1, 1);
    dim3 dimBlockSharedMemMerge(Constants::CUDA_THREADS_MERGE, 1, 1);

		sharedMemSort<<<dimGridSharedMemSort, dimBlockSharedMemSort>>>(((float2*)array) + ((padding / Constants::CUDA_ELEMENTS_SORT) * (Constants::CUDA_ELEMENTS_SORT / 2)), 
			((int2*)pos_indices) + ((padding / Constants::CUDA_ELEMENTS_SORT) * (Constants::CUDA_ELEMENTS_SORT / 2)));

    for (int k = (Constants::CUDA_ELEMENTS_SORT << 1), r = (int)(Constants::CUDA_SORT_NUM << 1); k <= numElements; k *= 2, r <<= 1) 
		{

        for (int j = k / 2, s = r >> 1; j > (Constants::CUDA_ELEMENTS_MERGE / 2); j >>= 1, s >>= 1) 
				{
            mergeArray<<<dimGridMergeArray, dimBlockMergeArray>>>(array, pos_indices, j, k, s);
        }

        sharedMemMerge<<<dimGridSharedMemMerge, dimBlockSharedMemMerge>>>(array, pos_indices, k);
    }		
		cudaThreadSynchronize();
		checkCUDAError("sortOnDevice");

		int num_threads = Constants::CUDA_BLOCK_SIZE_MAX;
		while (numElements < num_threads && num_threads > 1)
		{
			num_threads /= 2;
		};

		if (num_threads == 1) //this case should never happen
		{
			return (0);
		}; 

		dim3 dimGrid (numElements/num_threads);
		dim3 dimBlock (num_threads);

		void* dev_cut_off_index;
		cudaMalloc (&dev_cut_off_index, sizeof(int));
		cudaMemset (dev_cut_off_index, -1, sizeof(int));

		findCutOffIndex<<<dimGrid, dimBlock>>> (array, (int*)dev_cut_off_index);
		cudaThreadSynchronize();
		checkCUDAError("findCutoffIndex");
		int cut_off_index=-1;
		cudaMemcpy (&cut_off_index, dev_cut_off_index, sizeof(int), cudaMemcpyDeviceToHost);

		return (cut_off_index);
	}


	extern __shared__ float external_shared [];
	__global__ void scoreIndividuals (float* scores, const int overall_size, const int c, const int offset,  const int write_offset, const float ampl_cutoff)
	{		
		int v = threadIdx.x;
		int ref_index = tex1Dfetch (sorted_positions_indices_tex, blockIdx.x+offset);
		
		//printf ("my_index: %i\n", my_index);
		//printf ("ref_index: %i\n", ref_index); 	
	
		__shared__ int peak_cutoff, optimal_block_dim;
		__shared__ float seed_mz;
		float* c_scores = (float*) &external_shared[0]; 
		
		if (v==0)
		{	
			seed_mz = tex1Dfetch(pos_tex, ref_index);
			peak_cutoff = getNumPeakCutOff(seed_mz, c+1);	
			optimal_block_dim = 4*(peak_cutoff-1) -1;
			//optimal_block_dim = 2*(peak_cutoff-1);
		};

		__syncthreads();
		if (v < optimal_block_dim)
		{
			float my_mz, l_pos, l_intens; int l_index;
			my_mz = seed_mz-((peak_cutoff-1)*Constants::IW_NEUTRON_MASS-(v+1)*Constants::IW_HALF_NEUTRON_MASS)/((float)(c+1));
			//my_mz = seed_mz+v*Constants::IW_HALF_NEUTRON_MASS/((float)(c+1));

			l_index = ref_index;
			if (my_mz > seed_mz)
			{ 
				while (l_index < overall_size && tex1Dfetch(pos_tex, l_index++) < my_mz) 
				{ 
				};
				if (l_index<overall_size)
				{
					l_index -= 2;
				};
			}
			else
			{
				while (l_index >= 0 && tex1Dfetch(pos_tex,l_index--) > my_mz) 
				{							
				};
				if (l_index >=0)
				{
					++l_index;
				};					
			};

			if (l_index >=0  && l_index+1 < overall_size)
			{	
				l_pos = tex1Dfetch(pos_tex, l_index);
				l_intens = tex1Dfetch(trans_intensities_tex, l_index);
				c_scores[v] = l_intens + ( tex1Dfetch(trans_intensities_tex, l_index+1)-l_intens ) / (tex1Dfetch(pos_tex, l_index+1) - l_pos) * (my_mz - l_pos); 				
				//printf ("Scoring: %f\t\t%f\t%f\t%f\t%f\t%f\n",  seed_mz, my_mz, tex1Dfetch(pos_tex, l_index+1), l_pos, tex1Dfetch(trans_intensities_tex, l_index+1), l_intens);
			}
			else
			{
				c_scores[v]=INT_MIN;
			};
		};

		__syncthreads();

		//It has been test that an advanced reduction scheme does not offer
		//any performance advantages in our case; so we use the greedy way here ...	
		
		
		 //It has been test that an advanced reduction scheme does not offer
                //any performance advantages in our case; so we use the greedy way here ...
             /*   if (v==0)
                {
                        float final_score = 0, mid_val=0, l_score=0;
                        int minus = -1, i;
                        for (i=0; i<(int)ceil(optimal_block_dim/2.); ++i)
                        {
                                if (c_scores[i] != INT_MIN)
                                {
                                        final_score += minus*c_scores[i];
                                };
                                minus *=-1;
                        };

                        l_score = final_score;
                        mid_val = c_scores[i];

                        for (; i<optimal_block_dim && c_scores[i] != INT_MIN; ++i)
                        {
                                final_score += minus*c_scores[i];
                                minus *=-1;
                        };

                        if (!(l_score <=0 || final_score-l_score-mid_val <= 0 || final_score-mid_val <= ampl_cutoff))
                        {
                                scores[blockIdx.x+write_offset] = final_score;
                        };
                        //printf ("blockid: %i\t%i\n", blockIdx.x, write_offset);
                        //printf("final_score: %f\t\t%f\n", seed_mz, final_score);
                };*/

		
		
		__shared__ float mid_val, l_score, r_score;
		if (v==0)
		{
			l_score=0; mid_val=0;
			int minus = -1; int i;
			for (i=0; i<(int)floor(optimal_block_dim/2.); ++i)
			{
				if (c_scores[i] != INT_MIN)
				{
					//if (trunc(seed_mz*100) == 80442) printf("l: %f\n", c_scores[i]); 
					
					l_score += minus*c_scores[i];
				};
				minus *=-1;
			};
			mid_val = c_scores[i];
			//if (trunc(seed_mz*100) == 80442) printf("m: %f\n", c_scores[i]);
		};
		if(v==1)
		{
			r_score=0;
			int minus = -1;
			for (int i=(int)floor(optimal_block_dim/2.)+1; i<optimal_block_dim && c_scores[i] != INT_MIN; ++i)
			{
				//if (trunc(seed_mz*100) == 80442) printf("r: %f\n", c_scores[i]);

				r_score += minus*c_scores[i];
				minus *=-1;
			};
		};

		__syncthreads();

		if(v==0)
		{
			if (!(l_score <=0 || r_score <= 0 || l_score + r_score <= ampl_cutoff))
			{
				scores[blockIdx.x+write_offset] = l_score + r_score + mid_val;	
			};
			//printf ("blockid: %i\t%i\n", blockIdx.x, write_offset);
			/*if (trunc(seed_mz*100) == 80442) 
			{
				printf("final_score: %f\t\t%f\t%f\t%f\n", seed_mz,  l_score, mid_val, r_score);
			};*/
		};
	};

	
	void scoreOnDevice (int* sorted_positions_indices, float* trans_intensities, float* pos, float* scores, 
		const int c, const int num_of_scores, const int overall_size, const unsigned int max_peak_cutoff, const float ampl_cutoff)
	{
		int theo_block_dim = 4*(max_peak_cutoff-1)-1; //the number of scoring points per candidates, due to numerical reasons we increse max_peak_cutoff by one
		dim3 blockDim (theo_block_dim); //the number of scoring points per candidates

		cudaBindTexture(0, trans_intensities_tex, trans_intensities, overall_size*sizeof(float));
		cudaBindTexture(0, pos_tex, pos, overall_size*sizeof(float));
		cudaBindTexture(0, sorted_positions_indices_tex, sorted_positions_indices, overall_size*sizeof(int));	
		size_t offset = overall_size - num_of_scores;

		//printf ("num_of_scores: %i\n", num_of_scores);
		//printf ("overall_size: %i\n", overall_size);
		dim3 gridDim (Constants::CUDA_BLOCKS_PER_GRID_MAX);
		int counts=0, c_size = num_of_scores;

		while ((c_size -= Constants::CUDA_BLOCKS_PER_GRID_MAX) > 0)
		{		
			scoreIndividuals<<<gridDim, blockDim, blockDim.x*sizeof(float)>>> (scores, overall_size, c, 
				counts*Constants::CUDA_BLOCKS_PER_GRID_MAX+offset, counts*Constants::CUDA_BLOCKS_PER_GRID_MAX, ampl_cutoff);	
			++counts;
		};

		if ((c_size += Constants::CUDA_BLOCKS_PER_GRID_MAX) > 0)
		{
			gridDim = dim3 (c_size);
			scoreIndividuals<<<gridDim, blockDim, blockDim.x*sizeof(float)+2*sizeof(int)+sizeof(float)>>> (scores, overall_size, c, 
				counts*Constants::CUDA_BLOCKS_PER_GRID_MAX+offset, counts*Constants::CUDA_BLOCKS_PER_GRID_MAX, ampl_cutoff);		
		};
		
		cudaThreadSynchronize();
		checkCUDAError("scoreOnDevice");
		
		cudaUnbindTexture (trans_intensities_tex);
		cudaUnbindTexture (pos_tex);
		cudaUnbindTexture (sorted_positions_indices_tex);
	}

}	
