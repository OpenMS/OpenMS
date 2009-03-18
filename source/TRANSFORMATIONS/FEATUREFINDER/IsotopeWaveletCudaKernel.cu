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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletCudaKernel.h>


#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <iomanip>

#include <assert.h>
#include <cuda.h>


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


#define MAX_BLOCKS_PER_GRID 65535


texture<float,1> trans_intensities_tex, pos_tex;
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
		float fac (-(LAMBDA_Q_0 + LAMBDA_Q_1*mz + LAMBDA_Q_2*mz*mz));
		fac += (tz1-1)*__log2f(-fac)*ONEOLOG2E - lgammaf(tz1);
			
		return (__sinf((tz1-1)*WAVELET_PERIODICITY) * __expf(fac));
	}



	__global__ void ConvolutionIsotopeWaveletKernel(float* signal_pos, float* signal_int, const int from_max_to_left, const int from_max_to_right, float* result, 
		const unsigned int charge, const int block_size, const float peak_cutoff_intercept, const float peak_cutoff_slope)
	{
		// the device-shared memory storing one data block
		// This is currently hard-coded to 256 points, since we require two 4B floats for each
		// data point, leading to 2kB per block.
		__shared__ float signal_pos_block[BLOCK_SIZE_MAX];
		__shared__ float signal_int_block[BLOCK_SIZE_MAX];

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
			
		int problem_size_num_of_threads = block_size-(from_max_to_left+from_max_to_right/*+1-1*/);
		// the position in the original signal array that corresponds the data point computed by this thread
		//
		//                 left padding,                                                      position of thread
		//                 ignored in output    the points computed by the previous blocks    in block
		int my_data_pos  = from_max_to_left    +  blockIdx.x*problem_size_num_of_threads + threadIdx.x;
		int my_local_pos = threadIdx.x + from_max_to_left;

		//every thread loads its own (maximum) position into the shared memory
		signal_pos_block[my_local_pos] = signal_pos[my_data_pos];
		//signal_pos_block[my_local_pos] = tex1Dfetch(pos_tex, my_data_pos);
		signal_int_block[my_local_pos] = signal_int[my_data_pos];
		//printf ("%i\t\t%i\t\t %f\t\t%f\n",  blockIdx.x, threadIdx.x, signal_int[my_data_pos], tex1Dfetch(cuda_device_intens_texture_reference_, my_data_pos));
		//printf ("norm-loading: %f\t%i\n", signal_pos[my_data_pos], my_local_pos);

		//every thread with an ID smaller than the number of from_max_to_left loads the additional boundary points
		//at the left end
		if (threadIdx.x < from_max_to_left)
		{
			signal_pos_block[threadIdx.x] = signal_pos[my_data_pos-from_max_to_left];
			//signal_pos_block[threadIdx.x] = tex1Dfetch(pos_tex, my_data_pos-from_max_to_left);
			signal_int_block[threadIdx.x] = signal_int[my_data_pos-from_max_to_left];
			//printf ("pre-loading: %f\t%i\n",  signal_pos[my_data_pos-from_max_to_left], threadIdx.x);
			

			signal_pos_block[block_size-threadIdx.x-1] = signal_pos[my_data_pos+block_size-2*threadIdx.x-1-from_max_to_left];
			//signal_pos_block[block_size-threadIdx.x-1] = tex1Dfetch(pos_tex, my_data_pos+block_size-2*threadIdx.x-1-from_max_to_left);
			signal_int_block[block_size-threadIdx.x-1] = signal_int[my_data_pos+block_size-2*threadIdx.x-1-from_max_to_left];
			//printf ("extra-loading: %f\t%i\n", signal_pos[my_data_pos+block_size-2*threadIdx.x-1-from_max_to_left], block_size-threadIdx.x-1);

		}
		else
		{
			//int num_threads_with_only_one_load = problem_size_num_of_threads-from_max_to_left;
			//additional loads to be done at the right end: from_max_to_right
			int additional_right_end_loads=0;
			//printf ("crit: %i\t%i\t%i\n", (int)(from_max_to_right/num_threads_with_only_one_load), from_max_to_right, num_threads_with_only_one_load); 
			//while ((int)(from_max_to_right/num_threads_with_only_one_load) - additional_right_end_loads > 0)
			//while (my_local_pos + (additional_right_end_loads+1)*num_threads_with_only_one_load < block_size-2*from_max_to_left)
			while (my_local_pos + (additional_right_end_loads+1)*(problem_size_num_of_threads-from_max_to_left) < block_size-from_max_to_left)
			{
					++additional_right_end_loads;
					signal_pos_block[my_local_pos+additional_right_end_loads*(problem_size_num_of_threads-from_max_to_left)] = signal_pos[my_data_pos+additional_right_end_loads*(problem_size_num_of_threads-from_max_to_left)];
					//signal_pos_block[my_local_pos+additional_right_end_loads*(problem_size_num_of_threads-from_max_to_left)] = tex1Dfetch(pos_tex, my_data_pos+additional_right_end_loads*(problem_size_num_of_threads-from_max_to_left));

					signal_int_block[my_local_pos+additional_right_end_loads*(problem_size_num_of_threads-from_max_to_left)] = signal_int[my_data_pos+additional_right_end_loads*(problem_size_num_of_threads-from_max_to_left)];
					//printf ("post-loading: %f\t%i\n", signal_pos[my_data_pos+additional_right_end_loads*num_threads_with_only_one_load], my_local_pos+additional_right_end_loads*num_threads_with_only_one_load);
			};
			/*additional_right_end_loads=0;		
			//leave this loops separated from each other; this reduces the number of registers in use, which can have significant impact on the occupancy of the program
			while (my_local_pos + (additional_right_end_loads+1)*(problem_size_num_of_threads-from_max_to_left) < block_size-from_max_to_left)
			{
					++additional_right_end_loads;
					//signal_pos_block[my_local_pos+additional_right_end_loads*(problem_size_num_of_threads-from_max_to_left)] = signal_pos[my_data_pos+additional_right_end_loads*(problem_size_num_of_threads-from_max_to_left)];
					signal_int_block[my_local_pos+additional_right_end_loads*(problem_size_num_of_threads-from_max_to_left)] = signal_int[my_data_pos+additional_right_end_loads*(problem_size_num_of_threads-from_max_to_left)];
					//printf ("post-loading: %f\t%i\n", signal_pos[my_data_pos+additional_right_end_loads*num_threads_with_only_one_load], my_local_pos+additional_right_end_loads*num_threads_with_only_one_load);
			};*/

		};
		//wait until the shared data is loaded completely
		__syncthreads(); 

		//exit(-1);

		//my_local_pos = threadIdx.x*(from_max_to_left+from_max_to_right+1) % problem_size_num_of_threads + from_max_to_left;
		//my_data_pos = my_local_pos + blockIdx.x*problem_size_num_of_threads;

		//printf ("%i\t%i\n", my_local_pos, my_data_pos); 

		float value = 0, boundary=(ceil(peak_cutoff_intercept+peak_cutoff_slope*charge*signal_pos_block[my_local_pos])*NEUTRON_MASS)/charge, c_diff;
		for (int current_conv_pos = my_local_pos-from_max_to_left; 
								current_conv_pos < my_local_pos+from_max_to_right; 
							++current_conv_pos)
		{
			//printf ("%i\t", current_conv_pos);
			//current_value = signal_int_block[current_conv_pos];
			//current_pos   = signal_pos_block[current_conv_pos];
			//if ((signal_pos_block[current_conv_pos]-center_pos) < 0) continue;
			//c_val = isotope_wavelet((current_pos-center_pos)*charge+1., current_pos);
			c_diff = signal_pos_block[current_conv_pos]-signal_pos_block[my_local_pos]+QUARTER_NEUTRON_MASS/charge;

			value += c_diff > 0 && c_diff <= boundary ? isotope_wavelet(c_diff*charge+1., signal_pos_block[current_conv_pos]*charge)*signal_int_block[current_conv_pos] : 0;
		};

		result[my_data_pos] = value;
	}


	void getExternalCudaTransforms (dim3 dimGrid, dim3 dimBlock, float* positions_dev, float* intensities_dev, int from_max_to_left, int from_max_to_right, float* result_dev, 
		const int charge, const int block_size, const float peak_cutoff_intercept, const float peak_cutoff_slope) 
	{
		ConvolutionIsotopeWaveletKernel<<<dimGrid,dimBlock>>> (positions_dev, intensities_dev, from_max_to_left, from_max_to_right, result_dev, charge, block_size, peak_cutoff_intercept, peak_cutoff_slope);
		cudaThreadSynchronize();
		checkCUDAError("ConvolutionIsotopeWaveletKernel");
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

			__shared__ float shmem[ELEMENTS_MERGE];
			__shared__ int posshmem[ELEMENTS_MERGE];

			int tmp = blockIdx.x * ELEMENTS_MERGE + threadIdx.x;

			float data = array[tmp];
			float data2 = array[tmp + (ELEMENTS_MERGE / 2)];

			float data3 = array[tmp + THREADS_MERGE];
			float data4 = array[tmp + THREADS_MERGE + (ELEMENTS_MERGE / 2)];
			
			int posdata = pos[tmp];
			int posdata2 = pos[tmp + (ELEMENTS_MERGE / 2)];

			int posdata3 = pos[tmp + THREADS_MERGE];
			int posdata4 = pos[tmp + THREADS_MERGE + (ELEMENTS_MERGE / 2)];

			int dir = k & (blockIdx.x * (ELEMENTS_MERGE));


			if (dir == 0) {
					if (data > data2) {  // ascending
							shmem[threadIdx.x] = data2;
							shmem[threadIdx.x + (ELEMENTS_MERGE / 2)] = data;
							posshmem[threadIdx.x] = posdata2;
							posshmem[threadIdx.x + (ELEMENTS_MERGE / 2)] = posdata;
					} else {
							shmem[threadIdx.x] = data;
							shmem[threadIdx.x + (ELEMENTS_MERGE / 2)] = data2;
							posshmem[threadIdx.x] = posdata;
							posshmem[threadIdx.x + (ELEMENTS_MERGE / 2)] = posdata2;
					}

					if (data3 > data4) {  // ascending
							shmem[threadIdx.x + THREADS_MERGE] = data4;
							shmem[threadIdx.x + (ELEMENTS_MERGE / 2) + THREADS_MERGE] = data3;
							posshmem[threadIdx.x + THREADS_MERGE] = posdata4;
							posshmem[threadIdx.x + (ELEMENTS_MERGE / 2) + THREADS_MERGE] = posdata3;
					} else {
							shmem[threadIdx.x + THREADS_MERGE] = data3;
							shmem[threadIdx.x + (ELEMENTS_MERGE / 2) + THREADS_MERGE] = data4;
							posshmem[threadIdx.x + THREADS_MERGE] = posdata3;
							posshmem[threadIdx.x + (ELEMENTS_MERGE / 2) + THREADS_MERGE] = posdata4;
					}
			} else {
					if (data < data2) {  // descending
							shmem[threadIdx.x] = data2;
							shmem[threadIdx.x + (ELEMENTS_MERGE / 2)] = data;
							posshmem[threadIdx.x] = posdata2;
							posshmem[threadIdx.x + (ELEMENTS_MERGE / 2)] = posdata;
					} else {
							shmem[threadIdx.x] = data;
							shmem[threadIdx.x + (ELEMENTS_MERGE / 2)] = data2;							
							posshmem[threadIdx.x] = posdata;
							posshmem[threadIdx.x + (ELEMENTS_MERGE / 2)] = posdata2;
					}

					if (data3 < data4) {  // descending
							shmem[threadIdx.x + THREADS_MERGE] = data4;
							shmem[threadIdx.x + (ELEMENTS_MERGE / 2) + THREADS_MERGE] = data3;
							posshmem[threadIdx.x + THREADS_MERGE] = posdata4;
							posshmem[threadIdx.x + (ELEMENTS_MERGE / 2) + THREADS_MERGE] = posdata3;
					} else {
							shmem[threadIdx.x + THREADS_MERGE] = data3;
							shmem[threadIdx.x + (ELEMENTS_MERGE / 2) + THREADS_MERGE] = data4;							
							posshmem[threadIdx.x + THREADS_MERGE] = posdata3;
							posshmem[threadIdx.x + (ELEMENTS_MERGE / 2) + THREADS_MERGE] = posdata4;
					}
			}



			int j = 256, s = MERGE_NUM >> 2; 


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

			x = (threadIdx.x + THREADS_MERGE) + ((threadIdx.x + THREADS_MERGE) & s);
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

			x = (threadIdx.x + THREADS_MERGE) + ((threadIdx.x + THREADS_MERGE) & s);
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

			x = (threadIdx.x + THREADS_MERGE) + ((threadIdx.x + THREADS_MERGE) & s);
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

			x = (threadIdx.x + THREADS_MERGE) + ((threadIdx.x + THREADS_MERGE) & s);
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

			x = (threadIdx.x + THREADS_MERGE) + ((threadIdx.x + THREADS_MERGE) & s);
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

			x = (threadIdx.x + THREADS_MERGE) + ((threadIdx.x + THREADS_MERGE) & s);
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

			x = (threadIdx.x + THREADS_MERGE) + ((threadIdx.x + THREADS_MERGE) & s);
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

			x = (threadIdx.x + THREADS_MERGE) + ((threadIdx.x + THREADS_MERGE) & s);
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

			x = (threadIdx.x + THREADS_MERGE) + ((threadIdx.x + THREADS_MERGE) & s);
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

			int i = blockIdx.x * ELEMENTS_MERGE + threadIdx.x;
			array[i] = shmem[threadIdx.x];
			array[i + (ELEMENTS_MERGE / 4)] = shmem[(ELEMENTS_MERGE / 4) + threadIdx.x];
			array[i + (ELEMENTS_MERGE / 2)] = shmem[(ELEMENTS_MERGE / 2) + threadIdx.x];
			array[i + (3 * ELEMENTS_MERGE / 4)] = shmem[(3 * ELEMENTS_MERGE / 4) + threadIdx.x];

			pos[i] = posshmem[threadIdx.x];
			pos[i + (ELEMENTS_MERGE / 4)] = posshmem[(ELEMENTS_MERGE / 4) + threadIdx.x];
			pos[i + (ELEMENTS_MERGE / 2)] = posshmem[(ELEMENTS_MERGE / 2) + threadIdx.x];
			pos[i + (3 * ELEMENTS_MERGE / 4)] = posshmem[(3 * ELEMENTS_MERGE / 4) + threadIdx.x];
	}




	__global__ void mergeArray(float *array, int* pos, int j, int k, int s) {
			int tmp = (blockIdx.x * THREADS_GL);
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
			__shared__ float shmem[ELEMENTS_SORT];
			__shared__ int posshmem[ELEMENTS_SORT];

			float2 data = array[blockIdx.x * (ELEMENTS_SORT / 2) + threadIdx.x];
			int2 posdata = pos[blockIdx.x * (ELEMENTS_SORT / 2) + threadIdx.x];


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


			for (int k = 4, r = 0xFFFFFFFC; k <= (ELEMENTS_SORT / 2); k *= 2, r <<= 1) {

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


			data = array[blockIdx.x * (ELEMENTS_SORT / 2) + THREADS_SORT + threadIdx.x];
			posdata = pos[blockIdx.x * (ELEMENTS_SORT / 2) + THREADS_SORT + threadIdx.x];
			float* shmem2 = &shmem[ELEMENTS_SORT / 2];
			int* posshmem2 = &posshmem[ELEMENTS_SORT / 2];
			
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


			for (int k = 4, r = 0xFFFFFFFC; k <= (ELEMENTS_SORT / 2); k *= 2, r <<= 1) {

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

					for (int j = ELEMENTS_SORT / 2, s = SORT_NUM >> 1; j > 0; j >>= 1, s >>= 1) {

							int x = threadIdx.x + (threadIdx.x & s);
							int y = x + j;
							__syncthreads();

							if (shmem[x] > shmem[y]) {  // ascending
									swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
							}

							x = (threadIdx.x + THREADS_SORT) + ((threadIdx.x + THREADS_SORT) & s);
							y = x + j;
			
							if (shmem[x] > shmem[y]) {  // ascending
									swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
							}
					}
			
			} else {
			
					for (int j = ELEMENTS_SORT / 2, s = SORT_NUM >> 1; j > 0; j >>= 1, s >>= 1) {

							int x = threadIdx.x + (threadIdx.x & s);
							int y = x + j;
							__syncthreads();

							if (shmem[x] < shmem[y]) {  // descending
									swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
							}

							x = (threadIdx.x + THREADS_SORT) + ((threadIdx.x + THREADS_SORT) & s);
							y = x + j;
			
							if (shmem[x] < shmem[y]) {  // descending
									swap(shmem[x], shmem[y], posshmem[x], posshmem[y]);
							}
					}
			
			}


			__syncthreads();

			int i = blockIdx.x * ELEMENTS_SORT + threadIdx.x;
			((float*)array)[i] = shmem[threadIdx.x];
			((float*)array)[i + (ELEMENTS_SORT / 4)] = shmem[(ELEMENTS_SORT / 4) + threadIdx.x];
			((float*)array)[i + (ELEMENTS_SORT / 2)] = shmem[(ELEMENTS_SORT / 2) + threadIdx.x];
			((float*)array)[i + (3 * ELEMENTS_SORT / 4)] = shmem[(3 * ELEMENTS_SORT / 4) + threadIdx.x];

			((int*)pos)[i] = posshmem[threadIdx.x];
			((int*)pos)[i + (ELEMENTS_SORT / 4)] = posshmem[(ELEMENTS_SORT / 4) + threadIdx.x];
			((int*)pos)[i + (ELEMENTS_SORT / 2)] = posshmem[(ELEMENTS_SORT / 2) + threadIdx.x];
			((int*)pos)[i + (3 * ELEMENTS_SORT / 4)] = posshmem[(3 * ELEMENTS_SORT / 4) + threadIdx.x];
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
	  dim3 dimGridSharedMemSort((numElements / ELEMENTS_SORT) - (padding / ELEMENTS_SORT), 1, 1);
    dim3 dimBlockSharedMemSort(THREADS_SORT, 1, 1);

    dim3 dimGridMergeArray(numElements / (THREADS_GL * ELEMENTS_GL), 1, 1);
    dim3 dimBlockMergeArray(THREADS_GL, 1, 1);

    dim3 dimGridSharedMemMerge(numElements / ELEMENTS_MERGE, 1, 1);
    dim3 dimBlockSharedMemMerge(THREADS_MERGE, 1, 1);

		sharedMemSort<<<dimGridSharedMemSort, dimBlockSharedMemSort>>>(((float2*)array) + ((padding / ELEMENTS_SORT) * (ELEMENTS_SORT / 2)), ((int2*)pos_indices) + ((padding / ELEMENTS_SORT) * (ELEMENTS_SORT / 2)));

    for (int k = (ELEMENTS_SORT << 1), r = (int)(SORT_NUM << 1); k <= numElements; k *= 2, r <<= 1) 
		{

        for (int j = k / 2, s = r >> 1; j > (ELEMENTS_MERGE / 2); j >>= 1, s >>= 1) 
				{
            mergeArray<<<dimGridMergeArray, dimBlockMergeArray>>>(array, pos_indices, j, k, s);
        }

        sharedMemMerge<<<dimGridSharedMemMerge, dimBlockSharedMemMerge>>>(array, pos_indices, k);
    }		
		cudaThreadSynchronize();
		checkCUDAError("sortOnDevice");

		int num_threads = BLOCK_SIZE_MAX;
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

		//std::cout << "numElements:" << numElements << "\t" << dimGrid.x << "\t" << dimBlock.x << std::endl;

		findCutOffIndex<<<dimGrid, dimBlock>>> (array, (int*)dev_cut_off_index);
		cudaThreadSynchronize();
		checkCUDAError("findCutoffIndex");
		int cut_off_index=-1;
		cudaMemcpy (&cut_off_index, dev_cut_off_index, sizeof(int), cudaMemcpyDeviceToHost);

		//std::cout << "Found cutoff index as: " << cut_off_index << std::endl;
		return (cut_off_index);
	}


	

	extern __shared__ float c_scores [];
	__global__ void scoreIndividuals (int* sorted_positions_indices, float* pos, float* trans_intensities, float* scores, const int overall_size, 
		const int c, const int offset,  const int write_offset, const float peak_cutoff_intercept, const float peak_cutoff_slope)
	{		
		int v = threadIdx.x;
		//int ref_index = sorted_positions_indices[blockIdx.x+offset];
		int ref_index = tex1Dfetch (sorted_positions_indices_tex, blockIdx.x+offset);

		//printf ("my_index: %i\n", my_index);
		//printf ("ref_index: %i\n", ref_index); 	
	
		__shared__ int peak_cutoff, optimal_block_dim;
		__shared__ float seed_mz;

		if (v==0)
		{	
			seed_mz = tex1Dfetch(pos_tex, ref_index);
			peak_cutoff = (int) ceil(peak_cutoff_intercept+peak_cutoff_slope*seed_mz*(c+1));	
			optimal_block_dim = 4*(peak_cutoff-1) -1;
		};

		__syncthreads();
		if (v < optimal_block_dim)
		{
			float my_mz, l_pos, l_intens; int l_index;
			my_mz = seed_mz-((peak_cutoff-1)*NEUTRON_MASS-(v+1)*HALF_NEUTRON_MASS)/((float)(c+1));

			l_index = ref_index;
			if (my_mz > seed_mz)
			{ 
				while (l_index < overall_size && tex1Dfetch(pos_tex, l_index++) < my_mz) 
				{ 
				};
				l_index -= 2;
			}
			else
			{
				while (l_index >= 0 && tex1Dfetch(pos_tex,l_index--) > my_mz) 
				{							
				};
				++l_index;					
			};

			if (l_index >=0  && l_index+1 < overall_size)
			{	
				//l_pos = pos[l_index];
				l_pos = tex1Dfetch(pos_tex, l_index);
				//l_intens = trans_intensities[l_index];
				l_intens = tex1Dfetch(trans_intensities_tex, l_index);
				//c_scores[v] = l_intens + ( trans_intensities[l_index+1]-l_intens ) / (pos[l_index+1] - l_pos) * (my_mz - l_pos); 
				
				c_scores[v] = l_intens + ( tex1Dfetch(trans_intensities_tex, l_index+1)-l_intens ) / (tex1Dfetch(pos_tex, l_index+1) - l_pos) * (my_mz - l_pos); 				
				//printf ("Scoring: %f\t\t%f\t%f\t%f\t%f\t%f\n",  seed_mz, my_mz, tex1Dfetch(pos_tex, l_index+1), l_pos, tex1Dfetch(trans_intensities_tex, l_index+1), l_intens);
			}
			else
			{
				c_scores[v]=0;
			};
		};

		__syncthreads();	

		if (v==0)
		{
			float final_score = 0;
			int minus = -1;
			for (int i=0; i<optimal_block_dim; ++i)
			{
				final_score += minus*c_scores[i];
				minus *=-1;
			};
			scores[blockIdx.x+write_offset] = final_score;
		};			
	};


	void scoreOnDevice (int* sorted_positions_indices, float* trans_intensities, float* pos, float* scores, 
		const int c, const int num_of_scores, const int overall_size, const float peak_cutoff_intercept, const float peak_cutoff_slope, const unsigned int max_peak_cutoff)
	{
		//just to be sure, we will have 4 additional threads for numerical reasons that might trigger some additional scoring points
		dim3 blockDim (4*(max_peak_cutoff) -1); //the number of scoring points per candidates

		cudaBindTexture(0, trans_intensities_tex, trans_intensities, overall_size*sizeof(float));
		cudaBindTexture(0, pos_tex, pos, overall_size*sizeof(float));
		cudaBindTexture(0, sorted_positions_indices_tex, sorted_positions_indices, overall_size*sizeof(int));	
		size_t schrott = overall_size - num_of_scores;

		dim3 gridDim;
		int counts=0, c_size = num_of_scores;
		while ((c_size -= MAX_BLOCKS_PER_GRID) > 0)
		{		
			gridDim = dim3 (MAX_BLOCKS_PER_GRID);

			scoreIndividuals<<<gridDim, blockDim, blockDim.x*sizeof(float)>>> (sorted_positions_indices, pos, trans_intensities, scores, overall_size, c, counts*MAX_BLOCKS_PER_GRID+schrott, counts*MAX_BLOCKS_PER_GRID, 
				peak_cutoff_intercept, peak_cutoff_slope);	
			++counts;
		};


		if ((c_size += MAX_BLOCKS_PER_GRID) > 0)
		{
			gridDim = dim3 (c_size);

			scoreIndividuals<<<gridDim, blockDim, blockDim.x*sizeof(float)>>> (sorted_positions_indices, pos, trans_intensities, scores, overall_size, c, counts*MAX_BLOCKS_PER_GRID+schrott, counts*MAX_BLOCKS_PER_GRID, 
				peak_cutoff_intercept, peak_cutoff_slope);		
		};
		
//		CUDA_SAFE_CALL(cudaThreadSynchronize());
		cudaThreadSynchronize();
		checkCUDAError("scoreOnDevice");
		
		cudaUnbindTexture (trans_intensities_tex);
		cudaUnbindTexture (pos_tex);
		cudaUnbindTexture (sorted_positions_indices_tex);
	}

}	
