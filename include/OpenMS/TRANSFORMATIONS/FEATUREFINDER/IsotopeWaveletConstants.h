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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELETCONSTANTS_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELETCONSTANTS_H

namespace OpenMS 
{
	//We cannot define these constants as extern variable because of the cuda interface
	namespace Constants 
	{
		#undef OPENMS_DEBUG_ISOTOPE_WAVELET

		//the following case should never happen (due to configure)
		//but just to be sure ...
		#ifdef OPENMS_HAS_TBB_H
			#ifndef OPENMS_HAS_CUDA //the use threading building blocks requires a working GPU device
				#undef OPENMS_HAS_TBB_H
			#endif
		#endif 

		const unsigned int DEFAULT_NUM_OF_INTERPOLATION_POINTS = 3;

		const double MASS_EPSILON = 1e-4f;

		const double MARR_WAVELET_CUTOFF = 4.f;

		const double PEPTIDE_MASS_RULE_FACTOR = 0.000507f;
		const	double PEPTIDE_MASS_RULE_BOUND = 	1./PEPTIDE_MASS_RULE_FACTOR;
		const double PEPTIDE_MASS_RULE_THEO_PPM_BOUND = 200;

		//exact
		/*const double IW_NEUTRON_MASS = 1.00866491578f; 
		const double IW_HALF_NEUTRON_MASS = 0.5043325f;
		const double IW_QUARTER_NEUTRON_MASS = 0.252166228f;
		const double WAVELET_PERIODICITY = 6.229209734f;*/

		const double IW_NEUTRON_MASS = 1.00235f; 
		const double IW_HALF_NEUTRON_MASS = 0.501175f;
		const double IW_QUARTER_NEUTRON_MASS = 0.2505875f;
		const double WAVELET_PERIODICITY = 6.268454439f;
		const double ONEOLOG2E = 0.6931471806f;

		const double IW_PROTON_MASS = 1.00727646688f;

		//Linear Fit (standard)
		const double LAMBDA_L_0 = -0.472998839574110749e-1f;
		const double LAMBDA_L_1 = 0.743579753540513913e-3f;

		//Quadratic Fit (maybe a overkill, since the dependency looks quite linear, at least in moderate mass ranges)
		const double LAMBDA_Q_0 = -0.137152573151174711f;
		const double LAMBDA_Q_1 = 0.851289601785403817e-3f;
		const double LAMBDA_Q_2 = -0.2834469691e-7f;

		const double CUTOFF_FIT99_POLY_0 = 0.813303208695311253f;
		const double CUTOFF_FIT99_POLY_1 = 0.00198668967903372230f;
		const double CUTOFF_FIT99_POLY_2 = -1.17957442070674527e-7f;

		const double BORDER_MZ_FIT99 = 4200.109412;

		const double CUTOFF_FIT99_LOG_0 = -31.9263753269393682f;
		const double CUTOFF_FIT99_LOG_1 = 4.67502601032944209f;

		const int SHIFT23 = (1<<23);
		const double SHIFT23_00 = (1.0/(1<<23));
		const double LOG_CONST = 0.346607f;
		const double POW_CONST = 0.33971f;
		
		const int CUDA_INIT_SUCCESS = 1;
		const int CUDA_INIT_FAIL = -1; 

		const int CUDA_BLOCKS_PER_GRID_MAX = 65535;
		const int CUDA_BLOCK_SIZE_MAX = 256;//limited due to the shared memory
		const int CUDA_EXTENDED_BLOCK_SIZE_MAX = 2039;
		const int CUDA_TEXTURE_THREAD_LIMIT = 384;//limited due to the number of used registers

		const int CUDA_ELEMENTS_SORT = 512;
		const int CUDA_THREADS_SORT =  (CUDA_ELEMENTS_SORT >> 2);
		const int CUDA_SORT_NUM = 0xFFFFFE00;
		const int CUDA_ELEMENTS_MERGE = 1024;
		const int CUDA_MERGE_NUM = 0xFFFFFC00;
		const int CUDA_THREADS_MERGE = (CUDA_ELEMENTS_MERGE >> 2);
		const int CUDA_ELEMENTS_GL = 2;
		const int CUDA_THREADS_GL = 256;
		const int CUDA_MIN_SORT_SIZE = CUDA_ELEMENTS_SORT;

		const int TBB_NUM_OF_GPU_DEVICES = 2;

		#define CUDA_CHECK_ERROR
	}
}

#endif
