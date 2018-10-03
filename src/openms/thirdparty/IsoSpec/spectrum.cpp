/*
 *   Copyright (C) 2015-2016 Mateusz Łącki and Michał Startek.
 *
 *   This file is part of IsoSpec.
 *
 *   IsoSpec is free software: you can redistribute it and/or modify
 *   it under the terms of the Simplified ("2-clause") BSD licence.
 *
 *   IsoSpec is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
 *
 *   You should have received a copy of the Simplified BSD Licence
 *   along with IsoSpec.  If not, see <https://opensource.org/licenses/BSD-2-Clause>.
 */

#include <math.h>
#include <string.h>
#include <iostream>
#include "spectrum.h"
#include "isoMath.h"
#include "isoSpec++.h"

Kernel::Kernel(double _delta, double* _k, double _bucketsize, double _buckets) : 
delta(_delta),
k(_k),
bucketsize(_bucketsize),
buckets(_buckets)
{}

Kernel* Kernel::SinglePoint(double bucketsize)
{
    double* k = new double[1];
    k[0] = 1.0;
    return new Kernel(bucketsize/2.0, k, bucketsize, 1);
}

Kernel* Kernel::Gaussian(double stdev, double bucketsize, double prob)
{
	double rg_end = -NormalCDFInverse((1.0 - prob)/2.0, 0.0, stdev);
	double bucklen = ceil(rg_end - bucketsize/2.0);
	unsigned int buck_offset = static_cast<unsigned int>(bucklen);
	unsigned int buckets = 2 * buck_offset + 1;
	double* k = new double[buckets];
	for (unsigned int ii=0; ii<buckets; ii++)
		k[ii] = NormalPDF((ii - buck_offset) * bucketsize, 0.0, stdev);
	return new Kernel((static_cast<double>(buck_offset) + 0.5) * bucketsize, k, bucketsize, buckets);
}

void Kernel::print()
{
    for( unsigned int ii = 0; ii < buckets; ii++ )
    {
    	std::cout << k[ii] << std::endl;
    }
}


SinglePointFunctionalKernel::SinglePointFunctionalKernel() {}

double SinglePointFunctionalKernel::getMass(double bucketStart, double bucketEnd)
{
	if(bucketStart <= 0.0 and 0.0 < bucketEnd)
		return 1.0;
	return 0.0;
}

double SinglePointFunctionalKernel::getSupportMin()
{ return 0.0; }

double SinglePointFunctionalKernel::getSupportMax()
{ return 0.0; }


TruncatedGaussianFunctionalKernel::TruncatedGaussianFunctionalKernel(double _stdev, double _prob) :
stdev(_stdev), prob(_prob)
{
	support_min = NormalCDFInverse((1.0 - prob)/2.0, 0.0, stdev);
	support_max = -support_max;
//	support_len = support_max - support_min;
	correction = 1.0/prob;
}

double TruncatedGaussianFunctionalKernel::getMass(double bucketStart, double bucketEnd)
{
	double start = std::max(support_min, bucketStart);
	double end   = std::min(support_max, bucketEnd);
	return (NormalCDF(end, 0.0, stdev) - NormalCDF(start, 0.0, stdev)) * correction;
}

double TruncatedGaussianFunctionalKernel::getSupportMin()
{ return support_min; }

double TruncatedGaussianFunctionalKernel::getSupportMax()
{ return support_max; }


RectangularFunctionalKernel::RectangularFunctionalKernel(double start, double end) :
support_min(start), support_max(end)
{
	support_len = support_max - support_min;
}

double RectangularFunctionalKernel::getSupportMin()
{ return support_min; }

double RectangularFunctionalKernel::getSupportMax()
{ return support_max; }


Spectrum::Spectrum(double _start, double _bucketsize, int _buckets, bool _clear)
: start(_start), end(_start + _bucketsize * static_cast<float>(_buckets)), bucketsize(_bucketsize), buckets(_buckets)
{
	spectrum = new double[buckets];
	memset(spectrum, buckets, sizeof(double));
}
	

Spectrum::Spectrum(IsoSpec& iso, FunctionalKernel& kernel, double _bucketsize)
{
	double ker_supp_min = kernel.getSupportMin();
	double ker_supp_max = kernel.getSupportMax();

	double lpm = iso.getLightestPeakMass() + ker_supp_min;
	double hpm = iso.getHeaviestPeakMass() + ker_supp_max;

	start = floor(lpm) - bucketsize * 1.5;
	buckets = floor((hpm - start)/bucketsize) + 2;
	end = start + buckets * bucketsize;

	spectrum = new double[buckets];
	memset(spectrum, 0, buckets * sizeof *spectrum);

	iso.processConfigurationsUntilCutoff();

	double* masses = new double[iso.cnt];
	double* lprobs = new double[iso.cnt];

	iso.getCurrentProduct(masses, lprobs, nullptr);

	unsigned int kernel_bucketstart_offset = static_cast<unsigned int>(ceil(-ker_supp_min/bucketsize)) + 1;
	unsigned int buckets_needed = static_cast<unsigned int>(ceil((ker_supp_max - ker_supp_min)/bucketsize)) + 1;

	for(unsigned int ii = 0; ii < iso.cnt; ii++)
	{
		unsigned int pos = position(masses[ii]);
		unsigned int start_iter = pos - kernel_bucketstart_offset;
		unsigned int end_iter   = start_iter + buckets_needed;

		double current_bucket_start;
		double current_bucket_end   = mass_at_index_start(start_iter);

		double prob = exp(lprobs[ii]);
		
		for(unsigned int jj = start_iter; jj <= end_iter; jj++)
		{
			current_bucket_start = current_bucket_end;
			current_bucket_end += bucketsize;

			spectrum[jj] += kernel.getMass(current_bucket_start, current_bucket_end) * prob;
		}
	}
}




