// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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


#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletTransform.h>
#include <math.h>

namespace OpenMS
{
/** @brief Computes the discrete-time continuous wavelet transform simultaneously for several charges.
	* 
	* The function returns a pointer to the transformed versions of @p scan. Note that after having called this 
	* function the deletion of this pointer is up to you! The function computes the transform for several charge states
	* at the same time.
	* With the help of @p jiggle_bins it is possible to compute an "interpolated" transform, s.t. the maximum of
	* wavelet is not exactly aligned with some data point. In some (few) cases this could help to find the exact 
	* m/z position of some feature, since the wavelet even "samples" points not included in the scan. 
	* Usually, this parameter should be left unchanged, since it involves some expert knowledge to interpret 
	* the results. Note that the size of each MSSpectrum that the function returns increased of @p jiggle_bins 
	* is strictly positive.
	* This function is slightly less efficient than the non static version, since it has to do some preprocessing
	* on the MS scan. In particular the function computes the average MZ spacing of @scan. Hence, you might get 
	* (slightly) different results if you set up an instance of this class and provide this average
	* from yourself. In the latter case you would prefer to compute the average spacing over the whole map instead
	* on a single scan. 
	*
	* @param scan The MS scan you wish to transform.
	* @param max_charge The maximal charge state that is considered.
	* @param jiggle_bins The number of additional "virtual sampling points" between each pair of signal points.
	* @return The transformed spectra. Entry i in this vector corresponds to the "i+1"-charge-state-transform 
	* of @p scan. The length of each MSSpectrum contained in this vector equals either the length of @p scan.
	* If @param scan is smaller than the internally computed wavelet length, no transform can be computed and 
	* the original scan will be returend for each charge state. */
std::vector<MSSpectrum<RawDataPoint1D> >* IsotopeWaveletTransform::getTransforms 
	(const MSSpectrum<RawDataPoint1D>& scan, const unsigned int max_charge, const double jiggle_bins) throw ()
{	
	unsigned int scan_size = scan.size();
			
	double av_MZ_spacing=0;
	for (unsigned int i=0; i<scan_size-1; ++i)
		 av_MZ_spacing += scan[i+1].getMZ() - scan[i].getMZ();
	av_MZ_spacing /= (double) (scan_size-1);

	unsigned int peak_cutoff = IsotopeWavelet::getPeakCutoff();
	unsigned int wavelet_length = (unsigned int) trunc(peak_cutoff/av_MZ_spacing);	

	//Creating the result vector
	std::vector<MSSpectrum<RawDataPoint1D> >* res = new std::vector<MSSpectrum<RawDataPoint1D> > (max_charge, scan);
	std::vector<double> psi (wavelet_length, 0); //The wavelet

	if (scan_size < wavelet_length)
		return (res); //original scan will be returned

	
	double cum_spacing, c_spacing, //Helping variables
		max_w_monoi_intens=0.25*NEUTRON_MASS, //The position of the monoisotopic peak within the coordinate sys. of the wavelet 
		sums=0, //Helping variables
		max_position_scan=0, //The position of the data point (within the scan) we want to align with
		align_offset, //Correction term; shifts the wavelet to get the desired alignment
		last, distance; //Helping variables
	unsigned int c=0, k=0, j=0;
	double c_charge; //double, since we will oven devide by c_charge 
	
	//The upcoming variable is necessary to capture strange effects in special types of unequally spaced data sets.
	//Imagine some wholes in the m/z range (points the mass spectrometer did not sample). If they become larger than 
	//0.25*NEUTRON_MASS (considering the case of charge 1), several data points wll share the same max_position, 
	//causing the upcoming code to crash since suddenly some m/z pom/z positions willl occure twice. The interval of multiple 
	//occuring points is stored by multiple_s and implicitly by i.
	std::vector<int> multiple_s (max_charge,-1);
	std::vector<double> last_max_position_scan (max_charge, -1);
	bool repair=false;

	unsigned int c_res_index=0; //this is the index of the convolution function
	//Starting convolution
	for (unsigned int i=0; i<scan_size; ++i)
	{
		//Now, let's sample the wavelets
		for (c=0; c<max_charge; ++c)
		{	
			c_charge=c+1;
			cum_spacing=0;				
			max_w_monoi_intens=0.25*NEUTRON_MASS/c_charge; //This is the position of the monoistopic peak (centered)
			
			//Align the maximum monoisotopic peak of the wavelet with some scan point. This is step is critical, since
			//otherwise we might - especially in the case of badly resolved data - miss patterns, since scan maxima and
			//wavelet maxima are "anticorrelated"
			j=0; last=0;
			while (cum_spacing < max_w_monoi_intens)
			{
				c_spacing = scan[(i+j+1)%scan_size].getMZ() - scan[(i+j)%scan_size].getMZ();
			 	last=cum_spacing;	
				if (c_spacing < 0) //I.e. we are at the end of the scan
					cum_spacing += av_MZ_spacing;
				else //The "normal" case
					cum_spacing += c_spacing; 
				++j;
			};
				
			align_offset = max_w_monoi_intens - last; //I.e. we have to shift the wavelet by this amount to align the data
			--j;				

			//The upcoming variable holds the position of the spectrum that is aligned with the monoisotopic 
			//maximum of the wavelet. We do not add the overall correction term for the left shift at this point, 
			//since we will get trouble by the NEUTRON_MASS and the resulting numerical instabilities. 
			//We will add this correcting term at the end of the whole processing.
			max_position_scan = scan[(i+j)%scan_size].getMZ();
			if (max_position_scan == last_max_position_scan[c]) //Uuups, multiple times the some m/z coordinate
			{
				if (multiple_s[c] < 0) //This is the first entry where this artifact occured
					multiple_s[c] = i-1; 
				//Notice that the problematic case of multiple_s being at the end of the spectrum (this might happen for 
				//the overlapping part of the transform) can be ignored. 				
				//The special case if we are the boundary (exactly the last point in the spectrum).
				if (i == scan_size-1)
					repair = true;
			}
			else //Denotes the end of the multiple pos interval and triggers a repair.
			{
				if (multiple_s[c] >= 0)
					repair=true; //We cannot do this now. Just after the transform at the actual point is completed.
			};

			last_max_position_scan[c] = max_position_scan;

			//We devide the interval between i+j and i+j-1 in equally spaced bins and try to jiggle the wavelet
			//a little bit around in this region to get a better fit for the maximum position. 	
			//If this makes sense??? Could not find any improvements of the results by now. Therefore, jiggle_bins=0 by
			//default.	
			if (jiggle_bins > 0)
				distance = (scan[(i+j)%scan_size].getMZ() - scan[(i+j-1)%scan_size].getMZ()) / jiggle_bins;
			else
				distance = 0;

			

			//The loop for refinement; will be executed only once, if there is no refinement.
			int count = (jiggle_bins<=0) ? -1: 0; //We need this correction, to avoid double evalution in the refinement steps.
			cum_spacing = align_offset;
			for (double c_jig=0; count<jiggle_bins; c_jig += distance, ++count)
			{
				//std::cout << i << "\t before sampling" << std::endl;
				//Sampling the wavelet 
				sampleTheWavelet (scan, i, cum_spacing, (unsigned int) c_charge, av_MZ_spacing, psi);
				//std::cout << i << "\t after sampling" << std::endl;

				
				k=0; sums=0;
				for (unsigned int j=i; j<scan_size && k<wavelet_length; ++j, ++k)
					sums += scan[j].getIntensity()*psi[k];

				if (k< wavelet_length) // I.e. we have an overlapping wavelet
					sums=0; // => We can absolutely neglect this feature since it is too near at the boundary.

				//Store the current convolution result
				if (jiggle_bins <= 0) //Normal correlation, no interpolating step
				{
					(*res)[c][i].setIntensity(sums);
					(*res)[c][i].setMZ(max_position_scan);	
		
					if (repair) //Only necessary in this case, not in the refinement
					{		
						unsigned int noi2interpol = i - multiple_s[c]; //NOT +1
		
						//The special case if we are the boundary (exactly the last point in the spectrum)
						if (i == scan_size-1)
						{
							//We do not care about the intenities, since we will set them to zero anyway.
							//We would just like to avoid multiple positions to occur in the transform
							for (unsigned int ii=0; ii<=noi2interpol; ++ii) 
							//it must be "<=noi..." !!! not "<", since in this case we do not want to keep the last multiple	
							//the same holds for "ii=0"
							{
								(*res)[c][multiple_s[c]+ii].setMZ((*res)[c][multiple_s[c]-1].getMZ() + (ii+1)*av_MZ_spacing);		
							};
							
							last_max_position_scan[c] = max_position_scan; //Reset
							multiple_s[c]=-1; //Reset
							repair=false;
							continue;
						};
						
						double x1 = (*res)[c][multiple_s[c]].getMZ();
						double y1 = (*res)[c][multiple_s[c]].getIntensity();					
						double x2 = (*res)[c][i].getMZ();
						if (i >= scan_size) //Is still possible and ugly => reset x2
							x2 = ((*res)[c][i].getMZ() - (*res)[c][i-1].getMZ()) + (*res)[c][i].getMZ();
						double y2 = (*res)[c][i].getIntensity();
						if (i >= scan_size)
							y2 = (*res)[c][i].getIntensity(); //Do just anything, does not matter what, since we are at the boundary
						double dx = (x2-x1)/(double)(noi2interpol);
						for (unsigned int ii=1; ii<noi2interpol; ++ii) //ii=1, not 0, since we want to keep the first of the multiples
						{	
							(*res)[c][multiple_s[c]+ii].setMZ((*res)[c][multiple_s[c]].getMZ()+ii*dx);
							(*res)[c][multiple_s[c]+ii].setIntensity(y1 + (y2-y1)/(x2-x1)*((*res)[c][multiple_s[c]].getMZ()+ii*dx-x1));
						};
					
						last_max_position_scan[c] = max_position_scan; //Reset
						multiple_s[c]=-1; //Reset
						repair=false;
					};				
				}
				else //Refinement
				{	
					(*res)[c][c_res_index].setIntensity(sums);
					(*res)[c][c_res_index].setMZ(max_position_scan);
					++c_res_index;	
				};				
				cum_spacing = align_offset;
				cum_spacing -= c_jig;		
				max_position_scan += distance;	
			};
		};
	};
	return (res); //user carries over the responsibility for deleting this object!
}


/** @brief Samples the wavelet at discrete time points, s.t. they match automatically to the m/z positions provided
	* in @p scan and returns the discrete values of psi in @psi. Usually (unless you would like to write your own 
	* @see FeatureFinder), you do not need to call this function.
	*	@param mz_index The start index of @p scan for which the the wavelet should be adapted.  
	* @param offset Is the offset the wavelet function needs to be aligned with a scan point.
	* @param z The z the wavelet function should adapt (corresponds to z in the paper).
	* @param av_MZ_spacing The average spacing in the m/z domain.
	* @param psi The sampled values.
	* @param mode Indicates wheter positive mode (+1) or negative mode (-1) has been used for ionization. */ 
void IsotopeWaveletTransform::sampleTheWavelet (const MSSpectrum<RawDataPoint1D>& scan, const unsigned int mz_index, 
	const double offset, const unsigned int z, const double av_MZ_spacing, std::vector<double>& psi, const unsigned int mode)
	throw ()
{
	unsigned int scan_size = scan.size(), help;
	double c_pos, c_pos1, lambda, c_spacing;

	psi.resize (scan_size); //just to be sure; if psi is already scan_size large, this will is a simple test
	
	c_pos = scan[mz_index].getMZ();				
	lambda = IsotopeWavelet::getLambdaL(c_pos*z-mode*z*PROTON_MASS);
	unsigned int peak_cutoff = IsotopeWavelet::getPeakCutoff();
	unsigned int wavelet_length = (unsigned int) trunc(peak_cutoff/av_MZ_spacing);	

	double cum_spacing=offset;
	//Building up (sampling) the wavelet
	for (unsigned int j=0; j<wavelet_length; ++j)
	{
		c_pos = scan[(mz_index+j)%scan_size].getMZ();	
		c_pos1 = scan[(mz_index+j+1)%scan_size].getMZ();
		psi[j] = ((cum_spacing<=0) ? 0 : IsotopeWavelet::getValueByLambda (cum_spacing, lambda, z));

		c_spacing = c_pos1 - c_pos;
		//c_spacing might get negative, as soon as the wavelet approaches the end of the scan (if i=scan_size-1).
		//Since this case is only of theoretical interest (we do not expect any important scans at the very end of the 
		//spectrum), we simlply use the average spacing in that case.
		if (c_spacing < 0)
			cum_spacing += av_MZ_spacing;
		else //The "normal" case
			cum_spacing += c_spacing;			
	};

	double mean =0;
	for (unsigned int j=0; j<wavelet_length-1; ++j)
	{
		help = mz_index+j;
		mean += chordTrapezoidRule (scan[help%scan_size].getMZ(), scan[(help+1)%scan_size].getMZ(), 
			psi[j], psi[j+1]);
	};

	//Subtracting the av_MZ_spacing
	for (unsigned int j=0; j<wavelet_length; ++j)
		psi[j] -= mean/(double)peak_cutoff;
}

} //namespace
