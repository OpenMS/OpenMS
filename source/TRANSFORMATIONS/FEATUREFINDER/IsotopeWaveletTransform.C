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
#include <OpenMS/CHEMISTRY/Averagine.h>
#include <math.h>

namespace OpenMS
{
	
bool comparator (const RawDataPoint1D a, const RawDataPoint1D b)
{
	return (a.getIntensity() > b.getIntensity());
}		


IsotopeWaveletTransform::IsotopeWaveletTransform () throw() : hash_precision_ (DEFAULT_HASH_PRECISION)
{ }


IsotopeWaveletTransform::~IsotopeWaveletTransform () throw()
{ }

		
void IsotopeWaveletTransform::getTransforms (const MSSpectrum<RawDataPoint1D>& scan, 
	std::vector<MSSpectrum<RawDataPoint1D> > &transforms, const unsigned int max_charge) throw ()
{	
	unsigned int scan_size = scan.size();
	double av_MZ_spacing = getAvMZSpacing(scan);
	unsigned int peak_cutoff = IsotopeWavelet::getPeakCutoff();
	unsigned int wavelet_length = (unsigned int) trunc(peak_cutoff/av_MZ_spacing);	
	std::vector<double> psi (wavelet_length, 0); //The wavelet

	if (scan_size < wavelet_length)
		return;

	
	double cum_spacing, c_spacing, //Helping variables
		max_w_monoi_intens=0.25*NEUTRON_MASS, //The position of the monoisotopic peak within the coordinate sys. of the wavelet 
		sums=0, //Helping variables
		max_position_scan=0, //The position of the data point (within the scan) we want to align with
		align_offset, //Correction term; shifts the wavelet to get the desired alignment
		last;
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
			if (max_position_scan == last_max_position_scan[c]) //Uuups, multiple times the same m/z coordinate
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
			cum_spacing = align_offset;
				
			//Sampling the wavelet 
			sampleTheWavelet (scan, i, cum_spacing, (unsigned int) c_charge, av_MZ_spacing, psi);
				
			//The convolution
			k=0; sums=0;
			for (unsigned int j=i; j<scan_size && k<wavelet_length; ++j, ++k)
				sums += scan[j].getIntensity()*psi[k];

			if (k< wavelet_length) // I.e. we have an overlapping wavelet
			{
				sums=0; // => We can absolutely neglect this feature since it is too near at the boundary.
				max_position_scan = transforms[c][i-1].getMZ()+av_MZ_spacing;
			};

			//Store the current convolution result
			transforms[c][i].setIntensity(sums);
			transforms[c][i].setMZ(max_position_scan);	

			if (repair)
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
						transforms[c][multiple_s[c]+ii].setMZ(transforms[c][multiple_s[c]-1].getMZ() + (ii+1)*av_MZ_spacing);		
					};

					last_max_position_scan[c] = max_position_scan; //Reset
					multiple_s[c]=-1; //Reset
					repair=false;
					continue;
				}

				double x1 = transforms[c][multiple_s[c]].getMZ();
				double y1 = transforms[c][multiple_s[c]].getIntensity();					
				double x2 = transforms[c][i].getMZ();
				if (i >= scan_size) //Is still possible and ugly => reset x2
					x2 = (transforms[c][i].getMZ() - transforms[c][i-1].getMZ()) + transforms[c][i].getMZ();
				double y2 = transforms[c][i].getIntensity();
				if (i >= scan_size)
					y2 = transforms[c][i].getIntensity(); //Do just anything, does not matter what, since we are at the boundary
				double dx = (x2-x1)/(double)(noi2interpol);
				for (unsigned int ii=1; ii<noi2interpol; ++ii) //ii=1, not 0, since we want to keep the first of the multiples
				{	
					transforms[c][multiple_s[c]+ii].setMZ(transforms[c][multiple_s[c]].getMZ()+ii*dx);
					transforms[c][multiple_s[c]+ii].setIntensity(y1 + (y2-y1)/(x2-x1)*(transforms[c][multiple_s[c]].getMZ()+ii*dx-x1));
				};

				last_max_position_scan[c] = max_position_scan; //Reset
				multiple_s[c]=-1; //Reset
				repair=false;
			}				
		}
	}
	return;
}


void IsotopeWaveletTransform::identifyCharges (const std::vector<MSSpectrum<RawDataPoint1D> >& candidates, const unsigned int scan_index, 	
	const double ampl_cutoff) throw ()
{	
	double av_MZ_spacing = getAvMZSpacing(candidates[0]);
	unsigned int peak_cutoff = IsotopeWavelet::getPeakCutoff();
	unsigned int wavelet_length = (unsigned int) trunc(peak_cutoff/av_MZ_spacing);	
	unsigned int cands_size = candidates.size();
	unsigned int signal_size=candidates[0].size(), c_index, i_iter, end_index; 
	int start_index; //Do not change this to unsigned int
	MSSpectrum<RawDataPoint1D>::iterator iter;
	double seed_mz, c_av_intens, c_score, c_sd_intens, threshold;

	//For all charges do ...
	for (unsigned int c=0; c<cands_size; ++c)		
	{
		//Indicates wheter for some specific region of the scan an isotopic pattern has already been identified
		//i.e.: In the moment, we do not care about overlapping signals! 
		std::vector<bool> processed = std::vector<bool> (signal_size, false); 
		MSSpectrum<RawDataPoint1D> c_sorted_candidate = candidates[c]; 

		//The fowllong hash map allows a fast transform from m/z positions to m/z indices w.r.t. the transformed vector
		hash_multimap<int, unsigned int> index_hash;
		for (unsigned int i=0; i<signal_size; ++i)
		{
			int hash_key = (int) trunc(candidates[c][i].getMZ()*hash_precision_); 
			index_hash.insert (std::pair<int, unsigned int> (hash_key, i));
		};

		//Sort the transform in descending order according to the intensities present in the transform 	
		sort (c_sorted_candidate.begin(), c_sorted_candidate.end(), comparator); 
		c_av_intens = getAvIntens (candidates[c]);		
		c_sd_intens = getSdIntens (candidates[c], c_av_intens);

		//Eliminate uninteresting regions
		//In principle we should do that in a binary fashion for efficiency reasons ...  
		for (iter=c_sorted_candidate.begin(); iter != c_sorted_candidate.end(); ++iter)
		{
			if (iter->getIntensity() <= c_av_intens)
				break;
		};
		c_sorted_candidate.erase (iter, c_sorted_candidate.end());	

		if (ampl_cutoff < 0)
			threshold = 0;
		else
			threshold=ampl_cutoff*c_sd_intens + c_av_intens;

		i_iter=0;
		for (iter=c_sorted_candidate.begin(); iter != c_sorted_candidate.end(); ++iter, ++i_iter)
		{					
			seed_mz=iter->getMZ();
			c_index = index_hash.find((int)trunc(seed_mz*hash_precision_))->second;

			if (processed[c_index])
        continue;
				
			//In order to determine the start and end indices, we first need to know the width of the region one should consider 
			//to estimate the mean and the sd of the pattern candidate. 
			//That region is defined by the position of the heighst amplitude +/- wavelet_length_.
			start_index = c_index-wavelet_length-1;
			end_index = c_index+wavelet_length+1;

			if (isinf(start_index)) //Error check, should be removed after some intensive code tests 
				std::cout << "Ups:  start_index is inf" << std::endl; 

			if (start_index < 0)
				start_index = 0;
			if (end_index >= signal_size)
				end_index = signal_size-1;			
			
			//Mark as processed
			for (unsigned int z=start_index; z<=end_index; ++z)
				processed[z] = true;	
			
			c_score = scoreThis (candidates[c], start_index, seed_mz, end_index, c, iter->getIntensity(), threshold);
	
			if (c_score <= 0)
				continue;

			//Push the seed into its corresponding box (or create a new one, if necessary)
			push2Box (seed_mz, scan_index, c, c_score, iter->getIntensity(), candidates[0].getRT());
		};	
	};
}
	

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


double IsotopeWaveletTransform::scoreThis (const MSSpectrum<RawDataPoint1D>& candidate, unsigned int start_index, 
	const double seed_mz, const unsigned int end_index, const unsigned int c, const double intens, const double ampl_cutoff) throw ()
{	
	unsigned int cands_size = candidate.size(); 
	double c_av_MZ_spacing = getAvMZSpacing (candidate, start_index, end_index);
	unsigned int peak_cutoff = IsotopeWavelet::getPeakCutoff();
	double c_score=0, c_check_point=-1, c_val;
	std::pair<int, int> c_between_left, c_between_right; 				

	std::vector<double> xs, ys;		
	double leftBound, rightBound; unsigned int interpol_range;
	//p_h_ind indicates if we are looking for a whole or a peak
	int p_h_ind=1, end=4*peak_cutoff -1; //4 times and not 2 times, since we move by 0.5 m/z entities

	std::vector<double> xvec, yvec, weights;

	for (int v=1; v<end; ++v, ++p_h_ind)
	{
		c_check_point = seed_mz-(peak_cutoff*NEUTRON_MASS-v*0.5*NEUTRON_MASS)/((double)c+1);

		leftBound = c_check_point;
		rightBound = c_check_point;		

		do
		{ leftBound -= c_av_MZ_spacing;
			c_between_left = getNearBys (candidate, leftBound, start_index);
		}	while (leftBound >= candidate[0].getMZ() && (c_between_left.first < 0 || c_between_left.second < 0)); 

		do 
		{ rightBound += c_av_MZ_spacing;
			c_between_right = getNearBys (candidate, rightBound, start_index);
		} while (rightBound <=  candidate[cands_size-1].getMZ() && (c_between_right.first < 0 || c_between_right.second < 0)); 

		interpol_range = c_between_right.second - c_between_left.first +1; 


		if (interpol_range < 3) //We need at least 3 points for cubic interpolation
		{
			--c_between_left.first; 
			++c_between_right.first;
		};
		
		if (c_between_left.first <= 0 || c_between_left.second <= 0 || c_between_right.first <= 0 || c_between_right.second <= 0)
		{		
			--end; //I.e. we also drop a checkpoint at the right side of the interval
			continue;
		};
		
		interpol_range = c_between_right.second - c_between_left.first +1; 
		xs.resize(interpol_range); ys.resize(interpol_range);
		
		for (int i=c_between_left.first, j=0; i<=c_between_right.second; ++i, ++j)
		{			
			xs[j] = candidate[i].getMZ();	
			ys[j] = candidate[i].getIntensity();
		};

		c_val = getCubicInterpolatedValue (xs, c_check_point, ys);

		if (p_h_ind%2 == 1) //I.e. a whole
			c_score -= c_val;
		else
			c_score +=c_val;
	};

	if (c_score <= ampl_cutoff+intens)
		return(0);
	
	return (log(c_score)+ log(intens));
}


double IsotopeWaveletTransform::getAvMZSpacing (const MSSpectrum<RawDataPoint1D>& scan, int start_index, int end_index) throw ()
{ 
	double av_MZ_spacing=0;
	if (end_index < 0)
		end_index = scan.size();
	for (int i=start_index; i<end_index-1; ++i)
		 av_MZ_spacing += scan[i+1].getMZ() - scan[i].getMZ();
	return (av_MZ_spacing / (double) (end_index-1-start_index));
}


double IsotopeWaveletTransform::getAvIntens (const MSSpectrum<RawDataPoint1D>& scan) throw ()
{ 
	double av_intens=0;
	for (unsigned int i=0; i<scan.size(); ++i)
		 av_intens += scan[i].getIntensity();
	return (av_intens / (double) (scan.size()));
}


double IsotopeWaveletTransform::getSdIntens (const MSSpectrum<RawDataPoint1D>& scan, const double mean) throw ()
{
	double res=0, intens;
	for (unsigned int i=0; i<scan.size(); ++i)
	{
		intens = (scan[i].getIntensity() < 0) ? 0 : scan[i].getIntensity();
		res += (intens-mean)*(intens-mean);
	};

	return (sqrt(res/(double)(scan.size()-1)));	
}


double IsotopeWaveletTransform::getCubicInterpolatedValue (const std::vector<double>& x, const double xi, const std::vector<double>& y) throw ()
{
	gsl_interp_accel* acc = gsl_interp_accel_alloc ();
	gsl_spline* spline = gsl_spline_alloc (gsl_interp_cspline, x.size());

	gsl_spline_init (spline, &x[0], &y[0], x.size());
	double yi = gsl_spline_eval (spline, xi, acc);

	gsl_spline_free (spline);
	gsl_interp_accel_free (acc);
	return (yi);	
}


std::pair<int, int> IsotopeWaveletTransform::getNearBys (const MSSpectrum<RawDataPoint1D>& signal, const double mz, 
	const unsigned int start) const throw ()
{
	for (unsigned int i=start; i<signal.size(); ++i)
	{
		if (signal[i].getMZ() > mz)
		{
			if (i>start) //everything's fine
				return (std::pair<int, int> (i-1, i));
			else //wrong boundaries!
				break;
		};
	};

	//not found
	return (std::pair<int, int> (-1, -1));
}


void IsotopeWaveletTransform::push2Box (const double mz, const unsigned int scan, unsigned int charge, 
	const double score, const double intens, const double rt) throw ()
{	
	std::map<double, Box>::iterator upper_iter = openBoxes_.upper_bound(mz);
	std::map<double, Box>::iterator lower_iter; 
	if (openBoxes_.empty())
		lower_iter = openBoxes_.end();
	else
		lower_iter = openBoxes_.lower_bound(mz);

	//Ugly, but necessary due to the implementation of STL lower_bound
	if (mz != openBoxes_.lower_bound(mz)->first && lower_iter != openBoxes_.begin())
		lower_iter = --(openBoxes_.lower_bound(mz));

	std::map<double, Box>::iterator insert_iter;
	bool createNewBox=false;
	if (lower_iter == openBoxes_.end()) //I.e. there is no open Box for that mz position
		createNewBox=true;

	if (upper_iter == openBoxes_.end() && fabs(lower_iter->first - mz) < 0.5*NEUTRON_MASS) //Found matching Box
	{
		insert_iter = lower_iter;
		createNewBox=false;
	}
	else
		createNewBox=true;

	if (upper_iter != openBoxes_.end() && lower_iter != openBoxes_.end())
	{	
		//Figure out which entry is closer to m/z
		double dist_lower = fabs(lower_iter->first - mz);
		double dist_upper = fabs(upper_iter->first - mz);
		dist_lower = (dist_lower < 0.5*NEUTRON_MASS) ? dist_lower : INT_MAX;
		dist_upper = (dist_upper < 0.5*NEUTRON_MASS) ? dist_upper : INT_MAX;

		if (dist_lower>=0.5*NEUTRON_MASS && dist_upper>=0.5*NEUTRON_MASS) // they are both too far away
			createNewBox=true;
		else
		{
			insert_iter = (dist_lower < dist_upper) ? lower_iter : upper_iter;	
			createNewBox=false;
		};
	}; 

	BoxElement element; element.c = charge; element.mz = mz; element.score = score; element.RT = rt, element.intens=intens;
	std::pair<unsigned int, BoxElement> help2 (scan, element);

	if (createNewBox == false)
	{	
		insert_iter->second.insert (help2);	

		//Unfortunately, we need to change the m/z key to the average of all keys inserted in that box.
		Box replacement (insert_iter->second);	

		//We cannot devide both m/z by 2, since we already inserted some m/z's whose weight would be lowered.
		//Also note that we alread inserted the new entry, leading to size-1.
		double c_mz = insert_iter->first * (insert_iter->second.size()-1) + mz;	
		c_mz /= ((double) insert_iter->second.size());			
		
		//Now let's remove the old and insert the new one
		openBoxes_.erase (insert_iter);	
		std::pair<double, std::map<unsigned int, BoxElement> > help3 (c_mz, replacement);	
		openBoxes_.insert (help3);		
	}
	else
	{
		std::map<unsigned int, BoxElement> help3;
		help3.insert (help2);
		std::pair<double, std::map<unsigned int, BoxElement> > help4 (mz, help3);
		openBoxes_.insert (help4);
	};
}


void IsotopeWaveletTransform::updateBoxStates (const unsigned int c_scan_number, const unsigned int RT_interleave, 
	const unsigned int RT_votes_cutoff) throw ()
{
	std::map<double, Box>::iterator iter, iter2;
	for (iter=openBoxes_.begin(); iter!=openBoxes_.end(); )
	{
		//For each Box we need to figure out, if and when the last RT value has been inserted
		//If the Box his unchanged since RT_interleave_ scans, we will close the Box.
		unsigned int lastScan = (--(iter->second.end()))->first;
		if (c_scan_number - lastScan > RT_interleave) //I.e. close the box!
		{
			iter2 = iter;
			++iter2;
			if (iter->second.size() >= RT_votes_cutoff)
				closedBoxes_.insert (*iter);
			openBoxes_.erase (iter);
			iter=iter2;
		}
		else
			++iter;
	};
}


FeatureMap<Feature> IsotopeWaveletTransform::mapSeeds2Features 
	(const unsigned int max_charge, const unsigned int RT_votes_cutoff) throw ()
{
	FeatureMap<Feature> feature_map;
	std::list<RawDataPoint2D> thelist;
	std::map<double, Box>::iterator iter;
	Box::iterator box_iter;
	unsigned int bestChargeIndex; double bestChargeScore, c_mz, c_RT; unsigned int c_charge; 		
	ConvexHull2D c_conv_hull;
	Feature c_feature;

	std::cout << "Size of closed boxes: " << closedBoxes_.size() << std::endl;
	std::pair<double, double> c_extend;
	for (iter=closedBoxes_.begin(); iter!=closedBoxes_.end(); ++iter)
	{		
		Box& c_box = iter->second;
		std::vector<double> chargeVotes (max_charge, 0), chargeBinaryVotes (max_charge, 0);
	
		//Let's first determine the charge
		//Therefor, we can use two types of votes: qulitative ones (chargeBinaryVotes) or quantitaive ones (chargeVotes)
		//Collting the votes ...
		for (box_iter=iter->second.begin(); box_iter!=iter->second.end(); ++box_iter)
		{
			chargeVotes[box_iter->second.c] += box_iter->second.score;
			++chargeBinaryVotes[box_iter->second.c];
		};
		
		//... dertermining the best fitting charge
		bestChargeIndex=0; bestChargeScore=0; 
		for (unsigned int i=0; i<max_charge; ++i)
			if (chargeVotes[i] > bestChargeScore)
			{
				bestChargeIndex = i;
				bestChargeScore = chargeVotes[i];
			};			

		//Pattern found in too few RT scan 
		if (chargeBinaryVotes[bestChargeIndex] < RT_votes_cutoff)
			continue;

		c_charge = bestChargeIndex + 1; //that's the finally predicted charge state for the pattern

		double av_intens=0, av_mz=0;
		//Now, let's get the RT boundaries for the box
		std::vector<DPosition<2> > point_set;
		for (box_iter=c_box.begin(); box_iter!=c_box.end(); ++box_iter)
		{
			c_mz = box_iter->second.mz;
			c_RT = box_iter->second.RT;
			Averagine::getModel (c_mz, c_charge, &c_extend);
			std::cout << box_iter->first << " (" << c_RT << ")\t" << c_extend.first << "\t" 
				<< c_extend.second << "\t" << box_iter->second.c +1 << "#" << c_charge << std::endl;
			point_set.push_back (DPosition<2> (c_extend.first, c_RT)); //mz start, RT
			point_set.push_back (DPosition<2> (c_extend.second, c_RT)); //mz start, RT
			av_intens += box_iter->second.intens;
			av_mz += c_mz;
		};
		av_mz /= (double)c_box.size();
		av_intens /= (double)c_box.size();	
		std::cout << "**************************************************" << std::endl;

		c_conv_hull = point_set;
		c_feature.setConvexHulls (std::vector<ConvexHull2D> (1, c_conv_hull));
		c_feature.setMZ (av_mz);
		c_feature.setIntensity (av_intens);
		feature_map.push_back (c_feature);
	};		
	return (feature_map);
}

} //namespace
