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


template <typename MapType>
IsotopeWaveletFF<MapType>* IsotopeWaveletFF<MapType>::me_ = NULL;


template <typename MapType>		
IsotopeWaveletFF<MapType>::IsotopeWaveletFF() throw()
	: experiment_ (NULL), max_charge_ (1), threshold_ (0), RT_votes_cutoff_ (1), RT_interleave_(0), hash_precision_(1000)
{ }

template <typename MapType>	
IsotopeWaveletFF<MapType>::~IsotopeWaveletFF() throw()
{ }


template <typename MapType>		
IsotopeWaveletFF<MapType>::IsotopeWaveletFF (MapType& experiment, const unsigned int max_charge, const double threshold, 
	const unsigned int RT_votes_cutoff, const unsigned int RT_interleave, const unsigned int hash_precision) throw ()
	: experiment_ (experiment), max_charge_ (max_charge), threshold_ (threshold), RT_votes_cutoff_ (RT_votes_cutoff), 
	RT_interleave_(RT_interleave), hash_precision_(hash_precision)
{ }


template <typename MapType>		
FeatureMap<Feature> IsotopeWaveletFF<MapType>::runFF (int start_scan, int end_scan) throw ()
{
	IsotopeWaveletTransform iwt;

	if (end_scan < 0)
	{
		end_scan = experiment_.size();
	};

	for (unsigned int i=start_scan; i<end_scan; ++i)
	{	
		std::vector<MSSpectrum<RawDataPoint1D> > pwts (max_charge_, experiment_[i]);
		std::cout << "Spectrum " << i << " (" << ::std::fixed << ::std::setprecision(4) << 
			experiment_[i].getRT() << ") of " << end_scan-1 << " ... " ; 
		std::cout.flush();
		
		IsotopeWaveletTransform::getTransforms (experiment_[i], pwts, max_charge_);
	
		std::cout << "transform ok ... "; std::cout.flush();
		iwt.identifyCharges (pwts, i, threshold_);
		std::cout << "charge recognition ok ... "; std::cout.flush();
		iwt.updateBoxStates(i, RT_interleave_, RT_votes_cutoff_);
		std::cout << "updated box states." << std::endl;
	};
	
	//And now ... a cute hack ;-) 
	//Forces to empty OpenBoxes_ and to synchronize ClosedBoxes_ 
	iwt.updateBoxStates(INT_MAX, RT_interleave_, RT_votes_cutoff_); 

	std::cout << "Final mapping."; std::cout.flush();
	return (iwt.mapSeeds2Features (max_charge_, RT_votes_cutoff_));
}
