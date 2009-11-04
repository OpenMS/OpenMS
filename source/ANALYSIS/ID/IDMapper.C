// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDMapper.h>

using namespace std;

namespace OpenMS 
{

	IDMapper::IDMapper()
  	: DefaultParamHandler("IDMapper"),
  		rt_delta_(5.0),
  		mz_delta_(1),
  		measure_(MEASURE_PPM)
  {
		defaults_.setValue("rt_delta",rt_delta_, "allowed RT delta in seconds"); 
		defaults_.setMinFloat("rt_delta",0);
		defaults_.setValue("mz_delta", mz_delta_, "allowed m/z delta in ppm or Da"); 
		defaults_.setMinFloat("mz_delta",0);
		defaults_.setValue("mz_measure", "ppm", "unit of mz_delta (ppm or Da)"); 
		defaults_.setValidStrings("mz_measure", StringList::create("ppm,Da"));
		defaultsToParam_();  
  }

			
	IDMapper::IDMapper(const IDMapper& cp)
	: DefaultParamHandler(cp),
		rt_delta_(cp.rt_delta_),
  	mz_delta_(cp.mz_delta_),
  	measure_(cp.measure_)
	{
		updateMembers_();
	}

	IDMapper& IDMapper::operator = (const IDMapper& rhs)
	{
		if (this == &rhs) return *this;
		
		DefaultParamHandler::operator = (rhs);
		rt_delta_=rhs.rt_delta_;
  	mz_delta_=rhs.mz_delta_;
  	measure_=rhs.measure_;
		updateMembers_();
		
		return *this;
	}

	void IDMapper::updateMembers_()
	{
		rt_delta_ = param_.getValue("rt_delta");
		mz_delta_ = param_.getValue("mz_delta");
		measure_ = param_.getValue("mz_measure")=="ppm"? MEASURE_PPM : MEASURE_DA;
	}	
	
	void IDMapper::annotate(ConsensusMap& map, const std::vector<PeptideIdentification>& ids, const std::vector<ProteinIdentification>& protein_ids, bool measure_from_subelements)
	{
		checkHits_(ids);
				
		//append protein identifications to Map
		map.getProteinIdentifications().insert(map.getProteinIdentifications().end(),protein_ids.begin(),protein_ids.end());

		//keep track of assigned/unassigned peptide identifications
		std::set<Size> assigned;
					
		// store which peptides fit which feature (and avoid double entries)
		// consensusMap -> {peptide_index}
		std::vector < std::set< size_t> > mapping(map.size());
		
		//iterate over the peptide IDs
		for (Size i=0; i<ids.size(); ++i)
		{
			if (ids[i].getHits().size()==0) continue;

			DoubleReal rt_pep = ids[i].getMetaValue("RT");
			DoubleReal mz_pep = ids[i].getMetaValue("MZ");

			//iterate over the features
			for(Size cm_index = 0 ; cm_index<map.size(); ++cm_index)
			{
				//check if we compare distance from centroid or subelements
				if (!measure_from_subelements)
				{
					if ( isMatch_(rt_pep-map[cm_index].getRT(), mz_pep, map[cm_index].getMZ()) )
					{
						map[cm_index].getPeptideIdentifications().push_back(ids[i]);
						assigned.insert(i);
					}
				}
				else
				{
					for(ConsensusFeature::HandleSetType::const_iterator it_handle = map[cm_index].getFeatures().begin(); 
							it_handle != map[cm_index].getFeatures().end(); 
							++it_handle)
					{
						if (isMatch_(rt_pep - it_handle->getRT(), mz_pep,it_handle->getMZ()))
						{
							if (mapping[cm_index].count(i) == 0)
							{
								map[cm_index].getPeptideIdentifications().push_back(ids[i]);
								assigned.insert(i);
								mapping[cm_index].insert(i);
							}
							continue; // we added this peptide already.. no need to check further
						}
					}
				}
			}
		}

		//append unassigned peptide identifications
		for (Size i=0; i<ids.size(); ++i)
		{
			if (assigned.count(i)==0)
			{
				map.getUnassignedPeptideIdentifications().push_back(ids[i]);
			}
		}

	}

	const DoubleReal IDMapper::getAbsoluteMZDelta_(const DoubleReal mz) const
	{
		if (measure_==MEASURE_PPM)
		{
			return (mz * mz_delta_ / 1e6);
		} 
		else if (measure_==MEASURE_DA)
		{
			return mz_delta_;
		}
		throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "IDMapper::getAbsoluteDelta_(): illegal internal state of measure_!", String(measure_));
	}

	const bool IDMapper::isMatch_(const DoubleReal rt_distance, const DoubleReal mz_theoretical, const DoubleReal mz_observed) const
	{
		if (measure_==MEASURE_PPM)
		{
			return (fabs(rt_distance) <= rt_delta_) && (fabs((mz_theoretical*1e6-mz_observed*1e6)/mz_theoretical) <= mz_delta_);
		} 
		else if (measure_==MEASURE_DA)
		{
			return (fabs(rt_distance) <= rt_delta_) && (fabs(mz_theoretical-mz_observed) <= mz_delta_);
		}
		throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "IDMapper::getAbsoluteDelta_(): illegal internal state of measure_!", String(measure_));
	}

	void IDMapper::checkHits_(const std::vector<PeptideIdentification>& ids) const
	{
		for (Size i=0; i<ids.size(); ++i)
		{
			if (!ids[i].metaValueExists("RT"))
			{
				throw Exception::MissingInformation(__FILE__,__LINE__,__PRETTY_FUNCTION__, "IDMapper: meta data value 'RT' missing for peptide identification!"); 
			}
			if (!ids[i].metaValueExists("MZ"))
			{
				throw Exception::MissingInformation(__FILE__,__LINE__,__PRETTY_FUNCTION__, "IDMapper: meta data value 'MZ' missing for peptide identification!"); 
			}
		}
	}

} // namespace OpenMS
