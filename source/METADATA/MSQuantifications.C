// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Mathias Walzer$
// $Authors: Mathias Walzer$
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/MSQuantifications.h>

using namespace std;

namespace OpenMS
{
			/// Constructor
			MSQuantifications::MSQuantifications() :
				ExperimentalSettings(),
				PersistentObject()
			{
			}

			/// Copy constructor
			MSQuantifications::MSQuantifications(const MSQuantifications& source) :
				ExperimentalSettings(source),
				PersistentObject(source)
			{
			}
			
			MSQuantifications::~MSQuantifications()
			{  	
			}
			
			/// Assignment operator
			MSQuantifications& MSQuantifications::operator= (const MSQuantifications& source)
			{
				if (&source == this) return *this;

				ExperimentalSettings::operator=(source);
				PersistentObject::operator=(source);

				//~ reassign members
				
				return *this;
			}

			/// Equality operator
			bool MSQuantifications::operator== (const MSQuantifications& rhs) const
			{
				return ExperimentalSettings::operator==(rhs);
			}
			
			/// Equality operator
			bool MSQuantifications::operator!= (const MSQuantifications& rhs) const
			{
				return !(operator==(rhs));
			}
			
			std::vector<DataProcessing> MSQuantifications::getDataProcessingList() const
			{
				std::vector<DataProcessing> list; 
				for (std::vector<FeatureMap<> >::const_iterator fit = feature_maps_.begin(); fit != feature_maps_.end(); ++fit)
				{
					list.insert(list.end(),fit->getDataProcessing().begin(), fit->getDataProcessing().end());
				}

				for (std::vector<ConsensusMap>::const_iterator cit = consensus_maps_.begin(); cit != consensus_maps_.end(); ++cit)
				{
					list.insert(list.end(),cit->getDataProcessing().begin(), cit->getDataProcessing().end());
				}
				
				return list;				
			}
			
			std::vector<MSQuantifications::Assay> MSQuantifications::getAssays() const
			{
				return assays_;
			}

			std::vector<FeatureMap<> > MSQuantifications::getFeatureMaps() const
			{
				return feature_maps_;
			}			
			
			std::vector<ConsensusMap> MSQuantifications::getConsensusMaps() const
			{
				return consensus_maps_;
			}			

			MSQuantifications::AnalysisSummary MSQuantifications::getAnalysisSummary() const
			{
				return analysis_summary_;
			}

}//namespace OpenMS

