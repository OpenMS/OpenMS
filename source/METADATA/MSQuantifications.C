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
#include<set>
#include<iostream>

using namespace std;

namespace OpenMS
{
		const std::string MSQuantifications::NamesOfQuantTypes[] = {"MS1LABEL", "MS2LABEL", "LABELFREE"};

			/// Constructor
			MSQuantifications::MSQuantifications() :
				ExperimentalSettings()
			{
			}

			/// Copy constructor
			MSQuantifications::MSQuantifications(const MSQuantifications& source) :
				ExperimentalSettings(source)
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
				//~ PersistentObject::operator=(source);

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

			void MSQuantifications::setDataProcessingList(std::vector<DataProcessing>& dpl)
			{
				data_processings_ = dpl;
			}

			const std::vector<DataProcessing> MSQuantifications::getDataProcessingList() const
			{
				std::vector<DataProcessing> list = data_processings_;

				//This is one way street for dataprocessing - it probably wont get mapped back after writeout and readin
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

			const std::vector<MSQuantifications::Assay>& MSQuantifications::getAssays() const
			{
				return assays_;
			}

			std::vector<MSQuantifications::Assay>& MSQuantifications::getAssays()
			{
				return assays_;
			}
			
			//~ std::map<String,ConsensusFeature::Ratio>& MSQuantifications::getRatioCalculations()
			//~ {
				//~ return ratio_calculations_;
			//~ }

			const std::vector<FeatureMap<> >& MSQuantifications::getFeatureMaps() const
			{
				return feature_maps_;
			}

			const std::vector<ConsensusMap>& MSQuantifications::getConsensusMaps() const
			{
				return consensus_maps_;
			}
			
			std::vector<ConsensusMap>& MSQuantifications::getConsensusMaps()
			{
				return consensus_maps_;
			}

			const MSQuantifications::AnalysisSummary& MSQuantifications::getAnalysisSummary() const
			{
				return analysis_summary_;
			}

			MSQuantifications::AnalysisSummary& MSQuantifications::getAnalysisSummary()
			{
				return analysis_summary_;
			}

			void MSQuantifications::setAnalysisSummaryQuantType(MSQuantifications::QUANT_TYPES r)
			{
				analysis_summary_.quant_type_ = r;
			}

			void MSQuantifications::addConsensusMap(ConsensusMap& m)
			{
				consensus_maps_.push_back(m);
			}

			void MSQuantifications::assignUIDs()
			{
				for(std::vector<Assay>::iterator ait = assays_.begin(); ait != assays_.end(); ++ait)
				{
					ait->uid_ = String(UniqueIdGenerator::getUniqueId());
				}
			}

			void MSQuantifications::registerExperiment(MSExperiment<Peak1D> & exp, std::vector< std::vector< std::pair<String, DoubleReal> > > label)
			{
				for (std::vector< std::vector< std::pair<String, DoubleReal> > >::const_iterator lit = label.begin(); lit != label.end(); ++lit)
				{
					//TODO look for existing labels
					Assay a;
					a.mods_ = (*lit);
					a.raw_files_.push_back(exp.getExperimentalSettings());
					assays_.push_back(a);
				}

				data_processings_ = exp[0].getDataProcessing(); //todo overwrite MSExperiments inherited front method to work. [0] operator is ugly!
			}

}//namespace OpenMS

