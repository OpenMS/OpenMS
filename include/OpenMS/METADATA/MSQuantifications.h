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

#ifndef OPENMS_METADATA_MSQUANTIFICATIONS_H
#define OPENMS_METADATA_MSQUANTIFICATIONS_H

#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
//~ #include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/DataProcessing.h>

#include <vector>
#include <map>

namespace OpenMS
{
	class OPENMS_DLLAPI MSQuantifications
		:	public ExperimentalSettings
	{
		public:
			/// @name Base type definitions
			//@{
			/// typedef docu
			typedef CVTermList ParamGroupList; // userparams are exclusively inside the CVTermList's MetaInfoInterface

			enum QUANT_TYPES {MS1LABEL=0, MS2LABEL, LABELFREE, SIZE_OF_QUANT_TYPES}; // derived from processing applied
			static const std::string NamesOfQuantTypes[SIZE_OF_QUANT_TYPES];
			
			//@}
			//~ InputFiles: //~ searchdb abbildung version,releasedate,#entries,dbname über paramgrouplist
			//~ struct ParamGroupList
			//~ {
				//~ ParamGroupList()
				//~ {
				//~ }

				//~ ParamGroupList(const ParamGroupList& rhs)
					//~ :	cv_params(rhs.cv_params)
				//~ {
				//~ }

				//~ ~ParamGroupList()
				//~ {
				//~ }

				//~ ParamGroupList& operator = (const ParamGroupList& rhs)
				//~ {
					//~ if (&rhs != this)
					//~ {
						//~ cv_params = rhs.cv_params;
						//~ user_params = rhs.user_params;
					//~ }
					//~ return *this;
				//~ }

				//~ MetaInfoInterface user_params;
				//~ CVTermList cv_params;
			//~ };

			struct AnalysisSummary
			{
				AnalysisSummary()
				{
				}

				AnalysisSummary(const AnalysisSummary& rhs)
					:	cv_params_(rhs.cv_params_)
				{
					user_params_ = rhs.user_params_;
					quant_type_ = rhs.quant_type_;
				}

				virtual ~AnalysisSummary()
				{
				}

				AnalysisSummary& operator = (const AnalysisSummary& rhs)
				{
					if (&rhs != this)
					{
						cv_params_ = rhs.cv_params_;
						user_params_ = rhs.user_params_;
						quant_type_ = rhs.quant_type_;
					}
					return *this;
				}

				MetaInfo user_params_;
				CVTermList cv_params_;
				QUANT_TYPES quant_type_;
			};

			struct Assay
			{
				Assay()
				{
				}

				Assay(const Assay& rhs)
				{
					uid_ = rhs.uid_;
					mods_ = rhs.mods_;
					raw_files_ = rhs.raw_files_;
					feature_maps_ = rhs.feature_maps_;
				}

				virtual ~Assay()
				{
				}

				Assay& operator = (const Assay& rhs)
				{
					if (&rhs != this)
					{
						uid_ = rhs.uid_;
						mods_ = rhs.mods_;
						raw_files_ = rhs.raw_files_;
						feature_maps_ = rhs.feature_maps_;
					}
					return *this;
				}

				String uid_;
				std::vector< std::pair<String, DoubleReal> > mods_;
				std::vector<ExperimentalSettings> raw_files_;
				std::map<size_t, FeatureMap<> > feature_maps_; // iTRAQ needs no FeatureMaps so ExperimentalSettings are not directly mapped to FeatureMaps
			};
			
			// TODO handle referencing from consensusmaps to featuremaps/rawfiles
			// TODO add ContactPerson or something to (Consensus)FeatureMap or DataProcessing (see below)
			// TODO rewrite OpenMS::DataProcessing - data not yet linked in openms core formats - below should go in analysissummary of MSQuantifications - input/output not possible to be carried along
			//~ if(DataProcessing::NamesOfProcessingAction[*it] == String("Quantitation"))
					//~ {
					//~ if (processing.getSoftware().getName()==String("SILACAnalyzer"))
					//~ {
							//~ experiment_type = MS1LABEL;
					//~ }
					//~ else if (processing.getSoftware().getName()==String("ITRAQAnalyzer"))
					//~ {
							//~ experiment_type = MS2LABEL;
					//~ }
					//~ else
					//~ {
							//~ experiment_type = LABELFREE;
					//~ }
			//~ }
			//~ QUANT_TYPES experiment_type = MS1LABEL;

			/// Constructor
			MSQuantifications();

			/// Destructor
			~MSQuantifications();

			/// Copy constructor
			MSQuantifications(const MSQuantifications& source);

			/// Assignment operator
			MSQuantifications& operator= (const MSQuantifications& source);

		/// Equality operator
			bool operator== (const MSQuantifications& rhs) const;

			/// Equality operator
			bool operator!= (const MSQuantifications& rhs) const;

			/**
				@brief Loads data from a text file.

				@param filename The input file name.
				@param trim_lines Whether or not the lines are trimmed when reading them from file.
				@param first_n If set, only @p first_n lines the lines from the beginning of the file are read.

				@note this function uses unix-style linebreaks

				@exception Exception::FileNotFound is thrown if the file could not be opened.
			*/
			void load(const String& filename, bool trim_lines=false, Int first_n=-1);

			const std::vector<DataProcessing> getDataProcessingList() const;
			const std::vector<Assay>& getAssays() const;
			std::vector<Assay>& getAssays();
			std::map<String,ConsensusFeature::Ratio>& getRatios();
			const std::vector<ConsensusMap>& getConsensusMaps() const;
			std::vector<ConsensusMap>& getConsensusMaps();
			const std::vector<FeatureMap<> >& getFeatureMaps() const;
			const AnalysisSummary& getAnalysisSummary() const;
			AnalysisSummary& getAnalysisSummary();
			void setDataProcessingList(std::vector<DataProcessing>& dpl);
			void setAnalysisSummaryQuantType(QUANT_TYPES r);
			void addConsensusMap(ConsensusMap& m);
			void assignUIDs();
			void registerExperiment(MSExperiment<Peak1D> & exp, std::vector< std::vector< std::pair<String, DoubleReal> > > labels);

		private:
			AnalysisSummary analysis_summary_;
			std::vector<MetaInfo> bibliographic_reference_;
			std::vector<ConsensusMap> consensus_maps_;
			std::vector<FeatureMap<> > feature_maps_;
			std::vector<Assay> assays_;
			std::vector<DataProcessing> data_processings_;
			//~ std::map<String,ConsensusFeature::Ratio > ratio_calculations_;
		};

} // namespace OpenMS

#endif // OPENMS_METADATA_MSQUANTIFICATIONS_H
