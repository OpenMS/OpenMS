// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <map>
#include <utility>
#include <unordered_map>
#include <set>
#include <vector>

namespace OpenMS
{
  using IndProtGrp = OpenMS::ProteinIdentification::ProteinGroup;
  using IndProtGrps = std::vector<IndProtGrp>;

  /**
    @brief File adapter for Triqler files
    @ingroup FileIO
  */
  
    class OPENMS_DLLAPI TriqlerFile
    {
    public:
        /// Default constructor
        TriqlerFile() = default;
        /// Destructor
        ~TriqlerFile() = default;

        /// store label free experiment
        void storeLFQ(const std::string& filename, 
                      const ConsensusMap &consensus_map,
                      const ExperimentalDesign& design,
                      const StringList& reannotate_filenames,
                      const std::string& condition);
        
    private:
      typedef OpenMS::Peak2D::IntensityType Intensity;
      typedef OpenMS::Peak2D::CoordinateType Coordinate;

      static const std::string na_string_;
      static const char delim_ = ',';
      static const char accdelim_ = ';';
      static const char quote_ = '"';

      /*
        *  @brief: Struct to aggregate intermediate information from ConsensusFeature and ConsensusMap,
        *  such as filenames, intensities, retention times, labels and features (for further processing)
        */
      struct AggregatedConsensusInfo
      {
        std::vector< std::vector< std::string> > consensus_feature_filenames;           //< Filenames of ConsensusFeature
        std::vector< std::vector< Intensity > > consensus_feature_intensities;       //< Intensities of ConsensusFeature
        std::vector< std::vector< Coordinate > > consensus_feature_retention_times; //< Retention times of ConsensusFeature
        std::vector< std::vector< unsigned > > consensus_feature_labels;          //< Labels of ConsensusFeature
        std::vector<BaseFeature> features;                                        //<s Features of ConsensusMap
      };

      /*
        *  @brief: Aggregates information from ConsensusFeature and ConsensusMap,
        *  such as filenames, intensities, retention times, labels and features.
        *  Stores them in AggregatedConsensusInfo for later processing
        */
      TriqlerFile::AggregatedConsensusInfo aggregateInfo_(const ConsensusMap& consensus_map,
                                                          const std::vector<std::string>& spectra_paths);

      /*
        *  @brief: Internal function to check if condition exists in Experimental Design
        */
      static void checkConditionLFQ_(const ExperimentalDesign::SampleSection& sampleSection, const std::string& condition);

      /*
        *  In OpenMS, a run is split into multiple fractions.
        */
      static void assembleRunMap_(
              std::map< std::pair< std::string, unsigned>, unsigned> &run_map,
              const ExperimentalDesign &design);

      /*
        * @brief checks two vectors for same content
        */
      static bool checkUnorderedContent_(const std::vector< std::string> &first, const std::vector< std::string> &second);

      OpenMS::Peak2D::IntensityType sumIntensity_(const std::set< OpenMS::Peak2D::IntensityType > &intensities) const
      {
        OpenMS::Peak2D::IntensityType result = 0;
        for (const OpenMS::Peak2D::IntensityType &intensity : intensities)
        {
          result += intensity;
        }
        return result;
      }

      OpenMS::Peak2D::IntensityType meanIntensity_(const std::set< OpenMS::Peak2D::IntensityType > &intensities) const
      {
        return sumIntensity_(intensities) / intensities.size();
      }

      class TriqlerLine_
      {
      public :
        TriqlerLine_(
            const std::string& run,
            const std::string& condition,
            const std::string& precursor_charge,
            const std::string& search_score,           
            const std::string& intensity,
            const std::string& sequence,
            const std::string& accession
        ):  run_(run),
            condition_(condition),
            precursor_charge_(precursor_charge),
            search_score_(search_score),
            intensity_(intensity),
            sequence_(sequence),
            accession_(accession)
             {}
        
        TriqlerLine_(TriqlerLine_&& m) = default;

        TriqlerLine_(const TriqlerLine_& m) = default;

        /// as string
        std::string toString() const;

        friend bool operator<(const TriqlerLine_ &l,
                              const TriqlerLine_ &r) 
        {
          return std::tie(l.accession_, l.run_, l.condition_, l.precursor_charge_, l.intensity_, l.sequence_) <
                 std::tie(r.accession_, r.run_, r.condition_, r.precursor_charge_, r.intensity_, r.sequence_);
        }

      private:
        std::string run_;
        std::string condition_;
        std::string precursor_charge_;
        std::string search_score_;
        std::string intensity_;
        std::string sequence_;
        std::string accession_;
      };

      using MapSequenceToLines_ = std::map<std::string, std::set<TriqlerLine_>>;
      /*
        *  @brief Constructs the lines and adds them to the TextFile
        *  @param[out] peptideseq_quantifyable Has to be a set (only) for deterministic ordered output
        */
      void constructFile_(TextFile& csv_out,
                          const std::set<std::string>& peptideseq_quantifyable,
                          const MapSequenceToLines_& peptideseq_to_line) const;

      /*
      *  @brief Constructs the accession to indist. group mapping
      */
      static std::unordered_map<std::string, const IndProtGrp* > getAccessionToGroupMap_(const IndProtGrps& ind_prots);


      /*
       * @brief Based on the evidence accession set in a PeptideHit, checks if is unique and therefore quantifyable
       * in a group context.
       *
       */
      bool isQuantifyable_(
          const std::set<std::string>& accs,
          const std::unordered_map<std::string, const IndProtGrp*>& accession_to_group) const;

    };
} // namespace OpenMS
