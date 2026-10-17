// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Lukas Heumos $
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
    @brief File adapter for MSstats files
    @ingroup FileIO
  */

    class OPENMS_DLLAPI MSstatsFile
    {
    public:
        /// Default constructor
        MSstatsFile() = default;
        /// Destructor
        ~MSstatsFile() = default;

        /**
         * @brief Store label free experiment (MSstats)
         * @param[in] filename Output filename
         * @param[in] consensus_map ConsensusMap with quantification data
         * @param[in] design ExperimentalDesign file
         * @param[in] reannotate_filenames Optional filenames for reannotation
         * @param[in] is_isotope_label_type If true, IsotopeLabelType is 'H', else 'L'
         * @param[in] bioreplicate Column name for biological replicate in design
         * @param[in] condition Column name for condition in design
         * @param[in] retention_time_summarization_method Method for RT summarization
         * @param[in] remove_shared_peptides If true, shared peptides mapping to multiple indistinguishable protein groups are removed (default)
         */
        void storeLFQ(const std::string& filename,
                      const ConsensusMap &consensus_map, // we might add singleton protein groups
                      const ExperimentalDesign& design,
                      const StringList& reannotate_filenames,
                      const bool is_isotope_label_type,
                      const std::string& bioreplicate,
                      const std::string& condition,
                      const std::string& retention_time_summarization_method,
                      const bool remove_shared_peptides = true);

        /**
         * @brief Store isobaric experiment (MSstatsTMT)
         * @param[in] filename Output filename
         * @param[in] consensus_map ConsensusMap with quantification data
         * @param[in] design ExperimentalDesign file
         * @param[in] reannotate_filenames Optional filenames for reannotation
         * @param[in] bioreplicate Column name for biological replicate in design
         * @param[in] condition Column name for condition in design
         * @param[in] mixture Column name for mixture in design (used for TMT experiments)
         * @param[in] retention_time_summarization_method Method for RT summarization
         * @param[in] remove_shared_peptides If true, shared peptides mapping to multiple indistinguishable protein groups are removed (default)
         */
        void storeISO(const std::string& filename,
                      const ConsensusMap &consensus_map,
                      const ExperimentalDesign& design,
                      const StringList& reannotate_filenames,
                      const std::string& bioreplicate,
                      const std::string& condition,
                      const std::string& mixture,
                      const std::string& retention_time_summarization_method,
                      const bool remove_shared_peptides = true);

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
        std::vector< std::vector< std::string > > consensus_feature_filenames;           //< Filenames of ConsensusFeature
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
      MSstatsFile::AggregatedConsensusInfo aggregateInfo_(const ConsensusMap& consensus_map,
                                                          const std::vector<std::string>& spectra_paths);

      /*
        *  @brief: Internal function to check if MSstats_BioReplicate and MSstats_Condition exists in Experimental Design
        */
      static void checkConditionLFQ_(const ExperimentalDesign::SampleSection& sampleSection, const std::string& bioreplicate, const std::string& condition);

      /*
        *  @brief: Internal function to check if MSstats_BioReplicate, MSstats_Condition and MSstats_Mixture in Experimental Design
        */
      static void checkConditionISO_(const ExperimentalDesign::SampleSection& sampleSection, const std::string& bioreplicate, const std::string& condition, const std::string& mixture);

      /*
        *  @brief MSstats treats runs differently than OpenMS. In MSstats, runs are an enumeration of (SpectraFilePath, Fraction)
        *  In OpenMS, a run is split into multiple fractions.
        */
      static void assembleRunMap_(
              std::map< std::pair< std::string, unsigned>, unsigned> &run_map,
              const ExperimentalDesign &design);

      /*
        * @brief checks if the first vector is a subset of the second
        */
      static bool isSubsetOf_(const std::vector< std::string> &first, const std::vector< std::string > &second);
      static void warnOnSubsetFiles_(const std::vector<std::string>& spectra_paths, const std::vector<std::string>& design_filenames);

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

      class MSstatsLine_
      {
      public :
        MSstatsLine_(
            bool _has_fraction,
            const std::string& _accession,
            const std::string& _sequence,
            const std::string& _precursor_charge,
            const std::string& _fragment_ion,
            const std::string& _frag_charge,
            const std::string& _isotope_label_type,
            const std::string& _condition,
            const std::string& _bioreplicate,
            const std::string& _run,
            const std::string& _fraction
        ): has_fraction_(_has_fraction),
            accession_(_accession),
            sequence_(_sequence),
            precursor_charge_(_precursor_charge),
            fragment_ion_(_fragment_ion),
            frag_charge_(_frag_charge),
            isotope_label_type_(_isotope_label_type),
            condition_(_condition),
            bioreplicate_(_bioreplicate),
            run_(_run),
            fraction_(_fraction) {}

        const std::string& accession() const {return this->accession_;}
        const std::string& sequence() const {return this->sequence_;}
        const std::string& precursor_charge() const {return this->precursor_charge_;}
        const std::string& run() const {return this->run_;}

        std::string toString() const
        {
          const std::string delim(",");
          return  accession_
                  + delim + sequence_
                  + delim + precursor_charge_
                  + delim + fragment_ion_
                  + delim + frag_charge_
                  + delim + isotope_label_type_
                  + delim + condition_
                  + delim + bioreplicate_
                  + delim + run_
                  + (this->has_fraction_ ? delim + std::string(fraction_) : "");
        }

        friend bool operator<(const MSstatsLine_ &l,
                              const MSstatsLine_ &r) {

          return std::tie(l.accession_, l.run_, l.condition_, l.bioreplicate_, l.precursor_charge_, l.sequence_) <
                  std::tie(r.accession_, r.run_, r.condition_, r.bioreplicate_, r.precursor_charge_, r.sequence_);
        }


      private:
        bool has_fraction_;
        std::string accession_;
        std::string sequence_;
        std::string precursor_charge_;
        std::string fragment_ion_;
        std::string frag_charge_;
        std::string isotope_label_type_;
        std::string condition_;
        std::string bioreplicate_;
        std::string run_;
        std::string fraction_;
      };

      class MSstatsTMTLine_
      {
      public :
        MSstatsTMTLine_(
            const std::string& _accession,
            const std::string& _sequence,
            const std::string& _precursor_charge,
            const std::string& _channel,
            const std::string& _condition,
            const std::string& _bioreplicate,
            const std::string& _run,
            const std::string& _mixture,
            const std::string& _techrepmixture,
            const std::string& _fraction
        ): accession_(_accession),
            sequence_(_sequence),
            precursor_charge_(_precursor_charge),
            channel_(_channel),
            condition_(_condition),
            bioreplicate_(_bioreplicate),
            run_(_run),
            mixture_(_mixture),
            techrepmixture_(_techrepmixture),
            fraction_(_fraction) {}

        const std::string& accession() const {return this->accession_;}
        const std::string& sequence() const {return this->sequence_;}
        const std::string& precursor_charge() const {return this->precursor_charge_;}
        const std::string& run() const {return this->run_;}

        std::string toString() const
        {
          const std::string delim(",");
          return  accession_
                  + delim + sequence_
                  + delim + precursor_charge_
                  + delim + channel_
                  + delim + condition_
                  + delim + bioreplicate_
                  + delim + run_
                  + delim + mixture_
                  + delim + techrepmixture_
                  + delim + std::string(fraction_);
        }

        friend bool operator<(const MSstatsTMTLine_ &l,
                              const MSstatsTMTLine_ &r) {

          return std::tie(l.accession_, l.run_, l.condition_, l.bioreplicate_, l.mixture_, l.precursor_charge_, l.sequence_, l.channel_) <
                  std::tie(r.accession_, r.run_, r.condition_, r.bioreplicate_, r.mixture_, r.precursor_charge_, r.sequence_, r.channel_);
        }


      private:
        std::string accession_;
        std::string sequence_;
        std::string precursor_charge_;
        std::string channel_;
        std::string condition_;
        std::string bioreplicate_;
        std::string run_;
        std::string mixture_;
        std::string techrepmixture_;
        std::string fraction_;
      };

      /*
        *  @brief Constructs the lines and adds them to the TextFile
        *  @param[out] peptideseq_quantifyable Has to be a set (only) for deterministic  ordered output
        */
      template <class LineType>
      void constructFile_(const std::string& retention_time_summarization_method,
                          const bool rt_summarization_manual,
                          TextFile& csv_out,
                          const std::set<std::string>& peptideseq_quantifyable,
                          LineType & peptideseq_to_prefix_to_intensities) const;

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
