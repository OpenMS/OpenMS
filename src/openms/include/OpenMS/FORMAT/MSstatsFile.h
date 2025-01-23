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

        /// store label free experiment (MSstats)
        void storeLFQ(const String& filename, 
                      const ConsensusMap &consensus_map, // we might add singleton protein groups
                      const ExperimentalDesign& design,
                      const StringList& reannotate_filenames,
                      const bool is_isotope_label_type,
                      const String& bioreplicate,
                      const String& condition,
                      const String& retention_time_summarization_method);
        
        /// store isobaric experiment (MSstatsTMT)
        void storeISO(const String& filename, 
                      const ConsensusMap &consensus_map,
                      const ExperimentalDesign& design,
                      const StringList& reannotate_filenames,
                      const String& bioreplicate,
                      const String& condition,
                      const String& mixture,
                      const String& retention_time_summarization_method);

    private:
      typedef OpenMS::Peak2D::IntensityType Intensity;
      typedef OpenMS::Peak2D::CoordinateType Coordinate;

      static const String na_string_;
      static const char delim_ = ',';
      static const char accdelim_ = ';';
      static const char quote_ = '"';

      /*
        *  @brief: Struct to aggregate intermediate information from ConsensusFeature and ConsensusMap,
        *  such as filenames, intensities, retention times, labels and features (for further processing)
        */
      struct AggregatedConsensusInfo
      {
        std::vector< std::vector< String > > consensus_feature_filenames;           //< Filenames of ConsensusFeature
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
                                                          const std::vector<String>& spectra_paths);

      /*
        *  @brief: Internal function to check if MSstats_BioReplicate and MSstats_Condition exists in Experimental Design
        */
      static void checkConditionLFQ_(const ExperimentalDesign::SampleSection& sampleSection, const String& bioreplicate, const String& condition);

      /*
        *  @brief: Internal function to check if MSstats_BioReplicate, MSstats_Condition and MSstats_Mixture in Experimental Design
        */
      static void checkConditionISO_(const ExperimentalDesign::SampleSection& sampleSection, const String& bioreplicate, const String& condition, const String& mixture);

      /*
        *  @brief MSstats treats runs differently than OpenMS. In MSstats, runs are an enumeration of (SpectraFilePath, Fraction)
        *  In OpenMS, a run is split into multiple fractions.
        */
      static void assembleRunMap_(
              std::map< std::pair< String, unsigned>, unsigned> &run_map,
              const ExperimentalDesign &design);

      /*
        * @brief checks two vectors for same content
        */
      static bool checkUnorderedContent_(const std::vector< String> &first, const std::vector< String > &second);

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
            const String& _accession,
            const String& _sequence,
            const String& _precursor_charge,
            const String& _fragment_ion,
            const String& _frag_charge,
            const String& _isotope_label_type,
            const String& _condition,
            const String& _bioreplicate,
            const String& _run,
            const String& _fraction
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

        const String& accession() const {return this->accession_;}
        const String& sequence() const {return this->sequence_;}
        const String& precursor_charge() const {return this->precursor_charge_;}
        const String& run() const {return this->run_;}

        String toString() const
        {
          const String delim(",");
          return  accession_
                  + delim + sequence_
                  + delim + precursor_charge_
                  + delim + fragment_ion_
                  + delim + frag_charge_
                  + delim + isotope_label_type_
                  + delim + condition_
                  + delim + bioreplicate_
                  + delim + run_
                  + (this->has_fraction_ ? delim + String(fraction_) : "");
        }

        friend bool operator<(const MSstatsLine_ &l,
                              const MSstatsLine_ &r) {

          return std::tie(l.accession_, l.run_, l.condition_, l.bioreplicate_, l.precursor_charge_, l.sequence_) <
                  std::tie(r.accession_, r.run_, r.condition_, r.bioreplicate_, r.precursor_charge_, r.sequence_);
        }


      private:
        bool has_fraction_;
        String accession_;
        String sequence_;
        String precursor_charge_;
        String fragment_ion_;
        String frag_charge_;
        String isotope_label_type_;
        String condition_;
        String bioreplicate_;
        String run_;
        String fraction_;
      };

      class MSstatsTMTLine_
      {
      public :
        MSstatsTMTLine_(
            const String& _accession,
            const String& _sequence,
            const String& _precursor_charge,
            const String& _channel,
            const String& _condition,
            const String& _bioreplicate,
            const String& _run,
            const String& _mixture,
            const String& _techrepmixture,
            const String& _fraction
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

        const String& accession() const {return this->accession_;}
        const String& sequence() const {return this->sequence_;}
        const String& precursor_charge() const {return this->precursor_charge_;}
        const String& run() const {return this->run_;}

        String toString() const
        {
          const String delim(",");
          return  accession_
                  + delim + sequence_
                  + delim + precursor_charge_
                  + delim + channel_
                  + delim + condition_
                  + delim + bioreplicate_
                  + delim + run_
                  + delim + mixture_
                  + delim + techrepmixture_
                  + delim + String(fraction_);
        }

        friend bool operator<(const MSstatsTMTLine_ &l,
                              const MSstatsTMTLine_ &r) {

          return std::tie(l.accession_, l.run_, l.condition_, l.bioreplicate_, l.mixture_, l.precursor_charge_, l.sequence_, l.channel_) <
                  std::tie(r.accession_, r.run_, r.condition_, r.bioreplicate_, r.mixture_, r.precursor_charge_, r.sequence_, r.channel_);
        }


      private:
        String accession_;
        String sequence_;
        String precursor_charge_;
        String channel_;
        String condition_;
        String bioreplicate_;
        String run_;
        String mixture_;
        String techrepmixture_;
        String fraction_;
      };

      /*
        *  @brief Constructs the lines and adds them to the TextFile
        *  @param peptideseq_quantifyable Has to be a set (only) for deterministic  ordered output
        */
      template <class LineType>
      void constructFile_(const String& retention_time_summarization_method,
                          const bool rt_summarization_manual,
                          TextFile& csv_out,
                          const std::set<String>& peptideseq_quantifyable,
                          LineType & peptideseq_to_prefix_to_intensities) const;

      /*
      *  @brief Constructs the accession to indist. group mapping
      */
      static std::unordered_map<OpenMS::String, const IndProtGrp* > getAccessionToGroupMap_(const IndProtGrps& ind_prots);


      /*
       * @brief Based on the evidence accession set in a PeptideHit, checks if is unique and therefore quantifyable
       * in a group context.
       *
       */
      bool isQuantifyable_(
          const std::set<String>& accs,
          const std::unordered_map<String, const IndProtGrp*>& accession_to_group) const;

    };
} // namespace OpenMS
