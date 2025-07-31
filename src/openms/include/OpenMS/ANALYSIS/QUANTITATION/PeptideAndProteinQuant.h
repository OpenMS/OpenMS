// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>


namespace OpenMS
{
  /**
      @brief Helper class for peptide and protein quantification based on feature data annotated with IDs.

      This class is used by @ref TOPP_ProteinQuantifier. See there for further documentation.

      @htmlinclude OpenMS_PeptideAndProteinQuant.parameters
  */
  class OPENMS_DLLAPI PeptideAndProteinQuant :
    public DefaultParamHandler
  {
public:

    /// Mapping: sample ID -> abundance
    typedef std::map<UInt64, double> SampleAbundances;

    /// Quantitative and associated data for a peptide
    struct PeptideData
    {
      /// mapping: fraction -> filename -> charge -> channel/label -> abundance
      std::map<Int, std::map<String, std::map<Int, std::map<Int, double>>>> abundances;

      /// mapping: fraction -> filename -> charge -> abundance
      std::map<Int, std::map<String, std::map<Int, UInt64>>> psm_counts;

      /// mapping: sample -> total abundance
      SampleAbundances total_abundances;

      /// spectral counting-based abundances
      SampleAbundances total_psm_counts;

      /// protein accessions for this peptide
      std::set<String> accessions;

      /// total number of identifications
      Size psm_count = 0;

      /// constructor
      PeptideData() = default;
    };

    /// Mapping: peptide sequence (modified) -> peptide data
    typedef std::map<AASequence, PeptideData> PeptideQuant;

    /// Quantitative and associated data for a protein
    struct ProteinData
    {
      /// mapping: peptide (unmodified) -> sample -> abundance
      std::map<String, SampleAbundances> peptide_abundances;

      std::map<String, SampleAbundances> peptide_psm_counts;

      /// mapping: filename -> channel/label -> abundance
      std::map<String, std::map<Int, double>> channel_level_abundances;

      /// mapping: filename -> PSM counts
      std::map<String, UInt64> file_level_psm_counts;

      /// mapping: sample -> total abundance
      SampleAbundances total_abundances;

      /// spectral counting-based abundances
      SampleAbundances total_psm_counts;

      /// number of distinct peptide sequences
      SampleAbundances total_distinct_peptides;

      /// total number of PSMs mapping to this protein
      Size psm_count = 0;

      /// constructor
      ProteinData() = default;
    };

    /// Mapping: protein accession -> protein data
    typedef std::map<String, ProteinData> ProteinQuant;

    /// Statistics for processing summary
    struct Statistics
    {
      /// number of samples (or assays in mzTab terms)
      Size n_samples;

      /// number of fractions
      Size n_fractions;

      /// number of MS files
      Size n_ms_files;

      /// protein statistics
      Size quant_proteins, too_few_peptides;

      /// peptide statistics
      Size quant_peptides, total_peptides;

      /// feature statistics
      Size quant_features, total_features, blank_features, ambig_features;

      /// constructor
      Statistics() :
        n_samples(0), quant_proteins(0), too_few_peptides(0),
        quant_peptides(0), total_peptides(0), quant_features(0),
        total_features(0), blank_features(0), ambig_features(0) {}
    };

    /// Constructor
    PeptideAndProteinQuant();

    /// Destructor
    ~PeptideAndProteinQuant() override {}

    /**
         @brief Read quantitative data from a feature map.

         Parameters should be set before using this method, as setting parameters will clear all results.
    */
    void readQuantData(FeatureMap& features, const ExperimentalDesign& ed);

    /**
         @brief Read quantitative data from a consensus map.

         Parameters should be set before using this method, as setting parameters will clear all results.
    */
    void readQuantData(ConsensusMap& consensus, const ExperimentalDesign& ed);

    /**
         @brief Read quantitative data from identification results (for quantification via spectral counting).

         Parameters should be set before using this method, as setting parameters will clear all results.
    */
    void readQuantData(std::vector<ProteinIdentification>& proteins,
                       PeptideIdentificationList& peptides,
                       const ExperimentalDesign& ed);

    /**
         @brief Compute peptide abundances.

         Based on quantitative data for individual charge states (in member @p pep_quant_), overall abundances for peptides are computed (and stored again in @p pep_quant_).

         Quantitative data must first be read via readQuantData().

         Optional (peptide-level) protein inference information (e.g. from Fido or ProteinProphet) can be supplied via @p peptides. In that case, peptide-to-protein associations - the basis for protein-level quantification - will also be read from @p peptides!
    */
    void quantifyPeptides(const PeptideIdentificationList& peptides =
                          PeptideIdentificationList());


    /**
         @brief Compute protein abundances.

         Peptide abundances must be computed first with quantifyPeptides(). Optional protein inference information (e.g. BasicProteinInference or Epifany) can be supplied via @p proteins.
         
         @param proteins Optional protein inference information
    */
    void quantifyProteins(const ProteinIdentification& proteins = ProteinIdentification());


    std::map<OpenMS::String, OpenMS::String> mapAccessionToLeader(const OpenMS::ProteinIdentification& proteins) const;

    /// Get summary statistics
    const Statistics& getStatistics();

    /// Get peptide abundance data
    const PeptideQuant& getPeptideResults();

    /// Get protein abundance data
    const ProteinQuant& getProteinResults();

    /// Annotate protein quant results as meta data to protein ids
    void annotateQuantificationsToProteins(
      const ProteinQuant& protein_quants, 
      ProteinIdentification& proteins,
      bool remove_unquantified = true);

private:

    /// Processing statistics for output in the end
    Statistics stats_;

    /// Peptide quantification data
    PeptideQuant pep_quant_;

    /// Protein quantification data
    ProteinQuant prot_quant_;

    /// Experimental design for filename/channel to sample mapping
    ExperimentalDesign experimental_design_;


    /**
         @brief Get the "canonical" annotation (a single peptide hit) of a feature/consensus feature from the associated list of peptide identifications.

         Only the best-scoring peptide hit of each ID in @p peptides is taken into account. The hits of each ID must already be sorted! If there's more than one ID and the best hits are not identical by sequence, or if there's no peptide ID, an empty peptide hit (for "ambiguous/no annotation") is returned.
         Protein accessions from identical peptide hits are accumulated.
    */
    PeptideHit getAnnotation_(PeptideIdentificationList& peptides);

    /**
         @brief Gather quantitative information from a feature.

         Store quantitative information from @p feature in member @p pep_quant_, based on the peptide annotation in @p hit.
         @p fraction, use 0 for first fraction (or if no fractionation was performed)
         @p filename, the base filename (without path/extension) from which the feature originates
         @p channel_or_label, the channel/label identifier (e.g., TMT channel, typically 1 for LFQ)
         If @p hit is empty ("ambiguous/no annotation"), nothing is stored.
    */
    void quantifyFeature_(const FeatureHandle& feature,
      size_t fraction,
      const String& filename,
      const PeptideHit& hit,
      Int channel_or_label);

    /**
     *   @brief Determine fraction, filename, charge state, and channel of a peptide with the highest
     *   number of abundances.
     *   @param peptide_abundances Const input map fraction -> filename -> charge -> channel -> abundance
     *   @param best Will additionally return the best fraction, filename, charge state, and channel
     *   @return true if at least one abundance was found, false otherwise
     */
    bool getBest_(
      const std::map<Int, std::map<String, std::map<Int, std::map<Int, double>>>> & peptide_abundances,
      std::tuple<size_t, String, size_t, Int> & best);

    /**
         @brief Order keys (charges/peptides for peptide/protein quantification) according to how many samples they allow to quantify, breaking ties by total abundance.

         The keys of @p abundances are stored ordered in @p result, best first.
    */
    template <typename T>
    void orderBest_(const std::map<T, SampleAbundances> & abundances,
                    std::vector<T>& result)
    {
      typedef std::pair<Size, double> PairType;
      std::multimap<PairType, T, std::greater<PairType> > order;
      for (typename std::map<T, SampleAbundances>::const_iterator ab_it =
             abundances.begin(); ab_it != abundances.end(); ++ab_it)
      {
        double total = 0.0;
        for (SampleAbundances::const_iterator samp_it = ab_it->second.begin();
             samp_it != ab_it->second.end(); ++samp_it)
        {
          total += samp_it->second;
        }
        if (total <= 0.0) continue;         // not quantified
        PairType key = std::make_pair(ab_it->second.size(), total);
        order.insert(std::make_pair(key, ab_it->first));
      }
      result.clear();
      for (typename std::multimap<PairType, T, std::greater<PairType> >::
           iterator ord_it = order.begin(); ord_it != order.end(); ++ord_it)
      {
        result.push_back(ord_it->second);
      }
    }



    /**
         @brief Normalize peptide abundances across samples by (multiplicative) scaling to equal medians.
    */
    void normalizePeptides_();

    /**
         @brief Transfer peptide-level quantitative data to protein-level data structures.
         
         This method populates prot_quant_ with peptide abundance and PSM count data.
         
         @param proteins Protein identification information
    */
    void transferPeptideDataToProteins_(const ProteinIdentification& proteins);

    /**
         @brief Select peptides for protein quantification based on filtering criteria.
         
         @param protein_accession The protein accession to select peptides for
         @param top_n Maximum number of peptides to select (0 = no limit)
         @param fix_peptides Whether to use consistent peptides across samples
         @return Vector of selected peptide sequences
    */
    std::vector<String> selectPeptidesForQuantification_(const String& protein_accession,
                                                        Size top_n,
                                                        bool fix_peptides);

    /**
         @brief Aggregate abundances using the specified mathematical method.
         
         @param abundances Vector of abundance values to aggregate
         @param method Aggregation method ("median", "mean", "weighted_mean", "sum")
         @return Aggregated abundance value
    */
    double aggregateAbundances_(const std::vector<double>& abundances,
                               const String& method) const;

    /**
         @brief Calculate protein abundances for a single protein using selected peptides.
         
         @param protein_accession The protein accession
         @param selected_peptides Vector of peptide sequences to use for quantification
         @param aggregate_method Method to aggregate peptide abundances
         @param top_n Maximum number of peptides to use per sample
         @param include_all Whether to include proteins with insufficient peptides
    */
    void calculateProteinAbundances_(const String& protein_accession,
                                    const std::vector<String>& selected_peptides,
                                    const String& aggregate_method,
                                    Size top_n,
                                    bool include_all);

    /**
         @brief Calculate detailed protein abundances at channel level using selected peptides.
         
         @param protein_accession The protein accession
         @param selected_peptides Vector of peptide sequences to use for quantification
         @param aggregate_method Method to aggregate peptide abundances
         @param top_n Maximum number of peptides to use per sample
         @param include_all Whether to include proteins with insufficient peptides
         @param accession_to_leader Map for resolving protein group leaders
    */
    void calculateFileAndChannelLevelProteinAbundances_(const String& protein_accession,
                                           const std::vector<String>& selected_peptides,
                                           const String& aggregate_method,
                                           Size top_n,
                                           bool include_all,
                                           const std::map<String, String>& accession_to_leader);

    /**
         @brief Perform iBAQ normalization on protein abundances.
         
         @param proteins Protein identification information containing sequences
    */
    void performIbaqNormalization_(const ProteinIdentification& proteins);

    /**
         @brief Get the "canonical" protein accession from the list of protein accessions of a peptide.

         @param pep_accessions Protein accessions of a peptide
         @param accession_to_leader Captures information about indistinguishable proteins (maps accession to accession of group leader)

         If there is no information about indistinguishable proteins (from protXML) available, a canonical accession exists only for proteotypic peptides - it's the single accession for the respective peptide.

         Otherwise, a peptide has a canonical accession if it maps only to proteins of one indistinguishable group. In this case, the canonical accession is that of the group leader.

         If there is no canonical accession, the empty string is returned.
    */
    String getAccession_(const std::set<String>& pep_accessions,
                         const std::map<String, String>& accession_to_leader) const;

    /**
         @brief Count the number of identifications (best hits only) of each peptide sequence.

         The peptide hits in @p peptides are sorted by score in the process.
    */
    void countPeptides_(PeptideIdentificationList& peptides);

    /**
         @brief Map (filename, channel) to sample using ExperimentalDesign.
         
         @param filename The base filename (without path/extension)
         @param channel_or_label The channel/label identifier
         @param ed The experimental design containing the mapping information
         @return The sample ID corresponding to the filename and channel
    */
    size_t getSampleIDFromFilenameAndChannel_(const String& filename,
                                           Int channel_or_label,
                                           const ExperimentalDesign& ed) const;

    /// Clear all data when parameters are set
    void updateMembers_() override;

  };   // class

} // namespace
