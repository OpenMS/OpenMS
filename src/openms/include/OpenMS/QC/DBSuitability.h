// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Waschischeck $
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <cfloat>
#include <vector>

#include <boost/regex.hpp>

namespace OpenMS
{
  class ParamXMLFile;
  class PeptideIdentification;
  class PeptideHit;
  class MSExperiment;

  /**
  * @brief This class holds the functionality of calculating the database suitability.
  *
  * To calculate the suitability of a database for a specific mzML for identification search, it
  * is vital to perform a combined deNovo+database identification search. Meaning that the database
  * should be appended with an additional entry derived from concatenated deNovo sequences from said mzML.
  * Currently only Comet search is supported.
  *
  * This class will calculate q-values by itself and will throw an error if any q-value calculation was done beforehand.
  *
  * The algorithm parameters can be set using setParams().
  *
  * Allows for multiple usage of the compute function. The result of each call is stored internally in a vector.
  * Therefore old results will not be overridden by a new call. This vector then can be returned using getResults().
  *
  * This class serves as the library representation of @ref TOPP_DatabaseSuitability
  */
  class OPENMS_DLLAPI DBSuitability:
    public DefaultParamHandler
  {
  public:
    /// struct to store results
    struct OPENMS_DLLAPI SuitabilityData
    {
      /// number of times the top hit is considered to be a deNovo hit
      Size num_top_novo = 0;

      /// number of times the top hit is considered to be a database hit
      Size num_top_db = 0;
      
      /// number of times a deNovo hit scored on top of a database hit
      Size num_interest = 0;

      /// number of times a deNovo hit scored on top of a database hit,
      /// but their score difference was small enough, that it was still counted as a database hit
      Size num_re_ranked = 0;

      /// the cut-off that was used to determine when a score difference was "small enough"
      /// this is normalized by mw
      double cut_off = DBL_MAX;

      /// the suitability of the database used for identification search, calculated with:
      ///               \#db_hits / (\#db_hits + \#deNovo_hit)
      /// can reach from 0 -> the database was not at all suited to 1 -> the perfect database was used
      ///
      /// Preliminary tests have shown that databases of the right organism or close related organisms
      /// score around 0.9 to 0.95, organisms from the same class can still score around 0.8, organisms
      /// from the same phylum score around 0.5 to 0.6 and after that it quickly falls to suitabilities
      /// of 0.15 or even 0.05.
      /// Note that these test were only performed for one mzML and your results might differ.
      double suitability = 0;

      /// the suitability if re-ranking would have been turned off
      /// if re-ranking is actually turned off, this will be the same as the normal suitability
      double suitability_no_rerank = 0;

      /// the suitability after correcting the top deNovo hits, if re-ranking would have been disabled
      double suitability_corr_no_rerank = 0;

      // resets all members to their defaults
      void clear();

      /// apply a correction factor to the already calculated suitability
      /// only works if num_top_db and num_top_novo contain a non-zero value
      void setCorrectionFactor(double factor);

      double getCorrectionFactor() const;

      double getCorrectedNovoHits() const;

      double getCorrectedSuitability() const;

      /**
      * @brief Returns a SuitabilityData object containing the data if re-ranking didn't happen
      *
      * Cases that are re-ranked are already counted. To get the 'no re-ranking' data these cases need to be
      * subtracted from the number of top database hits and added to the number of top deNovo hits.
      *
      * @returns       simulated suitability data where re-ranking didn't happen
      */
      SuitabilityData simulateNoReRanking() const;

    private:
      /// \#IDs with only deNovo search / \#IDs with only database search
      /// used for correcting the number of deNovo hits
      /// worse databases will have less IDs than good databases
      /// this punishes worse databases more than good ones and will result in
      /// a worse suitability
      double corr_factor = -1;

      /// number of top deNovo hits multiplied by the correction factor
      double num_top_novo_corr = 0;

      /// the suitability after correcting the top deNovo hits to impact worse databases more
      ///
      /// The corrected suitability has a more linear behaviour. It basicly translates to the ratio
      /// of the theoretical perfect database the used database corresponds to. (i.e. a corrected
      /// suitability of 0.5 means the used database contains half the proteins of the 'perfect' database)
      double suitability_corr = 0;
    };

    /// Constructor
    /// Settings are initialized with their default values:
    /// no_rerank = false, reranking_cutoff_percentile = 1, FDR = 0.01
    DBSuitability();

    /// Destructor
    ~DBSuitability() override = default;

    /// To test private member functions
    friend class DBSuitability_friend;

    /**
    * @brief Computes suitability of a database used to search a mzML
    *
    * Top deNovo and top database hits from a combined deNovo+database search
    * are counted. The ratio of db hits vs all hits yields the suitability.
    * To re-rank cases, where a de novo peptide scores just higher than
    * the database peptide, a decoy cut-off is calculated. This functionality
    * can be turned off. This will result in an underestimated suitability,
    * but it can solve problems like different search engines or to few decoy hits.
    * 
    * Parameters can be set using the functionality of DefaultParamHandler.
    * Parameters are:
    *           no_rerank                   - re-ranking can be turned off with this (will be set automatically
    *                                         if no cross correlation score is found)
    *           reranking_cutoff_percentile - percentile that determines which cut-off will be returned
    *           FDR                         - q-value that should be filtered for
    *                                         Preliminary tests have shown that database suitability
    *                                         is rather stable across common FDR thresholds from 0 - 5 %
    *           keep_search_files           - temporary files created for and by the internal ID search are kept
    *           disable_correction          - disables corrected suitability calculations
    *           force                       - forces re-ranking to be done even without a cross correlation score,
    *                                         in which case the default main score is used
    *
    * The calculated suitability is then tried to be corrected. For this a correction factor for the number of found top
    * deNovo hits is calculated.
    * This is done by perfoming an additional combined identification search with a smaller sample of the database.
    * It was observed that the number of top deNovo and db hits behave linear according to the sampling ratio of the
    * database. This can be used to extrapolate the number of database hits that would be needed to get a suitability
    * of 1. This number in combination with the maximum number of deNovo hits (found with an identification search
    * where only deNovo is used as a database) can be used to calculate a correction factor like this:
    *                     \#database hits for suitability of 1 / \#maximum deNovo hits
    * This formula can be simplified in a way that the maximum number of deNovo hits isn't needed:
    *                     - (database hits slope) / deNovo hits slope
    * Both of these values can easily be calculated with the original suitability data in conjunction with the one sampled search.
    * 
    * Correcting the number of found top deNovo hits with this factor results in them being more comparable to the top
    * database hits. This in return results in a more linear behaviour of the suitability according to the sampling ratio.
    * The corrected suitability reflects what sampling ratio your database represents regarding to the theoretical 'perfect'
    * database. Or in other words: Your database needs to be (1 - corrected suitability) bigger to get a suitability of 1.
    *
    * Both the original suitability as well as the corrected one are reported in the result.
    *
    * Since q-values need to be calculated the identifications are taken by copy.
    * Since decoys need to be calculated for the fasta input those are taken by copy as well.
    *
    * Result is appended to the result member. This allows for multiple usage.
    *
    * @param pep_ids            vector containing pepIDs with target/decoy annotation coming from a deNovo+database
    *                           identification search without FDR
    *                           (Comet is recommended - to use other search engines either disable reranking or set the '-force' flag)
    *                           vector is modified internally, and is thus copied
    * @param exp                MSExperiment that was searched to produce the identifications
    *                           given in @p pep_ids
    * @param original_fasta     FASTAEntries of the database used for the ID search (without decoys)
    * @param novo_fasta         FASTAEntry derived from deNovo peptides
    * @param search_params      SearchParameters object containing information which adapter
    *                           was used with which settings for the identification search
    *                           that resulted in @p pep_ids
    * @throws                   MissingInformation if no target/decoy annotation is found on @p pep_ids
    * @throws                   MissingInformation if no xcorr is found,
    *                           this happens when another adapter than CometAdapter was used
    * @throws                   Precondition if a q-value is found in @p pep_ids
    */
    void compute(std::vector<PeptideIdentification>&& pep_ids, const MSExperiment& exp, const std::vector<FASTAFile::FASTAEntry>& original_fasta, const std::vector<FASTAFile::FASTAEntry>& novo_fasta, const ProteinIdentification::SearchParameters& search_params);

    /**
    * @brief Returns results calculated by this metric
    *
    * The returned vector contains one DBSuitabilityData object for each time compute was called.
    * Each of these objects contains the suitability information that was extracted from the
    * identifications used for the corresponding call of compute.
    *
    * @returns  DBSuitabilityData objects in a vector
    */
    const std::vector<SuitabilityData>& getResults() const;

  private:
    /// result vector
    std::vector<SuitabilityData> results_;

    /// pattern for finding a decoy string
    const boost::regex decoy_pattern_;

    /**
    * @brief Calculates the xcorr difference between the top two hits marked as decoy
    *
    * Searches for the top two decoys hits and returns their score difference.
    * By default the xcorr from Comet is used. If no xcorr can be found and the 'force' flag is set
    * the main score from the peptide hit is used, else an error is thrown.
    *
    * If there aren't two decoys, DBL_MAX is returned.
    *
    * @param pep_id     pepID from where the decoy difference will be calculated
    * @returns          xcorr difference
    * @throws           MissingInformation if no target/decoy annotation is found
    * @throws           MissingInformation if no xcorr is found
    */
    double getDecoyDiff_(const PeptideIdentification& pep_id) const;

    /**
    * @brief Calculates a xcorr cut-off based on decoy hits
    *
    * Decoy differences of all N pepIDs are calculated. The (1-reranking_cutoff_percentile)*N highest
    * one is returned.
    * It is assumed that this difference accounts for 'reranking_cutoff_percentile' of the re-ranking cases.
    *
    * @param pep_ids                      vector containing the pepIDs
    * @param reranking_cutoff_percentile  percentile that determines which cut-off will be returned
    * @returns                            xcorr cut-off
    * @throws                             IllegalArgument if reranking_cutoff_percentile isn't in range [0,1]
    * @throws                             IllegalArgument if reranking_cutoff_percentile is too low for a decoy cut-off to be calculated
    * @throws                             MissingInformation if no more than 20 % of the peptide IDs have two decoys in their top ten peptide hits
    */
    double getDecoyCutOff_(const std::vector<PeptideIdentification>& pep_ids, double reranking_cutoff_percentile) const;

    /**
    * @brief Tests if a PeptideHit is considered a deNovo hit
    *
    * To test this the function looks into the protein accessions.
    * If only the deNovo protein is found, 'true' is returned.
    * If at least one database protein is found, 'false' is returned.
    *
    * This function also uses boost::regex_search to make sure the deNovo accession doesn't contain a decoy string.
    * This is needed for 'target+decoy' hits.
    *
    * @param hit      PepHit in question
    * @returns        true/false
    */
    bool isNovoHit_(const PeptideHit& hit) const;

    /**
    * @brief Tests if a PeptideHit has a score better than the given threshold
    *
    * @param hit                    PepHit in question
    * @param threshold              threshold to check against
    * @param higher_score_better    true/false depending if a higher or a lower score is better
    * @returns                      true/false
    */
    bool checkScoreBetterThanThreshold_(const PeptideHit& hit, double threshold, bool higher_score_better) const;

    /**
    * @brief Looks through meta values of SearchParameters to find out which search adapter was used
    *
    * Checks for the following adapters:
    * CometAdapter, MSGFPlusAdapter, MSFraggerAdapter, MyriMatchAdapter, OMSSAAdapter and XTandemAdapter
    *
    * @param meta_values   SearchParameters object, since the adapters write their parameters here
    * @returns             A pair containing the name of the adapter and the parameters used to run it
    * @throws              MissingInformation if none of the adapters above is found in the meta values
    */
    std::pair<String, Param> extractSearchAdapterInfoFromMetaValues_(const ProteinIdentification::SearchParameters& meta_values) const;

    /**
    * @brief Writes parameters into a given file
    *
    * @param parameters    parameters to write
    * @param filename      name of the file where the parameters should be written to
    * @throws              UnableToCreateFile if filename isn't writable
    */
    void writeIniFile_(const Param& parameters, const String& filename) const;

    /**
    * @brief Executes the workflow from search adapter, followed by PeptideIndexer and finishes with FDR
    *
    * Which adapter should run with which parameters can be controlled.
    * Make sure the search adapter you wish to use is built on your system and the executable is on your PATH variable.
    *
    * Indexing and FDR are always done the same way.
    *
    * The inputs are stored in temporary files to execute the Adapter.
    * (MSExperiment -> .mzML, vector<FASTAEntry> -> .fasta, Param -> .INI)
    *
    * @param exp            MSExperiment that will be searched
    * @param fasta_data     represents the database that should be used to search
    * @param adapter_name   name of the adapter to search with
    * @param parameters     parameters for the adapter
    * @returns              peptide identifications with annotated q-values
    * @throws               MissingInformation if no adapter name is given
    * @throws               InvalidParameter if a not supported adapter name is given
    * @throws               InternalToolError if any error occures while running the adapter
    * @throws               InternalToolError if any error occures while running PeptideIndexer functionalities
    * @throws               InvalidParameter if the needed FDR parameters are not found
    */
    std::vector<PeptideIdentification> runIdentificationSearch_(const MSExperiment& exp, const std::vector<FASTAFile::FASTAEntry>& fasta_data, const String& adapter_name, Param& parameters) const;

    /**
    * @brief Creates a subsampled fasta with the given subsampling rate
    *
    * The subsampling is based on the number of amino acides and not on the number of fasta entries.
    *
    * @param fasta_data         fasta of which the subsampling should be done
    * @param subsampling_rate   subsampling rate to be used [0,1]
    * @returns                  fasta entries with total number of AA = original number of AA * subsampling_rate
    * @throws                   IllegalArgument if subsampling rate is not between 0 and 1
    */
    std::vector<FASTAFile::FASTAEntry> getSubsampledFasta_(const std::vector<FASTAFile::FASTAEntry>& fasta_data, double subsampling_rate) const;

    /**
    * @brief Calculates all suitability data from a combined deNovo+database search
    *
    * Counts top database and top deNovo hits.
    *
    * Calculates a decoy score cut-off to compare high scoring deNovo hits with lower scoring database hits.
    * If the score difference is smaller than the cut-off the database hit is counted and the deNovo hit ignored.
    *
    * Suitability is calculated: # database hits / # all hits
    *
    * @param pep_ids    peptide identifications coming from the combined search, each peptide identification should be sorted
    * @param data       SuitabilityData object where the result should be written into
    * @throws           MissingInformation if no target/decoy annotation is found on @p pep_ids
    * @throws           MissingInformation if no xcorr is found,
    *                   this happens when another adapter than CometAdapter was used
    */
    void calculateSuitability_(const std::vector<PeptideIdentification>& pep_ids, SuitabilityData& data) const;

    /**
    * @brief Calculates and appends decoys to a given vector of FASTAEntry
    *
    * Each sequence is digested with Trypsin. The resulting peptides are reversed and appended to one another.
    * This results in the decoy sequences.
    * The identifier is given a 'DECOY_' prefix.
    *
    * @param fasta     reference to fasta vector where the decoys are needed
    */
    void appendDecoys_(std::vector<FASTAFile::FASTAEntry>& fasta) const;

    /**
    * @brief Returns the cross correlation score normalized by MW (if existing), else if the 'force' flag is set the current main score is returned
    *
    * @param pep_hit    PeptideHit of which the score is needed
    * @returns          cross correlation score normalized by MW or current score
    * @throws           MissingInformation if no xcorr is found and 'force' flag isn't set
    */
    double extractScore_(const PeptideHit& pep_hit) const;

    /**
    * @brief Calculates the correction factor from two suitability calculations
    *
    * Two suitability calculations need to be done for this. One with the original data and one with data from a search with a sampled database.
    * The number of db hits and deNovo hits behaves linear. The two searches can than be used to calculate the
    * corresponding linear functions.
    * The factor is calculated with the negative ratio of the db slope and the deNovo slope.
    *
    * @param data            suitability data from the original search
    * @param data_sampled    vector of suitability data from the sampled search(s)
    * @param sampling_rate   the sampling rate used for sampled db [0,1)
    * @returns               correction factor
    */
    double calculateCorrectionFactor_(const SuitabilityData& data, const SuitabilityData& data_sampled, double sampling_rate) const;

    /**
    * @brief Determines the number of unique proteins found in the protein accessions of PeptideIdentifications
    *
    * @param peps               vector of PeptideIdentifications
    * @param number_of_hits     the number of hits to search in (if this is bigger than the actual number of hits all hits are looked at)
    * @returns                  number of unique protein accessions
    * @throws                   MissingInformation if no target/decoy annotation is found on @p peps
    */
    UInt numberOfUniqueProteins_(const std::vector<PeptideIdentification>& peps, UInt number_of_hits = 1) const;

    /**
    * @brief Finds the SuitabilityData object with the median number of de novo hits
    *
    *  If the median isn't distinct (e.g. two entries could be considered median) the upper one is chosen.
    *
    * @param data     vector of SuitabilityData objects
    * @returns        index to object with median number of de novo hits
    */
    Size getIndexWithMedianNovoHits_(const std::vector<SuitabilityData>& data) const;

    /**
    * @brief Extracts the worst score that still passes a FDR (q-value) threshold
    *
    * This can be used to 'convert' a FDR threshold to a threshold for the desired score (score and FDR need to be dependent)
    *
    * @param pep_ids              vector of PeptideIdentifications
    * @param FDR                  FDR threshold, hits with a worse q-value score aren't looked at
    * @param score_name           name of the score to search for
    *                             The score name doesn't need to be the exact metavalue name, but a metavalue key should contain it.
    *                             i.e. "e-value" as metavalue "e-value_score"
    * @param higher_score_better  true/false depending if a higher or lower score (@p score_name) is better
    * @returns                    the worst score that is still in the FDR threshold
    *
    * @throws                     IllegalArgument if @p score_name isn't found in the metavalues
    * @throws                     Precondition if main score of @p pep_ids isn't 'q-value'
    */
    double getScoreMatchingFDR_(const std::vector<PeptideIdentification>& pep_ids, double FDR, const String& score_name, bool higher_score_better) const;
  };

  // friend class to test private member functions
  class DBSuitability_friend
  {
  public:
    DBSuitability_friend() = default;

    ~DBSuitability_friend() = default;

    std::vector<FASTAFile::FASTAEntry> getSubsampledFasta(const std::vector<FASTAFile::FASTAEntry>& fasta_data, double subsampling_rate)
    {
      return suit_.getSubsampledFasta_(fasta_data, subsampling_rate);
    }

    void appendDecoys(std::vector<FASTAFile::FASTAEntry>& fasta)
    {
      suit_.appendDecoys_(fasta);
    }

    double calculateCorrectionFactor(const DBSuitability::SuitabilityData& data, const DBSuitability::SuitabilityData& data_sampled, double sampling_rate)
    {
      return suit_.calculateCorrectionFactor_(data, data_sampled, sampling_rate);
    }

    UInt numberOfUniqueProteins(const std::vector<PeptideIdentification>& peps, UInt number_of_hits = 1)
    {
      return suit_.numberOfUniqueProteins_(peps, number_of_hits);
    }

    Size getIndexWithMedianNovoHits(const std::vector<DBSuitability::SuitabilityData>& data)
    {
      return suit_.getIndexWithMedianNovoHits_(data);
    }

    double getScoreMatchingFDR(const std::vector<PeptideIdentification>& pep_ids, double FDR, String score_name, bool higher_score_better)
    {
      return suit_.getScoreMatchingFDR_(pep_ids, FDR, score_name, higher_score_better);
    }

    /* Not tested:
      getDecoyDiff_, getDecoyCutOff_, isNovoHit_, checkScoreBetterThanThreshold_
      Reason: These functions are essential to the normal suitability calculation and if something would not work, the test for 'compute' would fail.

      extractSearchAdapterInfoFromMetaValues_, writeIniFile_, extractScore_
      Reason: These functions are very straightforeward.

      runIdentificationSearch_
      Reason: This function simulates a whole workflow and testing it would be to complicated.
    */

  private:
    DBSuitability suit_;
  };
}

