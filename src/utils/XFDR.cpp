// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Lukas Zimmermann, Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/XQuestResultXMLFile.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/MATH/STATISTICS/Histogram.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XQuestResultXMLHandler.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <string>
#include <cmath>

#include <cassert>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_XFDR XFDR

    @brief Calculates false discovery rate estimates on crosslink identifications.

    This tool calculates and FDR estimate for crosslink identifications, which are produced by OpenPepXL.
    The method employed currently is identical to the target-decoy approach used by xProphet (Walzthoeni et al., 2012).
    Consequently, this tool can also consume xquest.xml files (produced either by OpenPepXL or xQuest). The tool supports
    output in the idXML and mzIdentML formats.

    @experimental This tool is work in progress and usage and input requirements might change.

    <center>
        <table>
            <tr>
                <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
                <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ XFDR \f$ \longrightarrow \f$</td>
                <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
            </tr>
            <tr>
                <td VALIGN="middle" ALIGN = "center" ROWSPAN=1>  OpenPepXL </td>
                <td VALIGN="middle" ALIGN = "center" ROWSPAN=1>  OpenPepXLLF </td>
                <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> - </td>
            </tr>
            <tr>
            </tr>
        </table>
    </center>

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_XFDR.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_XFDR.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPXFDR final :
public TOPPBase
{
public:

  static const String param_in;  // Parameter for the input file
  static const String param_in_type;
  static const String param_out_idXML;
  static const String param_out_mzid;
  static const String param_out_xquest;
  static const String param_decoy_string; // full prefix for decoy proteins
  static const String param_minborder;  // minborder  # filter for minimum precursor mass error (ppm)
  static const String param_maxborder;  // maxborder  # filter for maximum precursor mass error (ppm)
  static const String param_mindeltas;  // mindeltas  0.95 # filter for delta score, 0 is no filter, minimum delta score required, hits are rejected if larger or equal
  static const String param_minionsmatched; // minionsmatched 0 # Filter for minimum matched ions per peptide
  static const String param_uniquexl; // calculate statistics based on unique IDs
  static const String param_no_qvalues; // Do not transform to qvalues
  static const String param_minscore; // minscore 0 # minimum ld-score to be considered
  static const String param_binsize; // bin size for cumulative histograms

  // Number of ranks used
  static const UInt n_rank;

  // Constants related to particular crosslink classes
  static const String crosslink_class_intradecoys; // intradecoys
  static const String crosslink_class_fulldecoysintralinks; // fulldecoysintralinks
  static const String crosslink_class_interdecoys; // interdecoys
  static const String crosslink_class_fulldecoysinterlinks; // fulldecoysinterlinks
  static const String crosslink_class_intralinks; // intralinks
  static const String crosslink_class_interlinks; // interlinks
  static const String crosslink_class_monolinks;  // monolinks
  static const String crosslink_class_monodecoys; // monodecoys
  static const String crosslink_class_decoys; // decoys
  static const String crosslink_class_targets; // targets
  static const String crosslink_class_hybriddecoysintralinks; // hybriddecoysintralinks
  static const String crosslink_class_hybriddecoysinterlinks; // hybriddecoysintralinks

  // Meta values used to identify cross-links
  static const String crosslink_type;
  static const String crosslink_rank;
  static const String target_decoy;

  TOPPXFDR() :
    TOPPBase("XFDR", "Calculates false discovery rate estimates on crosslink identifications", false),
    min_score_(0),
    max_score_(0)
  {
  }

protected:

  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() final
  {
    StringList formats = ListUtils::create<String>("xml,idXML,mzid,xquest.xml");

    // File input
    registerInputFile_(TOPPXFDR::param_in, "<file>", "", "Crosslink Identifications in either xquest.xml, idXML, or mzIdentML format (as produced by OpenPepXL)", false);
    setValidFormats_(TOPPXFDR::param_in, formats);

    // File input type (if omitted, guessed from the file extension)
    registerStringOption_(TOPPXFDR::param_in_type, "<in_type>", "", "Type of input file provided with -in. If omitted, the file type is guessed from the file extension.", false, false);
    setValidStrings_(TOPPXFDR::param_in_type, formats);

    // idXML output
    registerOutputFile_(TOPPXFDR::param_out_idXML, "<idXML_file>", "", "Output as idXML file", false, false);
    setValidFormats_(TOPPXFDR::param_out_idXML, ListUtils::create<String>("idXML"));

    // mzIdentML output
    registerOutputFile_(TOPPXFDR::param_out_mzid, "<mzIdentML_file>", "", "Output as mzIdentML file", false, false);
    setValidFormats_(TOPPXFDR::param_out_mzid, ListUtils::create<String>("mzid"));

    // xquest.xml output
    registerOutputFile_(TOPPXFDR::param_out_xquest, "<xQuestXML_file>", "", "Output as xquest.xml file", false, false);
    setValidFormats_(TOPPXFDR::param_out_xquest, ListUtils::create<String>("xquest.xml"));

    // decoy prefix
    registerStringOption_(TOPPXFDR::param_decoy_string, "<string>", "DECOY_", "Prefix of decoy protein ids. The correspondig target protein id should be retrievable by deleting this prefix.", false);

    // Minborder
    registerDoubleOption_(TOPPXFDR::param_minborder, "<minborder>", -50, "Filter for minimum precursor mass error (ppm) before FDR estimation. Values outside of the tolerance window of the original search will effectively disable this filter.", false);

    // Maxborder
    registerDoubleOption_(TOPPXFDR::param_maxborder, "<maxborder>", 50, "Filter for maximum precursor mass error (ppm) before FDR estimation. Values outside of the tolerance window of the original search will effectively disable this filter.", false);

    // Mindeltas
    registerDoubleOption_(TOPPXFDR::param_mindeltas, "<mindeltas>", 0, "Filter for delta score, 0 disables the filter. Minimum delta score required, hits are rejected if larger or equal. The delta score is a ratio of the score of a hit and the score of the next best hit to the same spectrum, so the value range is between 0 and 1 with 1.0 meaning the scores are equal and 0.5 meaning the next best score is half as high as the current one.", false);
    setMinFloat_(TOPPXFDR::param_mindeltas, 0.0);
    setMaxFloat_(TOPPXFDR::param_mindeltas, 1.0);

    // Minionsmatched
    registerIntOption_(TOPPXFDR::param_minionsmatched, "<minionsmatched>", 0, "Filter for minimum matched ions per peptide.", false);
    setMinInt_(TOPPXFDR::param_minionsmatched, 0);

    // Uniquexl
    registerFlag_(TOPPXFDR::param_uniquexl, "Calculate statistics based only on unique IDs. For a set of IDs from equal candidates (same pair of peptides, modifications and cross-linked positions), only the highest scoring hit will be considered. By default the score distribution will be estimated using all 1st ranked candidates.");

    // Qtransform
    registerFlag_(TOPPXFDR::param_no_qvalues, "Do not transform simple FDR to q-values");

    // Minscore
    registerDoubleOption_(TOPPXFDR::param_minscore, "<minscore>", -10, "Minimum score to be considered for FDR calculation. A number lower than the lowest score will effectively disable this filter.", false);

    // Cumulative Histograms bin size
    registerDoubleOption_(TOPPXFDR::param_binsize, "<binsize>", 0.0001, "Bin size for the cumulative histograms for score distributions. Should be about the same size as the smallest expected difference between scores. Smaller numbers will make XFDR more robust, but much slower. Negative numbers are not allowed. Should only be changed if the range of the main score changes or another score than the OpenPepXL score is used.", false, true);
    setMinFloat_(TOPPXFDR::param_binsize, 1e-15);
  }

    /**
   * @brief Prepares vector of PeptideIdentification such that it can be processed downstream.
   * The encompassed steps are:
   *  * Set min_score_ and max_score_ encountered in the data
   *  * Ensure that crosslink_type and crosslink_rank are available in the PeptideIdentification
   *  * Define peptide_identification as decoy if at least one of the peptide_hits is decoy
   *    (this semantic is also used by the XQuestResultXMLHandler)
   *  * Define the crosslink as either inter/or intraprotein
   *  * Set the identifier of the Peptide Identification if there is only one protein identification
   *
   */
  void initDataStructures()
  {
    const String prot_identifier = prot_id_.getIdentifier();
    String decoy_string = getStringOption_(TOPPXFDR::param_decoy_string);

    // if the metaValue exists in search_params and the default value for XFDR was not changed, use the one in search_params
    ProteinIdentification::SearchParameters search_params = prot_id_.getSearchParameters();
    if (search_params.metaValueExists("decoy_string") && decoy_string == "DECOY_")
    {
      decoy_string = search_params.getMetaValue("decoy_string");
    }

    // Preprocess all peptide identifications and construct derived data structures necessary for XFDR
    for (Size i = 0; i < all_pep_ids_.size(); ++i)
    {
      PeptideIdentification &pep_id = all_pep_ids_[i];
      if (pep_id.getHits().size() < 1)
      {
        continue;
      }

      pep_id.setIdentifier(prot_identifier);

      vector< PeptideHit > &pep_hits = pep_id.getHits();

      for (PeptideHit& ph : pep_hits)
      {
        // Set the minScore and MaxScore attribute depending on the input data
        // const double score = getCrosslinkScore(pep_id);
        const double score = ph.getScore();

        // Set score boundaries
        if (score < this->min_score_)
        {
          this->min_score_ = std::floor(score);
        }
        if (score > this->max_score_)
        {
          this->max_score_ = std::ceil(score);
        }
        assert(this->min_score_ <= this->max_score_);

        // figure out if crosslink is inter- or intra protein
        // for cases with multiple proteins, count as true, if any one possible combination of proteins fits the criteria
        // so both can be true at the same time (or false for mono-links)
        setIntraProtein(ph, false);
        setInterProtein(ph, false);

        if (ph.metaValueExists(Constants::OPENPEPXL_XL_TYPE) && ph.getMetaValue(Constants::OPENPEPXL_XL_TYPE) == "cross-link")
        {
          // String alpha_prots_string;
          StringList alpha_prots;
          const std::vector<PeptideEvidence> pevs_alpha = ph.getPeptideEvidences();
          for (PeptideEvidence pev : pevs_alpha)
          {
            // alpha_prots_string = alpha_prots_string + "," + pev->getProteinAccession();
            alpha_prots.push_back(pev.getProteinAccession());
          }

          // StringList alpha_prots = ListUtils::create<String>(ph.getAccession());
          StringList beta_prots = ListUtils::create<String>(ph.getMetaValue(Constants::OPENPEPXL_BETA_ACCESSIONS).toString());


          for (String& alpha_prot : alpha_prots)
          {
            for (String& beta_prot : beta_prots)
            {
              if (isSameProtein(alpha_prot, beta_prot, decoy_string))
              {
                // cout << "TEST Protein Combo: " << alpha_prot << " | " << beta_prot << " == INTRA" << endl;
                setIntraProtein(ph, true);
              }
              else
              {
                // cout << "TEST Protein Combo: " << alpha_prot << " | " << beta_prot << " == INTER" << endl;
                setInterProtein(ph, true);
              }
            }
          }
        }

        String id = getId(ph);
        ph.setMetaValue("OpenPepXL:id", id);
        // candidates with the same ID will also have the same types
        if (this->cross_link_classes_.find(id) == this->cross_link_classes_.end())
        {
          assignTypes(ph, this->cross_link_classes_[id]);
        }
      }
    }
    if (arg_uniquex)
    {
      findTopUniqueHits();
    }
  }


  /**
   * @brief Inspects PeptideIdentification pep_id and assigns all cross-link types that this identification belongs to
   * @param pep_id Peptide ID to be assigned.
   * @param types Result vector containing the names of the crosslink classes
   */
  static void assignTypes(PeptideHit &ph, StringList &types)
  {
    types.clear();
    bool xl_is_decoy = ph.getMetaValue(Constants::TARGET_DECOY) == "decoy";

    // target or decoy
    if (xl_is_decoy)
    {
      types.push_back(TOPPXFDR::crosslink_class_decoys);
    }
    else
    {
      types.push_back(TOPPXFDR::crosslink_class_targets);
    }

    // intralinks
    if (ph.getMetaValue("XFDR:is_intraprotein").toBool() && (!xl_is_decoy))
    {
      types.push_back(TOPPXFDR::crosslink_class_intralinks);
    }

    // intradecoys
    if (ph.getMetaValue("XFDR:is_intraprotein").toBool() && xl_is_decoy)
    {
      types.push_back(TOPPXFDR::crosslink_class_intradecoys);
    }

    // interlinks
    if (ph.getMetaValue("XFDR:is_interprotein").toBool() && (!xl_is_decoy))
    {
      types.push_back(TOPPXFDR::crosslink_class_interlinks);
    }

    // interdecoys
    if (ph.getMetaValue("XFDR:is_interprotein").toBool() && xl_is_decoy)
    {
      types.push_back(TOPPXFDR::crosslink_class_interdecoys);
    }

    assert(ph.metaValueExists(Constants::OPENPEPXL_XL_TYPE));
    String current_crosslink_type = ph.getMetaValue(Constants::OPENPEPXL_XL_TYPE);

    // monolinks
    if ( (!xl_is_decoy) && (current_crosslink_type == "mono-link"
        ||  current_crosslink_type == "loop-link"))
    {
      types.push_back(TOPPXFDR::crosslink_class_monolinks);
    }

    // monodecoys
    if ( xl_is_decoy && (current_crosslink_type == "mono-link"
        ||  current_crosslink_type == "loop-link"))
    {
      types.push_back(TOPPXFDR::crosslink_class_monodecoys);
    }

    if (current_crosslink_type == "cross-link")
    {
      const bool alpha_is_decoy = ph.getMetaValue(Constants::OPENPEPXL_TARGET_DECOY_ALPHA).toString() == "decoy";
      const bool beta_is_decoy = ph.getMetaValue(Constants::OPENPEPXL_TARGET_DECOY_BETA).toString() == "decoy";

      // fulldecoysintralinks
      if (ph.getMetaValue("XFDR:is_intraprotein").toBool() && alpha_is_decoy && beta_is_decoy)
      {
        types.push_back(TOPPXFDR::crosslink_class_fulldecoysintralinks);
      }

      // fulldecoysinterlinks
      if (ph.getMetaValue("XFDR:is_interprotein").toBool() && alpha_is_decoy && beta_is_decoy)
      {
        types.push_back(TOPPXFDR::crosslink_class_fulldecoysinterlinks);
      }

      // hybriddecoysintralinks
      if (ph.getMetaValue("XFDR:is_intraprotein").toBool()
          && (( (!alpha_is_decoy) &&   beta_is_decoy)
          ||     (alpha_is_decoy  && (!beta_is_decoy))))
      {
        types.push_back(TOPPXFDR::crosslink_class_hybriddecoysintralinks);
      }

      // hybriddecoysinterlinks
      if (ph.getMetaValue("XFDR:is_interprotein").toBool()
          && (( (!alpha_is_decoy) &&   beta_is_decoy)
          ||     (alpha_is_decoy  && (!beta_is_decoy))))
      {
        types.push_back(TOPPXFDR::crosslink_class_hybriddecoysinterlinks);
      }
    }
  }

  /** Target counting as performed by the xProphet software package
   *
   * @brief xprophet  method for target hits counting as implemented in xProphet
   * @param cum_histograms Cumulative score distributions
   */
  void fdr_xprophet(std::map< String, Math::Histogram<> > & cum_histograms,
                    const String  & targetclass, const String & decoyclass, const String & fulldecoyclass,
                    vector< double > & fdr, bool mono)
  {
    // Determine whether targetclass, decoyclass, and fulldecoyclass are present in the histogram map
    bool targetclass_present = cum_histograms.find(targetclass) != cum_histograms.end();
    bool decoyclass_present = cum_histograms.find(decoyclass) != cum_histograms.end();
    bool fulldecoyclass_present = cum_histograms.find(fulldecoyclass) != cum_histograms.end();

    for (double current_score = this->min_score_ +  (arg_binsize/2);
        current_score <= this->max_score_ - (arg_binsize/2);
        current_score += arg_binsize)
    {
      double estimated_n_decoys = decoyclass_present ? cum_histograms[decoyclass].binValue(current_score) : 0;
      if ( ! mono)
      {
        estimated_n_decoys -= 2 * ( fulldecoyclass_present ? cum_histograms[fulldecoyclass].binValue(current_score) : 0);
      }
      double n_targets = targetclass_present ? cum_histograms[targetclass].binValue(current_score) : 0;
      fdr.push_back(n_targets > 0 ? estimated_n_decoys / (n_targets) : 0);
    }
  }

  /**
   * @brief Calculates the qFDR values for the provided FDR values, assuming that the FDRs are sorted by score in the input vector
   * @param fdr Vector with FDR values which should be used for qFDR calculation
   * @param qfdr Result qFDR values
   */
  void calc_qfdr(const vector< double > &fdr, vector< double > &qfdr)
  {
    qfdr.resize(fdr.size());
    for (Int i = fdr.size() - 1; i >= 0; --i)
    {
      double current_fdr = fdr[i];
      double smallest_fdr = current_fdr;
      for (Int j = i; j >= 0; j--)
      {
        double fdr_to_check = fdr[j];
        if (fdr_to_check < smallest_fdr)
        {
          smallest_fdr = fdr_to_check;
        }
      }
      qfdr[i] = smallest_fdr < current_fdr ? smallest_fdr : current_fdr;
    }
  }

  void findTopUniqueHits()
  {
    for (PeptideIdentification& pep_id : this->all_pep_ids_)
    {
      for (PeptideHit& ph : pep_id.getHits())
      {
        String id = ph.getMetaValue("OpenPepXL:id");
        auto uid_it = std::find(this->unique_ids_.begin(), this->unique_ids_.end(), id);
        // if an ID for this candidate already exists, check if the new score is higher than the last
        if (uid_it != this->unique_ids_.end())
        {
          int index = std::distance(this->unique_ids_.begin(), uid_it);
          if (this->unique_id_scores_[index] < ph.getScore())
          {
            this->unique_id_scores_[index] = ph.getScore();
          }
        }
        else
        {
          this->unique_ids_.push_back(id);
          this->unique_id_scores_.push_back(ph.getScore());
        }
      }
    }
  }

  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char **) final
  {
    // Tool Arguments
    loadArguments();
    validateArguments();
    writeArgumentsLog();

    writeLog_("Reading input file...");
    // Input File loading, initializes all_pep_ids_ vector
    ExitCodes load_result = loadInputFile();
    if (load_result != EXECUTION_OK)
    {
      logFatal("Loading of input file has failed");
      return load_result;
    }

    writeLog_("Initializing data structures...");
    // Initialize and validate data structures that are derived from the main peptide identification vector 'all_pep)ids'
    initDataStructures();

    // Maps the cross link class to the encountered scores
    map<String, vector<double>> scores;
    UInt num_flagged(0);

    writeLog_("Collecting scores for each class...");
    // Loop through the peptides, apply filter, and assign cross-link types
    // for (const Size &idx : this->rank_one_pep_idx)
    for (PeptideIdentification& pep_id : this->all_pep_ids_)
    {
      if (pep_id.getHits().size() < 1)
      {
        continue;
      }

      // usually the 1st ranked hit should be the first in the list
      // but we make sure here and use the first hit with rank 1 that we find
      int rank_one_hit_index(0);
      for (Size i = 0; i < pep_id.getHits().size(); ++i)
      {
        PeptideHit& ph = pep_id.getHits()[i];
        if (int(ph.getMetaValue(Constants::OPENPEPXL_XL_RANK)) == 1)
        {
          rank_one_hit_index = i;
          break;
        }
      }
      PeptideHit& ph = pep_id.getHits()[rank_one_hit_index];


      // if after the search above we don't have a rank 1 hit, skip this pep_id
      if (int(ph.getMetaValue(Constants::OPENPEPXL_XL_RANK)) != 1)
      {
        continue;
      }

      // Attributes of peptide identification that can be used for filtering
      // const double delta_score = calculateDeltaScore(pep_id, rank_one_hit_index, rank_two_hit_index);
      const double delta_score = ph.getMetaValue(Constants::DELTA_SCORE);
      const double score = ph.getScore();

      double error_rel(0);
      if (ph.metaValueExists(Constants::PRECURSOR_ERROR_PPM_USERPARAM))
      {
        error_rel = ph.getMetaValue(Constants::PRECURSOR_ERROR_PPM_USERPARAM);
      }

      const Size min_ions_matched = getMinIonsMatched(ph);
      const String id = ph.getMetaValue("OpenPepXL:id");

      // Only consider IDs which fullfill all filter criteria specified by the user
      if (   (arg_minborder <= error_rel)   // minborder fullfilled
          && (arg_maxborder >= error_rel)   // maxborder fullfilled
          && (arg_mindeltas == 0  || delta_score < arg_mindeltas)
          && (min_ions_matched  >= (Size)arg_minionsmatched)
          &&  score >= arg_minscore
         )
      {
        // check for the unique ID criterion
        if (arg_uniquex)
        {
          auto uid_it = std::find(this->unique_ids_.begin(), this->unique_ids_.end(), id);
          if (uid_it != this->unique_ids_.end())
          {
            int index = std::distance(this->unique_ids_.begin(), uid_it);
            if (this->unique_id_scores_[index] != ph.getScore())
            {
              // this is not the highest scoring ID for this candidate
              continue;
            }
          }
        }
        num_flagged++;

        ph.setMetaValue("XFDR:used_for_FDR", 1);

        for (const String &cross_link_class : this->cross_link_classes_[id])
        {
          scores[cross_link_class].push_back(score);
        }
      }
    }
    writeLog_(this->toolName_() + " has used " + num_flagged + " hits to calculate the FDR");

    // Log number of scores within each class
    writeLog_("Number of Scores for each class:");

    for (const auto &score : scores)
    {
      writeLog_(score.first + ": " + score.second.size());
    }

    // Generate Histograms of the scores for each class
    // Use cumulative histograms to count the number of scores above consecutive thresholds
    std::map< String, Math::Histogram<> >  cum_histograms;
    for (const auto &class_scores: scores)
    {
      vector< double > current_scores = class_scores.second;

      Math::Histogram<> histogram(this->min_score_, this->max_score_, arg_binsize);
      Math::Histogram<>::getCumulativeHistogram(current_scores.begin(), current_scores.end(), true, true, histogram);
      cum_histograms[class_scores.first] = histogram;
    }

    writeLog_("Calculating Score Distributions...");
    // Calculate FDR for interlinks
    vector< double > fdr_interlinks;
    this->fdr_xprophet(cum_histograms, TOPPXFDR::crosslink_class_interlinks, TOPPXFDR::crosslink_class_interdecoys, TOPPXFDR::crosslink_class_fulldecoysinterlinks, fdr_interlinks, false);

    // Calculate FDR for intralinks
    vector< double > fdr_intralinks;
    this->fdr_xprophet(cum_histograms, TOPPXFDR::crosslink_class_intralinks, TOPPXFDR::crosslink_class_intradecoys, TOPPXFDR::crosslink_class_fulldecoysintralinks, fdr_intralinks, false);

    // Calculate FDR for monolinks and looplinks
    vector< double > fdr_monolinks;
    this->fdr_xprophet(cum_histograms, TOPPXFDR::crosslink_class_monolinks, TOPPXFDR::crosslink_class_monodecoys, "", fdr_monolinks, true);

    // Determine whether qTransform should be performed (and consequently the score type)
    bool arg_no_qvalues = getFlag_(TOPPXFDR::param_no_qvalues);
    String score_type = arg_no_qvalues ? "FDR" : "q-value";

    if ( ! arg_no_qvalues)
    {
      writeLog_("Performing qFDR transformation...");

      vector< double > qfdr_interlinks;
      this->calc_qfdr(fdr_interlinks, qfdr_interlinks);

      vector< double > qfdr_intralinks;
      this->calc_qfdr(fdr_intralinks, qfdr_intralinks);

      vector< double > qfdr_monolinks;
      this->calc_qfdr(fdr_monolinks, qfdr_monolinks);

      fdr_interlinks = qfdr_interlinks;
      fdr_intralinks = qfdr_intralinks;
      fdr_monolinks = qfdr_monolinks;
    }

    writeLog_("Assigning FDRs...");
    // Assign FDR values to all identifications
    for (PeptideIdentification &pep_id : all_pep_ids_)
    {
      for (PeptideHit& ph : pep_id.getHits())
      {
        if ( ! ph.metaValueExists("XFDR:used_for_FDR"))
        {
          ph.setMetaValue("XFDR:used_for_FDR", 0);
        }
        double score = ph.getScore();

        StringList crosslink_types;
        assignTypes(ph, crosslink_types);

        ph.setMetaValue("XFDR:fdr_type", score_type);

        // Assign FDR value as meta value and also set as score
        bool assigned = false;
        double fdr = 1;
        for (StringList::const_iterator crosslink_types_it = crosslink_types.begin();
            crosslink_types_it != crosslink_types.end(); ++crosslink_types_it)
        {
          String current_crosslink_type = *crosslink_types_it;
          Size idx = std::floor((score - this->min_score_) / arg_binsize);
          if (   current_crosslink_type == TOPPXFDR::crosslink_class_fulldecoysinterlinks
              || current_crosslink_type == TOPPXFDR::crosslink_class_hybriddecoysinterlinks
              || current_crosslink_type == TOPPXFDR::crosslink_class_interdecoys
              || current_crosslink_type == TOPPXFDR::crosslink_class_interlinks)
          {
            fdr = fdr_interlinks[idx];
            assigned = true;
            break;
          }
          else if (   current_crosslink_type == TOPPXFDR::crosslink_class_fulldecoysintralinks
              || current_crosslink_type == TOPPXFDR::crosslink_class_hybriddecoysintralinks
              || current_crosslink_type == TOPPXFDR::crosslink_class_intradecoys
              || current_crosslink_type == TOPPXFDR::crosslink_class_intralinks)
          {
            fdr = fdr_intralinks[idx];
            assigned = true;
            break;
          }
          else if (   current_crosslink_type == TOPPXFDR::crosslink_class_monodecoys
              || current_crosslink_type == TOPPXFDR::crosslink_class_monolinks)
          {
            fdr = fdr_monolinks[idx];
            assigned = true;
            break;
          }
        }
        if ( assigned)
        {
          ph.setMetaValue("XFDR:FDR", fdr);
        }
        else
        {
          OPENMS_LOG_WARN << "WARNING: A Crosslink could not be identified as either interlink, intralink, or monolink, so no FDR will be available for it." << endl;
        }
      }
    }

    writeLog_("Writing output...");
    // write idXML
    if (! arg_out_idXML.empty())
    {
      IdXMLFile().store( arg_out_idXML, all_prot_ids_, all_pep_ids_);
    }

    // write mzid file
    if (! arg_out_mzid.empty())
    {
      MzIdentMLFile().store( arg_out_mzid, all_prot_ids_, all_pep_ids_);
    }

    // write xquest.xml file
    if (! arg_out_xquest.empty())
    {
      XQuestResultXMLFile().store(arg_out_xquest, all_prot_ids_, all_pep_ids_);
    }
    return EXECUTION_OK;
  }

private:

  // Score range for this of the tool
  Int min_score_;
  Int max_score_;

  // unique top hits
  vector<String> unique_ids_;
  vector<double> unique_id_scores_;

  // Data structures
  // Vector of peptideIdentifications and indizes for rank one hits
  vector<PeptideIdentification> all_pep_ids_;

  // maps index of peptide id all_pep_ids_ to vector of cross link class
  map<String, vector<String>> cross_link_classes_;

  // Protein Identification
  vector<ProteinIdentification> all_prot_ids_;
  ProteinIdentification prot_id_;

  // Whether input comes from xQuest
  bool is_xquest_input_;

  // Program arguments
  String arg_out_idXML;
  String arg_out_mzid;
  String arg_out_xquest;
  String arg_in;
  double arg_mindeltas;
  double arg_minborder;
  double arg_maxborder;
  Int arg_minionsmatched;
  double arg_minscore;
  bool arg_uniquex;
  double arg_binsize;


  void logFatal(const String &message) const
  {
    OPENMS_LOG_ERROR << "FATAL: " << message << " Terminating now!" << std::endl;
  }

  void loadArguments()
  {
    arg_out_idXML = getStringOption_(TOPPXFDR::param_out_idXML);
    arg_out_mzid = getStringOption_(TOPPXFDR::param_out_mzid);
    arg_out_xquest = getStringOption_(TOPPXFDR::param_out_xquest);
    arg_in = getStringOption_(TOPPXFDR::param_in);
    arg_mindeltas = getDoubleOption_(TOPPXFDR::param_mindeltas);
    arg_minborder = getDoubleOption_(TOPPXFDR::param_minborder);
    arg_maxborder = getDoubleOption_(TOPPXFDR::param_maxborder);
    arg_minionsmatched = getIntOption_(TOPPXFDR::param_minionsmatched);
    arg_minscore = getDoubleOption_(TOPPXFDR::param_minscore);
    arg_uniquex = getFlag_(TOPPXFDR::param_uniquexl);
    arg_binsize = getDoubleOption_(TOPPXFDR::param_binsize);
  }

  void writeArgumentsLog() const
  {
    //-------------------------------------------------------------
    // Printing parameters to log
    //-------------------------------------------------------------
    writeLog_(arg_minborder != -1 ? "Lower bound for precursor mass error for FDR calculation is " + String(arg_minborder) + " ppm"
                                  : "No lower bound for precursor mass error for FDR calculation");
    writeLog_(arg_maxborder != -1 ? "Upper bound for precursor mass error for FDR calculation is " + String(arg_maxborder) + " ppm"
                                  : "No upper bound for precursor mass error for FDR calculation");
    writeLog_(arg_mindeltas != 0  ? "Filtering of hits by a deltascore of " + String(arg_mindeltas) + " is used."
                                  : "No filtering of hits by deltascore");
    writeLog_(arg_minionsmatched > 0 ? "Filtering of hits by minimum ions matched: " + String(arg_minionsmatched) + " is used"
                                     : "No filtering of hits by minimum ions matched.");
    writeLog_(arg_minscore > 0 ? "Filtering of hits by minimum score of " + String(arg_minscore) + " is used."
                               : "No filtering of hits by minimum score.");
    writeLog_(arg_uniquex ? "Error model is generated based on unique cross-links."
                          : "Error model is generated based on redundant cross-links.");
    writeLog_("Bin size for cumulative histograms is " + String(arg_binsize));
  }

  ExitCodes validateArguments() const
  {
    if (this->arg_out_idXML.empty() && this->arg_out_mzid.empty() && this->arg_out_xquest.empty())
    {
      logFatal(
              "No output file specified. You must at least specify one output with -"
              + String(TOPPXFDR::param_out_idXML)
              + " or -" + String(TOPPXFDR::param_out_mzid)
              + " or -" + String(TOPPXFDR::param_out_xquest)
              + " or -" + String(TOPPXFDR::param_out_xquest)
      );
      return ILLEGAL_PARAMETERS;
    }

    if (arg_in.empty())
    {
      logFatal("Input file is empty");
      return ILLEGAL_PARAMETERS;
    }
    if (arg_minborder > arg_maxborder)
    {
      logFatal("Minborder cannot be larger than Maxboder!");
      return ILLEGAL_PARAMETERS;
    }
    return EXECUTION_OK;
  }


  /**
  * Loads the input file depending on the type. Returns 0 if the loading of the input was successful, error
  * code otherwise
  * @return 0 if the loading of the input was successful, error code otherwise
  */
  ExitCodes loadInputFile()
  {
    //------------------------------------------------------------
    // Determine type of input file
    //-------------------------------------------------------------
    const String arg_in_type = getStringOption_(TOPPXFDR::param_in_type);
    const FileTypes::Type in_type = arg_in_type.empty() ?
                                      FileHandler::getType(this->arg_in) : FileTypes::nameToType(arg_in_type);

    this->is_xquest_input_ = false;
    if (in_type == FileTypes::XQUESTXML)
    {
      this->is_xquest_input_ = true;
      XQuestResultXMLFile xquest_file;
      xquest_file.load(arg_in, all_pep_ids_, all_prot_ids_);

     writeLog_("\nTotal number of hits in xQuest input: " + String(xquest_file.getNumberOfHits()));
    }
    else if (in_type == FileTypes::MZIDENTML)
    {
       MzIdentMLFile().load(arg_in, all_prot_ids_, all_pep_ids_);
     }
     else if (in_type == FileTypes::IDXML)
     {
       // TODO I doubt that the filters are supported for IDXML Input
       IdXMLFile().load(arg_in, all_prot_ids_, all_pep_ids_);
     }
     else
     {
       logFatal("Input file type not recognized.");
       return ILLEGAL_PARAMETERS;
     }
     const Size n_pep_ids = all_pep_ids_.size();
     const Size n_prot_ids = all_prot_ids_.size();

     writeLog_("Number of Peptide IDs in input file: " + String(n_pep_ids));
     writeLog_("Number of Protein IDs in input file: " + String(n_prot_ids));

     // Terminate if no hits could be found
     if (n_pep_ids == 0)
     {
       logFatal("Input file does not contain any identifications.");
       return INPUT_FILE_EMPTY;
     }

     // Terminate if do not exactly encounter one protein id
     if (n_prot_ids != 1)
     {
       logFatal("There is not exactly one protein identification in the input file. This is unsupported!");
       return INPUT_FILE_CORRUPT;
     }
     this->prot_id_ = all_prot_ids_[0];

     return EXECUTION_OK;
  }

  String getId(const PeptideHit& ph) const
  {
    if (is_xquest_input_)
    {
      return ph.getMetaValue("OpenPepXL:id").toString();
    }

    if (ph.getMetaValue(Constants::OPENPEPXL_XL_TYPE) == "cross-link")
    {
      return   ph.getSequence().toUnmodifiedString()
               + "-" + AASequence::fromString(ph.getMetaValue(Constants::OPENPEPXL_BETA_SEQUENCE)).toUnmodifiedString()
               + "-a" + String(ph.getMetaValue(Constants::OPENPEPXL_XL_POS1))
               + "-b" + String(ph.getMetaValue(Constants::OPENPEPXL_XL_POS2));

    }
    else if (ph.getMetaValue(Constants::OPENPEPXL_XL_TYPE) == "loop-link")
    {
      return   ph.getSequence().toUnmodifiedString()
               + "-a" + String(ph.getMetaValue(Constants::OPENPEPXL_XL_POS1))
               + "-b" + String(ph.getMetaValue(Constants::OPENPEPXL_XL_POS2));

    }
    else if (ph.metaValueExists(Constants::OPENPEPXL_XL_MASS))
    {
      return   ph.getSequence().toUnmodifiedString()
               + "-" + String(ph.getMetaValue(Constants::OPENPEPXL_XL_POS1))
               + "-" + String(ph.getMetaValue(Constants::OPENPEPXL_XL_MASS));
    }
    else
    {
      return   ph.getSequence().toUnmodifiedString()
               + "-" + String(ph.getMetaValue(Constants::OPENPEPXL_XL_POS1));
    }
  }

  static Size getMinIonsMatched(const PeptideHit& ph)
  {
    Size alpha_ions = Size(ph.getMetaValue("matched_linear_alpha")) + Size(ph.getMetaValue("matched_xlink_alpha"));
    Size beta_ions = Size(ph.getMetaValue("matched_linear_beta")) + Size(ph.getMetaValue("matched_xlink_beta"));
    return std::min(alpha_ions, beta_ions);
  }


  static void setIntraProtein(PeptideHit& ph, const bool value)
  {
    ph.setMetaValue("XFDR:is_intraprotein", DataValue(value ? "true" : "false"));
  }

  static void setInterProtein(PeptideHit& ph, const bool value)
  {
    ph.setMetaValue("XFDR:is_interprotein", DataValue(value ? "true" : "false"));
  }

  /**
   *  @brief Determines whether the Petide Evidences belong to the same protein, modulo decoy
   */
  static bool isSameProtein(
          String prot1,
          String prot2,
          const String &decoy_string)
  {
    prot1.substitute(decoy_string, "");
    prot2.substitute(decoy_string, "");
    assert( ! prot1.hasSubstring(decoy_string));
    assert( ! prot2.hasSubstring(decoy_string));
    return prot1 == prot2;
  }
};


const String TOPPXFDR::param_in = "in";
const String TOPPXFDR::param_in_type = "in_type";
const String TOPPXFDR::param_out_idXML = "out_idXML";
const String TOPPXFDR::param_out_mzid = "out_mzIdentML";
const String TOPPXFDR::param_out_xquest = "out_xquest";
const String TOPPXFDR::param_decoy_string = "decoy_string";
const String TOPPXFDR::param_minborder = "minborder";
const String TOPPXFDR::param_maxborder = "maxborder";
const String TOPPXFDR::param_mindeltas = "mindeltas";
const String TOPPXFDR::param_minionsmatched = "minionsmatched";
const String TOPPXFDR::param_uniquexl = "uniquexl";
const String TOPPXFDR::param_no_qvalues = "no_qvalues";
const String TOPPXFDR::param_minscore = "minscore";
const String TOPPXFDR::param_binsize = "binsize";

const String TOPPXFDR::crosslink_class_intradecoys = "intradecoys";
const String TOPPXFDR::crosslink_class_fulldecoysintralinks = "fulldecoysintralinks";
const String TOPPXFDR::crosslink_class_interdecoys = "interdecoys";
const String TOPPXFDR::crosslink_class_fulldecoysinterlinks = "fulldecoysinterlinks";
const String TOPPXFDR::crosslink_class_monodecoys = "monodecoys";
const String TOPPXFDR::crosslink_class_intralinks = "intralinks";
const String TOPPXFDR::crosslink_class_interlinks = "interlinks";
const String TOPPXFDR::crosslink_class_monolinks  = "monolinks";
const String TOPPXFDR::crosslink_class_decoys = "decoys";
const String TOPPXFDR::crosslink_class_targets = "targets";
const String TOPPXFDR::crosslink_class_hybriddecoysintralinks = "hybriddecoysintralinks";
const String TOPPXFDR::crosslink_class_hybriddecoysinterlinks = "hybriddecoysinterlinks";

// the actual main function needed to create an executable
int main(int argc, const char ** argv)
{
  TOPPXFDR tool;
  return tool.main(argc, argv);
}
/// @endcond
