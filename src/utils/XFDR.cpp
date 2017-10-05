// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Lukas Zimmermann $
// $Authors: Lukas Zimmermann $
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

#include <boost/iterator/counting_iterator.hpp>

#include <string>
#include <math.h>

#include <assert.h>

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


class TOPPXFDR :
    public TOPPBase
{
  public:

    static const String param_in;  // Parameter for the input file
    static const String param_in_type;
    static const String param_out_idXML;
    static const String param_out_mzid;
    static const String param_minborder;  // minborder  # filter for minimum precursor mass error (ppm)
    static const String param_maxborder;  // maxborder  # filter for maximum precursor mass error (ppm)
    static const String param_mindeltas;  // mindeltas  0.95 # filter for delta score, 0 is no filter, minimum delta score required, hits are rejected if larger or equal
    static const String param_minionsmatched; // minionsmatched 0 # Filter for minimum matched ions per peptide
    static const String param_uniquexl; // calculate statistics based on unique IDs
    static const String param_no_qvalues; // Do not transform to qvalues
    static const String param_minscore; // minscore 0 # minimum ld-score to be considered

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
    static const String crosslink_class_hybriddecoysintralinks; // hybriddecoysintralinks
    static const String crosslink_class_hybriddecoysinterlinks; // hybriddecoysintralinks

    // Score range for calculating the FPs
    static const double fpnum_score_step;

    // Meta values used to identify cross-links
    static const String crosslink_type;
    static const String crosslink_rank;
    static const String target_decoy;

    TOPPXFDR() :
      TOPPBase("XFDR", "Calculates false discovery rate estimates on crosslink identifications", false)
    {
    }

  protected:

    // this function will be used to register the tool parameters
    // it gets automatically called on tool execution
    void registerOptionsAndFlags_()
    {
      StringList formats = ListUtils::create<String>("xml,idXML,mzid");

      // File input
      registerInputFile_(TOPPXFDR::param_in, "<file>", "", "Crosslink Identifications in either xquest.xml, idXML, or mzIdentML format (as produced by OpenPepXL)", false);
      setValidFormats_(TOPPXFDR::param_in, formats);

      // File input type (if omitted, guessed from the file extension)
      registerStringOption_(TOPPXFDR::param_in_type, "<in_type>", "", "Type of input file provided with -in", false, false);
      setValidStrings_(TOPPXFDR::param_in_type, formats);

      // idXML output
      registerOutputFile_(TOPPXFDR::param_out_idXML, "<idXML_file>", "", "Output as idXML file", false, false);
      setValidFormats_(TOPPXFDR::param_out_idXML, ListUtils::create<String>("idXML"));

      // mzIdentML output
      registerOutputFile_(TOPPXFDR::param_out_mzid, "<mzIdentML_file>", "", "Output as mzIdentML file", false, false);
      setValidFormats_(TOPPXFDR::param_out_mzid, ListUtils::create<String>("mzid"));

      // Minborder
      registerIntOption_(TOPPXFDR::param_minborder, "<minborder>", -1, "Filter for minimum precursor mass error (ppm). Can only be applied for xquest.xml input.", false);

      // Maxborder
      registerIntOption_(TOPPXFDR::param_maxborder, "<maxborder>", -1, "Filter for maximum precursor mass error (ppm). Can only be applied for xquest.xml input.", false);

      // Mindeltas
      registerDoubleOption_(TOPPXFDR::param_mindeltas, "<mindeltas>", 0, "Filter for delta score, 0 is no filter. Minimum delta score required, hits are rejected if larger or equal. Can only be applied for xquest.xml input.", false);

      // Minionsmatched
      registerIntOption_(TOPPXFDR::param_minionsmatched, "<minionsmatched>", 0, "Filter for minimum matched ions per peptide. Can only be applied for xquest.xml input.", false);

      // Uniquexl
      registerFlag_(TOPPXFDR::param_uniquexl, "Calculate statistics based on unique IDs. Can only be applied for xquest.xml input.");

      // Qtransform
      registerFlag_(TOPPXFDR::param_no_qvalues, "Do not transform simple FDR to q-values");

      // Minscore
      registerIntOption_(TOPPXFDR::param_minscore, "<minscore>", 0, "Minimum score to be considered for FDR calculation", false);
    }

    /**
     * @brief Used to define how PeptideIdentification are to be sorted based on some meta value (usually some score)
     * "greater_than" ensures that the sorting will be in descending order
     */
    struct greater_than_by_key
    {
      // Vector containing all crosslink peptide identifications of the input file
      const vector< PeptideIdentification > & all_ids;

      //Indices of the rank one elements within the all_ids vector
      const vector< UInt >  & rank_one_ids;

      inline bool operator()(const UInt & index1, const UInt & index2)
      {
        return (  TOPPXFDR::getCrosslinkScore(all_ids[rank_one_ids[index1]])
                > TOPPXFDR::getCrosslinkScore(all_ids[rank_one_ids[index2]]));
      }
    };

    /**
     * @brief Returns the score of a crosslink peptide identification.
     * @param pep_id Which peptide identification the score should be taken from
     * @return crosslink score of that peptide identification
     */
    static double getCrosslinkScore(const PeptideIdentification & pep_id)
    {
      const vector< PeptideHit > & pep_hits = pep_id.getHits();

#ifndef NDEBUG
      UInt n_hits = pep_hits.size();
      assert(n_hits == 1 || n_hits == 2);
      if (n_hits == 2)
      {
        assert(pep_hits[0].getScore() == pep_hits[1].getScore());
      }
#endif
      return pep_hits[0].getScore();
    }

    /**
     * @brief Moves the meta value denotes by @p meta_value from the first peptide hit to the peptide identification
     * @param pep_id The PeptideIdentification whose meta value will be adjusted
     * @param meta_value The MetaValue that will be moved from the PeptideHit to the PeptideIdentification @p pep_id
     * @return Whether the MetaValue @p meta_value could be found at all
     */
    static bool moveToPeptideIdentification(PeptideIdentification & pep_id, const String & meta_value)
    {
      const vector< PeptideHit > & pep_hits = pep_id.getHits();
      assert(pep_hits.size() == 1 || pep_hits.size() == 2);

      if ( ! pep_id.metaValueExists(meta_value))
      {
        if ( ! pep_hits[0].metaValueExists(meta_value))
        {
          // Meta value has not been found neither in the peptide identification nor in the first peptide hit
          return false;
        }
        pep_id.setMetaValue(meta_value, pep_hits[0].getMetaValue(meta_value));
        assert(pep_id.metaValueExists(meta_value));
      }
      return true;
    }

    /**
     * @brief Prepares vector of PeptideIdentification such that it can be processed downstream.
     * The encompassed steps are:
     *  * Set min_score and max_score encountered in the data
     *  * Ensure that crosslink_type and crosslink_rank are available in the PeptideIdentification
     *  * Define peptide_identification as decoy if at least one of the peptide_hits is decoy
     *    (this semantic is also used by the XQuestResultXMLHandler)
     *  * Define the crosslink as either inter/or intraprotein
     *  * Set the identifier of the Peptide Identification if there is only one protein identification
     *
     * @param pep_ids vector of PeptideIdentification to be prepared
     */
    bool prepareInput(vector< PeptideIdentification > & pep_ids, const vector< ProteinIdentification > & prot_ids)
    {
      bool set_prot_id = prot_ids.size() == 1;

      // Preprocess all peptide identifications
      for (vector< PeptideIdentification >::iterator pep_ids_it = pep_ids.begin();
           pep_ids_it != pep_ids.end(); ++pep_ids_it)
      {
        PeptideIdentification & pep_id = *pep_ids_it;

        if (set_prot_id)
        {
          pep_id.setIdentifier(prot_ids[0].getIdentifier());
        }

        // Set the minScore and MaxScore attribute depending on the input data
        double score = getCrosslinkScore(pep_id);
        if (score < this->min_score)
        {
          this->min_score = std::floor(score);
        }
        if (score > this->max_score)
        {
          this->max_score = std::ceil(score);
        }

        // Fetch the PeptideHits
        const vector< PeptideHit > & pep_hits = pep_id.getHits();
        UInt n_hits = pep_hits.size();
        assert(n_hits == 1 || n_hits == 2);

        // Pull crosslink_type and crosslink_rank to the peptide_identification if necessary
        // Return false if one of the methods fails

        if ( ! TOPPXFDR::moveToPeptideIdentification(pep_id, TOPPXFDR::crosslink_type) ||
             ! TOPPXFDR::moveToPeptideIdentification(pep_id, TOPPXFDR::crosslink_rank) )
        {
          return false;
        }

        // Set peptide identification to target or decoy
        if (    pep_hits[0].getMetaValue(TOPPXFDR::target_decoy).toString() == "decoy"
                || (n_hits == 2 && pep_hits[1].getMetaValue(TOPPXFDR::target_decoy).toString() == "decoy"))
        {
          pep_id.setMetaValue(TOPPXFDR::target_decoy, DataValue("decoy"));
        }
        else
        {
          pep_id.setMetaValue(TOPPXFDR::target_decoy, DataValue("target"));
        }

        // figure out if crosslink is inter- or intra protein
        if (n_hits == 2)
        {
          vector< PeptideEvidence > alpha_ev = pep_hits[0].getPeptideEvidences();
          vector< PeptideEvidence > beta_ev = pep_hits[1].getPeptideEvidences();

          for (vector< PeptideEvidence >::const_iterator alpha_ev_it  = alpha_ev.begin();
               alpha_ev_it != alpha_ev.end(); ++alpha_ev_it)
          {
            for (vector< PeptideEvidence >::const_iterator beta_ev_it = beta_ev.begin();
                 beta_ev_it != beta_ev.end(); ++beta_ev_it)
            {
              String alpha_prot = alpha_ev_it->getProteinAccession();
              String beta_prot = beta_ev_it->getProteinAccession();

              alpha_prot.substitute("reverse_", "");
              alpha_prot.substitute(Internal::XQuestResultXMLHandler::decoy_string, "");
              beta_prot.substitute("reverse_", "");
              beta_prot.substitute(Internal::XQuestResultXMLHandler::decoy_string, "");
              pep_id.setMetaValue( alpha_prot == beta_prot ? "OpenXQuest:is_intraprotein" : "OpenXQuest:is_interprotein" , DataValue());
            }
          }
        }
      }
      return true;
    }


    /**
     * @brief Inspects PeptideIdentification pep_id and assigns all cross-link types that this identification belongs to
     * @param pep_id Peptide ID to be assigned.
     * @param types Result vector containing the names of the crosslink classes
     */
    inline static void assignTypes(PeptideIdentification & pep_id, StringList & types)
    {
      types.clear();
      bool pep_is_decoy = pep_id.getMetaValue("target_decoy").toString() == "decoy";

      // Intradecoys
      if (pep_id.metaValueExists("OpenXQuest:is_intraprotein") && pep_is_decoy)
      {
        types.push_back(TOPPXFDR::crosslink_class_intradecoys);
      }

      // Decoys
      if (pep_is_decoy)
      {
        types.push_back(TOPPXFDR::crosslink_class_decoys);
      }

      // intralinks
      if (pep_id.metaValueExists("OpenXQuest:is_intraprotein") && ! pep_is_decoy)
      {
        types.push_back(TOPPXFDR::crosslink_class_intralinks);
      }

      // interdecoys
      if (pep_id.metaValueExists("OpenXQuest:is_interprotein") && pep_is_decoy)
      {
        types.push_back(TOPPXFDR::crosslink_class_interdecoys);
      }

      // interlinks
      if (pep_id.metaValueExists("OpenXQuest:is_interprotein") && ! pep_is_decoy)
      {
        types.push_back(TOPPXFDR::crosslink_class_interlinks);
      }

      assert(pep_id.metaValueExists(TOPPXFDR::crosslink_type));
      String crosslink_type = pep_id.getMetaValue(TOPPXFDR::crosslink_type);

      // monolinks
      if ( ! pep_is_decoy && (crosslink_type == "mono-link"
                              ||  crosslink_type == "loop-link"))
      {
        types.push_back(TOPPXFDR::crosslink_class_monolinks);
      }

      // monodecoys
      if ( pep_is_decoy && (crosslink_type == "mono-link"
                            ||  crosslink_type == "loop-link"))
      {
        types.push_back(TOPPXFDR::crosslink_class_monodecoys);
      }
      const vector< PeptideHit > & pep_hits = pep_id.getHits();
      if (pep_hits.size() == 2)
      {
        PeptideHit alpha = pep_hits[0];
        PeptideHit beta = pep_hits[1];
        
        bool alpha_is_decoy = alpha.getMetaValue("target_decoy").toString() == "decoy";
        bool beta_is_decoy = beta.getMetaValue("target_decoy").toString() == "decoy";

        // fulldecoysintralinks
        if (   pep_id.metaValueExists("OpenXQuest:is_intraprotein")
               && alpha_is_decoy
               && beta_is_decoy)
        {
          types.push_back(TOPPXFDR::crosslink_class_fulldecoysintralinks);
        }

        // fulldecoysinterlinks
        if (   pep_id.metaValueExists("OpenXQuest:is_interprotein")
               && alpha_is_decoy
               && beta_is_decoy)
        {
          types.push_back(TOPPXFDR::crosslink_class_fulldecoysinterlinks);

        }

        // hybriddecoysintralinks
        if (       pep_id.metaValueExists("OpenXQuest:is_intraprotein")
                   && (( ! alpha_is_decoy
                         &&     beta_is_decoy)
                       ||     (alpha_is_decoy
                               &&   ! beta_is_decoy)))
        {
          types.push_back(TOPPXFDR::crosslink_class_hybriddecoysintralinks);
        }

        // hybriddecoysinterlinks
        if (       pep_id.metaValueExists("OpenXQuest:is_interprotein")
                   && (( ! alpha_is_decoy
                         &&     beta_is_decoy)
                       ||     (alpha_is_decoy
                               &&   ! beta_is_decoy)))
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

      for (double current_score = this->min_score +  (TOPPXFDR::fpnum_score_step/2) ;
           current_score <= this->max_score - (TOPPXFDR::fpnum_score_step/2);
           current_score += TOPPXFDR::fpnum_score_step)
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
    void calc_qfdr(const vector< double > & fdr, vector< double > & qfdr)
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

    // the main_ function is called after all parameters are read
    ExitCodes main_(int, const char **)
    {
      //-------------------------------------------------------------
      // Initialize instance variables
      //-------------------------------------------------------------
      bool is_xquest_input = false;
      this->min_score = 0;
      this->max_score = 0;

      //----------------------------------------------------------------
      // parsing parameters, terminate if invalid values are encountered
      //----------------------------------------------------------------

      // Check whether at least one output file has been specified
      String arg_out_idXML = getStringOption_(TOPPXFDR::param_out_idXML);
      String arg_out_mzid = getStringOption_(TOPPXFDR::param_out_mzid);

      if (arg_out_idXML.empty() && arg_out_mzid.empty())
      {
        LOG_ERROR << "ERROR: No output file specified. You must at least specify one output with -"
                  <<  TOPPXFDR::param_out_idXML << " or -" << TOPPXFDR::param_out_mzid << ". Terminating." << endl;
        return ILLEGAL_PARAMETERS;
      }

      double arg_mindeltas = getDoubleOption_(TOPPXFDR::param_mindeltas);
      if (arg_mindeltas < 0)
      {
        LOG_ERROR << "ERROR: Negative values for parameter 'mindeltas' are not allowed." << endl;
        return ILLEGAL_PARAMETERS;
      }
      if (arg_mindeltas > 1)
      {
        LOG_ERROR << "ERROR: Values larger than 1 for parameter 'mindeltas' are not allowed." << endl;
        return ILLEGAL_PARAMETERS;
      }
      // mindelta if 0 disables this filter (according to the documentation of xProphet)
      bool mindelta_filter_disabled = (arg_mindeltas == 0);

      Int arg_minborder = getIntOption_(TOPPXFDR::param_minborder);
      Int arg_maxborder = getIntOption_(TOPPXFDR::param_maxborder);

      if (arg_minborder > arg_maxborder)
      {
        LOG_ERROR << "ERROR: Minborder cannot be larger than Maxborder. Terminating" << endl;
        return ILLEGAL_PARAMETERS;
      }
      Int arg_minionsmatched = getIntOption_(TOPPXFDR::param_minionsmatched);

      if (arg_minionsmatched < 0)
      {
        LOG_ERROR << "ERROR: Minionsmatched cannot be negative. Terminating." << endl;
        return ILLEGAL_PARAMETERS;
      }

      Int arg_minscore = getIntOption_(TOPPXFDR::param_minscore);
      bool arg_uniquex = getFlag_(TOPPXFDR::param_uniquexl);

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

      //-------------------------------------------------------------
      // Determine type of input file
      //-------------------------------------------------------------
      String arg_in = getStringOption_(TOPPXFDR::param_in);
      if (arg_in.empty())
      {
        LOG_ERROR << "ERROR: Input file is empty. Terminating." << endl;
        return ILLEGAL_PARAMETERS;
      }
      String arg_in_type = getStringOption_(TOPPXFDR::param_in_type);
      FileTypes::Type in_type = arg_in_type.empty() ? FileHandler::getType(arg_in)
                                                    : FileTypes::nameToType(arg_in_type);

      //-------------------------------------------------------------
      // Declare important variables
      //-------------------------------------------------------------
      Size pep_id_index = TOPPXFDR::n_rank - 1;

      // Variable is only used if the input is xQuest XML file
      Size n_spectra = 0;

      // Main data structures
      vector < PeptideIdentification > all_ids;
      vector < ProteinIdentification > prot_ids;
      vector < UInt > rank_one_ids; // Stores the indizes of the rank one hits within all_ids
      vector < vector < PeptideIdentification > > spectra;

      //-------------------------------------------------------------
      // Parse the input file
      //-------------------------------------------------------------

      if (in_type == FileTypes::XML)
      {
        is_xquest_input = true;
        
        XQuestResultXMLFile xquest_result_file;
        xquest_result_file.load(arg_in, spectra, prot_ids, 1, true);

        // currently, cross-link identifications are stored within one ProteinIdentification
        assert(prot_ids.size() == 1);
        n_spectra = spectra.size();
        writeLog_("Total number of spectra: " + String(n_spectra)  +
                  "\nTotal number of hits: " + String(xquest_result_file.getNumberOfHits()));

        Size rank_counter = 0;
        for (vector < vector < PeptideIdentification > >::const_iterator spectra_it = spectra.begin();
             spectra_it != spectra.end(); ++spectra_it)
        {
          vector< PeptideIdentification > spectrum = *spectra_it;
          for (vector< PeptideIdentification >::const_iterator spectrum_it = spectrum.begin(); spectrum_it != spectrum.end(); ++spectrum_it)
          {
            PeptideIdentification pep_id = *spectrum_it;
            all_ids.push_back(pep_id);

            if ( static_cast<int>(pep_id.getMetaValue(TOPPXFDR::crosslink_rank)) == 1)
            {
              rank_one_ids.push_back(rank_counter);
            }
            rank_counter++;
          }
        }
        writeLog_("Number of IDs in input file: " + String(all_ids.size()));

        // Terminate if no hits could be foud
        if (all_ids.size() == 0)
        {
          LOG_ERROR << "ERROR: Input file does not contain any identifications. Terminating." << endl;
          return INPUT_FILE_EMPTY;
        }

#ifndef NDEBUG
        for (vector< PeptideIdentification >::const_iterator all_ids_it = all_ids.begin();
             all_ids_it != all_ids.end(); ++ all_ids_it)
        {
          const PeptideIdentification & pep_id = *all_ids_it;
          const vector< PeptideHit > & pep_hits = pep_id.getHits();
          UInt n_hits = pep_hits.size();
          String pep_id_target_decoy = pep_id.getMetaValue("target_decoy").toString();
          String alpha_target_decoy = pep_hits[0].getMetaValue("target_decoy").toString();

          if (n_hits == 1)
          {
            assert(    (alpha_target_decoy == "decoy" && pep_id_target_decoy == "decoy")
                    || (alpha_target_decoy == "target" && pep_id_target_decoy == "target"));

          }
          else if (n_hits == 2)
          {
            String beta_target_decoy = pep_hits[1].getMetaValue("target_decoy").toString();
            if (alpha_target_decoy == "target" && beta_target_decoy == "target" )
            {
              assert(pep_id_target_decoy == "target");
            }
            else
            {
              assert(pep_id_target_decoy == "decoy");
            }
          }
        }
#endif
      }
      else if (in_type == FileTypes::MZIDENTML)
      {
        // Prevent filter options for this input (currently not supported)
        if (arg_uniquex || arg_minborder != -1 || arg_maxborder != -1 || arg_minionsmatched != 0 || arg_mindeltas != 0)
        {
          LOG_ERROR << "ERROR: The filters uniquexl min/maxborder, minionsmatched, and mindeltas are not supported for idXML. Terminating." << endl;
          return ILLEGAL_PARAMETERS;
        }

        MzIdentMLFile().load(arg_in, prot_ids, all_ids);
        writeLog_("Number of IDs in input file: " + String(all_ids.size()));

        // Terminate if no hits could be foud
        if (all_ids.size() == 0)
        {
          LOG_ERROR << "ERROR: Input file does not contain any identifications. Terminating." << endl;
          return INPUT_FILE_EMPTY;
        }

        Size rank_counter = 0;
        for (vector< PeptideIdentification >::const_iterator all_ids_it = all_ids.begin();
             all_ids_it != all_ids.end(); ++all_ids_it)
        {
          PeptideIdentification pep_id = *all_ids_it;

          if (static_cast<UInt>(pep_id.getMetaValue(TOPPXFDR::crosslink_rank)) == 1)
          {
            rank_one_ids.push_back(rank_counter);
          }
          rank_counter++;
        }
      }
      else if (in_type == FileTypes::IDXML)
      {
        // Prevent filter options for this input (currently not supported)
        if (arg_uniquex || arg_minborder != -1 || arg_maxborder != -1 || arg_minionsmatched != 0 || arg_mindeltas != 0)
        {
          LOG_ERROR << "ERROR: The filters uniquexl min/maxborder, minionsmatched, and mindeltas are not supported for idXML. Terminating." << endl;
          return ILLEGAL_PARAMETERS;
        }

        IdXMLFile().load(arg_in, prot_ids, all_ids);
        writeLog_("Number of IDs in input file: " + String(all_ids.size()));

        // Terminate if no hits could be foud
        if (all_ids.size() == 0)
        {
          LOG_ERROR << "ERROR: Input file does not contain any identifications. Terminating." << endl;
          return INPUT_FILE_EMPTY;
        }

        Size rank_counter = 0;
        for (vector< PeptideIdentification >::const_iterator all_ids_it = all_ids.begin();
             all_ids_it != all_ids.end(); ++all_ids_it)
        {
          PeptideIdentification pep_id = *all_ids_it;
          
          if ( static_cast<UInt>(pep_id.getMetaValue(TOPPXFDR::crosslink_rank)) == 1)
          {
            rank_one_ids.push_back(rank_counter);
          }
          rank_counter++;
        }
      }

      // Prepare input data
      if ( ! this->prepareInput(all_ids, prot_ids))
      {
        LOG_ERROR << "ERROR: Input data could not be prepared. Terminating." << endl;
        return INPUT_FILE_CORRUPT;
      }

      // Number of peptide identifications that need to be considered
      Size n_ids = rank_one_ids.size();

      //-------------------------------------------------------------
      // Calculate the delta score for each hit
      // Calculate n_min_ions_matched
      // Currently only for xQuest input files
      //-------------------------------------------------------------
      // The score is calculated for each hit h on the set of all hits of the spectrum that encompasses
      std::vector< std::vector< double > > delta_scores;
      std::vector< Size > n_min_ions_matched;

      // For xQuest input,calculate delta scores and min_ions_matched
      if (is_xquest_input)
      {
        writeLog_("Input is a xQuest result file. Compute the delta scores and the number of matched ions");
        delta_scores.reserve(n_spectra);
        n_min_ions_matched.resize(n_spectra);

        for (Size i = 0; i < n_spectra; ++i)
        {
          Size n_hits = spectra[i].size();
          //delta_scores[i] = new std::vector<double>(n_hits);
          vector<double> current(n_hits);

          assert(n_hits > 0); // because we initially do not load 'empty' spectra
          // calculate n_min_ions_matched
          PeptideIdentification * pep_id1 = &spectra[i][0];
          assert( static_cast<int>(pep_id1->getMetaValue(TOPPXFDR::crosslink_rank)) == 1); // because hits are sorted according to their rank within the spectrum
          const vector<PeptideHit> & pep_hits = pep_id1->getHits();

          if ( pep_id1->getMetaValue(TOPPXFDR::crosslink_type) == "cross-link")
          {
            n_min_ions_matched[i] = std::min( static_cast<int>(pep_hits[0].getMetaValue("OpenXQuest:num_of_matched_ions")),
                static_cast<int>(pep_hits[1].getMetaValue("OpenXQuest:num_of_matched_ions")));
          }
          else
          {
            n_min_ions_matched[i] = static_cast<int>(pep_hits[0].getMetaValue("OpenXQuest:num_of_matched_ions"));
          }
          // Calculate delta score
          if (n_hits > 1)
          {
            for (Size j = 0; j < n_hits - 1; ++j)
            {
              pep_id1 = &spectra[i][j];
              for (Size k = 1; j+k < n_hits; ++k )
              {
                PeptideIdentification * pep_id2 = &spectra[i][j+k];
                if (pep_id1->getMetaValue("OpenXQuest:structure") != pep_id2->getMetaValue("OpenXQuest:structure"))
                {
                  current[j] =   static_cast<double>(pep_id2->getMetaValue("OpenXQuest:score"))
                      / static_cast<double>(pep_id1->getMetaValue("OpenXQuest:score"));
                  break;
                }
              }
            }
          }
          delta_scores.push_back(current);
        }
      }
      else
      {
        writeLog_("Input is not xQuest. Omit computing delta score and min. number of matched ions");
      }

      /*
       * Sort rank one hits in descending order according to the score.
       */

      // Contains the indices of the rank_one_ids sorted by the score (Init from 0,1,2,...,rank_one_ids.size()
      vector< UInt > order_score ( boost::counting_iterator<Size>(0),
                                   boost::counting_iterator<Size>(rank_one_ids.size()));

      // Configure the sorting of the Peptide Identifications
      greater_than_by_key order_conf = {
        all_ids,
        rank_one_ids,
      };

      std::sort(order_score.begin(), order_score.end(), order_conf);
      assert(std::is_sorted(order_score.begin(), order_score.end(), order_conf));

      // For unique IDs
      std::set<String> unique_ids;

      // Applies user specified filters and aggregates the scores for the corresponding classes
      std::map< String, vector< double > > scores;
      UInt num_flagged = 0;

      //-------------------------------------------------------------
      // Sort peptide ID based on the crosslink class and apply filters
      //-------------------------------------------------------------
      for (Size i = 0; i < n_ids; ++i)
      {
        // Extract required attributes of the peptide_identification (filter criteria)
        PeptideIdentification & pep_id = all_ids[rank_one_ids[order_score[i]]];
        double error_rel = 0;
        double delta_score = 0;
        Int ions_matched = 0;
        String id = "";
        if (is_xquest_input)
        {
          id = pep_id.getMetaValue("OpenXQuest:id").toString();
          error_rel = static_cast<double>(pep_id.getMetaValue("OpenXQuest:error_rel"));
          delta_score = (delta_scores[order_score[i]])[pep_id_index];
          ions_matched = n_min_ions_matched[order_score[i]];
        }
        num_flagged++;
        double score = getCrosslinkScore(pep_id);

        // Only consider peptide identifications which  fullfill all filter criteria specified by the user
        if (        (is_xquest_input ? (    (arg_minborder <= error_rel || arg_minborder == -1)
                                            && (arg_maxborder >= error_rel || arg_maxborder == -1)) : true)
                    && (is_xquest_input ? (mindelta_filter_disabled || delta_score < arg_mindeltas) : true)
                    && (is_xquest_input ? ions_matched  >= arg_minionsmatched : true)
                    &&  score >= arg_minscore
                    && (is_xquest_input ? ( ! arg_uniquex || unique_ids.find(id) == unique_ids.end()) : true))
        {
          pep_id.setMetaValue("OpenXQuest:xprophet_f", 1);
          unique_ids.insert(id);
          StringList crosslink_types;
          assignTypes(pep_id, crosslink_types);

          for (StringList::const_iterator crosslink_types_it = crosslink_types.begin(); crosslink_types_it != crosslink_types.end();
               ++crosslink_types_it)
          {
            scores[*crosslink_types_it].push_back(score);
          }
        }
      }
      writeLog_(this->toolName_() + " has used " + num_flagged + " hits to calculate the FDR");

      // Log number of scores within each class
      writeLog_("Number of Scores for each class:");

      for (std::map< String, vector< double > >::const_iterator scores_it = scores.begin();
           scores_it != scores.end(); ++scores_it)
      {
        std::pair< String, vector< double > > pair = *scores_it;
        writeLog_(pair.first + ": " + pair.second.size());
      }

      // Generate Histograms of the scores for each class
      // Use cumulative histograms to count the number of scores above consecutive thresholds
      std::map< String, Math::Histogram<> >  cum_histograms;
      for (std::map< String, vector< double > >::const_iterator scores_it = scores.begin();
           scores_it != scores.end(); ++scores_it)
      {
        vector< double > current_scores = scores_it->second;
        String classname = scores_it->first;
        Math::Histogram<> histogram(this->min_score, this->max_score, TOPPXFDR::fpnum_score_step);
        Math::Histogram<>::getCumulativeHistogram(current_scores.begin(), current_scores.end(), true, true, histogram);
        cum_histograms[classname] = histogram;
      }

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
        writeLog_("Performing qFDR transformation");

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
      
      // Assign FDR values to all identifications
      for (vector< PeptideIdentification >::iterator all_ids_it = all_ids.begin();
           all_ids_it != all_ids.end(); ++all_ids_it)
      {
        PeptideIdentification & pep_id = *all_ids_it;

        if ( ! pep_id.metaValueExists("OpenXQuest:xprophet_f"))
        {
          pep_id.setMetaValue("OpenXQuest:xprophet_f", 0);
        }
        double score = getCrosslinkScore(pep_id);
        
        StringList crosslink_types;
        assignTypes(pep_id, crosslink_types);

        pep_id.setMetaValue("OpenXQuest:fdr_type", score_type);
        pep_id.setScoreType(score_type);

        // Get PeptideHits
        vector< PeptideHit > & pep_hits = pep_id.getHits();
        Size n_hits = pep_hits.size();
        assert(n_hits == 1 || n_hits == 2);
        // Assign FDR value as meta value and also set as score
        bool assigned = false;
        double fdr;
        for (StringList::const_iterator crosslink_types_it = crosslink_types.begin();
             crosslink_types_it != crosslink_types.end(); ++crosslink_types_it)
        {
          String crosslink_type = *crosslink_types_it;
          Size idx = std::floor((score - this->min_score) / TOPPXFDR::fpnum_score_step);
          if (   crosslink_type == TOPPXFDR::crosslink_class_fulldecoysinterlinks
                 || crosslink_type == TOPPXFDR::crosslink_class_hybriddecoysinterlinks
                 || crosslink_type == TOPPXFDR::crosslink_class_interdecoys
                 || crosslink_type == TOPPXFDR::crosslink_class_interlinks)
          {
            fdr = fdr_interlinks[idx];
            assigned = true;
            break;
          }
          else if (   crosslink_type == TOPPXFDR::crosslink_class_fulldecoysintralinks
                      || crosslink_type == TOPPXFDR::crosslink_class_hybriddecoysintralinks
                      || crosslink_type == TOPPXFDR::crosslink_class_intradecoys
                      || crosslink_type == TOPPXFDR::crosslink_class_intralinks)
          {
            fdr = fdr_intralinks[idx];
            assigned = true;
            break;
          }
          else if (   crosslink_type == TOPPXFDR::crosslink_class_monodecoys
                      || crosslink_type == TOPPXFDR::crosslink_class_monolinks)
          {
            fdr = fdr_monolinks[idx];
            assigned = true;
            break;
          }
        }
        if ( assigned)
        {
            pep_id.setMetaValue("OpenXQuest:fdr", fdr);

            // Set FDR as score for the peptide hits
            for (Size i = 0; i < n_hits; ++i)
            {
              pep_hits[i].setScore(fdr);
            }

            // For FDR score, higher score is no longer better
            pep_id.setHigherScoreBetter(false);
        }
        else
        {
           LOG_WARN << "WARNING: Crosslink could not be identified as either interlink, intralink, or monolink, so no FDR will be available." << endl;
        }
      }
      // Write idXML
      if ( ! arg_out_idXML.empty())
      {
        IdXMLFile().store( arg_out_idXML, prot_ids, all_ids);
      }

      // Write mzid file
      if (! arg_out_mzid.empty())
      {
        MzIdentMLFile().store( arg_out_mzid, prot_ids, all_ids);
      }

      return EXECUTION_OK;
    }

  private:

    // Score range for this exection of the tool
    Int min_score;
    Int max_score;
};
const String TOPPXFDR::param_in = "in";
const String TOPPXFDR::param_in_type = "in_type";
const String TOPPXFDR::param_out_idXML = "out_idXML";
const String TOPPXFDR::param_out_mzid = "out_mzIdentML";
const String TOPPXFDR::param_minborder = "minborder";
const String TOPPXFDR::param_maxborder = "maxborder";
const String TOPPXFDR::param_mindeltas = "mindeltas";
const String TOPPXFDR::param_minionsmatched = "minionsmatched";
const String TOPPXFDR::param_uniquexl = "uniquexl";
const String TOPPXFDR::param_no_qvalues = "no_qvalues";
const String TOPPXFDR::param_minscore = "minscore";

const UInt TOPPXFDR::n_rank = 1; //  Number of ranks used

const String TOPPXFDR::crosslink_class_intradecoys = "intradecoys";
const String TOPPXFDR::crosslink_class_fulldecoysintralinks = "fulldecoysintralinks";
const String TOPPXFDR::crosslink_class_interdecoys = "interdecoys";
const String TOPPXFDR::crosslink_class_fulldecoysinterlinks = "fulldecoysinterlinks";
const String TOPPXFDR::crosslink_class_monodecoys = "monodecoys";
const String TOPPXFDR::crosslink_class_intralinks = "intralinks";
const String TOPPXFDR::crosslink_class_interlinks = "interlinks";
const String TOPPXFDR::crosslink_class_monolinks  = "monolinks";
const String TOPPXFDR::crosslink_class_decoys = "decoys";
const String TOPPXFDR::crosslink_class_hybriddecoysintralinks = "hybriddecoysintralinks";
const String TOPPXFDR::crosslink_class_hybriddecoysinterlinks = "hybriddecoysinterlinks";

// Parameters for actually calculating the number of FPs
const double TOPPXFDR::fpnum_score_step = 0.1;

// meta values for crosslink identifications
const String TOPPXFDR::crosslink_type = "xl_type";
const String TOPPXFDR::crosslink_rank = "xl_rank";
const String TOPPXFDR::target_decoy = "target_decoy";

// the actual main function needed to create an executable
int main(int argc, const char ** argv)
{
  TOPPXFDR tool;
  return tool.main(argc, argv);
}
/// @endcond
