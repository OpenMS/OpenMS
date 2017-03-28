// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
#include <OpenMS/METADATA/XQuestResultMeta.h>
#include <OpenMS/FORMAT/XQuestResultXMLFile.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/MATH/STATISTICS/Histogram.h>
#include <OpenMS/MATH/STATISTICS/CumulativeHistogram.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>


#include <boost/iterator/counting_iterator.hpp>
#include <boost/function.hpp>

#include <string>
#include <math.h>

using namespace OpenMS;
using namespace std;


// TODO Switch to Debug mode in CMake and remove this undef
#undef NDEBUG
#include <assert.h>

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_NewTool

    @brief Template for a new Tool

    This tool can be used for scientific stuff.

    And more scientific applications.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_NewTool.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_NewTool.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


// Struct which is used for sorting the order vector by a meta value
// of the peptide identification (usually score)
/**
 * Struct which is used for sorting the order vector by a meta value of the peptide identification
 * containing the identified cross-link.
 * @brief The less_than_by_key struct
 */

/*
 *
 * Less than by key which used the vector < vector < PeptideIdentification data structure /
 *
struct less_than_by_key
{
    vector< vector< PeptideIdentification > > & elements;  // The elements according to which the rankm vector should be computed
    size_t idx;  // Which element in the internal vector should be used (normally determined by rank, and currently only rank 1 is supported)
    const String & meta_value;  // By which meta_value each peptide identification should be sorted (usually by score)

    inline bool operator()(const size_t & index1, const size_t & index2)
    {
      return ((double) elements[index1][idx].getMetaValue(meta_value)
              > (double) elements[index2][idx].getMetaValue(meta_value));
    }
};
*/


inline void addEmptyClass(std::map< String, vector< double> > & scores, const String & name) 
{
  if (scores.find(name) == scores.end())
  {
    scores.insert( std::pair<String, vector< double> >(name, vector< double>()));
  }
}

void removeSubstring(String &  large, const String & small)
{
  std::string::size_type i = large.find(small);
  if (i != std::string::npos)
  {
    large.erase(i, small.length());
  }
}


/**
   Look up Meta Value for cross-link identification first in PeptideIdentification, then in the PeptideHits 
*/
template<typename T> 
T getXLMetaValue(const String & name, const PeptideIdentification & pep_id, bool is_score)
{
  if( pep_id.metaValueExists(name))
  {
    return (T) pep_id.getMetaValue(name);   
  }
  vector< PeptideHit > pep_hits = pep_id.getHits();
  
  if (is_score)
  {
    return (T) pep_hits[0].getScore();
  }
      
  if (pep_hits[0].metaValueExists(name))
  {  
    return (T) pep_hits[0].getMetaValue(name);
  }
  if (pep_hits.size() == 2 && pep_hits[1].metaValueExists(name))
  {
    return (T) pep_hits[1].getMetaValue(name);
  }
  return DataValue::EMPTY;
}


/** 
   Prepares the peptide identifications of cross-links such that they can be further processed
 */
void prepareIDXML(vector< PeptideIdentification > & pep_ids, vector< ProteinIdentification > & prot_ids)
{
  for (vector< PeptideIdentification >::iterator pep_ids_it = pep_ids.begin();
       pep_ids_it != pep_ids.end(); ++pep_ids_it)
  {
   PeptideIdentification & pep_id = *pep_ids_it;
   
   // SetPeptideIdentification to Decoy is either alpha or beta or both are decoys
   vector< PeptideHit > pep_hits = pep_id.getHits();
    
   if (    pep_hits[0].getMetaValue("target_decoy").toString() == "decoy" 
       || (pep_hits.size() == 2 && pep_hits[1].getMetaValue("target_decoy").toString() == "decoy"))
   {
     pep_id.setMetaValue("target_decoy", DataValue("decoy"));
   }
   else
   {
     pep_id.setMetaValue("target_decoy", DataValue("target"));  
   } 
  
  if (pep_hits.size() == 2)
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
         
         removeSubstring(alpha_prot, "reverse_");
         removeSubstring(beta_prot, "decoy_");
         removeSubstring(alpha_prot, "reverse_");
         removeSubstring(beta_prot, "decoy_"); 
         pep_id.setMetaValue( alpha_prot == beta_prot ? "OpenXQuest:is_intraprotein" : "OpenXQuest:is_interprotein" , DataValue());    
      }            
    }     
  }
 
  // Determine whether inter or intra protein
   assert(pep_id.metaValueExists("target_decoy"));
  }
}



struct less_than_by_key
{
    //Indizes of the rank one elements within the all_ids vector
    vector < PeptideIdentification > & all_ids;
    vector< Size >  & rank_one_ids;
    const String & meta_value;  // By which meta_value each peptide identification should be sorted (usually by some score)

    inline bool operator()(const size_t & index1, const size_t & index2)
    {
      return (  getXLMetaValue<double>(meta_value, all_ids[rank_one_ids[index1]], true)
              > getXLMetaValue<double>(meta_value, all_ids[rank_one_ids[index2]], true));
    }
};




/**
 * Ensures that the order vector is sorted in descending order with respect to the Peptide Identifications (meta value: Score)
 * This is only used for assertions.
 *
 */
bool isSortedDescending(vector< Size > & order, vector< PeptideIdentification > & all_ids, vector< Size > & rank_one_ids )
{
  for (Size i = 0; i < order.size() - 1; ++i)
  {
    double a_score = getXLMetaValue<double>("OpenXQuest:score", all_ids[rank_one_ids[order[i]]], true);
    double b_score = getXLMetaValue<double>("OpenXQuest:score", all_ids[rank_one_ids[order[i+1]]], true);

    if(a_score < b_score)
    {
      return false;
    }
  }
  return true;
}


class TOPPXFDR :
    public TOPPBase
{
  public:

    static const String param_in;  // Parameter for the input file
    static const String param_out_idXML;
    static const String param_minborder;  // minborder -5 # filter for minimum precursor mass error (ppm)
    static const String param_maxborder;  // maxborder  5 # filter for maximum precursor mass error (ppm)
    static const String param_mindeltas;  // mindeltas  0.95 # filter for delta score, 0 is no filter, minimum delta score required, hits are rejected if larger or equal
    static const String param_minionsmatched; // minionsmatched 0 # Filter for minimum matched ions per peptide
    static const String param_uniquexl; // calculate statistics based on unique IDs
    static const String param_qtransform; // transform simple FDR to q-FDR values
    static const String param_minscore; // minscore 0 # minimum ld-score to be considered
    static const String param_verbose; // Whether or not the output of the tool should be verbose.

    // Number of ranks used
    static const Int n_rank;

    // Constants related to particular xl classes
    // decoy classes
    static const String xlclass_intradecoys; // intradecoys
    static const String xlclass_fulldecoysintralinks; // fulldecoysintralinks
    static const String xlclass_interdecoys; // interdecoys
    static const String xlclass_fulldecoysinterlinks; // fulldecoysinterlinks
    static const String xlclass_intralinks; // intralinks
    static const String xlclass_interlinks; // interlinks
    static const String xlclass_monolinks;  // monolinks
    static const String xlclass_monodecoys; // monodecoys
    static const String xlclass_decoys; // decoys
    static const String xlclass_hybriddecoysintralinks; // hybriddecoysintralinks
    static const String xlclass_hybriddecoysinterlinks; // hybriddecoysintralinks
    static const StringList xlclasses;

    // Parameters to actually calculate the number of FPs // TODO Needs to more dynamic in the future
    static const double fpnum_score_start;
    static const double fpnum_score_end;
    static const double fpnum_score_step;

    TOPPXFDR() :
      TOPPBase("XFDR", "Template for Tool creation", false)
    {
    }

  protected:

   inline void assign_types(PeptideIdentification &pep_id, StringList & types)
   {
      types.clear();
      bool pep_is_decoy = pep_id.getMetaValue("target_decoy").toString() == "decoy";
     
      // Intradecoys
      if (pep_id.metaValueExists("OpenXQuest:is_intraprotein") && pep_is_decoy)
      {
        types.push_back(TOPPXFDR::xlclass_intradecoys);
      }

      // Decoys
      if (pep_is_decoy)
      {
        types.push_back(TOPPXFDR::xlclass_decoys);
      }

      // intralinks
      if (pep_id.metaValueExists("OpenXQuest:is_intraprotein") && ! pep_is_decoy)
      {
        types.push_back(TOPPXFDR::xlclass_intralinks);
      }

      // interdecoys
      if (pep_id.metaValueExists("OpenXQuest:is_interprotein") && pep_is_decoy)
      {
        types.push_back(TOPPXFDR::xlclass_interdecoys);
      }

      // interlinks
      if (pep_id.metaValueExists("OpenXQuest:is_interprotein") && ! pep_is_decoy)
      {
        types.push_back(TOPPXFDR::xlclass_interlinks);
      }
      String xl_type = getXLMetaValue<String>("xl_type", pep_id, false);
     
      // monolinks
      if ( ! pep_is_decoy && (xl_type == "mono-link"
                          ||  xl_type == "loop-link"))
      {
        types.push_back(TOPPXFDR::xlclass_monolinks);
      }

      // monodecoys
      if ( pep_is_decoy && (xl_type == "mono-link"
                        ||  xl_type == "loop-link"))
      {
        types.push_back(TOPPXFDR::xlclass_monodecoys);
      }
      vector< PeptideHit > pep_hits = pep_id.getHits();
      if(pep_hits.size() == 2)
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
          types.push_back(TOPPXFDR::xlclass_fulldecoysintralinks);
        }

        // fulldecoysinterlinks
        if (   pep_id.metaValueExists("OpenXQuest:is_interprotein")
            && alpha_is_decoy
            && beta_is_decoy)
        {
          types.push_back(TOPPXFDR::xlclass_fulldecoysinterlinks);

        }

        // hybriddecoysintralinks
        if (       pep_id.metaValueExists("OpenXQuest:is_intraprotein")
            && ( ! alpha_is_decoy
            &&     beta_is_decoy
            ||     alpha_is_decoy
            &&   ! beta_is_decoy))
        {
          types.push_back(TOPPXFDR::xlclass_hybriddecoysintralinks);
        }

        // hybriddecoysinterlinks
        if (       pep_id.metaValueExists("OpenXQuest:is_interprotein")
            && ( ! alpha_is_decoy
            &&     beta_is_decoy
            ||     alpha_is_decoy
            &&   ! beta_is_decoy))
        {
          types.push_back(TOPPXFDR::xlclass_hybriddecoysinterlinks);
        }
      }
    }


    /** False positve counting as performed by xProphet software package.
      *
      * @brief xprophet  method for FP counting as implemented in xProphet
      * @param cum_histograms Cumulative score distributions
      * @param fp_counts Number of FPs for each score threshold.
      */

    /*
     void fp_xprophet(std::map< String, Math::CumulativeHistogram<>  * > & cum_histograms,
                       Math::Histogram<> & fp_counts)
     {

       // Required cumulative histograms for FPs
       Math::CumulativeHistogram<>  intradecoys_histogram = *cum_histograms[TOPPXFDR::xlclass_intradecoys];
       Math::CumulativeHistogram<>  fulldecoysintralinks_histogram = *cum_histograms[TOPPXFDR::xlclass_fulldecoysintralinks];
       Math::CumulativeHistogram<>  interdecoys_histogram = *cum_histograms[TOPPXFDR::xlclass_interdecoys];
       Math::CumulativeHistogram<>  fulldecoysinterlinks_histogram = *cum_histograms[TOPPXFDR::xlclass_fulldecoysinterlinks];
       Math::CumulativeHistogram<>  monodecoys_histogram = *cum_histograms[TOPPXFDR::xlclass_monodecoys];

       for (double current_score = TOPPXFDR::fpnum_score_start +  (TOPPXFDR::fpnum_score_step/2) ;
            current_score <= TOPPXFDR::fpnum_score_end - (TOPPXFDR::fpnum_score_step/2);
            current_score += TOPPXFDR::fpnum_score_step)
       {
           UInt n_intrafp = intradecoys_histogram.binValue(current_score) - 2 * fulldecoysintralinks_histogram.binValue(current_score);
           UInt n_interfp = interdecoys_histogram.binValue(current_score) - 2 * fulldecoysinterlinks_histogram.binValue(current_score);
           UInt n_monofp = monodecoys_histogram.binValue(current_score);
           // Aim for the center of the bin when inserting the score
           fp_counts.inc(current_score, n_interfp + n_intrafp + n_monofp);
       }
    }
    */

    /** Target counting as performed by the xProphet software package
       *
       * @brief xprophet  method for target hits counting as implemented in xProphet
       * @param cum_histograms Cumulative score distributions
       * @param target_counts Number of target hits
       */
    /*
      void target_xprophet(std::map< String, Math::CumulativeHistogram<>  * > & cum_histograms,
                        Math::Histogram<> & target_counts)
      {
        // Required cumulative histograms for target hits
        Math::CumulativeHistogram<>  intralinks_histogram = *cum_histograms[TOPPXFDR::xlclass_intralinks];
        Math::CumulativeHistogram<>  interlinks_histogram = *cum_histograms[TOPPXFDR::xlclass_interlinks];
        Math::CumulativeHistogram<>  monolinks_histogram = *cum_histograms[TOPPXFDR::xlclass_monolinks];

        for (double current_score = TOPPXFDR::fpnum_score_start +  (TOPPXFDR::fpnum_score_step/2) ;
             current_score <= TOPPXFDR::fpnum_score_end - (TOPPXFDR::fpnum_score_step/2);
             current_score += TOPPXFDR::fpnum_score_step)
        {
            UInt n_intralinks = intralinks_histogram.binValue(current_score);
            UInt n_interlinks = interlinks_histogram.binValue(current_score);
            UInt n_monolinks = monolinks_histogram.binValue(current_score);
            target_counts.inc(current_score, n_intralinks + n_interlinks + n_monolinks);
        }
     }
     */
    

    /** Target counting as performed by the xProphet software package
      *
      * @brief xprophet  method for target hits counting as implemented in xProphet
      * @param cum_histograms Cumulative score distributions
      */
    void fdr_xprophet(std::map< String, Math::CumulativeHistogram<>  * > & cum_histograms,
                      const String  & targetclass, const String & decoyclass, const String & fulldecoyclass,
                      vector< double > & fdr, bool mono)
    {
      for (double current_score = TOPPXFDR::fpnum_score_start +  (TOPPXFDR::fpnum_score_step/2) ;
           current_score <= TOPPXFDR::fpnum_score_end - (TOPPXFDR::fpnum_score_step/2);
           current_score += TOPPXFDR::fpnum_score_step)
      {
        double estimated_n_decoys = cum_histograms[decoyclass]->binValue(current_score);
        if ( ! mono)
        {
          estimated_n_decoys -= 2 * cum_histograms[fulldecoyclass]->binValue(current_score);
        }
        double n_targets = cum_histograms[targetclass]->binValue(current_score);
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
      for (int i = fdr.size(); i > -1; --i)
      {
        double current_fdr = fdr[i];
        double smallest_fdr = current_fdr;
        for (int j = i; j > -1; j--)
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

    // this function will be used to register the tool parameters
    // it gets automatically called on tool execution
    void registerOptionsAndFlags_()
    {
      // Verbose Flag
      registerFlag_(TOPPXFDR::param_verbose, "Whether the log of information will be loud and noisy");


      // File input
      registerInputFile_(TOPPXFDR::param_in, "<file>", "", "Results in the original xquest.xml format", false);
      setValidFormats_(TOPPXFDR::param_in, ListUtils::create<String>("xml,mzid,idXML"));
      
      // idXML output
      registerOutputFile_(TOPPXFDR::param_out_idXML, "<idXML_file>", "", "Output as idXML file", false, false);
      setValidFormats_(TOPPXFDR::param_out_idXML, ListUtils::create<String>("idXML"));

      // Minborder
      registerIntOption_(TOPPXFDR::param_minborder, "<minborder>", -5, "Filter for minimum precursor mass error (ppm)", false);

      // Maxborder
      registerIntOption_(TOPPXFDR::param_maxborder, "<maxborder>",  5, "Filter for maximum precursor mass error (ppm)", false);

      // Mindeltas
      registerDoubleOption_(TOPPXFDR::param_mindeltas, "<mindeltas>", 0.95, "Filter for delta score, 0 is no filter, minimum delta score required, hits are rejected if larger or equal", false);

      // Minionsmatched
      registerIntOption_(TOPPXFDR::param_minionsmatched, "<minionsmatched>", 0, "Filter for minimum matched ions per peptide", false);

      // Uniquexl
      registerFlag_(TOPPXFDR::param_uniquexl, "Calculate statistics based on unique IDs");

      // Qtransform
      registerFlag_(TOPPXFDR::param_qtransform, "Transform simple FDR to q-FDR values");

      // Minscore
      registerIntOption_(TOPPXFDR::param_minscore, "<minscore>", 0, "Minimum ld-score to be considered", false);

    }

    // the main_ function is called after all parameters are read
    ExitCodes main_(int, const char **)
    {
      bool arg_verbose = getFlag_(TOPPXFDR::param_verbose);
      bool is_xquest_input = false;

      //-------------------------------------------------------------
      // parsing parameters
      //-------------------------------------------------------------
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
      Int arg_minionsmatched = getIntOption_(TOPPXFDR::param_minionsmatched);
      Int arg_minscore = getIntOption_(TOPPXFDR::param_minscore);
      bool arg_uniquex = getFlag_(TOPPXFDR::param_uniquexl);

      //-------------------------------------------------------------
      // Printing parameters
      //-------------------------------------------------------------
      LOG_INFO << "Filtering of precursor mass error from " << arg_minborder << " to " << arg_maxborder << " ppm is used.\n"
               << "Filtering of hits by a deltascore of " << arg_mindeltas << " is used.\n";
      if (arg_minionsmatched > 0)
      {
        LOG_INFO << "Filtering of hits by minimum ions matched: " << arg_minionsmatched << " is used\n";
      }
      else
      {
        LOG_INFO << "No filtering of hits by minimum ions matched.\n";
      }
      if (arg_minscore > 0)
      {
        LOG_INFO << "Filtering of hits by minimum score of " << arg_minscore << " is used.\n";
      }
      else
      {
        LOG_INFO << "No filtering of hits by minimum score.\n";
      }
      if (arg_uniquex)
      {
        LOG_INFO << "Error model is generated based on unique cross-links.\n";
      }
      else
      {
        LOG_INFO << "Error model is generated based on redundant cross-links.\n";
      }
      if(arg_verbose)
      {
        LOG_INFO << "Output will be verbose.\n";
      }
      else
      {
        LOG_INFO << "Output will not be verbose\n";
      }
      LOG_INFO << "-----------------------------------------\n" << endl;

      //-------------------------------------------------------------
      // Parse the input file
      //-------------------------------------------------------------
      String arg_in = getStringOption_(TOPPXFDR::param_in);

      if (arg_verbose)
      {
        LOG_INFO << "INFO: Parsing input file: " << arg_in << endl;
      }
      Size pep_id_index = TOPPXFDR::n_rank - 1; // This is of course trash if more than 1 ranks are used.
      Size n_spectra;

      vector < PeptideIdentification > all_ids;
      vector < ProteinIdentification > prot_ids;
      
      vector < Size > rank_one_ids; // Stores the indizes of the rank one hits within all_ids

      vector < vector < PeptideIdentification > > spectra;

      // assume XQuestXML here
      if (arg_in.hasSuffix("xml"))
      {
        is_xquest_input = true;
        
        // Core data structures for the util
        vector< XQuestResultMeta > metas;
        // Parse XQuestResultXMLFile (TODO Also support idXML and mzIdentML)
        XQuestResultXMLFile xquest_result_file;
        xquest_result_file.load(arg_in, metas, spectra, false, 1); // We do not load 'empty' spectra here
        n_spectra = spectra.size();

        if (arg_verbose)
        {
          LOG_INFO << "Total number of spectra: " << n_spectra << "\n"
                   << "Total number of hits: "    << xquest_result_file.get_n_hits() << endl;
        }

        Size rank_counter = 0;
        for (vector < vector < PeptideIdentification > >::const_iterator spectra_it = spectra.begin();
             spectra_it != spectra.end(); ++spectra_it)
        {
          vector< PeptideIdentification > spectrum = *spectra_it;
          for (vector< PeptideIdentification >::const_iterator spectrum_it = spectrum.begin(); spectrum_it != spectrum.end(); ++spectrum_it)
          {
             PeptideIdentification pep_id = *spectrum_it;
             all_ids.push_back(pep_id);
             if( getXLMetaValue<int>("xl_rank", pep_id, false) == 1)
             {
                rank_one_ids.push_back(rank_counter);
             }
             rank_counter++;
          }
        }
      }
      
      // TODO Add support
      else if (arg_in.hasSuffix("mzid"))
      {
        vector< ProteinIdentification > prot_ids;
        vector< PeptideIdentification > pep_ids;
        MzIdentMLFile().load(arg_in, prot_ids, pep_ids);

        PeptideIdentification pep_id = pep_ids[0];
        StringList keys;
        pep_id.getKeys(keys);

        cout << "PRINT" << endl;
        for (StringList::const_iterator keys_it = keys.begin(); keys_it != keys.end(); ++keys_it)
        {
          cout << *keys_it << endl;
        }
      }
      else if (arg_in.hasSuffix("idXML"))
      {
        IdXMLFile().load(arg_in, prot_ids, all_ids); 
        prepareIDXML(all_ids, prot_ids);
       
        for (vector< PeptideIdentification >::const_iterator all_ids_it = all_ids.begin();
             all_ids_it != all_ids.end(); ++all_ids_it)
        {
           assert(all_ids_it->metaValueExists("target_decoy"));
        }
    
         
        Size rank_counter = 0;
        for (vector< PeptideIdentification >::const_iterator all_ids_it = all_ids.begin();
             all_ids_it != all_ids.end(); ++all_ids_it)
        {
          PeptideIdentification pep_id = *all_ids_it;
          
          if (getXLMetaValue<int>("xl_rank", pep_id, false) == 1)
          {
            rank_one_ids.push_back(rank_counter);
          }
          rank_counter++;
        }

      }
      
      // Number of peptide identifications that need to be considered
      Size n_ids = rank_one_ids.size();
      
   
      //      for(vector< vector < PeptideIdentification > >::const_iterator it = spectra.begin(); it != spectra.end(); ++it)
      //      {
      //        vector< PeptideIdentification > csm = *it;
      //        cout << csm.size() << endl;
      //        for(vector< PeptideIdentification >::const_iterator it2 = csm.begin(); it2 != csm.end(); ++it2)
      //        {
      //            cout << it2->getHits().size() << endl;
      //        }
      //      }

      //-------------------------------------------------------------
      // Calculate the delta score for each hit
      // Calculate n_min_ions_matched
      // Currently only for xQuest input files
      //-------------------------------------------------------------
      // The score is calculated for each hit h on the set of all hits of the spectrum that encompasses
      std::vector< std::vector< double >* > delta_scores;
      std::vector< size_t > n_min_ions_matched;


      // For xQuest input,calculate delta scores and min_ions_matched
      if (is_xquest_input)
      {
        if (arg_verbose)
        {
          LOG_INFO << "Input is XQuest. Compute the delta scores and the number of matched ions" << endl;
        }
        delta_scores.resize(n_spectra);
        n_min_ions_matched.resize(n_spectra);

        for(size_t i = 0; i < delta_scores.size(); ++i)
        {
          size_t n_hits = spectra[i].size();
          delta_scores[i] = new std::vector<double>(n_hits);
          vector<double> * current = delta_scores[i];

          assert(n_hits > 0); // because we initially do not load 'empty' spectra
          // calculate n_min_ions_matched
          PeptideIdentification * pep_id1 = &spectra[i][0];

          // Currently the correct rank order in the xQuest result file is assumed
          //assert((int) pep_id1->getMetaValue("xl_rank") == 1); // because hits are sorted according to their rank within the spectrum

          vector<PeptideHit> pep_hits = pep_id1->getHits();

          if( pep_id1->getMetaValue("xl_type") == "cross-link")
          {
            n_min_ions_matched[i] = std::min((int) pep_hits[0].getMetaValue("OpenXQuest:num_of_matched_ions"),
                (int) pep_hits[1].getMetaValue("OpenXQuest:num_of_matched_ions"));
          }
          else
          {
            n_min_ions_matched[i] = (int) pep_hits[0].getMetaValue("OpenXQuest:num_of_matched_ions");
          }
          // Calculate delta score
          if (n_hits > 1)
          {
            for (size_t j = 0; j < n_hits - 1; ++j)
            {
              pep_id1 = &spectra[i][j];
              for (size_t k = 1; j+k < n_hits; ++k )
              {
                PeptideIdentification * pep_id2 = &spectra[i][j+k];
                if(pep_id1->getMetaValue("OpenXQuest:structure") != pep_id2->getMetaValue("OpenXQuest:structure"))
                {
                  (*current)[j] =  ((double) pep_id2->getMetaValue("OpenXQuest:score"))
                      / ((double) pep_id1->getMetaValue("OpenXQuest:score"));
                  break;
                }
              }
            }
          }
        }
      }
      else if (arg_verbose)
      {
        LOG_INFO << "Input is not xQuest. Omit computing delta score and min. number of matched ions." << endl;
      }
     
      /*
       * Sort rank one hits in descending order according to the score
       */
      typedef std::vector<size_t> ranks;
    
      ranks order_score ( boost::counting_iterator<size_t>(0),
                          boost::counting_iterator<size_t>(rank_one_ids.size()));

      // Configure the sorting of the Peptide Identifications
      less_than_by_key order_conf = {
        all_ids,
        rank_one_ids,            // elements
        "OpenXQuest:score"  // key
      };
     
     std::sort(order_score.begin(), order_score.end(), order_conf);     
     assert(isSortedDescending(order_score, all_ids, rank_one_ids));
 
      // For unique IDs
      std::set<String> unique_ids;

      // Applies user specified filters and aggregates the scores for the corresponding classes
      std::map< String, vector< double > > scores;

      // Fetch all relevant peptide identifications in descending order of score
      for (Size i = 0; i < n_ids; ++i)
      {    
        // Extract required attributes of the peptide_identification (filter criteria)
        PeptideIdentification & pep_id = all_ids[rank_one_ids[order_score[i]]];
        
        double error_rel;
        double delta_score;
        Size ions_matched;
        String id;
        if (is_xquest_input)
        {
          id = (String) pep_id.getMetaValue("OpenXQuest:id");
          error_rel = (double) pep_id.getMetaValue("OpenXQuest:error_rel");
          delta_score = (*(delta_scores[order_score[i]]))[pep_id_index];
          ions_matched = n_min_ions_matched[order_score[i]];
        }     
        double score = getXLMetaValue<double>("OpenXQuest:score", pep_id, true); 
        
        // Only consider peptide identifications which  fullfill all filter criteria specified by the user
        if (        (is_xquest_input ? (    arg_minborder <= error_rel
                                         && arg_maxborder >= error_rel ) : true)
                 && (is_xquest_input ? (mindelta_filter_disabled || delta_score < arg_mindeltas) : true)  // Only apply for xQuest Input
                 && (is_xquest_input ? ions_matched  >= arg_minionsmatched : true)                       // Only apply for xQuest Input
                 &&  score >= arg_minscore
                 && (is_xquest_input ? ( ! arg_uniquex || unique_ids.find(id) == unique_ids.end()) : true))
        {
          pep_id.setMetaValue("OpenXQuest:xprophet_f", 1);
          unique_ids.insert(id);
          StringList xl_types;
          
          assign_types(pep_id, xl_types); // TODO Currently hard-coded (maybe similar to onthologies).
          
          for (StringList::const_iterator xl_types_it = xl_types.begin(); xl_types_it != xl_types.end(); ++xl_types_it)
          {
            // Assign score to this xl type
            scores[*xl_types_it].push_back(score);
          }
        }
      }
      // Push empty vector for all remaining empty classes  
      addEmptyClass(scores, TOPPXFDR::xlclass_intradecoys);
      addEmptyClass(scores, TOPPXFDR::xlclass_fulldecoysintralinks);
      addEmptyClass(scores, TOPPXFDR::xlclass_interdecoys);
      addEmptyClass(scores, TOPPXFDR::xlclass_fulldecoysinterlinks );
      addEmptyClass(scores, TOPPXFDR::xlclass_monodecoys );
      addEmptyClass(scores, TOPPXFDR::xlclass_intralinks );
      addEmptyClass(scores, TOPPXFDR::xlclass_interlinks );
      addEmptyClass(scores, TOPPXFDR::xlclass_monolinks );
      addEmptyClass(scores, TOPPXFDR::xlclass_decoys );
      addEmptyClass(scores, TOPPXFDR::xlclass_hybriddecoysintralinks );
      addEmptyClass(scores, TOPPXFDR::xlclass_hybriddecoysinterlinks );
      
     
      // Print number of scores within each class
      if (arg_verbose)
      {
        LOG_INFO << "\nINFO: Number of Scores for each class:" << endl;

        for (std::map< String, vector< double > >::const_iterator scores_it = scores.begin();
             scores_it != scores.end(); ++scores_it)
        {
          std::pair< String, vector< double > > pair = *scores_it;
          LOG_INFO << pair.first << ": " << pair.second.size() << endl;
        }
      }
     

      // Generate Histograms of the scores for each class
      // Use cumulative histograms to count the number of scores above consecutive thresholds

      //std::map< String, Math::Histogram<> * > histograms;
      std::map< String, Math::CumulativeHistogram<>  * >  cum_histograms;
      for (std::map< String, vector< double > >::const_iterator scores_it = scores.begin();
           scores_it != scores.end(); ++scores_it)
      {
        vector< double > current_scores = scores_it->second;
        String classname = scores_it->first;
        //histograms[classname] = new Math::Histogram<>(current_scores.begin(), current_scores.end(), 0, 100, 1);
        cum_histograms[classname] = new Math::CumulativeHistogram<>(current_scores.begin(), current_scores.end(),
                                                                    TOPPXFDR::fpnum_score_start,
                                                                    TOPPXFDR::fpnum_score_end,
                                                                    TOPPXFDR::fpnum_score_step, true, true);
      }
      
      // This is currently not needed for the FDR calculation
      //Math::Histogram<> fp_counts(TOPPXFDR::fpnum_score_start, TOPPXFDR::fpnum_score_end, TOPPXFDR::fpnum_score_step);
      //this->fp_xprophet(cum_histograms, fp_counts);

      //Math::Histogram<> target_counts(TOPPXFDR::fpnum_score_start, TOPPXFDR::fpnum_score_end, TOPPXFDR::fpnum_score_step);
      //this->target_xprophet(cum_histograms, target_counts);
      
      // Calculate FDR for interlinks
      vector< double > fdr_interlinks;   
      this->fdr_xprophet(cum_histograms, TOPPXFDR::xlclass_interlinks, TOPPXFDR::xlclass_interdecoys, TOPPXFDR::xlclass_fulldecoysinterlinks, fdr_interlinks, false);
     
      // Calculate FDR for intralinks
      vector< double > fdr_intralinks; 
      this->fdr_xprophet(cum_histograms, TOPPXFDR::xlclass_intralinks, TOPPXFDR::xlclass_intradecoys, TOPPXFDR::xlclass_fulldecoysintralinks, fdr_intralinks, false);
      
      // Calculate FDR for monolinks and looplinks
      vector< double > fdr_monolinks;
      this->fdr_xprophet(cum_histograms, TOPPXFDR::xlclass_monolinks, TOPPXFDR::xlclass_monodecoys, "", fdr_monolinks, true);

      // Determine whether qTransform should be performed
      bool arg_qtransform = getFlag_(TOPPXFDR::param_qtransform);


      if(arg_qtransform)
      {
        LOG_INFO << "Performing qFDR transformation" << endl;

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
           all_ids_it != all_ids.end(); all_ids_it++)
      {
        PeptideIdentification & pep_id = *all_ids_it;
       
        if ( ! pep_id.metaValueExists("OpenXQuest:xprophet_f"))
        {
          pep_id.setMetaValue("OpenXQuest:xprophet_f", 0);
        }       
        double score = getXLMetaValue<double>("OpenXQuest:score", pep_id, true);
        
        StringList xl_types;
        assign_types(pep_id, xl_types);
        pep_id.setMetaValue("OpenXQuest:fdr_type", arg_qtransform ? "qfdr" : "fdr");     
       
        // Assign FDR value
        bool assigned = false;
        for(StringList::const_iterator xl_types_it = xl_types.begin(); xl_types_it != xl_types.end(); xl_types_it++)
        {
          String xl_type = *xl_types_it;
          Size idx = std::floor((score - TOPPXFDR::fpnum_score_start) / TOPPXFDR::fpnum_score_step);
          if (   xl_type == TOPPXFDR::xlclass_fulldecoysinterlinks
              || xl_type == TOPPXFDR::xlclass_hybriddecoysinterlinks
              || xl_type == TOPPXFDR::xlclass_interdecoys
              || xl_type == TOPPXFDR::xlclass_interlinks)
          {
            pep_id.setMetaValue("OpenXQuest:fdr", fdr_interlinks[idx]);
            assigned = true;
            break;
          }
          else if (   xl_type == TOPPXFDR::xlclass_fulldecoysintralinks
                   || xl_type == TOPPXFDR::xlclass_hybriddecoysintralinks
                   || xl_type == TOPPXFDR::xlclass_intradecoys
                   || xl_type == TOPPXFDR::xlclass_intralinks)
          {
            pep_id.setMetaValue("OpenXQuest:fdr", fdr_intralinks[idx]);
            assigned = true;
            break;
          }
          else if (   xl_type == TOPPXFDR::xlclass_monodecoys
                   || xl_type == TOPPXFDR::xlclass_monolinks)
          {
            pep_id.setMetaValue("OpenXQuest:fdr", fdr_monolinks[idx]);
            assigned = true;
            break;
          }
        }
        if ( ! assigned)
        {
          LOG_WARN << "WARNING: Cross-link could not be identified as either interlink, intralink, or monolink, so no FDR will be available." << endl;
        }
      }     

      // Delete Delta Scores
      for (std::vector< std::vector< double >* >::const_iterator it = delta_scores.begin(); it != delta_scores.end(); ++it)
      {
        delete *it;
      }

      // Delete cumulative_histograms
      for (std::map< String, Math::CumulativeHistogram<> * >::iterator cum_histograms_it = cum_histograms.begin();
           cum_histograms_it != cum_histograms.end(); ++cum_histograms_it)
      {
        delete cum_histograms_it->second;
      }

     // Write output
     String arg_out_idXML = getStringOption_(TOPPXFDR::param_out_idXML);
     
     if ( ! arg_out_idXML.empty())
     {
       IdXMLFile().store( arg_out_idXML, prot_ids, all_ids);   
     }
    
      return EXECUTION_OK;
    }
};
const String TOPPXFDR::param_in = "in";
const String TOPPXFDR::param_out_idXML = "out_idXML";
const String TOPPXFDR::param_minborder = "minborder";
const String TOPPXFDR::param_maxborder = "maxborder";
const String TOPPXFDR::param_mindeltas = "mindeltas";
const String TOPPXFDR::param_minionsmatched = "minionsmatched";
const String TOPPXFDR::param_uniquexl = "uniquexl";
const String TOPPXFDR::param_qtransform = "qtransform";
const String TOPPXFDR::param_minscore = "minscore";
const String TOPPXFDR::param_verbose = "verbose";

const Int    TOPPXFDR::n_rank = 1; //  Number of ranks used

const String TOPPXFDR::xlclass_intradecoys = "intradecoys";
const String TOPPXFDR::xlclass_fulldecoysintralinks = "fulldecoysintralinks";
const String TOPPXFDR::xlclass_interdecoys = "interdecoys";
const String TOPPXFDR::xlclass_fulldecoysinterlinks = "fulldecoysinterlinks";
const String TOPPXFDR::xlclass_monodecoys = "monodecoys";
const String TOPPXFDR::xlclass_intralinks = "intralinks";
const String TOPPXFDR::xlclass_interlinks = "interlinks";
const String TOPPXFDR::xlclass_monolinks  = "monolinks";
const String TOPPXFDR::xlclass_decoys = "decoys";
const String TOPPXFDR::xlclass_hybriddecoysintralinks = "hybriddecoysintralinks";
const String TOPPXFDR::xlclass_hybriddecoysinterlinks = "hybriddecoysinterlinks";



// Parameters for actually calculating the number of FPs
const double TOPPXFDR::fpnum_score_start = 0;
const double TOPPXFDR::fpnum_score_end = 100;
const double TOPPXFDR::fpnum_score_step = 0.1;


// the actual main function needed to create an executable
int main(int argc, const char ** argv)
{
  TOPPXFDR tool;
  return tool.main(argc, argv);
}
/// @endcond
