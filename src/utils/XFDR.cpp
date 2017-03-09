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
#include <OpenMS/FORMAT/CrossLinkClassesFile.h>
#include <OpenMS/FORMAT/CsvFile.h>

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


/**
 * Ensures that the order vector is sorted in descending order with respect to the Peptide Identifications (meta value: Score)
 * This is only used for assertions.
 *
 */
bool isSortedDescending(vector< size_t > & order, vector< vector < PeptideIdentification > > & spectra, Int idx)
{
  for (size_t i = 0; i < spectra.size() - 1; ++i)
  {
      double a_score = (double) spectra[order[i]][idx].getMetaValue("OpenXQuest:score");
      double b_score = (double) spectra[order[i+1]][idx].getMetaValue("OpenXQuest:score");

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

    static const String param_in_xlclasses;  // Parameter for specifying the cross-link classes
    static const String param_in_xquestxml;  // Parameter for the original xQuest XML file
    static const String param_out_csvresults; // Reports peptide identifications with annotated FDR in the CSV file format
    static const String param_minborder;  // minborder -5 # filter for minimum precursor mass error (ppm)
    static const String param_maxborder;  // maxborder  5 # filter for maximum precursor mass error (ppm)
    static const String param_mindeltas;  // mindeltas  0.95 # filter for delta score, 0 is no filter, minimum delta score required, hits are rejected if larger or equal
    static const String param_minionsmatched; // minionsmatched 0 # Filter for minimum matched ions per peptide
    static const String param_uniquexl; // calculate statistics based on unique IDs
    static const String param_qtransform; // transform simple FDR to q-FDR values
    static const String param_minscore; // minscore 0 # minimum ld-score to be considered
    static const String param_verbose; // Whether or not the output of the tool should be verbose.

    // Parameters which are not directly related to xProphet
    static const String param_fdr_calc_method;

    // Number of ranks used
    static const Int n_rank;

    // Constants related to particular xl classes
    // decoy classes
    static const String xlclass_intradecoys; // intradecoys
    static const String xlclass_fulldecoysintralinks; // fulldecoysintralinks
    static const String xlclass_interdecoys; // interdecoys
    static const String xlclass_fulldecoysinterlinks; // fulldecoysinterlinks
    static const String xlclass_monodecoys; // monodecoys
    static const String xlclass_intralinks; // intralinks
    static const String xlclass_interlinks; // interlinks
    static const String xlclass_monolinks;  // monolinks

    // Parameters to actually calculate the number of FPs // TODO Needs to more dynamic in the future
    static const double fpnum_score_start;
    static const double fpnum_score_end;
    static const double fpnum_score_step;

    TOPPXFDR() :
            TOPPBase("XFDR", "Template for Tool creation", false)
    {
    }

protected:


    /** False positve counting as performed by xProphet software package.
      *
      * @brief xprophet  method for FP counting as implemented in xProphet
      * @param cum_histograms Cumulative score distributions
      * @param fp_counts Number of FPs for each score threshold.
      */
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


     /** Target counting as performed by the xProphet software package
       *
       * @brief xprophet  method for target hits counting as implemented in xProphet
       * @param cum_histograms Cumulative score distributions
       * @param target_counts Number of target hits
       */
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

     /** Target counting as performed by the xProphet software package
      *
      * @brief xprophet  method for target hits counting as implemented in xProphet
      * @param cum_histograms Cumulative score distributions
      */
     void fdr_xprophet(std::map< String, Math::CumulativeHistogram<>  * > & cum_histograms,
                       const String  & targetclass, const String & decoyclass, const String & fulldecoyclass,
                       vector< double > & fdr, bool adjusted)
     {
       for (double current_score = TOPPXFDR::fpnum_score_start +  (TOPPXFDR::fpnum_score_step/2) ;
            current_score <= TOPPXFDR::fpnum_score_end - (TOPPXFDR::fpnum_score_step/2);
            current_score += TOPPXFDR::fpnum_score_step)
       {
           double estimated_n_decoys = cum_histograms[decoyclass]->binValue(current_score) - 2 * cum_histograms[fulldecoyclass]->binValue(current_score);
           double n_targets = cum_histograms[targetclass]->binValue(current_score);
           double summand = adjusted ? estimated_n_decoys : 0;
           fdr.push_back(n_targets > 0 ? estimated_n_decoys / (n_targets + summand) : 0);
       }
    }


    // this function will be used to register the tool parameters
    // it gets automatically called on tool execution
    void registerOptionsAndFlags_()
    {
      // Verbose Flag
      registerFlag_(TOPPXFDR::param_verbose, "Whether the log of information will be loud and noisy");

      // File input (Cross-link classes)
      registerInputFile_(TOPPXFDR::param_in_xlclasses, "<classes_file>", "",
                         "Specification of cross-link classes to compute statistics on", true, false);

      // File input (XQuest result XML)
      registerInputFile_(TOPPXFDR::param_in_xquestxml, "<file>", "", "Results in the original xquest.xml format", false);
      setValidFormats_(TOPPXFDR::param_in_xquestxml, ListUtils::create<String>("xml"));

      registerOutputFile_(TOPPXFDR::param_out_csvresults, "<csv_file>", "", "Peptide Identification annotated with the FDR values", false, false);
      setValidFormats_(TOPPXFDR::param_out_csvresults, ListUtils::create<String>("csv"));

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

      // False positve Counting  method
      registerStringOption_(TOPPXFDR::param_fdr_calc_method, "<fp_count_method>", "xprophet", "Method how the false discovery rate is computed", false, true);
      setValidStrings_(TOPPXFDR::param_fdr_calc_method, ListUtils::create<String>("xprophet,xprophet_adjusted"));

    }

    // the main_ function is called after all parameters are read
    ExitCodes main_(int, const char **)
    {
      bool arg_verbose = getFlag_(TOPPXFDR::param_verbose);

      //-------------------------------------------------------------
      // Handling the cross-link classes specification
      //-------------------------------------------------------------
      String arg_in_xlclasses = getStringOption_(TOPPXFDR::param_in_xlclasses);

      // Container for rules of specifying classes
      CrossLinkClassesFile cross_link_classes_file;

      if (! cross_link_classes_file.load(arg_in_xlclasses))
      {
         LOG_ERROR << "ERROR: Reading of cross-link class specification file failed." << endl;
         return PARSE_ERROR;
      }

      //-------------------------------------------------------------
      // Determine FDR calculation method
      //-------------------------------------------------------------
      String arg_fdr_calc_method = getStringOption_(TOPPXFDR::param_fdr_calc_method);
      if (arg_fdr_calc_method == "xprophet" || arg_fdr_calc_method == "xprophet_adjusted")
      {
         StringList required_classes = ListUtils::create<String>("intralinks,interlinks,monolinks,intradecoys,fulldecoysintralinks,interdecoys,fulldecoysinterlinks,monodecoys");
         for (StringList::const_iterator required_classes_it = required_classes.begin();
              required_classes_it != required_classes.end(); ++required_classes_it)
         {
           String classname = *required_classes_it;
           if ( ! cross_link_classes_file.has((classname)))
           {
             LOG_ERROR << "ERROR: xProphet target counting selected, but the following xlink class has not been defined: " << classname << endl;
             return ILLEGAL_PARAMETERS;
           }
         }
      }
      else
      {
        LOG_ERROR << "ERROR: Unsupported method for FDR calculation. Aborting" << endl;
        return ILLEGAL_PARAMETERS;
      }

      if (arg_verbose)
      {
        for(CrossLinkClassesFile::ClassNameConstIterator it = cross_link_classes_file.begin(); it != cross_link_classes_file.end(); ++it)
        {
              cout << "Class defined: " <<  *it << '\n';
        }
        cout << "----------------------------" << endl;
      }

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
      // Parsing XML file of xQuest
      //-------------------------------------------------------------
      String arg_in_xquestxml = getStringOption_(TOPPXFDR::param_in_xquestxml);

      if (arg_verbose)
      {
        LOG_INFO << "INFO: Parsing xQuest input XML file: " << arg_in_xquestxml << endl;
      }
      // Core data structures for the util
      vector< XQuestResultMeta > metas;
      vector< vector < PeptideIdentification > > spectra;

      size_t  pep_id_index = TOPPXFDR::n_rank - 1; // This is of course trash if more than 1 ranks are used.

      // Parse XQuestResultXMLFile (TODO Also support idXML and mzID)
      XQuestResultXMLFile xquest_result_file;
      xquest_result_file.load(arg_in_xquestxml, metas, spectra, false, 1); // We do not load 'empty' spectra here
      size_t n_spectra = spectra.size();

      if (arg_verbose)
      {
        LOG_INFO << "Total number of spectra: " << n_spectra << "\n"
                 << "Total number of hits: "    << xquest_result_file.get_n_hits() << endl;
      }

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
      //-------------------------------------------------------------
      // The score is calculated for each hit h on the set of all hits of the spectrum that encompasses
      // TODO Maybe it makes sense to wrap the vectors into a struct (they belong together by rank)
      std::vector< std::vector< double >* > delta_scores(n_spectra);
      std::vector< size_t > n_min_ions_matched(n_spectra);

      for(size_t i = 0; i < delta_scores.size(); ++i)
      {
          size_t n_hits = spectra[i].size();
          delta_scores[i] = new std::vector<double>(n_hits);
          vector<double> * current = delta_scores[i];

          assert(n_hits > 0); // because we initially do not load 'empty' spectra
          // calculate n_min_ions_matched
          PeptideIdentification * pep_id1 = &spectra[i][0];
          assert((int) pep_id1->getMetaValue("xl_rank") == 1); // because hits are sorted according to their rank within the spectrum

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
      typedef std::vector<size_t> ranks;

      ranks order_score ( boost::counting_iterator<size_t>(0),
                          boost::counting_iterator<size_t>(n_spectra));

      // Configure the sorting of the Peptide Identifications
      less_than_by_key order_conf = {
        spectra,            // elements
        pep_id_index,       // idx
        "OpenXQuest:score"  // key
      };

      std::sort(order_score.begin(), order_score.end(), order_conf);
      assert(isSortedDescending(order_score, spectra, pep_id_index));

      // For unique IDs
      std::set<String> unique_ids;

      // Applies user specified filters and aggregates the scores for the corresponding classes
      // TODO Maybe put the stuff on the heap? Stack is getting full now.
      std::map< String, vector< double > > scores;

      for (size_t i = 0; i < n_spectra; ++i)
      {
       // Extract required attributes of the peptide_identification (filter criteria)
       PeptideIdentification * pep_id = &spectra[order_score[i]][pep_id_index];
       double error_rel = (double) pep_id->getMetaValue("OpenXQuest:error_rel");
       double delta_score = (*(delta_scores[order_score[i]]))[pep_id_index];    // Index 0 because we only consider rank 1 here
       size_t ions_matched = n_min_ions_matched[order_score[i]];
       double score = (double) pep_id->getMetaValue("OpenXQuest:score");
       String id = (String) pep_id->getMetaValue("OpenXQuest:id");

       // Only consider peptide identifications which  fullfill all filter criteria specified by the user
       if (     arg_minborder <= error_rel
            &&  arg_maxborder >= error_rel
            &&  (mindelta_filter_disabled || delta_score < arg_mindeltas)
            &&  ions_matched  >= arg_minionsmatched
            &&  score         >= arg_minscore
            &&  ( ! arg_uniquex || unique_ids.find(id) == unique_ids.end()))
       {
           unique_ids.insert(id);
           cross_link_classes_file.collect(*pep_id, scores, "OpenXQuest:score");
       }
      }

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


      bool adjusted;
      if (arg_fdr_calc_method == "xprophet")
      {
          adjusted = false;
      }
      else if (arg_fdr_calc_method == "xprophet_adjusted")
      {
        adjusted = true;
      }
      else
      {
          LOG_ERROR << "ERROR: Unsupported FDR calculation method. Aborting." << endl;
          return ILLEGAL_PARAMETERS;
      }

      vector< double > fdr_interlinks;
      this->fdr_xprophet(cum_histograms, TOPPXFDR::xlclass_interlinks, TOPPXFDR::xlclass_interdecoys, TOPPXFDR::xlclass_fulldecoysinterlinks, fdr_interlinks, adjusted);
      vector< double > fdr_intralinks;
      this->fdr_xprophet(cum_histograms, TOPPXFDR::xlclass_intralinks, TOPPXFDR::xlclass_intradecoys, TOPPXFDR::xlclass_fulldecoysintralinks, fdr_intralinks, adjusted);

      /*
      for (vector< double>::const_iterator it = fdr_intralinks.begin(); it != fdr_intralinks.end(); ++it)
      {
         cout << *it << endl;
      }
      */

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


      // Write the CSV output file if requested
      String arg_out_csvresults = getStringOption_(TOPPXFDR::param_out_csvresults);
      if ( ! arg_out_csvresults.empty())
      {
          CsvFile csv_file;

      }
      return EXECUTION_OK;
    }
};
const String TOPPXFDR::param_in_xlclasses = "in_xlclasses";
const String TOPPXFDR::param_in_xquestxml = "in_xquestxml";
const String TOPPXFDR::param_out_csvresults = "out_csvresults";
const String TOPPXFDR::param_minborder = "minborder";
const String TOPPXFDR::param_maxborder = "maxborder";
const String TOPPXFDR::param_mindeltas = "mindeltas";
const String TOPPXFDR::param_minionsmatched = "minionsmatched";
const String TOPPXFDR::param_uniquexl = "uniquexl";
const String TOPPXFDR::param_qtransform = "qtransform";
const String TOPPXFDR::param_minscore = "minscore";
const String TOPPXFDR::param_verbose = "verbose";
const String TOPPXFDR::param_fdr_calc_method = "fdr_calc_method";

const Int    TOPPXFDR::n_rank = 1; //  Number of ranks used

const String TOPPXFDR::xlclass_intradecoys = "intradecoys";
const String TOPPXFDR::xlclass_fulldecoysintralinks = "fulldecoysintralinks";
const String TOPPXFDR::xlclass_interdecoys = "interdecoys";
const String TOPPXFDR::xlclass_fulldecoysinterlinks = "fulldecoysinterlinks";
const String TOPPXFDR::xlclass_monodecoys = "monodecoys";
const String TOPPXFDR::xlclass_intralinks = "intralinks";
const String TOPPXFDR::xlclass_interlinks = "interlinks";
const String TOPPXFDR::xlclass_monolinks  = "monolinks";


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
