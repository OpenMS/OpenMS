// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MasstraceCorrelator.h>
#include <OpenMS/FORMAT/FileHandler.h>

#ifdef TESTING
#define DEBUG_MASSTRACES
#endif


//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_ClusterMassTracesByPrecursor ClusterMassTracesByPrecursor

@brief Identifies precursor mass traces and tries to correlate them with fragment ion mass traces in SWATH maps.

This algorithm will try to correlate the masstraces to find co-eluting traces and cluster them.

This program looks at mass traces in a precursor MS1 map and tries to
correlate them with features found in the corresponding MS2 map based on
their elution profile. It uses

 - the mass traces from the MS1 in consensusXML format [note this is an unintended use of the consesusXML format to also store intensities]
 - the mass traces from the MS2 (SWATH map)

 It does a separate correlation analysis on the MS1 and the MS2 map,
 both produces a set of pseudo spectra.
 In a second (optional) step, the MS2 pseudo spectra are correlated with
 the MS1 traces and the most likely precursor is assigned to the pseudo
 spectrum.
  
It is based on the following papers:
ETISEQ -- an algorithm for automated elution time ion sequencing of concurrently fragmented peptides for mass spectrometry-based proteomics
  BMC Bioinformatics 2009, 10:244 doi:10.1186/1471-2105-10-244 ; http://www.biomedcentral.com/1471-2105/10/244
  they use FFT to correlate and then use lag of at least 1 scan and pearson correlation of 0.7 to assign precursors to product ions
  If one fragment matches to multiple precursors, it is assigned to all of them. If it doesn't match any, it is assigned to all

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_ClusterMassTracesByPrecursor.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_ClusterMassTracesByPrecursor.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace std;
using namespace OpenMS;

class TOPPCorrelateMasstraces
  : public TOPPBase, 
    public ProgressLogger
{

 public:

  TOPPCorrelateMasstraces()
    : TOPPBase("ClusterMassTracesByPrecursor", "Correlate precursor masstraces with fragment ion masstraces in SWATH maps based on their elution profile.")
  {
  }

 protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in_ms1","<file>","","MS1 mass traces");
    setValidFormats_("in_ms1",ListUtils::create<String>("consensusXML"));

    registerInputFile_("in_swath","<file>","","MS2 / SWATH mass traces");
    setValidFormats_("in_swath",ListUtils::create<String>("consensusXML"));

    registerOutputFile_("out","<file>","","output file");
    setValidFormats_("out",ListUtils::create<String>("mzML"));

    // registerFlag_("ms1_centric","MS1 centric - find MS1 features first and then add MS2s (MSE like)");
    registerFlag_("assign_unassigned_to_all","Assign unassigned MS2 fragments to all precursors (only for ms1_centrif)");

    registerDoubleOption_("min_pearson_correlation", "<double>", 0.7, "Minimal pearson correlation score to match elution profiles to each other.", false); // try 0.3, 0.5 and 0.7
    registerIntOption_("max_lag", "<number>", 1, "Maximal lag (e.g. by how many spectra the peak may be shifted at most). This parameter will depend on your chromatographic setup but a number between 1 and 3 is usually sensible.", false);
    registerIntOption_("min_nr_ions", "<number>", 3, "Minimal number of ions to report a spectrum.", false);
    registerDoubleOption_("max_rt_apex_difference", "<double>", 5.0, "Maximal difference of the apex in retention time (in seconds). This is a hard parameter, all profiles further away will not be considered at all.", false);

    registerDoubleOption_("swath_lower", "<double>", 0.0, "Swath lower isolation window", false);
    registerDoubleOption_("swath_upper", "<double>", 0.0, "Swath upper isolation window", false);
  }

 public:

  ExitCodes main_(int , const char**) override
  {
    setLogType(log_type_); 

    String ms1 = getStringOption_("in_ms1");
    String in_swath =  getStringOption_("in_swath");
    String out = getStringOption_("out");

    // bool ms1_centric = getFlag_("ms1_centric");

    double swath_lower = getDoubleOption_("swath_lower");
    double swath_upper = getDoubleOption_("swath_upper");

    // Load input:
    // - MS1 feature map containing the MS1 mass traces
    // - MS2 feature map containing the MS2 (SWATH) mass traces
    ConsensusMap MS1_feature_map;
    ConsensusMap MS2_feature_map;
    FileHandler().loadConsensusFeatures(ms1, MS1_feature_map, {FileTypes::CONSENSUSXML}, log_type_);
    FileHandler().loadConsensusFeatures(in_swath, MS2_feature_map, {FileTypes::CONSENSUSXML}, log_type_);
    cout << "Loaded consensus maps" << endl;

#ifdef DEBUG_MASSTRACES
    for (Size i=0; i<MS1_feature_map.size(); ++i)
    {
        ConsensusFeature f1 = MS1_feature_map[i];
        cout << "MS1 mass trace " << i << " at " << f1.getMZ() << " and " << 
          f1.getRT() <<  " +/- " << f1.getWidth() << " with " << f1.getIntensity() << endl;
    }
#endif

    MSExperiment pseudo_spectra_ms1centric;
    MS1CentricClustering(MS1_feature_map, MS2_feature_map, 
        swath_lower, swath_upper, pseudo_spectra_ms1centric);
    FileHandler().storeExperiment(out,pseudo_spectra_ms1centric, {FileTypes::MZML}, log_type_);

    return EXECUTION_OK;
  }

  /** @brief Cluster fragments ions with their corresponding precursors 
   *
   * This is based on the ETISEQ algorithm and works as follows:
   *
   *  - Identify the precursor traces
   *  - For each precursor determine which are the most likely fragments and
   *    then assign those to the precursor
   *  - Assign unassigned fragments to scans
   *  - Create actual precursor spectra
   *
   * TODO: incorporate elements from DIAUmpire
   *  - allow ions to be assigned to multiple precursors
   *  - also generate mass traces from the unfragmented precursors
   *
  */
  void MS1CentricClustering(ConsensusMap& MS1_feature_map, ConsensusMap& MS2_feature_map, 
      double swath_lower, double swath_upper, 
      MSExperiment& pseudo_spectra_precursors1)
  {
    // -----------------------------------
    // Parameters 
    // -----------------------------------
    double min_pscore = getDoubleOption_("min_pearson_correlation");
    int max_lag = getIntOption_("max_lag");
    double rt_max_distance = getDoubleOption_("max_rt_apex_difference");
    Size min_nr_ions = (Size)getIntOption_("min_nr_ions");
    bool unassigned = getFlag_("assign_unassigned_to_all");
    // to consider all signals within 2 seconds equal makes sense with
    // 3.2 seconds between each recording => each swath will be within
    // +/- 2.0 seconds of a full scan
    double mindiff = 2.0;

    OpenMS::MasstraceCorrelator mtcorr;
    std::map< int, std::vector< std::vector<double> > >  feature_attributes; // temporary array storing the attributes for all the features
    std::vector<bool> ms2feature_used;
    std::map< int, std::vector<int> > ms1_assignment_map; // map MS1 feature ids to MS2 feature ids
    ms2feature_used.resize(MS2_feature_map.size());

    // -----------------------------------
    // Cache datastructures
    // -----------------------------------
    // We cache the RT and intensities of each feature
    std::vector< MasstraceCorrelator::MasstracePointsType > feature_points_ms2;
    std::vector< std::pair<double,double> > max_intensities_ms2; 
    std::vector< double > rt_cache_ms2;
    mtcorr.createConsensusMapCache(MS2_feature_map, feature_points_ms2, max_intensities_ms2, rt_cache_ms2);

    std::vector< MasstraceCorrelator::MasstracePointsType > feature_points_ms1;
    std::vector< std::pair<double,double> > max_intensities_ms1; 
    std::vector< double > rt_cache_ms1;
    mtcorr.createConsensusMapCache(MS1_feature_map, feature_points_ms1, max_intensities_ms1, rt_cache_ms1);

    // cache the m/z of each MS1 feature
    std::vector< double > mz_cache_ms1;
    for (Size i = 0; i < MS1_feature_map.size(); ++i)
    {
      mz_cache_ms1.push_back(MS1_feature_map[i].getMZ());
    }

    double* rt_cache_ptr;
    double current_rt;

    // -----------------------------------
    // Step 1 - assign fragment mass traces to precursors
    //
    // Go through all precursors and find suitable MS2 signals which could
    // potentially belong to this precursor.
    //
    startProgress(0, MS1_feature_map.size(), "assigning precursor to fragment ions");
    for (Size i=0; i<MS1_feature_map.size(); ++i)
    {
      setProgress(i);
      if (mz_cache_ms1[i] < swath_lower || mz_cache_ms1[i] > swath_upper) continue;
      ms1_assignment_map[i].clear();

      // Identify a given precursor and get its RT (current_rt) 
      // 
      // Obtain a pointer to the beginning of the RT vector of all MS2 features
      // (and decrement by one since in the loop we first increment the ptr)
      current_rt = rt_cache_ms1[i];
      rt_cache_ptr = &rt_cache_ms2[0];
      --rt_cache_ptr;

      for (Size j=0; j<MS2_feature_map.size(); ++j)
      {
        ++rt_cache_ptr;

        // First check whether this feature is within a suitable RT distance
        // and that is not already used.
        // Check whether the feature is already used
        //  TODO : this implies we can assign only one feature to one
        //         precursor, we might have to change that! See DIA Umpire!
        if (fabs(current_rt - (*rt_cache_ptr) ) > rt_max_distance ) continue;
        if (ms2feature_used[j]) continue;

#ifdef DEBUG_MASSTRACES
        for (Size kk=0; kk<f1_points.size(); kk++)
        { 
          cout << f1_points[kk].first << " f/s " << f1_points[kk].second << endl; 
        }
        cout << " above prec, below frag " << endl;
        for (Size kk=0; kk<f2_points.size(); kk++)
        { 
          cout << f2_points[kk].first << " f/s " << f2_points[kk].second << endl; 
        }
#endif

        // Score the MS1 mass trace against the MS2 mass trace
        int lag; double lag_intensity; double pearson_score;
        mtcorr.scoreHullpoints(feature_points_ms1[i], feature_points_ms2[j], 
            lag, lag_intensity, pearson_score, min_pscore, max_lag, mindiff);

        if (pearson_score > min_pscore && lag >= -max_lag && lag <= max_lag)
        {
#ifdef DEBUG_MASSTRACES
          cout <<  "assign fragment to precursor! " << f1.getMZ() << " -> " << f2.getMZ() << 
            " [scores " <<  lag << " " << pearson_score << "]" << endl;
#endif
          ms2feature_used[j] = true;
          ms1_assignment_map[i].push_back(j);
          std::vector< double > feature_arr;
          feature_arr.push_back(rt_cache_ms2[j]);  // MS2 retention time
          feature_arr.push_back(fabs(rt_cache_ms1[i] - rt_cache_ms2[j] ) ); // difference between MS1 and MS2 RT
          feature_arr.push_back(lag); // lag
          feature_arr.push_back(pearson_score); // pearson score
          feature_arr.push_back(lag_intensity); // lag intensity
          feature_attributes[i].push_back(feature_arr);
        }
      }

      // only keep those assignments which have enough ions
      if (ms1_assignment_map[i].size() <= min_nr_ions) 
      {
        ms1_assignment_map[i].clear();
      }

#ifdef DEBUG_MASSTRACES
      if (ms1_assignment_map[i].size() > 1)
      {
        cout << i << " idx " << " " << f1 << " size " <<  MS1_feature_map[i].size() << endl;
        cout << " to precursor " << i << " i assigned " << ms1_assignment_map[i].size() << " points" << endl;
      }
      cout << "MS1 mass trace " << i << " at " << f1.getMZ() << " and " << f1.getRT( ) << " with " << f1.getIntensity() << endl;
#endif

    }
    endProgress();

    // Stats
    Size cnt_ms2_used = 0;
    Size cnt_ms1_used = 0;
    for (Size i = 0; i < MS1_feature_map.size(); i++) 
    {
      if (!ms1_assignment_map[i].empty()) cnt_ms1_used++;
    }
    for (Size i = 0; i < ms2feature_used.size(); i++) {
      if (ms2feature_used[i]) cnt_ms2_used++;
    }

    std::cout <<"I have assigned " << cnt_ms2_used << " (out of " << MS2_feature_map.size() << 
      ") MS2 features to " << cnt_ms1_used << " (out of " << MS1_feature_map.size() << ") MS1 features " << std::endl;

    // -----------------------------------
    // Step 2 - assign the unused fragment ions (if requested)
    //
    // TODO : 
    // i) just assign them to all potentially matching spectra
    // ii) assign a fragment ion only to a single precursor
    int cnt = 0;
    startProgress(0, MS2_feature_map.size(), "assigning the unused fragments ");
    for (Size j=0; j<MS2_feature_map.size() && unassigned; ++j)
    {
      setProgress(j);
      if (ms2feature_used[j]) continue;
      cnt++;

      // find suitable MS1 spectra to assign these
      for (Size i=0; i<MS1_feature_map.size(); ++i)
      {
        if (mz_cache_ms1[i] < swath_lower || mz_cache_ms1[i] > swath_upper ) continue;
        if (ms1_assignment_map[i].empty()) continue;
        if (fabs(rt_cache_ms1[i] - rt_cache_ms2[j]) > rt_max_distance) continue;

        // Assign to all matching MS1 precursors
        ms1_assignment_map[i].push_back(j);
      }
    }
    endProgress();
    cout << "There were " << cnt << " (out of " << MS2_feature_map.size() << " ) unused fragment ions that were assigned to all spectra within RT range." << endl;

    // -----------------------------------
    // Step 3 - create spectra and assign precursor and fragments to spectra
    cnt = 0;
    startProgress(0, MS1_feature_map.size(), "create the spectra and assign the fragments ");
    for (Size i=0; i<MS1_feature_map.size(); ++i)
    {
      setProgress(i);
      if (mz_cache_ms1[i] < swath_lower || mz_cache_ms1[i] > swath_upper) continue;

      MSSpectrum spectrum;
      ConsensusFeature f2 = MS1_feature_map[i];
      spectrum.setRT(f2.getRT());
      spectrum.setMSLevel(2);
      Precursor p;
      p.setMZ(f2.getMZ());
      std::vector<Precursor> preclist;
      preclist.push_back(p);
      spectrum.setPrecursors(preclist);

      // fill meta data
      spectrum.getFloatDataArrays().clear();
      spectrum.getFloatDataArrays().resize(5);
      spectrum.getFloatDataArrays()[0].setName("RT_apex");
      spectrum.getFloatDataArrays()[1].setName("RT_diff");
      spectrum.getFloatDataArrays()[2].setName("lag");
      spectrum.getFloatDataArrays()[3].setName("pearson_score");
      spectrum.getFloatDataArrays()[4].setName("lag_intensity");
      int j = 0;
      for (std::vector<int>::iterator it = ms1_assignment_map[i].begin(); it != ms1_assignment_map[i].end(); ++it)
      {
        ConsensusFeature f1 = MS2_feature_map[*it];
        Peak1D peak;
        peak.setMZ(f1.getMZ());
        peak.setIntensity(f1.getIntensity());
        spectrum.push_back(peak);

        spectrum.getFloatDataArrays()[0].push_back(feature_attributes[i][j][0]);
        spectrum.getFloatDataArrays()[1].push_back(feature_attributes[i][j][1]);
        spectrum.getFloatDataArrays()[2].push_back(feature_attributes[i][j][2]);
        spectrum.getFloatDataArrays()[3].push_back(feature_attributes[i][j][3]);
        spectrum.getFloatDataArrays()[4].push_back(feature_attributes[i][j][4]);
        j++;
      }
 
      if (spectrum.size() > min_nr_ions) 
      {
        pseudo_spectra_precursors1.addSpectrum(spectrum);
        cnt++;
#ifdef DEBUG_MASSTRACES
        cout << "MS1 mass trace " << i << " was assigned " << ms1_assignment_map[i].size() << " " << f2.getRT() << " " << f2.getMZ() << " " << f2.getIntensity() << endl;
#endif
      }
    }
    endProgress();
    cout << "There were " << cnt << " precursor ions with more than " << min_nr_ions << " fragment ion assigned." << endl;
  }

};

int main( int argc, const char** argv )
{
  TOPPCorrelateMasstraces tool;
  return tool.main(argc,argv);
}

///@endcond
