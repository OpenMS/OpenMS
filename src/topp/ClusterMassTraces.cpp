// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MasstraceCorrelator.h>
#include <OpenMS/FORMAT/FileHandler.h>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**

@page TOPP_ClusterMassTraces ClusterMassTraces

@brief Cluster mass traces occurring in the same map together

Cluster mass traces together found in a mass spectrometric map (MS1 or MS2).
Input is a consensus map containing individual mass traces, the output may be
spectra containing all clustered features.

Mass traces are clustered independent of precursor traces in another map
(this is the more simple approach)  and pseudo spectra are created without
any precursors assigned. This is useful for 

 - clustering of features in an MS1 map (isotope traces, charge states etc)
 - clustering of features in an SWATH map (fragment ions from the same precursor, isotope traces, charge states etc)

On the clustered fragments in an MS2 map, one can then (optionally) do 

 - de novo searches 
 - calculate the most likely precursor(s) and DB-search

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_ClusterMassTraces.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_ClusterMassTraces.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


#include <OpenMS/APPLICATIONS/TOPPBase.h>

class TOPPClusterMassTraces
  : public TOPPBase, 
    public ProgressLogger

{

  // Docu
  //

 public:

  TOPPClusterMassTraces()
    : TOPPBase("ClusterMassTraces","Creates pseudo spectra.")
  {
  }

 protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in","<file>","","Mass traces");
    setValidFormats_("in",ListUtils::create<String>("consensusXML"));

    registerOutputFile_("out","<file>","","output file");
    setValidFormats_("out",ListUtils::create<String>("mzML"));

    registerDoubleOption_("min_pearson_correlation", "<double>", 0.7, "Minimal pearson correlation score", false);
    registerIntOption_("min_peak_nr", "<number>", 1, "Minimal peak nr to output pseudo spectra", false);
    registerIntOption_("max_lag", "<number>", 1, "Maximal lag", false);
    registerDoubleOption_("max_rt_apex_difference", "<double>", 5.0, "Maximal difference of the apex in retention time", false);
    registerDoubleOption_("max_intensity_cutoff", "<double>", 0.0, "Maximal intensity to be added to a spectrum", false);

    registerDoubleOption_("add_precursor", "<double>", 0.0, "Add a precursor mass", false);
  }

 public:

  ExitCodes main_(int , const char**) override
  {

    setLogType(log_type_); 

    String infile = getStringOption_("in");
    String out = getStringOption_("out");

    double min_pearson_correlation_ = getDoubleOption_("min_pearson_correlation");
    int max_lag_ = getIntOption_("max_lag");
    int min_peak_nr = getIntOption_("min_peak_nr");
    double max_rt_apex_difference_ = getDoubleOption_("max_rt_apex_difference");
    double add_precursor = getDoubleOption_("add_precursor");
    // double max_intensity_cutoff_ = getDoubleOption_("max_intensity_cutoff");

    ConsensusMap masstrace_map;
    FileHandler().loadConsensusFeatures(infile, masstrace_map, {FileTypes::CONSENSUSXML}, log_type_);

    MSExperiment pseudo_spectra;

    if (masstrace_map.empty())
    {
      // Error
    }

    std::cout << "Input map " << infile <<" has size: " << masstrace_map.size() << std::endl;

    masstrace_map.sortByIntensity(true);

    std::cout << "Input map " << infile <<" has size: " << masstrace_map.size() << std::endl;

    OpenMS::MasstraceCorrelator mtcorr;
    mtcorr.setLogType(log_type_); 
    mtcorr.createPseudoSpectra(masstrace_map, pseudo_spectra, min_peak_nr,
        min_pearson_correlation_, max_lag_, max_rt_apex_difference_/* , max_intensity_cutoff_ */);
    pseudo_spectra.sortSpectra();

    // If we want to set a specific precursor, do this now
    if (add_precursor > 0 )
    {
      for (Size i = 0; i < pseudo_spectra.size(); i++)
      {
        Precursor p;
        //p.setIsolationWindowLowerOffset(swath_lower);
        //p.setIsolationWindowUpperOffset(swath_upper);
        p.setMZ(add_precursor);
        std::vector<Precursor> preclist;
        preclist.push_back(p);
        pseudo_spectra[i].setPrecursors(preclist);
      }
    }
    FileHandler().storeExperiment(out,pseudo_spectra, {FileTypes::MZML});

    return EXECUTION_OK;
  }

};

int main( int argc, const char** argv )
{

  TOPPClusterMassTraces tool;
  return tool.main(argc,argv);
}

///@endcond
