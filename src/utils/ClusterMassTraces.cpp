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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MasstraceCorrelator.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_ClusterMassTraces ClusterMassTraces

  @brief Cluster mass traces occuring in the same map together

  Cluster mass traces together found in a mass spectrometric map (MS1 or MS2).
  Input is a consensus map containing individual mass traces, the output may be
  spectra containing all clustered features.

  Mass traces are clustered independend of precursor traces in another map
  (this is the more simple approach)  and pseudo spectra are created without
  any precursors assigned. This is useful for 

   - clustering of features in an MS1 map (isotope traces, charge states etc)
   - clustering of features in an SWATH map (fragment ions from the same precursor, isotope traces, charge states etc)

  On the clustered fragments in an MS2 map, one can then (optionally) do 

   - de novo searches 
   - calculate the most likely precursor(s) and DB-search

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
    : TOPPBase("ClusterMassTraces","Creates pseudo spectra.", false)
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

    ConsensusXMLFile consensus_f;
    consensus_f.setLogType(log_type_);
    ConsensusMap masstrace_map;
    consensus_f.load(infile, masstrace_map);

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
    MzMLFile().store(out,pseudo_spectra);

    return EXECUTION_OK;
  }

};

int main( int argc, const char** argv )
{

  TOPPClusterMassTraces tool;
  return tool.main(argc,argv);
}

