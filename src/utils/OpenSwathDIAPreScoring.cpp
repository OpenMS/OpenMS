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
// $Maintainer: Timo Sachsenberg $
// $Authors: Witold Wolski $
// --------------------------------------------------------------------------

#include <fstream>
#include <iostream>
//#include <boost/filesystem.hpp>

#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DIAPrescoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataFrameWriter.h>

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/MzMLFile.h>

//using namespace OpenMS;
//using namespace std;


/**
  @page UTILS_OpenSwathDIAPreScoring OpenSwathDIAPreScoring

  @brief ...

  SWATH specific parameters only apply if you have full MS2 spectra maps.

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_OpenSwathDIAPreScoring.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_OpenSwathDIAPreScoring.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
using namespace OpenMS;

class DIAPreScoring :
  public TOPPBase
{
public:

  DIAPreScoring() :
    TOPPBase("OpenSwathDIAPreScoring", "Scoring spectra using the DIA scores.", false)
  {
  }

protected:

  typedef PeakMap MapType;
  typedef boost::shared_ptr<PeakMap> MapTypePtr;

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("tr", "<file>", "", "transition file");
    setValidFormats_("tr", ListUtils::create<String>("TraML"));

    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", ListUtils::create<String>("tsv"));

    registerInputFileList_("swath_files", "<files>", StringList(),
                           "Swath files that were used to extract the transitions. If present, SWATH specific scoring will be applied.",
                           false);
    setValidFormats_("swath_files", ListUtils::create<String>("mzML"));

    registerDoubleOption_("min_upper_edge_dist", "<double>", 0.0,
                          "Minimal distance to the edge to still consider a precursor, in Thomson (only in SWATH)",
                          false);


  }

  Param getSubsectionDefaults_(const String&) const override
  {
    return OpenMS::DiaPrescore().getDefaults();
  }

  ExitCodes main_(int, const char**) override
  {
    OpenMS::StringList file_list = getStringList_("swath_files");
    std::string tr_file = getStringOption_("tr");
    std::cout << tr_file << std::endl;
    //std::string out = getStringOption_("out");
    //std::cout << out << std::endl;
    double min_upper_edge_dist = getDoubleOption_("min_upper_edge_dist");

    // If we have a transformation file, trafo will transform the RT in the
    // scoring according to the model. If we dont have one, it will apply the
    // null transformation.
    Param feature_finder_param = getParam_().copy("algorithm:", true);

    // Create the output map, load the input TraML file and the chromatograms
    MapType exp;
    OpenSwath::LightTargetedExperiment transition_exp;

    std::cout << "Loading TraML file" << std::endl;
    {
      OpenMS::TargetedExperiment transition_exp_;
      TraMLFile t;
      t.load(tr_file, transition_exp_);
      //int pept =  transition_exp_.getPeptides().size();
      //int prot = transition_exp_.getProteins().size();
      //int trans = transition_exp_.getTransitions().size();
      OpenSwathDataAccessHelper::convertTargetedExp(transition_exp_, transition_exp);
      int ltrans = transition_exp.transitions.size();

      std::cout << ltrans << std::endl;
    }
    // Here we deal with SWATH files (can be multiple files)

    for (Size i = 0; i < file_list.size(); ++i)
    {
      MzMLFile swath_file;
      MapTypePtr swath_map (new MapType);
      FeatureMap featureFile;
      std::cout << "Loading file " << file_list[i] << std::endl;

      // no progress log on the console in parallel

      std::string fileout = file_list[i];

      /// Returns the basename of the file (without the path).
      /// Returns the path of the file (without the file name).

      //boost::filesystem::path x(fileout);
      //boost::filesystem::path y = x.parent_path()  ;
      //std::string fname = x.stem().string();


      //std::string tmp = File.basename(fileout);
      std::string fname = File::removeExtension(fileout);
      fname += ".tsv";


      swath_file.setLogType(log_type_);
      swath_file.load(file_list[i], *swath_map);
      if (swath_map->size() == 0 || (*swath_map)[0].getPrecursors().size() == 0)
      {
        std::cerr << "WARNING: File " << swath_map->getLoadedFilePath()
                  << " does not have any experiments or any precursors. Is it a SWATH map?"
                  << std::endl;
        continue;
      }
      // Find the transitions to extract and extract them
      OpenSwath::LightTargetedExperiment transition_exp_used;
      double upper, lower;
      const std::vector<Precursor> prec = (*swath_map)[0].getPrecursors();
      lower = prec[0].getMZ() - prec[0].getIsolationWindowLowerOffset();
      upper = prec[0].getMZ() + prec[0].getIsolationWindowUpperOffset();
      OpenSwathHelper::selectSwathTransitions(transition_exp, transition_exp_used,
                                              min_upper_edge_dist, lower, upper);
      if (transition_exp_used.getTransitions().size() == 0)
      {
        std::cerr << "WARNING: For file " << swath_map->getLoadedFilePath()
                  << " there are no transitions to extract." << std::endl;
        continue;
      }
      //OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType transition_group_map;
      std::cout << "Using Spectrum Interface!" << std::endl;
      OpenSwath::SpectrumAccessPtr  spectrumAccess = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(
        swath_map);

      //std::cout << "using data frame writer for storing data. Outfile :" << out << std::endl;
      OpenSwath::IDataFrameWriter* dfw = new OpenSwath::CSVWriter(fname);
      OpenMS::DiaPrescore dp;
      dp.operator()(spectrumAccess, transition_exp_used, dfw);
      delete dfw;
      //featureFinder.pickExperiment(chromatogram_ptr, out_featureFile,
      //transition_exp_used, trafo, swath_ptr, transition_group_map);
      //FeatureXMLFile().store(out, out_featureFile);
    }         //end of for loop
    return EXECUTION_OK;
  }       //end of _main

};

int main(int argc, const char** argv)
{
  DIAPreScoring tool;
  int code = tool.main(argc, argv);
  return code;

}

/// @endcond
