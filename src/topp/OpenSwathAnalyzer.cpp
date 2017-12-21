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

#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <fstream>
#include <boost/shared_ptr.hpp>

using namespace OpenMS;
using namespace std;

//
//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
 @page TOPP_OpenSwathAnalyzer OpenSwathAnalyzer

 @brief  Executes a peak-picking and scoring algorithm on MRM/SRM data.

    <CENTER>
        <table>
            <tr>
                <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
                <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ OpenSwathAnalyzer \f$ \longrightarrow \f$</td>
                <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
            </tr>
            <tr>
                <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_OpenSwathChromatogramExtractor </td>
                <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_OpenSwathFeatureXMLToTSV </td>
            </tr>
            <tr>
                <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MRMMapper </td>
                <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_OpenSwathConfidenceScoring </td>
            </tr>
        </table>
    </CENTER>

 The idea of the OpenSwath Analyzer is to analyze a series of chromatograms
 together with the associated meta information (stored in TraML format) in
 order to determine likely places of elution of a peptide in targeted
 proteomics data (derived from SWATH-MS or MRM/SRM).

 <B>The command line parameters of this tool are:</B>
 @verbinclude TOPP_OpenSwathAnalyzer.cli

 <B>The algorithm parameters for the Analyzer filter are:</B>
 @htmlinclude TOPP_OpenSwathAnalyzer.html

 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPOpenSwathAnalyzer : public TOPPBase
{
public:

  TOPPOpenSwathAnalyzer() :
  TOPPBase("OpenSwathAnalyzer",
           "Picks peaks and finds features in an SWATH-MS or SRM experiment.", true)
  {
  }

protected:

  typedef PeakMap MapType;

  void registerModelOptions_(const String &default_model)
  {
    registerTOPPSubsection_("model", "Options to control the modeling of retention time transformations from data");
    registerStringOption_("model:type", "<name>", default_model, "Type of model", false, true);
    StringList model_types;
    TransformationDescription::getModelTypes(model_types);
    if (!ListUtils::contains(model_types, default_model))
    {
      model_types.insert(model_types.begin(), default_model);
    }
    setValidStrings_("model:type", model_types);
    registerFlag_("model:symmetric_regression", "Only for 'linear' model: Perform linear regression on 'y - x' vs. 'y + x', instead of on 'y' vs. 'x'.", true);
  }

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "",
                       "input file containing the chromatograms." /* , false */);
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFile_("tr", "<file>", "", "transition file");
    setValidFormats_("tr", ListUtils::create<String>("TraML"));

    registerInputFile_("rt_norm", "<file>", "",
                       "RT normalization file (how to map the RTs of this run to the ones stored in the library)",
                       false);
    setValidFormats_("rt_norm", ListUtils::create<String>("trafoXML"));

    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", ListUtils::create<String>("featureXML"));

    registerFlag_("no-strict",
                  "run in non-strict mode and allow some chromatograms to not be mapped.");

    addEmptyLine_();
    registerInputFileList_("swath_files", "<files>", StringList(),
                           "[applies only if you have full MS2 spectra maps] "
                           "Swath files that were used to extract the transitions. "
                           "If present, SWATH specific scoring will be used.",
                           false);
    setValidFormats_("swath_files", ListUtils::create<String>("mzML"));

    registerDoubleOption_("min_upper_edge_dist", "<double>", 0.0,
                          "[applies only if you have full MS2 spectra maps] "
                          "Minimal distance to the edge to still consider a precursor, in Thomson (only in SWATH)",
                          false);

    registerModelOptions_("linear");

    registerSubsection_("algorithm", "Algorithm parameters section");

  }

  Param getSubsectionDefaults_(const String &) const override
  {
    return MRMFeatureFinderScoring().getDefaults();
  }

  ExitCodes main_(int, const char **) override
  {

    StringList file_list = getStringList_("swath_files");
    String in = getStringOption_("in");
    String tr_file = getStringOption_("tr");
    String out = getStringOption_("out");
    double min_upper_edge_dist = getDoubleOption_("min_upper_edge_dist");
    bool nostrict = getFlag_("no-strict");

    // If we have a transformation file, trafo will transform the RT in the
    // scoring according to the model. If we dont have one, it will apply the
    // null transformation.
    String trafo_in = getStringOption_("rt_norm");
    TransformationDescription trafo;
    if (trafo_in.size() > 0)
    {
      TransformationXMLFile trafoxml;
      String model_type = getStringOption_("model:type");
      Param model_params = getParam_().copy("model:", true);
      trafoxml.load(trafo_in, trafo);
      trafo.fitModel(model_type, model_params);
    }

    Param feature_finder_param = getParam_().copy("algorithm:", true);

    // Create the output map, load the input TraML file and the chromatograms
    boost::shared_ptr<MapType> exp (new MapType());
    FeatureMap out_featureFile;
    OpenSwath::LightTargetedExperiment transition_exp;

    std::cout << "Loading TraML file" << std::endl;
    {
      TargetedExperiment *transition_exp__ = new TargetedExperiment();
      TargetedExperiment &transition_exp_ = *transition_exp__;
      {
        TraMLFile *t = new TraMLFile;
        t->load(tr_file, transition_exp_);
        delete t;
      }
      OpenSwathDataAccessHelper::convertTargetedExp(transition_exp_, transition_exp);
      delete transition_exp__;
    }

    MzMLFile mzmlfile;
    mzmlfile.setLogType(log_type_);
    mzmlfile.load(in, *exp.get());

    // If there are no SWATH files, it's just regular SRM/MRM Scoring
    if (file_list.size() == 0)
    {
      MRMFeatureFinderScoring featureFinder;
      featureFinder.setParameters(feature_finder_param);
      featureFinder.setLogType(log_type_);
      featureFinder.setStrictFlag(!nostrict);
      OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType transition_group_map;
      OpenSwath::SpectrumAccessPtr chromatogram_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
      std::vector< OpenSwath::SwathMap > empty_maps;
      featureFinder.pickExperiment(chromatogram_ptr, out_featureFile,
                                   transition_exp, trafo, empty_maps, transition_group_map);
      out_featureFile.ensureUniqueId();
      addDataProcessing_(out_featureFile, getProcessingInfo_(DataProcessing::QUANTITATION));
      FeatureXMLFile().store(out, out_featureFile);
      return EXECUTION_OK;
    }

    // Here we deal with SWATH files (can be multiple files)
    // Only in OpenMP 3.0 are unsigned loop variables allowed
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize i = 0; i < boost::numeric_cast<SignedSize>(file_list.size()); ++i)
    {
      MRMFeatureFinderScoring featureFinder;
      MzMLFile swath_file;
      boost::shared_ptr<MapType> swath_map (new MapType());
      FeatureMap featureFile;
      cout << "Loading file " << file_list[i] << endl;

////#ifndef _OPENMP
      // no progress log on the console in parallel
      swath_file.setLogType(log_type_);
      featureFinder.setLogType(log_type_);
//#endif

      swath_file.load(file_list[i], *swath_map.get());

      // Logging and output to the console
#ifdef _OPENMP
#pragma omp critical (featureFinder)
#endif
      {
        cout << "Doing file " << file_list[i]
#ifdef _OPENMP
             << " (" << i << " out of " << file_list.size() / omp_get_num_threads() << " -- total for all threads: " << file_list.size() << ")" << endl;
#else
        << " (" << i << " out of " << file_list.size() << ")" << endl;
#endif
      }

      OpenSwath::LightTargetedExperiment transition_exp_used;
      bool do_continue = OpenSwathHelper::checkSwathMapAndSelectTransitions(*swath_map.get(), transition_exp, transition_exp_used, min_upper_edge_dist);

      if (do_continue)
      {
        featureFinder.setParameters(feature_finder_param);
        featureFinder.setStrictFlag(!nostrict);
        OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType transition_group_map;
        OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);
        OpenSwath::SpectrumAccessPtr chromatogram_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
        std::vector< OpenSwath::SwathMap > swath_maps(1);
        swath_maps[0].sptr = swath_ptr;
        featureFinder.pickExperiment(chromatogram_ptr, featureFile,
                                     transition_exp_used, trafo, swath_maps, transition_group_map);

        // write all features and the protein identifications from tmp_featureFile into featureFile
#ifdef _OPENMP
#pragma omp critical (featureFinder)
#endif
        {
          for (FeatureMap::iterator feature_it = featureFile.begin();
               feature_it != featureFile.end(); ++feature_it)
          {
            out_featureFile.push_back(*feature_it);
          }
          for (std::vector<ProteinIdentification>::iterator protid_it =
                 featureFile.getProteinIdentifications().begin();
               protid_it != featureFile.getProteinIdentifications().end();
               ++protid_it)
          {
            out_featureFile.getProteinIdentifications().push_back(*protid_it);
          }

        }
      } // end of do_continue
    } // end of loop over all files / end of OpenMP

    addDataProcessing_(out_featureFile, getProcessingInfo_(DataProcessing::QUANTITATION));
    out_featureFile.ensureUniqueId();
    FeatureXMLFile().store(out, out_featureFile);

    return EXECUTION_OK;
  }

};

int main(int argc, const char **argv)
{
  TOPPOpenSwathAnalyzer tool;
  return tool.main(argc, argv);
}

/// @endcond
