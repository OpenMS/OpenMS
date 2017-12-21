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

#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>


using namespace std;

#ifdef _OPENMP
#include <omp.h>
#endif

// #ifdef _OPENMP
//   #define IF_MASTERTHREAD if (omp_get_thread_num() ==0)  
// #else
//   #define IF_MASTERTHREAD 
// #endif    

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_OpenSwathChromatogramExtractor OpenSwathChromatogramExtractor

  @brief Extracts chromatograms (XICs) from a file containing spectra.

  <CENTER>
      <table>
          <tr>
              <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
              <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ OpenSwathChromatogramExtractor \f$ \longrightarrow \f$</td>
              <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
          </tr>
          <tr>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileFilter </td>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_OpenSwathAnalyzer </td>
          </tr>
      </table>
  </CENTER>

  This module extracts ion traces (extracted ion chromatograms or XICs) from a
  file containing spectra.  The masses at which the chromatograms should be
  extracted are stored in a TraML file and the result is stored in a mzML file
  holding chromatograms. This tool is designed to extract chromatograms from
  SWATH (data independent acquisition) data (see ref[1]), thus it will extract
  the masses found in the product ion section of the TraML transitions,
  returning as many chromatograms as input transitions were provided.

  For SWATH data, the "is_swath" flag which will check the precursor
  isolation window of the first scan and assume all scans in that file were
  recorded with this precursor window (thus making it necessary to provide one
  input file per SWATH window). The module will then only extract transitions
  whose precursors fall into the corresponding isolation window.

  For the extraction method, two convolution functions are available: top-hat
  and bartlett. While top-hat will just sum up the signal within a quadratic
  window, bartlett will weigh the signal in the center of the window more than
  the signal on the edge.

  [1] Gillet LC, Navarro P, Tate S, Röst H, Selevsek N, Reiter L, Bonner R, Aebersold R. \n
  <a href="http://dx.doi.org/10.1074/mcp.O111.016717"> Targeted data extraction of the MS/MS spectra generated by data-independent
  acquisition: a new concept for consistent and accurate proteome analysis. </a> \n
  Mol Cell Proteomics. 2012 Jun;11(6):O111.016717. 

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_OpenSwathChromatogramExtractor.cli
 
  <B>The algorithm parameters for the Analyzer filter are:</B>
  @htmlinclude TOPP_OpenSwathChromatogramExtractor.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPOpenSwathChromatogramExtractor 
  : public TOPPBase
{
public:

  TOPPOpenSwathChromatogramExtractor() 
    : TOPPBase("OpenSwathChromatogramExtractor", "Extract chromatograms (XIC) from a MS2 map file.", true)
  {
  }

protected:

  typedef PeakMap MapType;

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", StringList(), "Input files separated by blank");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFile_("tr", "<file>", "", "transition file ('TraML' or 'csv')");
    setValidFormats_("tr", ListUtils::create<String>("csv,traML"));
    
    registerInputFile_("rt_norm", "<file>", "", "RT normalization file (how to map the RTs of this run to the ones stored in the library)", false);
    setValidFormats_("rt_norm", ListUtils::create<String>("trafoXML"));

    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", ListUtils::create<String>("mzML"));

    registerDoubleOption_("min_upper_edge_dist", "<double>", 0.0, "Minimal distance to the edge to still consider a precursor, in Thomson", false);
    registerDoubleOption_("mz_window", "<double>", 0.05, "Extraction window in m/z dimension (in Thomson, to use ppm see -ppm flag). This is the full window size, e.g. 100 ppm would extract 50 ppm on either side.", false);
    registerDoubleOption_("rt_window", "<double>", -1, "Extraction window in RT dimension (-1 means extract over the whole range). This is the full window size, e.g. a value of 1000 seconds would extract 500 seconds on either side.", false);
    setMinFloat_("mz_window", 0.0);

    registerFlag_("is_swath", "Set this flag if the data is SWATH data");
    registerFlag_("ppm", "m/z extraction_window is in ppm");

    registerFlag_("extract_MS1", "Extract the MS1 transitions based on the precursor values in the TraML file");

    registerStringOption_("extraction_function", "<name>", "tophat", "Function used to extract the signal", false);
    StringList model_types;
    model_types.push_back("tophat");
    model_types.push_back("bartlett"); // bartlett if we use zeros at the end
    setValidStrings_("extraction_function", model_types);

    registerModelOptions_("linear");
  }

  void registerModelOptions_(const String & default_model)
  {
    registerTOPPSubsection_("model", "Options to control the modeling of retention time transformations from data");
    registerStringOption_("model:type", "<name>", default_model, "Type of model", false);
    StringList model_types;
    TransformationDescription::getModelTypes(model_types);
    if (!ListUtils::contains(model_types, default_model)) {
      model_types.insert(model_types.begin(), default_model);
    }
    setValidStrings_("model:type", model_types);
    registerFlag_("model:symmetric_regression", "Only for 'linear' model: Perform linear regression on 'y - x' vs. 'y + x', instead of on 'y' vs. 'x'.");
  }

  ExitCodes main_(int, const char **) override
  {
    StringList file_list = getStringList_("in");
    String tr_file_str = getStringOption_("tr");
    String out = getStringOption_("out");
    bool is_swath = getFlag_("is_swath");
    bool ppm = getFlag_("ppm");
    bool extract_MS1 = getFlag_("extract_MS1");
    double min_upper_edge_dist = getDoubleOption_("min_upper_edge_dist");
    double mz_extraction_window = getDoubleOption_("mz_window");
    double rt_extraction_window = getDoubleOption_("rt_window");

    String extraction_function = getStringOption_("extraction_function");

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
    TransformationDescription trafo_inverse = trafo;
    trafo_inverse.invert();

    const char * tr_file = tr_file_str.c_str();

    MapType out_exp;
    std::vector< OpenMS::MSChromatogram > chromatograms;
    TraMLFile traml;
    OpenMS::TargetedExperiment targeted_exp;

    std::cout << "Loading TraML file" << std::endl;
    traml.load(tr_file, targeted_exp);
    std::cout << "Loaded TraML file" << std::endl;

    // Do parallelization over the different input files
    // Only in OpenMP 3.0 are unsigned loop variables allowed
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize i = 0; i < boost::numeric_cast<SignedSize>(file_list.size()); ++i)
    {
      boost::shared_ptr<PeakMap > exp(new PeakMap);
      MzMLFile f;
      // Logging and output to the console
      // IF_MASTERTHREAD f.setLogType(log_type_); 

      // Find the transitions to extract and extract them
      MapType tmp_out;
      OpenMS::TargetedExperiment transition_exp_used;
      f.load(file_list[i], *exp);
      if (exp->empty() ) { continue; } // if empty, go on
      OpenSwath::SpectrumAccessPtr expptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
      bool do_continue = true;
      if (is_swath)
      {
        do_continue = OpenSwathHelper::checkSwathMapAndSelectTransitions(*exp, targeted_exp, transition_exp_used, min_upper_edge_dist);  
      }
      else
      {
        transition_exp_used = targeted_exp;
      }

#ifdef _OPENMP
#pragma omp critical (OpenSwathChromatogramExtractor_metadata)
#endif
      // after loading the first file, copy the meta data from that experiment
      // this may happen *after* chromatograms were already added to the
      // output, thus we do NOT fill the experiment here but rather store all
      // the chromatograms in the "chromatograms" array and store them in
      // out_exp afterwards.
      if (i == 0) 
      {
        out_exp = *exp;
        out_exp.clear(false);
      }

      std::cout << "Extracting " << transition_exp_used.getTransitions().size() << " transitions" << std::endl;
      std::vector< OpenSwath::ChromatogramPtr > chromatogram_ptrs;
      std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinates;

      // continue if the map is not empty
      if (do_continue)
      {

        // Prepare the coordinates (with or without rt extraction) and then extract the chromatograms
        ChromatogramExtractor extractor;
        if (rt_extraction_window < 0)
        {
          extractor.prepare_coordinates(chromatogram_ptrs, coordinates, transition_exp_used, rt_extraction_window, extract_MS1);
        }
        else
        {
          // Use an rt extraction window of 0.0 which will just write the retention time in start / end positions
          extractor.prepare_coordinates(chromatogram_ptrs, coordinates, transition_exp_used, 0.0, extract_MS1);
          for (std::vector< ChromatogramExtractor::ExtractionCoordinates >::iterator it = coordinates.begin(); it != coordinates.end(); ++it)
          {
            it->rt_start = trafo_inverse.apply(it->rt_start) - rt_extraction_window / 2.0;
            it->rt_end = trafo_inverse.apply(it->rt_end) + rt_extraction_window / 2.0;
          }
        }
        extractor.extractChromatograms(expptr, chromatogram_ptrs, coordinates, 
            mz_extraction_window, ppm, extraction_function);

#ifdef _OPENMP
#pragma omp critical (OpenSwathChromatogramExtractor_insertMS1)
#endif
        {
          // Remove potential meta value indicating cached data
          SpectrumSettings exp_settings = (*exp)[0];
          for (Size j = 0; j < exp_settings.getDataProcessing().size(); j++)
          {
            if (exp_settings.getDataProcessing()[j]->metaValueExists("cached_data"))
            { exp_settings.getDataProcessing()[j]->removeMetaValue("cached_data"); }
          }
          extractor.return_chromatogram(chromatogram_ptrs, coordinates, transition_exp_used, exp_settings, chromatograms, extract_MS1);
        }

      } // end of do_continue
    } // end of loop over all files / end of OpenMP

    // TODO check that no chromatogram IDs occur multiple times !
    
    // store the output
    out_exp.setChromatograms(chromatograms);
    MzMLFile mzf;
    mzf.setLogType(log_type_); 
    addDataProcessing_(out_exp, getProcessingInfo_(DataProcessing::SMOOTHING));
    mzf.store(out, out_exp);

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPOpenSwathChromatogramExtractor tool;
  return tool.main(argc, argv);
}

/// @endcond
