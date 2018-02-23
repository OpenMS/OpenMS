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
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderIdentificationAlgorithm.h>

#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_ProteomicLFQ
 **/

class UTILProteomicLFQ :
  public TOPPBase
{
public:
  UTILProteomicLFQ() :
    TOPPBase("ProteomicLFQ", "A standard proteomics LFQ pipeline.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<file list>", StringList(), "input files");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerInputFileList_("in_ids", "<file list>", StringList(), "unfiltered identifications");
    setValidFormats_("in_ids", ListUtils::create<String>("idXML,mzId"));
    registerInputFile_("design", "<file>", "", "design file");
    setValidFormats_("design", ListUtils::create<String>("tsv"));

    registerOutputFile_("out", "<file>", "", "output mzTab file");
    setValidFormats_("out", ListUtils::create<String>("tsv")); // TODO: add file extension for mzTab

    Param pp_defaults = PeakPickerHiRes().getDefaults();
    Param ff_defaults = FeatureFinderIdentificationAlgorithm().getDefaults();

    Param combined;
    combined.insert("Centroiding:", pp_defaults);
    combined.insert("Quantitation:", ff_defaults);
    // TODO: add params of other steps as well

    registerFullParam_(combined);
  }

/*
  ExperimentalDesign getExperimentalDesignIds_(
    const String & design_file, 
    const vector<ProteinIdentification> & proteins)
  {
    if (!design_file.empty()) // load experimental design file
    {
      return ExperimentalDesign::load(design_file);
      // TODO FRACTIONS: check if ed sane
    }
    else  // no design file provided
    {
      return ExperimentalDesign::fromIdentifications(proteins);
    }
  }
*/
  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    // read tool parameters
    StringList in = getStringList_("in");
    String out = getStringOption_("out");
    StringList in_ids = getStringList_("in_ids");
    String design_file = getStringOption_("design");
    
    //TODO: if no experimental design file is given, create a trivial von from protein ids (no fractions)
    //ExperimentalDesign design = getExperimentalDesignIds_(design_file, proteins);
    //
    ExperimentalDesign design = ExperimentalDesign::load(design_file);
    std::map<unsigned int, std::vector<String> > frac2ms = design.getFractionToMSFilesMapping();

    Param pp_param = getParam_().copy("Centroiding:", true);
    writeDebug_("Parameters passed to PeakPickerHiRes algorithm", pp_param, 3);
    PeakPickerHiRes pp;
    pp.setLogType(log_type_);
    pp.setParameters(pp_param);

    Param ff_param = getParam_().copy("Quantitation:", true);
    writeDebug_("Parameters passed to FeatureFinderIdentification algorithm", ff_param, 3);
    FeatureFinderIdentificationAlgorithm ff;
    ff.getProgressLogger().setLogType(log_type_);
    ff.setParameters(ff_param);

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    for (auto const ms_files : frac2ms) // for each fraction->ms file(s)
    {
      for (String const & mz_file : ms_files.second) // for each MS file
      {
        // TODO: check if s is part of in 
      
        // load raw file
        MzMLFile mzML_file;
        mzML_file.setLogType(log_type_);

        PeakMap ms_raw;
        mzML_file.load(mz_file, ms_raw);

        if (ms_raw.empty())
        {
          LOG_WARN << "The given file does not contain any spectra.";
          return INCOMPATIBLE_INPUT_DATA;
        }

        // check if spectra are sorted
        for (Size i = 0; i < ms_raw.size(); ++i)
        {
          if (!ms_raw[i].isSorted())
          {
            ms_raw[i].sortByPosition();
            writeLog_("Info: Sorte peaks by m/z.");
          }
        }

        //-------------------------------------------------------------
        // pick
        //-------------------------------------------------------------
        // TODO: only peak if not already picked (add auto mode that skips already picked ones)
        PeakMap ms_centroided;
        pp.pickExperiment(ms_raw, ms_centroided, true);
        // TODO: free memory of profile PeakMaps (if we needed to pick sth.), otherwise pass through

        // writing picked mzML files for data submission
        //annotate output with data processing info
        // TODO: how to store picked files? by specifying a folder? or by output files that match in number to input files
        // TODO: overwrite primaryMSRun with picked mzML name (for submission)
        // mzML_file.store(OUTPUTFILENAME, ms_centroided);
        // TODO: free all MS2 spectra (release memory!)
      }

      //-------------------------------------------------------------
      // feature detection
      //-------------------------------------------------------------
      // TODO: call FeatureFinderIdentification on all maps of this fraction 
      // TODO: free memory of centroided PeakMaps of this fraction

      //-------------------------------------------------------------
      // align
      //-------------------------------------------------------------
      // TODO: MapAlignerIdentification on all FeatureXMLs of this fraction
      
      //-------------------------------------------------------------
      // link
      //-------------------------------------------------------------
      // TODO: FeatureLinkerUnlabeledKD on all FeatureXMLs of this fraction
    }
    // TODO: FileMerger merge ids (here? or already earlier? filtered?)

    //-------------------------------------------------------------
    // Protein inference
    //-------------------------------------------------------------
    // TODO: ProteinInference on merged ids (how merged?)

    //-------------------------------------------------------------
    // Protein quantification and export to mzTab
    //-------------------------------------------------------------
    // TODO: ProteinQuantifier on (merged?) consensusXML (with 1% FDR?) + inference ids (unfiltered?)? 
    // export of MzTab file as final output

    return EXECUTION_OK;
  }
};


int main(int argc, const char ** argv)
{
  UTILProteomicLFQ tool;
  return tool.main(argc, argv);
}

/// @endcond
