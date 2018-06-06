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
// $Authors: Sandro Andreotti $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGeneratorTrainer.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>


using namespace std;
using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------
//\link CHEMISTRY_SvmTheoreticalSpectrumGeneratorTrainer.
/**
  @page UTILS_SvmTheoreticalSpectrumGeneratorTrainer SvmTheoreticalSpectrumGeneratorTrainer

  @brief Trainer for SVM model as input for SvmTheoreticalSpectrumGenerator.

  This application requires mzML file with ms2 spectra and annotations in an idXml file and trains a SVM model usable by
  SvmTheoreticalSpectrumGenerator. Please refer to the documentation of the corresponding class @ref OpenMS::SvmTheoreticalSpectrumGeneratorTrainer

  @note This tool is experimental!

  @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_SvmTheoreticalSpectrumGeneratorTrainer.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_SvmTheoreticalSpectrumGeneratorTrainer.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class SvmTheoreticalSpectrumGeneratorTrainerTOPP :
  public TOPPBase
{
  typedef SvmTheoreticalSpectrumGenerator::IonType IonType;

public:
  SvmTheoreticalSpectrumGeneratorTrainerTOPP() :
    TOPPBase("SvmTheoreticalSpectrumGeneratorTrainer", "Trainer for SVM models as input for SvmTheoreticalSpectrumGenerator", false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    // I/O settings
    registerInputFile_("in_spectra", "<file>", "", "Input Training Spectra in mzML", true);
    setValidFormats_("in_spectra",  ListUtils::create<String>("mzML"));
    registerInputFile_("in_identifications", "<file>", "", "Input file with corresponding sequences in idXML", true);
    setValidFormats_("in_identifications",  ListUtils::create<String>("idXML"));
    registerOutputFile_("model_output_file", "<file>", "",
                        "Name for output files. For each ion_type one file <filename>_residue_loss_charge.svm and one <filename>.info which has to be passed to the SvmTheoretical SpectrumGenerator", true);
    //TODO: check how to handle file prefix properly for TOPPAS/KNIME/CTD
    registerIntOption_("precursor_charge", "<Int>", 2, "Precursor charge state used for model training", false);
    setMinInt_("precursor_charge", 1);
    setMaxInt_("precursor_charge", 3);
    registerFlag_("write_training_files", "No models are trained but input training files for libSVM command line tools are produced");

    registerSubsection_("algorithm", "");
  }

  Param getSubsectionDefaults_(const String & /* section*/) const override
  {
    Param tmp = SvmTheoreticalSpectrumGeneratorTrainer().getDefaults();
    tmp.remove("write_training_files");
    return tmp;
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String in_spectra = getStringOption_("in_spectra");
    String in_identifications = getStringOption_("in_identifications");
    String outfile = getStringOption_("model_output_file");
    Int precursor_charge = getIntOption_("precursor_charge");

    //-------------------------------------------------------------
    // init SvmTheoreticalSpectrumGeneratorTrainer
    //-------------------------------------------------------------
    SvmTheoreticalSpectrumGeneratorTrainer trainer;

    Param param = getParam_().copy("algorithm:", true);
    String write_files = getFlag_("write_training_files") ? "true" : "false";
    param.setValue("write_training_files", write_files);
    trainer.setParameters(param);

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    PeakMap map;
    MzMLFile().load(in_spectra, map);

    std::vector<PeptideIdentification> pep_ids;
    std::vector<ProteinIdentification> prot_ids;
    String tmp_str;
    IdXMLFile().load(in_identifications, prot_ids, pep_ids, tmp_str);

    IDMapper idmapper;
    Param par;
    par.setValue("rt_tolerance", 0.001);
    par.setValue("mz_tolerance", 0.001);
    idmapper.setParameters(par);
    idmapper.annotate(map, pep_ids, prot_ids);

    //generate vector of annotations
    std::vector<AASequence> annotations;
    PeakMap::iterator it;
    for (it = map.begin(); it != map.end(); ++it)
    {
      annotations.push_back(it->getPeptideIdentifications()[0].getHits()[0].getSequence());
    }

    trainer.trainModel(map, annotations, outfile, precursor_charge);
    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  SvmTheoreticalSpectrumGeneratorTrainerTOPP tool;
  return tool.main(argc, argv);
}

/// @endcond
