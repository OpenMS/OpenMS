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
// $Maintainer: Chris Bielow $
// $Authors: Stephan Aiche, Chris Bielow, Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <iostream>
#include <string>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/SYSTEM/StopWatch.h>

#include <OpenMS/SIMULATION/MSSim.h>
#include <OpenMS/SIMULATION/SimTypes.h>

// file types
#include <OpenMS/FORMAT/DTA2DFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
using namespace OpenMS;
using namespace std;


/**
  @page UTILS_MSSimulator MSSimulator

  @brief A highly configurable simulator for mass spectrometry experiments.

  This implementation is described in
  <p>
  Bielow C, Aiche S, Andreotti S, Reinert K<br>
  MSSimulator: Simulation of Mass Spectrometry Data<br>
  Journal of Proteome Research (2011), DOI: 10.1021/pr200155f<br>
  </p>

  The most important features are:
  <ul>
    <li>Simulation of Capillary electrophoresis and HPLC as separation step</li>
    <li>Simulation of MS spectra</li>
    <li>Simulation of MS/MS spectra with configurable precursor-selection strategy</li>
    <li>Simulation of iTRAQ labels</li>
    <li>Simulation of different noise models and instrument types (resolution, peak shape)</li>
  </ul>

  Look at the INI file (via "MSSimulator -write_ini myini.ini") to see the available parameters and more functionality.

  <h3>Input: FASTA files</h3>
  Protein sequences (including amino acid modifications) can be provided as FASTA file.
  We allow a special tag in the description of each entry to specify protein abundance.
  If you want to create a complex FASTA file with a Gaussian protein abundance model in log space,
  see our Python script shipping with your OpenMS installation (e.g., <OpenMS-dir>/share/OpenMS/examples/simulation/FASTAProteinAbundanceSampling.py).
  It supports (random) sampling from a large FASTA file, protein weight filtering and adds an
  intensity tag to each entry.

  If multiplexed data is simulated (like SILAC or iTRAQ) you need to supply multiple FASTA input files.
  For the label-free setting, all FASTA input files will be merged into one, before simulation.

  <p>
   For MS/MS simulation only a test model is shipped with %OpenMS.<br>
   Please find trained models at: http://sourceforge.net/projects/open-ms/files/Supplementary/Simulation/.
  </p>

  To specify intensity values for certain proteins, add an abundance tag for the corresponding protein in the FASTA input file:<br>
    - add '[# &lt;key&gt;=&lt;value&gt; #]' at the end of the &gt; line to specify intensity
  For RT control (disable digestion, to make this work!)
  <ul>
    <li> rt (subjected to small local error by randomization)
    <li> RT (used as is without local error)
  </ul>

  For amino acid modifications, insert their name at the respective amino acid residues. The modifications are fixed. If you need variable modifications,
  you have to add the desired combinatorial variants (presence/absence of one or all modifications) to the FASTA file.
  Valid modification names are listed in many TOPP/UTILS, e.g @ref TOPP_MSGFPlusAdapter 's @em -fixed_modifications parameter.

e.g.
@code
>seq1 optional comment [# intensity=567.4 #]
(Acetyl).M(Oxidation)ASQYLATARHGC(Carbamidomethyl)FLPRHRDTGILP
>seq2 optional comment [# intensity=117.4, RT=405.3 #]
QKRPSQRHGLATAC(Carbamidomethyl)RHGTGGGDRAT.(Dehydrated)
@endcode

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_MSSimulator.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_MSSimulator.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMSSimulator :
  public TOPPBase
{
public:
  TOPPMSSimulator() :
    TOPPBase("MSSimulator", "A highly configurable simulator for mass spectrometry experiments.", false)
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    // I/O settings
    registerInputFileList_("in", "<files>", ListUtils::create<String>(""), "Input protein sequences", true, false);
    setValidFormats_("in", ListUtils::create<String>("FASTA"));
    registerOutputFile_("out", "<file>", "", "output: simulated MS raw (profile) data", false);
    setValidFormats_("out", ListUtils::create<String>("mzML"));
    registerOutputFile_("out_pm", "<file>", "", "output: ground-truth picked (centroided) MS data", false);
    setValidFormats_("out_pm", ListUtils::create<String>("mzML"));
    registerOutputFile_("out_fm", "<file>", "", "output: ground-truth features", false);
    setValidFormats_("out_fm", ListUtils::create<String>("featureXML"));
    registerOutputFile_("out_cm", "<file>", "", "output: ground-truth features, grouping ESI charge variants of each parent peptide", false);
    setValidFormats_("out_cm", ListUtils::create<String>("consensusXML"));
    registerOutputFile_("out_lcm", "<file>", "", "output: ground-truth features, grouping labeled variants", false);
    setValidFormats_("out_lcm", ListUtils::create<String>("consensusXML"));
    registerOutputFile_("out_cntm", "<file>", "", "output: ground-truth features caused by contaminants", false);
    setValidFormats_("out_cntm", ListUtils::create<String>("featureXML"));
    registerOutputFile_("out_id", "<file>", "", "output: ground-truth MS2 peptide identifications", false);
    setValidFormats_("out_id", ListUtils::create<String>("idXML"));

    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String& /*section*/) const
  {
    Param tmp;
    tmp.insert("MSSim:", MSSim().getParameters());

    // set parameters for the different types of random number generators
    // we support one for the technical and one for the biological variability
    tmp.setValue("RandomNumberGenerators:biological", "random", "Controls the 'biological' randomness of the generated data (e.g. systematic effects like deviations in RT). If set to 'random' each experiment will look different. If set to 'reproducible' each experiment will have the same outcome (given that the input data is the same).");
    tmp.setValidStrings("RandomNumberGenerators:biological", ListUtils::create<String>("reproducible,random"));
    tmp.setValue("RandomNumberGenerators:technical", "random", "Controls the 'technical' randomness of the generated data (e.g. noise in the raw signal). If set to 'random' each experiment will look different. If set to 'reproducible' each experiment will have the same outcome (given that the input data is the same).");
    tmp.setValidStrings("RandomNumberGenerators:technical", ListUtils::create<String>("reproducible,random"));
    tmp.setSectionDescription("RandomNumberGenerators", "Parameters for generating the random aspects (e.g. noise) in the simulated data. The generation is separated into two parts, the technical part, like noise in the raw signal, and the biological part, like systematic deviations in the predicted retention times.");
    return tmp;
  }

  // Load proteins from FASTA file
  void loadFASTA_(const String& filename, SimTypes::SampleProteins& proteins)
  {
    writeLog_(String("Loading sequence data from ") + filename +  String(" ..."));

    FASTAFile fastafile;
    typedef std::vector<FASTAFile::FASTAEntry> FASTAdata;
    FASTAdata fastadata;

    // load FASTA file contents
    fastafile.load(filename, fastadata);

    // add data from file to protein storage
    String::size_type index;

    StringList valid_meta_values = ListUtils::create<String>("intensity,RT,rt");
    // re-parse FASTA description to obtain quantitation info
    for (FASTAdata::iterator it = fastadata.begin(); it != fastadata.end(); ++it)
    {
      // remove all ambiguous characters from FASTA entry
      // TODO: this is somehow problematic since we modify user input
      it->sequence.remove('X');
      it->sequence.remove('B');
      it->sequence.remove('Z');

      // parsed abundance
      MetaInfoInterface data;
      data.setMetaValue("intensity", 10000.0);

      // Look for a relative quantity given in the comment line of a FASTA entry
      // e.g. >BSA [#120]
      index = (it->description).find("[#");
      // if found, extract and set relative quantity accordingly
      if (index != String::npos)
      {
        String::size_type index_end = (it->description).find(']', index);
        if (index_end == String::npos) throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MSSimulator: Invalid entry (" + it->identifier + ") in FASTA file; abundance section has open tag '[#' but missing close tag ']'.");

        //std::cout << (it->description).substr(index+2,index_end-index-2) << std::endl;
        StringList meta_values = ListUtils::create<String>((it->description).substr(index + 2, index_end - index - 3).removeWhitespaces(), ',');
        for (Size i = 0; i < meta_values.size(); ++i)
        {
          StringList components;
          meta_values[i].split('=', components);
          if (components.size() != 2) throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MSSimulator: Invalid entry (" + it->identifier + ") in FASTA file; the component '" + meta_values[i] + "' is missing an assignment ('=').");
          // check if component is known
          if (!ListUtils::contains(valid_meta_values, components[0])) throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MSSimulator: Invalid entry (" + it->identifier + ") in FASTA file; the component '" + meta_values[i] + "' has an unsupported meta value.");

          if (components[0] == "intensity" || String(components[0]).toUpper() == "RT")
          {
            data.setMetaValue(components[0], components[1].toDouble());
          }
          else
          {
            data.setMetaValue(components[0], components[1]);
          }
        }
      }

      proteins.push_back(SimTypes::SimProtein(*it, data));
    }

    writeLog_(String("done (") + fastadata.size() + String(" protein(s) loaded)"));
  }

  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    // check if at least one output file is
    if (getStringOption_("out") == "" &&
        getStringOption_("out_pm") == "" &&
        getStringOption_("out_fm") == "" &&
        getStringOption_("out_cm") == "" &&
        getStringOption_("out_lcm") == "" &&
        getStringOption_("out_cntm") == "" &&
        getStringOption_("out_id") == "")
    {
      LOG_ERROR << "Error: At least one output file needs to specified!" << std::endl;
      return MISSING_PARAMETERS;
    }

    MSSim ms_simulation;
    ms_simulation.setParameters(getParam_().copy("algorithm:MSSim:", true));

    // read proteins
    SimTypes::SampleChannels channels;
    StringList input_files = getStringList_("in");
    for (Size i = 0; i < input_files.size(); ++i)
    {
      SimTypes::SampleProteins proteins;
      loadFASTA_(input_files[i], proteins);
      channels.push_back(proteins);
    }

    // initialize the random number generators
    bool biological_random = getParam_().getValue("algorithm:RandomNumberGenerators:biological") == "random";
    bool technical_random = getParam_().getValue("algorithm:RandomNumberGenerators:technical") == "random";
    SimTypes::MutableSimRandomNumberGeneratorPtr rnd_gen(new SimTypes::SimRandomNumberGenerator);
    rnd_gen->initialize(biological_random, technical_random);

    ms_simulation.setLogType(this->log_type_);

    // start simulation
    writeLog_("Starting simulation");
    StopWatch w;

    w.start();
    ms_simulation.simulate(rnd_gen, channels);
    w.stop();
    writeLog_(String("Simulation took ") + String(w.getClockTime()) + String(" seconds"));

    String outputfile_name = getStringOption_("out");
    if (outputfile_name != "")
    {
      writeLog_(String("Storing simulated raw data in: ") + outputfile_name);
      MzMLFile().store(outputfile_name, ms_simulation.getExperiment());
    }

    String pxml_out = getStringOption_("out_pm");
    if (pxml_out != "")
    {
      writeLog_(String("Storing simulated peak/centroided data in: ") + pxml_out);
      MzMLFile().store(pxml_out, ms_simulation.getPeakMap());
    }

    String fxml_out = getStringOption_("out_fm");
    if (fxml_out != "")
    {
      writeLog_(String("Storing simulated features in: ") + fxml_out);
      FeatureXMLFile().store(fxml_out, ms_simulation.getSimulatedFeatures());
    }

    String cxml_out = getStringOption_("out_cm");
    if (cxml_out != "")
    {
      writeLog_(String("Storing charged consensus features in: ") + cxml_out);

      ConsensusMap& charge_consensus = ms_simulation.getChargeConsensus();
      charge_consensus.getFileDescriptions()[0].filename = fxml_out;
      charge_consensus.getFileDescriptions()[0].size = ms_simulation.getSimulatedFeatures().size();
      charge_consensus.getFileDescriptions()[0].unique_id = ms_simulation.getSimulatedFeatures().getUniqueId();

      ConsensusXMLFile().store(cxml_out, charge_consensus);
    }

    String lcxml_out = getStringOption_("out_lcm");
    if (lcxml_out != "")
    {
      writeLog_(String("Storing labeling consensus features in: ") + lcxml_out);

      // set file name for all (sub)feature maps
      ConsensusMap& labeling_consensus = ms_simulation.getLabelingConsensus();
      for (ConsensusMap::FileDescriptions::iterator fdI = labeling_consensus.getFileDescriptions().begin();
           fdI != labeling_consensus.getFileDescriptions().end();
           ++fdI)
      {
        fdI->second.filename = fxml_out;
      }

      ConsensusXMLFile().store(lcxml_out, labeling_consensus);
    }

    String cntxml_out = getStringOption_("out_cntm");
    if (cntxml_out != "")
    {
      writeLog_(String("Storing simulated contaminant features in: ") + cntxml_out);
      FeatureXMLFile().store(cntxml_out, ms_simulation.getContaminants());
    }

    String id_out = getStringOption_("out_id");
    if (id_out != "")
    {
      writeLog_(String("Storing ground-truth peptide IDs in: ") + id_out);
      vector<ProteinIdentification> proteins;
      vector<PeptideIdentification> peptides;
      ms_simulation.getIdentifications(proteins, peptides);
      IdXMLFile().store(id_out, proteins, peptides);
    }

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPMSSimulator tool;
  return tool.main(argc, argv);
}

/// @endcond
