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
// $Authors: George Rosenberger, Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_TargetedFileConverter TargetedFileConverter

  @brief Converts different spectral libraries / transition files for targeted proteomics and metabolomics analysis.
  
  Can convert multiple formats to and from TraML (standardized transition format). The following formats are supported:

        <ul>
          <li> TraML </li>
          <li> OpenSWATH TSV transition lists </li>
          <li> OpenSWATH PQP SQLite files </li>
          <li> SpectraST MRM transition lists </li>
          <li> Skyline transition lists </li>
          <li> Spectronaut transition lists </li>
        </ul>

  Transition lists can be either comma- or tab-separated. Modifications should be provided in UniMod format<sup>1</sup>, but can also be provided in TPP format. The following columns are required:

        <ul>
          <li> PrecursorMz* (float) </li>
          <li> ProductMz* (float; synonyms: FragmentMz) </li>
          <li> LibraryIntensity* (float; synonyms: RelativeFragmentIntensity) </li>
          <li> NormalizedRetentionTime* (float; synonyms: RetentionTime, Tr_recalibrated, iRT, RetentionTimeCalculatorScore) (normalized retention time) </li>
        </ul>

  For targeted proteomics files, the following additional columns should be provided:
        <ul>
          <li> ProteinId** (free text; synonyms: ProteinName) </li>
          <li> PeptideSequence** (free text, sequence only (no modifications); synonyms: Sequence, StrippedSequence) </li>
          <li> ModifiedPeptideSequence** (free text, should contain modifications<sup>1</sup>; synonyms: FullUniModPeptideName, FullPeptideName, ModifiedSequence)  </li>
          <li> PrecursorCharge** (integer, contains the charge of the precursorl synonyms: Charge) </li>
          <li> ProductCharge** (integer, contains the fragment charge; synonyms: FragmentCharge) </li>
          <li> FragmentType (free text, contains the type of the fragment, e.g. "b" or "y") </li>
          <li> FragmentSeriesNumber (integer, e.g. for y7 use "7" here; synonyms: FragmentNumber) </li>
        </ul>

  OpenSWATH uses grouped transitions to detect candidate analyte signals. These groups are by default generated based on the input, but can also be manually specified:

        <ul>
          <li> TransitionGroupId (free text, designates the transition group [e.g. peptide] to which this transition belongs; synomys: TransitionGroupName, transition_group_id) </li>
          <li> TransitionId (free text, needs to be unique for each transition [in this file]; synonyms: TransitionName, transition_name) </li>
          <li> Decoy (1: decoy, 0: target, i.e. no decoy; determines whether the transition is a decoy transition or not, synomys: decoy, isDecoy) </li>
          <li> PeptideGroupLabel (free text, designates to which peptide label group (as defined in MS:1000893) the peptide belongs to<sup>2</sup>) </li>
          <li> DetectingTransition (1: use transition to detect peak group, 0: don't use transition for detection; synonyms: detecting_transition) </li>
          <li> IdentifyingTransition (1: use transition for peptidoform inference using IPF, 0: don't use transition for identification; synonyms: identifying_transition) </li>
          <li> QuantifyingTransition (1: use transition to quantify peak group, 0: don't use transition for quantification; synonyms: quantifying_transition) </li>
        </ul>

  Optionally, the following columns can be specified but they are not actively used by OpenSWATH:
        <ul>
          <li> CollisionEnergy (float; synonyms: CE) </li>
          <li> Annotation (free text, e.g. y7) </li>
          <li> UniprotId (free text; synonyms: UniprotID) </li>
          <li> LabelType (free text, optional description of which label was used, e.g. heavy or light) </li>
        </ul>

  For targeted metabolomics files, the following fields are also supported:

        <ul>
          <li> CompoundName** (synonyms: CompoundId) </li>
          <li> SMILES </li>
          <li> SumFormula </li>
        </ul>

  Fields indicated with * are strictly required while fields indicated with **
  are only required in the specific context (proteomics or metabolomics).

<p>
Remarks:
</p>
<ul>
  <li>
    1. modifications should be supplied inside the sequence using UniMod
      identifiers or freetext identifiers that are understood by %OpenMS. <br/>
      example: PEPT(Phosphorylation)IDE(UniMod:27)A )
  </li>
  <li>
    2. peptide label groups designate groups of peptides that are isotopically
    modified forms of the same peptide species. For example, the heavy and
    light forms of the same peptide will both be assigned the same peptide
    group label.  <br/>
      example: <br/>
      PEPTIDEAK -> gets label "PEPTIDEAK_gr1" <br/>
      PEPTIDEAK[+8] -> gets label "PEPTIDEAK_gr1" <br/>
      PEPT(Phosphorylation)IDEAK -> gets label "PEPTIDEAK_gr2" <br/>
      PEPT(Phosphorylation)IDEAK[+8] -> gets label "PEPTIDEAK_gr2" <br/>
  </li>
</ul>
</p>

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_TargetedFileConverter.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_TargetedFileConverter.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPTargetedFileConverter :
  public TOPPBase
{
public:

  TOPPTargetedFileConverter() :
    TOPPBase("TargetedFileConverter", "Converts different transition files for targeted proteomics / metabolomics analysis.", false)
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file to convert.\n "
                                           "See http://www.openms.de/current_doxygen/html/UTILS_TargetedFileConverter.html for format of OpenSWATH transition TSV file or SpectraST MRM file.");
    registerStringOption_("in_type", "<type>", "", "input file type -- default: determined from file extension or content\n", false);
    String formats("tsv,mrm,pqp,TraML");
    setValidFormats_("in", ListUtils::create<String>(formats));
    setValidStrings_("in_type", ListUtils::create<String>(formats));

    formats = "tsv,pqp,TraML";
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>(formats));
    registerStringOption_("out_type", "<type>", "", "Output file type -- default: determined from file extension or content\nNote: that not all conversion paths work or make sense.", false);
    setValidStrings_("out_type", ListUtils::create<String>(formats));

    registerSubsection_("algorithm", "Algorithm parameters section");
    registerFlag_("legacy_traml_id", "PQP to TraML: Should legacy TraML IDs be used?", true);

  }

  Param getSubsectionDefaults_(const String&) const override
  {
    return TransitionTSVFile().getDefaults();
  }

  ExitCodes main_(int, const char**) override
  {
    FileHandler fh;

    //input file type
    String in = getStringOption_("in");
    FileTypes::Type in_type = FileTypes::nameToType(getStringOption_("in_type"));

    if (in_type == FileTypes::UNKNOWN)
    {
      in_type = fh.getType(in);
      writeDebug_(String("Input file type: ") + FileTypes::typeToName(in_type), 2);
    }

    if (in_type == FileTypes::UNKNOWN)
    {
      writeLog_("Error: Could not determine input file type!");
      return PARSE_ERROR;
    }

    //output file names and types
    String out = getStringOption_("out");
    FileTypes::Type out_type = FileTypes::nameToType(getStringOption_("out_type"));

    if (out_type == FileTypes::UNKNOWN)
    {
      out_type = fh.getTypeByFileName(out);
    }

    if (out_type == FileTypes::UNKNOWN)
    {
      writeLog_("Error: Could not determine output file type!");
      return PARSE_ERROR;
    }

    bool legacy_traml_id = getFlag_("legacy_traml_id");

    //--------------------------------------------------------------------------- 
    // Start Conversion
    //--------------------------------------------------------------------------- 
    TargetedExperiment targeted_exp;
    if (in_type == FileTypes::TSV || in_type == FileTypes::MRM)
    {
      const char* tr_file = in.c_str();
      Param reader_parameters = getParam_().copy("algorithm:", true);
      TransitionTSVFile tsv_reader = TransitionTSVFile();
      tsv_reader.setLogType(log_type_);
      tsv_reader.setParameters(reader_parameters);
      tsv_reader.convertTSVToTargetedExperiment(tr_file, in_type, targeted_exp);
      tsv_reader.validateTargetedExperiment(targeted_exp);
    }
    else if (in_type == FileTypes::PQP)
    {
      const char* tr_file = in.c_str();
      TransitionPQPFile pqp_reader = TransitionPQPFile();
      Param reader_parameters = getParam_().copy("algorithm:", true);
      pqp_reader.setLogType(log_type_);
      pqp_reader.setParameters(reader_parameters);
      pqp_reader.convertPQPToTargetedExperiment(tr_file, targeted_exp, legacy_traml_id);
      pqp_reader.validateTargetedExperiment(targeted_exp);
    }
    else if (in_type == FileTypes::TRAML)
    {
      TraMLFile traml;
      traml.load(in, targeted_exp);
    }

    if (out_type == FileTypes::TSV)
    {
      const char* tr_file = out.c_str();
      TransitionTSVFile tsv_reader = TransitionTSVFile();
      tsv_reader.setLogType(log_type_);
      tsv_reader.convertTargetedExperimentToTSV(tr_file, targeted_exp);
    }
    if (out_type == FileTypes::PQP)
    {
      const char * tr_file = out.c_str();
      TransitionPQPFile pqp_reader = TransitionPQPFile();
      pqp_reader.setLogType(log_type_);
      pqp_reader.convertTargetedExperimentToPQP(tr_file, targeted_exp);
    }
    else if (out_type == FileTypes::TRAML)
    {
      TraMLFile traml;
      traml.store(out, targeted_exp);
    }

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{

  TOPPTargetedFileConverter tool;
  return tool.main(argc, argv);
}

/// @endcond
