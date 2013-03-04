// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVReader.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_ConvertTSVToTraML ConvertTSVToTraML

  @brief Converts OpenSWATH transition TSV files to TraML files

  The OpenSWATH transition TSV files need to have the following headers, all fields need to be separated by tabs:

    <CENTER>
        <table>

          <tr> PrecursorMz (float) </tr>
          <tr> ProductMz (float) </tr>
          <tr> Tr_calibrated (float) </tr>
          <tr> transition_name (free text, needs to be unique for each transition [in this file]) </tr>
          <tr> CE (float) </tr>
          <tr> LibraryIntensity (float) </tr>
          <tr> transition_group_id (free text, designates the transition group [e.g. peptide] to which this transition belongs) </tr>
          <tr> decoy (1==decoy, 0== no decoy; determines whether the transition is a decoy transition or not) </tr>
          <tr> PeptideSequence  (free text, sequence only (no modifications) ) </tr>
          <tr> ProteinName  (free text) </tr>
          <tr> Annotation  (free text, e.g. y7) </tr>
          <tr> FullUniModPeptideName  (free text, should contain modifications*)  </tr>
          <tr> PrecursorCharge (integer, contains the charge of the precursor) </tr>
          <tr> Labelgroup (free text, e.g. heavy or light) </tr>
          <tr> UniprotID (free text) </tr>
          <tr> FragmentType (free text, contains the type of the fragment, e.g. "b" or "y") </tr>
          <tr> FragmentCharge (integer, contains the fragment charge) </tr>
          <tr> FragmentSeriesNumber (integer, e.g. for y7 use "7" here) </tr>

        </table>
    </CENTER>


* modifications should be supplied inside the sequence using UniMod
  identifiers or freetext identifiers that are understood by OpenMS. Please do
  not use the ambiguous bracket notation (e.g. PEPT[+80]IDE or PEPT[181]IDE)
  since this is ambiguous and will NOT be interpreted correctly!
  example: PEPT(Phosphorylation)IDE(UniMod:27)A )

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPConvertTSVToTraML : public TOPPBase
{
public:

  TOPPConvertTSVToTraML() :
  TOPPBase("ConvertTSVToTraML", "Converts an OpenSWATH transition TSV file to a TraML file")
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input OpenSWATH transition TSV file");
    setValidFormats_("in", StringList::create("csv"));

    registerOutputFile_("out", "<file>", "", "Output TraML file");
    setValidFormats_("out", StringList::create("TraML"));

  }

  ExitCodes main_(int, const char **)
  {
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    const char * tr_file = in.c_str();

    TraMLFile traml;
    TargetedExperiment targeted_exp;

    TransitionTSVReader tsv_reader = TransitionTSVReader();
    std::cout << "Reading " << in << std::endl;
    tsv_reader.setLogType(log_type_);
    tsv_reader.convertTSVToTargetedExperiment(tr_file, targeted_exp);
    tsv_reader.validateTargetedExperiment(targeted_exp);

    std::cout << "Writing " << out << std::endl;
    traml.store(out, targeted_exp);

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{

  TOPPConvertTSVToTraML tool;
  return tool.main(argc, argv);
}

/// @endcond
