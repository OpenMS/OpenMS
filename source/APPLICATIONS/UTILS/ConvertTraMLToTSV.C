// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
  @page UTILS_ConvertTraMLToTSV ConvertTraMLToTSV

  @brief Converts TraML files to OpenSWATH transition TSV files

  The OpenSWATH transition TSV files will have the following headers, all fields are separated by tabs:

        <ul>
          <li> PrecursorMz (float) </li>
          <li> ProductMz (float) </li>
          <li> Tr_calibrated (float) </li>
          <li> transition_name (free text, needs to be unique for each transition [in this file]) </li>
          <li> Collision Energy (float) </li>
          <li> LibraryIntensity (float) </li>
          <li> transition_group_id (free text, designates the transition group [e.g. peptide] to which this transition belongs) </li>
          <li> decoy (1==decoy, 0== no decoy; determines whether the transition is a decoy transition or not) </li>
          <li> PeptideSequence  (free text, sequence only (no modifications) ) </li>
          <li> ProteinName  (free text) </li>
          <li> Annotation  (free text, e.g. y7) </li>
          <li> FullUniModPeptideName  (free text, should contain modifications*)  </li>
          <li> PrecursorCharge (integer, contains the charge of the precursor) </li>
          <li> GroupLabel (free text, e.g. heavy or light) </li>
          <li> UniprotID (free text) </li>
          <li> FragmentType (free text, contains the type of the fragment, e.g. "b" or "y") </li>
          <li> FragmentCharge (integer, contains the fragment charge) </li>
          <li> FragmentSeriesNumber (integer, e.g. for y7 use "7" here) </li>
        </ul>

* modifications are returned in UniMod annotation.

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPConvertTraMLToTSV : public TOPPBase
{
public:

  TOPPConvertTraMLToTSV() :
  TOPPBase("ConvertTraMLToTSV", "Converts a TraML file to an OpenSWATH transition TSV file")
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input TraML file");
    setValidFormats_("in", StringList::create("TraML"));

    registerOutputFile_("out", "<file>", "", "Output OpenSWATH transition TSV file");
    setValidFormats_("out", StringList::create("csv"));
  }

  ExitCodes main_(int, const char **)
  {
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    const char * tr_file = out.c_str();

    TraMLFile traml;
    TargetedExperiment targeted_exp;

    std::cout << "Reading " << in << std::endl;
    traml.load(in, targeted_exp);
    TransitionTSVReader tsv_reader = TransitionTSVReader();
    tsv_reader.setLogType(log_type_);
    tsv_reader.convertTargetedExperimentToTSV(tr_file, targeted_exp);
    std::cout << "Writing " << out << std::endl;

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{

  TOPPConvertTraMLToTSV tool;
  return tool.main(argc, argv);
}

/// @endcond
