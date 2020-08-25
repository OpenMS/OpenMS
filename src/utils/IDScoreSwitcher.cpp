// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_IDScoreSwitcher IDScoreSwitcher

    @brief Switches between different scores of peptide hits (PSMs) or protein hits in identification data.

    In the idXML file format and in OpenMS' internal representation of identification data, every peptide spectrum match (PSM, "peptide hit") and every protein hit can be associated with a single numeric (quality) score of an arbitrary type. However, database search engines that generate PSMs or tools for post-processing of identification data may assign multiple scores of different types to each PSM/protein. These scores can be captured as meta data associated with the PSMs/proteins (in idXML: "UserParam" elements), but they are typically not considered by TOPP tools that utilize the scores. This utility allows to switch between "primary" scores and scores stored as meta values.

    By default this tool operates on PSM scores; to consider protein scores instead, set the @p proteins flag. The meta value that is supposed to replace the PSM/protein score - given by parameter @p new_score - has to be numeric (type "float") and exist for every peptide or protein hit, respectively. The old score will be stored as a meta value, the name for which is given by the parameter @p old_score. It is an error if a meta value with this name already exists for any hit, unless that meta value already stores the same score.

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_IDScoreSwitcher.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_IDScoreSwitcher.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDScoreSwitcher :
  public TOPPBase
{
public:

  TOPPIDScoreSwitcher() :
    TOPPBase("IDScoreSwitcher", "Switches between different scores of peptide or protein hits in identification data", false)
  {
  }

protected:

  IDScoreSwitcherAlgorithm switcher_{};

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerFlag_("proteins", "Apply to protein scores instead of PSM scores");
    registerFullParam_(switcher_.getParameters());
  }

  ExitCodes main_(int, const char**) override
  {
    switcher_.setParameters(getParam_().copySubset(switcher_.getParameters()));
    String in = getStringOption_("in"), out = getStringOption_("out");
    bool do_proteins_ = getFlag_("proteins");

    vector<ProteinIdentification> proteins;
    vector<PeptideIdentification> peptides;

    IdXMLFile().load(in, proteins, peptides);

    Size counter = 0;
    if (do_proteins_)
    {
      for (auto& pid : proteins)
      {
        switcher_.switchScores<ProteinIdentification>(pid, counter);
      }
    }
    else
    {
      for (auto& pepid : peptides)
      {
        switcher_.switchScores<PeptideIdentification>(pepid, counter);
      }
    }

    IdXMLFile().store(out, proteins, peptides);

    OPENMS_LOG_INFO << "Successfully switched " << counter << " "
             << (do_proteins_ ? "protein" : "PSM") << " scores." << endl;

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPIDScoreSwitcher tool;
  return tool.main(argc, argv);
}

/// @endcond
