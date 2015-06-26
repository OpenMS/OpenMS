// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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

    @brief Switches between different scores of PSMs (peptide hits) in identification data.

    In the idXML file format and in OpenMS' internal representation of identification data, every peptide spectrum match (PSM, "peptide hit") can be associated with a single numeric (quality) score of an arbitrary type. However, database search engines that generate PSMs or tools for post-processing of identification data may assign multiple scores of different types to each PSM. These scores can be captured as meta data associated with the PSMs, but they are typically not considered by TOPP tools that utilize the PSM scores. This utility allows to switch between "primary" scores and scores stored as meta values.

    The meta value that is supposed to replace the PSM score - given by parameter @p new_score - has to be numeric (type "float") and exist for every PSM. The old PSM score will be stored as a meta value, the name for which is given by the parameter @p old_score. It is an error if a meta value with this name already exists for any PSM, unless that meta value already stores the same score.

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
    TOPPBase("IDScoreSwitcher", "Switches between different scores of PSMs (peptide hits) in identification data", false)
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input file");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>("idXML"));

    registerStringOption_("new_score", "<name>", "", "Name of the meta value to use as the new score");
    registerStringOption_("new_score_orientation", "<choice>", "", "Orientation of the new score (are higher or lower values better?)");
    setValidStrings_("new_score_orientation", ListUtils::create<String>("lower_better,higher_better"));
    registerStringOption_("new_score_type", "<name>", "", "Name to use as the type of the new score (default: same as 'new_score')", false);
    registerStringOption_("old_score", "<name>", "", "Name to use for the meta value storing the old score (default: old score type)", false);
  }

  String describePeptideHit(const PeptideHit& hit)
  {
    return "peptide hit with sequence '" + hit.getSequence().toString() +
      "', charge " + String(hit.getCharge()) + ", score " + 
      String(hit.getScore());
  }


  ExitCodes main_(int, const char**)
  {
    String in = getStringOption_("in"), out = getStringOption_("out"),
      new_score = getStringOption_("new_score"),
      new_score_type = getStringOption_("new_score_type"),
      old_score = getStringOption_("old_score");
    bool higher_better = (getStringOption_("new_score_orientation") == 
                          "higher_better");

    if (new_score_type.empty()) new_score_type = new_score;

    const double tolerance = 1e-6; // relative tolerance

    vector<ProteinIdentification> proteins;
    vector<PeptideIdentification> peptides;

    IdXMLFile().load(in, proteins, peptides);

    Size counter = 0;
    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      for (vector<PeptideHit>::iterator hit_it = pep_it->getHits().begin();
           hit_it != pep_it->getHits().end(); ++hit_it, ++counter)
      {
        if (!hit_it->metaValueExists(new_score))
        {
          String msg = "Meta value '" + new_score + "' not found for " + 
            describePeptideHit(*hit_it);
          throw Exception::MissingInformation(__FILE__, __LINE__, __PRETTY_FUNCTION__, msg);
        }

        String old_score_meta = (old_score.empty() ? pep_it->getScoreType() :
                                 old_score);
        DataValue dv = hit_it->getMetaValue(old_score_meta);
        if (!dv.isEmpty()) // meta value for old score already exists
        {
          if (fabs((double(dv) - hit_it->getScore()) * 2.0 /
                   (double(dv) + hit_it->getScore())) > tolerance)
          {
            String msg = "Meta value '" + old_score_meta + "' already exists "
              "with a conflicting value for " + describePeptideHit(*hit_it);
            throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, msg, dv.toString());
          } // else: values match, nothing to do
        }
        else
        {
          hit_it->setMetaValue(old_score_meta, hit_it->getScore());
        }
        hit_it->setScore(hit_it->getMetaValue(new_score));
      }
      pep_it->setScoreType(new_score_type);
      pep_it->setHigherScoreBetter(higher_better);
    }

    IdXMLFile().store(out, proteins, peptides);

    LOG_INFO << "Successfully switched " << counter << " PSM scores." << endl;

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPIDScoreSwitcher tool;
  return tool.main(argc, argv);
}

/// @endcond
