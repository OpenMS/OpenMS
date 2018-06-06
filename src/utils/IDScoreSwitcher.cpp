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
    TOPPBase("IDScoreSwitcher", "Switches between different scores of peptide or protein hits in identification data", false), tolerance_(1e-6)
  {
  }

protected:

  /// relative tolerance for score comparisons:
  const double tolerance_;

  String new_score_, new_score_type_, old_score_; // tool parameters
  bool higher_better_; // for the new scores, are higher ones better?


  void registerOptionsAndFlags_() override
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
    registerFlag_("proteins", "Apply to protein scores instead of PSM scores");
  }


  String describeHit_(const PeptideHit& hit)
  {
    return "peptide hit with sequence '" + hit.getSequence().toString() +
      "', charge " + String(hit.getCharge()) + ", score " + 
      String(hit.getScore());
  }


  String describeHit_(const ProteinHit& hit)
  {
    return "protein hit with accession '" + hit.getAccession() + "', score " +
      String(hit.getScore());
  }


  template <typename IDType>
  void switchScores_(IDType& id, Size& counter)
  {
    for (typename vector<typename IDType::HitType>::iterator hit_it = id.getHits().begin();
         hit_it != id.getHits().end(); ++hit_it, ++counter)
    {
      if (!hit_it->metaValueExists(new_score_))
      {
        String msg = "Meta value '" + new_score_ + "' not found for " + 
          describeHit_(*hit_it);
        throw Exception::MissingInformation(__FILE__, __LINE__,
                                            OPENMS_PRETTY_FUNCTION, msg);
      }

      String old_score_meta = (old_score_.empty() ? id.getScoreType() : 
                               old_score_);
      DataValue dv = hit_it->getMetaValue(old_score_meta);
      if (!dv.isEmpty()) // meta value for old score already exists
      {
        if (fabs((double(dv) - hit_it->getScore()) * 2.0 /
                 (double(dv) + hit_it->getScore())) > tolerance_)
        {
          String msg = "Meta value '" + old_score_meta + "' already exists "
            "with a conflicting value for " + describeHit_(*hit_it);
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        msg, dv.toString());
        } // else: values match, nothing to do
      }
      else
      {
        hit_it->setMetaValue(old_score_meta, hit_it->getScore());
      }
      hit_it->setScore(hit_it->getMetaValue(new_score_));
    }
    id.setScoreType(new_score_type_);
    id.setHigherScoreBetter(higher_better_);
  }


  ExitCodes main_(int, const char**) override
  {
    String in = getStringOption_("in"), out = getStringOption_("out");
    bool do_proteins = getFlag_("proteins");
    new_score_ = getStringOption_("new_score");
    new_score_type_ = getStringOption_("new_score_type");
    old_score_ = getStringOption_("old_score");
    higher_better_ = (getStringOption_("new_score_orientation") == 
                      "higher_better");

    if (new_score_type_.empty()) new_score_type_ = new_score_;

    vector<ProteinIdentification> proteins;
    vector<PeptideIdentification> peptides;

    IdXMLFile().load(in, proteins, peptides);

    Size counter = 0;
    if (do_proteins)
    {
      for (vector<ProteinIdentification>::iterator prot_it = proteins.begin();
           prot_it != proteins.end(); ++prot_it)
      {
        switchScores_<ProteinIdentification>(*prot_it, counter);
      }
    }
    else
    {
      for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
           pep_it != peptides.end(); ++pep_it)
      {
        switchScores_<PeptideIdentification>(*pep_it, counter);
      }
    }

    IdXMLFile().store(out, proteins, peptides);

    LOG_INFO << "Successfully switched " << counter << " "
             << (do_proteins ? "protein" : "PSM") << " scores." << endl;

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPIDScoreSwitcher tool;
  return tool.main(argc, argv);
}

/// @endcond
