// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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

#include <OpenMS/FORMAT/PercolatorOutfile.h>

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/CsvFile.h>

#include <boost/math/special_functions/fpclassify.hpp> // for "isnan"
#include <boost/regex.hpp>

namespace OpenMS
{
  using namespace std;

  // initialize static variable:
  const std::string PercolatorOutfile::score_type_names[] =
    {"qvalue", "PEP", "score"};


  PercolatorOutfile::PercolatorOutfile()
  {
  }


  enum PercolatorOutfile::ScoreType PercolatorOutfile::getScoreType(
    String score_type_name)
  {
    score_type_name.toLower();
    if ((score_type_name == "q-value") || (score_type_name == "qvalue") ||
        (score_type_name == "q value"))
    {
      return QVALUE;
    }
    if ((score_type_name == "pep") ||
        (score_type_name == "posterior error probability"))
    {
      return POSTERRPROB;
    }
    if (score_type_name == "score")
    {
      return SCORE;
    }
    throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                  "Not a valid Percolator score type",
                                  score_type_name);
  }


  void PercolatorOutfile::getPeptideSequence_(String peptide, AASequence& seq)
    const
  {
    // 'peptide' includes neighboring amino acids, e.g.: K.AAAR.A
    // but unclear to which protein neighboring AAs belong, so we ignore them:
    size_t len = peptide.size(), start = 0, count = std::string::npos;
    if (peptide[1] == '.') start = 2;
    if (peptide[len - 2] == '.') count = len - start - 2;
    peptide = peptide.substr(start, count);

    // re-format modifications:
    String unknown_mod = "[unknown]";
    if (peptide.hasSubstring(unknown_mod))
    {
      LOG_WARN << "Removing unknown modification(s) from peptide '" << peptide
               << "'" << endl;
      peptide.substitute(unknown_mod, "");
    }
    boost::regex re("\\[UNIMOD:(\\d+)\\]");
    std::string replacement = "(UniMod:$1)";
    peptide = boost::regex_replace(peptide, re, replacement);
    seq = AASequence::fromString(peptide);
  }


  String PercolatorOutfile::getFullModName_(const String& residue,
                                            const String& mod) const
  {
    if (residue == "N-term")
    {
      try
      {
        const ResidueModification& n_term_mod =
          ModificationsDB::getInstance()->
          getModification(mod, "", ResidueModification::N_TERM);
        return n_term_mod.getFullId();
      }
      catch (Exception::ElementNotFound) {};
    }
    else if (residue == "C-term")
    {
      try
      {
        const ResidueModification& c_term_mod =
          ModificationsDB::getInstance()->
          getModification(mod, "", ResidueModification::C_TERM);
        return c_term_mod.getFullId();
      }
      catch (Exception::ElementNotFound) {};
    }

    return mod + " (" + residue + ")";
  }


  void PercolatorOutfile::addModsToSearchParams_(
    const map<String, set<String> >::const_iterator& map_it,
    ProteinIdentification::SearchParameters& params) const
  {
    set<String>::const_iterator set_it = map_it->second.begin();

    // if there's only one mod, it's probably a fixed one:
    if ((map_it->second.size() == 1) && (!set_it->empty()))
    {
      params.fixed_modifications.push_back(getFullModName_(map_it->first,
                                                           *set_it));
    }
    else // variable mod(s)
    {
      for (; set_it != map_it->second.end(); ++set_it)
      {
        if (!set_it->empty())
        {
          params.variable_modifications.push_back(getFullModName_(map_it->first,
                                                                  *set_it));
        }
      }
    }
  }


  void PercolatorOutfile::getSearchModifications_(
    const vector<PeptideIdentification>& peptides,
    ProteinIdentification::SearchParameters& params) const
  {
    // amino acid (or terminus) -> set of modifications (incl. no mod.):
    map<String, set<String> > mod_map;

    for (vector<PeptideIdentification>::const_iterator pep_it =
           peptides.begin(); pep_it != peptides.end(); ++pep_it)
    {
      for (vector<PeptideHit>::const_iterator hit_it =
             pep_it->getHits().begin(); hit_it != pep_it->getHits().end();
           ++hit_it)
      {
        const String& n_term_mod =
          hit_it->getSequence().getNTerminalModificationName();
        mod_map["N-term"].insert(n_term_mod);
        const String& c_term_mod =
          hit_it->getSequence().getCTerminalModificationName();
        mod_map["C-term"].insert(c_term_mod);

        for (AASequence::ConstIterator seq_it = hit_it->getSequence().begin();
             seq_it != hit_it->getSequence().end(); ++seq_it)
        {
          const String& mod = seq_it->getModificationName();
          mod_map[seq_it->getOneLetterCode()].insert(mod);
        }
      }
    }

    for (map<String, set<String> >::const_iterator map_it =
           mod_map.begin(); map_it != mod_map.end(); ++map_it)
    {
      addModsToSearchParams_(map_it, params);
    }
    sort(params.fixed_modifications.begin(), params.fixed_modifications.end());
    sort(params.variable_modifications.begin(),
         params.variable_modifications.end());
  }


  void PercolatorOutfile::load(const String& filename,
                               ProteinIdentification& proteins,
                               vector<PeptideIdentification>& peptides,
                               SpectrumMetaDataLookup& lookup,
                               enum ScoreType output_score)
  {
    SpectrumMetaDataLookup::MetaDataFlags lookup_flags =
      (SpectrumMetaDataLookup::MDF_RT |
       SpectrumMetaDataLookup::MDF_PRECURSORMZ |
       SpectrumMetaDataLookup::MDF_PRECURSORCHARGE);

    if (lookup.reference_formats.empty())
    {
      // MS-GF+ Percolator (mzid?) format:
      lookup.addReferenceFormat("_SII_(?<INDEX1>\\d+)_\\d+_\\d+_(?<CHARGE>\\d+)_\\d+");
      // Mascot Percolator format (RT may be missing, e.g. for searches via
      // ProteomeDiscoverer):
      lookup.addReferenceFormat("spectrum:[^;]+[(scans:)(scan=)(spectrum=)](?<INDEX0>\\d+)[^;]+;rt:(?<RT>\\d*(\\.\\d+)?);mz:(?<MZ>\\d+(\\.\\d+)?);charge:(?<CHARGE>-?\\d+)");
    }

    vector<String> items;
    CsvFile source(filename, '\t');
    source.getRow(0, items);
    String header = ListUtils::concatenate<String>(items, '\t');
    if (header != 
        "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds")
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                  header, "Not a valid header for Percolator (PSM level) output");
    }

    set<String> accessions;
    Size no_charge = 0, no_rt = 0, no_mz = 0; // counters for missing data
    peptides.clear();

    for (Size row = 1; row < source.rowCount(); ++row)
    {
      source.getRow(row, items);

      SpectrumMetaDataLookup::SpectrumMetaData meta_data;
      try
      {
        lookup.getSpectrumMetaData(items[0], meta_data, lookup_flags);
      }
      catch (...)
      {
        String msg = "Error: Could not extract data for spectrum reference '" +
          items[0] + "' from row " + String(row);
        LOG_ERROR << msg << endl;
      }

      PeptideHit hit;
      if (meta_data.precursor_charge != 0)
      {
        hit.setCharge(meta_data.precursor_charge);
      }
      else
      {
        ++no_charge; // maybe TODO: calculate charge from m/z and peptide mass?
      }

      PeptideIdentification peptide;
      peptide.setIdentifier("id");
      if (!boost::math::isnan(meta_data.rt))
      {
        peptide.setRT(meta_data.rt);
      }
      else
      {
        ++no_rt;
      }
      if (!boost::math::isnan(meta_data.precursor_mz))
      {
        peptide.setMZ(meta_data.precursor_mz);
      }
      else
      {
        ++no_mz;
      }

      double score = items[1].toDouble();
      double qvalue = items[2].toDouble();
      double posterrprob = items[3].toDouble();
      hit.setMetaValue("Percolator_score", score);
      hit.setMetaValue("Percolator_qvalue", qvalue);
      hit.setMetaValue("Percolator_PEP", posterrprob);
      switch (output_score)
      {
        case SCORE:
          hit.setScore(score);
          peptide.setScoreType("Percolator_score");
          peptide.setHigherScoreBetter(true);
          break;
        case QVALUE:
          hit.setScore(qvalue);
          peptide.setScoreType("q-value");
          peptide.setHigherScoreBetter(false);
          break;
        case POSTERRPROB:
          hit.setScore(posterrprob);
          peptide.setScoreType("Posterior Error Probability");
          peptide.setHigherScoreBetter(false);
          break;
        case SIZE_OF_SCORETYPE:
          throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "'output_score' must not be 'SIZE_OF_SCORETYPE'!");
      }

      AASequence seq;
      getPeptideSequence_(items[4], seq);
      hit.setSequence(seq);

      for (Size pos = 5; pos < items.size(); ++pos)
      {
        accessions.insert(items[pos]);
        PeptideEvidence evidence;
        evidence.setProteinAccession(items[pos]);
        hit.addPeptideEvidence(evidence);
      }

      peptide.insertHit(hit);
      peptides.push_back(peptide);
    }

    proteins = ProteinIdentification();
    proteins.setIdentifier("id");
    proteins.setDateTime(DateTime::now());
    proteins.setSearchEngine("Percolator");

    for (set<String>::const_iterator it = accessions.begin();
         it != accessions.end(); ++it)
    {
      ProteinHit hit;
      hit.setAccession(*it);
      proteins.insertHit(hit);
    }

    // add info about allowed modifications:
    ProteinIdentification::SearchParameters params;
    getSearchModifications_(peptides, params);
    proteins.setSearchParameters(params);

    LOG_INFO << "Created " << proteins.getHits().size() << " protein hits.\n"
             << "Created " << peptides.size() << " peptide hits (PSMs)."
             << endl;
    if (no_charge > 0)
    {
      LOG_WARN << no_charge << " peptide hits without charge state information."
               << endl;
    }
    if (no_rt > 0)
    {
      LOG_WARN << no_rt << " peptide hits without retention time information." 
               << endl;
    }
    if (no_mz > 0)
    {
      LOG_WARN << no_mz << " peptide hits without mass-to-charge information." 
               << endl;
    }
  }

} // namespace OpenMS
