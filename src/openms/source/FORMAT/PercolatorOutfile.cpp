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

#include <OpenMS/FORMAT/PercolatorOutfile.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/CsvFile.h>

#include <boost/math/special_functions/fpclassify.hpp> // for "isnan"


using namespace std;

namespace OpenMS
{

  PercolatorOutfile::PercolatorOutfile()
  {
    PSMInfoExtractor extractor;
    // MS-GF+ Percolator (mzid?) format:
    extractor.re.assign("_SII_(?<SCAN>\\d+)_\\d+_\\d+_(?<CHARGE>\\d+)_\\d+");
    extractor.count_from_zero = false;
    extractors_.push_back(extractor);
    // Mascot Percolator format (RT may be missing, e.g. for searches via
    // ProteomeDiscoverer):
    extractor.re.assign("spectrum:[^;]+[(scans:)(scan=)(spectrum=)](?<SCAN>\\d+)[^;]+;rt:(?<RT>\\d*(\\.\\d+)?);mz:(?<MZ>\\d+(\\.\\d+)?);charge:(?<CHARGE>-?\\d+)");
    extractor.count_from_zero = true;
    extractors_.push_back(extractor);
  }


  bool PercolatorOutfile::getPSMInfo_(const String& PSM_ID, 
                                      const vector<struct PSMInfoExtractor>& 
                                      extractors, Int& scan_number, Int& charge,
                                      double& rt, double& mz)
  {
    scan_number = -1;
    charge = 0;
    rt = numeric_limits<double>::quiet_NaN();
    mz = 0.0;

    vector<struct PSMInfoExtractor>::const_iterator ex_it = extractors.begin();
    try
    {
      for (; ex_it != extractors.end(); ++ex_it)
      {
        boost::smatch match;
        bool found = boost::regex_search(PSM_ID, match, ex_it->re);
        if (found)
        {
          if (match["RT"].matched)
          {
            String rt_val = match["RT"].str();
            if (!rt_val.empty()) rt = rt_val.toDouble();
          }
          if (match["MZ"].matched)
          {
            String mz_val = match["MZ"].str();
            if (!mz_val.empty()) mz = mz_val.toDouble();
          }
          if (match["CHARGE"].matched)
          {
            String charge_val = match["CHARGE"].str();
            if (!charge_val.empty()) charge = charge_val.toInt();
          }
          if (match["SCAN"].matched)
          {
            String scan_val = match["SCAN"].str();
            if (!scan_val.empty()) scan_number = scan_val.toInt();
            if (!ex_it->count_from_zero) --scan_number;
          }
          return true;
        }
      }
    }
    catch (Exception::ConversionError&)
    {
      String msg = "Error: PSM ID has unexpected format '" + PSM_ID + "'. The "
        "regular expression '" + ex_it->re.str() + "' matched, but the "
        "extracted information could not be converted to a number.";
      LOG_ERROR << msg << endl;
    }
    return false;
  }


  void PercolatorOutfile::getPeptideSequence_(String peptide, AASequence& seq)
    const
  {
    // 'peptide' includes neighboring amino acids, e.g.: K.AAAR.A
    // but unclear to which protein neighboring AAs belong, so we ignore them:
    size_t len = peptide.size(), start = 0, count = string::npos;
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
    string replacement = "(UniMod:$1)";
    peptide = boost::regex_replace(peptide, re, replacement);
    seq = AASequence::fromString(peptide);
  }


  void PercolatorOutfile::preprocessExperiment_(
    const MSExperiment<>& experiment, ScanInfoMap& scan_map)
  {
    Size index = 0;
    for (MSExperiment<>::ConstIterator it = experiment.begin(); 
         it != experiment.end(); ++it, ++index)
    {
      if (it->getMSLevel() > 1)
      {
        const Precursor& precursor = it->getPrecursors()[0];
        ScanInfo scan_info;
        scan_info.charge = precursor.getCharge();
        scan_info.rt = it->getRT();
        scan_info.mz = precursor.getMZ();
        scan_map[index] = scan_info;
      }
    }
  }


  void PercolatorOutfile::load(const String& filename,
                               ProteinIdentification& proteins, 
                               vector<PeptideIdentification>& peptides,
                               enum ScoreType output_score,
                               const String& psm_regex, bool count_from_zero,
                               const MSExperiment<>* experiment_p)
  {
    vector<struct PSMInfoExtractor> extractors;
    if (psm_regex.empty()) // use default regular expressions
    {
      extractors = extractors_;
    }
    else // use only user-supplied regular expression
    {
      struct PSMInfoExtractor extractor;
      extractor.re.assign(psm_regex);
      extractor.count_from_zero = count_from_zero;
      extractors.push_back(extractor);
    }

    ScanInfoMap scan_map;
    if (experiment_p) preprocessExperiment_(*experiment_p, scan_map);

    vector<String> items;
    CsvFile source(filename, '\t');
    source.getRow(0, items);
    String header = ListUtils::concatenate<String>(items, '\t');
    if (header != 
        "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds")
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                  header,
                                  "Not a valid header for Percolator output");
    }

    set<String> accessions;
    Size no_charge = 0, no_rt = 0, no_mz = 0; // counters for missing data
    peptides.clear();

    for (Size row = 1; row < source.rowCount(); ++row)
    {
      source.getRow(row, items);
      Int scan_number, charge;
      double rt, mz;
      bool success = getPSMInfo_(items[0], extractors, scan_number, charge,
                                 rt, mz);

      // can we look up additional information in the raw data?
      if (success && (scan_number >= 0) && !scan_map.empty())
      {
        ScanInfoMap::const_iterator pos = scan_map.find(Size(scan_number));
        if (pos == scan_map.end())
        {
          String msg = "Error: Could not find spectrum for scan number " + 
            String(scan_number) + ", extracted from row " + String(row);
          LOG_ERROR << msg << endl;
        }
        else
        {
          if (charge == 0) charge = pos->second.charge;
          if (boost::math::isnan(rt)) rt = pos->second.rt;
          if (mz == 0.0) mz = pos->second.mz;
        }
      }
      
      PeptideHit hit;
      if (charge != 0)
      {
        hit.setCharge(charge);
      }
      else
      {
        no_charge++; // maybe TODO: calculate charge from m/z and peptide mass?
      }

      PeptideIdentification peptide;
      peptide.setIdentifier("id");
      if (!boost::math::isnan(rt))
      {
        peptide.setRT(rt);
      }
      else
      {
        no_rt++;
      }
      if (mz != 0.0)
      {
        peptide.setMZ(mz);
      }
      else
      {
        no_mz++;
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
