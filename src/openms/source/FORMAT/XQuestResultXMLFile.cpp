// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Eugen Netz $
// $Authors: Lukas Zimmermann, Eugen Netz $
// --------------------------------------------------------------------------
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/XQuestResultXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XQuestResultXMLHandler.h>
#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <fstream>
#include <OpenMS/ANALYSIS/XLMS/OPXLHelper.h>

namespace OpenMS
{
  XQuestResultXMLFile::XQuestResultXMLFile() :
    XMLFile("/SCHEMAS/xQuest_1_0.xsd", "1.0"),
    n_hits_(-1)
  {
  }
  XQuestResultXMLFile::~XQuestResultXMLFile()
  {
  }

  void XQuestResultXMLFile::load(const String & filename,
                                 std::vector < PeptideIdentification > & pep_ids,
                                 std::vector< ProteinIdentification > & prot_ids
                                )
  {
   Internal::XQuestResultXMLHandler handler(filename, pep_ids, prot_ids);
   this->parse_(filename, &handler);

   this->n_hits_ = handler.getNumberOfHits();
   this->min_score_ = handler.getMinScore();
   this->max_score_ = handler.getMaxScore();

   // this helper function adds additional explicit "xl_target_decoy" meta values derived from parsed data
   OPXLHelper::addXLTargetDecoyMV(pep_ids);
   // this helper function adds beta peptide accessions
   OPXLHelper::addBetaAccessions(pep_ids);
   // this helper function bases the ranked lists of labeled XLMS searches on each light spectrum instead of pairs
   // the second parameter here should be the maximal number of hits per spectrum,
   // but using the total number of hits we will just keep everything contained in the file
   // (just reassigned to single spectra and reranked by score)
   pep_ids = OPXLHelper::combineTopRanksFromPairs(pep_ids, this->n_hits_);
  }

  int XQuestResultXMLFile::getNumberOfHits() const
  {
    return this->n_hits_;
  }

  double XQuestResultXMLFile::getMinScore() const
  {
    return this->min_score_;
  }

  double XQuestResultXMLFile::getMaxScore() const
  {
    return this->max_score_;
  }

  void XQuestResultXMLFile::store(const String& filename, const std::vector<ProteinIdentification>& poid, const std::vector<PeptideIdentification>& peid) const
  {
    if (!FileHandler::hasValidExtension(filename, FileTypes::XQUESTXML))
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "invalid file extension, expected '" + FileTypes::typeToName(FileTypes::XQUESTXML) + "'");
    }

    Internal::XQuestResultXMLHandler handler(poid, peid, filename, schema_version_);
    save_(filename, &handler);
  }

  // version for labeled linkers
  void XQuestResultXMLFile::writeXQuestXMLSpec(const String& out_file, const String& base_name, const OPXLDataStructs::PreprocessedPairSpectra& preprocessed_pair_spectra, const std::vector< std::pair<Size, Size> >& spectrum_pairs, const std::vector< std::vector< OPXLDataStructs::CrossLinkSpectrumMatch > >& all_top_csms, const PeakMap& spectra)
  {
    // XML Header
    std::ofstream spec_xml_file;
    std::cout << "Writing spec.xml to " << out_file << std::endl;
    spec_xml_file.open(out_file.c_str(), std::ios::trunc); // ios::app = append to file, ios::trunc = overwrites file
    // TODO write actual data
    spec_xml_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?><xquest_spectra author=\"Eugen Netz\" deffile=\"xquest.def\" >" << std::endl;

    // collect indices of spectra, that need to be written out
    std::vector <std::pair <Size, Size> > spectrum_indices;

    for (Size i = 0; i < all_top_csms.size(); ++i)
    {
      if (!all_top_csms[i].empty())
      {
        if (all_top_csms[i][0].scan_index_light < spectra.size() && all_top_csms[i][0].scan_index_heavy < spectra.size())
        {
          spectrum_indices.push_back( std::make_pair(all_top_csms[i][0].scan_index_light, all_top_csms[i][0].scan_index_heavy) );
        }
      }
    }

    // loop over list of indices and write out spectra
    for (Size i = 0; i < spectrum_indices.size(); ++i)
    {
      Size scan_index_light = spectrum_indices[i].first;
      Size scan_index_heavy = spectrum_indices[i].second;
      // TODO more correct alternative
      String spectrum_light_name = base_name + ".light." + scan_index_light;
      String spectrum_heavy_name = base_name + ".heavy." + scan_index_heavy;
      String spectrum_name = spectrum_light_name + String("_") + spectrum_heavy_name;

      if (scan_index_light < spectra.size() && scan_index_heavy < spectra.size() && i < preprocessed_pair_spectra.spectra_linear_peaks.size() && i < preprocessed_pair_spectra.spectra_xlink_peaks.size())
      {
        // 4 Spectra resulting from a light/heavy spectra pair.  Write for each spectrum, that is written to xquest.xml (should be all considered pairs, or better only those with at least one sensible Hit, meaning a score was computed)
        spec_xml_file << "<spectrum filename=\"" << spectrum_light_name << ".dta" << "\" type=\"light\">" << std::endl;
        spec_xml_file << getxQuestBase64EncodedSpectrum_(spectra[scan_index_light], String(""));
        spec_xml_file << "</spectrum>" << std::endl;

        spec_xml_file << "<spectrum filename=\"" << spectrum_heavy_name << ".dta" << "\" type=\"heavy\">" << std::endl;
        spec_xml_file << getxQuestBase64EncodedSpectrum_(spectra[scan_index_heavy], String(""));
        spec_xml_file << "</spectrum>" << std::endl;

        // the preprocessed pair spectra are sorted by another index
        // because some pairs do not yield any hits worth reporting (e.g. no matching peaks), the index from the spectrum matches or spectrum_indices does not address the right pair anymore
        // use find with the pair of spectrum indices to find the correct index for the preprocessed linear and cross-linked ion spectra
        std::vector<std::pair <Size, Size> >::const_iterator pair_it = std::find(spectrum_pairs.begin(), spectrum_pairs.end(), spectrum_indices[i]);
        Size pair_index = std::distance(spectrum_pairs.begin(), pair_it);

        String spectrum_common_name = spectrum_name + String("_common.txt");
        spec_xml_file << "<spectrum filename=\"" << spectrum_common_name << "\" type=\"common\">" << std::endl;
        spec_xml_file << getxQuestBase64EncodedSpectrum_(preprocessed_pair_spectra.spectra_linear_peaks[pair_index], spectrum_light_name + ".dta," + spectrum_heavy_name + ".dta");
        spec_xml_file << "</spectrum>" << std::endl;

        String spectrum_xlink_name = spectrum_name + String("_xlinker.txt");
        spec_xml_file << "<spectrum filename=\"" << spectrum_xlink_name << "\" type=\"xlinker\">" << std::endl;
        spec_xml_file << getxQuestBase64EncodedSpectrum_(preprocessed_pair_spectra.spectra_xlink_peaks[pair_index], spectrum_light_name + ".dta," + spectrum_heavy_name + ".dta");
        spec_xml_file << "</spectrum>" << std::endl;
      }
    }

    spec_xml_file << "</xquest_spectra>" << std::endl;
    spec_xml_file.close();

    return;
  }

  // version for label-free linkers
  void XQuestResultXMLFile::writeXQuestXMLSpec(const String& out_file, const String& base_name, const std::vector< std::vector< OPXLDataStructs::CrossLinkSpectrumMatch > >& all_top_csms, const PeakMap& spectra)
  {
    // String spec_xml_filename = base_name + "_matched.spec.xml";
    // XML Header
    std::ofstream spec_xml_file;
    std::cout << "Writing spec.xml to " << out_file << std::endl;
    spec_xml_file.open(out_file.c_str(), std::ios::trunc); // ios::app = append to file, ios::trunc = overwrites file
    // TODO write actual data
    spec_xml_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?><xquest_spectra author=\"Eugen Netz\" deffile=\"xquest.def\" >" << std::endl;

    // collect indices of spectra, that need to be written out
    std::vector <Size> spectrum_indices;

    for (Size i = 0; i < all_top_csms.size(); ++i)
    {
      if (!all_top_csms[i].empty())
      {
        if (all_top_csms[i][0].scan_index_light < spectra.size())
        {
          spectrum_indices.push_back(all_top_csms[i][0].scan_index_light);
        }
      }
    }

    // loop over list of indices and write out spectra
    for (Size i = 0; i < spectrum_indices.size(); ++i)
    {
      String spectrum_light_name = base_name + ".light." + spectrum_indices[i];
      String spectrum_heavy_name = base_name + ".heavy." + spectrum_indices[i];

      String spectrum_name = spectrum_light_name + String("_") + spectrum_heavy_name;

      // 4 Spectra resulting from a light/heavy spectra pair.  Write for each spectrum, that is written to xquest.xml (should be all considered pairs, or better only those with at least one sensible Hit, meaning a score was computed)
      spec_xml_file << "<spectrum filename=\"" << spectrum_light_name << ".dta" << "\" type=\"light\">" << std::endl;
      spec_xml_file << getxQuestBase64EncodedSpectrum_(spectra[spectrum_indices[i]], String(""));
      spec_xml_file << "</spectrum>" << std::endl;

      spec_xml_file << "<spectrum filename=\"" << spectrum_heavy_name << ".dta" << "\" type=\"heavy\">" << std::endl;
      spec_xml_file << getxQuestBase64EncodedSpectrum_(spectra[spectrum_indices[i]], String(""));
      spec_xml_file << "</spectrum>" << std::endl;

      String spectrum_common_name = spectrum_name + String("_common.txt");
      spec_xml_file << "<spectrum filename=\"" << spectrum_common_name << "\" type=\"common\">" << std::endl;
      spec_xml_file << getxQuestBase64EncodedSpectrum_(spectra[spectrum_indices[i]], spectrum_light_name + ".dta," + spectrum_heavy_name + ".dta");
      spec_xml_file << "</spectrum>" << std::endl;

      String spectrum_xlink_name = spectrum_name + String("_xlinker.txt");
      spec_xml_file << "<spectrum filename=\"" << spectrum_xlink_name << "\" type=\"xlinker\">" << std::endl;
      spec_xml_file << getxQuestBase64EncodedSpectrum_(spectra[spectrum_indices[i]], spectrum_light_name + ".dta," + spectrum_heavy_name + ".dta");
      spec_xml_file << "</spectrum>" << std::endl;
    }

    spec_xml_file << "</xquest_spectra>" << std::endl;
    spec_xml_file.close();

    return;
  }

  String XQuestResultXMLFile::getxQuestBase64EncodedSpectrum_(const PeakSpectrum& spec, String header)
  {
    std::vector<String> in_strings;
    StringList sl;

    double precursor_mz = 0;
    double precursor_z = 0;
    if (spec.getPrecursors().size() > 0)
    {
      precursor_mz = Math::roundDecimal(spec.getPrecursors()[0].getMZ(), -9);
      precursor_z = spec.getPrecursors()[0].getCharge();
    }

    // header lines
    if (!header.empty()) // common or xlinker spectrum will be reported
    {
      sl.push_back(header + "\n"); // e.g. GUA1372-S14-A-LRRK2_DSS_1A3.03873.03873.3.dta,GUA1372-S14-A-LRRK2_DSS_1A3.03863.03863.3.dta
      sl.push_back(String(precursor_mz) + "\n");
      sl.push_back(String(precursor_z) + "\n");
    }
    else // light or heavy spectrum will be reported
    {
      sl.push_back(String(precursor_mz) + "\t" + String(precursor_z) + "\n");
    }

    PeakSpectrum::IntegerDataArray charges;
    if (spec.getIntegerDataArrays().size() > 0)
    {
      charges = spec.getIntegerDataArrays()[0];
    }

    // write peaks
    for (Size i = 0; i != spec.size(); ++i)
    {
      String s;
      s += String(Math::roundDecimal(spec[i].getMZ(), -9)) + "\t";
      s += String(spec[i].getIntensity()) + "\t";

      if (charges.size() > 0)
      {
        s += String(charges[i]);
      }
      else
      {
        s += "0";
      }

      s += "\n";

      sl.push_back(s);
    }

    String out;
    out.concatenate(sl.begin(), sl.end(), "");
    in_strings.push_back(out);

    String out_encoded;
    Base64().encodeStrings(in_strings, out_encoded, false, false);
    String out_wrapped;
    wrap_(out_encoded, 76, out_wrapped);
    return out_wrapped;
  }

  void XQuestResultXMLFile::wrap_(const String& input, Size width, String & output)
  {
    Size start = 0;

    while (start + width < input.size())
    {
      output += input.substr(start, width) + "\n";
      start += width;
    }

    if (start < input.size())
    {
      output += input.substr(start, input.size() - start) + "\n";
    }
  }
}
