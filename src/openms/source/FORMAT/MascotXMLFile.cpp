// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MascotXMLFile.h>

using namespace xercesc;
using namespace std;

namespace OpenMS
{

  MascotXMLFile::MascotXMLFile() :
    Internal::XMLFile()
  {
  }

  void MascotXMLFile::load(const String& filename,
                           ProteinIdentification& protein_identification,
                           vector<PeptideIdentification>& id_data,
                           const SpectrumMetaDataLookup& lookup)
  {
    map<String, vector<AASequence> > peptides;

    load(filename, protein_identification, id_data, peptides, lookup);
  }

  void MascotXMLFile::load(const String& filename,
                           ProteinIdentification& protein_identification,
                           vector<PeptideIdentification>& id_data,
                           map<String, vector<AASequence> >& peptides,
                           const SpectrumMetaDataLookup& lookup)
  {
    //clear
    protein_identification = ProteinIdentification();
    id_data.clear();

    Internal::MascotXMLHandler handler(protein_identification, id_data, 
                                       filename, peptides, lookup);
    parse_(filename, &handler);

    // since the Mascot XML can contain "peptides" without sequences,
    // the identifications without any real peptide hit are removed
    vector<PeptideIdentification> filtered_hits;
    filtered_hits.reserve(id_data.size());
    Size missing_sequence = 0; // counter

    for (PeptideIdentification& id_it : id_data)
    {
      const vector<PeptideHit>& peptide_hits = id_it.getHits();
      if (!peptide_hits.empty() && 
          (peptide_hits.size() > 1 || !peptide_hits[0].getSequence().empty()))
      {
        filtered_hits.push_back(id_it);
      }
      else if (!id_it.empty()) ++missing_sequence;
    }
    if (missing_sequence) 
    {
      OPENMS_LOG_WARN << "Warning: Removed " << missing_sequence 
               << " peptide identifications without sequence." << endl;
    }
    id_data.swap(filtered_hits);

    // check if we have (some) RT information:
    Size no_rt_count = 0;
    for (PeptideIdentification& id_it : id_data)
    {
      if (!id_it.hasRT())
      {
        ++no_rt_count;
      }
    }
    if (no_rt_count)
    {
      OPENMS_LOG_WARN << "Warning: " << no_rt_count << " (of " << id_data.size() 
               << ") peptide identifications have no retention time value."
               << endl;
    }
    // if we have a mapping, but couldn't find any RT values, that's an error:
    if (!lookup.empty() && (no_rt_count == id_data.size()))
    {
      throw Exception::MissingInformation(
        __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
        "No retention time information for peptide identifications found");
    }

    // argh! Mascot 2.2 tends to repeat the first hit (yes it appears twice),
    // so we delete one of them
    for (PeptideIdentification& pip : id_data)
    {
      vector<PeptideHit> peptide_hits = pip.getHits();
      // check if equal, except for rank
      if (peptide_hits.size() > 1 &&
          peptide_hits[0].getScore() == peptide_hits[1].getScore() &&
          peptide_hits[0].getSequence() == peptide_hits[1].getSequence() &&
          peptide_hits[0].getCharge() == peptide_hits[1].getCharge())
      {
        // erase first hit
        peptide_hits.erase(peptide_hits.begin() + 1);
        pip.setHits(peptide_hits);
      }
    }
  }


  void MascotXMLFile::initializeLookup(SpectrumMetaDataLookup& lookup, const PeakMap& exp, const String& scan_regex)
  {
    // load spectra and extract scan numbers from the native IDs
    // (expected format: "... scan=#"):
    lookup.readSpectra(exp.getSpectra());
    if (scan_regex.empty()) // use default formats
    {
      if (!lookup.empty()) // raw data given -> spectrum look-up possible
      {
        // possible formats and resulting scan numbers:
        // - Mascot 2.3 (?):
        // <pep_scan_title>scan=818</pep_scan_title> -> 818
        // - ProteomeDiscoverer/Mascot 2.3 or 2.4:
        // <pep_scan_title>Spectrum136 scans:712,</pep_scan_title> -> 712
        // - other variants:
        // <pep_scan_title>Spectrum3411 scans: 2975,</pep_scan_title> -> 2975
        // <...>File773 Spectrum198145 scans: 6094</...> -> 6094
        // <...>6860: Scan 10668 (rt=5380.57)</...> -> 10668
        // <pep_scan_title>Scan Number: 1460</pep_scan_title> -> 1460
        lookup.addReferenceFormat("[Ss]can( [Nn]umber)?s?[=:]? *(?<SCAN>\\d+)");
        // - with .dta input to Mascot:
        // <...>/path/to/FTAC05_13.673.673.2.dta</...> -> 673
        lookup.addReferenceFormat(R"(\.(?<SCAN>\d+)\.\d+\.(?<CHARGE>\d+)(\.dta)?)");
      }
      // title containing RT and MZ instead of scan number:
      // <...>575.848571777344_5018.0811_controllerType=0 controllerNumber=1 scan=11515_EcoliMS2small</...>
      lookup.addReferenceFormat(R"(^(?<MZ>\d+(\.\d+)?)_(?<RT>\d+(\.\d+)?))");
    }
    else // use only user-defined format
    {
      lookup.addReferenceFormat(scan_regex);
    }
  }

} // namespace OpenMS
