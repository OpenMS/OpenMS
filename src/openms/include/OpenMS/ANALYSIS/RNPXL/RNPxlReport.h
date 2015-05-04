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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_RNPXL_RNPXLREPORT
#define OPENMS_ANALYSIS_RNPXL_RNPXLREPORT

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/ANALYSIS/RNPXL/RNPxlMarkerIonExtractor.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>

namespace OpenMS
{

struct OPENMS_DLLAPI RNPxlReportRow
{
  bool no_id;
  double rt;
  double original_mz;
  String accessions;
  String RNA;
  String peptide;
  Int charge;
  double score;
  double peptide_weight;
  double RNA_weight;
  double xl_weight;
  double abs_prec_error;
  double rel_prec_error;
  RNPxlMarkerIonExtractor::MarkerIonsType marker_ions;
  double m_H;
  double m_2H;
  double m_3H;
  double m_4H;

  String getString(String separator) const;

};

struct OPENMS_DLLAPI RNPxlReportRowHeader
{
  String getString(String separator)
  {
    StringList sl;
    sl << "#RT" << "original m/z" << "proteins" << "RNA" << "peptide" << "charge" << "score"
       << "peptide weight" << "RNA weight" << "cross-link weight";

    // marker ion fields
    RNPxlMarkerIonExtractor::MarkerIonsType marker_ions = RNPxlMarkerIonExtractor::extractMarkerIons(PeakSpectrum(), 0.0); // call only to generate header entries
    for (RNPxlMarkerIonExtractor::MarkerIonsType::const_iterator it = marker_ions.begin(); it != marker_ions.end(); ++it)
    {
      for (Size i = 0; i != it->second.size(); ++i)
      {
        sl << String(it->first + "_" + it->second[i].first);
      }
    }
    sl << "abs prec. error Da" << "rel. prec. error ppm" << "M+H" << "M+2H" << "M+3H" << "M+4H";
    return ListUtils::concatenate(sl, separator);
  }
};

// create report
struct OPENMS_DLLAPI RNPxlReport
{

static std::vector<RNPxlReportRow> annotate(const PeakMap& spectra, std::vector<PeptideIdentification>& peptide_ids, double marker_ions_tolerance)
{
  std::map<Size, Size> map_spectra_to_id;
  for (Size i = 0; i != peptide_ids.size(); ++i)
  {
    OPENMS_PRECONDITION(!peptide_ids[i].getHits().empty(), "Error: no empty peptide ids allowed.");
    Size scan_index = (unsigned int)peptide_ids[i].getMetaValue("scan_index");
    map_spectra_to_id[scan_index] = i;
  }

  std::vector<RNPxlReportRow> csv_rows;

  for (PeakMap::ConstIterator s_it = spectra.begin(); s_it != spectra.end(); ++s_it)
  {
    int scan_index = s_it - spectra.begin();
    std::vector<Precursor> precursor = s_it->getPrecursors();

    // there should only one precursor and MS2 should contain at least a few peaks to be considered (e.g. at least for every AA in the peptide)
    if (s_it->getMSLevel() == 2 && precursor.size() == 1)
    {
      Size charge = precursor[0].getCharge();
      double mz = precursor[0].getMZ();
      RNPxlMarkerIonExtractor::MarkerIonsType marker_ions = RNPxlMarkerIonExtractor::extractMarkerIons(*s_it, marker_ions_tolerance);

      double rt = s_it->getRT();

      RNPxlReportRow row;

      // case 1: no peptide identification: store rt, mz, charge and marker ion intensities
      if (map_spectra_to_id.find(scan_index) == map_spectra_to_id.end())
      {
        row.no_id = true;
        row.rt = rt;
        row.original_mz = mz;
        row.charge = charge;
        row.marker_ions = marker_ions;
        csv_rows.push_back(row);
        continue;
      }

      PeptideIdentification& pi = peptide_ids[map_spectra_to_id[scan_index]];
      std::vector<PeptideHit>& phs = pi.getHits();

      // case 2: identification data present for spectrum
      PeptideHit& ph = phs[0];
      const AASequence& sequence = ph.getSequence();
      double peptide_weight = sequence.getMonoWeight();
      String rna_name = ph.getMetaValue("RNPxl:RNA");
      double rna_weight = ph.getMetaValue("RNPxl:RNA_MASS_z0");

      // crosslink weight for different charge states
      double weight_z1 = (peptide_weight + rna_weight + 1.0 * Constants::PROTON_MASS_U);
      double weight_z2 = (peptide_weight + rna_weight + 2.0 * Constants::PROTON_MASS_U) / 2.0;
      double weight_z3 = (peptide_weight + rna_weight + 3.0 * Constants::PROTON_MASS_U) / 3.0;
      double weight_z4 = (peptide_weight + rna_weight + 4.0 * Constants::PROTON_MASS_U) / 4.0;

      double xl_weight = peptide_weight + rna_weight;
      double theo_mz = (xl_weight + static_cast<double>(charge) * Constants::PROTON_MASS_U) / (double)charge;
      double absolute_difference = theo_mz - mz;
      double ppm_difference =  absolute_difference / theo_mz * 1e6;

      String protein_accessions;
      std::set<String> accs = ph.extractProteinAccessions();

      // concatenate set into String
      for (std::set<String>::const_iterator a_it = accs.begin(); a_it != accs.end(); ++a_it)
      {
        if (a_it != accs.begin())
        {
          protein_accessions += ",";
        }
        protein_accessions += *a_it;
      }

      row.no_id = false;
      row.rt = rt;
      row.original_mz = mz;
      row.accessions = protein_accessions;
      row.RNA = rna_name;
      row.peptide = sequence.toString();
      row.charge = charge;
      row.score = ph.getScore();
      row.peptide_weight = peptide_weight;
      row.RNA_weight = rna_weight;
      row.xl_weight = peptide_weight + rna_weight;

      ph.setMetaValue("RNPxl:peptide_mass_z0", DataValue(peptide_weight));
      ph.setMetaValue("RNPxl:xl_mass_z0", xl_weight);

      for (RNPxlMarkerIonExtractor::MarkerIonsType::const_iterator it = marker_ions.begin(); it != marker_ions.end(); ++it)
      {
        for (Size i = 0; i != it->second.size(); ++i)
        {
          ph.setMetaValue(it->first + "_" + it->second[i].first, static_cast<double>(it->second[i].second * 100.0));
        }
      }

      row.marker_ions = marker_ions;
      row.abs_prec_error = absolute_difference;
      row.rel_prec_error = ppm_difference;
      row.m_H = weight_z1;
      row.m_2H = weight_z2;
      row.m_3H = weight_z3;
      row.m_4H = weight_z4;

      ph.setMetaValue("RNPxl:Da difference", (double)absolute_difference);
      ph.setMetaValue("RNPxl:ppm difference", (double)ppm_difference);
      ph.setMetaValue("RNPxl:z1 mass", (double)weight_z1);
      ph.setMetaValue("RNPxl:z2 mass", (double)weight_z2);
      ph.setMetaValue("RNPxl:z3 mass", (double)weight_z3);
      ph.setMetaValue("RNPxl:z4 mass", (double)weight_z4);

      csv_rows.push_back(row);

    }
  }
  return csv_rows;
}
};

}

#endif

