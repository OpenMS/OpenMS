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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/ANALYSIS/RNPXL/RNPxlReport.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

namespace OpenMS
{
  String RNPxlReportRow::getString(const String& separator) const
  {
    StringList sl;

    // rt mz
    sl << String::number(rt, 3) << String::number(original_mz, 4);

    // id if available
    if (no_id)
    {
      sl << "" << "" << "" << "" << "" << "" << "" << "" << "" << "" << "" << "";
    }
    else
    {
      sl << accessions << RNA << peptide << String(charge) << String(score)
         << best_localization_score  << localization_scores << best_localization
         << String::number(peptide_weight, 4) << String::number(RNA_weight, 4) 
         << String::number(peptide_weight + RNA_weight, 4);
    }

    // marker ions
    for (RNPxlMarkerIonExtractor::MarkerIonsType::const_iterator it = marker_ions.begin(); it != marker_ions.end(); ++it)
    {
      for (Size i = 0; i != it->second.size(); ++i)
      {
        sl << String::number(it->second[i].second * 100.0, 2);
      }
    }

    // id error and multiple charged mass
    if (no_id)
    {
      sl << "" << "" << ""
         << "" << "" << "" << "";
    }
    else
    {
      // error
      sl << String::number(abs_prec_error, 4)
         << String::number(rel_prec_error, 1);

      // weight
      sl << String::number(m_H, 4)
         << String::number(m_2H, 4)
         << String::number(m_3H, 4)
         << String::number(m_4H, 4);
      // rank
      sl << String(rank);
    }

    return ListUtils::concatenate(sl, separator);
  }

  String RNPxlReportRowHeader::getString(const String& separator)
  {
    StringList sl;
    sl << "#RT" << "original m/z" << "proteins" << "RNA" << "peptide" << "charge" << "score"
       << "best localization score" << "localization scores" << "best localization(s)"
       << "peptide weight" << "RNA weight" << "cross-link weight";

    // marker ion fields
    RNPxlMarkerIonExtractor::MarkerIonsType marker_ions = RNPxlMarkerIonExtractor::extractMarkerIons(PeakSpectrum(), 0.0); // call only to generate header entries
    for (auto const & ma : marker_ions)
    {
      for (Size i = 0; i != ma.second.size(); ++i)
      {
        sl << String(ma.first + "_" + ma.second[i].first);
      }
    }
    sl << "abs prec. error Da" << "rel. prec. error ppm" << "M+H" << "M+2H" << "M+3H" << "M+4H" << "rank";
    return ListUtils::concatenate(sl, separator);
  }

  std::vector<RNPxlReportRow> RNPxlReport::annotate(const PeakMap& spectra, std::vector<PeptideIdentification>& peptide_ids, double marker_ions_tolerance)
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
        Size rank(0);
        for (PeptideHit& ph : phs)
        {
          ++rank;
          // total weight = precursor RNA weight + peptide weight
          // this ensures that sequences with additional reported partial loss match the total weight
          // Note that the partial loss is only relevent on the MS2 and would otherwise be added to the totalweight
          String sequence_string = ph.getSequence().toString();
          sequence_string.substitute("(RNA:U_prime-H2O)", "");
          sequence_string.substitute("(RNA:U_prime)", "");
          sequence_string.substitute("(RNA:U-H3PO4)", "");
          sequence_string.substitute("(RNA:U-H2O)", "");
          sequence_string.substitute("(RNA:U)", "");
          sequence_string.substitute("(RNA:C3O)", "");

          const AASequence sequence = AASequence::fromString(sequence_string);

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
          double ppm_difference =  Math::getPPM(mz, theo_mz);

          String protein_accessions;
          std::set<String> accs = ph.extractProteinAccessionsSet();

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
          row.peptide = ph.getSequence().toString();
          row.charge = charge;
          row.score = ph.getScore();
          row.peptide_weight = peptide_weight;
          row.RNA_weight = rna_weight;
          row.xl_weight = peptide_weight + rna_weight;
          row.rank = rank;

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

          if (ph.metaValueExists("RNPxl:best_localization_score") && 
              ph.metaValueExists("RNPxl:localization_scores") && 
              ph.metaValueExists("RNPxl:best_localization"))
          {
            row.best_localization_score = ph.getMetaValue("RNPxl:best_localization_score");
            row.localization_scores = ph.getMetaValue("RNPxl:localization_scores");
            row.best_localization = ph.getMetaValue("RNPxl:best_localization");;
          }

          ph.setMetaValue("RNPxl:Da difference", (double)absolute_difference);
          ph.setMetaValue(Constants::PRECURSOR_ERROR_PPM_USERPARAM, (double)ppm_difference);
          ph.setMetaValue("RNPxl:z1 mass", (double)weight_z1);
          ph.setMetaValue("RNPxl:z2 mass", (double)weight_z2);
          ph.setMetaValue("RNPxl:z3 mass", (double)weight_z3);
          ph.setMetaValue("RNPxl:z4 mass", (double)weight_z4);
          csv_rows.push_back(row);

          // In the last annotation step we add the oligo as delta mass modification
          // to get the proper theoretical mass annotated in the PeptideHit
          // Try to add it to the C- then N-terminus. 
          // If already modified search for an unmodified amino acid and add it there
          if (rna_weight > 0)
          {
            const AASequence& aa = ph.getSequence();
            const String seq = ph.getSequence().toString();

            if (!aa.hasCTerminalModification()) 
            {
              AASequence new_aa = AASequence::fromString(seq + ".[" + String(rna_weight) + "]");
              ph.setSequence(new_aa);
            }
            else if (!aa.hasNTerminalModification()) 
            {
              AASequence new_aa = AASequence::fromString("[" + String(rna_weight) + "]." + seq);
              ph.setSequence(new_aa);
            }
            else // place it anywhere
            {
              for (auto a : aa) // NOTE: performs copy, this cannot possibly work!
              {
                if (!a.isModified())
                {
                  a.setModification("[" + String(rna_weight) + "]");
                  break;
                }
              }             
            }
          }          
      }
    }
  }
  return csv_rows;
}
}
