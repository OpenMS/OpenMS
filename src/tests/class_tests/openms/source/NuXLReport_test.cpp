// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/NUXL/NuXLModificationsGenerator.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLReport.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CONCEPT/Constants.h>

using namespace OpenMS;
using namespace std;

START_TEST(NuXLReport, "$Id$")

START_SECTION((static std::vector<NuXLReportRow> annotate(const PeakMap&, PeptideIdentificationList&, const StringList&, double)))
{
  const double adduct_mass = 306.025304840900048;
  const Int charge = 2;
  const AASequence peptide = AASequence::fromString("M(Oxidation)PEPTIDE");
  const double peptide_weight = peptide.getMonoWeight();
  const double xl_weight = peptide_weight + adduct_mass;
  const double precursor_mz = (xl_weight + charge * Constants::PROTON_MASS_U) / charge;

  AASequence peptide_with_adduct = peptide;
  peptide_with_adduct.setModificationByDiffMonoMass(0, adduct_mass);

  // The named definition on a clean residue, and folded into the oxidised one.
  const EmpiricalFormula adduct_formula("C9H11N2O8P1");
  AASequence peptide_with_definition = peptide;
  peptide_with_definition.setModification(4, NuXLModificationsGenerator::registerPrecursorAdduct("U-H2O1", adduct_formula, peptide[4]));
  TEST_STRING_EQUAL(peptide_with_definition.toString(), "M(Oxidation)PEPT(NuXL:U-H2O1)IDE")
  AASequence peptide_with_combined_definition = peptide;
  peptide_with_combined_definition.setModification(0, NuXLModificationsGenerator::registerPrecursorAdduct("U-H2O1", adduct_formula, peptide[0]));
  TEST_STRING_EQUAL(peptide_with_combined_definition.toString(), "M(NuXL:U-H2O1~Oxidation)PEPTIDE")

  PeakMap spectra;
  PeptideIdentificationList peptide_ids;

  auto add_psm = [&](const AASequence& sequence, const Int localization_position)
  {
    MSSpectrum spectrum;
    spectrum.setMSLevel(2);
    Precursor precursor;
    precursor.setCharge(charge);
    precursor.setMZ(precursor_mz);
    spectrum.setPrecursors({precursor});
    spectra.addSpectrum(spectrum);

    PeptideHit hit;
    hit.setCharge(charge);
    hit.setSequence(sequence);
    hit.setMetaValue("CalcMass", precursor_mz);
    hit.setMetaValue("NuXL:NA", "U-H2O1");
    hit.setMetaValue("NuXL:NA_MASS_z0", adduct_mass);
    hit.setMetaValue("NuXL:isXL", 1);
    hit.setMetaValue("NuXL:best_localization_position", localization_position);
    hit.setMetaValue("isotope_error", 0);

    PeptideIdentification peptide_id;
    peptide_id.setMetaValue("scan_index", static_cast<unsigned int>(peptide_ids.size()));
    peptide_id.setHits({hit});
    peptide_ids.push_back(peptide_id);
  };

  // Legacy localized results carry the adduct only in CalcMass.
  add_psm(peptide, 0);
  // New localized results additionally carry the adduct in the AASequence.
  add_psm(peptide_with_adduct, 0);
  // Unlocalized results intentionally retain the legacy representation.
  add_psm(peptide, -1);
  // Named definitions must yield the same masses as the mass-delta form.
  add_psm(peptide_with_definition, 4);
  add_psm(peptide_with_combined_definition, 0);

  const vector<NuXLReportRow> rows = NuXLReport::annotate(spectra, peptide_ids, {}, 0.05);
  TEST_EQUAL(rows.size(), 5)

  for (Size i = 0; i != rows.size(); ++i)
  {
    TEST_REAL_SIMILAR(rows[i].peptide_weight, peptide_weight)
    TEST_REAL_SIMILAR(rows[i].xl_weight, xl_weight)
    TEST_REAL_SIMILAR(rows[i].m_2H, precursor_mz)
    TEST_REAL_SIMILAR(rows[i].abs_prec_error, 0.0)

    const PeptideHit& annotated_hit = peptide_ids[i].getHits()[0];
    TEST_REAL_SIMILAR(static_cast<double>(annotated_hit.getMetaValue("NuXL:peptide_mass_z0")), peptide_weight)
    TEST_REAL_SIMILAR(static_cast<double>(annotated_hit.getMetaValue("NuXL:xl_mass_z0")), xl_weight)
    TEST_REAL_SIMILAR(static_cast<double>(annotated_hit.getMetaValue("NuXL:z2 mass")), precursor_mz)
  }

  TEST_STRING_EQUAL(rows[0].peptide, peptide.toString())
  TEST_STRING_EQUAL(rows[1].peptide, peptide_with_adduct.toString())
  TEST_STRING_EQUAL(rows[2].peptide, peptide.toString())
  TEST_STRING_EQUAL(rows[3].peptide, peptide_with_definition.toString())
  TEST_STRING_EQUAL(rows[4].peptide, peptide_with_combined_definition.toString())
}
END_SECTION

END_TEST
