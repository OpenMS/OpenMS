// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: OpenMS Team $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/IsoelectricPoint.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <cmath>

namespace OpenMS
{
  const IsoelectricPoint::PkaValues& IsoelectricPoint::getPkaValues_(ProteomicsPkaScale scale)
  {
    // Lehninger: Nelson DL, Cox MM. Lehninger Principles of Biochemistry.
    static const PkaValues lehninger = {9.69, 2.34, 3.65, 4.25, 8.18, 10.07, 6.00, 10.53, 12.48};
    // EMBOSS: used by pepstats (EMBOSS suite)
    static const PkaValues emboss = {8.6, 3.6, 3.9, 4.1, 8.5, 10.1, 6.5, 10.8, 12.5};
    // Sillero: Sillero & Ribeiro, Anal Biochem 1989
    static const PkaValues sillero = {8.2, 3.2, 4.0, 4.5, 9.0, 10.0, 6.4, 10.4, 12.0};

    switch (scale)
    {
      case ProteomicsPkaScale::LEHNINGER:
        return lehninger;
      case ProteomicsPkaScale::EMBOSS:
        return emboss;
      case ProteomicsPkaScale::SILLERO:
        return sillero;
      default:
        return lehninger;
    }
  }

  double IsoelectricPoint::chargeAcidic_(double pH, double pKa)
  {
    // For an acidic group: charge = -1 / (1 + 10^(pKa - pH))
    return -1.0 / (1.0 + std::pow(10.0, pKa - pH));
  }

  double IsoelectricPoint::chargeBasic_(double pH, double pKa)
  {
    // For a basic group: charge = +1 / (1 + 10^(pH - pKa))
    return 1.0 / (1.0 + std::pow(10.0, pH - pKa));
  }

  double IsoelectricPoint::computeCharge(
    const AASequence& seq,
    double pH,
    ProteomicsPkaScale scale)
  {
    if (seq.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Cannot compute charge of an empty sequence", "");
    }

    const PkaValues& pka = getPkaValues_(scale);

    // Start with terminal contributions
    double charge = chargeBasic_(pH, pka.nterm) + chargeAcidic_(pH, pka.cterm);

    // Add side chain contributions
    for (const auto& residue : seq)
    {
      const String& olc = residue.getOneLetterCode();
      if (olc == "D")
      {
        charge += chargeAcidic_(pH, pka.D);
      }
      else if (olc == "E")
      {
        charge += chargeAcidic_(pH, pka.E);
      }
      else if (olc == "C")
      {
        charge += chargeAcidic_(pH, pka.C);
      }
      else if (olc == "Y")
      {
        charge += chargeAcidic_(pH, pka.Y);
      }
      else if (olc == "H")
      {
        charge += chargeBasic_(pH, pka.H);
      }
      else if (olc == "K")
      {
        charge += chargeBasic_(pH, pka.K);
      }
      else if (olc == "R")
      {
        charge += chargeBasic_(pH, pka.R);
      }
    }

    return charge;
  }

  double IsoelectricPoint::computePI(
    const AASequence& seq,
    ProteomicsPkaScale scale,
    double tolerance)
  {
    if (seq.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Cannot compute pI of an empty sequence", "");
    }

    // Bisection method over pH range [0, 14]
    double pH_low = 0.0;
    double pH_high = 14.0;
    double pH_mid = 7.0;

    while ((pH_high - pH_low) > tolerance)
    {
      pH_mid = (pH_low + pH_high) / 2.0;
      double charge = computeCharge(seq, pH_mid, scale);

      if (charge > 0.0)
      {
        pH_low = pH_mid;
      }
      else
      {
        pH_high = pH_mid;
      }
    }

    return (pH_low + pH_high) / 2.0;
  }

} // namespace OpenMS
