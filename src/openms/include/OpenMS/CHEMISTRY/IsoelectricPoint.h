// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: OpenMS Team $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{
  /**
  @ingroup Chemistry

  @brief Utility class for computing the isoelectric point (pI) and net charge of peptides.

  This class provides methods to compute:
  - The net charge of an amino acid sequence at a given pH using the Henderson-Hasselbalch equation
  - The isoelectric point (pI), i.e. the pH at which the net charge is zero, found via bisection

  The calculation considers:
  - The N-terminal amino group (pKa), unless the peptide has an N-terminal modification
  - The C-terminal carboxyl group (pKa), unless the peptide has a C-terminal modification
  - Ionizable side chains of D, E, C, U, Y, H, K, R

  The default pK values follow the Lehninger scale (Nelson & Cox, Lehninger Principles of Biochemistry).
  Other supported scales (EMBOSS, Sillero) can be selected via the ProteomicsPkaScale enum.

  Side-chain PTM-specific pKa shifts are currently not modeled; modified residues are evaluated using the
  parent residue one-letter code. Pyrrolysine (O) is not supported and raises Exception::InvalidValue.
  If the zero crossing lies outside the searched pH interval [0, 14], computePI() returns the closer boundary
  and logs a warning.

  References:
  - Nelson DL, Cox MM. Lehninger Principles of Biochemistry. 6th ed. (2013).
  - Sillero A, Ribeiro JM. Isoelectric points of proteins... Anal Biochem. 1989;179:319-325.

  @todo Timo Sachsenberg Add the Bjellqvist scale with residue-dependent N-terminal pKa values.
  */

  /// @brief Enum for pKa scales used in isoelectric point calculation
  enum class ProteomicsPkaScale
  {
    LEHNINGER = 0, ///< Lehninger (Nelson & Cox) scale
    EMBOSS = 1,    ///< EMBOSS scale (used by pepstats)
    SILLERO = 2,   ///< Sillero & Ribeiro scale
    SIZE_OF_PROTEOMICS_PKA_SCALES
  };

  class OPENMS_DLLAPI IsoelectricPoint
  {

  public:

    /// @brief Computes the net charge of an amino acid sequence at a given pH
    /// @param seq The amino acid sequence
    /// @param pH The pH value at which to compute the charge
    /// @param scale The pKa scale to use (default: Lehninger)
    /// @return The net charge of the peptide at the given pH
    /// @throw Exception::InvalidValue if the sequence is empty
    static double computeCharge(
      const AASequence& seq,
      double pH,
      ProteomicsPkaScale scale = ProteomicsPkaScale::LEHNINGER);

    /// @brief Computes the isoelectric point (pI) of an amino acid sequence via bisection
    /// @param seq The amino acid sequence
    /// @param scale The pKa scale to use (default: Lehninger)
    /// @param tolerance The convergence tolerance for the bisection (default: 1e-4)
    /// @return The isoelectric point (pH where net charge ≈ 0)
    /// @throw Exception::InvalidValue if the sequence is empty
    static double computePI(
      const AASequence& seq,
      ProteomicsPkaScale scale = ProteomicsPkaScale::LEHNINGER,
      double tolerance = 1e-4);

  private:
    /// Internal struct holding pKa values for a given scale
    struct PkaValues
    {
      double nterm;  ///< pKa of the N-terminal amino group
      double cterm;  ///< pKa of the C-terminal carboxyl group
      double D;      ///< Aspartate side chain
      double E;      ///< Glutamate side chain
      double C;      ///< Cysteine side chain
      double Y;      ///< Tyrosine side chain
      double H;      ///< Histidine side chain
      double K;      ///< Lysine side chain
      double R;      ///< Arginine side chain
    };

    /// Returns the pKa values for the specified scale
    static const PkaValues& getPkaValues_(ProteomicsPkaScale scale);

    /// Charge contribution from an acidic group (COOH, side chains of D, E, C, Y)
    static double chargeAcidic_(double pH, double pKa);

    /// Charge contribution from a basic group (NH2, side chains of H, K, R)
    static double chargeBasic_(double pH, double pKa);
  };

} // namespace OpenMS
