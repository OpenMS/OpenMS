// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Timo Sachsenberg, Eugen Netz $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <map>

namespace OpenMS
{
  class MSSpectrum;

  /**
    *  @brief Calculate sequence tags from m/z values
    *
    *  Current restrictions and potential extensions:
    *  - first prefix/suffix ion don't give rise to first character in tags (as currently ion types are not specified)
    *  - gaps are not supported
    *  - I/L are treated redundantly, all possible combinations are written out
    **/
  class OPENMS_DLLAPI Tagger
  {
    public:
      /**
            @brief Constructor for Tagger

            The parameter @p max_charge_ should be >= @p min_charge_.
            Also @p max_tag_length should be >= @p min_tag_length.

            @param[in] min_tag_length the minimal sequence tag length.
            @param[in] tolerance the tolerance for matching residue masses to peak delta masses.
            @param[in] max_tag_length the maximal sequence tag length.
            @param[in] min_charge minimal fragment charge considered for each sequence tag.
            @param[in] max_charge maximal fragment charge considered for each sequence tag.
            @param[in,out] fixed_mods a list of modification names. The modified residues replace the unmodified versions.
            @param[in,out] var_mods a list of modification names. The modified residues are added as additional entries to the list of residues.
            @param[in] tol_is_ppm if set to true, the tolerance is interpreted as ppm, otherwise as absolute value. Defaults to true.
          */
      Tagger(size_t min_tag_length, double tolerance, size_t max_tag_length = 65535, size_t min_charge = 1, size_t max_charge = 1, const StringList& fixed_mods = StringList(), const StringList& var_mods = StringList(), bool tol_is_ppm = true);

      /**
            @brief Generate tags from mass vector @p mzs

            The parameter @p tags is filled with one string per sequence tag.
            It uses the standard residues from ResidueDB including
            the fixed and variable modifications given to the constructor.

            @param[in] mzs a vector of mz values, containing the mz values from a centroided fragment spectrum.
            @param[out] tags the vector of tags, that is filled with this function.
          */
      void getTag(const std::vector<double>& mzs, std::vector<std::string>& tags) const;

      /**
            @brief Generate tags from an MSSpectrum

            The parameter @p tags is filled with one string per sequence tag.
            It uses the standard residues from ResidueDB including
            the fixed and variable modifications given to the constructor.

            @param[in] spec a centroided fragment spectrum.
            @param[out] tags the vector of tags, that is filled with this function.
          */
      void getTag(const MSSpectrum& spec, std::vector<std::string>& tags) const;

      /**
            @brief Change the maximal charge considered by the tagger

            Allows to change the maximal considered charge e.g. based on a spectra
            precursor charge without calling the constructor multiple times.

            @param[in] max_charge the new maximal charge.
          */
      void setMaxCharge(size_t max_charge);


    private:
      double min_gap_; ///< will be set to smallest residue mass in ResidueDB
      double max_gap_; ///< will be set to highest residue mass in ResidueDB
      double tolerance_; ///< tolerance
      bool tol_is_ppm_; ///< is the tolerance in ppm
      size_t min_tag_length_; ///< minimum tag length
      size_t max_tag_length_; ///< maximum tag length
      size_t min_charge_; ///< minimal fragment charge
      size_t max_charge_; ///< maximal fragment charge
      std::map<double, char> mass2aa_; ///< mapping of residue masses to their one letter codes

      /// get a residue one letter code by matching the mass @p m to the map of residues mass2aa_, returns ' ' if there is no match.
      char getAAByMass_(double m) const;

      /// start searching for tags starting from peak @p i of the mz vector @p mzs
      void getTag_(std::string& tag, const std::vector<double>& mzs, const size_t i, std::vector<std::string>& tags, const size_t charge) const;
  };
}

