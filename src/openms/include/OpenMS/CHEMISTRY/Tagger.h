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
// $Authors: Timo Sachsenberg, Eugen Netz $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <set>
#include <vector>
#include <string>
#include <map>

namespace OpenMS
{
  /**
    *  @brief Calculate sequence tags from m/z values
    *
    *  Current restrictions and potential extensions:
    *  - first prefix/suffix ion don't give rise to first character in tags (as currently ion types are not specified)
    *  - gaps are not supported
    **/
  class OPENMS_DLLAPI Tagger
  {
    public:
      // initalize tagger with minimum/maximum tag length, +- tolerance ppm, min/max charge, fixed and variable modifications
      Tagger(size_t min_tag_length, double ppm, size_t max_tag_length = 65535, size_t min_charge = 1, size_t max_charge = 1, const StringList& fixed_mods = StringList(), const StringList& var_mods = StringList());

      // generate tags from mass vector @p mzs using the standard residues in ResidueDB
      void getTag(const std::vector<double>& mzs, std::vector<std::string>& tags) const;
      void getTag(const MSSpectrum& spec, std::vector<std::string>& tags) const;
      void setMaxCharge(size_t max_charge);

    private:
      double min_gap_; // will be set to smallest residue mass in ResidueDB
      double max_gap_; // will be set to highest residue mass in ResidueDB
      double ppm_; // < tolerance
      size_t min_tag_length_; // < minimum tag length
      size_t max_tag_length_; // < maximum tag length
      size_t min_charge_; // < minimal fragment charge
      size_t max_charge_; // < maximal fragment charge
      std::map<double, char> mass2aa;
      char getAAByMass_(double m) const;
      void getTag_(std::string& tag, const std::vector<double>& mzs, const size_t i, std::vector<std::string>& tags) const;
  };
}
