// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <string_view>

namespace OpenMS
{
  // add common enums here to avoid big includes of large classes and break circular dependencies

  /// Enum for different units which can be displayed on a plotting axis
  /// The order is arbitrary.
  enum class DIM_UNIT
  {
    RT = 0,   ///< RT in seconds
    MZ,       ///< m/z
    INT,      ///< intensity
    IM_MS,    ///< ion mobility milliseconds
    IM_VSSC,  ///< volt-second per square centimeter (i.e. 1/K_0)
    FAIMS_CV, ///< FAIMS compensation voltage
    SIZE_OF_DIM_UNITS
  };
  inline std::string_view DIM_NAMES[(int)DIM_UNIT::SIZE_OF_DIM_UNITS] = {"RT [s]", "m/z [Th]", "intensity", "IM [milliseconds]", "IM [vs / cm2]", "FAIMS CV"};
  inline std::string_view DIM_NAMES_SHORT[(int)DIM_UNIT::SIZE_OF_DIM_UNITS] = {"RT", "m/z", "int", "IM", "IM", "FCV"};

} // namespace OpenMS
