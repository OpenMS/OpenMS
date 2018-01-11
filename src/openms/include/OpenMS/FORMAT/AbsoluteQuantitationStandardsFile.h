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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_ABSOLUTEQUANTITATIONSTANDARDSFILE_H
#define OPENMS_FORMAT_ABSOLUTEQUANTITATIONSTANDARDSFILE_H

#include <OpenMS/config.h> // OPENMS_DLLAPI
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/METADATA/AbsoluteQuantitationStandards.h>

namespace OpenMS
{
  /**
    @brief Load files containing runConcentration data.
  */
  class OPENMS_DLLAPI AbsoluteQuantitationStandardsFile
  {
public:
    AbsoluteQuantitationStandardsFile() = default;
    virtual ~AbsoluteQuantitationStandardsFile() = default;

    /**
      @brief Load runConcentration data from a file and save it in memory.

      An example of the format expected:
      > sample_name,component_name,IS_component_name,actual_concentration,IS_actual_concentration,concentration_units,dilution_factor
      > 150516_CM1_Level1,23dpg.23dpg_1.Light,23dpg.23dpg_1.Heavy,0,1,uM,1
      > 150516_CM1_Level1,2mcit.2mcit_1.Light,2mcit.2mcit_1.Heavy,0,1,uM,1
      > 150516_CM1_Level1,2obut.2obut_1.Light,2obut.2obut_1.Heavy,0,1,uM,1

      @param[in] filename The file path from which the method is reading the data
      @param[out] run_concentrations Where the runConcentration data is going to be saved
    */
    void load(
      const String& filename,
      std::vector<AbsoluteQuantitationStandards::runConcentration>& run_concentrations
    ) const;

protected:
    /**
      @brief Extract one runConcentration from a single line

      Any missing information is going to be filled with default data:
      - an empty string for String data
      - the value 0.0 for concentration values
      - the value 1.0 for dilution factor

      The `headers` argument makes sure that the data is taken from the correct column/position in `line`.

      @param[in] line A list of strings each containing a column's info
      @param[in] headers A mapping from header name to position in the StringList given in input
    */
    AbsoluteQuantitationStandards::runConcentration extractRunFromLine_(
      const StringList& line,
      const std::map<String, Size>& headers
    ) const;
  };
}

#endif // OPENMS_FORMAT_ABSOLUTEQUANTITATIONSTANDARDSFILE_H
