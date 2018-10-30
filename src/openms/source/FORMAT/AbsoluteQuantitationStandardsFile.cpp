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

#include <OpenMS/FORMAT/AbsoluteQuantitationStandardsFile.h>

namespace OpenMS
{
  void AbsoluteQuantitationStandardsFile::load(
    const String& filename,
    std::vector<AbsoluteQuantitationStandards::runConcentration>& run_concentrations
  ) const
  {
    CsvFile csv(filename);
    StringList sl;
    std::map<String, Size> headers;
    if (csv.rowCount() > 0) // avoid accessing a row in an empty file
    {
      csv.getRow(0, sl);
    }
    for (Size i = 0; i < sl.size(); ++i)
    {
      headers[sl[i]] = i; // for each header found, assign an index value to it
    }
    run_concentrations.clear();
    for (Size i = 1; i < csv.rowCount(); ++i)
    {
      csv.getRow(i, sl);
      run_concentrations.push_back(extractRunFromLine_(sl, headers));
    }
  }

  AbsoluteQuantitationStandards::runConcentration AbsoluteQuantitationStandardsFile::extractRunFromLine_(
    const StringList& line,
    const std::map<String, Size>& headers
  ) const
  {
    AbsoluteQuantitationStandards::runConcentration rc;
    std::map<String, Size>::const_iterator it;
    it = headers.find("sample_name");
    rc.sample_name = it != headers.end() ? line[it->second] : "";
    it = headers.find("component_name");
    rc.component_name = it != headers.end() ? line[it->second] : "";
    it = headers.find("IS_component_name");
    rc.IS_component_name = it != headers.end() ? line[it->second] : "";
    it = headers.find("actual_concentration");
    rc.actual_concentration = it != headers.end() ? line[it->second].toDouble() : 0.0;
    it = headers.find("IS_actual_concentration");
    rc.IS_actual_concentration = it != headers.end() ? line[it->second].toDouble() : 0.0;
    it = headers.find("concentration_units");
    rc.concentration_units = it != headers.end() ? line[it->second] : "";
    it = headers.find("dilution_factor");
    rc.dilution_factor = it != headers.end() ? line[it->second].toDouble() : 1.0;
    return rc;
  }
}
