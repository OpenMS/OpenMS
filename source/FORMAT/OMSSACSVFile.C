// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/OMSSACSVFile.h>

#include <fstream>
#include <algorithm>

using namespace std;

namespace OpenMS
{

  OMSSACSVFile::OMSSACSVFile()
  {
  }

  OMSSACSVFile::~OMSSACSVFile()
  {
  }

  void OMSSACSVFile::load(const String & filename, ProteinIdentification & /* protein_identification */, vector<PeptideIdentification> & id_data) const
  {
    ifstream is(filename.c_str());
    if (!is)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }

    String line;
    getline(is, line, '\n');
    if (!line.hasPrefix("Spectrum"))
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "first line does not contain a description", "");
    }

    // line number counter
    Size line_number = 0;

    // ignore first line
    Int actual_spectrum_number(-1);
    while (getline(is, line, '\n'))
    {
      ++line_number;

      // Spectrum number, Filename/id, Peptide, E-value, Mass, gi, Accession, Start, Stop, Defline, Mods, Charge, Theo Mass, P-value
      // 1,,MSHHWGYGK,0.00336754,1101.49,0,6599,1,9,CRHU2 carbonate dehydratase (EC 4.2.1.1) II [validated] - human,,1,1101.48,1.30819e-08
      line.trim();

      // replace ',' in protein name
      String::ConstIterator it = find(line.begin(), line.end(), '"');
      UInt offset(0);
      if (it != line.end())
      {
        while (*(++it) != '"')
        {
          if (*it == ',')
          {
            offset++;
          }
        }
      }
      vector<String> split;
      line.split(',', split);
      if (split.size() != 14 && split.size() != 14 + offset)
      {
        throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, line, "number of columns should be 14 in line " + String(line_number));
      }
      PeptideHit p;
      p.setSequence(AASequence(split[2].trim()));
      p.setScore(split[13 + offset].trim().toDouble());
      p.setCharge(split[11 + offset].trim().toInt());

      if (actual_spectrum_number != split[0].trim().toInt())
      {
        // new id
        //id_data.push_back(IdentificationData());
        id_data.push_back(PeptideIdentification());
        id_data.back().setScoreType("OMSSA");
        actual_spectrum_number = (UInt)split[0].trim().toInt();
      }

      id_data.back().insertHit(p);
    }

  }

} // namespace OpenMS
