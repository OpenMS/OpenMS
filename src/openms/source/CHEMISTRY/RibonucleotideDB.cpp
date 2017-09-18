// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>

#include <QFile>
#include <QTextStream>

using namespace std;

namespace OpenMS
{
  RibonucleotideDB::RibonucleotideDB():
    max_code_length_(0)
  {
    readFromFile_("CHEMISTRY/Modomics.tsv");
  }


  RibonucleotideDB::~RibonucleotideDB()
  {
  }


  void RibonucleotideDB::readFromFile_(const String& path)
  {
    String full_path = File::find(path);

    const String header = "name\tshort_name\tnew_nomenclature\toriginating_base\trnamods_abbrev\thtml_abbrev\tformula\tmonoisotopic_mass\taverage_mass";

    // the input file is Unicode encoded, so we need Qt to read it:
    QFile file(full_path.toQString());
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, full_path);
    }

    QTextStream source(&file);
    source.setCodec("UTF-8");
    String line = source.readLine();
    if (line != header)
    {
      String msg = "expected header line starting with: '" + header + "'";
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  line, msg);
    }
    Size line_count = 1;

    QChar prime(0x2032); // Unicode "prime" character
    while (!source.atEnd())
    {
      line_count++;
      QString row = source.readLine();
      // replace all "prime" characters with apostrophes (e.g. in "5'", "3'"):
      row.replace(prime, '\'');
      try
      {
        Ribonucleotide ribo = parseRow_(row, line_count);
        code_map_[ribo.code] = ribonucleotides_.size();
        ribonucleotides_.push_back(ribo);
        max_code_length_ = max(max_code_length_, ribo.code.size());
      }
      catch (...)
      {
        LOG_ERROR << "Error: Failed to parse input line " << line_count
                  << " - skipping this line." << endl;
      }
    }
  }


  Ribonucleotide RibonucleotideDB::parseRow_(const String& row,
                                             Size line_count)
  {
    vector<String> parts;
    row.split('\t', parts);
    if (parts.size() != 9)
    {
      String msg = "9 tab-separated fields expected, found " +
        String(parts.size()) + " in line " + String(line_count);
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  row, msg);
    }
    Ribonucleotide ribo;
    ribo.name = parts[0];
    ribo.code = parts[1];
    ribo.new_code = parts[2];
    if (parts[3] == "preQ0base")
    {
      ribo.origin = '0';
    }
    else if (parts[3].size() == 1) // A, C, G, U
    {
      ribo.origin = parts[3][0];
    }
    // "parts[4]" is the Unicode equivalent to "parts[5]", so we can skip it
    ribo.html_code = parts[5];
    ribo.formula = parts[6];
    if (!parts[7].empty() && (parts[7] != "None"))
    {
      ribo.mono_mass = parts[7].toDouble();
    }
    if (!parts[8].empty() && (parts[8] != "None"))
    {
      ribo.avg_mass = parts[8].toDouble();
    }
    return ribo;
  }


  const Ribonucleotide& RibonucleotideDB::getRibonucleotide(const String& code)
  {
    boost::unordered_map<String, Size>::const_iterator pos =
      code_map_.find(code);
    if (pos == code_map_.end())
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, code);
    }
    return ribonucleotides_[pos->second];
  }


  const Ribonucleotide& RibonucleotideDB::getRibonucleotidePrefix(const String&
                                                                  seq)
  {
    String prefix = seq.substr(0, max_code_length_);
    while (!prefix.empty())
    {
      boost::unordered_map<String, Size>::const_iterator pos =
        code_map_.find(prefix);
      if (pos != code_map_.end())
      {
        return ribonucleotides_[pos->second];
      }
      prefix = prefix.substr(0, prefix.size() - 1);
    }
    throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     seq);
  }
}
