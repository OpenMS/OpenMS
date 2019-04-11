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
    readFromFile_("CHEMISTRY/Custom_RNA_modifications.tsv");
  }


  RibonucleotideDB::~RibonucleotideDB()
  {
    for (auto& r : ribonucleotides_)
    {
      delete(r);
    }
  }


  void RibonucleotideDB::readFromFile_(const std::string& path)
  {
    String full_path = File::find(path);

    String header = "name\tshort_name\tnew_nomenclature\toriginating_base\trnamods_abbrev\thtml_abbrev\tformula\tmonoisotopic_mass\taverage_mass";

    // the input file is Unicode encoded, so we need Qt to read it:
    QFile file(full_path.toQString());
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, full_path);
    }

    QTextStream source(&file);
    source.setCodec("UTF-8");
    Size line_count = 1;
    String line = source.readLine();
    while (line[0] == '#') // skip leading comments
    {
      line = source.readLine();
      ++line_count;
    }
    if (!line.hasPrefix(header)) // additional columns are allowed
    {
      String msg = "expected header line starting with: '" + header + "'";
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  line, msg);
    }

    QChar prime(0x2032); // Unicode "prime" character
    while (!source.atEnd())
    {
      line_count++;
      QString row = source.readLine();

      // replace all "prime" characters with apostrophes (e.g. in "5'", "3'"):
      row.replace(prime, '\'');
      try
      {
        ConstRibonucleotidePtr ribo = parseRow_(row.toStdString(), line_count);
        code_map_[ribo->getCode()] = ribonucleotides_.size();
        ribonucleotides_.push_back(ribo);
        max_code_length_ = max(max_code_length_, ribo->getCode().size());
      }
      catch (Exception::BaseException& e)
      {
        LOG_ERROR << "Error: Failed to parse input line " << line_count
                  << ". Reason:\n" << e.getName()
                  << " - " << e.getMessage() << "\nSkipping this line." << endl;
      }
    }
  }


  RibonucleotideDB::ConstRibonucleotidePtr RibonucleotideDB::parseRow_(
    const std::string& row, Size line_count)
  {
    vector<String> parts;
    String(row).split('\t', parts);
    if (parts.size() < 9)
    {
      String msg = "9 tab-separated fields expected, found " +
        String(parts.size()) + " in line " + String(line_count);
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  row, msg);
    }
    Ribonucleotide* ribo = new Ribonucleotide();
    ribo->setName(parts[0]);
    if (parts[1].hasSuffix("QtRNA")) // use just "Q" instead of "QtRNA"
    {
      ribo->setCode(parts[1].chop(4));
    }
    else
    {
      ribo->setCode(parts[1]);
    }
    ribo->setNewCode(parts[2]);
    if (parts[3] == "preQ0base") // queuosine and its derivatives
    {
      ribo->setOrigin('G'); // queuosine replaces "G" in tRNA-Asp/Asn
    }
    else if (parts[3].size() == 1) // A, C, G, U
    {
      ribo->setOrigin(parts[3][0]);
    }
    // "parts[4]" is the Unicode equivalent to "parts[5]", so we can skip it
    ribo->setHTMLCode(parts[5]);
    if (!parts[6].empty() && (parts[6] != "-"))
    {
      ribo->setFormula(EmpiricalFormula(parts[6]));
    }
    if (!parts[7].empty() && (parts[7] != "None"))
    {
      ribo->setMonoMass(parts[7].toDouble());
      if ((ribo->getMonoMass() == 0.0) && (!ribo->getFormula().isEmpty()))
      {
        ribo->setMonoMass(ribo->getFormula().getMonoWeight());
      }
    }
    if (!parts[8].empty() && (parts[8] != "None"))
    {
      ribo->setAvgMass(parts[8].toDouble());
      if ((ribo->getAvgMass() == 0.0) && (!ribo->getFormula().isEmpty()))
      {
        ribo->setAvgMass(ribo->getFormula().getAverageWeight());
      }
    }
    // Modomics' "new code" contains information on terminal specificity:
    if (parts[2].back() == 'N') // terminal mod., exception: "GN"
    {
      if (parts[2].hasSubstring("55") || (parts[2] == "N"))
      {
        ribo->setTermSpecificity(Ribonucleotide::FIVE_PRIME);
      }
      else if (parts[2].hasSubstring("33"))
      {
        ribo->setTermSpecificity(Ribonucleotide::THREE_PRIME);
      }
    }
    else // default specificity is "ANYWHERE"; now set formula after base loss:
    {
      if (parts[1].back() == 'm') // mod. attached to the ribose, not base
      {
        ribo->setBaselossFormula(EmpiricalFormula("C6H12O5"));
      }
      else if (parts[1].back() == '?') // ambiguity code -> fill the map
      {
        if (parts.size() < 10)
        {
          String msg =
            "10th field expected for ambiguous modification in line " +
            String(line_count);
          throw Exception::ParseError(__FILE__, __LINE__,
                                      OPENMS_PRETTY_FUNCTION, row, msg);
        }
        String code1 = parts[9].prefix(' '), code2 = parts[9].suffix(' ');
        ambiguity_map_[parts[1]] = make_pair(getRibonucleotide(code1),
                                             getRibonucleotide(code2));
      }
      else if ((parts[1] == "Ar(p)") || (parts[1] == "Gr(p)"))
      {
        ribo->setBaselossFormula(EmpiricalFormula("C10H19O21P"));
      }
    }

    return ribo;
  }


  RibonucleotideDB::ConstRibonucleotidePtr
  RibonucleotideDB::getRibonucleotide(const std::string& code)
  {
    std::unordered_map<std::string, Size>::const_iterator pos =
      code_map_.find(code);
    if (pos == code_map_.end())
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, code);
    }
    return ribonucleotides_[pos->second];
  }


  RibonucleotideDB::ConstRibonucleotidePtr
  RibonucleotideDB::getRibonucleotidePrefix(const std::string& seq)
  {
    std::string prefix = seq.substr(0, max_code_length_);
    while (!prefix.empty())
    {
      auto pos = code_map_.find(prefix);
      if (pos != code_map_.end())
      {
        return ribonucleotides_[pos->second];
      }
      prefix = prefix.substr(0, prefix.size() - 1);
    }
    throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     seq);
  }


  pair<RibonucleotideDB::ConstRibonucleotidePtr,
       RibonucleotideDB::ConstRibonucleotidePtr>
  RibonucleotideDB::getRibonucleotideAlternatives(const std::string& code)
  {
    auto pos = ambiguity_map_.find(code);
    if (pos == ambiguity_map_.end())
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, code);
    }
    return pos->second;
  }
}

