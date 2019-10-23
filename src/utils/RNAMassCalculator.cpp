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

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/FORMAT/SVOutStream.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <iomanip> // for "setprecision"
#include <iostream>
#include <ostream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_RNAMassCalculator RNAMassCalculator

    @brief Calculates masses and mass-to-charge ratios of RNA sequences.

    Given an RNA (oligonucleotide) sequence and a charge state, the charged mass (including H+ adducts or losses) and the mass-to-charge ratio are computed.
    The sequence can include modifications (for information on valid notation see the @ref OpenMS::NASequence "NASequence" class documentation).
    Neutral masses can be computed by using "0" as charge state.

    Input can be given directly as values of the parameters: @p in_seq for RNA sequences and @p charge for charge states.
    Alternatively, it can be read from a file (see parameter @p in) with the following format: An RNA sequence at the beginning of each line, optionally followed by any number of charge states.
    Whitespace, commas or semicolons can de used to delimit the different items.
    Parts of the input that cannot be understood will be skipped.
    If charge states are given in the input file as well as via the @p charge parameter, results are returned for the union of both sets of charge states.

    Output can be written to a file or to the screen (see parameter @p out).
    Results for different charge states are always ordered from lowest to highest charge.
    A number of different output formats are available via the parameter @p format:
    - @p list writes a human-readable list of the form "ABCDEF: z=1 m=566.192 m/z=566.192, z=2 m=567.199 m/z=283.599";
    - @p table produces a CSV-like table (using parameter @p separator to delimit fields) with the columns "sequence", "charge", "mass", and "mass-to-charge", and with one row per sequence and charge state;
    - @p mass_only writes only mass values (one line per sequence, values for different charge states separated by spaces);
    - @p mz_only writes only mass-to-charge ratios (one line per sequence, values for different charge states separated by spaces).


    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_RNAMassCalculator.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_RNAMassCalculator.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPRNAMassCalculator :
  public TOPPBase
{
public:

  TOPPRNAMassCalculator() :
    TOPPBase("RNAMassCalculator", "Calculates masses and mass-to-charge ratios of RNA sequences", false), use_avg_mass_(false), output_(nullptr), format_(), frag_type_(NASequence::NASFragmentType::Full)
  {
    frag_type_names_["full"] = NASequence::NASFragmentType::Full;
    frag_type_names_["internal"] = NASequence::NASFragmentType::Internal;
    frag_type_names_["5-prime"] = NASequence::NASFragmentType::FivePrime;
    frag_type_names_["3-prime"] = NASequence::NASFragmentType::ThreePrime;
    frag_type_names_["a-B-ion"] = NASequence::NASFragmentType::AminusB;
    frag_type_names_["a-ion"] = NASequence::NASFragmentType::AIon;
    frag_type_names_["b-ion"] = NASequence::NASFragmentType::BIon;
    frag_type_names_["c-ion"] = NASequence::NASFragmentType::CIon;
    frag_type_names_["d-ion"] = NASequence::NASFragmentType::DIon;
    frag_type_names_["w-ion"] = NASequence::NASFragmentType::WIon;
    frag_type_names_["x-ion"] = NASequence::NASFragmentType::XIon;
    frag_type_names_["y-ion"] = NASequence::NASFragmentType::YIon;
    frag_type_names_["z-ion"] = NASequence::NASFragmentType::ZIon;
  }

protected:

  bool use_avg_mass_;
  ostream* output_;  // pointer to output stream (stdout or file)
  String format_, separator_;
  NASequence::NASFragmentType frag_type_;
  map<String, NASequence::NASFragmentType> frag_type_names_;

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file with RNA sequences and optionally charge numbers (mutually exclusive to 'in_seq')", false);
    setValidFormats_("in", vector<String>(1, "txt"));

    registerStringList_("in_seq", "<sequences>", StringList(), "List of RNA sequences (mutually exclusive to 'in')", false, false);

    registerOutputFile_("out", "<file>", "", "Output file; if empty, output is written to the screen", false);
    setValidFormats_("out", vector<String>(1, "txt"));

    registerIntList_("charge", "<numbers>", ListUtils::create<Int>("0"), "List of charge states; required if 'in_seq' is given", false);
    registerStringOption_("format", "<choice>", "list", "Output format ('list': human-readable list, 'table': CSV-like table, 'mass_only': mass values only, 'mz_only': m/z values only)\n", false);
    setValidStrings_("format", ListUtils::create<String>("list,table,mass_only,mz_only"));
    registerFlag_("average_mass", "Compute average (instead of monoisotopic) oligonucleotide masses");
    registerStringOption_("fragment_type", "<choice>", "full", "For what type of sequence/fragment the mass should be computed\n", false);
    setValidStrings_("fragment_type", ListUtils::create<String>("full,internal,5-prime,3-prime,a-B-ion,a-ion,b-ion,c-ion,d-ion,w-ion,x-ion,y-ion,z-ion"));
    registerStringOption_("separator", "<sep>", "", "Field separator for 'table' output format; by default, the 'tab' character is used", false);
  }

  double computeMass_(const NASequence& seq, Int charge) const
  {
    if (use_avg_mass_) return seq.getAverageWeight(frag_type_, charge);
    else return seq.getMonoWeight(frag_type_, charge);
  }

  void writeTable_(const NASequence& seq, const set<Int>& charges)
  {
    SVOutStream sv_out(*output_, separator_);
    for (set<Int>::const_iterator it = charges.begin(); it != charges.end();
         ++it)
    {
      double mass = computeMass_(seq, *it);
      sv_out << seq.toString() << *it << mass;
      sv_out.writeValueOrNan(abs(mass / *it));
      sv_out << endl;
    }
  }

  void writeList_(const NASequence& seq, const set<Int>& charges)
  {
    *output_ << seq.toString() << ": ";
    for (set<Int>::const_iterator it = charges.begin(); it != charges.end();
         ++it)
    {
      double mass = computeMass_(seq, *it);
      if (it != charges.begin()) *output_ << ", ";
      *output_ << "z=" << *it << " m=" << mass << " m/z=";
      if (*it != 0) *output_ << abs(mass / *it);
      else *output_ << "inf";
    }
    *output_ << endl;
  }

  void writeMassOnly_(const NASequence& seq, const set<Int>& charges,
                      bool mz = false)
  {
    for (set<Int>::const_iterator it = charges.begin(); it != charges.end();
         ++it)
    {
      double mass = computeMass_(seq, *it);
      if (it != charges.begin()) *output_ << " ";
      if (!mz) *output_ << mass;
      else if (*it == 0) *output_ << "inf";
      else *output_ << abs(mass / *it);
    }
    *output_ << endl;
  }

  void writeLine_(const NASequence& seq, const set<Int>& charges)
  {
    if (format_ == "list") writeList_(seq, charges);
    else if (format_ == "table") writeTable_(seq, charges);
    else if (format_ == "mass_only") writeMassOnly_(seq, charges);
    else writeMassOnly_(seq, charges, true); // "mz_only"
  }

  String getItem_(String& line, const String& skip = " \t,;")
  {
    Size pos = line.find_first_of(skip);
    String prefix = line.substr(0, pos);
    pos = line.find_first_not_of(skip, pos);
    if (pos == String::npos) line = "";
    else line = line.substr(pos);
    return prefix.trim();
  }

  void readFile_(const String& filename, const set<Int>& charges)
  {
    ifstream input(filename.c_str());
    String line;
    Size line_count(0);
    while (getline(input, line))
    {
      ++line_count;
      String item = getItem_(line);
      if ((item[0] == '"') && (item[item.size() - 1] == '"'))
      {
        item.unquote();
      }

      NASequence seq;
      try
      {
        seq = NASequence::fromString(item);
      }
      catch (Exception::ParseError& /*e*/)
      {
        OPENMS_LOG_WARN << "Warning: '" << item << "' is not a valid RNA sequence - skipping\n";
        continue;
      }

      set<Int> local_charges(charges);
      Size conversion_failed_count(0);
      while (!line.empty())
      {
        item = getItem_(line);
        try
        {
          local_charges.insert(item.toInt());
        }
        catch (Exception::ConversionError& /*e*/)
        {
          ++conversion_failed_count;
        }
      }
      if (conversion_failed_count)
      {
        OPENMS_LOG_WARN << "Warning: Invalid charge state specified in line:" << line_count << ".\n";
      }
      if (local_charges.empty())
      {
        OPENMS_LOG_WARN << "Warning: No charge state specified - skipping (line:" << line_count << ")\n";
        continue;
      }
      writeLine_(seq, local_charges);
    }
    input.close();
  }


  ExitCodes main_(int, const char**) override
  {
    String in = getStringOption_("in");
    StringList in_seq = getStringList_("in_seq");
    String out = getStringOption_("out");
    IntList charge_list = getIntList_("charge");
    set<Int> charges(charge_list.begin(), charge_list.end());
    use_avg_mass_ = getFlag_("average_mass");
    frag_type_ = frag_type_names_[getStringOption_("fragment_type")];

    ofstream outfile;
    if (out.empty())
    {
      output_ = &cout;
    }
    else
    {
      outfile.open(out.c_str());
      output_ = &outfile;
    }
    // use 4 decimal places:
    *output_ << std::fixed << std::setprecision(4);

    format_ = getStringOption_("format");
    if (format_ == "table")
    {
      separator_ = getStringOption_("separator");
      if (separator_.empty()) separator_ = "\t";
      // write header:
      SVOutStream sv_out(*output_, separator_);
      sv_out << "sequence" << "charge" << "mass" << "mass-to-charge" << endl;
    }

    if ((in.size() > 0) && (in_seq.size() > 0))
    {
      OPENMS_LOG_ERROR << "Specifying an input file and input sequences at the same time is not allowed!";
      return ILLEGAL_PARAMETERS;
    }

    if (in.size() > 0)
    {
      readFile_(in, charges);
    }
    else
    {
      if (charges.empty())
      {
        OPENMS_LOG_ERROR << "Error: No charge state specified";
        return ILLEGAL_PARAMETERS;
      }
      for (StringList::iterator it = in_seq.begin(); it != in_seq.end(); ++it)
      {
        NASequence seq;
        try
        {
         seq = NASequence::fromString(*it);
        }
        catch (Exception::ParseError& /*e*/)
        {
          OPENMS_LOG_WARN << "Warning: '" << *it << "' is not a valid RNA sequence - skipping\n";
          continue;
        }

        writeLine_(seq, charges);
      }
    }

    if (!out.empty()) outfile.close();

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPRNAMassCalculator tool;
  return tool.main(argc, argv);
}

/// @endcond
