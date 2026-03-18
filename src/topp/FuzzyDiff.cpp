// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <fstream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_FuzzyDiff FuzzyDiff

@brief Compares two files, tolerating numeric differences.

In the diff output, \"position\" refers to the characters in the string, whereas \"column\" is meant for the text editor.

Only one of 'ratio' or 'absdiff' has to be satisfied.  Use \"absdiff\" to deal with cases like \"zero vs. epsilon\".

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_FuzzyDiff.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_FuzzyDiff.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFuzzyDiff :
  public TOPPBase
{
public:
  TOPPFuzzyDiff() :
    TOPPBase("FuzzyDiff", "Compares two files, tolerating numeric differences.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    addEmptyLine_();
    registerInputFile_("in1", "<file>", "", "first input file", true, false);
    registerInputFile_("in2", "<file>", "", "second input file", true, false);
    addEmptyLine_();

    registerDoubleOption_("ratio", "<double>", 1, R"(acceptable relative error. Only one of 'ratio' or 'absdiff' has to be satisfied.  Use "absdiff" to deal with cases like "zero vs. epsilon".)", false, false);
    setMinFloat_("ratio", 1);
    registerDoubleOption_("absdiff", "<double>", 0, "acceptable absolute difference. Only one of 'ratio' or 'absdiff' has to be satisfied. ", false, false);
    setMinFloat_("absdiff", 0);
    addEmptyLine_();

    registerStringList_("whitelist", "<string list>", ListUtils::create<String>("<?xml-stylesheet"), "Lines containing one of these strings are skipped", false, true);

    registerStringList_("matched_whitelist", "<string list>", ListUtils::create<String>(""), "Lines where one file contains one string and the other file another string are skipped. Input is given as list of colon separated tuples, e.g. String1:String2 String3:String4", false, true);

    registerIntOption_("verbose", "<int>", 2, "set verbose level:\n"
                                              "0 = very quiet mode (absolutely no output)\n"
                                              "1 = quiet mode (no output unless differences detected)\n"
                                              "2 = default (include summary at end)\n"
                                              "3 = continue after errors\n",
                       false, false
                       );
    setMinInt_("verbose", 0);
    setMaxInt_("verbose", 3);
    registerIntOption_("tab_width", "<int>", 8, "tabulator width, used for calculation of column numbers", false, false);
    setMinInt_("tab_width", 1);
    registerIntOption_("first_column", "<int>", 1, "number of first column, used for calculation of column numbers", false, false);
    setMinInt_("first_column", 0);
    addEmptyLine_();
    registerFlag_("sort", "sort the input files before comparison (useful for tabular files where row order may vary). The first line of each file is assumed to be a header and is not sorted.");
  }

  ExitCodes main_(int, const char **) override
  {

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    String in1 = getStringOption_("in1");
    String in2 = getStringOption_("in2");
    double acceptable_ratio = getDoubleOption_("ratio");
    double acceptable_absdiff = getDoubleOption_("absdiff");
    StringList whitelist = getStringList_("whitelist");
    StringList raw_matched_whitelist = getStringList_("matched_whitelist");
    int verbose_level = getIntOption_("verbose");
    int tab_width = getIntOption_("tab_width");
    int first_column = getIntOption_("first_column");
    bool do_sort = getFlag_("sort");

    // This is for debugging the parsing of whitelist_ from cmdline or ini file.  Converting StringList back to String is intentional.
    writeDebug_(String("whitelist: ") + String(whitelist) + " (size: " + whitelist.size() + ")", 1);
    writeDebug_(String("matched_whitelist: ") + String(raw_matched_whitelist) + " (size: " + raw_matched_whitelist.size() + ")", 1);

    OpenMS::FuzzyStringComparator fsc;

    std::vector< std::pair<std::string, std::string> > parsed_matched_whitelist; 
    for (Size i = 0; i < raw_matched_whitelist.size(); i++)
    {

      // Split each entry at the colon to produce a pair of strings
      std::vector<String> tmp;
      raw_matched_whitelist[i].split(":", tmp);
      if (tmp.size() != 2)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          String(raw_matched_whitelist[i]) + " does not have the format String1:String2");
      }

      parsed_matched_whitelist.emplace_back(tmp[0], tmp[1]);
    }

    fsc.setAcceptableRelative(acceptable_ratio);
    fsc.setAcceptableAbsolute(acceptable_absdiff);
    fsc.setWhitelist(whitelist);
    fsc.setMatchedWhitelist(parsed_matched_whitelist);
    fsc.setVerboseLevel(verbose_level);
    fsc.setTabWidth(tab_width);
    fsc.setFirstColumn(first_column);

    // If sorting is requested, sort both files (keeping header) and write to temp files
    String compare_in1 = in1;
    String compare_in2 = in2;
    String temp_file1, temp_file2;

    if (do_sort)
    {
      auto sortFile = [](const String& input_file, String& output_file)
      {
        std::ifstream infile(input_file);
        if (!infile)
        {
          throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, input_file);
        }

        std::string header;
        std::vector<std::string> lines;

        // Read header (first line)
        if (std::getline(infile, header))
        {
          // Read remaining lines
          std::string line;
          while (std::getline(infile, line))
          {
            lines.push_back(line);
          }
        }
        infile.close();

        // Sort data lines
        std::sort(lines.begin(), lines.end());

        // Write to temp file
        output_file = File::getTempDirectory() + "/" + File::basename(input_file) + ".sorted." + File::getUniqueName() + ".tmp";
        std::ofstream outfile(output_file);
        if (!outfile)
        {
          throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, output_file);
        }
        outfile << header << "\n";
        for (const auto& l : lines)
        {
          outfile << l << "\n";
        }
        outfile.close();
      };

      sortFile(in1, temp_file1);
      sortFile(in2, temp_file2);
      compare_in1 = temp_file1;
      compare_in2 = temp_file2;
    }

    bool result = fsc.compareFiles(compare_in1, compare_in2);

    // Clean up temp files
    if (do_sort)
    {
      File::remove(temp_file1);
      File::remove(temp_file2);
    }

    if (result)
    {
      return EXECUTION_OK;
    }
    else
    {
      // TODO think about better exit codes.
      return PARSE_ERROR;
    }
  }

};

int main(int argc, const char ** argv)
{
  TOPPFuzzyDiff tool;
  return tool.main(argc, argv);
}

/// @endcond
