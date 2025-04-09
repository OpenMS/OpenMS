// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ChromeleonFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <boost/regex.hpp>
#include <fstream>

namespace OpenMS
{
  void ChromeleonFile::load(const String& filename, MSExperiment& experiment) const
  {
    experiment.clear(true);
    std::ifstream ifs(filename, std::ifstream::in);
    if (!ifs.is_open())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    String line;
    MSChromatogram chromatogram;
    boost::smatch m;
    boost::regex re_channel("^Channel\t(.*)", boost::regex::no_mod_s);
    boost::regex re_injection("^Injection\t(.*)", boost::regex::no_mod_s);
    boost::regex re_processing_method("^Processing Method\t(.*)", boost::regex::no_mod_s);
    boost::regex re_instrument_method("^Instrument Method\t(.*)", boost::regex::no_mod_s);
    boost::regex re_injection_date("^Injection Date\t(.*)", boost::regex::no_mod_s);
    boost::regex re_injection_time("^Injection Time\t(.*)", boost::regex::no_mod_s);
    boost::regex re_detector("^Detector\t(.*)", boost::regex::no_mod_s);
    boost::regex re_signal_quantity("^Signal Quantity\t(.*)", boost::regex::no_mod_s);
    boost::regex re_signal_unit("^Signal Unit\t(.*)", boost::regex::no_mod_s);
    boost::regex re_signal_info("^Signal Info\t(.*)", boost::regex::no_mod_s);
    boost::regex re_raw_data("^Raw Data:", boost::regex::no_mod_s);
    boost::regex re_chromatogram_data("^Chromatogram Data:", boost::regex::no_mod_s);
    while (!ifs.eof())
    {
      TextFile::getLine(ifs, line);
      if (boost::regex_match(line, m, re_injection))
      {
        experiment.setMetaValue("mzml_id", m.str(1));
      }
      else if (boost::regex_match(line, m, re_channel))
      {
        experiment.setMetaValue("acq_method_name", m.str(1));
      }
      else if (boost::regex_match(line, m, re_processing_method))
      {
        experiment.getExperimentalSettings().getInstrument().getSoftware().setName(m.str(1));
      }
      else if (boost::regex_match(line, m, re_instrument_method))
      {
        experiment.getExperimentalSettings().getInstrument().setName(m.str(1));
      }
      else if (boost::regex_match(line, m, re_injection_date))
      {
        experiment.setMetaValue("injection_date", m.str(1));
      }
      else if (boost::regex_match(line, m, re_injection_time))
      {
        experiment.setMetaValue("injection_time", m.str(1));
      }
      else if (boost::regex_match(line, m, re_detector))
      {
        experiment.setMetaValue("detector", m.str(1));
      }
      else if (boost::regex_match(line, m, re_signal_quantity))
      {
        experiment.setMetaValue("signal_quantity", m.str(1));
      }
      else if (boost::regex_match(line, m, re_signal_unit))
      {
        experiment.setMetaValue("signal_unit", m.str(1));
      }
      else if (boost::regex_match(line, m, re_signal_info))
      {
        experiment.setMetaValue("signal_info", m.str(1));
      }
      else if (boost::regex_match(line, m, re_raw_data) ||
               boost::regex_match(line, m, re_chromatogram_data))
      {
        TextFile::getLine(ifs, line); // remove the subsequent line, right before the raw data
        break;
      }
    }
    while (!ifs.eof())
    {
      TextFile::getLine(ifs, line);
      std::vector<String> substrings;
      line.split('\t', substrings);
      if (substrings.size() == 3)
      {
        chromatogram.push_back(ChromatogramPeak(
          removeCommasAndParseDouble(substrings[0]),
          removeCommasAndParseDouble(substrings[2])));
      }
      else if (line.empty())
      {
        continue; // skips eventual empty lines, eg. the last before EOF
      }
      else
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, line, "Couldn't parse the raw data.");
      }
    }
    ifs.close();
    experiment.addChromatogram(chromatogram);
  }

  double ChromeleonFile::removeCommasAndParseDouble(String& number) const
  {
    return number.remove(',').toDouble();
  }
}
