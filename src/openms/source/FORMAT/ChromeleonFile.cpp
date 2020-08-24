// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
