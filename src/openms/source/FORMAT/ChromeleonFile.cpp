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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ChromeleonFile.h>

namespace OpenMS
{
  void ChromeleonFile::load(const String& filename, MSExperiment& experiment) const
  {
    std::ifstream ifs(filename, std::ifstream::in);
    if (!ifs.is_open())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    const Size BUFSIZE = 65536;
    char line[BUFSIZE];
    experiment.clear(true);
    MSChromatogram chromatogram;
    boost::cmatch m;
    boost::regex re_injection("^Injection\t(.+)", boost::regex::no_mod_s);
    boost::regex re_processing_method("^Processing Method\t(.+)", boost::regex::no_mod_s);
    boost::regex re_instrument_method("^Instrument Method\t(.+)", boost::regex::no_mod_s);
    boost::regex re_injection_date("^Injection Date\t(.+)", boost::regex::no_mod_s);
    boost::regex re_injection_time("^Injection Time\t(.+)", boost::regex::no_mod_s);
    boost::regex re_detector("^Detector\t(.+)", boost::regex::no_mod_s);
    boost::regex re_signal_quantity("^Signal Quantity\t(.+)", boost::regex::no_mod_s);
    boost::regex re_signal_unit("^Signal Unit\t(.+)", boost::regex::no_mod_s);
    boost::regex re_signal_info("^Signal Info\t(.+)", boost::regex::no_mod_s);
    boost::regex re_raw_data("^Raw Data:", boost::regex::no_mod_s);
    while (!ifs.eof())
    {
      ifs.getline(line, BUFSIZE);
      if (boost::regex_search(line, m, re_injection))
      {
        experiment.setMetaValue("mzml_id", std::string(m[1]));
      }
      else if (boost::regex_search(line, m, re_processing_method))
      {
        experiment.getExperimentalSettings().getInstrument().getSoftware().setName(String(m[1]));
      }
      else if (boost::regex_search(line, m, re_instrument_method))
      {
        experiment.getExperimentalSettings().getInstrument().setName(String(m[1]));
      }
      else if (boost::regex_search(line, m, re_injection_date))
      {
        experiment.setMetaValue("injection_date", std::string(m[1]));
      }
      else if (boost::regex_search(line, m, re_injection_time))
      {
        experiment.setMetaValue("injection_time", std::string(m[1]));
      }
      else if (boost::regex_search(line, m, re_detector))
      {
        experiment.setMetaValue("detector", std::string(m[1]));
      }
      else if (boost::regex_search(line, m, re_signal_quantity))
      {
        experiment.setMetaValue("signal_quantity", std::string(m[1]));
      }
      else if (boost::regex_search(line, m, re_signal_unit))
      {
        experiment.setMetaValue("signal_unit", std::string(m[1]));
      }
      else if (boost::regex_search(line, m, re_signal_info))
      {
        experiment.setMetaValue("signal_info", std::string(m[1]));
      }
      else if (boost::regex_search(line, m, re_raw_data))
      {
        ifs.getline(line, BUFSIZE); // remove the subsequent line, right before the raw data
        break;
      }
    }
    while (!ifs.eof())
    {
      ifs.getline(line, BUFSIZE);
      double rt, intensity;
      int ret = std::sscanf(line, "%lf\t%*f\t%lf", &rt, &intensity);
      if (ret == 2)
      {
        chromatogram.push_back(ChromatogramPeak(rt, intensity));
      }
      else if (!strcmp(line, "\r") || !strcmp(line, ""))
      {
        continue; // skips eventual empty lines, eg. the last before EOF
      }
      else
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, std::string(line), "Couldn't parse the raw data.");
      }
    }
    ifs.close();
    experiment.addChromatogram(chromatogram);
  }
}
