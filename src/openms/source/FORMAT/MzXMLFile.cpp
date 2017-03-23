// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzXMLFile.h>

#include <fstream>

using namespace std;

namespace OpenMS
{
  MzXMLFile::MzXMLFile() :
    XMLFile("/SCHEMAS/mzXML_2_1.xsd", "2.1")
  {
  }

  MzXMLFile::~MzXMLFile()
  {
  }

  PeakFileOptions & MzXMLFile::getOptions()
  {
    return options_;
  }

  const PeakFileOptions & MzXMLFile::getOptions() const
  {
    return options_;
  }

  void MzXMLFile::setOptions(const PeakFileOptions & options)
  {
      options_ = options;
  }

  void MzXMLFile::load(const String & filename, MapType & map)
  {
    map.reset();

    //set DocumentIdentifier
    map.setLoadedFileType(filename);
    map.setLoadedFilePath(filename);

    Internal::MzXMLHandler handler(map, filename, schema_version_, *this);
    handler.setOptions(options_);
    parse_(filename, &handler);
  }

  void MzXMLFile::load(const String & filename, MSExperiment<RichPeak1D> & map)
  {
#ifdef OPENMS_ASSERTIONS
    std::cout << "===========================================================================" << std::endl;
    std::cout << "WARNING: you are using a deprecated interface MzXMLFile::load with MSExperiment<RichPeak1D>." 
              << "\nPlease consider switching to MSExperiment<Peak1D>" << std::endl;
    std::cout << "===========================================================================" << std::endl;
#endif
    MSExperiment<Peak1D> map_;

    //set DocumentIdentifier
    map.setLoadedFileType(filename);
    map.setLoadedFilePath(filename);

    Internal::MzXMLHandler handler(map_, filename, schema_version_, *this);
    handler.setOptions(options_);
    parse_(filename, &handler);

    map.reset();
    map = (ExperimentalSettings)map_;
    map.setChromatograms(map_.getChromatograms());

    // convert regular Spectra to RichPeakSpectra
    for (Size k = 0; k < map_.getNrSpectra(); k++)
    {
      // copy each spectrum over
      MSSpectrum<RichPeak1D> s;
      const MSSpectrum<Peak1D> s2 = map_.getSpectrum(k);
      s = (SpectrumSettings)map_.getSpectrum(k);
      s.setRT(s2.getRT());
      s.setDriftTime(s2.getDriftTime());
      s.setMSLevel(s2.getMSLevel());
      s.setName(s2.getName());
      s.setFloatDataArrays(s2.getFloatDataArrays());
      s.setStringDataArrays(s2.getStringDataArrays());
      s.setIntegerDataArrays(s2.getIntegerDataArrays());
      for (MSSpectrum<Peak1D>::iterator it = map_.getSpectrum(k).begin(); it != map_.getSpectrum(k).end(); ++it)
      {
        RichPeak1D p;
        p.setMZ(it->getMZ());
        p.setIntensity(it->getIntensity());
      }
      map.addSpectrum(s);
    }
  }

} // namespace OpenMS

