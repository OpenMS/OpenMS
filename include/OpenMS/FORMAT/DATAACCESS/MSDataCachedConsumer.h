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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_DATAACCESS_MSDATACACHEDCONSUMER_H
#define OPENMS_FORMAT_DATAACCESS_MSDATACACHEDCONSUMER_H

#include <OpenMS/INTERFACES/IMSDataConsumer.h>

#include <OpenMS/ANALYSIS/OPENSWATH/CachedmzML.h>

namespace OpenMS
{
    /**
      @brief Transforming and cached writing consumer of MS data

      Is able to transform a spectrum on the fly while it is read using a
      function pointer that can be set on the object. The spectra is then
      cached to disk using the functions provided in CachedmzML.
    */
    class OPENMS_DLLAPI CachedMzMLConsumer :
      public CachedmzML,
      public Interfaces::IMSDataConsumer<>
    {
      typedef MSExperiment<> MapType;
      typedef MapType::SpectrumType SpectrumType;
      typedef MapType::ChromatogramType ChromatogramType;

    public:
      /// Default constructor
      CachedMzMLConsumer(String filename, bool clearData=true) :
        ofs(filename.c_str(), std::ios::binary),
        clearData_(clearData),
        spectra_written(0),
        chromatograms_written(0),
        spectra_expected(0),
        chromatograms_expected(0)
      {
      }

      /// Default destructor
      ~CachedMzMLConsumer()
      {
        // Close file stream: close() _should_ call flush() but it might not in
        // all cases. To be sure call flush() first.
        ofs.flush();
        ofs.close();
      }

      /// Write a spectrum
      void consumeSpectrum(SpectrumType & s)
      {
        if (spectra_written >= spectra_expected || chromatograms_written > 0)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                  "Cannot write spectra, reached expected spectra or have already written chromatograms.");
        }
        writeSpectrum_(s, ofs);
        spectra_written++;
        if (clearData_) {s.clear(false);}
      }

      /// Write a chromatogram
      void consumeChromatogram(ChromatogramType & c)
      {
        if (chromatograms_written >= chromatograms_expected || spectra_written != spectra_expected)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                  "Cannot write spectra, reached expected spectra or have already written chromatograms.");
        }
        writeChromatogram_(c, ofs);
        chromatograms_written++;
        if (clearData_) {c.clear(false);}
      }

      /// Write the header of a file to disk
      void setExpectedSize(Size expectedSpectra, Size expectedChromatograms)
      {
        if (spectra_expected != 0)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                  "Can only set expected size of the experiment once since this will open the file.");
        }

        spectra_expected = expectedSpectra;
        chromatograms_expected = expectedChromatograms;

        int magic_number = MAGIC_NUMBER;
        ofs.write((char*)&magic_number, sizeof(magic_number));
        ofs.write((char*)&spectra_expected, sizeof(spectra_expected));
        ofs.write((char*)&chromatograms_expected, sizeof(chromatograms_expected));
      }

      void setExperimentalSettings(const ExperimentalSettings& /* exp */) {;}

    protected:
      std::ofstream ofs;
      bool clearData_;
      Size spectra_written;
      Size chromatograms_written;
      Size spectra_expected;
      Size chromatograms_expected;

    };

} //end namespace OpenMS

#endif
