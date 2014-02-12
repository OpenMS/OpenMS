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

#include <OpenMS/FORMAT/CachedMzML.h>

namespace OpenMS
{
    /**
      @brief Transforming and cached writing consumer of MS data

      Is able to transform a spectrum on the fly while it is read using a
      function pointer that can be set on the object. The spectra is then
      cached to disk using the functions provided in CachedmzML.
    */
    class OPENMS_DLLAPI MSDataCachedConsumer :
      public CachedmzML,
      public Interfaces::IMSDataConsumer<>
    {
      typedef MSExperiment<> MapType;
      typedef MapType::SpectrumType SpectrumType;
      typedef MapType::ChromatogramType ChromatogramType;

    public:
      /// Default constructor
      MSDataCachedConsumer(String filename, bool clearData=true) :
        ofs_(filename.c_str(), std::ios::binary),
        clearData_(clearData),
        spectra_written_(0),
        chromatograms_written_(0),
        spectra_expected_(0),
        chromatograms_expexted_(0)
      {
      }

      /// Default destructor
      ~MSDataCachedConsumer()
      {
        // Close file stream: close() _should_ call flush() but it might not in
        // all cases. To be sure call flush() first.
        ofs_.flush();
        ofs_.close();
      }

      /// Write a spectrum
      void consumeSpectrum(SpectrumType & s)
      {
        if (spectra_written_ >= spectra_expected_ || chromatograms_written_ > 0)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                  "Cannot write spectra, reached expected spectra or have already written chromatograms.");
        }
        writeSpectrum_(s, ofs_);
        spectra_written_++;
        if (clearData_) {s.clear(false);}
      }

      /// Write a chromatogram
      void consumeChromatogram(ChromatogramType & c)
      {
        if (chromatograms_written_ >= chromatograms_expexted_ || spectra_written_ != spectra_expected_)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                  "Cannot write spectra, reached expected spectra or have already written chromatograms.");
        }
        writeChromatogram_(c, ofs_);
        chromatograms_written_++;
        if (clearData_) {c.clear(false);}
      }

      /// Write the header of a file to disk
      void setExpectedSize(Size expectedSpectra, Size expectedChromatograms)
      {
        if (spectra_expected_ != 0)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                  "Can only set expected size of the experiment once since this will open the file.");
        }

        spectra_expected_ = expectedSpectra;
        chromatograms_expexted_ = expectedChromatograms;

        int file_identifier = CACHED_MZML_FILE_IDENTIFIER;
        ofs_.write((char*)&file_identifier, sizeof(file_identifier));
        ofs_.write((char*)&spectra_expected_, sizeof(spectra_expected_));
        ofs_.write((char*)&chromatograms_expexted_, sizeof(chromatograms_expexted_));
      }

      void setExperimentalSettings(const ExperimentalSettings& /* exp */) {;}

    protected:
      std::ofstream ofs_;
      bool clearData_;
      Size spectra_written_;
      Size chromatograms_written_;
      Size spectra_expected_;
      Size chromatograms_expexted_;

    };

} //end namespace OpenMS

#endif
