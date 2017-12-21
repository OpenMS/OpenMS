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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_DATAACCESS_MSDATACACHEDCONSUMER_H
#define OPENMS_FORMAT_DATAACCESS_MSDATACACHEDCONSUMER_H

#include <OpenMS/INTERFACES/IMSDataConsumer.h>

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>

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
      public Interfaces::IMSDataConsumer
    {
      typedef MSSpectrum SpectrumType;
      typedef MSChromatogram ChromatogramType;

    public:

      /**
        @brief Constructor
  
        Opens the output file and writes the header.
      */
      MSDataCachedConsumer(String filename, bool clearData=true) :
        ofs_(filename.c_str(), std::ios::binary),
        clearData_(clearData),
        spectra_written_(0),
        chromatograms_written_(0)
      {
        int file_identifier = CACHED_MZML_FILE_IDENTIFIER;
        ofs_.write((char*)&file_identifier, sizeof(file_identifier));
      }

      /**
        @brief Destructor
  
        Closes the output file and writes the footer.
      */
      ~MSDataCachedConsumer() override
      {
        // Write size of file (to the end of the file)
        ofs_.write((char*)&spectra_written_, sizeof(spectra_written_));
        ofs_.write((char*)&chromatograms_written_, sizeof(chromatograms_written_));

        // Close file stream: close() _should_ call flush() but it might not in
        // all cases. To be sure call flush() first.
        ofs_.flush();
        ofs_.close();
      }

      /**
        @brief Write a spectrum to the output file
      */
      void consumeSpectrum(SpectrumType & s) override
      {
        if (chromatograms_written_ > 0)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Cannot write spectra after writing chromatograms.");
        }
        writeSpectrum_(s, ofs_);
        spectra_written_++;
        if (clearData_) {s.clear(false);}
      }

      /**
        @brief Write a chromatogram to the output file
      */
      void consumeChromatogram(ChromatogramType & c) override
      {
        writeChromatogram_(c, ofs_);
        chromatograms_written_++;
        if (clearData_) {c.clear(false);}
      }

      void setExpectedSize(Size /* expectedSpectra */, Size /* expectedChromatograms */) override {;}

      void setExperimentalSettings(const ExperimentalSettings& /* exp */) override {;}

    protected:
      std::ofstream ofs_;
      bool clearData_;
      Size spectra_written_;
      Size chromatograms_written_;

    };

} //end namespace OpenMS

#endif
