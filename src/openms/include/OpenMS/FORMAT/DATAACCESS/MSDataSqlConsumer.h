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

#ifndef OPENMS_FORMAT_DATAACCESS_MSDATASQLCONSUMER_H
#define OPENMS_FORMAT_DATAACCESS_MSDATASQLCONSUMER_H

#include <OpenMS/INTERFACES/IMSDataConsumer.h>

#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{

    namespace Internal
    {
      class MzMLSqliteHandler;
    }

    /**
      @brief A data consumer that inserts MS data into a SQLite database

      Consumes spectra and chromatograms and inserts them into an file-based
      SQL database using SQLite. As SQLite is highly inefficient when inserting
      one spectrum/chromatogram at a time, the consumer collects the data in an
      internal buffer and then flushes them all together to disk.

      It uses MzMLSqliteHandler internally to write batches of data to disk.

    */
    class OPENMS_DLLAPI MSDataSqlConsumer :
      public Interfaces::IMSDataConsumer
    {
      typedef MSExperiment MapType;
      typedef MapType::SpectrumType SpectrumType;
      typedef MapType::ChromatogramType ChromatogramType;

    public:

      /**
        @brief Constructor

        Opens the SQLite file and writes the tables.

        @param filename The filename of the SQLite database
        @param buffer_size How large the internal buffer size should be (defaults to 500 spectra / chromatograms)
        @param full_meta Whether to write the full meta-data in the SQLite header
        @param lossy_compression Whether to use lossy compression (numpress)
        @param linear_mass_acc Desired mass accuracy for RT or m/z space (absolute value)
      */
      MSDataSqlConsumer(String filename, int buffer_size = 500, bool full_meta = true, bool lossy_compression=false, double linear_mass_acc=1e-4);

      /**
        @brief Destructor
  
        Flushes the data for good.
      */
      ~MSDataSqlConsumer() override;

      /**
        @brief Flushes the data for good.

        After calling this function, no more data is held in the buffer but the
        class is still able to receive new data.
      */
      void flush();

      /**
        @brief Write a spectrum to the output file
      */
      void consumeSpectrum(SpectrumType & s) override;

      /**
        @brief Write a chromatogram to the output file
      */
      void consumeChromatogram(ChromatogramType & c) override;

      void setExpectedSize(Size /* expectedSpectra */, Size /* expectedChromatograms */) override;

      void setExperimentalSettings(const ExperimentalSettings& /* exp */) override;

    protected:

      String filename_;
      OpenMS::Internal::MzMLSqliteHandler * handler_;

      size_t flush_after_;
      bool full_meta_;
      std::vector<SpectrumType> spectra_;
      std::vector<ChromatogramType> chromatograms_;

      MSExperiment peak_meta_;
    };

} //end namespace OpenMS

#endif // OPENMS_FORMAT_DATAACCESS_MSDATASQLCONSUMER_H

