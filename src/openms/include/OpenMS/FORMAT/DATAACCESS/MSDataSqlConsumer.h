// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

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

        @param sql_filename The filename of the SQLite database
        @param run_id Unique identifier which links the sqMass and OSW file
        @param buffer_size How large the internal buffer size should be (defaults to 500 spectra / chromatograms)
        @param full_meta Whether to write the full meta-data in the SQLite header
        @param lossy_compression Whether to use lossy compression (numpress)
        @param linear_mass_acc Desired mass accuracy for RT or m/z space (absolute value)
      */
      MSDataSqlConsumer(const String& sql_filename, UInt64 run_id, int buffer_size = 500, bool full_meta = true, bool lossy_compression=false, double linear_mass_acc=1e-4);

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


