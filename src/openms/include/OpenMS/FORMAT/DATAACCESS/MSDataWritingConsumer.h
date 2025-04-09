// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/INTERFACES/IMSDataConsumer.h>

#include <OpenMS/FORMAT/HANDLERS/MzMLHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <vector>
#include <string>
#include <fstream>
#include <boost/shared_ptr.hpp>

namespace OpenMS
{
    /**
      @brief Consumer class that writes MS data to disk using the mzML format.

      The MSDataWritingConsumer is able to write spectra and chromatograms to
      disk on the fly (as soon as they are consumed). This class is abstract
      and allows the derived class to define how spectra and chromatograms are
      processed before being written to disk.
      
      If you are looking for class that simply takes spectra and chromatograms
      and writes them to disk, please use the PlainMSDataWritingConsumer.

      Example usage:

      @code
      MSDataWritingConsumer * consumer = new MSDataWritingConsumer(outfile); // some implementation
      consumer->setExpectedSize(specsize, chromsize);
      consumer->setExperimentalSettings(exp_settings);
      consumer->addDataProcessing(dp); // optional, will be added to all spectra and chromatograms
      [...]
      // multiple times ... 
      consumer->consumeSpectrum(spec);
      consumer->consumeChromatogram(chrom);
      [...]
      delete consumer;
      @endcode

      @note The first usage of consumeChromatogram or consumeSpectrum will start
      writing of the mzML header to disk (and the first element).

      @note Currently it is not possible to add spectra after having already added
      chromatograms since this could lead to a situation with multiple
      @a spectrumList tags appear in an mzML file.

      @note The expected size will @a not be enforced but it will lead to an
      inconsistent mzML if the count attribute of spectrumList or
      chromatogramList is incorrect.

    */
    class OPENMS_DLLAPI MSDataWritingConsumer : 
      public Internal::MzMLHandler,
      public Interfaces::IMSDataConsumer
    {

    public:
      typedef PeakMap MapType;
      typedef MapType::SpectrumType SpectrumType;
      typedef MapType::ChromatogramType ChromatogramType;

      /**
        @brief Constructor

        @param filename Filename for the output mzML
      */
      explicit MSDataWritingConsumer(const String& filename);

      /// Destructor
      ~MSDataWritingConsumer() override;

      /// @name IMSDataConsumer interface
      //@{
      /**
        @brief Set experimental settings for the whole file

        @param exp Experimental settings to be used for this file (from this
          and the first spectrum/chromatogram, the class will deduce most of
          the header of the mzML file)
      */
      void setExperimentalSettings(const ExperimentalSettings& exp) override;

      /**
        @brief Set expected size of spectra and chromatograms to be written.

        These numbers will be written in the spectrumList and chromatogramList
        tag in the mzML file. Therefore, these will contain wrong numbers if
        the expected size is not set correctly.

        @param expectedSpectra Number of spectra expected
        @param expectedChromatograms Number of chromatograms expected
      */
      void setExpectedSize(Size expectedSpectra, Size expectedChromatograms) override;

      /**
        @brief Consume a spectrum

        The spectrum will be processed using the processSpectrum_ method of the
        current implementation and then written to the mzML file.

        @param s The spectrum to be written to mzML
      */
      void consumeSpectrum(SpectrumType & s) override;

      /**
        @brief Consume a chromatogram

        The chromatogram will be processed using the processChromatogram_
        method of the current implementation and then written to the mzML file.

        @param c The chromatogram to be written to mzML
      */
      void consumeChromatogram(ChromatogramType & c) override;
      //@}

      /**
        @brief Optionally add a data processing method to each chromatogram and spectrum.

        The provided DataProcessing object will be added to each chromatogram
        and spectrum written to to the mzML file.

        @param d The DataProcessing object to be added
      */
      virtual void addDataProcessing(DataProcessing d);

      /**
        @brief Return the number of spectra written.
      */
      virtual Size getNrSpectraWritten();

      /**
        @brief Return the number of chromatograms written.
      */
      virtual Size getNrChromatogramsWritten();

    private:

      /// @name Data Processing using the template method pattern
      //@{
      /**
        @brief Process a spectrum before storing to disk

        Redefine this function to determine spectra processing before writing to disk.
      */
      virtual void processSpectrum_(SpectrumType & s) = 0;

      /**
        @brief Process a chromatogram before storing to disk

        Redefine this function to determine chromatogram processing before writing to disk.
      */
      virtual void processChromatogram_(ChromatogramType & c) = 0;
      //@}

      /**
        @brief Cleanup function called by the destructor.

        Will write the last tags to the file and close the file stream.
      */
      virtual void doCleanup_();

    protected:

      /// File stream (to write mzML)
      std::ofstream ofs_;

      /// Stores whether we have already started writing any data
      bool started_writing_;
      /// Stores whether we are currently writing spectra
      bool writing_spectra_;
      /// Stores whether we are currently writing chromatograms
      bool writing_chromatograms_;
      /// Number of spectra written
      Size spectra_written_;
      /// Number of chromatograms written
      Size chromatograms_written_;
      /// Number of spectra expected
      Size spectra_expected_;
      /// Number of chromatograms expected
      Size chromatograms_expected_;
      /// Whether to add dataprocessing term to the data before writing
      bool add_dataprocessing_;

      /// Validator that knows about CV terms
      Internal::MzMLValidator * validator_;

      /// Experimental settings to use for the whole file
      ExperimentalSettings settings_;
      /// Vector of data processing objects -> will be filled by writeHeader_
      std::vector<std::vector< ConstDataProcessingPtr > > dps_;
      /// The dataprocessing to be added to each spectrum/chromatogram
      DataProcessingPtr additional_dataprocessing_;
    };

    /**
      @brief Consumer class that writes MS data to disk using the mzML format.

      A very simple implementation of MSDataWritingConsumer, not offering the
      ability to modify the spectra before writing. This is probably the class
      you want if you want to write mzML files to disk by providing spectra and
      chromatograms sequentially.

    */
    class OPENMS_DLLAPI PlainMSDataWritingConsumer :
      public MSDataWritingConsumer 
    {
      void processSpectrum_(MapType::SpectrumType & /* s */) override {}
      void processChromatogram_(MapType::ChromatogramType & /* c */) override {}

    public:

      explicit PlainMSDataWritingConsumer(String filename) : MSDataWritingConsumer(filename) {}
    };

    /**
      @brief Consumer class that perform no operation.

      This is sometimes necessary to fulfill the requirement of passing an
      valid MSDataWritingConsumer object or pointer but no operation is
      required.
    */
    class OPENMS_DLLAPI NoopMSDataWritingConsumer :
      public MSDataWritingConsumer 
    {
    public:

      explicit NoopMSDataWritingConsumer(String filename) : MSDataWritingConsumer(filename) {}
      void setExperimentalSettings(const ExperimentalSettings& /* exp */) override {}
      void consumeSpectrum(SpectrumType & /* s */) override {}
      void consumeChromatogram(ChromatogramType & /* c */) override {}

    private:

      void doCleanup_() override {}
      void processSpectrum_(MapType::SpectrumType & /* s */) override {}
      void processChromatogram_(MapType::ChromatogramType & /* c */) override {}
    };


} //end namespace OpenMS


