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

#ifndef OPENMS_FORMAT_DATAACCESS_MSDATAWRITINGCONSUMER_H
#define OPENMS_FORMAT_DATAACCESS_MSDATAWRITINGCONSUMER_H

#include <OpenMS/INTERFACES/IMSDataConsumer.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/MzMLHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <vector>
#include <string>
#include <fstream>
#include <boost/shared_ptr.hpp>

namespace OpenMS
{
    /**
      @brief Transforming and writing consumer of MS data

      Is able to transform a spectra on the fly while it is read using a
      function pointer that can be set on the object. The spectra is then
      written to disk using the functions provided in MzMLHandler.

      Example usage:

      @code
      MSDataWritingConsumer * ppConsumer = new MSDataWritingConsumer(outfile); // some implementation
      ppConsumer->setExpectedSize(specsize, chromsize);
      ppConsumer->setExperimentalSettings(exp_settings);
      ppConsumer->addDataProcessing(dp); // optional, will be added to all spectra and chromatograms
      [...]
      ppConsumer->consumeSpectrum(spec);
      ppConsumer->consumeChromatogram(chrom);
      [...]
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
    class OPENMS_DLLAPI MSDataWritingConsumer  : 
      public Internal::MzMLHandler< MSExperiment<> >,
      public Interfaces::IMSDataConsumer< MSExperiment<> >
    {

    public:
      typedef MSExperiment<> MapType;
      typedef MapType::SpectrumType SpectrumType;
      typedef MapType::ChromatogramType ChromatogramType;

      /**
        @brief Process a spectrum or chromatogram before storing to disk
      */
      virtual void processSpectrum_(SpectrumType & s) = 0;

      virtual void processChromatogram_(ChromatogramType & c) = 0;

      /**
        @brief Constructor

        MzMLHandler and the std::ofstream require us to know the output filename in advance.

      */
      MSDataWritingConsumer(String filename) :
        Internal::MzMLHandler<MapType>(MapType(), filename, MzMLFile().getVersion(), ProgressLogger()),
        ofs(filename.c_str()), 
        started_writing_(false),
        writing_spectra_(false),
        writing_chromatograms_(false),
        spectra_written_(0),
        chromatograms_written_(0),
        spectra_expected_(0),
        chromatograms_expected_(0),
        add_dataprocessing_(false)
      {
        validator_ = new Internal::MzMLValidator(this->mapping_, this->cv_);
        ofs.precision(writtenDigits(DoubleReal()));
      }

      /// Destructor
      virtual ~MSDataWritingConsumer()
      {
        doCleanup();
      }

      /// Set experimental settings for the whole file
      void setExperimentalSettings(const ExperimentalSettings& exp)
      {
        settings_ = exp;
      }

      /// Set expected size -> these numbers will be written in the spectrumList and chromatogramList tag
      void setExpectedSize(Size expectedSpectra, Size expectedChromatograms)
      {
        spectra_expected_ = expectedSpectra;
        chromatograms_expected_ = expectedChromatograms;
      }

      void consumeSpectrum(SpectrumType & s)
      {
        if (writing_chromatograms_)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
              "Cannot write spectra after writing chromatograms.");
        }

        // Create copy and add dataprocessing if required
        SpectrumType scpy = s;
        processSpectrum_(scpy);

        if (add_dataprocessing_)
        {
          scpy.getDataProcessing().push_back(additional_dataprocessing_);
        }

        if (!started_writing_)
        {
          // this is the first data to be written -> start writing the header
          // We also need to modify the map and add this dummy spectrum in
          // order to write the header correctly
          MapType dummy;
          dummy = settings_;
          dummy.addSpectrum(scpy);

          //--------------------------------------------------------------------
          //header
          //--------------------------------------------------------------------
          Internal::MzMLHandler<MapType>::writeHeader_(ofs, dummy, dps_, *validator_);
          started_writing_ = true;
        }
        if (!writing_spectra_)
        {
          ofs << "\t\t<spectrumList count=\"" << spectra_expected_ << "\" defaultDataProcessingRef=\"dp_sp_0\">\n";
          writing_spectra_ = true;
        }
        bool renew_native_ids = false;
        // TODO writeSpectrum assumes that dps_ has at least one value -> assert
        // this here ...
        Internal::MzMLHandler<MapType>::writeSpectrum_(ofs, scpy,
                spectra_written_++, *validator_, renew_native_ids, dps_);
      }

      void addDataProcessing(DataProcessing d)
      {
        additional_dataprocessing_ = d;
        add_dataprocessing_ = true;
      }

      void consumeChromatogram(ChromatogramType & c)
      {
        // make sure to close an open List tag
        if (writing_spectra_)
        {
          ofs << "\t\t</spectrumList>\n";
        }

        // Create copy and add dataprocessing if required
        ChromatogramType ccpy = c;
        processChromatogram_(ccpy);

        if (add_dataprocessing_)
        {
          ccpy.getDataProcessing().push_back(additional_dataprocessing_);
        }

        if (!started_writing_)
        {
          // this is the first data to be written -> start writing the header
          // We also need to modify the map and add this dummy chromatogram in
          // order to write the header correctly
          MapType dummy;
          dummy = settings_;
          dummy.addChromatogram(ccpy);

          //--------------------------------------------------------------------
          //header (fill also dps_ variable)
          //--------------------------------------------------------------------
          Internal::MzMLHandler<MapType>::writeHeader_(ofs, dummy, dps_, *validator_);
          started_writing_ = true;
        }
        if (!writing_chromatograms_)
        {
          ofs << "\t\t<chromatogramList count=\"" << chromatograms_expected_ << "\" defaultDataProcessingRef=\"dp_sp_0\">\n";
          writing_chromatograms_ = true;
          writing_spectra_ = false;
        }
        Internal::MzMLHandler<MapType>::writeChromatogram_(ofs, ccpy,
                chromatograms_written_++, *validator_);
      }

      Size getNrSpectraWritten() {return spectra_written_;}
      Size getNrChromatogramsWritten() {return chromatograms_written_;}

    protected:

      void doCleanup()
      {
        //--------------------------------------------------------------------------------------------
        //cleanup
        //--------------------------------------------------------------------------------------------
        // make sure to close an open List tag
        if (writing_spectra_)
        {
          ofs << "\t\t</spectrumList>\n";
        }
        else if (writing_chromatograms_)
        {
          ofs << "\t\t</chromatogramList>\n";
        }

        // Only write the footer if we actually did start writing ... 
        if (started_writing_) 
          Internal::MzMLHandlerHelper::writeFooter_(ofs, options_, spectra_offsets, chromatograms_offsets);

        delete validator_;
        ofs.close();
      }

      std::ofstream ofs;

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
      std::vector<std::vector<DataProcessing> > dps_;
      /// The dataprocessing to be added to each spectrum/chromatogram
      DataProcessing additional_dataprocessing_;
    };

    class OPENMS_DLLAPI PlainMSDataWritingConsumer :
      public MSDataWritingConsumer 
    {
    public:

      PlainMSDataWritingConsumer(String filename) : MSDataWritingConsumer(filename) {}
      void processSpectrum_(MapType::SpectrumType & /* s */) {}
      void processChromatogram_(MapType::ChromatogramType & /* c */) {}
    };

    class OPENMS_DLLAPI NoopMSDataWritingConsumer :
      public MSDataWritingConsumer 
    {
    public:

      NoopMSDataWritingConsumer(String filename) : MSDataWritingConsumer(filename) {}
      void processSpectrum_(MapType::SpectrumType & /* s */) {}
      void processChromatogram_(MapType::ChromatogramType & /* c */) {}
      void setExperimentalSettings(const ExperimentalSettings& /* exp */) {}
      void setExpectedSize(Size /* expectedSpectra */, Size /* expectedChromatograms */) {}
      void consumeSpectrum(SpectrumType & /* s */) {}
      void addDataProcessing(DataProcessing /* d */) {}
      void consumeChromatogram(ChromatogramType & /* c */) {}

    protected:

      void doCleanup() {}
    };

} //end namespace OpenMS

#endif
