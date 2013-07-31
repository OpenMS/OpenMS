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
        @brief Proecess a spectrum before storing to disk
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
        started_writing(false),
        writing_spectra(false),
        writing_chromatograms(false),
        spectra_written(0),
        chromatograms_written(0),
        spectra_expected(0),
        chromatograms_expected(0),
        add_dataprocessing_(false)
      {
        validator_ = new Internal::MzMLValidator(this->mapping_, this->cv_);
        ofs.precision(writtenDigits(DoubleReal()));
      }

      /// Destructor
      virtual ~MSDataWritingConsumer()
      {
        //--------------------------------------------------------------------------------------------
        //cleanup
        //--------------------------------------------------------------------------------------------
        // make sure to close an open List tag
        if (writing_spectra)
        {
          ofs << "\t\t</spectrumList>\n";
        }
        else if (writing_chromatograms)
        {
          ofs << "\t\t</chromatogramList>\n";
        }


        Internal::MzMLHandler<MapType>::writeFooter_(ofs);
        delete validator_;
        ofs.close();
      }

      void setExperimentalSettings(ExperimentalSettings& exp)
      {
        settings = exp;
      }

      void setExpectedSize(Size expectedSpectra, Size expectedChromatograms)
      {
        spectra_expected = expectedSpectra;
        chromatograms_expected = expectedChromatograms;
      }

      void consumeSpectrum(SpectrumType & s)
      {
        if (writing_chromatograms)
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

        if (!started_writing)
        {
          // this is the first spectrum -> start writing the header
          // We also need to modify the map and add this dummy spectrum in
          // order to write the header correctly
          MapType dummy;
          dummy = settings;
          dummy.addSpectrum(scpy);

          //--------------------------------------------------------------------
          //header
          //--------------------------------------------------------------------
          Internal::MzMLHandler<MapType>::writeHeader_(ofs, dummy, dps, *validator_);
          started_writing = true;
        }
        if (!writing_spectra)
        {
          ofs << "\t\t<spectrumList count=\"" << spectra_expected << "\" defaultDataProcessingRef=\"dp_sp_0\">\n";
          writing_spectra = true;
        }
        bool renew_native_ids = false;
        // TODO writeSpectrum assumes that dps has at least one value -> assert
        // this here ...
        Internal::MzMLHandler<MapType>::writeSpectrum_(ofs, scpy,
                spectra_written++, *validator_, renew_native_ids, dps);
      }

      void addDataProcessing(DataProcessing d)
      {
        additional_dataprocessing_ = d;
        add_dataprocessing_ = true;
      }

      void consumeChromatogram(ChromatogramType & c)
      {
        // make sure to close an open List tag
        if (writing_spectra)
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

        if (!started_writing)
        {
          // this is the first chromatogram -> start writing the header
          // We also need to modify the map and add this dummy chromatogram in
          // order to write the header correctly
          MapType dummy;
          dummy = settings;
          dummy.addChromatogram(ccpy);

          //--------------------------------------------------------------------
          //header
          //--------------------------------------------------------------------
          Internal::MzMLHandler<MapType>::writeHeader_(ofs, dummy, dps, *validator_);
          started_writing = true;
        }
        if (!writing_chromatograms)
        {
          ofs << "\t\t<chromatogramList count=\"" << chromatograms_expected << "\" defaultDataProcessingRef=\"dp_sp_0\">\n";
          writing_chromatograms = true;
          writing_spectra = false;
        }
        Internal::MzMLHandler<MapType>::writeChromatogram_(ofs, ccpy,
                chromatograms_written++, *validator_);
      }

    protected:

      std::ofstream ofs;

      bool started_writing;
      bool writing_spectra;
      bool writing_chromatograms;
      Size spectra_written;
      Size chromatograms_written;
      Size spectra_expected;
      Size chromatograms_expected;
      bool add_dataprocessing_;

      /*
      ControlledVocabulary cv_;
      CVMappings mapping_;
      */
      Internal::MzMLValidator * validator_;

      ExperimentalSettings settings;
      std::vector<std::vector<DataProcessing> > dps;
      DataProcessing additional_dataprocessing_;
    };

    class OPENMS_DLLAPI PlainMSDataWritingConsumer :
      public MSDataWritingConsumer 
    {
    public:

      PlainMSDataWritingConsumer(String filename) : MSDataWritingConsumer(filename) {}
      void processSpectrum_(typename MapType::SpectrumType & /* s */) {}
      void processChromatogram_(typename MapType::ChromatogramType & /* c */) {}
    };

} //end namespace OpenMS

#endif
