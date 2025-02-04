// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Chris Bielow, Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/MzMLHandler.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/MSNumpressCoder.h>
#include <OpenMS/FORMAT/VALIDATORS/MzMLValidator.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/SYSTEM/File.h>

#include <map>

namespace OpenMS::Internal
{

    thread_local ProgressLogger pg_outer; ///< an extra logger for nested logging

    /// Constructor for a read-only handler
    MzMLHandler::MzMLHandler(MapType& exp, const String& filename, const String& version, const ProgressLogger& logger)
      : MzMLHandler(filename, version, logger)
    {
      exp_ = &exp;
    }

    /// Constructor for a write-only handler
    MzMLHandler::MzMLHandler(const MapType& exp, const String& filename, const String& version, const ProgressLogger& logger)
      : MzMLHandler(filename, version, logger)
    {
      cexp_ = &exp;
    }

    /// delegated c'tor for the common things
    MzMLHandler::MzMLHandler(const String& filename, const String& version, const ProgressLogger& logger)
      : XMLHandler(filename, version),
        logger_(logger),
        cv_(ControlledVocabulary::getPSIMSCV())
    {
      CVMappingFile().load(File::find("/MAPPING/ms-mapping.xml"), mapping_);

      // check the version number of the mzML handler
      if (VersionInfo::VersionDetails::create(version_) == VersionInfo::VersionDetails::EMPTY)
      {
        OPENMS_LOG_ERROR << "MzMLHandler was initialized with an invalid version number: " << version_ << std::endl;
      }
      pg_outer = logger; // inherit the logtype etc
    }


    /// Destructor
    MzMLHandler::~MzMLHandler() = default;

    /// Set the peak file options
    void MzMLHandler::setOptions(const PeakFileOptions& opt)
    {
      options_ = opt;
      spectrum_data_.reserve(options_.getMaxDataPoolSize());
      chromatogram_data_.reserve(options_.getMaxDataPoolSize());
    }

    /// Get the peak file options
    PeakFileOptions& MzMLHandler::getOptions()
    {
      return options_;
    }


    /// handler which support partial loading, implement this method
    XMLHandler::LOADDETAIL MzMLHandler::getLoadDetail() const
    {
      return load_detail_;
    }

    /// handler which support partial loading, implement this method
    void MzMLHandler::setLoadDetail(const XMLHandler::LOADDETAIL d)
    {
      load_detail_ = d;
    }

    //@}

    /// Get the spectra and chromatogram counts of a file
    void MzMLHandler::getCounts(Size& spectra_counts, Size& chromatogram_counts)
    {
      if (load_detail_ == XMLHandler::LD_RAWCOUNTS)
      {
        spectra_counts = std::max(scan_count_total_, 0); // default is -1; if no specs were found, report 0
        chromatogram_counts = std::max(chrom_count_total_, 0);
      }
      else
      {
        spectra_counts = scan_count_;
        chromatogram_counts = chromatogram_count_;
      }
    }

    /// Set the IMSDataConsumer consumer which will consume the read data
    void MzMLHandler::setMSDataConsumer(Interfaces::IMSDataConsumer* consumer)
    {
      consumer_ = consumer;
    }

    void MzMLHandler::populateSpectraWithData_()
    {

      // Whether spectrum should be populated with data
      if (options_.getFillData())
      {
        size_t errCount = 0;
        String error_message;
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (SignedSize i = 0; i < (SignedSize)spectrum_data_.size(); i++)
        {
          // parallel exception catching and re-throwing business
          if (!errCount) // no need to parse further if already an error was encountered
          {
            try
            {
              populateSpectraWithData_(spectrum_data_[i].data,
                                       spectrum_data_[i].default_array_length,
                                       options_,
                                       spectrum_data_[i].spectrum);
              if (options_.getSortSpectraByMZ() && !spectrum_data_[i].spectrum.isSorted())
              {
                spectrum_data_[i].spectrum.sortByPosition();
              }
            }

            catch (OpenMS::Exception::BaseException& e)
            {
#pragma omp critical(MZMLErrorHandling)
              {
                ++errCount;
                error_message = e.what();
              }
            }
            catch (...)
            {
#pragma omp atomic
              ++errCount;
            }
          }
        }
        if (errCount != 0)
        {
          std::cerr << "  Parsing error: '" << error_message  << "'" << std::endl;
          std::cerr << "  You could try to disable sorting spectra while loading." << std::endl;
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file_, "Error during parsing of binary data: '" + error_message + "'");
        }
      }

      // Append all spectra to experiment / consumer
      for (Size i = 0; i < spectrum_data_.size(); i++)
      {
        if (consumer_ != nullptr)
        {
          consumer_->consumeSpectrum(spectrum_data_[i].spectrum);
          if (options_.getAlwaysAppendData())
          {
            exp_->addSpectrum(std::move(spectrum_data_[i].spectrum));
          }
        }
        else
        {
          exp_->addSpectrum(std::move(spectrum_data_[i].spectrum));
        }
      }

      // Delete batch
      spectrum_data_.clear();
    }

    void MzMLHandler::populateChromatogramsWithData_()
    {
      // Whether chromatogram should be populated with data
      if (options_.getFillData())
      {
        size_t errCount = 0;
        String error_message;
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (SignedSize i = 0; i < (SignedSize)chromatogram_data_.size(); i++)
        {
          // parallel exception catching and re-throwing business
          try
          {
            populateChromatogramsWithData_(chromatogram_data_[i].data,
                                           chromatogram_data_[i].default_array_length,
                                           options_,
                                           chromatogram_data_[i].chromatogram);
            if (options_.getSortChromatogramsByRT() && !chromatogram_data_[i].chromatogram.isSorted())
            {
              chromatogram_data_[i].chromatogram.sortByPosition();
            }
          }
          catch (OpenMS::Exception::BaseException& e)
          {
#pragma omp critical
            {
              ++errCount;
              error_message = e.what();
            }
          }
          catch (...)
          {
#pragma omp atomic
            ++errCount;
          }
        }
        if (errCount != 0)
        {
          // throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file_, "Error during parsing of binary data.");
          std::cerr << "  Parsing error: '" << error_message  << "'" << std::endl;
          std::cerr << "  You could try to disable sorting spectra while loading." << std::endl;
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file_, "Error during parsing of binary data: '" + error_message + "'");
        }

      }

      // Append all chromatograms to experiment / consumer
      for (Size i = 0; i < chromatogram_data_.size(); i++)
      {
        if (consumer_ != nullptr)
        {
          consumer_->consumeChromatogram(chromatogram_data_[i].chromatogram);
          if (options_.getAlwaysAppendData())
          {
            exp_->addChromatogram(std::move(chromatogram_data_[i].chromatogram));
          }
        }
        else
        {
          exp_->addChromatogram(std::move(chromatogram_data_[i].chromatogram));
        }
      }

      // Delete batch
      chromatogram_data_.clear();
    }

    void MzMLHandler::addSpectrumMetaData_(const std::vector<MzMLHandlerHelper::BinaryData>& input_data,
                                           const Size n,
                                           SpectrumType& spectrum) const
    {

      // add meta data
      UInt meta_float_array_index = 0;
      UInt meta_int_array_index = 0;
      UInt meta_string_array_index = 0;
      for (Size i = 0; i < input_data.size(); i++) //loop over all binary data arrays
      {
        if (input_data[i].meta.getName() != "m/z array" && input_data[i].meta.getName() != "intensity array") // is meta data array?
        {
          if (input_data[i].data_type == MzMLHandlerHelper::BinaryData::DT_FLOAT)
          {
            if (n < input_data[i].size)
            {
              double value = (input_data[i].precision == MzMLHandlerHelper::BinaryData::PRE_64) ? input_data[i].floats_64[n] : input_data[i].floats_32[n];
              spectrum.getFloatDataArrays()[meta_float_array_index].push_back(value);
            }
            ++meta_float_array_index;
          }
          else if (input_data[i].data_type == MzMLHandlerHelper::BinaryData::DT_INT)
          {
            if (n < input_data[i].size)
            {
              Int64 value = (input_data[i].precision == MzMLHandlerHelper::BinaryData::PRE_64) ? input_data[i].ints_64[n] : input_data[i].ints_32[n];
              spectrum.getIntegerDataArrays()[meta_int_array_index].push_back(value);
            }
            ++meta_int_array_index;
          }
          else if (input_data[i].data_type == MzMLHandlerHelper::BinaryData::DT_STRING)
          {
            if (n < input_data[i].decoded_char.size())
            {
              String value = input_data[i].decoded_char[n];
              spectrum.getStringDataArrays()[meta_string_array_index].push_back(value);
            }
            ++meta_string_array_index;
          }
        }
      }
    }

    void MzMLHandler::populateSpectraWithData_(std::vector<MzMLHandlerHelper::BinaryData>& input_data,
                                               Size& default_arr_length,
                                               const PeakFileOptions& peak_file_options,
                                               SpectrumType& spectrum)
    {
      typedef SpectrumType::PeakType PeakType;

      // decode all base64 arrays
      MzMLHandlerHelper::decodeBase64Arrays(input_data, options_.getSkipXMLChecks());

      //look up the precision and the index of the intensity and m/z array
      bool mz_precision_64 = true;
      bool int_precision_64 = true;
      SignedSize mz_index = -1;
      SignedSize int_index = -1;
      MzMLHandlerHelper::computeDataProperties_(input_data, mz_precision_64, mz_index, "m/z array");
      MzMLHandlerHelper::computeDataProperties_(input_data, int_precision_64, int_index, "intensity array");

      //Abort if no m/z or intensity array is present
      if (int_index == -1 || mz_index == -1)
      {
        //if defaultArrayLength > 0 : warn that no m/z or int arrays is present
        if (default_arr_length != 0)
        {
          warning(LOAD, String("The m/z or intensity array of spectrum '") + spectrum.getNativeID() + "' is missing and default_arr_length is " + default_arr_length + ".");
        }
        return;
      }

      // Error if intensity or m/z is encoded as int32|64 - they should be float32|64!
      if ((!input_data[mz_index].ints_32.empty()) || (!input_data[mz_index].ints_64.empty()))
      {
        fatalError(LOAD, "Encoding m/z array as integer is not allowed!");
      }
      if ((!input_data[int_index].ints_32.empty()) || (!input_data[int_index].ints_64.empty()))
      {
        fatalError(LOAD, "Encoding intensity array as integer is not allowed!");
      }

      // Warn if the decoded data has a different size than the defaultArrayLength
      Size mz_size = mz_precision_64 ? input_data[mz_index].floats_64.size() : input_data[mz_index].floats_32.size();
      Size int_size = int_precision_64 ? input_data[int_index].floats_64.size() : input_data[int_index].floats_32.size();
      // Check if int-size and mz-size are equal
      if (mz_size != int_size)
      {
        fatalError(LOAD, String("The length of m/z and integer values of spectrum '") + spectrum.getNativeID() + "' differ (mz-size: " + mz_size + ", int-size: " + int_size + "! Not reading spectrum!");
      }
      bool repair_array_length = false;
      if (default_arr_length != mz_size)
      {
        warning(LOAD, String("The m/z array of spectrum '") + spectrum.getNativeID() + "' has the size " + mz_size + ", but it should have size " + default_arr_length + " (defaultArrayLength).");
        repair_array_length = true;
      }
      if (default_arr_length != int_size)
      {
        warning(LOAD, String("The intensity array of spectrum '") + spectrum.getNativeID() + "' has the size " + int_size + ", but it should have size " + default_arr_length + " (defaultArrayLength).");
        repair_array_length = true;
      }
      if (repair_array_length)
      {
        default_arr_length = int_size;
        warning(LOAD, String("Fixing faulty defaultArrayLength to ") + default_arr_length + ".");
      }

      //create meta data arrays and reserve enough space for the content
      if (input_data.size() > 2)
      {
        for (Size i = 0; i < input_data.size(); i++)
        {
          if (input_data[i].meta.getName() != "m/z array" && input_data[i].meta.getName() != "intensity array")
          {
            if (input_data[i].data_type == MzMLHandlerHelper::BinaryData::DT_FLOAT)
            {
              //create new array
              spectrum.getFloatDataArrays().resize(spectrum.getFloatDataArrays().size() + 1);
              //reserve space in the array
              spectrum.getFloatDataArrays().back().reserve(input_data[i].size);
              //copy meta info into MetaInfoDescription
              spectrum.getFloatDataArrays().back().MetaInfoDescription::operator=(input_data[i].meta);
            }
            else if (input_data[i].data_type == MzMLHandlerHelper::BinaryData::DT_INT)
            {
              //create new array
              spectrum.getIntegerDataArrays().resize(spectrum.getIntegerDataArrays().size() + 1);
              //reserve space in the array
              spectrum.getIntegerDataArrays().back().reserve(input_data[i].size);
              //copy meta info into MetaInfoDescription
              spectrum.getIntegerDataArrays().back().MetaInfoDescription::operator=(input_data[i].meta);
            }
            else if (input_data[i].data_type == MzMLHandlerHelper::BinaryData::DT_STRING)
            {
              //create new array
              spectrum.getStringDataArrays().resize(spectrum.getStringDataArrays().size() + 1);
              //reserve space in the array
              spectrum.getStringDataArrays().back().reserve(input_data[i].decoded_char.size());
              //copy meta info into MetaInfoDescription
              spectrum.getStringDataArrays().back().MetaInfoDescription::operator=(input_data[i].meta);
            }
          }
        }
      }

      // Copy meta data from m/z and intensity binary
      // We don't have this as a separate location => store it in spectrum
      for (Size i = 0; i < input_data.size(); i++)
      {
        if (input_data[i].meta.getName() == "m/z array" || input_data[i].meta.getName() == "intensity array")
        {
          std::vector<UInt> keys;
          input_data[i].meta.getKeys(keys);
          for (Size k = 0; k < keys.size(); ++k)
          {
            spectrum.setMetaValue(keys[k], input_data[i].meta.getMetaValue(keys[k]));
          }
        }
      }

      // We found that the push back approach is about 5% faster than using
      // iterators (e.g. spectrum iterator that gets updated)

      //add the peaks and the meta data to the container (if they pass the restrictions)
      PeakType tmp;
      spectrum.reserve(default_arr_length);

      // the most common case: no ranges, 64 / 32 precision
      //  -> this saves about 10 % load time
      if ( mz_precision_64 && !int_precision_64 &&
           input_data.size() == 2 &&
           !peak_file_options.hasMZRange() &&
           !peak_file_options.hasIntensityRange()
         )
      {
        std::vector< double >::const_iterator mz_it = input_data[mz_index].floats_64.begin();
        std::vector< float >::const_iterator int_it = input_data[int_index].floats_32.begin();
        for (Size n = 0; n < default_arr_length; n++)
        {
          //add peak
          tmp.setIntensity(*int_it);
          tmp.setMZ(*mz_it);
          ++mz_it;
          ++int_it;
          spectrum.push_back(tmp);
        }
        return;
      }

      for (Size n = 0; n < default_arr_length; n++)
      {
        double mz = mz_precision_64 ? input_data[mz_index].floats_64[n] : input_data[mz_index].floats_32[n];
        double intensity = int_precision_64 ? input_data[int_index].floats_64[n] : input_data[int_index].floats_32[n];
        if ((!peak_file_options.hasMZRange() || peak_file_options.getMZRange().encloses(DPosition<1>(mz)))
           && (!peak_file_options.hasIntensityRange() || peak_file_options.getIntensityRange().encloses(DPosition<1>(intensity))))
        {
          //add peak
          tmp.setIntensity(intensity);
          tmp.setMZ(mz);
          spectrum.push_back(tmp);

          // Only if there are more than 2 data arrays, we need to check
          // for meta data (as there will always be an m/z and intensity
          // array)
          if (input_data.size() > 2)
          {
            addSpectrumMetaData_(input_data, n, spectrum);
          }
        }
      }
    }

    void MzMLHandler::populateChromatogramsWithData_(std::vector<MzMLHandlerHelper::BinaryData>& input_data,
                                                     Size& default_arr_length,
                                                     const PeakFileOptions& peak_file_options,
                                                     ChromatogramType& inp_chromatogram)
    {
      typedef ChromatogramType::PeakType ChromatogramPeakType;

      //decode all base64 arrays
      MzMLHandlerHelper::decodeBase64Arrays(input_data, options_.getSkipXMLChecks());

      //look up the precision and the index of the intensity and m/z array
      bool int_precision_64 = true;
      bool rt_precision_64 = true;
      SignedSize int_index = -1;
      SignedSize rt_index = -1;
      MzMLHandlerHelper::computeDataProperties_(input_data, rt_precision_64, rt_index, "time array");
      MzMLHandlerHelper::computeDataProperties_(input_data, int_precision_64, int_index, "intensity array");

      //Abort if no m/z or intensity array is present
      if (int_index == -1 || rt_index == -1)
      {
        //if defaultArrayLength > 0 : warn that no time or int arrays is present
        if (default_arr_length != 0)
        {
          warning(LOAD, String("The time or intensity array of chromatogram '") +
              inp_chromatogram.getNativeID() + "' is missing and default_arr_length is " + default_arr_length + ".");
        }
        return;
      }

      // Warn if the decoded data has a different size than the defaultArrayLength
      Size rt_size = rt_precision_64 ? input_data[rt_index].floats_64.size() : input_data[rt_index].floats_32.size();
      Size int_size = int_precision_64 ? input_data[int_index].floats_64.size() : input_data[int_index].floats_32.size();
      // Check if int-size and rt-size are equal
      if (rt_size != int_size)
      {
        fatalError(LOAD, String("The length of RT and intensity values of chromatogram '") + inp_chromatogram.getNativeID() + "' differ (rt-size: " + rt_size + ", int-size: " + int_size + "! Not reading chromatogram!");
      }
      bool repair_array_length = false;
      if (default_arr_length != rt_size)
      {
        warning(LOAD, String("The base64-decoded rt array of chromatogram '") + inp_chromatogram.getNativeID() + "' has the size " + rt_size + ", but it should have size " + default_arr_length + " (defaultArrayLength).");
        repair_array_length = true;
      }
      if (default_arr_length != int_size)
      {
        warning(LOAD, String("The base64-decoded intensity array of chromatogram '") + inp_chromatogram.getNativeID() + "' has the size " + int_size + ", but it should have size " + default_arr_length + " (defaultArrayLength).");
        repair_array_length = true;
      }
      // repair size of array, accessing memory that is beyond int_size will lead to segfaults later
      if (repair_array_length)
      {
        default_arr_length = int_size; // set to length of actual data (int_size and rt_size are equal, s.a.)
        warning(LOAD, String("Fixing faulty defaultArrayLength to ") + default_arr_length + ".");
      }

      // Create meta data arrays and reserve enough space for the content
      if (input_data.size() > 2)
      {
        for (Size i = 0; i < input_data.size(); i++)
        {
          if (input_data[i].meta.getName() != "intensity array" && input_data[i].meta.getName() != "time array")
          {
            if (input_data[i].data_type == MzMLHandlerHelper::BinaryData::DT_FLOAT)
            {
              //create new array
              inp_chromatogram.getFloatDataArrays().resize(inp_chromatogram.getFloatDataArrays().size() + 1);
              //reserve space in the array
              inp_chromatogram.getFloatDataArrays().back().reserve(input_data[i].size);
              //copy meta info into MetaInfoDescription
              inp_chromatogram.getFloatDataArrays().back().MetaInfoDescription::operator=(input_data[i].meta);
            }
            else if (input_data[i].data_type == MzMLHandlerHelper::BinaryData::DT_INT)
            {
              //create new array
              inp_chromatogram.getIntegerDataArrays().resize(inp_chromatogram.getIntegerDataArrays().size() + 1);
              //reserve space in the array
              inp_chromatogram.getIntegerDataArrays().back().reserve(input_data[i].size);
              //copy meta info into MetaInfoDescription
              inp_chromatogram.getIntegerDataArrays().back().MetaInfoDescription::operator=(input_data[i].meta);
            }
            else if (input_data[i].data_type == MzMLHandlerHelper::BinaryData::DT_STRING)
            {
              //create new array
              inp_chromatogram.getStringDataArrays().resize(inp_chromatogram.getStringDataArrays().size() + 1);
              //reserve space in the array
              inp_chromatogram.getStringDataArrays().back().reserve(input_data[i].decoded_char.size());
              //copy meta info into MetaInfoDescription
              inp_chromatogram.getStringDataArrays().back().MetaInfoDescription::operator=(input_data[i].meta);
            }
          }
        }
      }

      // Copy meta data from time and intensity binary
      // We don't have this as a separate location => store it directly in spectrum meta data
      for (Size i = 0; i < input_data.size(); i++)
      {
        if (input_data[i].meta.getName() == "time array" || input_data[i].meta.getName() == "intensity array")
        {
          std::vector<UInt> keys;
          input_data[i].meta.getKeys(keys);
          for (Size k = 0; k < keys.size(); ++k)
          {
            inp_chromatogram.setMetaValue(keys[k], input_data[i].meta.getMetaValue(keys[k]));
          }
        }
      }

      // Add the peaks and the meta data to the container (if they pass the restrictions)
      inp_chromatogram.reserve(default_arr_length);
      ChromatogramPeakType tmp;
      for (Size n = 0; n < default_arr_length; n++)
      {
        double rt = rt_precision_64 ? input_data[rt_index].floats_64[n] : input_data[rt_index].floats_32[n];
        double intensity = int_precision_64 ? input_data[int_index].floats_64[n] : input_data[int_index].floats_32[n];
        if ((!peak_file_options.hasRTRange() || peak_file_options.getRTRange().encloses(DPosition<1>(rt)))
           && (!peak_file_options.hasIntensityRange() || peak_file_options.getIntensityRange().encloses(DPosition<1>(intensity))))
        {
          //add peak
          tmp.setIntensity(intensity);
          tmp.setRT(rt);
          inp_chromatogram.push_back(tmp);

          //add meta data
          UInt meta_float_array_index = 0;
          UInt meta_int_array_index = 0;
          UInt meta_string_array_index = 0;
          for (Size i = 0; i < input_data.size(); i++) //loop over all binary data arrays
          {
            if (input_data[i].meta.getName() != "intensity array" && input_data[i].meta.getName() != "time array") // is meta data array?
            {
              if (input_data[i].data_type == MzMLHandlerHelper::BinaryData::DT_FLOAT)
              {
                if (n < input_data[i].size)
                {
                  double value = (input_data[i].precision == MzMLHandlerHelper::BinaryData::PRE_64) ? input_data[i].floats_64[n] : input_data[i].floats_32[n];
                  inp_chromatogram.getFloatDataArrays()[meta_float_array_index].push_back(value);
                }
                ++meta_float_array_index;
              }
              else if (input_data[i].data_type == MzMLHandlerHelper::BinaryData::DT_INT)
              {
                if (n < input_data[i].size)
                {
                  Int64 value = (input_data[i].precision == MzMLHandlerHelper::BinaryData::PRE_64) ? input_data[i].ints_64[n] : input_data[i].ints_32[n];
                  inp_chromatogram.getIntegerDataArrays()[meta_int_array_index].push_back(value);
                }
                ++meta_int_array_index;
              }
              else if (input_data[i].data_type == MzMLHandlerHelper::BinaryData::DT_STRING)
              {
                if (n < input_data[i].decoded_char.size())
                {
                  String value = input_data[i].decoded_char[n];
                  inp_chromatogram.getStringDataArrays()[meta_string_array_index].push_back(value);
                }
                ++meta_string_array_index;
              }
            }
          }
        }
      }
    }


    void MzMLHandler::characters(const XMLCh* const chars, const XMLSize_t length)
    {
      if (skip_spectrum_ || skip_chromatogram_)
      {
        return;
      }

      const String& current_tag = open_tags_.back();

      if (current_tag == "binary")
      {
        // Since we convert a Base64 string here, it can only contain plain ASCII
        sm_.appendASCII(chars, length, bin_data_.back().base64);
      }
      else if (current_tag == "offset" || current_tag == "indexListOffset" || current_tag == "fileChecksum")
      {
        //do nothing for
        // - index
        // - checksum
        // - binary chromatogram data
      }
      else
      {
        /*String transcoded_chars2 = sm_.convert(chars);
        transcoded_chars2.trim();
        if (transcoded_chars2 != "")
        {
          warning(LOAD, String("Unhandled character content in tag '") + current_tag + "': " +
              transcoded_chars2);
        }
        */
      }
    }


    void MzMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
    {
      constexpr XMLCh s_count[] = {'c','o','u','n','t', 0};
      constexpr XMLCh s_default_array_length[] = { 'd','e','f','a','u','l','t','A','r','r','a','y','L','e','n','g','t','h' , 0};
      constexpr XMLCh s_array_length[] = { 'a','r','r','a','y','L','e','n','g','t','h' , 0};
      constexpr XMLCh s_accession[] = { 'a','c','c','e','s','s','i','o','n' , 0};
      constexpr XMLCh s_name[] = { 'n','a','m','e' , 0};
      constexpr XMLCh s_type[] = { 't','y','p','e' , 0};
      constexpr XMLCh s_value[] = { 'v','a','l','u','e' , 0};
      constexpr XMLCh s_unit_accession[] = { 'u','n','i','t','A','c','c','e','s','s','i','o','n' , 0};
      constexpr XMLCh s_id[] = { 'i','d' , 0};
      constexpr XMLCh s_ref[] = { 'r','e','f' , 0};
      constexpr XMLCh s_version[] = { 'v','e','r','s','i','o','n' , 0};
      constexpr XMLCh s_version_mzml[] = { 'm','z','M','L',':','v','e','r','s','i','o','n' , 0};
      constexpr XMLCh s_order[] = { 'o','r','d','e','r' , 0};
      constexpr XMLCh s_location[] = { 'l','o','c','a','t','i','o','n' , 0};
      constexpr XMLCh s_sample_ref[] = { 's','a','m','p','l','e','R','e','f' , 0};
      constexpr XMLCh s_software_ref[] = { 's','o','f','t','w','a','r','e','R','e','f' , 0};
      constexpr XMLCh s_source_file_ref[] = { 's','o','u','r','c','e','F','i','l','e','R','e','f' , 0};
      constexpr XMLCh s_spectrum_ref[] = { 's','p','e','c','t','r','u','m','R','e','f' , 0};
      constexpr XMLCh s_default_instrument_configuration_ref[] = { 'd','e','f','a','u','l','t','I','n','s','t','r','u','m','e','n','t','C','o','n','f','i','g','u','r','a','t','i','o','n','R','e','f' , 0};
      constexpr XMLCh s_instrument_configuration_ref[] = { 'i','n','s','t','r','u','m','e','n','t','C','o','n','f','i','g','u','r','a','t','i','o','n','R','e','f' , 0};
      constexpr XMLCh s_default_data_processing_ref[] = { 'd','e','f','a','u','l','t','D','a','t','a','P','r','o','c','e','s','s','i','n','g','R','e','f' , 0};
      constexpr XMLCh s_data_processing_ref[] = { 'd','a','t','a','P','r','o','c','e','s','s','i','n','g','R','e','f' , 0};
      constexpr XMLCh s_start_time_stamp[] = { 's','t','a','r','t','T','i','m','e','S','t','a','m','p' , 0};
      constexpr XMLCh s_external_spectrum_id[] = { 'e','x','t','e','r','n','a','l','S','p','e','c','t','r','u','m','I','D' , 0};
      // constexpr XMLCh s_default_source_file_ref[] = { 'd','e','f','a','u','l','t','S','o','u','r','c','e','F','i','l','e','R','e','f' , 0};
      constexpr XMLCh s_scan_settings_ref[] = { 's','c','a','n','S','e','t','t','i','n','g','s','R','e','f' , 0};
      String tag = sm_.convert(qname);
      open_tags_.push_back(tag);

      // do nothing until a spectrum/chromatogram/spectrumList ends
      if (skip_spectrum_ || skip_chromatogram_)
      {
        return;
      }

      //determine parent tag
      String parent_tag;
      if (open_tags_.size() > 1)
      {
        parent_tag = *(open_tags_.end() - 2);
      }
      String parent_parent_tag;
      if (open_tags_.size() > 2)
      {
        parent_parent_tag = *(open_tags_.end() - 3);
      }

      if (tag == "spectrum")
      {
        // for cppcheck
        constexpr XMLCh s_spot_id[] = { 's','p','o','t','I','D', 0 };

        //number of peaks
        spec_ = SpectrumType();
        default_array_length_ = attributeAsInt_(attributes, s_default_array_length);
        //spectrum source file
        String source_file_ref;
        if (optionalAttributeAsString_(source_file_ref, attributes, s_source_file_ref))
        {
          if (source_files_.find(source_file_ref) != source_files_.end())
          {
            spec_.setSourceFile(source_files_[source_file_ref]);
          }
          else
          {
            OPENMS_LOG_WARN << "Error: unregistered source file reference " << source_file_ref << "." << std::endl;
          }
        }
        //native id
        spec_.setNativeID(attributeAsString_(attributes, s_id));
        //maldi spot id
        String maldi_spot_id;
        if (optionalAttributeAsString_(maldi_spot_id, attributes, s_spot_id))
        {
          spec_.setMetaValue("maldi_spot_id", maldi_spot_id);
        }
        //data processing
        String data_processing_ref;
        if (optionalAttributeAsString_(data_processing_ref, attributes, s_data_processing_ref))
        {
          spec_.setDataProcessing(processing_[data_processing_ref]);
        }
        else
        {
          spec_.setDataProcessing(processing_[default_processing_]);
        }
      }
      else if (tag == "chromatogram")
      {
        if (load_detail_ == XMLHandler::LD_COUNTS_WITHOPTIONS)
        { //, but we only want to count
          skip_chromatogram_ = true; // skip the remaining chrom, until endElement(chromatogram)
          ++chromatogram_count_;
        }

        chromatogram_ = ChromatogramType();
        default_array_length_ = attributeAsInt_(attributes, s_default_array_length);
        String source_file_ref;
        if (optionalAttributeAsString_(source_file_ref, attributes, s_source_file_ref))
        {
          chromatogram_.setSourceFile(source_files_[source_file_ref]);
        }
        // native id
        chromatogram_.setNativeID(attributeAsString_(attributes, s_id));
        // data processing
        String data_processing_ref;
        if (optionalAttributeAsString_(data_processing_ref, attributes, s_data_processing_ref))
        {
          chromatogram_.setDataProcessing(processing_[data_processing_ref]);
        }
        else
        {
          chromatogram_.setDataProcessing(processing_[default_processing_]);
        }
      }
      else if (tag == "spectrumList")
      {
        //default data processing
        default_processing_ = attributeAsString_(attributes, s_default_data_processing_ref);

        //Abort if we need meta data only
        if (options_.getMetadataOnly())
        {
          throw EndParsingSoftly(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
        }
        scan_count_total_ = attributeAsInt_(attributes, s_count);
        logger_.startProgress(0, scan_count_total_, "loading spectra list");
        in_spectrum_list_ = true;
        // we only want total scan count and chrom count
        if (load_detail_ == XMLHandler::LD_RAWCOUNTS)
        { // in case chromatograms came before spectra, we have all information --> end parsing
          if (chrom_count_total_ != -1) throw EndParsingSoftly(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
          // or skip the remaining spectra until </spectrumList>
          skip_spectrum_ = true;
        }
        else
        {
          exp_->reserveSpaceSpectra(scan_count_total_);
        }
      }
      else if (tag == "chromatogramList")
      {
        // default data processing
        default_processing_ = attributeAsString_(attributes, s_default_data_processing_ref);

        //Abort if we need meta data only
        if (options_.getMetadataOnly())
        {
          throw EndParsingSoftly(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
        }
        chrom_count_total_ = attributeAsInt_(attributes, s_count);
        logger_.startProgress(0, chrom_count_total_, "loading chromatogram list");
        in_spectrum_list_ = false;

        // we only want total scan count and chrom count
        if (load_detail_ == XMLHandler::LD_RAWCOUNTS)
        { // in case spectra came before chroms, we have all information --> end parsing
          if (scan_count_total_ != -1)
          {
            throw EndParsingSoftly(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
          }
          // or skip the remaining chroms until </chromatogramList>
          skip_chromatogram_ = true;
        }
        else
        {
          exp_->reserveSpaceChromatograms(chrom_count_total_);
        }
      }
      else if (tag == "binaryDataArrayList" /* && in_spectrum_list_*/)
      {
        bin_data_.reserve(attributeAsInt_(attributes, s_count));
      }
      else if (tag == "binaryDataArray" /* && in_spectrum_list_*/)
      {
        bin_data_.emplace_back();
        bin_data_.back().np_compression = MSNumpressCoder::NONE; // ensure that numpress compression is initially set to none ...
        bin_data_.back().compression = false; // ensure that zlib compression is initially set to none ...

        // array length
        Int array_length = (Int) default_array_length_;
        optionalAttributeAsInt_(array_length, attributes, s_array_length);
        bin_data_.back().size = array_length;

        // data processing
        String data_processing_ref;
        if (optionalAttributeAsString_(data_processing_ref, attributes, s_data_processing_ref))
        {
          bin_data_.back().meta.setDataProcessing(processing_[data_processing_ref]);
        }
      }
      else if (tag == "cvParam")
      {
        String value = "";
        optionalAttributeAsString_(value, attributes, s_value);
        String unit_accession = "";
        optionalAttributeAsString_(unit_accession, attributes, s_unit_accession);
        handleCVParam_(parent_parent_tag, parent_tag, attributeAsString_(attributes, s_accession), attributeAsString_(attributes, s_name), value, unit_accession);
      }
      else if (tag == "userParam")
      {
        String type = "";
        optionalAttributeAsString_(type, attributes, s_type);
        String value = "";
        optionalAttributeAsString_(value, attributes, s_value);
        String unit_accession = "";
        optionalAttributeAsString_(unit_accession, attributes, s_unit_accession);
        handleUserParam_(parent_parent_tag, parent_tag, attributeAsString_(attributes, s_name), type, value, unit_accession);
      }
      else if (tag == "referenceableParamGroup")
      {
        current_id_ = attributeAsString_(attributes, s_id);
      }
      else if (tag == "sourceFile")
      {
        current_id_ = attributeAsString_(attributes, s_id);
        // Name of the source file, without reference to location (either URI or local path). e.g. "control.mzML"
        String name_of_file = attributeAsString_(attributes, s_name);

        //URI-formatted location where the file was retrieved.
        String path_to_file = attributeAsString_(attributes, s_location);

        // mzML files often deviate from the specification by storing e.g. the full path in the name attribute etc.
        // error: whole path is stored in file name. fix: split into path and file name
        if (path_to_file.empty() && !name_of_file.empty())
        {
          path_to_file = File::path(name_of_file);
          name_of_file = File::basename(name_of_file);
          if (path_to_file == ".")
          {
            path_to_file = "file://./";
          }
        }

        // format URI prefix as in mzML spec.
        if (path_to_file.hasPrefix("File://"))
        {
          path_to_file.substitute("File://", "file://");
        }
        if (path_to_file.hasPrefix("FILE://"))
        {
          path_to_file.substitute("FILE://", "file://");
        }
        if (path_to_file.hasPrefix("file:///."))
        {
          path_to_file.substitute("file:///.", "file://./");
        }

        bool is_relative_path = path_to_file.hasPrefix("file://./") || path_to_file.hasPrefix("file://../");

        // ill formed absolute or relative path
        if (!is_relative_path && path_to_file.hasPrefix("file://") && !path_to_file.hasPrefix("file:///"))
        {
          warning(LOAD, "Ill formed absolute or relative sourceFile path: " + path_to_file);
        }

        // if possible convert relative path to absolute path
        if (is_relative_path && File::isDirectory(path_to_file))
        {
          String normal_path = String(path_to_file).substitute("file://", ""); // remove URI prefix
          path_to_file = String("file://") + File::absolutePath(normal_path); // on linux this e.g. file:///home... on win: file://C:/...
        }

        // absolute path to the root: remove additional / otherwise we will get file://// on concatenation
        if (!is_relative_path && path_to_file == "file:///")
        {
          path_to_file = "file://";
        }

        source_files_[current_id_].setNameOfFile(name_of_file);
        source_files_[current_id_].setPathToFile(path_to_file);
      }
      else if (tag == "referenceableParamGroupRef")
      {
        //call handleCVParam_ with the parent tag for each parameter in the group
        String ref = attributeAsString_(attributes, s_ref);
        for (Size i = 0; i < ref_param_[ref].size(); ++i)
        {
          handleCVParam_(parent_parent_tag, parent_tag, ref_param_[ref][i].accession, ref_param_[ref][i].name, ref_param_[ref][i].value, ref_param_[ref][i].unit_accession);
        }
      }
      else if (tag == "scan")
      {
        Acquisition tmp;
        //source file => meta data
        String source_file_ref;
        if (optionalAttributeAsString_(source_file_ref, attributes, s_source_file_ref))
        {
          tmp.setMetaValue("source_file_name", source_files_[source_file_ref].getNameOfFile());
          tmp.setMetaValue("source_file_path", source_files_[source_file_ref].getPathToFile());
        }
        //external spectrum id => meta data
        String external_spectrum_id;
        if (optionalAttributeAsString_(external_spectrum_id, attributes, s_external_spectrum_id))
        {
          tmp.setIdentifier(external_spectrum_id);
        }

        //spectrumRef - not really needed

        //instrumentConfigurationRef - not really needed: why should a scan have a different instrument?
        String instrument_configuration_ref;
        if (optionalAttributeAsString_(instrument_configuration_ref, attributes, s_instrument_configuration_ref))
        {
          warning(LOAD, "Unhandled attribute 'instrumentConfigurationRef' in 'scan' tag.");
        }

        spec_.getAcquisitionInfo().push_back(std::move(tmp));
      }
      else if (tag == "mzML")
      {
        scan_count_ = 0;
        chromatogram_count_ = 0;
        scan_count_total_ = -1;
        chrom_count_total_ = -1;


        //check file version against schema version
        String file_version;
        if (!(optionalAttributeAsString_(file_version, attributes, s_version) || optionalAttributeAsString_(file_version, attributes, s_version_mzml)) )
        {
          warning(LOAD, "No version attribute in mzML");
        }

        VersionInfo::VersionDetails current_version = VersionInfo::VersionDetails::create(file_version);
        static VersionInfo::VersionDetails mzML_min_version = VersionInfo::VersionDetails::create("1.1.0");

        if (current_version == VersionInfo::VersionDetails::EMPTY)
        {
          warning(LOAD, String("Invalid mzML version string '") + file_version + "'. Assuming mzML version " + version_ + "!");
        }
        else
        {
          if (current_version < mzML_min_version)
          {
            fatalError(LOAD, String("Only mzML 1.1.0 or higher is supported! This file has version '") + file_version + "'.");
          }
          else if (current_version > VersionInfo::VersionDetails::create(version_))
          {
            warning(LOAD, "The mzML file version (" + file_version + ") is newer than the parser version (" + version_ + "). This might lead to undefined behavior.");
          }
        }

        //handle file accession
        String accession;
        if (optionalAttributeAsString_(accession, attributes, s_accession))
        {
          exp_->setIdentifier(accession);
        }
        //handle file id
        String id;
        if (optionalAttributeAsString_(id, attributes, s_id))
        {
          exp_->setMetaValue("mzml_id", id);
        }
        pg_outer.startProgress(0, 1, "loading mzML");
      }
      else if (tag == "contact")
      {
        exp_->getContacts().emplace_back();
      }
      else if (tag == "sample")
      {
        current_id_ = attributeAsString_(attributes, s_id);
        String name;
        if (optionalAttributeAsString_(name, attributes, s_name))
        {
          samples_[current_id_].setName(name);
        }
      }
      else if (tag == "run")
      {
        //sample
        String sample_ref;
        if (optionalAttributeAsString_(sample_ref, attributes, s_sample_ref))
        {
          exp_->setSample(samples_[sample_ref]);
        }
        //instrument
        String instrument_ref = attributeAsString_(attributes, s_default_instrument_configuration_ref);
        exp_->setInstrument(instruments_[instrument_ref]);
        //start time
        String start_time;
        if (optionalAttributeAsString_(start_time, attributes, s_start_time_stamp))
        {
          exp_->setDateTime(asDateTime_(start_time));
        }
        /*
        //defaultSourceFileRef
        String default_source_file_ref;
        if (optionalAttributeAsString_(default_source_file_ref, attributes, s_default_source_file_ref))
        {
          exp_->getSourceFiles().push_back(source_files_[default_source_file_ref]);
        } 
        */       
      }
      else if (tag == "software")
      {
        current_id_ = attributeAsString_(attributes, s_id);
        software_[current_id_].setVersion(attributeAsString_(attributes, s_version));
      }
      else if (tag == "dataProcessing")
      {
        current_id_ = attributeAsString_(attributes, s_id);
      }
      else if (tag == "processingMethod")
      {
        DataProcessingPtr dp(new DataProcessing);
        // See ticket 452: Do NOT remove this try/catch block until foreign
        // software (e.g. ProteoWizard msconvert.exe) produces valid mzML.
        try
        {
          dp->setSoftware(software_[attributeAsString_(attributes, s_software_ref)]);
        }
        catch (Exception::ParseError& /*e*/)
        {
          OPENMS_LOG_ERROR << "Warning: Parsing error, \"processingMethod\" is missing the required attribute \"softwareRef\".\n" <<
          "The software tool which generated this mzML should be fixed. Please notify the maintainers." << std::endl;
        }
        processing_[current_id_].push_back(dp);
        //The order of processing methods is currently ignored
      }
      else if (tag == "instrumentConfiguration")
      {
        current_id_ = attributeAsString_(attributes, s_id);

        //scan settings
        String scan_settings_ref;
        if (optionalAttributeAsString_(scan_settings_ref, attributes, s_scan_settings_ref))
        {
          warning(LOAD, "Unhandled attribute 'scanSettingsRef' in 'instrumentConfiguration' tag.");
        }
      }
      else if (tag == "softwareRef")
      {
        //Set the software of the instrument
        instruments_[current_id_].setSoftware(software_[attributeAsString_(attributes, s_ref)]);
      }
      else if (tag == "source")
      {
        instruments_[current_id_].getIonSources().emplace_back();
        instruments_[current_id_].getIonSources().back().setOrder(attributeAsInt_(attributes, s_order));
      }
      else if (tag == "analyzer")
      {
        instruments_[current_id_].getMassAnalyzers().emplace_back();
        instruments_[current_id_].getMassAnalyzers().back().setOrder(attributeAsInt_(attributes, s_order));
      }
      else if (tag == "detector")
      {
        instruments_[current_id_].getIonDetectors().emplace_back();
        instruments_[current_id_].getIonDetectors().back().setOrder(attributeAsInt_(attributes, s_order));
      }
      else if (tag == "precursor")
      {
        if (in_spectrum_list_)
        {
          //initialize
          spec_.getPrecursors().emplace_back();

          //source file => meta data
          String source_file_ref;
          if (optionalAttributeAsString_(source_file_ref, attributes, s_source_file_ref))
          {
            spec_.getPrecursors().back().setMetaValue("source_file_name", source_files_[source_file_ref].getNameOfFile());
            spec_.getPrecursors().back().setMetaValue("source_file_path", source_files_[source_file_ref].getPathToFile());
          }
          //external spectrum id => meta data
          String external_spectrum_id;
          if (optionalAttributeAsString_(external_spectrum_id, attributes, s_external_spectrum_id))
          {
            spec_.getPrecursors().back().setMetaValue("external_spectrum_id", external_spectrum_id);
          }

          //spectrum_ref => meta data
          String spectrum_ref;
          if (optionalAttributeAsString_(spectrum_ref, attributes, s_spectrum_ref))
          {
            spec_.getPrecursors().back().setMetaValue("spectrum_ref",  spectrum_ref);
          }
          //reset selected ion count
          selected_ion_count_ = 0;
        }
        else
        {
          chromatogram_.setPrecursor(Precursor());

          String source_file_ref;
          if (optionalAttributeAsString_(source_file_ref, attributes, s_source_file_ref))
          {
            chromatogram_.getPrecursor().setMetaValue("source_file_name", source_files_[source_file_ref].getNameOfFile());
            chromatogram_.getPrecursor().setMetaValue("source_file_path", source_files_[source_file_ref].getPathToFile());
          }

          String external_spectrum_id;
          if (optionalAttributeAsString_(external_spectrum_id, attributes, s_external_spectrum_id))
          {
            chromatogram_.getPrecursor().setMetaValue("external_spectrum_id", external_spectrum_id);
          }
          selected_ion_count_ = 0;
        }
      }
      else if (tag == "product")
      {
        //initialize
        if (in_spectrum_list_)
        {
          spec_.getProducts().emplace_back();
        }
        else
        {
          chromatogram_.setProduct(Product());
        }
      }
      else if (tag == "selectedIon")
      {
        //increase selected ion count
        ++selected_ion_count_;
      }
      else if (tag == "selectedIonList")
      {
        //Warn if more than one selected ion is present
        if (attributeAsInt_(attributes, s_count) > 1)
        {
          warning(LOAD, "OpenMS can currently handle only one selection ion per precursor! Only the first ion is loaded!");
        }
      }
      else if (tag == "scanWindow")
      {
        spec_.getInstrumentSettings().getScanWindows().emplace_back();
      }
    }

    void MzMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
    {
      constexpr XMLCh s_spectrum[] = { 's','p','e','c','t','r','u','m' , 0};
      constexpr XMLCh s_chromatogram[] = { 'c','h','r','o','m','a','t','o','g','r','a','m' , 0};
      constexpr XMLCh s_spectrum_list[] = { 's','p','e','c','t','r','u','m','L','i','s','t' , 0};
      constexpr XMLCh s_chromatogram_list[] = { 'c','h','r','o','m','a','t','o','g','r','a','m','L','i','s','t' , 0};
      constexpr XMLCh s_mzml[] = { 'm','z','M','L' , 0};
      constexpr XMLCh s_sourceFileList[] = { 's','o','u','r','c','e','F','i','l','e','L','i','s','t', 0};

      open_tags_.pop_back();

      if (equal_(qname, s_spectrum))
      {
        if (!skip_spectrum_)
        {
          // catch errors stemming from confusion about elution time and scan time
          if (!rt_set_ && spec_.metaValueExists("elution time (seconds)"))
          {
            spec_.setRT(spec_.getMetaValue("elution time (seconds)"));
          }
          /* this is too hot (could be SRM as well? -- check!):
          // correct spectrum type if possible (i.e., make it more specific)
          if (spec_.getInstrumentSettings().getScanMode() == InstrumentSettings::MASSSPECTRUM)
          {
          if (spec_.getMSLevel() <= 1) spec_.getInstrumentSettings().setScanMode(InstrumentSettings::MS1SPECTRUM);
          else                         spec_.getInstrumentSettings().setScanMode(InstrumentSettings::MSNSPECTRUM);
          }
          */

          // Move current data to (temporary) spectral data object
          SpectrumData tmp;
          tmp.spectrum = std::move(spec_);
          tmp.default_array_length = default_array_length_;
          if (options_.getFillData())
          {
            tmp.data = std::move(bin_data_);
          }
          // append current spectral data to buffer
          spectrum_data_.push_back(std::move(tmp));

          if (spectrum_data_.size() >= options_.getMaxDataPoolSize())
          {
            populateSpectraWithData_();
          }
        }

        switch (load_detail_)
        {
          case XMLHandler::LD_ALLDATA:
          case XMLHandler::LD_COUNTS_WITHOPTIONS:
            skip_spectrum_ = false; // don't skip the next spectrum (unless via options later)
            break;
          case XMLHandler::LD_RAWCOUNTS:
            skip_spectrum_ = true; // we always skip spectra; we only need the outer <spectrumList/chromatogramList count=...>
            break;
        }

        rt_set_ = false;
        logger_.nextProgress();
        bin_data_.clear();
        default_array_length_ = 0;
      }
      else if (equal_(qname, s_chromatogram))
      {
        if (!skip_chromatogram_)
        {

          // Move current data to (temporary) spectral data object
          ChromatogramData tmp;
          tmp.default_array_length = default_array_length_;
          tmp.chromatogram = std::move(chromatogram_);
          if (options_.getFillData())
          {
            tmp.data = std::move(bin_data_);
          }
          // append current spectral data to buffer
          chromatogram_data_.push_back(std::move(tmp));

          if (chromatogram_data_.size() >= options_.getMaxDataPoolSize())
          {
            populateChromatogramsWithData_();
          }
        }

        switch (load_detail_)
        {
          case XMLHandler::LD_ALLDATA:
          case XMLHandler::LD_COUNTS_WITHOPTIONS:
            skip_chromatogram_ = false; // don't skip the next chrom
            break;
          case XMLHandler::LD_RAWCOUNTS:
            skip_chromatogram_ = true; // we always skip chroms; we only need the outer <spectrumList/chromatogramList count=...>
            break;
        }

        logger_.nextProgress();
        bin_data_.clear();
        default_array_length_ = 0;
      }
      else if (equal_(qname, s_spectrum_list))
      {
        skip_spectrum_ = false; // no more spectra to come, so stop skipping (for the LD_RAWCOUNTS case)
        in_spectrum_list_ = false;
        logger_.endProgress();
      }
      else if (equal_(qname, s_chromatogram_list))
      {
        skip_chromatogram_ = false; // no more chromatograms to come, so stop skipping
        in_spectrum_list_ = false;
        logger_.endProgress();
      }
      else if (equal_(qname, s_sourceFileList ))
      {        
        for (auto const& ref_sourcefile : source_files_)
        {
          auto& sfs = exp_->getSourceFiles();
          // only store source files once
          if (std::find(sfs.begin(), sfs.end(), ref_sourcefile.second) == sfs.end())
          {
            exp_->getSourceFiles().push_back(ref_sourcefile.second);
          }
        }
      }
      else if (equal_(qname, s_mzml))
      {
        ref_param_.clear();
        current_id_ = "";
        source_files_.clear();
        samples_.clear();
        software_.clear();
        instruments_.clear();
        processing_.clear();

        // Flush the remaining data
        populateSpectraWithData_();
        populateChromatogramsWithData_();
        pg_outer.endProgress(File::fileSize(file_)); // we cannot query the offset within the file when SAX'ing it (Xerces does not support that)
                                                     // , so we can only report I/O at the very end
      }
    }

    void MzMLHandler::handleCVParam_(const String& parent_parent_tag,
                                     const String& parent_tag,
                                     const String& accession,
                                     const String& name,
                                     const String& value,
                                     const String& unit_accession)
    {
      // the actual value stored in the CVParam
      DataValue termValue = XMLHandler::cvParamToValue(cv_, parent_tag, accession, name, value, unit_accession);
      
      if (termValue == DataValue::EMPTY) return; // conversion failed (warning message was emitted in cvParamToValue())

      //------------------------- run ----------------------------
      if (parent_tag == "run")
      {
        //MS:1000857 ! run attribute
        if (accession == "MS:1000858") //fraction identifier
        {
          exp_->setFractionIdentifier(value);
        }
        else
        {
          warning(LOAD, String("Unhandled cvParam '") + accession + "' in tag '" + parent_tag + "'.");
        }
      }
      //------------------------- binaryDataArray ----------------------------
      else if (parent_tag == "binaryDataArray")
      {
        // store name for all non-default arrays
        if (cv_.isChildOf(accession, "MS:1000513")) // other array names as string
        {
          bin_data_.back().meta.setName(cv_.getTerm(accession).name);
        }

        if (!MzMLHandlerHelper::handleBinaryDataArrayCVParam(bin_data_, accession, value, name, unit_accession))
        {
          if (!cv_.isChildOf(accession, "MS:1000513")) //other array names as string
          {
            warning(LOAD, String("Unhandled cvParam '") + accession + "' in tag '" + parent_tag + "'.");
          }
        }
      }
      //------------------------- spectrum ----------------------------
      else if (parent_tag == "spectrum")
      {
        //spectrum type
        if (accession == "MS:1000294") //mass spectrum
        {
          spec_.getInstrumentSettings().setScanMode(InstrumentSettings::MASSSPECTRUM);
        }
        else if (accession == "MS:1000579") //MS1 spectrum
        {
          spec_.getInstrumentSettings().setScanMode(InstrumentSettings::MS1SPECTRUM);
        }
        else if (accession == "MS:1000580") //MSn spectrum
        {
          spec_.getInstrumentSettings().setScanMode(InstrumentSettings::MSNSPECTRUM);
        }
        else if (accession == "MS:1000581") //CRM spectrum
        {
          spec_.getInstrumentSettings().setScanMode(InstrumentSettings::CRM);
        }
        else if (accession == "MS:1000582") //SIM spectrum
        {
          spec_.getInstrumentSettings().setScanMode(InstrumentSettings::SIM);
        }
        else if (accession == "MS:1000583") //SRM spectrum
        {
          spec_.getInstrumentSettings().setScanMode(InstrumentSettings::SRM);
        }
        else if (accession == "MS:1000804") //electromagnetic radiation spectrum
        {
          spec_.getInstrumentSettings().setScanMode(InstrumentSettings::EMR);
        }
        else if (accession == "MS:1000805") //emission spectrum
        {
          spec_.getInstrumentSettings().setScanMode(InstrumentSettings::EMISSION);
        }
        else if (accession == "MS:1000806") //absorption spectrum
        {
          spec_.getInstrumentSettings().setScanMode(InstrumentSettings::ABSORPTION);
        }
        else if (accession == "MS:1000325") //constant neutral gain spectrum
        {
          spec_.getInstrumentSettings().setScanMode(InstrumentSettings::CNG);
        }
        else if (accession == "MS:1000326") //constant neutral loss spectrum
        {
          spec_.getInstrumentSettings().setScanMode(InstrumentSettings::CNL);
        }
        else if (accession == "MS:1000341") //precursor ion spectrum
        {
          spec_.getInstrumentSettings().setScanMode(InstrumentSettings::PRECURSOR);
        }
        else if (accession == "MS:1000789") //enhanced multiply charged spectrum
        {
          spec_.getInstrumentSettings().setScanMode(InstrumentSettings::EMC);
        }
        else if (accession == "MS:1000790") //time-delayed fragmentation spectrum
        {
          spec_.getInstrumentSettings().setScanMode(InstrumentSettings::TDF);
        }
        //spectrum representation
        else if (accession == "MS:1000127") //centroid spectrum
        {
          spec_.setType(SpectrumSettings::CENTROID);
        }
        else if (accession == "MS:1000128") //profile spectrum
        {
          spec_.setType(SpectrumSettings::PROFILE);
        }
        else if (accession == "MS:1000525") //spectrum representation
        {
          spec_.setType(SpectrumSettings::UNKNOWN);
        }
        // spectrum attribute
        else if (accession == "MS:1000511") //ms level
        {
          spec_.setMSLevel(value.toInt());

          if (options_.hasMSLevels() && !options_.containsMSLevel(spec_.getMSLevel()))
          {
            skip_spectrum_ = true;
          }
          else
          { // MS level is ok
            if (load_detail_ == XMLHandler::LD_COUNTS_WITHOPTIONS)
            { //, and we only want to count
              // , but do not skip the spectrum yet if (load_detail_ == XMLHandler::LD_COUNTS_WITHOPTIONS), since it might be outside the RT range (so should not count)
              //skip_spectrum_ = false; // it is false right now... keep it that way
            }
          }

        }
        else if (accession == "MS:1000497") // deprecated: zoom scan is now a scan attribute
        {
          OPENMS_LOG_DEBUG << "MS:1000497 - zoom scan is now a scan attribute. Reading it for backwards compatibility reasons as spectrum attribute." 
                           << " You can make this warning go away by converting this file using FileConverter to a newer version of the PSI ontology."
                           << " Or by using a recent converter that supports the newest PSI ontology."
                           << std::endl;
          spec_.getInstrumentSettings().setZoomScan(true);
        }
        else if (accession == "MS:1000285") //total ion current
        {
          //No member => meta data
          spec_.setMetaValue("total ion current", termValue);
        }
        else if (accession == "MS:1000504") //base peak m/z
        {
          //No member => meta data
          spec_.setMetaValue("base peak m/z", termValue);
        }
        else if (accession == "MS:1000505") //base peak intensity
        {
          //No member => meta data
          spec_.setMetaValue("base peak intensity", termValue);
        }
        else if (accession == "MS:1000527") //highest observed m/z
        {
          //No member => meta data
          spec_.setMetaValue("highest observed m/z", termValue);
        }
        else if (accession == "MS:1000528") //lowest observed m/z
        {
          //No member => meta data
          spec_.setMetaValue("lowest observed m/z", termValue);
        }
        else if (accession == "MS:1000618") //highest observed wavelength
        {
          //No member => meta data
          spec_.setMetaValue("highest observed wavelength", termValue);
        }
        else if (accession == "MS:1000619") //lowest observed wavelength
        {
          //No member => meta data
          spec_.setMetaValue("lowest observed wavelength", termValue);
        }
        else if (accession == "MS:1000796") //spectrum title
        {
          //No member => meta data
          spec_.setMetaValue("spectrum title", termValue);
        }
        else if (accession == "MS:1000797") //peak list scans
        {
          //No member => meta data
          spec_.setMetaValue("peak list scans", termValue);
        }
        else if (accession == "MS:1000798") //peak list raw scans
        {
          //No member => meta data
          spec_.setMetaValue("peak list raw scans", termValue);
        }
        else if (accession == "MS:1001581") //FAIMS compensation voltage
        {
          // According to the PSI-MS ontology this term should be stored below the "scan" and not "spectrum" parent.
          // Some pwiz version put this term on the "spectrum" level so we also read it here.
          //TODO CV term is wrongly annotated without an xref data type -> cast to double
          spec_.setDriftTime(value.toDouble());
          spec_.setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
        }
        //scan polarity
        else if (accession == "MS:1000129") //negative scan
        {
          spec_.getInstrumentSettings().setPolarity(IonSource::NEGATIVE);
        }
        else if (accession == "MS:1000130") //positive scan
        {
          spec_.getInstrumentSettings().setPolarity(IonSource::POSITIVE);
        }
        else
          warning(LOAD, String("Unhandled cvParam '") + accession + "' in tag '" + parent_tag + "'.");
      }
      //------------------------- scanWindow ----------------------------
      else if (parent_tag == "scanWindow")
      {
        if (accession == "MS:1000501") //scan window lower limit
        {
          spec_.getInstrumentSettings().getScanWindows().back().begin = value.toDouble();
        }
        else if (accession == "MS:1000500") //scan window upper limit
        {
          spec_.getInstrumentSettings().getScanWindows().back().end = value.toDouble();
        }
        else
          warning(LOAD, String("Unhandled cvParam '") + accession + "' in tag '" + parent_tag + "'.");
      }
      //------------------------- referenceableParamGroup ----------------------------
      else if (parent_tag == "referenceableParamGroup")
      {
        SemanticValidator::CVTerm term;
        term.accession = accession;
        term.name = name;
        term.value = value;
        term.unit_accession = unit_accession;
        ref_param_[current_id_].push_back(std::move(term));
      }
      //------------------------- selectedIon ----------------------------
      else if (parent_tag == "selectedIon")
      {
        //parse only the first selected ion
        if (selected_ion_count_ > 1)
        {
          return;
        }
        if (accession == "MS:1000744") //selected ion m/z
        {
          double this_mz = value.toDouble();
          Precursor& precursor = in_spectrum_list_ ?
            spec_.getPrecursors().back() : chromatogram_.getPrecursor();
          if (this_mz != precursor.getMZ())
          {
            if (options_.getPrecursorMZSelectedIon())
            {
              // overwrite the m/z of the isolation window:
              precursor.setMetaValue("isolation window target m/z",
                                     precursor.getMZ());
              precursor.setMZ(this_mz);
            }
            else // keep precursor m/z from isolation window
            {
              precursor.setMetaValue("selected ion m/z", this_mz);
            }
          }
          // don't need to do anything if the two m/z values are the same
        }
        else if (accession == "MS:1000041") //charge state
        {
          if (in_spectrum_list_)
          {
            spec_.getPrecursors().back().setCharge(value.toInt());
          }
          else
          {
            chromatogram_.getPrecursor().setCharge(value.toInt());
          }
        }
        else if (accession == "MS:1000042") //peak intensity
        {
          if (in_spectrum_list_)
          {
            spec_.getPrecursors().back().setIntensity(value.toDouble());
          }
          else
          {
            chromatogram_.getPrecursor().setIntensity(value.toDouble());
          }
        }
        else if (accession == "MS:1000633") //possible charge state
        {
          if (in_spectrum_list_)
          {
            spec_.getPrecursors().back().getPossibleChargeStates().push_back(value.toInt());
          }
          else
          {
            chromatogram_.getPrecursor().getPossibleChargeStates().push_back(value.toInt());
          }
        }
        else if (accession == "MS:1002476" || accession == "MS:1002815" || accession == "MS:1001581") //ion mobility drift time or FAIM compensation voltage
        {
          // Drift time may be a property of the precursor (in case we are
          // acquiring a fragment ion spectrum) or of the spectrum itself.
          // According to the updated OBO, it can be a precursor or a scan
          // attribute.
          //
          // If we find here, this relates to a particular precursor. We still
          // also store it in MSSpectrum in case a client only checks there.
          // In most cases, there is a single precursor with a single drift
          // time.
          //
          // Note that only milliseconds and VSSC are valid units

          auto unit = DriftTimeUnit::MILLISECOND;
          if (accession == "MS:1002476")
          {
            unit = DriftTimeUnit::MILLISECOND;
          }
          else if (accession == "MS:1002815")
          {
            unit = DriftTimeUnit::VSSC;          
          }
          else if (accession == "MS:1001581")
          {
            unit = DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE;          
          }

          if (in_spectrum_list_)
          {
            spec_.getPrecursors().back().setDriftTime(value.toDouble());
            spec_.setDriftTime(value.toDouble());
            spec_.setDriftTimeUnit(unit);
            spec_.getPrecursors().back().setDriftTimeUnit(unit);
          }
          else
          {
            chromatogram_.getPrecursor().setDriftTime(value.toDouble());
            chromatogram_.getPrecursor().setDriftTimeUnit(unit);
          }
        }
        else
        {
          warning(LOAD, String("Unhandled cvParam '") + accession + "' in tag '" + parent_tag + "'.");
        }
      }
      //------------------------- activation ----------------------------
      else if (parent_tag == "activation")
      {
        //precursor activation attribute
        if (in_spectrum_list_)
        {
          if (accession == "MS:1000245") //charge stripping
          {
            //No member => meta data
            spec_.getPrecursors().back().setMetaValue("charge stripping", String("true"));
          }
          else if (accession == "MS:1000045") //collision energy (ev)
          {
            //No member => meta data
            spec_.getPrecursors().back().setMetaValue("collision energy", termValue);
          }
          else if (accession == "MS:1000412") //buffer gas
          {
            //No member => meta data
            spec_.getPrecursors().back().setMetaValue("buffer gas", termValue);
          }
          else if (accession == "MS:1000419") //collision gas
          {
            //No member => meta data
            spec_.getPrecursors().back().setMetaValue("collision gas", termValue);
          }
          else if (accession == "MS:1000509") //activation energy (ev)
          {
            spec_.getPrecursors().back().setActivationEnergy(value.toDouble());
          }
          else if (accession == "MS:1000138") //percent collision energy
          {
            //No member => meta data
            spec_.getPrecursors().back().setMetaValue("percent collision energy", termValue);
          }
          else if (accession == "MS:1000869") //collision gas pressure
          {
            //No member => meta data
            spec_.getPrecursors().back().setMetaValue("collision gas pressure", termValue);
          }
          //dissociation method
          else if (accession == "MS:1000044") //dissociation method
          {
            //nothing to do here
          }
          else if (accession == "MS:1000133") //collision-induced dissociation
          {
            spec_.getPrecursors().back().getActivationMethods().insert(Precursor::CID);
          }
          else if (accession == "MS:1000134") //plasma desorption
          {
            spec_.getPrecursors().back().getActivationMethods().insert(Precursor::PD);
          }
          else if (accession == "MS:1000135") //post-source decay
          {
            spec_.getPrecursors().back().getActivationMethods().insert(Precursor::PSD);
          }
          else if (accession == "MS:1000136") //surface-induced dissociation
          {
            spec_.getPrecursors().back().getActivationMethods().insert(Precursor::SID);
          }
          else if (accession == "MS:1000242") //blackbody infrared radiative dissociation
          {
            spec_.getPrecursors().back().getActivationMethods().insert(Precursor::BIRD);
          }
          else if (accession == "MS:1000250") //electron capture dissociation
          {
            spec_.getPrecursors().back().getActivationMethods().insert(Precursor::ECD);
          }
          else if (accession == "MS:1000262") //infrared multiphoton dissociation
          {
            spec_.getPrecursors().back().getActivationMethods().insert(Precursor::IMD);
          }
          else if (accession == "MS:1000282") //sustained off-resonance irradiation
          {
            spec_.getPrecursors().back().getActivationMethods().insert(Precursor::SORI);
          }
          else if (accession == "MS:1000422") //beam-type collision-induced dissociation / HCD
          {
            spec_.getPrecursors().back().getActivationMethods().insert(Precursor::HCD);
          }
          else if (accession == "MS:1002472") //trap-type collision-induced dissociation
          {
            spec_.getPrecursors().back().getActivationMethods().insert(Precursor::TRAP);
          }          
          else if (accession == "MS:1002481") //high-energy collision-induced dissociation
          {
            spec_.getPrecursors().back().getActivationMethods().insert(Precursor::HCID);
          }
          else if (accession == "MS:1000433") //low-energy collision-induced dissociation
          {
            spec_.getPrecursors().back().getActivationMethods().insert(Precursor::LCID);
          }
          else if (accession == "MS:1000435") //photodissociation
          {
            spec_.getPrecursors().back().getActivationMethods().insert(Precursor::PHD);
          }
          else if (accession == "MS:1000598") //electron transfer dissociation
          {
            spec_.getPrecursors().back().getActivationMethods().insert(Precursor::ETD);
          }
          else if (accession == "MS:1003182"  //electron transfer and collision-induced dissociation
            || accession == "MS:1002679")  // workaround: supplemental collision-induced dissociation (see https://github.com/compomics/ThermoRawFileParser/issues/182)
          {
            spec_.getPrecursors().back().getActivationMethods().insert(Precursor::ETciD);
          }
          else if (accession == "MS:1002631" //electron transfer and higher-energy collision dissociation
            || accession == "MS:1002678") // workaround: supplemental beam-type collision-induced dissociation (see https://github.com/compomics/ThermoRawFileParser/issues/182)
          {
            spec_.getPrecursors().back().getActivationMethods().insert(Precursor::EThcD);
          }
          else if (accession == "MS:1000599") //pulsed q dissociation
          {
            spec_.getPrecursors().back().getActivationMethods().insert(Precursor::PQD);
          }
          else if (accession == "MS:1001880") //in-source collision-induced dissociation
          {
            spec_.getPrecursors().back().getActivationMethods().insert(Precursor::INSOURCE);
          }
          else if (accession == "MS:1002000") //LIFT
          {
            spec_.getPrecursors().back().getActivationMethods().insert(Precursor::LIFT);
          }          
          else
            warning(LOAD, String("Unhandled cvParam '") + accession + "' in tag '" + parent_tag + "'.");
        }
        else
        {
          if (accession == "MS:1000245") //charge stripping
          {
            //No member => meta data
            chromatogram_.getPrecursor().setMetaValue("charge stripping", String("true"));
          }
          else if (accession == "MS:1000045") //collision energy (ev)
          {
            //No member => meta data
            chromatogram_.getPrecursor().setMetaValue("collision energy", termValue);
          }
          else if (accession == "MS:1000412") //buffer gas
          {
            //No member => meta data
            chromatogram_.getPrecursor().setMetaValue("buffer gas", termValue);
          }
          else if (accession == "MS:1000419") //collision gas
          {
            //No member => meta data
            chromatogram_.getPrecursor().setMetaValue("collision gas", termValue);
          }
          else if (accession == "MS:1000509") //activation energy (ev)
          {
            chromatogram_.getPrecursor().setActivationEnergy(value.toDouble());
          }
          else if (accession == "MS:1000138") //percent collision energy
          {
            //No member => meta data
            chromatogram_.getPrecursor().setMetaValue("percent collision energy", termValue);
          }
          else if (accession == "MS:1000869") //collision gas pressure
          {
            //No member => meta data
            chromatogram_.getPrecursor().setMetaValue("collision gas pressure", termValue);
          }
          //dissociation method
          else if (accession == "MS:1000044") //dissociation method
          {
            //nothing to do here
          }
          else if (accession == "MS:1000133") //collision-induced dissociation
          {
            chromatogram_.getPrecursor().getActivationMethods().insert(Precursor::CID);
          }
          else if (accession == "MS:1000134") //plasma desorption
          {
            chromatogram_.getPrecursor().getActivationMethods().insert(Precursor::PD);
          }
          else if (accession == "MS:1000135") //post-source decay
          {
            chromatogram_.getPrecursor().getActivationMethods().insert(Precursor::PSD);
          }
          else if (accession == "MS:1000136") //surface-induced dissociation
          {
            chromatogram_.getPrecursor().getActivationMethods().insert(Precursor::SID);
          }
          else if (accession == "MS:1000242") //blackbody infrared radiative dissociation
          {
            chromatogram_.getPrecursor().getActivationMethods().insert(Precursor::BIRD);
          }
          else if (accession == "MS:1000250") //electron capture dissociation
          {
            chromatogram_.getPrecursor().getActivationMethods().insert(Precursor::ECD);
          }
          else if (accession == "MS:1000262") //infrared multiphoton dissociation
          {
            chromatogram_.getPrecursor().getActivationMethods().insert(Precursor::IMD);
          }
          else if (accession == "MS:1000282") //sustained off-resonance irradiation
          {
            chromatogram_.getPrecursor().getActivationMethods().insert(Precursor::SORI);
          }
          else if (accession == "MS:1000422") //beam-type collision-induced dissociation / HCD
          {
            chromatogram_.getPrecursor().getActivationMethods().insert(Precursor::HCD);
          }
          else if (accession == "MS:1002472") //trap-type collision-induced dissociation
          {
            chromatogram_.getPrecursor().getActivationMethods().insert(Precursor::TRAP);
          }
          else if (accession == "MS:1002481") //high-energy collision-induced dissociation          
          {
            chromatogram_.getPrecursor().getActivationMethods().insert(Precursor::HCID);
          }
          else if (accession == "MS:1000433") //low-energy collision-induced dissociation
          {
            chromatogram_.getPrecursor().getActivationMethods().insert(Precursor::LCID);
          }
          else if (accession == "MS:1000435") //photodissociation
          {
            chromatogram_.getPrecursor().getActivationMethods().insert(Precursor::PHD);
          }
          else if (accession == "MS:1000598") //electron transfer dissociation
          {
            chromatogram_.getPrecursor().getActivationMethods().insert(Precursor::ETD);
          }
          else if (accession == "MS:1003182") //electron transfer and collision-induced dissociation
          {
            chromatogram_.getPrecursor().getActivationMethods().insert(Precursor::ETciD);
          }
          else if (accession == "MS:1002631") //electron transfer and higher-energy collision dissociation
          {
            chromatogram_.getPrecursor().getActivationMethods().insert(Precursor::EThcD);
          }
          else if (accession == "MS:1000599") //pulsed q dissociation
          {
            chromatogram_.getPrecursor().getActivationMethods().insert(Precursor::PQD);
          }
          else if (accession == "MS:1001880") //in-source collision-induced dissociation
          {
            chromatogram_.getPrecursor().getActivationMethods().insert(Precursor::INSOURCE);
          }
          else if (accession == "MS:1002000") //LIFT
          {
            chromatogram_.getPrecursor().getActivationMethods().insert(Precursor::LIFT);
          }          
          else
          {
            warning(LOAD, String("Unhandled cvParam '") + accession + "' in tag '" + parent_tag + "'.");
          }
        }
      }
      //------------------------- isolationWindow ----------------------------
      else if (parent_tag == "isolationWindow")
      {
        if (parent_parent_tag == "precursor")
        {
          if (accession == "MS:1000827") //isolation window target m/z
          {
            if (in_spectrum_list_)
            {
              spec_.getPrecursors().back().setMZ(value.toDouble());
            }
            else
            {
              chromatogram_.getPrecursor().setMZ(value.toDouble());
            }
          }
          else if (accession == "MS:1000828") //isolation window lower offset
          {
            if (in_spectrum_list_)
            {
              spec_.getPrecursors().back().setIsolationWindowLowerOffset(value.toDouble());
            }
            else
            {
              chromatogram_.getPrecursor().setIsolationWindowLowerOffset(value.toDouble());
            }
          }
          else if (accession == "MS:1000829") //isolation window upper offset
          {
            if (in_spectrum_list_)
            {
              spec_.getPrecursors().back().setIsolationWindowUpperOffset(value.toDouble());
            }
            else
            {
              chromatogram_.getPrecursor().setIsolationWindowUpperOffset(value.toDouble());
            }
          }
          else
            warning(LOAD, String("Unhandled cvParam '") + accession + "' in tag '" + parent_tag + "'.");
        }
        else if (parent_parent_tag == "product")
        {
          if (accession == "MS:1000827") //isolation window target m/z
          {
            if (in_spectrum_list_)
            {
              spec_.getProducts().back().setMZ(value.toDouble());
            }
            else
            {
              chromatogram_.getProduct().setMZ(value.toDouble());
            }
          }
          else if (accession == "MS:1000829") //isolation window upper offset
          {
            if (in_spectrum_list_)
            {
              spec_.getProducts().back().setIsolationWindowUpperOffset(value.toDouble());
            }
            else
            {
              chromatogram_.getProduct().setIsolationWindowUpperOffset(value.toDouble());
            }
          }
          else if (accession == "MS:1000828") //isolation window lower offset
          {
            if (in_spectrum_list_)
            {
              spec_.getProducts().back().setIsolationWindowLowerOffset(value.toDouble());
            }
            else
            {
              chromatogram_.getProduct().setIsolationWindowLowerOffset(value.toDouble());
            }
          }
          else
            warning(LOAD, String("Unhandled cvParam '") + accession + "' in tag '" + parent_tag + "'.");
        }
      }
      //------------------------- scanList ----------------------------
      else if (parent_tag == "scanList")
      {
        if (cv_.isChildOf(accession, "MS:1000570")) //method of combination as string
        {
          spec_.getAcquisitionInfo().setMethodOfCombination(cv_.getTerm(accession).name);
        }
        else
        {
          warning(LOAD, String("Unhandled cvParam '") + accession + "' in tag '" + parent_tag + "'.");
        }
      }
      //------------------------- scan ----------------------------
      else if (parent_tag == "scan")
      {
        //scan attributes
        if (accession == "MS:1000502") //dwell time
        {
          //No member => meta data
          spec_.setMetaValue("dwell time", termValue);
        }
        else if (accession == "MS:1002476" || accession == "MS:1002815" || accession == "MS:1001581") //ion mobility drift time or FAIMS compensation voltage
        {
          // Drift time may be a property of the precursor (in case we are
          // acquiring a fragment ion spectrum) or of the spectrum itself.
          // According to the updated OBO, it can be a precursor or a scan
          // attribute.
          //
          // If we find it here, it relates to the scan or spectrum itself and
          // not to a particular precursor.
          //
          // Note: this is where pwiz stores the ion mobility for a spectrum

          auto unit = DriftTimeUnit::MILLISECOND;
          if (accession == "MS:1002476")
          {
            unit = DriftTimeUnit::MILLISECOND;
          }
          else if (accession == "MS:1002815")
          {
            unit = DriftTimeUnit::VSSC;
          }
          else if (accession == "MS:1001581")
          {
            unit = DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE;
          }

          spec_.setDriftTime(value.toDouble());
          spec_.setDriftTimeUnit(unit);
        }
        else if (accession == "MS:1000011") //mass resolution
        {
          //No member => meta data
          spec_.setMetaValue("mass resolution", termValue);
        }
        else if (accession == "MS:1000015") //scan rate
        {
          //No member => meta data
          spec_.setMetaValue("scan rate", termValue);
        }
        else if (accession == "MS:1000016") //scan start time
        {
          if (unit_accession == "UO:0000031") //minutes
          {
            spec_.setRT(60.0 * value.toDouble());
          }
          else //seconds
          {
            spec_.setRT(value.toDouble());
          }
          rt_set_ = true;
          if (options_.hasRTRange())
          {
            if (!options_.getRTRange().encloses(DPosition<1>(spec_.getRT())))
            {
              skip_spectrum_ = true;
            }
            else
            { // we are within RT range
              if (load_detail_ == XMLHandler::LD_COUNTS_WITHOPTIONS)
              { //, but we only want to count
                skip_spectrum_ = true;
                ++scan_count_;
              }
            }
          }
          else if (load_detail_ == XMLHandler::LD_COUNTS_WITHOPTIONS)
          { // all RTs are valid, and the MS level of the current spectrum is in our MSLevels (otherwise we would not be here)
            skip_spectrum_ = true;
            ++scan_count_;
          }
        }
        else if (accession == "MS:1000826") //elution time
        {
          if (unit_accession == "UO:0000031") //minutes
          {
            spec_.setMetaValue("elution time (seconds)", 60.0 * value.toDouble());
          }
          else //seconds
          {
            spec_.setMetaValue("elution time (seconds)", value.toDouble());
          }
        }
        else if (accession == "MS:1000512") //filter string
        {
          //No member => meta data
          spec_.setMetaValue("filter string", termValue);
        }
        else if (accession == "MS:1000803") //analyzer scan offset
        {
          //No member => meta data
          spec_.setMetaValue("analyzer scan offset", termValue); // used in SpectraIDViewTab()
        }
        else if (accession == "MS:1000616") //preset scan configuration
        {
          //No member => meta data
          spec_.setMetaValue("preset scan configuration", termValue);
        }
        else if (accession == "MS:1000800") //mass resolving power
        {
          //No member => meta data
          spec_.setMetaValue("mass resolving power", termValue);
        }
        else if (accession == "MS:1000880") //interchannel delay
        {
          //No member => meta data
          spec_.setMetaValue("interchannel delay", termValue);
        }
        //scan direction
        else if (accession == "MS:1000092") //decreasing m/z scan
        {
          //No member => meta data
          spec_.setMetaValue("scan direction", String("decreasing"));
        }
        else if (accession == "MS:1000093") //increasing m/z scan
        {
          //No member => meta data
          spec_.setMetaValue("scan direction", String("increasing"));
        }
        //scan law
        else if (accession == "MS:1000094") //scan law: exponential
        {
          //No member => meta data
          spec_.setMetaValue("scan law", String("exponential"));
        }
        else if (accession == "MS:1000095") //scan law: linear
        {
          //No member => meta data
          spec_.setMetaValue("scan law", String("linear"));
        }
        else if (accession == "MS:1000096") //scan law: quadratic
        {
          //No member => meta data
          spec_.setMetaValue("scan law", String("quadratic"));
        }
        else if (accession == "MS:1000497") // zoom scan
        {
          spec_.getInstrumentSettings().setZoomScan(true);
        }
        else
        {
          //warning(LOAD, String("Unhandled cvParam '") + accession + "' in tag '" + parent_tag + "'."); //of course just pops up with debug flag set ...
          spec_.getAcquisitionInfo().back().setMetaValue(accession, termValue);
        }
      }
      //------------------------- contact ----------------------------
      else if (parent_tag == "contact")
      {
        if (accession == "MS:1000586") //contact name
        {
          exp_->getContacts().back().setName(value);
        }
        else if (accession == "MS:1000587") //contact address
        {
          exp_->getContacts().back().setAddress(value);
        }
        else if (accession == "MS:1000588") //contact URL
        {
          exp_->getContacts().back().setURL(value);
        }
        else if (accession == "MS:1000589") //contact email
        {
          exp_->getContacts().back().setEmail(value);
        }
        else if (accession == "MS:1000590") //contact organization
        {
          exp_->getContacts().back().setInstitution(value);
        }
        else
          warning(LOAD, String("Unhandled cvParam '") + accession + "' in tag '" + parent_tag + "'.");
      }
      //------------------------- sourceFile ----------------------------
      else if (parent_tag == "sourceFile")
      {
        if (accession == "MS:1000569") //SHA-1 checksum
        {
          source_files_[current_id_].setChecksum(value, SourceFile::SHA1);
        }
        else if (accession == "MS:1000568") //MD5 checksum
        {
          source_files_[current_id_].setChecksum(value, SourceFile::MD5);
        }
        else if (cv_.isChildOf(accession, "MS:1000560")) //source file type as string
        {
          source_files_[current_id_].setFileType(cv_.getTerm(accession).name);
        }
        else if (cv_.isChildOf(accession, "MS:1000767")) //native spectrum identifier format as string
        {
          source_files_[current_id_].setNativeIDType(cv_.getTerm(accession).name);
          source_files_[current_id_].setNativeIDTypeAccession(cv_.getTerm(accession).id);
        }
        else
          warning(LOAD, String("Unhandled cvParam '") + accession + "' in tag '" + parent_tag + "'.");
      }
      //------------------------- sample ----------------------------
      else if (parent_tag == "sample")
      {
        if (accession == "MS:1000004") //sample mass (gram)
        {
          samples_[current_id_].setMass(value.toDouble());
        }
        else if (accession == "MS:1000001") //sample number
        {
          samples_[current_id_].setNumber(value);
        }
        else if (accession == "MS:1000005") //sample volume (milliliter)
        {
          samples_[current_id_].setVolume(value.toDouble());
        }
        else if (accession == "MS:1000006") //sample concentration (gram per liter)
        {
          samples_[current_id_].setConcentration(value.toDouble());
        }
        else if (accession == "MS:1000053") //sample batch
        {
          //No member => meta data
          samples_[current_id_].setMetaValue("sample batch", termValue);
        }
        else if (accession == "MS:1000047") //emulsion
        {
          samples_[current_id_].setState(Sample::EMULSION);
        }
        else if (accession == "MS:1000048") //gas
        {
          samples_[current_id_].setState(Sample::GAS);
        }
        else if (accession == "MS:1000049") //liquid
        {
          samples_[current_id_].setState(Sample::LIQUID);
        }
        else if (accession == "MS:1000050") //solid
        {
          samples_[current_id_].setState(Sample::SOLID);
        }
        else if (accession == "MS:1000051") //solution
        {
          samples_[current_id_].setState(Sample::SOLUTION);
        }
        else if (accession == "MS:1000052") //suspension
        {
          samples_[current_id_].setState(Sample::SUSPENSION);
        }
        else if (accession.hasPrefix("PATO:")) //quality of an object
        {
          //No member => meta data
          samples_[current_id_].setMetaValue(String(name), termValue);
        }
        else if (accession.hasPrefix("GO:")) //cellular_component
        {
          //No member => meta data
          samples_[current_id_].setMetaValue("GO cellular component", String(name));
        }
        else if (accession.hasPrefix("BTO:")) //brenda source tissue ontology
        {
          //No member => meta data
          samples_[current_id_].setMetaValue("brenda source tissue", String(name));
        }
        else
          warning(LOAD, String("Unhandled cvParam '") + accession + "' in tag '" + parent_tag + "'.");
      }
      //------------------------- instrumentConfiguration ----------------------------
      else if (parent_tag == "instrumentConfiguration")
      {
        //instrument model
        if (accession == "MS:1000031")
        {
          //unknown instrument => nothing to do
        }
        else if (cv_.isChildOf(accession, "MS:1000031")) //instrument name as string
        {
          instruments_[current_id_].setName(cv_.getTerm(accession).name);
        }
        //instrument attribute
        else if (accession == "MS:1000529") //instrument serial number
        {
          //No member => meta data
          instruments_[current_id_].setMetaValue("instrument serial number", termValue);
        }
        else if (accession == "MS:1000032") //customization
        {
          instruments_[current_id_].setCustomizations(value);
        }
        else if (accession == "MS:1000236") //transmission
        {
          //No member => metadata
          instruments_[current_id_].setMetaValue("transmission", termValue);
        }
        //ion optics type
        else if (accession == "MS:1000246") //delayed extraction
        {
          instruments_[current_id_].setIonOptics(Instrument::DELAYED_EXTRACTION);
        }
        else if (accession == "MS:1000221") //magnetic deflection
        {
          instruments_[current_id_].setIonOptics(Instrument::MAGNETIC_DEFLECTION);
        }
        else if (accession == "MS:1000275") //collision quadrupole
        {
          instruments_[current_id_].setIonOptics(Instrument::COLLISION_QUADRUPOLE);
        }
        else if (accession == "MS:1000281") //selected ion flow tube
        {
          instruments_[current_id_].setIonOptics(Instrument::SELECTED_ION_FLOW_TUBE);
        }
        else if (accession == "MS:1000286") //time lag focusing
        {
          instruments_[current_id_].setIonOptics(Instrument::TIME_LAG_FOCUSING);
        }
        else if (accession == "MS:1000300") //reflectron
        {
          instruments_[current_id_].setIonOptics(Instrument::REFLECTRON);
        }
        else if (accession == "MS:1000307") //einzel lens
        {
          instruments_[current_id_].setIonOptics(Instrument::EINZEL_LENS);
        }
        else if (accession == "MS:1000309") //first stability region
        {
          instruments_[current_id_].setIonOptics(Instrument::FIRST_STABILITY_REGION);
        }
        else if (accession == "MS:1000310") //fringing field
        {
          instruments_[current_id_].setIonOptics(Instrument::FRINGING_FIELD);
        }
        else if (accession == "MS:1000311") //kinetic energy analyzer
        {
          instruments_[current_id_].setIonOptics(Instrument::KINETIC_ENERGY_ANALYZER);
        }
        else if (accession == "MS:1000320") //static field
        {
          instruments_[current_id_].setIonOptics(Instrument::STATIC_FIELD);
        }
        //ion optics attribute
        else if (accession == "MS:1000304") //accelerating voltage
        {
          //No member => metadata
          instruments_[current_id_].setMetaValue("accelerating voltage", termValue);
        }
        else if (accession == "MS:1000216") //field-free region
        {
          //No member => metadata
          instruments_[current_id_].setMetaValue("field-free region", String("true"));
        }
        else if (accession == "MS:1000308") //electric field strength
        {
          //No member => metadata
          instruments_[current_id_].setMetaValue("electric field strength", termValue);
        }
        else if (accession == "MS:1000319") //space charge effect
        {
          //No member => metadata
          instruments_[current_id_].setMetaValue("space charge effect", String("true"));
        }
        else
          warning(LOAD, String("Unhandled cvParam '") + accession + "' in tag '" + parent_tag + "'.");
      }
      else if (parent_tag == "source")
      {
        //inlet type
        if (accession == "MS:1000055") //continuous flow fast atom bombardment
        {
          instruments_[current_id_].getIonSources().back().setInletType(IonSource::CONTINUOUSFLOWFASTATOMBOMBARDMENT);
        }
        else if (accession == "MS:1000056") //direct inlet
        {
          instruments_[current_id_].getIonSources().back().setInletType(IonSource::DIRECT);
        }
        else if (accession == "MS:1000057") //electrospray inlet
        {
          instruments_[current_id_].getIonSources().back().setInletType(IonSource::ELECTROSPRAYINLET);
        }
        else if (accession == "MS:1000058") //flow injection analysis
        {
          instruments_[current_id_].getIonSources().back().setInletType(IonSource::FLOWINJECTIONANALYSIS);
        }
        else if (accession == "MS:1000059") //inductively coupled plasma
        {
          instruments_[current_id_].getIonSources().back().setInletType(IonSource::INDUCTIVELYCOUPLEDPLASMA);
        }
        else if (accession == "MS:1000060") //infusion
        {
          instruments_[current_id_].getIonSources().back().setInletType(IonSource::INFUSION);
        }
        else if (accession == "MS:1000061") //jet separator
        {
          instruments_[current_id_].getIonSources().back().setInletType(IonSource::JETSEPARATOR);
        }
        else if (accession == "MS:1000062") //membrane separator
        {
          instruments_[current_id_].getIonSources().back().setInletType(IonSource::MEMBRANESEPARATOR);
        }
        else if (accession == "MS:1000063") //moving belt
        {
          instruments_[current_id_].getIonSources().back().setInletType(IonSource::MOVINGBELT);
        }
        else if (accession == "MS:1000064") //moving wire
        {
          instruments_[current_id_].getIonSources().back().setInletType(IonSource::MOVINGWIRE);
        }
        else if (accession == "MS:1000065") //open split
        {
          instruments_[current_id_].getIonSources().back().setInletType(IonSource::OPENSPLIT);
        }
        else if (accession == "MS:1000066") //particle beam
        {
          instruments_[current_id_].getIonSources().back().setInletType(IonSource::PARTICLEBEAM);
        }
        else if (accession == "MS:1000067") //reservoir
        {
          instruments_[current_id_].getIonSources().back().setInletType(IonSource::RESERVOIR);
        }
        else if (accession == "MS:1000068") //septum
        {
          instruments_[current_id_].getIonSources().back().setInletType(IonSource::SEPTUM);
        }
        else if (accession == "MS:1000069") //thermospray inlet
        {
          instruments_[current_id_].getIonSources().back().setInletType(IonSource::THERMOSPRAYINLET);
        }
        else if (accession == "MS:1000248") //direct insertion probe
        {
          instruments_[current_id_].getIonSources().back().setInletType(IonSource::BATCH);
        }
        else if (accession == "MS:1000249") //direct liquid introduction
        {
          instruments_[current_id_].getIonSources().back().setInletType(IonSource::CHROMATOGRAPHY);
        }
        else if (accession == "MS:1000396") //membrane inlet
        {
          instruments_[current_id_].getIonSources().back().setInletType(IonSource::MEMBRANE);
        }
        else if (accession == "MS:1000485") //nanospray inlet
        {
          instruments_[current_id_].getIonSources().back().setInletType(IonSource::NANOSPRAY);
        }
        //ionization type
        else if (accession == "MS:1000071") //chemical ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::CI);
        }
        else if (accession == "MS:1000073") //electrospray ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::ESI);
        }
        else if (accession == "MS:1000074") //fast atom bombardment ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::FAB);
        }
        else if (accession == "MS:1000227") //multiphoton ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::MPI);
        }
        else if (accession == "MS:1000240") //atmospheric pressure ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::API);
        }
        else if (accession == "MS:1000247") //desorption ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::DI);
        }
        else if (accession == "MS:1000255") //flowing afterglow
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::FA);
        }
        else if (accession == "MS:1000258") //field ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::FII);
        }
        else if (accession == "MS:1000259") //glow discharge ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::GD_MS);
        }
        else if (accession == "MS:1000271") //Negative ion chemical ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::NICI);
        }
        else if (accession == "MS:1000272") //neutralization reionization mass spectrometry
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::NRMS);
        }
        else if (accession == "MS:1000273") //photoionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::PI);
        }
        else if (accession == "MS:1000274") //pyrolysis mass spectrometry
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::PYMS);
        }
        else if (accession == "MS:1000276") //resonance enhanced multiphoton ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::REMPI);
        }
        else if (accession == "MS:1000380") //adiabatic ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::AI);
        }
        else if (accession == "MS:1000381") //associative ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::ASI);
        }
        else if (accession == "MS:1000383") //autodetachment
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::AD);
        }
        else if (accession == "MS:1000384") //autoionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::AUI);
        }
        else if (accession == "MS:1000385") //charge exchange ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::CEI);
        }
        else if (accession == "MS:1000386") //chemi-ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::CHEMI);
        }
        else if (accession == "MS:1000388") //dissociative ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::DISSI);
        }
        else if (accession == "MS:1000389") //electron ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::EI);
        }
        else if (accession == "MS:1000395") //liquid secondary ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::LSI);
        }
        else if (accession == "MS:1000399") //penning ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::PEI);
        }
        else if (accession == "MS:1000400") //plasma desorption ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::PD);
        }
        else if (accession == "MS:1000402") //secondary ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::SI);
        }
        else if (accession == "MS:1000403") //soft ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::SOI);
        }
        else if (accession == "MS:1000404") //spark ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::SPI);
        }
        else if (accession == "MS:1000406") //surface ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::SUI);
        }
        else if (accession == "MS:1000407") //thermal ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::TI);
        }
        else if (accession == "MS:1000408") //vertical ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::VI);
        }
        else if (accession == "MS:1000446") //fast ion bombardment
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::FIB);
        }
        else if (accession == "MS:1000070") //atmospheric pressure chemical ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::APCI);
        }
        else if (accession == "MS:1000239") //atmospheric pressure matrix-assisted laser desorption ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::AP_MALDI);
        }
        else if (accession == "MS:1000382") //atmospheric pressure photoionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::APPI);
        }
        else if (accession == "MS:1000075") //matrix-assisted laser desorption ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::MALDI);
        }
        else if (accession == "MS:1000257") //field desorption
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::FD);
        }
        else if (accession == "MS:1000387") //desorption/ionization on silicon
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::SILI);
        }
        else if (accession == "MS:1000393") //laser desorption ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::LD);
        }
        else if (accession == "MS:1000405") //surface-assisted laser desorption ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::SALDI);
        }
        else if (accession == "MS:1000397") //microelectrospray
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::MESI);
        }
        else if (accession == "MS:1000398") //nanoelectrospray
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::NESI);
        }
        else if (accession == "MS:1000278") //surface enhanced laser desorption ionization
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::SELDI);
        }
        else if (accession == "MS:1000279") //surface enhanced neat desorption
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::SEND);
        }
        else if (accession == "MS:1000008") //ionization type (base term)
        {
          instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::IONMETHODNULL);
        }
        //source attribute
        else if (accession == "MS:1000392") //ionization efficiency
        {
          //No member => meta data
          instruments_[current_id_].getIonSources().back().setMetaValue("ionization efficiency", termValue);
        }
        else if (accession == "MS:1000486") //source potential
        {
          //No member => meta data
          instruments_[current_id_].getIonSources().back().setMetaValue("source potential", termValue);
        }
        else if (accession == "MS:1000875") // declustering potential
        {
          //No member => meta data
          instruments_[current_id_].getIonSources().back().setMetaValue("declustering potential", termValue);
        }
        else if (accession == "MS:1000876") // cone voltage
        {
          //No member => meta data
          instruments_[current_id_].getIonSources().back().setMetaValue("cone voltage", termValue);
        }
        else if (accession == "MS:1000877") // tube lens
        {
          //No member => meta data
          instruments_[current_id_].getIonSources().back().setMetaValue("tube lens", termValue);
        }
        //laser attribute
        else if (accession == "MS:1000843") // wavelength
        {
          //No member => meta data
          instruments_[current_id_].getIonSources().back().setMetaValue("wavelength", termValue);
        }
        else if (accession == "MS:1000844") // focus diameter x
        {
          //No member => meta data
          instruments_[current_id_].getIonSources().back().setMetaValue("focus diameter x", termValue);
        }
        else if (accession == "MS:1000845") // focus diameter y
        {
          //No member => meta data
          instruments_[current_id_].getIonSources().back().setMetaValue("focus diameter y", termValue);
        }
        else if (accession == "MS:1000846") // pulse energy
        {
          //No member => meta data
          instruments_[current_id_].getIonSources().back().setMetaValue("pulse energy", termValue);
        }
        else if (accession == "MS:1000847") // pulse duration
        {
          //No member => meta data
          instruments_[current_id_].getIonSources().back().setMetaValue("pulse duration", termValue);
        }
        else if (accession == "MS:1000848") // attenuation
        {
          //No member => meta data
          instruments_[current_id_].getIonSources().back().setMetaValue("attenuation", termValue);
        }
        else if (accession == "MS:1000849") // impact angle
        {
          //No member => meta data
          instruments_[current_id_].getIonSources().back().setMetaValue("impact angle", termValue);
        }
        //laser type
        else if (accession == "MS:1000850") // gas laser
        {
          //No member => meta data
          instruments_[current_id_].getIonSources().back().setMetaValue("laser type", "gas laser");
        }
        else if (accession == "MS:1000851") // solid-state laser
        {
          //No member => meta data
          instruments_[current_id_].getIonSources().back().setMetaValue("laser type", "solid-state laser");
        }
        else if (accession == "MS:1000852") // dye-laser
        {
          //No member => meta data
          instruments_[current_id_].getIonSources().back().setMetaValue("laser type", "dye-laser");
        }
        else if (accession == "MS:1000853") // free electron laser
        {
          //No member => meta data
          instruments_[current_id_].getIonSources().back().setMetaValue("laser type", "free electron laser");
        }
        //MALDI matrix application
        else if (accession == "MS:1000834") // matrix solution
        {
          //No member => meta data
          instruments_[current_id_].getIonSources().back().setMetaValue("matrix solution", termValue);
        }
        else if (accession == "MS:1000835") // matrix solution concentration
        {
          //No member => meta data
          instruments_[current_id_].getIonSources().back().setMetaValue("matrix solution concentration", termValue);
        }
        // matrix application type
        else if (accession == "MS:1000836") // dried dropplet
        {
          //No member => meta data
          instruments_[current_id_].getIonSources().back().setMetaValue("matrix application type", "dried dropplet");
        }
        else if (accession == "MS:1000837") // printed
        {
          //No member => meta data
          instruments_[current_id_].getIonSources().back().setMetaValue("matrix application type", "printed");
        }
        else if (accession == "MS:1000838") // sprayed
        {
          //No member => meta data
          instruments_[current_id_].getIonSources().back().setMetaValue("matrix application type", "sprayed");
        }
        else if (accession == "MS:1000839") //  precoated plate
        {
          //No member => meta data
          instruments_[current_id_].getIonSources().back().setMetaValue("matrix application type", " precoated plate");
        }
        else
          warning(LOAD, String("Unhandled cvParam '") + accession + "' in tag '" + parent_tag + "'.");
      }
      else if (parent_tag == "analyzer")
      {
        //mass analyzer type
        if (accession == "MS:1000079") //fourier transform ion cyclotron resonance mass spectrometer
        {
          instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::FOURIERTRANSFORM);
        }
        else if (accession == "MS:1000080") //magnetic sector
        {
          instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::SECTOR);
        }
        else if (accession == "MS:1000081") //quadrupole
        {
          instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::QUADRUPOLE);
        }
        else if (accession == "MS:1000084") //time-of-flight
        {
          instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::TOF);
        }
        else if (accession == "MS:1000254") //electrostatic energy analyzer
        {
          instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::ESA);
        }
        else if (accession == "MS:1000264") //ion trap
        {
          instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::IT);
        }
        else if (accession == "MS:1000284") //stored waveform inverse fourier transform
        {
          instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::SWIFT);
        }
        else if (accession == "MS:1000288") //cyclotron
        {
          instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::CYCLOTRON);
        }
        else if (accession == "MS:1000484") //orbitrap
        {
          instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::ORBITRAP);
        }
        else if (accession == "MS:1000078") //axial ejection linear ion trap
        {
          instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::AXIALEJECTIONLINEARIONTRAP);
        }
        else if (accession == "MS:1000082") //quadrupole ion trap
        {
          instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::PAULIONTRAP);
        }
        else if (accession == "MS:1000083") //radial ejection linear ion trap
        {
          instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::RADIALEJECTIONLINEARIONTRAP);
        }
        else if (accession == "MS:1000291") //linear ion trap
        {
          instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::LIT);
        }
        else if (accession == "MS:1000443") //mass analyzer type (base term)
        {
          instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::ANALYZERNULL);
        }
        //mass analyzer attribute
        else if (accession == "MS:1000014") //accuracy (ppm)
        {
          instruments_[current_id_].getMassAnalyzers().back().setAccuracy(value.toDouble());
        }
        else if (accession == "MS:1000022") //TOF Total Path Length (meter)
        {
          instruments_[current_id_].getMassAnalyzers().back().setTOFTotalPathLength(value.toDouble());
        }
        else if (accession == "MS:1000024") //final MS exponent
        {
          instruments_[current_id_].getMassAnalyzers().back().setFinalMSExponent(value.toInt());
        }
        else if (accession == "MS:1000025") //magnetic field strength (tesla)
        {
          instruments_[current_id_].getMassAnalyzers().back().setMagneticFieldStrength(value.toDouble());
        }
        else if (accession == "MS:1000105") //reflectron off
        {
          instruments_[current_id_].getMassAnalyzers().back().setReflectronState(MassAnalyzer::OFF);
        }
        else if (accession == "MS:1000106") //reflectron on
        {
          instruments_[current_id_].getMassAnalyzers().back().setReflectronState(MassAnalyzer::ON);
        }
        else
          warning(LOAD, String("Unhandled cvParam '") + accession + "' in tag '" + parent_tag + "'.");
      }
      else if (parent_tag == "detector")
      {
        //detector type
        if (accession == "MS:1000107") //channeltron
        {
          instruments_[current_id_].getIonDetectors().back().setType(IonDetector::CHANNELTRON);
        }
        else if (accession == "MS:1000110") //daly detector
        {
          instruments_[current_id_].getIonDetectors().back().setType(IonDetector::DALYDETECTOR);
        }
        else if (accession == "MS:1000112") //faraday cup
        {
          instruments_[current_id_].getIonDetectors().back().setType(IonDetector::FARADAYCUP);
        }
        else if (accession == "MS:1000114") //microchannel plate detector
        {
          instruments_[current_id_].getIonDetectors().back().setType(IonDetector::MICROCHANNELPLATEDETECTOR);
        }
        else if (accession == "MS:1000115") //multi-collector
        {
          instruments_[current_id_].getIonDetectors().back().setType(IonDetector::MULTICOLLECTOR);
        }
        else if (accession == "MS:1000116") //photomultiplier
        {
          instruments_[current_id_].getIonDetectors().back().setType(IonDetector::PHOTOMULTIPLIER);
        }
        else if (accession == "MS:1000253") //electron multiplier
        {
          instruments_[current_id_].getIonDetectors().back().setType(IonDetector::ELECTRONMULTIPLIER);
        }
        else if (accession == "MS:1000345") //array detector
        {
          instruments_[current_id_].getIonDetectors().back().setType(IonDetector::ARRAYDETECTOR);
        }
        else if (accession == "MS:1000346") //conversion dynode
        {
          instruments_[current_id_].getIonDetectors().back().setType(IonDetector::CONVERSIONDYNODE);
        }
        else if (accession == "MS:1000347") //dynode
        {
          instruments_[current_id_].getIonDetectors().back().setType(IonDetector::DYNODE);
        }
        else if (accession == "MS:1000348") //focal plane collector
        {
          instruments_[current_id_].getIonDetectors().back().setType(IonDetector::FOCALPLANECOLLECTOR);
        }
        else if (accession == "MS:1000349") //ion-to-photon detector
        {
          instruments_[current_id_].getIonDetectors().back().setType(IonDetector::IONTOPHOTONDETECTOR);
        }
        else if (accession == "MS:1000350") //point collector
        {
          instruments_[current_id_].getIonDetectors().back().setType(IonDetector::POINTCOLLECTOR);
        }
        else if (accession == "MS:1000351") //postacceleration detector
        {
          instruments_[current_id_].getIonDetectors().back().setType(IonDetector::POSTACCELERATIONDETECTOR);
        }
        else if (accession == "MS:1000621") //photodiode array detector
        {
          instruments_[current_id_].getIonDetectors().back().setType(IonDetector::PHOTODIODEARRAYDETECTOR);
        }
        else if (accession == "MS:1000624") //inductive detector
        {
          instruments_[current_id_].getIonDetectors().back().setType(IonDetector::INDUCTIVEDETECTOR);
        }
        else if (accession == "MS:1000108") //conversion dynode electron multiplier
        {
          instruments_[current_id_].getIonDetectors().back().setType(IonDetector::CONVERSIONDYNODEELECTRONMULTIPLIER);
        }
        else if (accession == "MS:1000109") //conversion dynode photomultiplier
        {
          instruments_[current_id_].getIonDetectors().back().setType(IonDetector::CONVERSIONDYNODEPHOTOMULTIPLIER);
        }
        else if (accession == "MS:1000111") //electron multiplier tube
        {
          instruments_[current_id_].getIonDetectors().back().setType(IonDetector::ELECTRONMULTIPLIERTUBE);
        }
        else if (accession == "MS:1000113") //focal plane array
        {
          instruments_[current_id_].getIonDetectors().back().setType(IonDetector::FOCALPLANEARRAY);
        }
        else if (accession == "MS:1000026") //detector type (base term)
        {
          instruments_[current_id_].getIonDetectors().back().setType(IonDetector::TYPENULL);
        }
        //detector attribute
        else if (accession == "MS:1000028") //detector resolution
        {
          instruments_[current_id_].getIonDetectors().back().setResolution(value.toDouble());
        }
        else if (accession == "MS:1000029") //sampling frequency
        {
          instruments_[current_id_].getIonDetectors().back().setADCSamplingFrequency(value.toDouble());
        }
        //detector acquisition mode
        else if (accession == "MS:1000117") //analog-digital converter
        {
          instruments_[current_id_].getIonDetectors().back().setAcquisitionMode(IonDetector::ADC);
        }
        else if (accession == "MS:1000118") //pulse counting
        {
          instruments_[current_id_].getIonDetectors().back().setAcquisitionMode(IonDetector::PULSECOUNTING);
        }
        else if (accession == "MS:1000119") //time-digital converter
        {
          instruments_[current_id_].getIonDetectors().back().setAcquisitionMode(IonDetector::TDC);
        }
        else if (accession == "MS:1000120") //transient recorder
        {
          instruments_[current_id_].getIonDetectors().back().setAcquisitionMode(IonDetector::TRANSIENTRECORDER);
        }
        else
          warning(LOAD, String("Unhandled cvParam '") + accession + "' in tag '" + parent_tag + "'.");
      }
      else if (parent_tag == "processingMethod")
      {
        //data processing parameter
        if (accession == "MS:1000629") //low intensity threshold (ion count)
        {
          processing_[current_id_].back()->setMetaValue("low_intensity_threshold", termValue);
        }
        else if (accession == "MS:1000631") //high intensity threshold (ion count)
        {
          processing_[current_id_].back()->setMetaValue("high_intensity_threshold", termValue);
        }
        else if (accession == "MS:1000787") //inclusive low intensity threshold
        {
          processing_[current_id_].back()->setMetaValue("inclusive_low_intensity_threshold", termValue);
        }
        else if (accession == "MS:1000788") //inclusive high intensity threshold
        {
          processing_[current_id_].back()->setMetaValue("inclusive_high_intensity_threshold", termValue);
        }
        else if (accession == "MS:1000747") //completion time
        {
          processing_[current_id_].back()->setCompletionTime(asDateTime_(value));
        }
        //file format conversion
        else if (accession == "MS:1000530") //file format conversion
        {
          processing_[current_id_].back()->getProcessingActions().insert(DataProcessing::FORMAT_CONVERSION);
        }
        else if (accession == "MS:1000544") //Conversion to mzML
        {
          processing_[current_id_].back()->getProcessingActions().insert(DataProcessing::CONVERSION_MZML);
        }
        else if (accession == "MS:1000545") //Conversion to mzXML
        {
          processing_[current_id_].back()->getProcessingActions().insert(DataProcessing::CONVERSION_MZXML);
        }
        else if (accession == "MS:1000546") //Conversion to mzData
        {
          processing_[current_id_].back()->getProcessingActions().insert(DataProcessing::CONVERSION_MZDATA);
        }
        else if (accession == "MS:1000741") //Conversion to DTA
        {
          processing_[current_id_].back()->getProcessingActions().insert(DataProcessing::CONVERSION_DTA);
        }
        //data processing action
        else if (accession == "MS:1000543") //data processing action
        {
          processing_[current_id_].back()->getProcessingActions().insert(DataProcessing::DATA_PROCESSING);
        }
        else if (accession == "MS:1000033") //deisotoping
        {
          processing_[current_id_].back()->getProcessingActions().insert(DataProcessing::DEISOTOPING);
        }
        else if (accession == "MS:1000034") //charge deconvolution
        {
          processing_[current_id_].back()->getProcessingActions().insert(DataProcessing::CHARGE_DECONVOLUTION);
        }
        else if (accession == "MS:1000035" || cv_.isChildOf(accession, "MS:1000035")) //peak picking (or child terms, we make no difference)
        {
          processing_[current_id_].back()->getProcessingActions().insert(DataProcessing::PEAK_PICKING);
        }
        else if (accession == "MS:1000592" || cv_.isChildOf(accession, "MS:1000592")) //smoothing (or child terms, we make no difference)
        {
          processing_[current_id_].back()->getProcessingActions().insert(DataProcessing::SMOOTHING);
        }
        else if (accession == "MS:1000778" || cv_.isChildOf(accession, "MS:1000778")) //charge state calculation (or child terms, we make no difference)
        {
          processing_[current_id_].back()->getProcessingActions().insert(DataProcessing::CHARGE_CALCULATION);
        }
        else if (accession == "MS:1000780" || cv_.isChildOf(accession, "MS:1000780")) //precursor recalculation (or child terms, we make no difference)
        {
          processing_[current_id_].back()->getProcessingActions().insert(DataProcessing::PRECURSOR_RECALCULATION);
        }
        else if (accession == "MS:1000593") //baseline reduction
        {
          processing_[current_id_].back()->getProcessingActions().insert(DataProcessing::BASELINE_REDUCTION);
        }
        else if (accession == "MS:1000745") //retention time alignment
        {
          processing_[current_id_].back()->getProcessingActions().insert(DataProcessing::ALIGNMENT);
        }
        else if (accession == "MS:1001484") //intensity normalization
        {
          processing_[current_id_].back()->getProcessingActions().insert(DataProcessing::NORMALIZATION);
        }
        else if (accession == "MS:1001485") //m/z calibration
        {
          processing_[current_id_].back()->getProcessingActions().insert(DataProcessing::CALIBRATION);
        }
        else if (accession == "MS:1001486" || cv_.isChildOf(accession, "MS:1001486")) //data filtering (or child terms, we make no difference)
        {
          processing_[current_id_].back()->getProcessingActions().insert(DataProcessing::FILTERING);
        }
        else
          warning(LOAD, String("Unhandled cvParam '") + accession + "' in tag '" + parent_tag + "'.");
      }
      else if (parent_tag == "fileContent")
      {
        if (cv_.isChildOf(accession, "MS:1000524")) //data file content
        {
          //ignored
          //exp_->setMetaValue(name, termValue);
        }
        else if (cv_.isChildOf(accession, "MS:1000525")) //spectrum representation
        {
          //ignored
          //exp_->setMetaValue(name, termValue);
        }
        else
          warning(LOAD, String("Unhandled cvParam '") + accession + "' in tag '" + parent_tag + "'.");
      }
      else if (parent_tag == "software")
      {
        if (cv_.isChildOf(accession, "MS:1000531")) //software as string
        {
          if (accession == "MS:1000799") //custom unreleased software tool => use value as name
          {
            software_[current_id_].setName(value);
          }
          else //use name as name
          {
            software_[current_id_].setName(name);
          }
        }
        else
        {
          warning(LOAD, String("Unhandled cvParam '") + accession + "' in tag '" + parent_tag + "'.");
        }
        //~ software_[current_id_].addCVTerm(   CVTerm (accession, value, const String &cv_identifier_ref, const String &value, const Unit &unit)   ); TODO somthing like that
      }
      else if (parent_tag == "chromatogram")
      {
        if (accession == "MS:1000810")
        {
          chromatogram_.setChromatogramType(ChromatogramSettings::MASS_CHROMATOGRAM);
        }
        else if (accession == "MS:1000235")
        {
          chromatogram_.setChromatogramType(ChromatogramSettings::TOTAL_ION_CURRENT_CHROMATOGRAM);
        }
        else if (accession == "MS:1000627")
        {
          chromatogram_.setChromatogramType(ChromatogramSettings::SELECTED_ION_CURRENT_CHROMATOGRAM);
        }
        else if (accession == "MS:1000628")
        {
          chromatogram_.setChromatogramType(ChromatogramSettings::BASEPEAK_CHROMATOGRAM);
        }
        else if (accession == "MS:1001472")
        {
          chromatogram_.setChromatogramType(ChromatogramSettings::SELECTED_ION_MONITORING_CHROMATOGRAM);
        }
        else if (accession == "MS:1001473")
        {
          chromatogram_.setChromatogramType(ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
        }
        else if (accession == "MS:1001474")
        {
          chromatogram_.setChromatogramType(ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
        }
        else if (accession == "MS:1000811")
        {
          chromatogram_.setChromatogramType(ChromatogramSettings::ELECTROMAGNETIC_RADIATION_CHROMATOGRAM);
        }
        else if (accession == "MS:1000812")
        {
          chromatogram_.setChromatogramType(ChromatogramSettings::ABSORPTION_CHROMATOGRAM);
        }
        else if (accession == "MS:1000813")
        {
          chromatogram_.setChromatogramType(ChromatogramSettings::EMISSION_CHROMATOGRAM);
        }
        else if (accession == "MS:1000809")
        {
          chromatogram_.setName(value);
        }
        else
          warning(LOAD, String("Unhandled cvParam '") + accession + "' in tag '" + parent_tag + "'.");
      }
      else if (parent_tag == "target")
      {
        //allowed but, not needed
      }
      else
        warning(LOAD, String("Unhandled cvParam '") + accession + "' in tag '" + parent_tag + "'.");
    }

    void MzMLHandler::handleUserParam_(const String& parent_parent_tag,
                                       const String& parent_tag,
                                       const String& name,
                                       const String& type,
                                       const String& value,
                                       const String& unit_accession)
    {
      // create a DataValue that contains the data in the right type
      DataValue data_value = fromXSDString(type, value);

      if (!unit_accession.empty())
      {
        if (unit_accession.hasPrefix("UO:"))
        {
          data_value.setUnit(unit_accession.suffix(unit_accession.size() - 3).toInt());
          data_value.setUnitType(DataValue::UnitType::UNIT_ONTOLOGY);
        }
        else if (unit_accession.hasPrefix("MS:"))
        {
          data_value.setUnit(unit_accession.suffix(unit_accession.size() - 3).toInt());
          data_value.setUnitType(DataValue::UnitType::MS_ONTOLOGY);
        }
        else
        {
          warning(LOAD, String("Unhandled unit '") + unit_accession + "' in tag '" + parent_tag + "'.");
        }
      }

      //find the right MetaInfoInterface
      if (parent_tag == "run")
      {
        exp_->setMetaValue(name, data_value);
      }
      else if (parent_tag == "instrumentConfiguration")
      {
        instruments_[current_id_].setMetaValue(name, data_value);
      }
      else if (parent_tag == "source")
      {
        instruments_[current_id_].getIonSources().back().setMetaValue(name, data_value);
      }
      else if (parent_tag == "analyzer")
      {
        instruments_[current_id_].getMassAnalyzers().back().setMetaValue(name, data_value);
      }
      else if (parent_tag == "detector")
      {
        instruments_[current_id_].getIonDetectors().back().setMetaValue(name, data_value);
      }
      else if (parent_tag == "sample")
      {
        samples_[current_id_].setMetaValue(name, data_value);
      }
      else if (parent_tag == "software")
      {
        software_[current_id_].setMetaValue(name, data_value);
      }
      else if (parent_tag == "contact")
      {
        exp_->getContacts().back().setMetaValue(name, data_value);
      }
      else if (parent_tag == "sourceFile")
      {
        source_files_[current_id_].setMetaValue(name, data_value);
      }
      else if (parent_tag == "binaryDataArray")
      {
        bin_data_.back().meta.setMetaValue(name, data_value);
      }
      else if (parent_tag == "spectrum")
      {
        spec_.setMetaValue(name, data_value);
      }
      else if (parent_tag == "chromatogram")
      {
        chromatogram_.setMetaValue(name, data_value);
      }
      else if (parent_tag == "scanList")
      {
        spec_.getAcquisitionInfo().setMetaValue(name, data_value);
      }
      else if (parent_tag == "scan")
      {
        spec_.getAcquisitionInfo().back().setMetaValue(name, data_value);
      }
      else if (parent_tag == "scanWindow")
      {
        spec_.getInstrumentSettings().getScanWindows().back().setMetaValue(name, data_value);
      }
      else if (parent_tag == "isolationWindow")
      {
        //We don't have this as a separate location => store it in the precursor
        if (parent_parent_tag == "precursor")
        {
          if (in_spectrum_list_)
          {
            spec_.getPrecursors().back().setMetaValue(name, data_value);
          }
          else
          {
            chromatogram_.getPrecursor().setMetaValue(name, data_value);
          }
        }
        else if (parent_parent_tag == "product")
        {
          if (in_spectrum_list_)
          {
            spec_.getProducts().back().setMetaValue(name, data_value);
          }
          else
          {
            chromatogram_.getProduct().setMetaValue(name, data_value);
          }
        }
      }
      else if (parent_tag == "selectedIon")
      {
        //parse only the first selected ion
        if (selected_ion_count_ > 1)
          return;

        //We don't have this as a separate location => store it in the precursor
        if (in_spectrum_list_)
        {
          spec_.getPrecursors().back().setMetaValue(name, data_value);
        }
        else
        {
          chromatogram_.getPrecursor().setMetaValue(name, data_value);
        }
      }
      else if (parent_tag == "activation")
      {
        //We don't have this as a separate location => store it in the precursor
        if (in_spectrum_list_)
        {
          spec_.getPrecursors().back().setMetaValue(name, data_value);
        }
        else
        {
          chromatogram_.getPrecursor().setMetaValue(name, data_value);
        }
      }
      else if (parent_tag == "processingMethod")
      {
        processing_[current_id_].back()->setMetaValue(name, data_value);
      }
      else if (parent_tag == "fileContent")
      {
        //exp_->setMetaValue(name, data_value);
      }
      else
      {
        warning(LOAD, String("Unhandled userParam '") + name + "' in tag '" + parent_tag + "'.");
      }
    }

    bool MzMLHandler::validateCV_(const ControlledVocabulary::CVTerm& c, const String& path, const Internal::MzMLValidator& validator) const
    {
      // We remember already validated path-term-combinations in cached_terms_
      // This avoids recomputing SemanticValidator::locateTerm() multiple times for the same terms and paths
      // validateCV_() is called very often for the same path-term-combinations, so we save lots of repetitive computations
      // By caching these combinations we save about 99% of the runtime of validateCV_()

      const auto it = cached_terms_.find(std::make_pair(path, c.id));
      if (it != cached_terms_.end())
      {
        return it->second;
      }

      SemanticValidator::CVTerm sc;
      sc.accession = c.id;
      sc.name = c.name;
      sc.has_unit_accession = false;
      sc.has_unit_name = false;

      bool isValid = validator.SemanticValidator::locateTerm(path, sc);
      cached_terms_[std::make_pair(path, c.id)] = isValid;
      return isValid;
    }

    String MzMLHandler::writeCV_(const ControlledVocabulary::CVTerm& c, const DataValue& metaValue) const
    {
      String cvTerm = "<cvParam cvRef=\"" + c.id.prefix(':') + "\" accession=\"" + c.id + "\" name=\"" + c.name;
      if (!metaValue.isEmpty())
      {
        cvTerm += "\" value=\"" + writeXMLEscape(metaValue.toString());
        if (metaValue.hasUnit())
        {
          //  unitAccession="UO:0000021" unitName="gram" unitCvRef="UO"
          //
          // We need to identify the correct CV term for the *unit* by
          // retrieving the identifier and looking up the term within the
          // correct ontology in our cv_ object.
          char s[8];
          snprintf(s, sizeof(s), "%07d", metaValue.getUnit()); // all CV use 7 digit identifiers padded with zeros
          String unitstring = String(s);
          if (metaValue.getUnitType() == DataValue::UnitType::UNIT_ONTOLOGY)
          {
            unitstring = "UO:" + unitstring;
          }
          else if (metaValue.getUnitType() == DataValue::UnitType::MS_ONTOLOGY)
          {
            unitstring = "MS:" + unitstring;
          }
          else
          {
            warning(LOAD, String("Unhandled unit ontology '") );
          }

          ControlledVocabulary::CVTerm unit = cv_.getTerm(unitstring);
          cvTerm += "\" unitAccession=\"" + unit.id + "\" unitName=\"" + unit.name + "\" unitCvRef=\"" + unit.id.prefix(2);
        }
      }
      cvTerm += "\"/>\n";
      return cvTerm;
    }

    void MzMLHandler::writeUserParam_(std::ostream& os, const MetaInfoInterface& meta, UInt indent, const String& path, const Internal::MzMLValidator& validator, const std::set<String>& exclude) const
    {
      std::vector<String> cvParams;
      std::vector<String> userParams;

      std::vector<String> keys;
      meta.getKeys(keys);

      for (std::vector<String>::iterator key = keys.begin(); key != keys.end(); ++key)
      {
        if (exclude.count(*key)) continue; // skip excluded entries

        // special treatment of GO and BTO terms
        // <cvParam cvRef="BTO" accession="BTO:0000199" name="cardiac muscle"/>

        if (*key == "GO cellular component" || *key == "brenda source tissue")
        {
          // the CVTerm info is in the meta value
          const ControlledVocabulary::CVTerm* c = cv_.checkAndGetTermByName(meta.getMetaValue(*key));

          if (c != nullptr)
          {
            // TODO: validate CV, we currently cannot do this as the relations in the BTO and GO are not captured by our CV impl
            cvParams.push_back(writeCV_(*c, DataValue::EMPTY));
          }
        }
        else
        {
          bool writtenAsCVTerm = false;
          const ControlledVocabulary::CVTerm* c = cv_.checkAndGetTermByName(*key);
          if (c != nullptr)
          {
            if (validateCV_(*c, path, validator))
            {
              // write CV
              cvParams.push_back(writeCV_(*c, meta.getMetaValue(*key)));
              writtenAsCVTerm = true;
            }
          }

          // if we could not write it as CVTerm we will store it at least as userParam
          if (!writtenAsCVTerm)
          {
            String userParam = "<userParam name=\"" + *key + "\" type=\"";

            const DataValue& d = meta.getMetaValue(*key);
            //determine type
            if (d.valueType() == DataValue::INT_VALUE)
            {
              userParam += "xsd:integer";
            }
            else if (d.valueType() == DataValue::DOUBLE_VALUE)
            {
              userParam += "xsd:double";
            }
            else //string or lists are converted to string
            {
              userParam += "xsd:string";
            }

            userParam += "\" value=\"" + writeXMLEscape(d.toString());

            if (d.hasUnit())
            {
              //  unitAccession="UO:0000021" unitName="gram" unitCvRef="UO"
              //
              // We need to identify the correct CV term for the *unit* by
              // retrieving the identifier and looking up the term within the
              // correct ontology in our cv_ object.
              char s[8];
              snprintf(s, sizeof(s), "%07d", d.getUnit()); // all CV use 7 digit identifiers padded with zeros
              String unitstring = String(s);
              if (d.getUnitType() == DataValue::UnitType::UNIT_ONTOLOGY)
              {
                unitstring = "UO:" + unitstring;
              }
              else if (d.getUnitType() == DataValue::UnitType::MS_ONTOLOGY)
              {
                unitstring = "MS:" + unitstring;
              }
              else
              {
                warning(LOAD, String("Unhandled unit ontology '") );
              }

              ControlledVocabulary::CVTerm unit = cv_.getTerm(unitstring);
              userParam += "\" unitAccession=\"" + unit.id + "\" unitName=\"" + unit.name + "\" unitCvRef=\"" + unit.id.prefix(2);
            }

            userParam += "\"/>\n";


            userParams.push_back(std::move(userParam));
          }
        }
      }

      // write out all the cvParams and userParams in correct order
      for (const auto& p : cvParams)
      {
        os << String(indent, '\t') << p;
      }

      for (const auto& p : userParams)
      {
        os << String(indent, '\t') << p;
      }
    }

    ControlledVocabulary::CVTerm MzMLHandler::getChildWithName_(const String& parent_accession, const String& name) const
    {
      ControlledVocabulary::CVTerm res;
      auto searcher = [&res, &name, this] (const String& child)
      {
        const ControlledVocabulary::CVTerm& current = this->cv_.getTerm(child);
        if (current.name == name)
        {
          res = current;
          return true;
        }
        return false;
      };
      cv_.iterateAllChildren(parent_accession, searcher);
      return res;
    }

    void MzMLHandler::writeSoftware_(std::ostream& os, const String& id, const Software& software, const Internal::MzMLValidator& validator)
    {
      os << "\t\t<software id=\"" << id << "\" version=\"" << software.getVersion() << "\" >\n";
      ControlledVocabulary::CVTerm so_term = getChildWithName_("MS:1000531", software.getName());
      if (so_term.id.empty())
      {
        so_term = getChildWithName_("MS:1000531", software.getName() + " software"); //act of desperation to find the right cv and keep compatible with older cv mzmls
      }
      if (so_term.id.empty())
      {
        so_term = getChildWithName_("MS:1000531", "TOPP " + software.getName()); //act of desperation to find the right cv and keep compatible with older cv mzmls
      }
      if (so_term.id == "MS:1000799")
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000799\" name=\"custom unreleased software tool\" value=\"\" />\n";
      }
      else if (!so_term.id.empty())
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"" << so_term.id << "\" name=\"" << writeXMLEscape(so_term.name) << "\" />\n";
      }
      else
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000799\" name=\"custom unreleased software tool\" value=\"" << writeXMLEscape(software.getName()) << "\" />\n";
      }
      writeUserParam_(os, software, 3, "/mzML/Software/cvParam/@accession", validator);
      os << "\t\t</software>\n";
    }

    void MzMLHandler::writeSourceFile_(std::ostream& os, const String& id, const SourceFile& source_file, const Internal::MzMLValidator& validator)
    {
      os << "\t\t\t<sourceFile id=\"" << id << "\" name=\"" << writeXMLEscape(source_file.getNameOfFile()) << "\" location=\"" << writeXMLEscape(source_file.getPathToFile()) << "\">\n";
      //checksum
      if (source_file.getChecksumType() == SourceFile::SHA1)
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000569\" name=\"SHA-1\" value=\"" << source_file.getChecksum() << "\" />\n";
      }
      else if (source_file.getChecksumType() == SourceFile::MD5)
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000568\" name=\"MD5\" value=\"" << source_file.getChecksum() << "\" />\n";
      }
      else //FORCED
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000569\" name=\"SHA-1\" value=\"\" />\n";
      }
      //file type
      ControlledVocabulary::CVTerm ft_term = getChildWithName_("MS:1000560", source_file.getFileType());
      if (ft_term.id.empty() && source_file.getFileType().hasSuffix("file"))
      {
        ft_term = getChildWithName_("MS:1000560", source_file.getFileType().chop(4) + "format");   // this is born out of desperation that sourcefile has a string interface for its filetype and not the enum, which could have been easily manipulated to the updated cv
      }
      if (!ft_term.id.empty())
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"" << ft_term.id << "\" name=\"" << ft_term.name << "\" />\n";
      }
      else //FORCED
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000564\" name=\"PSI mzData format\" />\n";
      }
      //native ID format
      ControlledVocabulary::CVTerm id_term = getChildWithName_("MS:1000767", source_file.getNativeIDType());
      if (!id_term.id.empty())
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"" << id_term.id << "\" name=\"" << id_term.name << "\" />\n";
      }
      else //FORCED
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000777\" name=\"spectrum identifier nativeID format\" />\n";
      }
      writeUserParam_(os, source_file, 4, "/mzML/fileDescription/sourceFileList/sourceFile/cvParam/@accession", validator);
      os << "\t\t\t</sourceFile>\n";
    }

    void MzMLHandler::writeDataProcessing_(std::ostream& os, const String& id, const std::vector< ConstDataProcessingPtr >& dps, const Internal::MzMLValidator& validator)
    {
      os << "\t\t<dataProcessing id=\"" << id << "\">\n";

      //FORCED
      if (dps.empty())
      {
        os << "\t\t\t<processingMethod order=\"0\" softwareRef=\"so_default\">\n";
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000544\" name=\"Conversion to mzML\" />\n";
        os << "\t\t\t\t<userParam name=\"warning\" type=\"xsd:string\" value=\"fictional processing method used to fulfill format requirements\" />\n";
        os << "\t\t\t</processingMethod>\n";
      }

      bool written = false;
      for (Size i = 0; i < dps.size(); ++i)
      {
        //data processing action
        os << "\t\t\t<processingMethod order=\"0\" softwareRef=\"so_" << id << "_pm_" << i << "\">\n";
        if (dps[i]->getProcessingActions().count(DataProcessing::DATA_PROCESSING) == 1)
        {
          os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000543\" name=\"data processing action\" />\n";
          written = true;
        }
        if (dps[i]->getProcessingActions().count(DataProcessing::CHARGE_DECONVOLUTION) == 1)
        {
          os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000034\" name=\"charge deconvolution\" />\n";
          written = true;
        }
        if (dps[i]->getProcessingActions().count(DataProcessing::DEISOTOPING) == 1)
        {
          os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000033\" name=\"deisotoping\" />\n";
          written = true;
        }
        if (dps[i]->getProcessingActions().count(DataProcessing::SMOOTHING) == 1)
        {
          os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000592\" name=\"smoothing\" />\n";
          written = true;
        }
        if (dps[i]->getProcessingActions().count(DataProcessing::CHARGE_CALCULATION) == 1)
        {
          os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000778\" name=\"charge state calculation\" />\n";
          written = true;
        }
        if (dps[i]->getProcessingActions().count(DataProcessing::PRECURSOR_RECALCULATION) == 1)
        {
          os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000780\" name=\"precursor recalculation\" />\n";
          written = true;
        }
        if (dps[i]->getProcessingActions().count(DataProcessing::BASELINE_REDUCTION) == 1)
        {
          os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000593\" name=\"baseline reduction\" />\n";
          written = true;
        }
        if (dps[i]->getProcessingActions().count(DataProcessing::PEAK_PICKING) == 1)
        {
          os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000035\" name=\"peak picking\" />\n";
          written = true;
        }
        if (dps[i]->getProcessingActions().count(DataProcessing::ALIGNMENT) == 1)
        {
          os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000745\" name=\"retention time alignment\" />\n";
          written = true;
        }
        if (dps[i]->getProcessingActions().count(DataProcessing::CALIBRATION) == 1)
        {
          os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1001485\" name=\"m/z calibration\" />\n";
          written = true;
        }
        if (dps[i]->getProcessingActions().count(DataProcessing::NORMALIZATION) == 1)
        {
          os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1001484\" name=\"intensity normalization\" />\n";
          written = true;
        }
        if (dps[i]->getProcessingActions().count(DataProcessing::FILTERING) == 1)
        {
          os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1001486\" name=\"data filtering\" />\n";
          written = true;
        }
        //file format conversion
        if (dps[i]->getProcessingActions().count(DataProcessing::FORMAT_CONVERSION) == 1)
        {
          os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000530\" name=\"file format conversion\" />\n";
          written = true;
        }
        if (dps[i]->getProcessingActions().count(DataProcessing::CONVERSION_MZDATA) == 1)
        {
          os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000546\" name=\"Conversion to mzData\" />\n";
          written = true;
        }
        if (dps[i]->getProcessingActions().count(DataProcessing::CONVERSION_MZML) == 1)
        {
          os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000544\" name=\"Conversion to mzML\" />\n";
          written = true;
        }
        if (dps[i]->getProcessingActions().count(DataProcessing::CONVERSION_MZXML) == 1)
        {
          os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000545\" name=\"Conversion to mzXML\" />\n";
          written = true;
        }
        if (dps[i]->getProcessingActions().count(DataProcessing::CONVERSION_DTA) == 1)
        {
          os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000741\" name=\"Conversion to dta\" />\n";
          written = true;
        }
        if (!written)
        {
          os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000543\" name=\"data processing action\" />\n";
        }

        //data processing attribute
        if (dps[i]->getCompletionTime().isValid())
        {
          os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000747\" name=\"completion time\" value=\"" << dps[i]->getCompletionTime().toString("yyyy-MM-dd+hh:mm") << "\" />\n";
        }

        writeUserParam_(os, *(dps[i].get()), 4, "/mzML/dataProcessingList/dataProcessing/processingMethod/cvParam/@accession", validator);
        os << "\t\t\t</processingMethod>\n";
      }

      os << "\t\t</dataProcessing>\n";
    }

    void MzMLHandler::writePrecursor_(std::ostream& os, const Precursor& precursor, const Internal::MzMLValidator& validator)
    {
      // optional attributes
      String external_spectrum_id =
          precursor.metaValueExists("external_spectrum_id") ?
          " externalSpectrumID=\"" + precursor.getMetaValue("external_spectrum_id").toString() + "\"" :
          "";
      String spectrum_ref =
          precursor.metaValueExists("spectrum_ref") ?
          " spectrumRef=\"" + precursor.getMetaValue("spectrum_ref").toString() + "\"":
          "";

      os << "\t\t\t\t\t<precursor" + external_spectrum_id + spectrum_ref + ">\n";
      //--------------------------------------------------------------------------------------------
      //isolation window (optional)
      //--------------------------------------------------------------------------------------------

      // precursor m/z may come from "selected ion":
      double mz = precursor.getMetaValue("isolation window target m/z",
                                         precursor.getMZ());
      // Note that TPP parsers break when the isolation window is written out
      // in mzML files and the precursorMZ gets set to zero.
      if (mz > 0.0 && !options_.getForceTPPCompatability())
      {
        os << "\t\t\t\t\t\t<isolationWindow>\n";
        os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000827\" name=\"isolation window target m/z\" value=\"" << mz << "\" unitAccession=\"MS:1000040\" unitName=\"m/z\" unitCvRef=\"MS\" />\n";
        if (precursor.getIsolationWindowLowerOffset() > 0.0)
        {
          os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000828\" name=\"isolation window lower offset\" value=\"" << precursor.getIsolationWindowLowerOffset() << "\" unitAccession=\"MS:1000040\" unitName=\"m/z\" unitCvRef=\"MS\" />\n";
        }
        if (precursor.getIsolationWindowUpperOffset() > 0.0)
        {
          os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000829\" name=\"isolation window upper offset\" value=\"" << precursor.getIsolationWindowUpperOffset() << "\" unitAccession=\"MS:1000040\" unitName=\"m/z\" unitCvRef=\"MS\" />\n";
        }
        os << "\t\t\t\t\t\t</isolationWindow>\n";
      }
      //userParam: no extra object for it => no user parameters

      //--------------------------------------------------------------------------------------------
      //selected ion list (optional)
      //--------------------------------------------------------------------------------------------
      //

      if (options_.getForceTPPCompatability() ||
          precursor.getCharge() != 0 ||
          precursor.getIntensity() > 0.0 ||
          precursor.getDriftTime() >= 0.0 ||
          precursor.getDriftTimeUnit() == DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE ||
          !precursor.getPossibleChargeStates().empty() ||
          precursor.getMZ() > 0.0)
      {
        // precursor m/z may come from "isolation window":
        mz = precursor.getMetaValue("selected ion m/z",
                                    precursor.getMZ());
        os << "\t\t\t\t\t\t<selectedIonList count=\"1\">\n";
        os << "\t\t\t\t\t\t\t<selectedIon>\n";
        os << "\t\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000744\" name=\"selected ion m/z\" value=\"" << mz << "\" unitAccession=\"MS:1000040\" unitName=\"m/z\" unitCvRef=\"MS\" />\n";
        if (options_.getForceTPPCompatability() || precursor.getCharge() != 0)
        {
          os << "\t\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000041\" name=\"charge state\" value=\"" << precursor.getCharge() << "\" />\n";
        }
        if ( precursor.getIntensity() > 0.0)
        {
          os << "\t\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000042\" name=\"peak intensity\" value=\"" << precursor.getIntensity() << "\" unitAccession=\"MS:1000132\" unitName=\"percent of base peak\" unitCvRef=\"MS\" />\n";
        }
        for (Size j = 0; j < precursor.getPossibleChargeStates().size(); ++j)
        {
          os << "\t\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000633\" name=\"possible charge state\" value=\"" << precursor.getPossibleChargeStates()[j] << "\" />\n";
        }

        if (precursor.getDriftTime() != IMTypes::DRIFTTIME_NOT_SET)
        {
          switch (precursor.getDriftTimeUnit())
          {
            default:
              // assume milliseconds, but warn
              warning(STORE, String("Precursor drift time unit not set, assume milliseconds"));
              [[fallthrough]];
            case DriftTimeUnit::MILLISECOND:
              os << "\t\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1002476\" name=\"ion mobility drift time\" value=\"" << precursor.getDriftTime()
                  << "\" unitAccession=\"UO:0000028\" unitName=\"millisecond\" unitCvRef=\"UO\" />\n";
              break;
            case DriftTimeUnit::VSSC:
              os << "\t\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1002815\" name=\"inverse reduced ion mobility\" value=\"" << precursor.getDriftTime()
                  << "\" unitAccession=\"MS:1002814\" unitName=\"volt-second per square centimeter\" unitCvRef=\"MS\" />\n";
              break;
          }
        }
        //userParam: no extra object for it => no user parameters
        os << "\t\t\t\t\t\t\t</selectedIon>\n";
        os << "\t\t\t\t\t\t</selectedIonList>\n";
      }

      //--------------------------------------------------------------------------------------------
      //activation (mandatory)
      //--------------------------------------------------------------------------------------------
      os << "\t\t\t\t\t\t<activation>\n";
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-equal"
      if (precursor.getActivationEnergy() != 0)
#pragma clang diagnostic pop
      {
        os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000509\" name=\"activation energy\" value=\"" << precursor.getActivationEnergy() << "\" unitAccession=\"UO:0000266\" unitName=\"electronvolt\" unitCvRef=\"UO\" />\n";
      }
      if (precursor.getActivationMethods().count(Precursor::CID) != 0)
      {
        os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000133\" name=\"collision-induced dissociation\" />\n";
      }
      if (precursor.getActivationMethods().count(Precursor::PD) != 0)
      {
        os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000134\" name=\"plasma desorption\" />\n";
      }
      if (precursor.getActivationMethods().count(Precursor::PSD) != 0)
      {
        os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000135\" name=\"post-source decay\" />\n";
      }
      if (precursor.getActivationMethods().count(Precursor::SID) != 0)
      {
        os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000136\" name=\"surface-induced dissociation\" />\n";
      }
      if (precursor.getActivationMethods().count(Precursor::BIRD) != 0)
      {
        os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000242\" name=\"blackbody infrared radiative dissociation\" />\n";
      }
      if (precursor.getActivationMethods().count(Precursor::ECD) != 0)
      {
        os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000250\" name=\"electron capture dissociation\" />\n";
      }
      if (precursor.getActivationMethods().count(Precursor::IMD) != 0)
      {
        os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000262\" name=\"infrared multiphoton dissociation\" />\n";
      }
      if (precursor.getActivationMethods().count(Precursor::SORI) != 0)
      {
        os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000282\" name=\"sustained off-resonance irradiation\" />\n";
      }
      if (precursor.getActivationMethods().count(Precursor::HCID) != 0)
      {
        os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1002481\" name=\"high-energy collision-induced dissociation\" />\n";
      }
      if (precursor.getActivationMethods().count(Precursor::HCD) != 0)
      {
        os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000422\" name=\"beam-type collision-induced dissociation\" />\n";
      }
      if (precursor.getActivationMethods().count(Precursor::TRAP) != 0)
      {
        os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1002472\" name=\"trap-type collision-induced dissociation\" />\n";
      }
      if (precursor.getActivationMethods().count(Precursor::LCID) != 0)
      {
        os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000433\" name=\"low-energy collision-induced dissociation\" />\n";
      }
      if (precursor.getActivationMethods().count(Precursor::PHD) != 0)
      {
        os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000435\" name=\"photodissociation\" />\n";
      }
      if (precursor.getActivationMethods().count(Precursor::ETD) != 0)
      {
        os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000598\" name=\"electron transfer dissociation\" />\n";
      }
      if (precursor.getActivationMethods().count(Precursor::ETciD) != 0)
      {
        os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1003182\" name=\"electron transfer and collision-induced dissociation\" />\n";
      }
      if (precursor.getActivationMethods().count(Precursor::EThcD) != 0)
      {
        os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1002631\" name=\"electron transfer and higher-energy collision dissociation\" />\n";
      }
      if (precursor.getActivationMethods().count(Precursor::PQD) != 0)
      {
        os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000599\" name=\"pulsed q dissociation\" />\n";
      }
      if (precursor.getActivationMethods().count(Precursor::INSOURCE) != 0)
      {
        os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1001880\" name=\"in-source collision-induced dissociation\" />\n";
      }
      if (precursor.getActivationMethods().count(Precursor::LIFT) != 0)
      {
        os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1002000\" name=\"LIFT\" />\n";
      }      
      if (precursor.getActivationMethods().empty())
      {
        os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000044\" name=\"dissociation method\" />\n";
      }
      // as "precursor" has no own user param its userParam is stored here;
      // don't write out parameters that are used internally to distinguish
      // between precursor m/z values from different sources:
      writeUserParam_(os, precursor, 7, "/mzML/run/spectrumList/spectrum/precursorList/precursor/activation/cvParam/@accession", validator, {"isolation window target m/z", "selected ion m/z", "external_spectrum_id", "spectrum_ref"});
      os << "\t\t\t\t\t\t</activation>\n";
      os << "\t\t\t\t\t</precursor>\n";

    }

    void MzMLHandler::writeProduct_(std::ostream& os, const Product& product, const Internal::MzMLValidator& validator)
    {
      os << "\t\t\t\t\t<product>\n";
      os << "\t\t\t\t\t\t<isolationWindow>\n";
      os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000827\" name=\"isolation window target m/z\" value=\"" << product.getMZ() << "\" unitAccession=\"MS:1000040\" unitName=\"m/z\" unitCvRef=\"MS\" />\n";
      if ( product.getIsolationWindowLowerOffset() > 0.0)
      {
          os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000828\" name=\"isolation window lower offset\" value=\"" << product.getIsolationWindowLowerOffset() << "\" unitAccession=\"MS:1000040\" unitName=\"m/z\" unitCvRef=\"MS\" />\n";
      }
      if ( product.getIsolationWindowUpperOffset() > 0.0)
      {
          os << "\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000829\" name=\"isolation window upper offset\" value=\"" << product.getIsolationWindowUpperOffset() << "\" unitAccession=\"MS:1000040\" unitName=\"m/z\" unitCvRef=\"MS\" />\n";
      }
      writeUserParam_(os, product, 7, "/mzML/run/spectrumList/spectrum/productList/product/isolationWindow/cvParam/@accession", validator);
      os << "\t\t\t\t\t\t</isolationWindow>\n";
      os << "\t\t\t\t\t</product>\n";
    }

    void MzMLHandler::writeTo(std::ostream& os)
    {
      const MapType& exp = *(cexp_);
      logger_.startProgress(0, exp.size() + exp.getChromatograms().size(), "storing mzML file");
      int progress = 0;
      UInt stored_spectra = 0;
      UInt stored_chromatograms = 0;
      Internal::MzMLValidator validator(mapping_, cv_);

      std::vector<std::vector< ConstDataProcessingPtr > > dps;
      //--------------------------------------------------------------------------------------------
      //header
      //--------------------------------------------------------------------------------------------
      writeHeader_(os, exp, dps, validator);

      //--------------------------------------------------------------------------------------------
      // spectra
      //--------------------------------------------------------------------------------------------
      if (!exp.empty())
      {
        // INFO : do not try to be smart and skip empty spectra or
        // chromatograms. There can be very good reasons for this (e.g. if the
        // meta information needs to be stored here but the actual data is
        // stored somewhere else).
        os << "\t\t<spectrumList count=\"" << exp.size() << "\" defaultDataProcessingRef=\"dp_sp_0\">\n";

        // check native ids
        bool renew_native_ids = false;
        for (Size s_idx = 0; s_idx < exp.size(); ++s_idx)
        {
          if (!exp[s_idx].getNativeID().has('='))
          {
            renew_native_ids = true;
            break;
          }
        }

        // issue warning if something is wrong
        if (renew_native_ids)
        {
          warning(STORE, String("Invalid native IDs detected. Using spectrum identifier nativeID format (spectrum=xsd:nonNegativeInteger) for all spectra."));
        }

        // write actual data
        for (Size s_idx = 0; s_idx < exp.size(); ++s_idx)
        {
          logger_.setProgress(progress++);
          const SpectrumType& spec = exp[s_idx];
          writeSpectrum_(os, spec, s_idx, validator, renew_native_ids, dps);
          ++stored_spectra;
        }
        os << "\t\t</spectrumList>\n";
      }

      //--------------------------------------------------------------------------------------------
      // chromatograms
      //--------------------------------------------------------------------------------------------
      if (!exp.getChromatograms().empty())
      {
        // INFO : do not try to be smart and skip empty spectra or
        // chromatograms. There can be very good reasons for this (e.g. if the
        // meta information needs to be stored here but the actual data is
        // stored somewhere else).
        os << "\t\t<chromatogramList count=\"" << exp.getChromatograms().size() << "\" defaultDataProcessingRef=\"dp_sp_0\">\n";
        for (Size c_idx = 0; c_idx != exp.getChromatograms().size(); ++c_idx)
        {
          logger_.setProgress(progress++);
          const ChromatogramType& chromatogram = exp.getChromatograms()[c_idx];
          writeChromatogram_(os, chromatogram, c_idx, validator);
          ++stored_chromatograms;
        }
        os << "\t\t</chromatogramList>" << "\n";
      }

      MzMLHandlerHelper::writeFooter_(os, options_, spectra_offsets_, chromatograms_offsets_);

      OPENMS_LOG_INFO << stored_spectra << " spectra and " << stored_chromatograms << " chromatograms stored.\n";

      logger_.endProgress(os.tellp());
    }

    void MzMLHandler::writeHeader_(std::ostream& os,
                                   const MapType& exp,
                                   std::vector<std::vector< ConstDataProcessingPtr > >& dps,
                                   const Internal::MzMLValidator& validator)
    {
      os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";

      if (options_.getWriteIndex())
      {
        os << "<indexedmzML xmlns=\"http://psi.hupo.org/ms/mzml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0_idx.xsd\">\n";
      }
      os << R"(<mzML xmlns="http://psi.hupo.org/ms/mzml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd" accession=")" << writeXMLEscape(exp.getIdentifier()) << "\" version=\"" << version_ << "\">\n";
      //--------------------------------------------------------------------------------------------
      // CV list
      //--------------------------------------------------------------------------------------------
      os << "\t<cvList count=\"5\">\n"
         << "\t\t<cv id=\"MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Ontology\" URI=\"http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo\"/>\n"
         << "\t\t<cv id=\"UO\" fullName=\"Unit Ontology\" URI=\"http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo\"/>\n"
         << "\t\t<cv id=\"BTO\" fullName=\"BrendaTissue545\" version=\"unknown\" URI=\"http://www.brenda-enzymes.info/ontology/tissue/tree/update/update_files/BrendaTissueOBO\"/>\n"
         << "\t\t<cv id=\"GO\" fullName=\"Gene Ontology - Slim Versions\" version=\"unknown\" URI=\"http://www.geneontology.org/GO_slims/goslim_goa.obo\"/>\n"
         << "\t\t<cv id=\"PATO\" fullName=\"Quality ontology\" version=\"unknown\" URI=\"http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/quality.obo\"/>\n"
         << "\t</cvList>\n";
      //--------------------------------------------------------------------------------------------
      // file content
      //--------------------------------------------------------------------------------------------
      os << "\t<fileDescription>\n";
      os << "\t\t<fileContent>\n";
      std::map<InstrumentSettings::ScanMode, UInt> file_content;
      for (Size i = 0; i < exp.size(); ++i)
      {
        ++file_content[exp[i].getInstrumentSettings().getScanMode()];
      }
      if (file_content.find(InstrumentSettings::MASSSPECTRUM) != file_content.end())
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000294\" name=\"mass spectrum\" />\n";
      }
      if (file_content.find(InstrumentSettings::MS1SPECTRUM) != file_content.end())
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000579\" name=\"MS1 spectrum\" />\n";
      }
      if (file_content.find(InstrumentSettings::MSNSPECTRUM) != file_content.end())
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" />\n";
      }
      if (file_content.find(InstrumentSettings::SIM) != file_content.end())
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000582\" name=\"SIM spectrum\" />\n";
      }
      if (file_content.find(InstrumentSettings::SRM) != file_content.end())
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000583\" name=\"SRM spectrum\" />\n";
      }
      if (file_content.find(InstrumentSettings::CRM) != file_content.end())
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000581\" name=\"CRM spectrum\" />\n";
      }
      if (file_content.find(InstrumentSettings::PRECURSOR) != file_content.end())
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000341\" name=\"precursor ion spectrum\" />\n";
      }
      if (file_content.find(InstrumentSettings::CNG) != file_content.end())
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000325\" name=\"constant neutral gain spectrum\" />\n";
      }
      if (file_content.find(InstrumentSettings::CNL) != file_content.end())
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000326\" name=\"constant neutral loss spectrum\" />\n";
      }
      if (file_content.find(InstrumentSettings::EMR) != file_content.end())
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000804\" name=\"electromagnetic radiation spectrum\" />\n";
      }
      if (file_content.find(InstrumentSettings::EMISSION) != file_content.end())
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000805\" name=\"emission spectrum\" />\n";
      }
      if (file_content.find(InstrumentSettings::ABSORPTION) != file_content.end())
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000806\" name=\"absorption spectrum\" />\n";
      }
      if (file_content.find(InstrumentSettings::EMC) != file_content.end())
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000789\" name=\"enhanced multiply charged spectrum\" />\n";
      }
      if (file_content.find(InstrumentSettings::TDF) != file_content.end())
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000789\" name=\"time-delayed fragmentation spectrum\" />\n";
      }
      if (file_content.find(InstrumentSettings::UNKNOWN) != file_content.end() || file_content.empty())
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000294\" name=\"mass spectrum\" />\n";
      }
      // writeUserParam_(os, exp, 3, "/mzML/fileDescription/fileContent/cvParam/@accession", validator);
      os << "\t\t</fileContent>\n";

      //--------------------------------------------------------------------------------------------
      // source file list
      //--------------------------------------------------------------------------------------------
      //find out how many spectra source files need to be written
      UInt sf_sp_count = 0;
      for (Size i = 0; i < exp.size(); ++i)
      {
        if (exp[i].getSourceFile() != SourceFile())
        {
          ++sf_sp_count;
        }
      }
      if (!exp.getSourceFiles().empty() || sf_sp_count > 0)
      {
        os << "\t\t<sourceFileList count=\"" << exp.getSourceFiles().size() + sf_sp_count << "\">\n";

        //write source file of run
        for (Size i = 0; i < exp.getSourceFiles().size(); ++i)
        {
          writeSourceFile_(os, String("sf_ru_") + String(i), exp.getSourceFiles()[i], validator);
        }

        // write source files of spectra
        if (sf_sp_count > 0)
        {
          const SourceFile sf_default;
          for (Size i = 0; i < exp.size(); ++i)
          {
            if (exp[i].getSourceFile() != sf_default)
            {
              writeSourceFile_(os, String("sf_sp_") + i, exp[i].getSourceFile(), validator);
            }
          }
        }

        os << "\t\t</sourceFileList>\n";
      }

      //--------------------------------------------------------------------------------------------
      // contacts
      //--------------------------------------------------------------------------------------------
      for (Size i = 0; i < exp.getContacts().size(); ++i)
      {
        const ContactPerson& cp = exp.getContacts()[i];
        os << "\t\t<contact>\n";
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000586\" name=\"contact name\" value=\"" << writeXMLEscape(cp.getLastName()) << ", " << writeXMLEscape(cp.getFirstName()) << "\" />\n";
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000590\" name=\"contact affiliation\" value=\"" << writeXMLEscape(cp.getInstitution()) << "\" />\n";

        if (!cp.getAddress().empty())
        {
          os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000587\" name=\"contact address\" value=\"" << writeXMLEscape(cp.getAddress()) << "\" />\n";
        }
        if (!cp.getURL().empty())
        {
          os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000588\" name=\"contact URL\" value=\"" << writeXMLEscape(cp.getURL()) << "\" />\n";
        }
        if (!cp.getEmail().empty())
        {
          os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000589\" name=\"contact email\" value=\"" << writeXMLEscape(cp.getEmail()) << "\" />\n";
        }
        if (!cp.getContactInfo().empty())
        {
          os << "\t\t\t<userParam name=\"contact_info\" type=\"xsd:string\" value=\"" << writeXMLEscape(cp.getContactInfo()) << "\" />\n";
        }
        writeUserParam_(os, cp, 3, "/mzML/fileDescription/contact/cvParam/@accession", validator);
        os << "\t\t</contact>\n";
      }
      os << "\t</fileDescription>\n";

      //--------------------------------------------------------------------------------------------
      // sample
      //--------------------------------------------------------------------------------------------
      const Sample& sa = exp.getSample();
      os << "\t<sampleList count=\"1\">\n";
      os << "\t\t<sample id=\"sa_0\" name=\"" << writeXMLEscape(sa.getName()) << "\">\n";
      if (!sa.getNumber().empty())
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000001\" name=\"sample number\" value=\"" << writeXMLEscape(sa.getNumber()) << "\" />\n";
      }
      os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000004\" name=\"sample mass\" value=\"" << sa.getMass() << "\" unitAccession=\"UO:0000021\" unitName=\"gram\" unitCvRef=\"UO\" />\n";
      os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000005\" name=\"sample volume\" value=\"" << sa.getVolume() << "\" unitAccession=\"UO:0000098\" unitName=\"milliliter\" unitCvRef=\"UO\" />\n";
      os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000006\" name=\"sample concentration\" value=\"" << sa.getConcentration() << "\" unitAccession=\"UO:0000175\" unitName=\"gram per liter\" unitCvRef=\"UO\" />\n";
      if (sa.getState() == Sample::EMULSION)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000047\" name=\"emulsion\" />\n";
      }
      else if (sa.getState() == Sample::GAS)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000048\" name=\"gas\" />\n";
      }
      else if (sa.getState() == Sample::LIQUID)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000049\" name=\"liquid\" />\n";
      }
      else if (sa.getState() == Sample::SOLID)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000050\" name=\"solid\" />\n";
      }
      else if (sa.getState() == Sample::SOLUTION)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000051\" name=\"solution\" />\n";
      }
      else if (sa.getState() == Sample::SUSPENSION)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000052\" name=\"suspension\" />\n";
      }
      if (!sa.getComment().empty())
      {
        os << "\t\t\t<userParam name=\"comment\" type=\"xsd:string\" value=\"" << writeXMLEscape(sa.getComment()) << "\" />\n";
      }
      writeUserParam_(os, sa, 3, "/mzML/sampleList/sample/cvParam/@accession", validator);
      os << "\t\t</sample>\n";
      os << "\t</sampleList>\n";

      //--------------------------------------------------------------------------------------------
      // Software
      //--------------------------------------------------------------------------------------------

      // instrument software and fallback software is always written (see below)
      Size num_software(2);

      // Create a list of all different data processings: check if the
      // DataProcessing of the current spectra/chromatogram is already present
      // and if not, append it to the dps vector
      for (Size s = 0; s < exp.size(); ++s)
      {
        bool already_present = false;
        for (Size j = 0; j < dps.size(); j++)
        {
          already_present = OpenMS::Helpers::cmpPtrContainer(
              exp[s].getDataProcessing(), dps[j]);
          if (already_present) break;
        }
        if (!already_present)
        {
          dps.push_back(exp[s].getDataProcessing());
          num_software += exp[s].getDataProcessing().size();
        }
      }
      for (Size s = 0; s < exp.getChromatograms().size(); ++s)
      {
        bool already_present = false;
        for (Size j = 0; j < dps.size(); j++)
        {
          already_present = OpenMS::Helpers::cmpPtrContainer(
              exp.getChromatograms()[s].getDataProcessing(), dps[j]);
          if (already_present) break;
        }
        if (!already_present)
        {
          dps.push_back(exp.getChromatograms()[s].getDataProcessing());
          num_software += exp.getChromatograms()[s].getDataProcessing().size();
        }
      }

      // count binary data array software
      Size num_bi_software(0);

      for (Size s = 0; s < exp.size(); ++s)
      {
        for (Size m = 0; m < exp[s].getFloatDataArrays().size(); ++m)
        {
          for (Size i = 0; i < exp[s].getFloatDataArrays()[m].getDataProcessing().size(); ++i)
          {
            ++num_bi_software;
          }
        }
      }

      os << "\t<softwareList count=\"" << num_software + num_bi_software << "\">\n";

      // write instrument software
      writeSoftware_(os, "so_in_0", exp.getInstrument().getSoftware(), validator);

      // write fallback software
      writeSoftware_(os, "so_default", Software(), validator);

      // write the software of the dps
      for (Size s1 = 0; s1 != dps.size(); ++s1)
      {
        for (Size s2 = 0; s2 != dps[s1].size(); ++s2)
        {
          writeSoftware_(os, String("so_dp_sp_") + s1 + "_pm_" + s2, dps[s1][s2]->getSoftware(), validator);
        }
      }

      //write data processing (for each binary data array)
      for (Size s = 0; s < exp.size(); ++s)
      {
        for (Size m = 0; m < exp[s].getFloatDataArrays().size(); ++m)
        {
          for (Size i = 0; i < exp[s].getFloatDataArrays()[m].getDataProcessing().size(); ++i)
          {
            writeSoftware_(os, String("so_dp_sp_") + s + "_bi_" + m + "_pm_" + i,
                exp[s].getFloatDataArrays()[m].getDataProcessing()[i]->getSoftware(), validator);
          }
        }
      }
      os << "\t</softwareList>\n";

      //--------------------------------------------------------------------------------------------
      // instrument configuration (enclosing ion source, mass analyzer and detector)
      //--------------------------------------------------------------------------------------------
      const Instrument& in = exp.getInstrument();
      os << "\t<instrumentConfigurationList count=\"1\">\n";
      os << "\t\t<instrumentConfiguration id=\"ic_0\">\n";
      ControlledVocabulary::CVTerm in_term = getChildWithName_("MS:1000031", in.getName());
      if (!in_term.id.empty())
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"" << in_term.id << "\" name=\"" << writeXMLEscape(in_term.name) << "\" />\n";
      }
      else
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000031\" name=\"instrument model\" />\n";
      }

      if (!in.getCustomizations().empty())
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000032\" name=\"customization\" value=\"" << writeXMLEscape(in.getCustomizations()) << "\" />\n";
      }

      //ion optics
      if (in.getIonOptics() == Instrument::MAGNETIC_DEFLECTION)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000221\" name=\"magnetic deflection\" />\n";
      }
      else if (in.getIonOptics() == Instrument::DELAYED_EXTRACTION)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000246\" name=\"delayed extraction\" />\n";
      }
      else if (in.getIonOptics() == Instrument::COLLISION_QUADRUPOLE)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000275\" name=\"collision quadrupole\" />\n";
      }
      else if (in.getIonOptics() == Instrument::SELECTED_ION_FLOW_TUBE)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000281\" name=\"selected ion flow tube\" />\n";
      }
      else if (in.getIonOptics() == Instrument::TIME_LAG_FOCUSING)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000286\" name=\"time lag focusing\" />\n";
      }
      else if (in.getIonOptics() == Instrument::REFLECTRON)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000300\" name=\"reflectron\" />\n";
      }
      else if (in.getIonOptics() == Instrument::EINZEL_LENS)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000307\" name=\"einzel lens\" />\n";
      }
      else if (in.getIonOptics() == Instrument::FIRST_STABILITY_REGION)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000309\" name=\"first stability region\" />\n";
      }
      else if (in.getIonOptics() == Instrument::FRINGING_FIELD)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000310\" name=\"fringing field\" />\n";
      }
      else if (in.getIonOptics() == Instrument::KINETIC_ENERGY_ANALYZER)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000311\" name=\"kinetic energy analyzer\" />\n";
      }
      else if (in.getIonOptics() == Instrument::STATIC_FIELD)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000320\" name=\"static field\" />\n";
      }

      writeUserParam_(os, in, 3, "/mzML/instrumentConfigurationList/instrumentConfiguration/cvParam/@accession", validator);
      Size component_count = in.getIonSources().size() + in.getMassAnalyzers().size() + in.getIonDetectors().size();
      if (component_count != 0)
      {
        os << "\t\t\t<componentList count=\"" << (std::max)((Size)3, component_count) << "\">\n";
        //--------------------------------------------------------------------------------------------
        // ion source
        //--------------------------------------------------------------------------------------------
        for (Size i = 0; i < in.getIonSources().size(); ++i)
        {
          const IonSource& so = in.getIonSources()[i];
          os << "\t\t\t\t<source order=\"" << so.getOrder() << "\">\n";

          if (so.getInletType() == IonSource::CONTINUOUSFLOWFASTATOMBOMBARDMENT)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000055\" name=\"continuous flow fast atom bombardment\" />\n";
          }
          else if (so.getInletType() == IonSource::DIRECT)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000056\" name=\"direct inlet\" />\n";
          }
          else if (so.getInletType() == IonSource::ELECTROSPRAYINLET)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000057\" name=\"electrospray inlet\" />\n";
          }
          else if (so.getInletType() == IonSource::FLOWINJECTIONANALYSIS)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000058\" name=\"flow injection analysis\" />\n";
          }
          else if (so.getInletType() == IonSource::INDUCTIVELYCOUPLEDPLASMA)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000059\" name=\"inductively coupled plasma\" />\n";
          }
          else if (so.getInletType() == IonSource::INFUSION)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000060\" name=\"infusion\" />\n";
          }
          else if (so.getInletType() == IonSource::JETSEPARATOR)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000061\" name=\"jet separator\" />\n";
          }
          else if (so.getInletType() == IonSource::MEMBRANESEPARATOR)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000062\" name=\"membrane separator\" />\n";
          }
          else if (so.getInletType() == IonSource::MOVINGBELT)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000063\" name=\"moving belt\" />\n";
          }
          else if (so.getInletType() == IonSource::MOVINGWIRE)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000064\" name=\"moving wire\" />\n";
          }
          else if (so.getInletType() == IonSource::OPENSPLIT)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000065\" name=\"open split\" />\n";
          }
          else if (so.getInletType() == IonSource::PARTICLEBEAM)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000066\" name=\"particle beam\" />\n";
          }
          else if (so.getInletType() == IonSource::RESERVOIR)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000067\" name=\"reservoir\" />\n";
          }
          else if (so.getInletType() == IonSource::SEPTUM)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000068\" name=\"septum\" />\n";
          }
          else if (so.getInletType() == IonSource::THERMOSPRAYINLET)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000069\" name=\"thermospray inlet\" />\n";
          }
          else if (so.getInletType() == IonSource::BATCH)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000248\" name=\"direct insertion probe\" />\n";
          }
          else if (so.getInletType() == IonSource::CHROMATOGRAPHY)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000249\" name=\"direct liquid introduction\" />\n";
          }
          else if (so.getInletType() == IonSource::MEMBRANE)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000396\" name=\"membrane inlet\" />\n";
          }
          else if (so.getInletType() == IonSource::NANOSPRAY)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000485\" name=\"nanospray inlet\" />\n";
          }

          if (so.getIonizationMethod() == IonSource::APCI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000070\" name=\"atmospheric pressure chemical ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::CI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000071\" name=\"chemical ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::ESI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000073\" name=\"electrospray ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::FAB)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000074\" name=\"fast atom bombardment ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::MALDI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000075\" name=\"matrix-assisted laser desorption ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::MPI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000227\" name=\"multiphoton ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::AP_MALDI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000239\" name=\"atmospheric pressure matrix-assisted laser desorption ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::API)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000240\" name=\"atmospheric pressure ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::DI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000247\" name=\"desorption ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::FA)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000255\" name=\"flowing afterglow\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::FD)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000257\" name=\"field desorption\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::FI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000258\" name=\"field ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::GD_MS)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000259\" name=\"glow discharge ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::NICI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000271\" name=\"Negative ion chemical ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::NRMS)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000272\" name=\"neutralization reionization mass spectrometry\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::PI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000273\" name=\"photoionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::PYMS)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000274\" name=\"pyrolysis mass spectrometry\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::REMPI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000276\" name=\"resonance enhanced multiphoton ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::SELDI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000278\" name=\"surface enhanced laser desorption ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::SEND)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000279\" name=\"surface enhanced neat desorption\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::AI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000380\" name=\"adiabatic ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::ASI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000381\" name=\"associative ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::APPI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000382\" name=\"atmospheric pressure photoionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::AD)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000383\" name=\"autodetachment\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::AUI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000384\" name=\"autoionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::CEI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000385\" name=\"charge exchange ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::CHEMI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000386\" name=\"chemi-ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::SILI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000387\" name=\"desorption/ionization on silicon\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::DISSI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000388\" name=\"dissociative ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::EI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000389\" name=\"electron ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::LD)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000393\" name=\"laser desorption ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::LSI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000395\" name=\"liquid secondary ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::MESI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000397\" name=\"microelectrospray\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::NESI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000398\" name=\"nanoelectrospray\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::PEI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000399\" name=\"penning ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::PD)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000400\" name=\"plasma desorption ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::SI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000402\" name=\"secondary ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::SOI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000403\" name=\"soft ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::SPI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000404\" name=\"spark ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::SALDI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000405\" name=\"surface-assisted laser desorption ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::SUI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000406\" name=\"surface ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::TI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000407\" name=\"thermal ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::VI)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000408\" name=\"vertical ionization\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::FIB)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000446\" name=\"fast ion bombardment\" />\n";
          }
          else if (so.getIonizationMethod() == IonSource::IONMETHODNULL)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000008\" name=\"ionization type\" />\n";
          }

          writeUserParam_(os, so, 5, "/mzML/instrumentConfigurationList/instrumentConfiguration/componentList/source/cvParam/@accession", validator);
          os << "\t\t\t\t</source>\n";
        }
        //FORCED
        if (component_count < 3 && in.getIonSources().empty())
        {
          os << "\t\t\t\t<source order=\"1234\">\n";
          os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000446\" name=\"fast ion bombardment\" />\n";
          os << "\t\t\t\t\t<userParam name=\"warning\" type=\"xsd:string\" value=\"invented ion source, to fulfill mzML schema\" />\n";
          os << "\t\t\t\t</source>\n";
        }
        //--------------------------------------------------------------------------------------------
        // mass analyzer
        //--------------------------------------------------------------------------------------------
        for (Size i = 0; i < in.getMassAnalyzers().size(); ++i)
        {
          const MassAnalyzer& ma = in.getMassAnalyzers()[i];
          os << "\t\t\t\t<analyzer order=\"" << ma.getOrder() << "\">\n";

          os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000014\" name=\"accuracy\" value=\"" << ma.getAccuracy() << "\" unitAccession=\"UO:0000169\" unitName=\"parts per million\" unitCvRef=\"UO\" />\n";
          // @todo: the parameters below are instrument specific and should not be written every time
          os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000022\" name=\"TOF Total Path Length\" value=\"" << ma.getTOFTotalPathLength() << "\" unitAccession=\"UO:0000008\" unitName=\"meter\" unitCvRef=\"UO\" />\n";
          os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000024\" name=\"final MS exponent\" value=\"" << ma.getFinalMSExponent() << "\" />\n";
          os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000025\" name=\"magnetic field strength\" value=\"" << ma.getMagneticFieldStrength() << "\" unitAccession=\"UO:0000228\" unitName=\"tesla\" unitCvRef=\"UO\" />\n";

          if (ma.getReflectronState() == MassAnalyzer::ON)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000106\" name=\"reflectron on\" />\n";

          }
          else if (ma.getReflectronState() == MassAnalyzer::OFF)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000105\" name=\"reflectron off\" />\n";
          }

          if (ma.getType() == MassAnalyzer::FOURIERTRANSFORM)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000079\" name=\"fourier transform ion cyclotron resonance mass spectrometer\" />\n";
          }
          else if (ma.getType() == MassAnalyzer::SECTOR)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000080\" name=\"magnetic sector\" />\n";
          }
          else if (ma.getType() == MassAnalyzer::QUADRUPOLE)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000081\" name=\"quadrupole\" />\n";
          }
          else if (ma.getType() == MassAnalyzer::TOF)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000084\" name=\"time-of-flight\" />\n";
          }
          else if (ma.getType() == MassAnalyzer::ESA)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000254\" name=\"electrostatic energy analyzer\" />\n";
          }
          else if (ma.getType() == MassAnalyzer::IT)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000264\" name=\"ion trap\" />\n";
          }
          else if (ma.getType() == MassAnalyzer::SWIFT)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000284\" name=\"stored waveform inverse fourier transform\" />\n";
          }
          else if (ma.getType() == MassAnalyzer::CYCLOTRON)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000288\" name=\"cyclotron\" />\n";
          }
          else if (ma.getType() == MassAnalyzer::ORBITRAP)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000484\" name=\"orbitrap\" />\n";
          }
          else if (ma.getType() == MassAnalyzer::AXIALEJECTIONLINEARIONTRAP)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000078\" name=\"axial ejection linear ion trap\" />\n";
          }
          else if (ma.getType() == MassAnalyzer::PAULIONTRAP)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000082\" name=\"quadrupole ion trap\" />\n";
          }
          else if (ma.getType() == MassAnalyzer::RADIALEJECTIONLINEARIONTRAP)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000083\" name=\"radial ejection linear ion trap\" />\n";
          }
          else if (ma.getType() == MassAnalyzer::LIT)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000291\" name=\"linear ion trap\" />\n";
          }
          else if (ma.getType() == MassAnalyzer::ANALYZERNULL)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000443\" name=\"mass analyzer type\" />\n";
          }

          writeUserParam_(os, ma, 5, "/mzML/instrumentConfigurationList/instrumentConfiguration/componentList/analyzer/cvParam/@accession", validator);
          os << "\t\t\t\t</analyzer>\n";
        }
        //FORCED
        if (component_count < 3 && in.getMassAnalyzers().empty())
        {
          os << "\t\t\t\t<analyzer order=\"1234\">\n";
          os << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000288\" name=\"cyclotron\" />\n";
          os << "\t\t\t\t\t<userParam name=\"warning\" type=\"xsd:string\" value=\"invented mass analyzer, to fulfill mzML schema\" />\n";
          os << "\t\t\t\t</analyzer>\n";
        }
        //--------------------------------------------------------------------------------------------
        // ion detector
        //--------------------------------------------------------------------------------------------
        for (Size i = 0; i < in.getIonDetectors().size(); ++i)
        {
          const IonDetector& id = in.getIonDetectors()[i];
          os << "\t\t\t\t<detector order=\"" << id.getOrder() << "\">\n";

          os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000028\" name=\"detector resolution\" value=\"" << id.getResolution() << "\" />\n";
          os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000029\" name=\"sampling frequency\" value=\"" << id.getADCSamplingFrequency() << "\" unitAccession=\"UO:0000106\" unitName=\"hertz\" unitCvRef=\"UO\" />\n";

          if (id.getAcquisitionMode() == IonDetector::ADC)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000117\" name=\"analog-digital converter\" />\n";
          }
          else if (id.getAcquisitionMode() == IonDetector::PULSECOUNTING)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000118\" name=\"pulse counting\" />\n";
          }
          else if (id.getAcquisitionMode() == IonDetector::TDC)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000119\" name=\"time-digital converter\" />\n";
          }
          else if (id.getAcquisitionMode() == IonDetector::TRANSIENTRECORDER)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000120\" name=\"transient recorder\" />\n";
          }

          if (id.getType() == IonDetector::CHANNELTRON)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000107\" name=\"channeltron\" />\n";
          }
          else if (id.getType() == IonDetector::DALYDETECTOR)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000110\" name=\"daly detector\" />\n";
          }
          else if (id.getType() == IonDetector::FARADAYCUP)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000112\" name=\"faraday cup\" />\n";
          }
          else if (id.getType() == IonDetector::MICROCHANNELPLATEDETECTOR)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000114\" name=\"microchannel plate detector\" />\n";
          }
          else if (id.getType() == IonDetector::MULTICOLLECTOR)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000115\" name=\"multi-collector\" />\n";
          }
          else if (id.getType() == IonDetector::PHOTOMULTIPLIER)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000116\" name=\"photomultiplier\" />\n";
          }
          else if (id.getType() == IonDetector::ELECTRONMULTIPLIER)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000253\" name=\"electron multiplier\" />\n";
          }
          else if (id.getType() == IonDetector::ARRAYDETECTOR)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000345\" name=\"array detector\" />\n";
          }
          else if (id.getType() == IonDetector::CONVERSIONDYNODE)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000346\" name=\"conversion dynode\" />\n";
          }
          else if (id.getType() == IonDetector::DYNODE)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000347\" name=\"dynode\" />\n";
          }
          else if (id.getType() == IonDetector::FOCALPLANECOLLECTOR)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000348\" name=\"focal plane collector\" />\n";
          }
          else if (id.getType() == IonDetector::IONTOPHOTONDETECTOR)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000349\" name=\"ion-to-photon detector\" />\n";
          }
          else if (id.getType() == IonDetector::POINTCOLLECTOR)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000350\" name=\"point collector\" />\n";
          }
          else if (id.getType() == IonDetector::POSTACCELERATIONDETECTOR)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000351\" name=\"postacceleration detector\" />\n";
          }
          else if (id.getType() == IonDetector::PHOTODIODEARRAYDETECTOR)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000621\" name=\"photodiode array detector\" />\n";
          }
          else if (id.getType() == IonDetector::INDUCTIVEDETECTOR)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000624\" name=\"inductive detector\" />\n";
          }
          else if (id.getType() == IonDetector::CONVERSIONDYNODEELECTRONMULTIPLIER)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000108\" name=\"conversion dynode electron multiplier\" />\n";
          }
          else if (id.getType() == IonDetector::CONVERSIONDYNODEPHOTOMULTIPLIER)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000109\" name=\"conversion dynode photomultiplier\" />\n";
          }
          else if (id.getType() == IonDetector::ELECTRONMULTIPLIERTUBE)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000111\" name=\"electron multiplier tube\" />\n";
          }
          else if (id.getType() == IonDetector::FOCALPLANEARRAY)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000113\" name=\"focal plane array\" />\n";
          }
          else if (id.getType() == IonDetector::TYPENULL)
          {
            os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000026\" name=\"detector type\" />\n";
          }

          writeUserParam_(os, id, 5, "/mzML/instrumentConfigurationList/instrumentConfiguration/componentList/detector/cvParam/@accession", validator);
          os << "\t\t\t\t</detector>\n";
        }
        //FORCED
        if (component_count < 3 && in.getIonDetectors().empty())
        {
          os << "\t\t\t\t<detector order=\"1234\">\n";
          os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000107\" name=\"channeltron\" />\n";
          os << "\t\t\t\t\t<userParam name=\"warning\" type=\"xsd:string\" value=\"invented ion detector, to fulfill mzML schema\" />\n";
          os << "\t\t\t\t</detector>\n";
        }
        os << "\t\t\t</componentList>\n";
      }
      os << "\t\t\t<softwareRef ref=\"so_in_0\" />\n";
      os << "\t\t</instrumentConfiguration>\n";
      os << "\t</instrumentConfigurationList>\n";

      //--------------------------------------------------------------------------------------------
      // data processing
      //--------------------------------------------------------------------------------------------

      // count number of float data array dps
      Size num_bi_dps(0);
      for (Size s = 0; s < exp.size(); ++s)
      {
        for (Size m = 0; m < exp[s].getFloatDataArrays().size(); ++m)
        {
          ++num_bi_dps;
        }
      }

      os << "\t<dataProcessingList count=\"" << (std::max)((Size)1, dps.size() + num_bi_dps) << "\">\n";

      // default (if experiment is empty and no actual data processing is here)
      if (dps.size() + num_bi_dps == 0)
      {
        std::vector< ConstDataProcessingPtr > dummy;
        writeDataProcessing_(os, "dp_sp_0", dummy, validator);
      }

      for (Size s = 0; s < dps.size(); ++s)
      {
        writeDataProcessing_(os, String("dp_sp_") + s, dps[s], validator);
      }

      //for each binary data array
      for (Size s = 0; s < exp.size(); ++s)
      {
        for (Size m = 0; m < exp[s].getFloatDataArrays().size(); ++m)
        {
          // if a DataArray has dataProcessing information, write it, otherwise we assume it has the
          // same processing as the rest of the spectra and use the implicit referencing of mzML
          // to the first entry (which is a dummy if none exists; see above)
          if (!exp[s].getFloatDataArrays()[m].getDataProcessing().empty())
          {
            writeDataProcessing_(os, String("dp_sp_") + s + "_bi_" + m,
              exp[s].getFloatDataArrays()[m].getDataProcessing(), validator);
          }
        }
      }

      os << "\t</dataProcessingList>\n";
      //--------------------------------------------------------------------------------------------
      // acquisitionSettings
      //--------------------------------------------------------------------------------------------

      //--------------------------------------------------------------------------------------------
      // run
      //--------------------------------------------------------------------------------------------
      os << "\t<run id=\"ru_0\" defaultInstrumentConfigurationRef=\"ic_0\" sampleRef=\"sa_0\"";
      if (exp.getDateTime().isValid())
      {
        os << " startTimeStamp=\"" << exp.getDateTime().get().substitute(' ', 'T') << "\"";
      }
      if (!exp.getSourceFiles().empty())
      {
        os << " defaultSourceFileRef=\"sf_ru_0\"";
      }
      os << ">\n";

      //run attributes
      if (!exp.getFractionIdentifier().empty())
      {
        os << "\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000858\" name=\"fraction identifier\" value=\"" << exp.getFractionIdentifier() << "\" />\n";
      }

      writeUserParam_(os, exp, 2, "/mzML/run/cvParam/@accession", validator);

    }

    void MzMLHandler::writeSpectrum_(std::ostream& os,
                                     const SpectrumType& spec,
                                     Size s,
                                     const Internal::MzMLValidator& validator,
                                     bool renew_native_ids,
                                     std::vector<std::vector< ConstDataProcessingPtr > >& dps)
    {
      //native id
      String native_id = spec.getNativeID();
      if (renew_native_ids)
      {
        native_id = String("spectrum=") + s;
      }

      Int64 offset = os.tellp();
      spectra_offsets_.emplace_back(native_id, offset + 3);

      // IMPORTANT make sure the offset (above) corresponds to the start of the <spectrum tag
      os << "\t\t\t<spectrum id=\"" << writeXMLEscape(native_id) << "\" index=\"" << s << "\" defaultArrayLength=\"" << spec.size() << "\"";
      if (spec.getSourceFile() != SourceFile())
      {
        os << " sourceFileRef=\"sf_sp_" << s << "\"";
      }
      //the data processing info of the first spectrum is the default
      //if (s==0 || spec.getDataProcessing()!=exp[0].getDataProcessing())
      if (s == 0 || spec.getDataProcessing() != dps[0])
      {
        Size dp_ref_num = s;
        if (s != 0)
        {
          for (Size i = 0; i < dps.size(); ++i)
          {
            if (spec.getDataProcessing() == dps[i])
            {
              dp_ref_num = i;
              break;
            }
          }
        }
        os << " dataProcessingRef=\"dp_sp_" << dp_ref_num << "\"";
      }
      os << ">\n";

      //spectrum representation
      if (spec.getType() == SpectrumSettings::CENTROID)
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000127\" name=\"centroid spectrum\" />\n";
      }
      else if (spec.getType() == SpectrumSettings::PROFILE)
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000128\" name=\"profile spectrum\" />\n";
      }
      else
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000525\" name=\"spectrum representation\" />\n";
      }

      //spectrum attributes
      if (spec.getMSLevel() != 0)
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000511\" name=\"ms level\" value=\"" << spec.getMSLevel() << "\" />\n";
      }

      //spectrum type
      if (spec.getInstrumentSettings().getScanMode() == InstrumentSettings::MASSSPECTRUM)
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000294\" name=\"mass spectrum\" />\n";
      }
      else if (spec.getInstrumentSettings().getScanMode() == InstrumentSettings::MS1SPECTRUM)
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000579\" name=\"MS1 spectrum\" />\n";
      }
      else if (spec.getInstrumentSettings().getScanMode() == InstrumentSettings::MSNSPECTRUM)
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" />\n";
      }
      else if (spec.getInstrumentSettings().getScanMode() == InstrumentSettings::SIM)
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000582\" name=\"SIM spectrum\" />\n";
      }
      else if (spec.getInstrumentSettings().getScanMode() == InstrumentSettings::SRM)
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000583\" name=\"SRM spectrum\" />\n";
      }
      else if (spec.getInstrumentSettings().getScanMode() == InstrumentSettings::CRM)
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000581\" name=\"CRM spectrum\" />\n";
      }
      else if (spec.getInstrumentSettings().getScanMode() == InstrumentSettings::PRECURSOR)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000341\" name=\"precursor ion spectrum\" />\n";
      }
      else if (spec.getInstrumentSettings().getScanMode() == InstrumentSettings::CNG)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000325\" name=\"constant neutral gain spectrum\" />\n";
      }
      else if (spec.getInstrumentSettings().getScanMode() == InstrumentSettings::CNL)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000326\" name=\"constant neutral loss spectrum\" />\n";
      }
      else if (spec.getInstrumentSettings().getScanMode() == InstrumentSettings::EMR)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000804\" name=\"electromagnetic radiation spectrum\" />\n";
      }
      else if (spec.getInstrumentSettings().getScanMode() == InstrumentSettings::EMISSION)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000805\" name=\"emission spectrum\" />\n";
      }
      else if (spec.getInstrumentSettings().getScanMode() == InstrumentSettings::ABSORPTION)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000806\" name=\"absorption spectrum\" />\n";
      }
      else if (spec.getInstrumentSettings().getScanMode() == InstrumentSettings::EMC)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000789\" name=\"enhanced multiply charged spectrum\" />\n";
      }
      else if (spec.getInstrumentSettings().getScanMode() == InstrumentSettings::TDF)
      {
        os << "\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000789\" name=\"time-delayed fragmentation spectrum\" />\n";
      }
      else   //FORCED
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000294\" name=\"mass spectrum\" />\n";
      }

      //scan polarity
      if (spec.getInstrumentSettings().getPolarity() == IonSource::NEGATIVE)
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000129\" name=\"negative scan\" />\n";
      }
      else if (spec.getInstrumentSettings().getPolarity() == IonSource::POSITIVE)
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000130\" name=\"positive scan\" />\n";
      }

      writeUserParam_(os, spec, 4, "/mzML/run/spectrumList/spectrum/cvParam/@accession", validator);
      //--------------------------------------------------------------------------------------------
      //scan list
      //--------------------------------------------------------------------------------------------
      os << "\t\t\t\t<scanList count=\"" << (std::max)((Size)1, spec.getAcquisitionInfo().size()) << "\">\n";
      ControlledVocabulary::CVTerm ai_term = getChildWithName_("MS:1000570", spec.getAcquisitionInfo().getMethodOfCombination());
      if (!ai_term.id.empty())
      {
        os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"" << ai_term.id << "\" name=\"" << ai_term.name << "\" />\n";
      }
      else
      {
        os << "\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000795\" name=\"no combination\" />\n";
      }
      writeUserParam_(os, spec.getAcquisitionInfo(), 5, "/mzML/run/spectrumList/spectrum/scanList/cvParam/@accession", validator);

      //--------------------------------------------------------------------------------------------
      //scan
      //--------------------------------------------------------------------------------------------
      for (Size j = 0; j < spec.getAcquisitionInfo().size(); ++j)
      {
        const Acquisition& ac = spec.getAcquisitionInfo()[j];
        os << "\t\t\t\t\t<scan "; // TODO
        if (!ac.getIdentifier().empty())
        {
          os << "externalSpectrumID=\"" << ac.getIdentifier() << "\"";
        }
        os << ">\n";
        if (j == 0)
        {
          os << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000016\" name=\"scan start time\" value=\"" << spec.getRT()
             << "\" unitAccession=\"UO:0000010\" unitName=\"second\" unitCvRef=\"UO\" />\n";

          if (spec.getDriftTimeUnit() == DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE)
          {
            os << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1001581\" name=\"FAIMS compensation voltage\" value=\"" << spec.getDriftTime()
                << "\" unitAccession=\"UO:000218\" unitName=\"volt\" unitCvRef=\"UO\" />\n";
          }          
          else if (spec.getDriftTime() != IMTypes::DRIFTTIME_NOT_SET)// if drift time was never set, don't report it
          {
            if (spec.getDriftTimeUnit() == DriftTimeUnit::MILLISECOND)
            {
              os << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1002476\" name=\"ion mobility drift time\" value=\"" << spec.getDriftTime()
                 << "\" unitAccession=\"UO:0000028\" unitName=\"millisecond\" unitCvRef=\"UO\" />\n";
            }
            else if (spec.getDriftTimeUnit() == DriftTimeUnit::VSSC)
            {
              os << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1002815\" name=\"inverse reduced ion mobility\" value=\"" << spec.getDriftTime()
                 << "\" unitAccession=\"MS:1002814\" unitName=\"volt-second per square centimeter\" unitCvRef=\"MS\" />\n";
            }
            else
            {
              // assume milliseconds, but warn
              warning(STORE, String("Spectrum drift time unit not set, assume milliseconds"));
              os << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1002476\" name=\"ion mobility drift time\" value=\"" << spec.getDriftTime()
                 << "\" unitAccession=\"UO:0000028\" unitName=\"millisecond\" unitCvRef=\"UO\" />\n";
            }
          }
        }
        writeUserParam_(os, ac, 6, "/mzML/run/spectrumList/spectrum/scanList/scan/cvParam/@accession", validator);

        if (spec.getInstrumentSettings().getZoomScan())
        {
          os << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000497\" name=\"zoom scan\" />\n";
        }

        //scan windows
        if (j == 0 && !spec.getInstrumentSettings().getScanWindows().empty())
        {
          os << "\t\t\t\t\t\t<scanWindowList count=\"" << spec.getInstrumentSettings().getScanWindows().size() << "\">\n";
          for (Size k = 0; k < spec.getInstrumentSettings().getScanWindows().size(); ++k)
          {
            os << "\t\t\t\t\t\t\t<scanWindow>\n";
            os << "\t\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000501\" name=\"scan window lower limit\" value=\"" << spec.getInstrumentSettings().getScanWindows()[k].begin << "\" unitAccession=\"MS:1000040\" unitName=\"m/z\" unitCvRef=\"MS\" />\n";
            os << "\t\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000500\" name=\"scan window upper limit\" value=\"" << spec.getInstrumentSettings().getScanWindows()[k].end << "\" unitAccession=\"MS:1000040\" unitName=\"m/z\" unitCvRef=\"MS\" />\n";
            writeUserParam_(os, spec.getInstrumentSettings().getScanWindows()[k], 8, "/mzML/run/spectrumList/spectrum/scanList/scan/scanWindowList/scanWindow/cvParam/@accession", validator);
            os << "\t\t\t\t\t\t\t</scanWindow>\n";
          }
          os << "\t\t\t\t\t\t</scanWindowList>\n";
        }
        os << "\t\t\t\t\t</scan>\n";
      }
      //fallback if we have no acquisition information (a dummy scan is created for RT and so on)
      if (spec.getAcquisitionInfo().empty())
      {
        os << "\t\t\t\t\t<scan>\n";
        os << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000016\" name=\"scan start time\" value=\"" << spec.getRT() << "\" unitAccession=\"UO:0000010\" unitName=\"second\" unitCvRef=\"UO\" />\n";

        if (spec.getInstrumentSettings().getZoomScan())
        {
          os << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000497\" name=\"zoom scan\" />\n";
        }
        //scan windows
        if (!spec.getInstrumentSettings().getScanWindows().empty())
        {
          os << "\t\t\t\t\t\t<scanWindowList count=\"" << spec.getInstrumentSettings().getScanWindows().size() << "\">\n";
          for (Size j = 0; j < spec.getInstrumentSettings().getScanWindows().size(); ++j)
          {
            os << "\t\t\t\t\t\t\t<scanWindow>\n";
            os << "\t\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000501\" name=\"scan window lower limit\" value=\"" << spec.getInstrumentSettings().getScanWindows()[j].begin << "\" unitAccession=\"MS:1000040\" unitName=\"m/z\" unitCvRef=\"MS\" />\n";
            os << "\t\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000500\" name=\"scan window upper limit\" value=\"" << spec.getInstrumentSettings().getScanWindows()[j].end << "\" unitAccession=\"MS:1000040\" unitName=\"m/z\" unitCvRef=\"MS\" />\n";
            writeUserParam_(os, spec.getInstrumentSettings().getScanWindows()[j], 8, "/mzML/run/spectrumList/spectrum/scanList/scan/scanWindowList/scanWindow/cvParam/@accession", validator);
            os << "\t\t\t\t\t\t\t</scanWindow>\n";
          }
          os << "\t\t\t\t\t\t</scanWindowList>\n";
        }
        os << "\t\t\t\t\t</scan>\n";
      }
      os << "\t\t\t\t</scanList>\n";

      //--------------------------------------------------------------------------------------------
      //precursor list
      //--------------------------------------------------------------------------------------------
      if (!spec.getPrecursors().empty())
      {
        os << "\t\t\t\t<precursorList count=\"" << spec.getPrecursors().size() << "\">\n";
        for (Size p = 0; p != spec.getPrecursors().size(); ++p)
        {
          writePrecursor_(os, spec.getPrecursors()[p], validator);
        }
        os << "\t\t\t\t</precursorList>\n";
      }

      //--------------------------------------------------------------------------------------------
      //product list
      //--------------------------------------------------------------------------------------------
      if (!spec.getProducts().empty())
      {
        os << "\t\t\t\t<productList count=\"" << spec.getProducts().size() << "\">\n";
        for (Size p = 0; p < spec.getProducts().size(); ++p)
        {
          writeProduct_(os, spec.getProducts()[p], validator);
        }
        os << "\t\t\t\t</productList>\n";
      }

      //--------------------------------------------------------------------------------------------
      //binary data array list
      //--------------------------------------------------------------------------------------------
      if (!spec.empty())
      {
        String encoded_string;
        os << "\t\t\t\t<binaryDataArrayList count=\"" << (2 + spec.getFloatDataArrays().size() + spec.getStringDataArrays().size() + spec.getIntegerDataArrays().size()) << "\">\n";

        writeContainerData_<SpectrumType>(os, options_, spec, "mz");
        writeContainerData_<SpectrumType>(os, options_, spec, "intensity");

        String compression_term = MzMLHandlerHelper::getCompressionTerm_(options_, options_.getNumpressConfigurationIntensity(), "\t\t\t\t\t\t", false);
        // write float data array
        for (Size m = 0; m < spec.getFloatDataArrays().size(); ++m)
        {
          const SpectrumType::FloatDataArray& array = spec.getFloatDataArrays()[m];
          writeBinaryFloatDataArray_(os, options_, array, s, m, true, validator);
        }
        // write integer data array
        for (Size m = 0; m < spec.getIntegerDataArrays().size(); ++m)
        {
          const SpectrumType::IntegerDataArray& array = spec.getIntegerDataArrays()[m];
          std::vector<Int64> data64_to_encode(array.size());
          for (Size p = 0; p < array.size(); ++p)
          {
            data64_to_encode[p] = array[p];
          }
          Base64::encodeIntegers(data64_to_encode, Base64::BYTEORDER_LITTLEENDIAN, encoded_string, options_.getCompression());

          String data_processing_ref_string = "";
          if (!array.getDataProcessing().empty())
          {
            data_processing_ref_string = String("dataProcessingRef=\"dp_sp_") + s + "_bi_" + m + "\"";
          }
          os << "\t\t\t\t\t<binaryDataArray arrayLength=\"" << array.size() << "\" encodedLength=\"" << encoded_string.size() << "\" " << data_processing_ref_string << ">\n";
          os << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000522\" name=\"64-bit integer\" />\n";
          os << "\t\t\t\t\t\t" << compression_term << "\n";
          ControlledVocabulary::CVTerm bi_term = getChildWithName_("MS:1000513", array.getName());
          if (!bi_term.id.empty())
          {
            os << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"" << bi_term.id << "\" name=\"" << bi_term.name << "\" />\n";
          }
          else
          {
            os << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000786\" name=\"non-standard data array\" value=\"" << array.getName() << "\" />\n";
          }
          writeUserParam_(os, array, 6, "/mzML/run/spectrumList/spectrum/binaryDataArrayList/binaryDataArray/cvParam/@accession", validator);
          os << "\t\t\t\t\t\t<binary>" << encoded_string << "</binary>\n";
          os << "\t\t\t\t\t</binaryDataArray>\n";
        }
        // write string data arrays
        for (Size m = 0; m < spec.getStringDataArrays().size(); ++m)
        {
          const SpectrumType::StringDataArray& array = spec.getStringDataArrays()[m];
          std::vector<String> data_to_encode;
          data_to_encode.resize(array.size());
          for (Size p = 0; p < array.size(); ++p)
            data_to_encode[p] = array[p];
          Base64::encodeStrings(data_to_encode, encoded_string, options_.getCompression());
          String data_processing_ref_string = "";
          if (!array.getDataProcessing().empty())
          {
            data_processing_ref_string = String("dataProcessingRef=\"dp_sp_") + s + "_bi_" + m + "\"";
          }
          os << "\t\t\t\t\t<binaryDataArray arrayLength=\"" << array.size() << "\" encodedLength=\"" << encoded_string.size() << "\" " << data_processing_ref_string << ">\n";
          os << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1001479\" name=\"null-terminated ASCII string\" />\n";
          os << "\t\t\t\t\t\t" << compression_term << "\n";
          os << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000786\" name=\"non-standard data array\" value=\"" << array.getName() << "\" />\n";
          writeUserParam_(os, array, 6, "/mzML/run/spectrumList/spectrum/binaryDataArrayList/binaryDataArray/cvParam/@accession", validator);
          os << "\t\t\t\t\t\t<binary>" << encoded_string << "</binary>\n";
          os << "\t\t\t\t\t</binaryDataArray>\n";
        }
        os << "\t\t\t\t</binaryDataArrayList>\n";
      }

      os << "\t\t\t</spectrum>\n";
    }

    template <typename ContainerT>
    void MzMLHandler::writeContainerData_(std::ostream& os, const PeakFileOptions& pf_options_, const ContainerT& container, const String& array_type)
    {
      // Intensity is the same for chromatograms and spectra, the second
      // dimension is either "time" or "mz" (both of these are controlled by
      // getMz32Bit)
      bool is32Bit = ((array_type == "intensity" && pf_options_.getIntensity32Bit()) || pf_options_.getMz32Bit());
      if (!is32Bit || pf_options_.getNumpressConfigurationMassTime().np_compression != MSNumpressCoder::NONE)
      {
        std::vector<double> data_to_encode(container.size());
        if (array_type == "intensity")
        {
          for (Size p = 0; p < container.size(); ++p)
          {
            data_to_encode[p] = container[p].getIntensity();
          }
        }
        else
        {
          for (Size p = 0; p < container.size(); ++p)
          {
            data_to_encode[p] = container[p].getPos();
          }
        }
        writeBinaryDataArray_(os, pf_options_, data_to_encode, false, array_type);
      }
      else
      {
        std::vector<float> data_to_encode(container.size());

        if (array_type == "intensity")
        {
          for (Size p = 0; p < container.size(); ++p)
          {
            data_to_encode[p] = container[p].getIntensity();
          }
        }
        else
        {
          for (Size p = 0; p < container.size(); ++p)
          {
            data_to_encode[p] = container[p].getPos();
          }
        }
        writeBinaryDataArray_(os, pf_options_, data_to_encode, true, array_type);
      }

    }

    template <typename DataType>
    void MzMLHandler::writeBinaryDataArray_(std::ostream& os,
                                            const PeakFileOptions& pf_options_,
                                            std::vector<DataType>& data_to_encode,
                                            bool is32bit,
                                            String array_type)
    {
      String encoded_string;
      bool no_numpress = true;

      // Compute the array-type and the compression CV term
      String cv_term_type;
      String compression_term;
      String compression_term_no_np;
      MSNumpressCoder::NumpressConfig np_config;
      if (array_type == "mz")
      {
        cv_term_type = "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000514\" name=\"m/z array\" unitAccession=\"MS:1000040\" unitName=\"m/z\" unitCvRef=\"MS\" />\n";
        compression_term = MzMLHandlerHelper::getCompressionTerm_(pf_options_, pf_options_.getNumpressConfigurationMassTime(), "\t\t\t\t\t\t", true);
        compression_term_no_np = MzMLHandlerHelper::getCompressionTerm_(pf_options_, pf_options_.getNumpressConfigurationMassTime(), "\t\t\t\t\t\t", false);
        np_config = pf_options_.getNumpressConfigurationMassTime();
      }
      else if (array_type == "time")
      {
        cv_term_type = "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000595\" name=\"time array\" unitAccession=\"UO:0000010\" unitName=\"second\" unitCvRef=\"MS\" />\n";
        compression_term = MzMLHandlerHelper::getCompressionTerm_(pf_options_, pf_options_.getNumpressConfigurationMassTime(), "\t\t\t\t\t\t", true);
        compression_term_no_np = MzMLHandlerHelper::getCompressionTerm_(pf_options_, pf_options_.getNumpressConfigurationMassTime(), "\t\t\t\t\t\t", false);
        np_config = pf_options_.getNumpressConfigurationMassTime();
      }
      else if (array_type == "intensity")
      {
        cv_term_type = "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\" unitAccession=\"MS:1000131\" unitName=\"number of detector counts\" unitCvRef=\"MS\"/>\n";
        compression_term = MzMLHandlerHelper::getCompressionTerm_(pf_options_, pf_options_.getNumpressConfigurationIntensity(), "\t\t\t\t\t\t", true);
        compression_term_no_np = MzMLHandlerHelper::getCompressionTerm_(pf_options_, pf_options_.getNumpressConfigurationIntensity(), "\t\t\t\t\t\t", false);
        np_config = pf_options_.getNumpressConfigurationIntensity();
      }
      else
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown array type", array_type);
      }

      // Try numpress encoding (if it is enabled) and fall back to regular encoding if it fails
      if (np_config.np_compression != MSNumpressCoder::NONE)
      {
        MSNumpressCoder().encodeNP(data_to_encode, encoded_string, pf_options_.getCompression(), np_config);
        if (!encoded_string.empty())
        {
          // numpress succeeded
          no_numpress = false;
          os << "\t\t\t\t\t<binaryDataArray encodedLength=\"" << encoded_string.size() << "\">\n";
          os << cv_term_type;
          os << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\" />\n";
        }
      }

      // Regular DataArray without numpress (either 32 or 64 bit encoded)
      if (is32bit && no_numpress)
      {
        compression_term = compression_term_no_np; // select the no-numpress term
        Base64::encode(data_to_encode, Base64::BYTEORDER_LITTLEENDIAN, encoded_string, pf_options_.getCompression());
        os << "\t\t\t\t\t<binaryDataArray encodedLength=\"" << encoded_string.size() << "\">\n";
        os << cv_term_type;
        os << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000521\" name=\"32-bit float\" />\n";
      }
      else if (!is32bit && no_numpress)
      {
        compression_term = compression_term_no_np; // select the no-numpress term
        Base64::encode(data_to_encode, Base64::BYTEORDER_LITTLEENDIAN, encoded_string, pf_options_.getCompression());
        os << "\t\t\t\t\t<binaryDataArray encodedLength=\"" << encoded_string.size() << "\">\n";
        os << cv_term_type;
        os << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\" />\n";
      }

      os << compression_term << "\n";
      os << "\t\t\t\t\t\t<binary>" << encoded_string << "</binary>\n";
      os << "\t\t\t\t\t</binaryDataArray>\n";
    }

    void MzMLHandler::writeBinaryFloatDataArray_(std::ostream& os,
                                                 const PeakFileOptions& pf_options_,
                                                 const OpenMS::DataArrays::FloatDataArray& array,
                                                 const Size spec_chrom_idx,
                                                 const Size array_idx,
                                                 bool isSpectrum,
                                                 const Internal::MzMLValidator& validator)
    {
      String encoded_string;
      bool no_numpress = true;
      std::vector<float> data_to_encode = array;
      MetaInfoDescription array_metadata = array;
      // bool is32bit = true;

      // Compute the array-type and the compression CV term
      String cv_term_type;
      String compression_term;
      String compression_term_no_np;
      MSNumpressCoder::NumpressConfig np_config;
      // if (array_type == "float_data")
      {
        // Try and identify whether we have a CV term for this particular array (otherwise write the array name itself)
        ControlledVocabulary::CVTerm bi_term = getChildWithName_("MS:1000513", array.getName()); // name: binary data array

        String unit_cv_term = "";
        if (array_metadata.metaValueExists("unit_accession"))
        {
          ControlledVocabulary::CVTerm unit = cv_.getTerm(array_metadata.getMetaValue("unit_accession"));
          unit_cv_term = " unitAccession=\"" + unit.id + "\" unitName=\"" + unit.name + "\" unitCvRef=\"" + unit.id.prefix(2) + "\"";
          array_metadata.removeMetaValue("unit_accession"); // prevent this from being written as userParam
        }

        if (!bi_term.id.empty())
        {
          cv_term_type = "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"" + bi_term.id + "\" name=\"" + bi_term.name + "\"" + unit_cv_term + " />\n";
        }
        else
        {
          cv_term_type = "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000786\" name=\"non-standard data array\" value=\"" +
            array.getName() + "\"" + unit_cv_term + " />\n";
        }

        compression_term = MzMLHandlerHelper::getCompressionTerm_(pf_options_, pf_options_.getNumpressConfigurationFloatDataArray(), "\t\t\t\t\t\t", true);
        compression_term_no_np = MzMLHandlerHelper::getCompressionTerm_(pf_options_, pf_options_.getNumpressConfigurationFloatDataArray(), "\t\t\t\t\t\t", false);
        np_config = pf_options_.getNumpressConfigurationFloatDataArray();
      }

      String data_processing_ref_string = "";
      if (!array.getDataProcessing().empty())
      {
        data_processing_ref_string = String("dataProcessingRef=\"dp_sp_") + spec_chrom_idx + "_bi_" + array_idx + "\"";
      }

      // Try numpress encoding (if it is enabled) and fall back to regular encoding if it fails
      if (np_config.np_compression != MSNumpressCoder::NONE)
      {
        MSNumpressCoder().encodeNP(data_to_encode, encoded_string, pf_options_.getCompression(), np_config);
        if (!encoded_string.empty())
        {
          // numpress succeeded
          no_numpress = false;
          os << "\t\t\t\t\t<binaryDataArray arrayLength=\"" << array.size() << "\" encodedLength=\"" << encoded_string.size() << "\" " << data_processing_ref_string << ">\n";
          os << cv_term_type;
          os << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\" />\n";
        }
      }

      // Regular DataArray without numpress (here: only 32 bit encoded)
      if (no_numpress)
      {
        compression_term = compression_term_no_np; // select the no-numpress term
        Base64::encode(data_to_encode, Base64::BYTEORDER_LITTLEENDIAN, encoded_string, pf_options_.getCompression());
        os << "\t\t\t\t\t<binaryDataArray arrayLength=\"" << array.size() << "\" encodedLength=\"" << encoded_string.size() << "\" " << data_processing_ref_string << ">\n";
        os << cv_term_type;
        os << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000521\" name=\"32-bit float\" />\n";
      }

      os << compression_term << "\n";
      if (isSpectrum)
      {
        writeUserParam_(os, array_metadata, 6, "/mzML/run/spectrumList/spectrum/binaryDataArrayList/binaryDataArray/cvParam/@accession", validator);
      }
      else
      {
        writeUserParam_(os, array_metadata, 6, "/mzML/run/chromatogramList/chromatogram/binaryDataArrayList/binaryDataArray/cvParam/@accession", validator);
      }
      os << "\t\t\t\t\t\t<binary>" << encoded_string << "</binary>\n";
      os << "\t\t\t\t\t</binaryDataArray>\n";
    }

    // We only ever need 2 instances for the following functions: one for Spectra / Chromatograms and one for floats / doubles
    template void MzMLHandler::writeContainerData_<SpectrumType>(std::ostream& os,
                                                                 const PeakFileOptions& pf_options_,
                                                                 const SpectrumType& container,
                                                                 const String& array_type);

    template void MzMLHandler::writeContainerData_<ChromatogramType>(std::ostream& os,
                                                                     const PeakFileOptions& pf_options_,
                                                                     const ChromatogramType& container,
                                                                     const String& array_type);

    template void MzMLHandler::writeBinaryDataArray_<float>(std::ostream& os,
                                                            const PeakFileOptions& pf_options_,
                                                            std::vector<float>& data_to_encode,
                                                            bool is32bit,
                                                            String array_type);

    template void MzMLHandler::writeBinaryDataArray_<double>(std::ostream& os,
                                                             const PeakFileOptions& pf_options_,
                                                             std::vector<double>& data_to_encode,
                                                             bool is32bit,
                                                             String array_type);

    void MzMLHandler::writeChromatogram_(std::ostream& os,
                                         const ChromatogramType& chromatogram,
                                         Size c,
                                         const Internal::MzMLValidator& validator)
    {
      Int64 offset = os.tellp();
      chromatograms_offsets_.emplace_back(chromatogram.getNativeID(), offset + 3);

      // TODO native id with chromatogram=?? prefix?
      // IMPORTANT make sure the offset (above) corresponds to the start of the <chromatogram tag
      os << "\t\t\t<chromatogram id=\"" << writeXMLEscape(chromatogram.getNativeID()) << "\" index=\"" << c << "\" defaultArrayLength=\"" << chromatogram.size() << "\">" << "\n";

      // write cvParams (chromatogram type)
      if (chromatogram.getChromatogramType() == ChromatogramSettings::MASS_CHROMATOGRAM)
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000810\" name=\"ion current chromatogram\" />\n";
      }
      else if (chromatogram.getChromatogramType() == ChromatogramSettings::TOTAL_ION_CURRENT_CHROMATOGRAM)
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000235\" name=\"total ion current chromatogram\" />\n";
      }
      else if (chromatogram.getChromatogramType() == ChromatogramSettings::SELECTED_ION_CURRENT_CHROMATOGRAM)
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000627\" name=\"selected ion current chromatogram\" />\n";
      }
      else if (chromatogram.getChromatogramType() == ChromatogramSettings::BASEPEAK_CHROMATOGRAM)
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000628\" name=\"basepeak chromatogram\" />\n";
      }
      else if (chromatogram.getChromatogramType() == ChromatogramSettings::SELECTED_ION_MONITORING_CHROMATOGRAM)
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1001472\" name=\"selected ion monitoring chromatogram\" />\n";
      }
      else if (chromatogram.getChromatogramType() == ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM)
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1001473\" name=\"selected reaction monitoring chromatogram\" />\n";
      }
      else if (chromatogram.getChromatogramType() == ChromatogramSettings::ELECTROMAGNETIC_RADIATION_CHROMATOGRAM)
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000811\" name=\"electromagnetic radiation chromatogram\" />\n";
      }
      else if (chromatogram.getChromatogramType() == ChromatogramSettings::ABSORPTION_CHROMATOGRAM)
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000812\" name=\"absorption chromatogram\" />\n";
      }
      else if (chromatogram.getChromatogramType() == ChromatogramSettings::EMISSION_CHROMATOGRAM)
      {
        os << "\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000813\" name=\"emission chromatogram\" />\n";
      }
      else
      {
        // TODO
      }
      writePrecursor_(os, chromatogram.getPrecursor(), validator);
      writeProduct_(os, chromatogram.getProduct(), validator);

      //--------------------------------------------------------------------------------------------
      //binary data array list
      //--------------------------------------------------------------------------------------------
      String compression_term;
      String encoded_string;
      os << "\t\t\t\t<binaryDataArrayList count=\"" << (2 + chromatogram.getFloatDataArrays().size() + chromatogram.getStringDataArrays().size() + chromatogram.getIntegerDataArrays().size()) << "\">\n";

      writeContainerData_<ChromatogramType>(os, options_, chromatogram, "time");
      writeContainerData_<ChromatogramType>(os, options_, chromatogram, "intensity");

      compression_term = MzMLHandlerHelper::getCompressionTerm_(options_, options_.getNumpressConfigurationIntensity(), "\t\t\t\t\t\t", false);
      // write float data array
      for (Size m = 0; m < chromatogram.getFloatDataArrays().size(); ++m)
      {
        const ChromatogramType::FloatDataArray& array = chromatogram.getFloatDataArrays()[m];
        writeBinaryFloatDataArray_(os, options_, array, c, m, false, validator);
      }
      //write integer data array
      for (Size m = 0; m < chromatogram.getIntegerDataArrays().size(); ++m)
      {
        const ChromatogramType::IntegerDataArray& array = chromatogram.getIntegerDataArrays()[m];
        std::vector<Int64> data64_to_encode(array.size());
        for (Size p = 0; p < array.size(); ++p)
        {
          data64_to_encode[p] = array[p];
        }
        Base64::encodeIntegers(data64_to_encode, Base64::BYTEORDER_LITTLEENDIAN, encoded_string, options_.getCompression());
        String data_processing_ref_string = "";
        if (!array.getDataProcessing().empty())
        {
          data_processing_ref_string = String("dataProcessingRef=\"dp_sp_") + c + "_bi_" + m + "\"";
        }
        os << "\t\t\t\t\t<binaryDataArray arrayLength=\"" << array.size() << "\" encodedLength=\"" << encoded_string.size() << "\" " << data_processing_ref_string << ">\n";
        os << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000522\" name=\"64-bit integer\" />\n";
        os << "\t\t\t\t\t\t" << compression_term << "\n";
        ControlledVocabulary::CVTerm bi_term = getChildWithName_("MS:1000513", array.getName());
        if (!bi_term.id.empty())
        {
          os << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"" << bi_term.id << "\" name=\"" << bi_term.name << "\" />\n";
        }
        else
        {
          os << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000786\" name=\"non-standard data array\" value=\"" << array.getName() << "\" />\n";
        }
        writeUserParam_(os, array, 6, "/mzML/run/chromatogramList/chromatogram/binaryDataArrayList/binaryDataArray/cvParam/@accession", validator);
        os << "\t\t\t\t\t\t<binary>" << encoded_string << "</binary>\n";
        os << "\t\t\t\t\t</binaryDataArray>\n";
      }
      //write string data arrays
      for (Size m = 0; m < chromatogram.getStringDataArrays().size(); ++m)
      {
        const ChromatogramType::StringDataArray& array = chromatogram.getStringDataArrays()[m];
        std::vector<String> data_to_encode;
        data_to_encode.resize(array.size());
        for (Size p = 0; p < array.size(); ++p)
        {
          data_to_encode[p] = array[p];
        }
        Base64::encodeStrings(data_to_encode, encoded_string, options_.getCompression());
        String data_processing_ref_string = "";
        if (!array.getDataProcessing().empty())
        {
          data_processing_ref_string = String("dataProcessingRef=\"dp_sp_") + c + "_bi_" + m + "\"";
        }
        os << "\t\t\t\t\t<binaryDataArray arrayLength=\"" << array.size() << "\" encodedLength=\"" << encoded_string.size() << "\" " << data_processing_ref_string << ">\n";
        os << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1001479\" name=\"null-terminated ASCII string\" />\n";
        os << "\t\t\t\t\t\t" << compression_term << "\n";
        os << "\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000786\" name=\"non-standard data array\" value=\"" << array.getName() << "\" />\n";
        writeUserParam_(os, array, 6, "/mzML/run/chromatogramList/chromatogram/binaryDataArrayList/binaryDataArray/cvParam/@accession", validator);
        os << "\t\t\t\t\t\t<binary>" << encoded_string << "</binary>\n";
        os << "\t\t\t\t\t</binaryDataArray>\n";
      }
      os << "\t\t\t\t</binaryDataArrayList>\n";
      os << "\t\t\t</chromatogram>" << "\n";
    }

} // namespace OpenMS   // namespace Internal
