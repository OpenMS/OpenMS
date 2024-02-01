// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once
#include <OpenMS/config.h>
#include <OpenMS/FORMAT/MzMLbSeekableDevice.h>
#include <OpenMS/FORMAT/HANDLERS/MzMLHandlerHelper.h>

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/concepts.hpp>  // seekable_device

namespace OpenMS
{
  namespace HDF5
  {
    // custom binary data array loader for HDF5 data (as opposed to the base64 extraction from XML)
    class OPENMS_DLLAPI MzMLbBinaryDataArrayLoader
  {
      using MzMLbInputStream = boost::iostreams::stream<MzMLbSeekableDevice>;
      MzMLbInputStream is_; // input stream from seekable device

    public:
      explicit MzMLbBinaryDataArrayLoader(const MzMLbSeekableDevice& device) : 
        is_(device)
        {          
        }
      MzMLbBinaryDataArrayLoader() = default;
      MzMLbBinaryDataArrayLoader(MzMLbBinaryDataArrayLoader&& other) = default;
      MzMLbBinaryDataArrayLoader(const MzMLbBinaryDataArrayLoader& other) = default;
      ~MzMLbBinaryDataArrayLoader() = default;

    private:
      // see https://github.com/ProteoWizard/pwiz/blob/072abd4c764157e922e1c6cb0cecf94166993d30/pwiz/data/msdata/IO.cpp#L2277-L2278
      enum class PredictionType
      {
        None,
        Delta,
        Linear
      };
  
      // only on float data arrays
      static void predict_(OpenMS::Internal::MzMLHandlerHelper::BinaryData& bda, PredictionType pred)
      {    
        using OpenMS::Internal::MzMLHandlerHelper;

        if (bda.data_type != OpenMS::Internal::MzMLHandlerHelper::BinaryData::DT_FLOAT) return; // TODO: check how pwiz prevents this

        switch (pred)
        {
          case PredictionType::Delta:
            switch (bda.precision)
            {
              case MzMLHandlerHelper::BinaryData::PRE_32:
                for (size_t i = 2; i < bda.floats_32.size(); i++) 
                  bda.floats_32[i] = bda.floats_32[i] + bda.floats_32[i - 1] - bda.floats_32[0];
                break;
              case MzMLHandlerHelper::BinaryData::PRE_64:
                  for (size_t i = 2; i < bda.floats_64.size(); i++)
                    bda.floats_64[i] = bda.floats_64[i] + bda.floats_64[i - 1] - bda.floats_64[0];
                  break;
              default:
                // TODO: error
              break;
            }
            break;
          case PredictionType::Linear:
            switch (bda.precision)
            {
                case MzMLHandlerHelper::BinaryData::PRE_32:
                    for (size_t i = 2; i < bda.floats_32.size(); i++) 
                      bda.floats_32[i] = bda.floats_32[i] + 2.0f * bda.floats_32[i - 1] - bda.floats_32[i - 2] - bda.floats_32[1];
                    break;
                case MzMLHandlerHelper::BinaryData::PRE_64:
                    for (size_t i = 2; i < bda.floats_64.size(); i++)
                      bda.floats_64[i] = bda.floats_64[i] + 2.0 * bda.floats_64[i - 1] - bda.floats_64[i - 2] - bda.floats_64[1];
                    break;
                default:
                  // TODO: error
                break;
            }
            break;
          case PredictionType::None:
            break;
        }
      }


    // read HDF5 dataset referenced in XML part (with external_offset) into target. Target could be the m/z or intensity dimension of a spectrum, 
    // int. or rt dim. of a chromatogram, or an openms data array.
    // The offset is needed to find the actual data in the HDF5 dataset item. mzMLb allows to store blocks of data for better compression.
    static void readMzMLbBinaryDataArray_(MzMLbInputStream& is, OpenMS::Internal::MzMLHandlerHelper::BinaryData& target)
    {
      const std::string& external_dataset = target.mzMLb_dataset;
      const size_t external_offset = target.mzMLb_offset;
      const size_t external_array_length = target.mzMLb_array_length;
      size_t array_length = external_array_length != 0 ? external_array_length : 0;

      using BD = OpenMS::Internal::MzMLHandlerHelper::BinaryData;
      auto precisionToString = [](auto dt) -> std::string {
          switch (dt) {
              case BD::PRE_NONE:
                  return "PRE_NONE";
              case BD::PRE_32:
                  return "PRE_32";
              case BD::PRE_64:
                  return "PRE_64";
              default:
                  return "Unknown DataType";
          }
      };

      //PredictionType pred;
      auto dataTypeToString = [](auto dt) -> std::string {
          switch (dt) {
              case BD::DT_NONE:
                  return "DT_NONE";
              case BD::DT_FLOAT:
                  return "DT_FLOAT";
              case BD::DT_INT:
                  return "DT_INT";
              case BD::DT_STRING:
                  return "DT_STRING";
              default:
                  return "Unknown DataType";
          }
      };

      std::string data_type_string = dataTypeToString(target.data_type);
      std::string precision_string = precisionToString(target.precision);
      
      
      std::cout << "Reading dataset/offset/length/type/precision: " 
        << external_dataset << " " 
        << external_offset << " " 
        << array_length  << " "
        << data_type_string << " "        
        << precision_string << " "
        << std::endl;
      
/*
 TODO: what is this?                   
      // primary array types so set the default array length
      if (binaryDataArray->hasCVParam(MS_m_z_array) ||
          binaryDataArray->hasCVParam(MS_time_array) ||
          binaryDataArray->hasCVParam(MS_intensity_array) ||
          binaryDataArray->hasCVParam(MS_wavelength_array))
      {
          if (defaultArrayLength)
              *defaultArrayLength = arrayLength_;
      } */

      if (!external_dataset.empty())
      {
        // jump to start of data we want to extract
        is->seek(external_dataset, external_offset, std::ios_base::beg);

/*
        // MSNumpress? then we extract raw bytes and decode them
        if (target.numpress != BinaryDataEncoder::Numpress_None)
        {
            vector<char> buf(encodedLength_);
            (*mzMLb_is)->read_opaque(external_dataset_, &buf[0], encodedLength_);
            config.format = BinaryDataEncoder::Format_MzMLb;
            BinaryDataEncoder encoder(config);
            encoder.decode(&buf[0], buf.size(), binaryDataArray->data);
        }
        else */       
        {
          // load the binary data at the given offset into the target
          if (external_array_length > 0)
          {
            // TODO: basically implement all conrner cases  void MzMLHandlerHelper::decodeBase64Arrays(std::vector<BinaryData>& data, const bool skipXMLCheck) without the base64 part
            // TODO: this is currently not thread safe. works with export OMP_NUM_THREADS=1. maybe just put a critical section here but this would prevent parallel parsing of different files

            if (target.precision == BD::PRE_64 && target.data_type == BD::DT_FLOAT)
            {
              target.floats_64.resize(external_array_length);
              is->read(external_dataset, &target.floats_64[0], external_array_length);
            }
            /* probably not supported because of MSNumpress always being 64 bit
            else if (target.precision == BD::PRE_32 && target.data_type == BD::DT_FLOAT)
            {
              target.floats_32.resize(external_array_length);
              is->read(external_dataset, &target.floats_32[0], external_array_length);
            }
            */
            else if (target.precision == BD::PRE_64 && target.data_type == BD::DT_INT)
            {
              target.ints_64.resize(external_array_length);
              is->read(external_dataset, &target.ints_64[0], external_array_length);
            }
            /*
            else if (target.precision == BD::PRE_32 && target.data_type == BD::DT_INT)
            {
              target.ints_32.resize(external_array_length);
              is->read(external_dataset, &target.ints_32[0], external_array_length);
            }
            */
// TODO: DT_STRING?
          }
        }
        PredictionType pt;
        predict_(target, pt);
      }
    
  } 


   public:
      // input_data BinaryData objects contain the HDF5 dataset, the offset and the array length
      // as well as precision etc.
      // use this information to extract the actual binary data from the HDF5
      void fill(std::vector<OpenMS::Internal::MzMLHandlerHelper::BinaryData>& input_data)
      {
        for (auto& bda : input_data)
        {
          readMzMLbBinaryDataArray_(is_, bda);          
        }        
      }

  };
  }
}
