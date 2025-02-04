// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MSNumpressCoder.h>

#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/FORMAT/MSNUMPRESS/MSNumpress.h>
#include <boost/math/special_functions/fpclassify.hpp> // std::isfinite
// #define NUMPRESS_DEBUG

#include <iostream>

namespace OpenMS
{
  const std::string MSNumpressCoder::NamesOfNumpressCompression[] = {"none", "linear", "pic", "slof"};

  using namespace ms; // numpress namespace

  void MSNumpressCoder::encodeNP(const std::vector<double> & in, String & result,
      bool zlib_compression, const NumpressConfig & config)
  {
    result.clear();
    encodeNPRaw(in, result, config);
    if (result.empty())
    {
      return;
    }

    // Encode in base64 and compress
    std::vector<String> tmp;
    tmp.push_back(result);
    Base64::encodeStrings(tmp, result, zlib_compression, false);
  }

  void MSNumpressCoder::encodeNP(const std::vector<float> & in, String & result,
      bool zlib_compression, const NumpressConfig & config)
  {
    std::vector<double> dvector(in.begin(), in.end());
    encodeNP(dvector, result, zlib_compression, config);
  }

  void MSNumpressCoder::decodeNP(const String & in, std::vector<double> & out,
      bool zlib_compression, const NumpressConfig & config)
  {
    QByteArray base64_uncompressed;
    Base64::decodeSingleString(in, base64_uncompressed, zlib_compression);

    // Create a temporary string (*not* null-terminated) to hold the data
    std::string tmpstring(base64_uncompressed.constData(), base64_uncompressed.size());
    decodeNPRaw(tmpstring, out, config);

    // NOTE: it is possible (and likely faster) to call directly the const
    // unsigned char * function but this would make it necessary to deal with
    // reinterpret_cast ugliness here ... 
    //
    // decodeNP_internal_(reinterpret_cast<const unsigned char*>(base64_uncompressed.constData()), base64_uncompressed.size(), out, config);
  }

  void MSNumpressCoder::encodeNPRaw(const std::vector<double>& in, String& result, const NumpressConfig & config)
  {
    if (in.empty())
    {
      return;
    }
    if (config.np_compression == NONE)
    {
      return;
    }
    Size dataSize = in.size();

    // using MSNumpress, from johan.teleman@immun.lth.se
    std::vector<unsigned char> numpressed;

    double fixedPoint  = config.numpressFixedPoint;

    try
    {
      size_t byteCount = 0;

      // 1. Resize the data
      switch (config.np_compression)
      {
      case LINEAR:
        numpressed.resize(dataSize * sizeof(std::vector<double>::value_type) + 8);
        break;

      case PIC:
        numpressed.resize(dataSize * sizeof(std::vector<double>::value_type));
        break;

      case SLOF:
      {
        numpressed.resize(dataSize * 2 + 8);
        break;
      }

      case NONE:
        break;

      default:
        break;
      }

      // 2. Convert the data
      std::vector<double> unpressed; // for checking excessive accuracy loss
      switch (config.np_compression)
      {
      case LINEAR:
      {
        if (config.estimate_fixed_point)
        {
          // estimate fixed point either by mass accuracy or by using maximal permissible value
          if (config.linear_fp_mass_acc > 0)
          {
            fixedPoint = numpress::MSNumpress::optimalLinearFixedPointMass(&in[0], dataSize, config.linear_fp_mass_acc);
            // catch failure
            if (fixedPoint < 0.0)
            {
              fixedPoint = numpress::MSNumpress::optimalLinearFixedPoint(&in[0], dataSize);
            }
          }
          else
          {
            fixedPoint = numpress::MSNumpress::optimalLinearFixedPoint(&in[0], dataSize);
          }
        }
        byteCount = numpress::MSNumpress::encodeLinear(&in[0], dataSize, &numpressed[0], fixedPoint);
        numpressed.resize(byteCount);
        if (config.numpressErrorTolerance > 0.0)   // decompress to check accuracy loss
        {
          numpress::MSNumpress::decodeLinear(numpressed, unpressed);
        }
        break;
      }

      case PIC:
      {
        byteCount = numpress::MSNumpress::encodePic(&in[0], dataSize, &numpressed[0]);
        numpressed.resize(byteCount);
        if (config.numpressErrorTolerance > 0.0)   // decompress to check accuracy loss
        {
          numpress::MSNumpress::decodePic(numpressed, unpressed);
        }
        break;
      }

      case SLOF:
      {
        if (config.estimate_fixed_point)
        {
          fixedPoint = numpress::MSNumpress::optimalSlofFixedPoint(&in[0], dataSize);
        }
        byteCount = numpress::MSNumpress::encodeSlof(&in[0], dataSize, &numpressed[0], fixedPoint);
        numpressed.resize(byteCount);
        if (config.numpressErrorTolerance > 0.0)   // decompress to check accuracy loss
        {
          numpress::MSNumpress::decodeSlof(numpressed, unpressed);
        }
        break;
      }

      case NONE:
      {
        // already checked above ...
        break;
      }

      default:
        break;
      }

#ifdef NUMPRESS_DEBUG
      std::cout << "encodeNPRaw: numpressed array with with length " << numpressed.size() << std::endl;
      for (int i = 0; i < byteCount; i++)
      {
        std::cout << "array[" << i << "] : " << (int)numpressed[i] << std::endl;
      }
#endif


      // 3. Now check to see if encoding introduces excessive error
      int n = -1;
      if (config.numpressErrorTolerance > 0.0)
      {
        if (PIC == config.np_compression) // integer rounding, abs accuracy is +- 0.5
        {
          for (n = static_cast<int>(dataSize)-1; n>=0; n-- ) // check for overflow, strange rounding
          {
            if ((!std::isfinite(unpressed[n])) || (fabs(in[n] - unpressed[n]) >= 1.0))
            {
              break;
            }
          }
        }
        else // check for tolerance as well as overflow
        {
          for (n=static_cast<int>(dataSize)-1; n>=0; n--)
          {
            double u = unpressed[n];
            double d = in[n];
            if (!std::isfinite(u) || !std::isfinite(d))
            {
#ifdef NUMPRESS_DEBUG
              std::cout << "infinite u: " << u << " d: " << d << std::endl;
#endif
              break;
            }
            if (!d)
            {
              if (fabs(u) > config.numpressErrorTolerance)
              {
#ifdef NUMPRESS_DEBUG
                std::cout << "fabs(u): " << fabs(u) << " > config.numpressErrorTolerance: " << config.numpressErrorTolerance << std::endl;
#endif
                break;
              }
            }
            else if (!u)
            {
              if (fabs(d) > config.numpressErrorTolerance)
              {
#ifdef NUMPRESS_DEBUG
                std::cout << "fabs(d): " << fabs(d) << " > config.numpressErrorTolerance: " << config.numpressErrorTolerance << std::endl;
#endif
                break;
              }
            }
            else if (fabs(1.0 - (d / u)) > config.numpressErrorTolerance)
            {
#ifdef NUMPRESS_DEBUG
              std::cout << "d: " << d << " u: " << u << std::endl;
              std::cout << "fabs(1.0 - (d / u)): " << fabs(1.0 - (d / u)) << " > config.numpressErrorTolerance: " << config.numpressErrorTolerance << std::endl;
#endif
              break;
            }
          }
        }
      }
      if (n >= 0)
      {
        //Comment: throw?
        std::cerr << "Error occurred at position n = " << n << ". Enable NUMPRESS_DEBUG to get more info." << std::endl;
      }
      else
      {
        result = String(std::string(reinterpret_cast<const char*>(&numpressed[0]), byteCount));
        // Other solution:
        // http://stackoverflow.com/questions/2840835/way-to-get-unsigned-char-into-a-stdstring-without-reinterpret-cast
        // result = String( std::string(&numpressed[0], &numpressed[0] + byteCount) );
      }
    }
    catch (int e)
    {
      std::cerr << "MSNumpress encoder threw exception: " << e << std::endl;
    }
    catch (char const * e)
    {
      std::cerr << "MSNumpress encoder threw exception: " << e << std::endl;
    }
    catch (...)
    {
      std::cerr << "Unknown exception while encoding " << dataSize << " doubles" << std::endl;
    }
  }

  void MSNumpressCoder::decodeNPRaw(const std::string & in, std::vector<double>& out, const NumpressConfig & config)
  {
    decodeNPInternal_(reinterpret_cast<const unsigned char*>(in.c_str()), in.size(), out, config);
  }

  void MSNumpressCoder::decodeNPInternal_(const unsigned char* in, size_t in_size, std::vector<double>& out, const NumpressConfig & config)
  {
    out.clear();
    if (in_size == 0) return;

    size_t byteCount = in_size;

#ifdef NUMPRESS_DEBUG
    std::cout << "decodeNPInternal_: array input with length " << in_size << std::endl;
    for (int i = 0; i < in_size; i++)
    {
      std::cout << "array[" << i << "] : " << (int)in[i] << std::endl;
    }
#endif

    try
    {
      size_t initialSize;

      switch (config.np_compression)
      {
      case LINEAR:
      {
        initialSize = byteCount * 2;
        if (out.size() < initialSize)
        { 
          out.resize(initialSize);
        }
        size_t count = numpress::MSNumpress::decodeLinear(in, byteCount, &out[0]);
        out.resize(count);
        break;
      }

      case PIC:
      {
        initialSize = byteCount * 2;
        if (out.size() < initialSize)
        { 
          out.resize(initialSize);
        }
        size_t count = numpress::MSNumpress::decodePic(in, byteCount, &out[0]);
        out.resize(count);
        break;
      }

      case SLOF:
      {
        initialSize = byteCount / 2;
        if (out.size() < initialSize)
        { 
          out.resize(initialSize);
        }
        size_t count = numpress::MSNumpress::decodeSlof(in, byteCount, &out[0]);
        out.resize(count);
        break;
      }

      case NONE:
      {
        return;
      }

      default:
        break;
      }

    }
    catch (...)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error in Numpress decompression");
    }

#ifdef NUMPRESS_DEBUG
    std::cout << "decodeNPInternal_: output size " << out.size() << std::endl;
    for (int i = 0; i < out.size(); i++)
    {
      std::cout << "array[" << i << "] : " << out[i] << std::endl;
    }
#endif


  }

} //namespace OpenMS
