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

#include <OpenMS/FORMAT/MSNumpressCoder.h>

#include <OpenMS/MATH/MISC/MSNumpress.h>
#include <boost/math/special_functions/fpclassify.hpp> // boost::math::isfinite

// #define NUMPRESS_DEBUG

namespace OpenMS
{
  using namespace ms; // numpress namespace

  void MSNumpressCoder::encodeNP_(std::vector<double> & in, String & result, NumpressConfig config)
  {
    if (in.empty()) return;
    if (config.np_compression == NONE) return;

    Size dataSize = in.size();

    // using MSNumpress, from johan.teleman@immun.lth.se
    size_t byteCount;
    std::vector<unsigned char> compressed;
    std::vector<unsigned char> numpressed;
    std::vector<float> data32;
    std::vector<double> data64endianized;

    double fixedPoint  = config.numpressFixedPoint;

    try 
    {
      // 1. Resize the data
      switch(config.np_compression)
      {
        case LINEAR:
          numpressed.resize(dataSize * sizeof(std::vector<double>::value_type) + 8);
          break;
        case PIC:
          numpressed.resize(dataSize * sizeof(std::vector<double>::value_type));
          break;
        case SLOF:
          numpressed.resize(dataSize * 2 + 8);
          break;
        case NONE:
          break;
      }

      // 2. Convert the data
      std::vector<double> unpressed; // for checking excessive accurary loss
      switch (config.np_compression) 
      {
        case LINEAR:
        {
          if (config.estimate_fixed_point) {fixedPoint = numpress::MSNumpress::optimalLinearFixedPoint(&in[0], dataSize);}
          byteCount = numpress::MSNumpress::encodeLinear(&in[0], dataSize, &numpressed[0], fixedPoint);
          if (config.numpressErrorTolerance) // decompress to check accuracy loss
          {
            numpress::MSNumpress::decodeLinear(numpressed,unpressed); 
          }
          break;
        }

        case PIC:
        {
          byteCount = numpress::MSNumpress::encodePic(&in[0], dataSize, &numpressed[0]);
          if (config.numpressErrorTolerance) // decompress to check accuracy loss
          {
            numpress::MSNumpress::decodePic(numpressed,unpressed); 
          }
          break;
        }

        case SLOF:
        {
          if (config.estimate_fixed_point) {fixedPoint = numpress::MSNumpress::optimalSlofFixedPoint(&in[0], dataSize);}
          byteCount = numpress::MSNumpress::encodeSlof(&in[0], dataSize, &numpressed[0], fixedPoint);
          if (config.numpressErrorTolerance) // decompress to check accuracy loss
          {
            numpress::MSNumpress::decodeSlof(numpressed,unpressed); 
          }
          break;
        }

        case NONE:
        {
          // already checked above ... 
          break;
        }
      }

#ifdef NUMPRESS_DEBUG
          std::cout << "encodeNP_: numpressed array with with length " << numpressed.size() << std::endl;
          for (int i = 0; i < numpressed.size(); i++)
          {
            std::cout << "array[" << i << "] : " << (int)numpressed[i] << std::endl;
          }
#endif

      // 3. Now check to see if encoding introduces excessive error
      int n=-1;
      if (config.numpressErrorTolerance) 
      {
        if (PIC == config.np_compression)  // integer rounding, abs accuracy is +- 0.5
        {
          for (n=(int)dataSize;n--;) // check for overflow, strange rounding
          {
            if ((!boost::math::isfinite(unpressed[n])) || (fabs(in[n]-unpressed[n])>=1.0)) 
            {
              break;
            }
          }
        }
        else // check for tolerance as well as overflow
        {
          for (n=(int)dataSize;n--;)
          {
            double d,u;
            if (!boost::math::isfinite(u = unpressed[n])||!boost::math::isfinite(d = in[n]))
            {
              break;
            }
            if (!d)
            {
              if (fabs(u) > config.numpressErrorTolerance)
              {
                break;
              }
            }
            else if (!u)
            {
              if (fabs(d) > config.numpressErrorTolerance)
              {
                break;
              }
            }
            else if (fabs(1.0-(d/u)) > config.numpressErrorTolerance)
            {
              break;
            }
          }
        }
      }
      if (n>=0)
      {
#ifdef NUMPRESS_DEBUG
        std::cout << "Error occured with n = " << n <<std::endl;
#endif
      }
      else
      {
        result = String(std::string(reinterpret_cast<const char *>(&numpressed[0]), byteCount));
        // Other solution:
        // http://stackoverflow.com/questions/2840835/way-to-get-unsigned-char-into-a-stdstring-without-reinterpret-cast
        // result = String( std::string(&numpressed[0], &numpressed[0] + byteCount) );
      }
    } catch (int e) 
    {
      std::cerr << "MZNumpress encoder threw exception: " << e << std::endl;
    } catch (...) 
    {
      std::cerr << "Unknown exception while encoding " << dataSize << " doubles" << std::endl;
    }
  }

  void MSNumpressCoder::decodeNP_(std::string in, std::vector<double> & out, NumpressConfig config)
  {
    decodeNP_internal_(reinterpret_cast<const unsigned char*>(in.c_str()), in.size(), out, config);
  }

  void MSNumpressCoder::decodeNP_internal_(const unsigned char* in, size_t in_size, std::vector<double> & out, NumpressConfig config)
  {
    out.clear();
    if (in_size == 0) return;

    size_t byteCount = in_size;
    size_t initialSize;

#ifdef NUMPRESS_DEBUG
    std::cout << "decodeNP_: array input with length " << in_size << std::endl;
    for (int i = 0; i < in_size; i++)
    {
      std::cout << "array[" << i << "] : " << (int)in[i] << std::endl;
    }
#endif

    try 
    {
      switch(config.np_compression)
      {
        case LINEAR:
        {
          initialSize = byteCount * 2;
          if (out.size() < initialSize) { out.resize(initialSize); }
          size_t count = numpress::MSNumpress::decodeLinear(in, byteCount, &out[0]);
          out.resize(count);
          break;
        }
        case PIC:
        {
          initialSize = byteCount * 2;
          if (out.size() < initialSize) { out.resize(initialSize); }
          size_t count = numpress::MSNumpress::decodePic(in, byteCount, &out[0]);
          out.resize(count);
          break;
        }
        case SLOF:
        {
          initialSize = byteCount / 2;
          if (out.size() < initialSize) { out.resize(initialSize); }
          size_t count = numpress::MSNumpress::decodeSlof(in, byteCount, &out[0]);
          out.resize(count);
          break;
        }
        case NONE:
        {
          return;
          break;
        }
      }

    } catch (...) 
    {
      throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Error in Numpress decompression");
    }
  }

} //namespace OpenMS

