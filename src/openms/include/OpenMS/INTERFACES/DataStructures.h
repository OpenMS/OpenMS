// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Authors: Darren Kessner, Hannes Roest, Witold Wolski$
// --------------------------------------------------------------------------

#pragma once

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>

#include <OpenMS/config.h>

namespace OpenMS
{
namespace Interfaces
{

  /**
    @brief The datastructures used by the OpenSwath interfaces

    Many of them are closely related to Proteowizard data structures,
    originally written by Darren Kessner and released under the Apache 2.0
    licence and can be found in the file pwiz/data/msdata/MSData.hpp.

    Original author: Darren Kessner <darren@proteowizard.org>

    Copyright 2007 Spielberg Family Center for Applied Proteomics
      Cedars-Sinai Medical Center, Los Angeles, California  90048


    The following datastructures are used :
    - BinaryDataArray : a struct that holds a std::vector<double> with the data
    - ChromatogramMeta : meta information of a chromatogram (index)
    - Chromatogram : chromatogram data. Contains a vector of pointers to BinaryDataArray,
                     the first one is time array (RT), the second one is intensity
    - SpectrumMeta : meta information of a spectrum (index, identifier, RT, ms_level)
    - Spectrum :     spectrum data. Contains a vector of pointers to BinaryDataArray,
                     the first one is mz array, the second one is intensity
  */

  /// The structure into which encoded binary data goes.
  struct OPENMS_DLLAPI BinaryDataArray
  {
    /// this optional attribute may reference the 'id' attribute of the appropriate dataProcessing.
    //DataProcessingPtr dataProcessingPtr;

    /// the binary data.
    std::vector<double> data;
  };
  typedef boost::shared_ptr<BinaryDataArray> BinaryDataArrayPtr;

  /// Identifying information for a chromatogram
  struct OPENMS_DLLAPI ChromatogramMeta
  {
    /// the zero-based, consecutive index of the chromatogram in the ChromatogramList.
    std::size_t index;
    /// a unique identifier for this chromatogram.
    std::string id;
    /// precursor and product m/z
    double precursor_isolation_target;
    double product_isolation_target;

    ChromatogramMeta() :
      index()
    {
    }

  };
  typedef boost::shared_ptr<ChromatogramMeta> ChromatogramMetaPtr;

  /// A single chromatogram.
  struct OPENMS_DLLAPI Chromatogram
  {
    /// default length of binary data arrays contained in this element.
    std::size_t defaultArrayLength;

private:
    /// list of binary data arrays.
    std::vector<BinaryDataArrayPtr> binaryDataArrayPtrs;

public:
    Chromatogram() :
      defaultArrayLength(2),
      binaryDataArrayPtrs(defaultArrayLength)
    {
      initvec();
    }

private:

    void initvec()
    {
      for (std::size_t i = 0; i < defaultArrayLength; ++i)
      {
        BinaryDataArrayPtr empty(new BinaryDataArray);
        binaryDataArrayPtrs[i] = empty;
      }
    }

public:
    /// get time array (may be null)
    BinaryDataArrayPtr getTimeArray() const
    {
      return binaryDataArrayPtrs[0];
    }

    /// set time array
    void setTimeArray(BinaryDataArrayPtr data)
    {
      binaryDataArrayPtrs[0] = data;
    }

    /// get intensity array (may be null)
    BinaryDataArrayPtr getIntensityArray() const
    {
      return binaryDataArrayPtrs[1];
    }

    /// set intensity array
    void setIntensityArray(BinaryDataArrayPtr data)
    {
      binaryDataArrayPtrs[1] = data;
    }

  };
  typedef boost::shared_ptr<Chromatogram> ChromatogramPtr;

  /// Identifying information for a spectrum
  struct OPENMS_DLLAPI SpectrumMeta
  {
    /// the zero-based, consecutive index of the spectrum in the SpectrumList.
    size_t index;

    /// a unique identifier for this spectrum.
    std::string id;

    /// retention time information
    double RT;

    /// ms level
    int ms_level;

    SpectrumMeta() :
      index(0)
    {
    }

  };
  typedef boost::shared_ptr<SpectrumMeta> SpectrumMetaPtr;

  /// The structure that captures the generation of a peak list (including the underlying acquisitions)
  struct OPENMS_DLLAPI Spectrum
  {
    /// default length of binary data arrays contained in this element.
    std::size_t defaultArrayLength;

private:
    /// list of binary data arrays.
    std::vector<BinaryDataArrayPtr> binaryDataArrayPtrs;

public:
    Spectrum() :
      defaultArrayLength(2),
      binaryDataArrayPtrs(defaultArrayLength)
    {
      initvec();
    }

private:

    void initvec()
    {
      for (std::size_t i = 0; i < defaultArrayLength; ++i)
      {
        BinaryDataArrayPtr empty(new BinaryDataArray);
        binaryDataArrayPtrs[i] = empty;
      }
    }

public:
    /// get m/z array (may be null)
    BinaryDataArrayPtr getMZArray() const
    {
      return binaryDataArrayPtrs[0];
    }

    /// set mz array
    void setMZArray(BinaryDataArrayPtr data)
    {
      binaryDataArrayPtrs[0] = data;
    }

    /// get intensity array (may be null)
    BinaryDataArrayPtr getIntensityArray() const
    {
      return binaryDataArrayPtrs[1];
    }

    /// set intensity array
    void setIntensityArray(BinaryDataArrayPtr data)
    {
      binaryDataArrayPtrs[1] = data;
    }

  };
  typedef boost::shared_ptr<Spectrum> SpectrumPtr;

} //end namespace Interfaces
} //end namespace OpenMS
