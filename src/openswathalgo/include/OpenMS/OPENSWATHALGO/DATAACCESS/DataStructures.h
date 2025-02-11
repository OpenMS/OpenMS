// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Darren Kessner, Hannes Roest, Witold Wolski$
// --------------------------------------------------------------------------

#pragma once

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>

#include <OpenMS/OPENSWATHALGO/OpenSwathAlgoConfig.h>

namespace OpenSwath
{
  /**
    @brief The datastructures used by the OpenSwath interfaces

    Many of them are closely related to Proteowizard data structures,
    originally written by Darren Kessner and released under the Apache 2.0 licence.

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
  struct OPENSWATHALGO_DLLAPI OSBinaryDataArray
  {
    /// this optional attribute may reference the 'id' attribute of the appropriate dataProcessing.
    //DataProcessingPtr dataProcessingPtr;

    /// the binary data.
    std::vector<double> data;

    /// (optional) data description for non-standard arrays.
    std::string description;
  };
  typedef OSBinaryDataArray BinaryDataArray;
  typedef boost::shared_ptr<BinaryDataArray> BinaryDataArrayPtr;

  /// Identifying information for a chromatogram
  struct OPENSWATHALGO_DLLAPI OSChromatogramMeta
  {
    /// the zero-based, consecutive index of the chromatogram in the ChromatogramList.
    std::size_t index;
    /// a unique identifier for this chromatogram.
    std::string id;
    OSChromatogramMeta() :
      index()
    {
    }

  };
  typedef OSChromatogramMeta ChromatogramMeta;
  typedef boost::shared_ptr<ChromatogramMeta> ChromatogramMetaPtr;

  /// A single chromatogram.
  struct OPENSWATHALGO_DLLAPI OSChromatogram
  {
private:
    /// default length of binary data arrays contained in this element.
    std::size_t defaultArrayLength;

    /// this attribute can optionally reference the 'id' of the appropriate dataProcessing.
    //DataProcessingPtr dataProcessingPtr;
    /// description of precursor ion information (i.e. Q1 settings)
    //Precursor precursor;
    /// description of product ion information (i.e. Q3 settings)
    //Product product;

    /// list of binary data arrays.
    std::vector<BinaryDataArrayPtr> binaryDataArrayPtrs;
public:

    OSChromatogram() :
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
    BinaryDataArrayPtr getTimeArray()
    {
      return binaryDataArrayPtrs[0];
    }

    /// set time array
    void setTimeArray(BinaryDataArrayPtr data)
    {
      binaryDataArrayPtrs[0] = data;
    }

    /// get intensity array (may be null)
    BinaryDataArrayPtr getIntensityArray()
    {
      return binaryDataArrayPtrs[1];
    }

    /// set intensity array
    void setIntensityArray(BinaryDataArrayPtr data)
    {
      binaryDataArrayPtrs[1] = data;
    }

    /// non-mutable access to the underlying data arrays
    const std::vector<BinaryDataArrayPtr> & getDataArrays() const
    {
      return binaryDataArrayPtrs;
    }

    /// mutable access to the underlying data arrays
    std::vector<BinaryDataArrayPtr> & getDataArrays()
    {
      return binaryDataArrayPtrs;
    }

  };
  typedef OSChromatogram Chromatogram;
  typedef boost::shared_ptr<Chromatogram> ChromatogramPtr;

  /// Identifying information for a spectrum
  struct OPENSWATHALGO_DLLAPI OSSpectrumMeta
  {
    /// the zero-based, consecutive index of the spectrum in the SpectrumList.
    size_t index;

    /// a unique identifier for this spectrum.
    std::string id;

    double RT;

    int ms_level;

    OSSpectrumMeta() :
      index(0)
    {
    }

    ///Comparator for the retention time.
    struct RTLess
    {
      inline bool operator()(const OSSpectrumMeta& a, const OSSpectrumMeta& b) const
      {
        return a.RT < b.RT;
      }
    };

  };
  typedef OSSpectrumMeta SpectrumMeta;
  typedef boost::shared_ptr<SpectrumMeta> SpectrumMetaPtr;

  /// The structure that captures the generation of a peak list (including the underlying acquisitions)
  struct OPENSWATHALGO_DLLAPI OSSpectrum
  {
private:
    /// default length of binary data arrays contained in this element.
    std::size_t defaultArrayLength;

    /// list of binary data arrays.
    std::vector<BinaryDataArrayPtr> binaryDataArrayPtrs;

public:
    OSSpectrum() :
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

    /// set m/z array
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

    void setDriftTimeArray(BinaryDataArrayPtr data)
    {
      data->description = "Ion Mobility";
      binaryDataArrayPtrs.push_back(data);
    }

    /// get drift time array (may be null)
    BinaryDataArrayPtr getDriftTimeArray() const
    {
      // The array name starts with "Ion Mobility", but may carry additional
      // information such as the actual unit in which it was measured (seconds,
      // milliseconds, volt-second per square centimeter). We currently ignore
      // the unit but return the correct array.
      // For diaPASEF data converted with proteowizard ion mobility arrays are stored in "inverse reduced ion mobility"
      for (auto & bda : binaryDataArrayPtrs)
      {
        if (bda->description.find("Ion Mobility") == 0)
        {
          return bda;
        }
        else if (bda->description.find("mean inverse reduced ion mobility array") == 0)
        {
          return bda;
        }
      }
      return BinaryDataArrayPtr(); // return null
    }

    /// non-mutable access to the underlying data arrays
    const std::vector<BinaryDataArrayPtr> & getDataArrays() const
    {
      return binaryDataArrayPtrs;
    }

    /// mutable access to the underlying data arrays
    std::vector<BinaryDataArrayPtr> & getDataArrays()
    {
      return binaryDataArrayPtrs;
    }

  };
  typedef OSSpectrum Spectrum;
  typedef boost::shared_ptr<Spectrum> SpectrumPtr;
} //end Namespace OpenSwath

