
// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#ifdef WITH_HDF5

///////////////////////////

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <OpenMS/FORMAT/MzMLbSeekableDevice.h>
#include <OpenMS/FORMAT/HANDLERS/MzMLbBinaryDataArrayLoader.h>

#include "H5Cpp.h"
#include "blosc_filter.h"

namespace OpenMS
{
class MzMLbFile
{
  public:
    MzMLbFile()
    {
      // load blosc plugin (could be part of a HDF5 singleton if we use it somewhere else)
      char *version, *date;
      auto return_code = register_blosc(&version, &date);
      //if (return_code < 0) throw ; // TODO
      std::cout << "Blosc version info: " << version << " " << date << std::endl;
    }

    MSExperiment load(const std::string& file_name)
    {
      // open mzMLb file
      auto mzMLb = OpenMS::HDF5::MzMLbSeekableDevice(file_name);
      std::streamsize xml_size = mzMLb.size("mzML");
      std::cout << xml_size << std::endl; // size of XML part?
      
      // Allocate the buffer (plus one for the null terminator)
      std::string xml_buffer(xml_size, '\0');

      // Read the XML blob
      mzMLb.read(&xml_buffer[0], xml_size);
      std::cout << xml_buffer << std::endl;
  
      // Create MSExperiment with all meta data but no peak or chromatogram and binary array data
      MzMLFile mzfile;

      // create experiment from XML buffer. 
      // setting the filename will use the MzMLbBinaryDataArrayLoader to fill spectra and chromatograms from the HDF5
      MSExperiment exp;
      mzfile.loadBuffer(xml_buffer, exp, file_name); //TODO: check if this also works if root element is "indexedMzML" (default: "mzML")
      std::cout << "chromatograms: " << exp.getNrChromatograms() << "\tspectra: " << exp.getNrSpectra() << std::endl;
      return exp; // RVOP
    }
};
}
#endif
