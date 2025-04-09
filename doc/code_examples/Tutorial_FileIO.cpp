// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/openms_data_path.h> // exotic header for path to tutorial data
#include <iostream>

using namespace OpenMS;
using namespace std;

int main(int argc, const char** argv)
{
  auto file_mzXML = OPENMS_DOC_PATH + String("/code_examples/data/Tutorial_FileIO.mzXML");

  // temporary data storage
  PeakMap map;

  // convert MzXML to MzML. Internally we use FileHandler to do the actual work.
  // Here we limit the input type to be MzXML only
  FileHandler().loadExperiment(file_mzXML, map, {FileTypes::MZXML});
  FileHandler().storeExperiment("Tutorial_FileIO.mzML", map, {FileTypes::MZML});

  // The FileHandler object can also hold options for how to load the file
  FileHandler f;
  PeakFileOptions opts;
  // Here we set the MZ range to load to 100-200
  opts.setMZRange({100, 200});
  f.setOptions(opts);
  f.loadExperiment(file_mzXML, map, {FileTypes::MZXML});


  // we can also load an experiment from a file without any restrictions on the file type:
  FileHandler().loadExperiment(File::path(file_mzXML) + "/Tutorial_Spectrum1D.dta", map);

  // if we want to allow all types that can store MS2 data we can do the following:
  FileHandler().loadExperiment(file_mzXML, map,
                               FileTypeList::typesWithProperties({FileTypes::FileProperties::PROVIDES_EXPERIMENT}));
  // The curly braces can contain multiple file properties. The FileTypeList that is created is the intersection of these properties
  // so: FileTypeList::typesWithProperties({FileTypes::FileProperties::PROVIDES_EXPERIMENT, FileTypes::FileProperties::READABLE})
  // returns only fileTypes which can store both MS1 and MS2 spectra

  // We use various FileHandler functions to load other types.
  FeatureMap feat;
  FileHandler().loadFeatures(File::path(file_mzXML) + "/Tutorial_Labeled.featureXML", feat);

  // If we try to load something from a file that can't store that info (for example trying to get an experiment from an idXML file)
  // An error gets thrown at run time. Check out @p FileHandler class for more info

} // end of main
