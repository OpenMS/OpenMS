// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

int main(int argc, const char** argv)
{
  if (argc < 2) return 1;
  // the path to the data should be given on the command line
  String tutorial_data_path(argv[1]);

  // temporary data storage
  PeakMap map;

  // convert MzXML to MzML. Internally we use FileHandler to do the actual work.
  // Here we limit the input type to be MZXML only
  FileHandler().loadExperiment(tutorial_data_path + "/data/Tutorial_FileIO.mzXML", map, {FileTypes::MZXML});
  FileHandler().storeExperiment("Tutorial_FileIO.mzML", map, {FileTypes::MZML});

  // The FileHandler object can also hold options for how to load the file
  FileHandler f = FileHandler();
  PeakFileOptions opts = PeakFileOptions();
  // Here we set the MZ range to load to 100-200
  opts.setMZRange( {100, 200} );
  f.setOptions(opts);
  f.loadExperiment(tutorial_data_path + "/data/Tutorial_FileIO.mzXML", map, {FileTypes::MZXML});


  // we can also load an experiment from a file without any restrictions on the file type:
  FileHandler().loadExperiment(tutorial_data_path + "/data/Tutorial_Spectrum1D.dta", map);

  // if we want to allow all types that can store MS2 data we can do the following:
  FileHandler().loadExperiment(tutorial_data_path + "/data/Tutorial_FileIO.mzXML", map, FileTypeList::typesWithProperties({FileTypes::FileProperties::PROVIDES_EXPERIMENT}));
  // The curly braces can contain multiple file properties. The FileTypeList that is created is the intersection of these properties
  // so: FileTypeList::typesWithProperties({FileTypes::FileProperties::PROVIDES_EXPERIMENT, FileTypes::FileProperties::READABLE})
  // returns only fileTypes which can store both MS1 and MS2 spectra

  // We use various FileHandler functions to load other types.
  FeatureMap feat;
  FileHandler().loadFeatures(tutorial_data_path + "/data/Tutorial_Labeled.featureXML", feat);

  // If we try to load something from a file that can't store that info (for example trying to get an experiment from an idXML file)
  // An error gets thrown at run time. Check out @p FileHandler class for more info

  return 0;
} //end of main
