// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#ifdef WITH_HDF5

///////////////////////////

///////////////////////////

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLbFile.h>

using namespace std;

#include <string>
#include <iostream>

START_TEST(MzMLb, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

//using mzMLbInputStream = boost::iostreams::stream<OpenMS::HDF5::MzMLbSeekableDevice>;


/*
  // fill the empty spectrum with HDF5 dataset data
  void fillSpectrum(mzMLbInputStream& is, size_t i, MzMLbMapping mz_mapping, MzMLbMapping int_mapping)
  {
    // TODO: in the XML blob these metavalues are present but we probably don't store them for m/z and int dim separately - need to check if with some prefix...
    // e.g., in the XML we have two entries - one for each dimension:
    
    // <binaryDataArray encodedLength="0">
    //           <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
    //           <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
    //           <cvParam cvRef="MS" accession="MS:1002841" name="external HDF5 dataset" value="spectrum_MS_1000514_double"/>
    //           <cvParam cvRef="MS" accession="MS:1002843" name="external array length" value="636"/>
    //           <cvParam cvRef="MS" accession="MS:1002842" name="external offset" value="304577"/>
    //           <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value="" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
    //           <binary/>
    //         </binaryDataArray>
    //         <binaryDataArray encodedLength="0">
    //           <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
    //           <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
    //           <cvParam cvRef="MS" accession="MS:1002841" name="external HDF5 dataset" value="spectrum_MS_1000515_double"/>
    //           <cvParam cvRef="MS" accession="MS:1002843" name="external array length" value="636"/>
    //           <cvParam cvRef="MS" accession="MS:1002842" name="external offset" value="304577"/>
    //           <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitCvRef="MS" unitAccession="MS:1000131" unitName="number of detector counts"/>
    //           <binary/>
    //       </binaryDataArray>

    //       Easiest is probably to keep track of these meta values during XML part parsing and build a map.
    

    const std::string mz_external_dataset = mz_mapping.external_dataset;
    const size_t mz_external_array_length = mz_mapping.external_array_length;
    const size_t mz_external_offset = mz_mapping.external_offset;
    const size_t mz_precision_32 = mz_mapping.precision_32;

    const std::string int_external_dataset = int_mapping.external_dataset;
    const size_t int_external_array_length = int_mapping.external_array_length;
    const size_t int_external_offset = int_mapping.external_offset;
    const size_t int_precision_32 = int_mapping.precision_32;

    // check if we can use populateSpectrumWithDate or similar... we basically need to
    // retrieve the data, put it into an OpenMS BinaryData object with annotated information (precision etc.)
    // and then do e.g. numpress decoding (optional) and predict on the decoded data with our existing functions

    // depends on precision...
    FLOAT64_DATA_ARRAY mz_array;
    mz_array.reserve(mz_external_array_length);

    FLOAT64_DATA_ARRAY int_array;
    int_array.reserve(int_external_array_length);

    readMzMLbBinaryDataArray<FLOAT64_DATA_ARRAY>(is, external_dataset, external_array_length, external_offset, mz_array)
    readMzMLbBinaryDataArray<FLOAT64_DATA_ARRAY>(is, external_dataset, external_array_length, external_offset, int_array)

    // build spectrum
    spectrum.reserve(external_array_length);
    for (Size j = 0; j < external_array_length; j++)
    {
      spectrum.emplace_back(mz_array[j], int_array[j]);
    }

    // TODO: check if there are additional data arrays and do the same
    
    //for (Size j = 2; j < data.size(); j++)
    //{
//      spectrum.getFloatDataArrays().push_back(MSSpectrum::FloatDataArray());
  //    spectrum.getFloatDataArrays().back().reserve(data[j]->data.size());
    //  spectrum.getFloatDataArrays().back().setName(data[j]->description);
    //  for (const auto& k : data[j]->data) spectrum.getFloatDataArrays().back().push_back(k);
    //}
  }

  void fillChromatogram(ChromatogramType& chromatogram, std::ifstream& ifs)
  {
    // TODO:...
  }
*/


START_SECTION((MzMLb()))
{
  const std::string filename( OPENMS_GET_TEST_DATA_PATH("msconvert.0.24017-6a003b2.mzMLb") ); // file converted with pwiz
  auto mzmlb = MzMLbFile();
  MSExperiment exp = mzmlb.load(filename);
  for (auto s : exp.getSpectrum(0))
  {
    std::cout << "mz: " << s.getMZ() << " int: " << s.getIntensity() << std::endl;
  }
}
END_SECTION

END_TEST

#endif