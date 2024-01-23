// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

///////////////////////////

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include "H5Cpp.h"
#include "blosc_filter.h"

using namespace std;

#include <OpenMS/FORMAT/Connection_mzMLb.h>

//using MzMLb = pwiz::msdata::mzmlb::Connection_mzMLb;

#include <string>
#include <iostream>

START_TEST(MzMLb, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

using MzMLb = OpenMS::Connection_mzMLb; // TODO: maybe rename? MzMLbStream? just MzMLb

using mzMLbInputStream = boost::iostreams::stream<MzMLb>;

// see https://github.com/ProteoWizard/pwiz/blob/072abd4c764157e922e1c6cb0cecf94166993d30/pwiz/data/msdata/IO.cpp#L2277-L2278
enum class PredictionType
{
  Prediction_None,
  Prediction_Delta,
  Prediction_Linear
};


/* TODO add
  // the pwiz code works on single data arrays (e.g., mz dimension) while we have peak structs. TODO: needs adaptation
  void predict(PredictionType pred, PrecisionType prec)
  {
      using OpenMS::Internal::MzMLHandlerHelper;
      
      switch (pred)
      {
          case PredictionType::Prediction_Delta:
              switch (p)
              {
                  case BinaryData::PRE_32:
                      for (size_t i = 2; i < binaryDataArray->data.size(); i++) 
                        binaryDataArray->data[i] = (float)binaryDataArray->data[i] + (float)binaryDataArray->data[i - 1] - (float)binaryDataArray->data[0];
                      break;
                  case BinaryData::PRE_64:
                      for (size_t i = 2; i < binaryDataArray->data.size(); i++)
                        binaryDataArray->data[i] = binaryDataArray->data[i] + binaryDataArray->data[i - 1] - binaryDataArray->data[0];
                      break;
              }
              break;
          case PredictionType::Prediction_Linear:
              switch (prec)
              {
                  case BinaryData::PRE_32:
                      for (size_t i = 2; i < binaryDataArray->data.size(); i++)
                        binaryDataArray->data[i] = (float)binaryDataArray->data[i] + 2.0f * (float)binaryDataArray->data[i - 1] - (float)binaryDataArray->data[i - 2] - (float)binaryDataArray->data[1];
                      break;
                  case BinaryData::PRE_64:
                      for (size_t i = 2; i < binaryDataArray->data.size(); i++)
                        binaryDataArray->data[i] = binaryDataArray->data[i] + 2.0 * binaryDataArray->data[i - 1] - binaryDataArray->data[i - 2] - binaryDataArray->data[1];
                      break;
              }
              break;
          case PredictionType::Prediction_None:
              break;
      }
  }

*/



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

class MzMLbFile
{
  public:
    MzMLbFile()
    {
      // load blosc plugin (could be part of a HDF5 singleton if we use it somewhere else)
      char *version, *date;
      auto return_code = register_blosc(&version, &date);
      TEST_EQUAL(return_code >= 0, true);
      std::cout << "Blosc version info: " << version << " " << date << std::endl;
    }

    MSExperiment load(const std::string& file_name)
    {
      // open mzMLb file
      auto mzMLb = MzMLb(file_name);
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
      return exp;
    }
};

START_SECTION((MzMLb()))
{
  const std::string filename( OPENMS_GET_TEST_DATA_PATH("msconvert.0.24017-6a003b2.mzMLb") ); // file converted with pwiz
  auto mzmlb = MzMLbFile();
  MSExperiment exp = mzmlb.load(filename);

/*
  for (Size i = 0; i != exp.getNrSpectra(); ++i)
  {
    // retrieve datasets, offsets, for data extraction
    auto mzs_in_hdf5 = mzfile.getMzMLbMapping(MSSpectrum, MZ_DIMENSION, i); // TODO:
    auto ints_in_hdf5 = mzfile.getMzMLbMapping(MSSpectrum, INT_DIMENSION, i);
    // TODO: other data arrays
    fillSpectrum(is, i, mzs_in_hdf5, ints_in_hdf5);
  }

  for (Size i = 0; i != exp.getChromatograms(); ++i)
  {
    // retrieve datasets, offsets, for data extraction
    auto rts_in_hdf5 = mzfile.getMzMLbMapping(MSCHROMATOGRAM, RT_DIMENSION, i);
    auto ints_in_hdf5 = mzfile.getMzMLbMapping(MSCHROMATOGRAM, INT_DIMENSION, i);
    // TODO: other data arrays
    fillChromatogram(is, i, rts_in_hdf5, ints_in_hdf5);
  }
*/

}
END_SECTION

END_TEST
