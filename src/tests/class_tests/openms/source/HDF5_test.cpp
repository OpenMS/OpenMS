// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

///////////////////////////

#include <OpenMS/CONCEPT/Types.h>

#include <memory>

using namespace std;


// Example from https://support.hdfgroup.org/HDF5/doc/cpplus_RM/create_8cpp-example.html

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the COPYING file, which can be found at the root of the source code       *
 * distribution tree, or in https://support.hdfgroup.org/ftp/HDF5/releases.  *
 * If you do not have access to either file, you may request a copy from     *
 * help@hdfgroup.org.                                                        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <string>
#include "H5Cpp.h"
#include "blosc_filter.h"


using namespace H5;
const H5std_string  FILE_NAME( "SDS.h5" );
const H5std_string  DATASET_NAME( "IntArray" );
const int   NX = 5;                    // dataset dimensions
const int   NY = 6;
const int   RANK = 2;

START_TEST(HDF5, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

START_SECTION((HDF5()))
   /*
    * Data initialization.
    */
   int i, j;
   int data[NX][NY];          // buffer for data to write
   for (j = 0; j < NX; j++)
   {
     for (i = 0; i < NY; i++)
       data[j][i] = i + j;
   }
   /*
    * 0 1 2 3 4 5
    * 1 2 3 4 5 6
    * 2 3 4 5 6 7
    * 3 4 5 6 7 8
    * 4 5 6 7 8 9
    */
   // Try block to detect exceptions raised by any of the calls inside it
   // try
   {
      /*
       * Turn off the auto-printing when failure occurs so that we can
       * handle the errors appropriately
       */
      H5::Exception::dontPrint();
      /*
       * Create a new file using H5F_ACC_TRUNC access,
       * default file creation properties, and default file
       * access properties.
       */
      H5File file( FILE_NAME, H5F_ACC_TRUNC );
      /*
       * Define the size of the array and create the data space for fixed
       * size dataset.
       */
      hsize_t     dimsf[2];              // dataset dimensions
      dimsf[0] = NX;
      dimsf[1] = NY;
      DataSpace dataspace( RANK, dimsf );
      /*
       * Define datatype for the data in the file.
       * We will store little endian INT numbers.
       */
      IntType datatype( PredType::NATIVE_INT );
      datatype.setOrder( H5T_ORDER_LE );
      /*
       * Create a new dataset within the file using defined dataspace and
       * datatype and default dataset creation properties.
       */
      DataSet dataset = file.createDataSet( DATASET_NAME, datatype, dataspace );
      /*
       * Write the data to the dataset using default memory space, file
       * space, and transfer properties.
       */
      dataset.write( data, PredType::NATIVE_INT );
   }  // end of try block
   // // catch failure caused by the H5File operations
   // catch( FileIException error )
   // {
   //    error.printError();
   //    return -1;
   // }
   // // catch failure caused by the DataSet operations
   // catch( DataSetIException error )
   // {
   //    error.printError();
   //    return -1;
   // }
   // // catch failure caused by the DataSpace operations
   // catch( DataSpaceIException error )
   // {
   //    error.printError();
   //    return -1;
   // }
   // // catch failure caused by the DataSpace operations
   // catch( DataTypeIException error )
   // {
   //    error.printError();
   //    return -1;
   // }

END_SECTION


START_SECTION((HDF5_BLOSC()))
{
  char *version, *date;
  auto return_code = register_blosc(&version, &date);
  TEST_EQUAL(return_code >= 0, true);
  std::cout << "Blosc version info: " << version << " " << date << std::endl;
  
  // try to load a BLOSC compressed HDF5 file (mzMLB) in read/write mode
  unsigned int openFlag = H5F_ACC_RDONLY;
  const H5std_string  filename( OPENMS_GET_TEST_DATA_PATH("test.mzMLb") );

  // thread-safe way to release file handle and free pointer
  std::unique_ptr<H5File, void(*)(H5::H5File*)> hdf5_file(new H5File(filename, openFlag /*, fcparm, faparm*/), 
    [](H5::H5File* f) {
      if (f)
      {
        f->flush(H5F_SCOPE_LOCAL);
        f->close();
        delete f;
        f = nullptr;
      }
    });

    // get a list of all objects in the hdf5_file
    std::vector<std::string> object_names;
    hdf5_file->iterateElems("/", NULL,
      // a c-style call-back function that takes takes an (&object_names), and makes it available as op_data.
      // the lambda is called for each step during iteration
      [](hid_t /*loc_id*/, const char* name, void* op_data) -> herr_t
      {
        std::vector<std::string>* object_names = static_cast<std::vector<std::string>*>(op_data);
        object_names->push_back(name);
        return 0;
      }, 
      &object_names);

    // print the object names and types
    for (const std::string& name : object_names) 
    {
      H5O_info_t object_info;
      // note: called without "H5O_INFO_ALL" because we are currently using the hdf5 compatibility macros for 1.10
      H5Oget_info_by_name(hdf5_file->getId(), name.c_str(), &object_info, H5P_DEFAULT/*, H5O_INFO_ALL*/); 

      if (object_info.type == H5O_TYPE_GROUP) 
      {
        std::cout << name << " (group)" << std::endl;
      } 
      else if (object_info.type == H5O_TYPE_DATASET) 
      {
        hid_t dataset_id = H5Dopen(hdf5_file->getId(), name.c_str(), H5P_DEFAULT);
        hid_t datatype_id = H5Dget_type(dataset_id);
        H5T_class_t datatype_class = H5Tget_class(datatype_id);

        if (datatype_class == H5T_COMPOUND) 
        {
          std::cout << name << " (compound)" << std::endl;
          int nmembers = H5Tget_nmembers(datatype_id);

          for (int i = 0; i < nmembers; ++i) 
          {
            std::string member_name = H5Tget_member_name(datatype_id, i);
            std::cout << "    " << member_name << std::endl;
          }
        } 
        else if (datatype_class == H5T_INTEGER) 
        {
          std::cout << name << " (integer)" << std::endl;
        } 
        else if (datatype_class == H5T_FLOAT) 
        {
          std::cout << name << " (float)" << std::endl;
        } 
        else if (datatype_class == H5T_STRING) 
        {
          std::cout << name << " (string)" << std::endl;
        }

        H5Tclose(datatype_id);
        H5Dclose(dataset_id);
      } 
      else 
      {
        std::cout << name << " (unknown)" << std::endl;
      }
    }  
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
