// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest$
// $Authors: Hannes Roest$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/SwathWindowLoader.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/SYSTEM/SysInfo.h>

#include <fstream>
#include <iostream>
#include <sstream>

namespace OpenMS
{

  void SwathWindowLoader::annotateSwathMapsFromFile(const std::string & filename,
    std::vector< OpenSwath::SwathMap >& swath_maps, bool do_sort, bool force)
  {
    std::vector<double> swath_prec_lower_, swath_prec_upper_;
    readSwathWindows(filename, swath_prec_lower_, swath_prec_upper_);

    // Sort the windows by the start of the lower window
    if (do_sort)
    {
      std::sort(swath_maps.begin(), swath_maps.end(), [](const OpenSwath::SwathMap& left, const OpenSwath::SwathMap& right)
        {
          return left.upper < right.upper;
        });
    }

    Size i = 0, j = 0;
    for (; i < swath_maps.size(); i++)
    {
      if (swath_maps[i].ms1)
      { // skip to next map (only increase i)
        continue;
      }

      if (j >= swath_prec_lower_.size())
      {
        OPENMS_LOG_FATAL_ERROR << "Trying to access annotation for SWATH map " << j
                  << " but there are only " << swath_prec_lower_.size() << " windows in the"
                  << " swath_windows_file. Please check your input." << std::endl;
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "The number of SWATH maps read from the raw data and from the annotation file do not match.");
      }

      OPENMS_LOG_INFO << "Re-annotate from file: SWATH " <<
        swath_maps[i].lower << " / " << swath_maps[i].upper << " (raw data) is annotated via swath_windows_file with " <<
        swath_prec_lower_[j] << " / " << swath_prec_upper_[j] << std::endl;
      
      // new boundaries should be smaller/equal than the original ones from the data
      if (!(swath_maps[i].lower <= swath_prec_lower_[j] && swath_prec_upper_[j] <= swath_maps[i].upper))
      { 
        String err = "SWATH window #" + String(j+1) + " from swath_windows_file extends beyond the Swath window of the data."
                     " Did you forget to apply the sort_swath_maps flag? (override with -force)";
        if (force)
        {
          OPENMS_LOG_ERROR << err << "\nOverridden with -force.\n";
        }
        else
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, err);
        }
      }

      swath_maps[i].lower = swath_prec_lower_[j];
      swath_maps[i].upper = swath_prec_upper_[j];
      j++;
    }

    if (j != swath_prec_upper_.size())
    {
      OPENMS_LOG_FATAL_ERROR << "The number of SWATH maps read from the raw data (" <<
        j << ") and from the annotation file (" << swath_prec_upper_.size() << ") do not match." << std::endl;
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "The number of SWATH maps read from the raw data and from the annotation file do not match.");
    }
  }

  void SwathWindowLoader::readSwathWindows(const std::string & filename,
    std::vector<double> & swath_prec_lower_,
    std::vector<double> & swath_prec_upper_ )
  {
    std::ifstream data(filename.c_str());
    String line;
    std::vector<String> headerSubstrings;
    double lower, upper;

    // Check for presence of header
    std::getline(data, line);
    try 
    { // If string can be successfully converted to double (excluding initial spaces) then the first line is not a header
      StringUtils::split(line.trim().substitute('\t', ' '), ' ', headerSubstrings);
      StringUtils::toDouble(headerSubstrings[0]);
      OPENMS_LOG_INFO << "Swath Header not found" << std::endl;
    }
    catch (const Exception::ConversionError &)
    {
      OPENMS_LOG_INFO << "Read Swath window header: '" << line << std::endl;
      std::getline(data, line);
    }

    // read the rest of the SWATH window file
    do
    {
      std::stringstream lineStream(line);

      lineStream >> lower;
      lineStream >> upper;

      swath_prec_lower_.push_back(lower);
      swath_prec_upper_.push_back(upper);
      if (!(lower < upper))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Swath window file contains illegal ranges", line);
      }
    } while (std::getline(data, line));

    assert(swath_prec_lower_.size() == swath_prec_upper_.size());
    OPENMS_LOG_INFO << "Read Swath window file with " << swath_prec_lower_.size() << " SWATH windows." << std::endl;
  }

}
