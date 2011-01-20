// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
 @page UTILS_UniqueIdAssigner UniqueIdAssigner

 @brief Assign new unique ids to FeatureXML or ConsensusXML files.

 <B>The command line parameters of this tool are:</B>
 @verbinclude UTILS_UniqueIdAssigner.cli
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

int
main(int argc, const char** argv)
{
  if ( argc != 4 )
  {
    std::cout << "Usage:  " << argv[0] << " method input [output]\n"
      "\n"
      "where:\n"
      "  method   indicates the method to be applied (none,ensure,reassign)\n"
      "           where:\n"
      "             none      =  do nothing, just read and write\n"
      "             ensure    =  if current unique id is invalid, assign a valid one\n"
      "             reassign  =  assign new valid unique ids unconditionally\n"
      "  input    is the input file, which must be in FeatureXML or ConsensusXML format\n"
      "  output   is the output file, which is written in the same format as input\n"
      "\n"
      "  WARNING!!!  If output is the special name '--overwrite', then input will be OVERWRITTEN!\n"
      "\n";
    return 1;
  }

  const char * const argv_method = argv[1];
  const char * const argv_input = argv[2];
  const char * argv_output = 0;
  if ( OpenMS::String (argv[3]) == "--overwrite" )
  {
    argv_output = argv[2];
  }
  else
  {
    argv_output = argv[3];
  }

  OpenMS::FileHandler file_handler;
  OpenMS::FileTypes::Type in_type = file_handler.getType(argv_input);

  if ( in_type == OpenMS::FileTypes::UNKNOWN )
  {
    std::cout << "Error: Could not determine input file type!" << std::endl;
    return 1;
  }

  try
  {
    if ( in_type == OpenMS::FileTypes::FEATUREXML )
    {

      OpenMS::FeatureXMLFile feature_file;
      OpenMS::FeatureMap<> feature_map;
      feature_file.load(argv_input, feature_map);

      std::string method(argv_method);
      if ( method == "reassign" )
      {
        feature_map.applyMemberFunction(&OpenMS::UniqueIdInterface::setUniqueId);
      }
      else if ( method == "ensure" )
      {
        feature_map.applyMemberFunction(&OpenMS::UniqueIdInterface::ensureUniqueId);
      }
      else if ( method == "none" )
      {
        // do nothing, but the output may be different nevertheless
      }
      else
      {
        std::cout << "unsupported method: " << method << std::endl;
      }

      feature_file.store(argv_output, feature_map);

    }
    else if ( in_type == OpenMS::FileTypes::CONSENSUSXML )
    {

      OpenMS::ConsensusXMLFile consensus_file;
      OpenMS::ConsensusMap consensus_map;
      consensus_file.load(argv_input, consensus_map);

      std::string method(argv_method);
      if ( method == "reassign" )
      {
        consensus_map.applyMemberFunction(&OpenMS::UniqueIdInterface::setUniqueId);
      }
      else if ( method == "ensure" )
      {
        consensus_map.applyMemberFunction(&OpenMS::UniqueIdInterface::ensureUniqueId);
      }
      else if ( method == "none" )
      {
        // do nothing, but the output may be different nevertheless
      }
      else
      {
        std::cout << "unsupported method: " << method << std::endl;
      }

      consensus_file.store(argv_output, consensus_map);

    }
    else
    {
      std::cout << "Error: unsupported input file type: " << file_handler.typeToName(in_type) << std::endl;
      return 1;
    }

  }
  catch ( ... )
  {
    std::cout << argv[0] << ' ' << argv[1] << ' ' << argv[2] << "  :  Something went wrong...n" << std::endl;
    return 2;
  }

  return 0;
}

/// @endcond
