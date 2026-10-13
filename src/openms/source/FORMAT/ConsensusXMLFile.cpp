// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Marc Sturm, Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ConsensusXMLFile.h>

#include <OpenMS/FORMAT/HANDLERS/ConsensusXMLHandler.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/METADATA/ProteinRunReferences.h>
#include <fstream>
#include <optional>

using namespace std;

namespace OpenMS
{
  ConsensusXMLFile::ConsensusXMLFile() :
    XMLFile("/SCHEMAS/ConsensusXML_1_7.xsd", "1.7")
  {
  }

  ConsensusXMLFile::~ConsensusXMLFile() = default;

  PeakFileOptions& ConsensusXMLFile::getOptions()
  {
    return options_;
  }

  const PeakFileOptions& ConsensusXMLFile::getOptions() const
  {
    return options_;
  }

  namespace
  {
    // The identifications of a file as peptide and protein identifications of the map
    void setPeptideIdentifications(ConsensusMap& map, MapPeptideIdentifications&& identifications)
    {
      map.setProteinIdentifications(std::move(identifications.proteins));
      map.setUnassignedPeptideIdentifications(std::move(identifications.unassigned));
      for (auto& [path, peptides] : identifications.features)
      {
        map.at(path[0]).setPeptideIdentifications(std::move(peptides));
      }
    }
  } // namespace

  void ConsensusXMLFile::store(const std::string& filename, const ConsensusMap& consensus_map)
  {
    // consensusXML holds peptide and protein identifications: those of the map, or its identification data exported
    const MapPeptideIdentifications identifications = IdentificationDataConverter::exportMapIdentifications(consensus_map);
    if (!FileHandler::hasValidExtension(filename, FileTypes::CONSENSUSXML))
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "invalid file extension, expected '" + FileTypes::typeToName(FileTypes::CONSENSUSXML) + "'");
    }

    // every peptide identification needs its protein identification run, and the proteins of its peptide evidences need
    // to be protein hits of that run (they are referenced there): fail before the file is opened (one pass for both)
    ProteinRunReferences::checkProteinAccessions(identifications.proteins, identifications.peptides());

    if (!consensus_map.isMapConsistent(&getGlobalLogWarn()))
    {
      // Currently it is possible that FeatureLinkerUnlabeledQT triggers this exception
      // throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The ConsensusXML file contains invalid maps or references thereof. No data was written! Please fix the file or notify the maintainer of this tool if you did not provide a consensusXML file!");
      std::cerr << "The ConsensusXML file contains invalid maps or references thereof. Please fix the file or notify the maintainer of this tool if you did not provide a consensusXML file! Note that this warning will be a fatal error in the next version of OpenMS!" << std::endl;
    }

    if (Size invalid_unique_ids = consensus_map.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId))
    {
      // TODO Take care *outside* that this does not happen.
      // We can detect this here but it is too late to fix the problem;
      // there is no straightforward action to be taken in all cases.
      // Note also that we are given a const reference.
      OPENMS_LOG_INFO << std::string("ConsensusXMLFile::store():  found ") + invalid_unique_ids + " invalid unique ids" << std::endl;
    }

    // This will throw if the unique ids are not unique,
    // so we never create bad files in this respect.
    try
    {
      consensus_map.updateUniqueIdToIndex();
    }
    catch (Exception::Postcondition& e)
    {
      OPENMS_LOG_FATAL_ERROR << e.getName() << ' ' << e.what() << std::endl;
      throw;
    }


    Internal::ConsensusXMLHandler handler(consensus_map, identifications, filename);
    handler.setOptions(options_);
    handler.setLogType(getLogType());
    save_(filename, &handler);
  }

  void ConsensusXMLFile::load(const std::string& filename, ConsensusMap& consensus_map)
  {
    consensus_map.clear(true); // clear map

    //set DocumentIdentifier
    consensus_map.setLoadedFileType(filename);
    consensus_map.setLoadedFilePath(filename);

    Internal::ConsensusXMLHandler handler(consensus_map, filename);
    handler.setOptions(options_);
    handler.setLogType(getLogType());
    parse_(filename, &handler);
    setPeptideIdentifications(consensus_map, std::move(handler.getIdentifications()));

    if (!consensus_map.isMapConsistent(&getGlobalLogWarn())) // a warning is printed to LOG_WARN during isMapConsistent()
    {
      // don't throw exception for now, since this would prevent us from reading old files...
      // throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The ConsensusXML file contains invalid maps or references thereof. Please fix the file!");

    }
    consensus_map.updateRanges();

  }

} // namespace OpenMS
