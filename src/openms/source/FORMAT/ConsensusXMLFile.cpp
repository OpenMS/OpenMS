// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Marc Sturm, Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ConsensusXMLFile.h>

#include <OpenMS/FORMAT/HANDLERS/ConsensusXMLHandler.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/SYSTEM/File.h>
#include <fstream>

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

  void ConsensusXMLFile::store(const String& filename, const ConsensusMap& consensus_map)
  {
    if (!FileHandler::hasValidExtension(filename, FileTypes::CONSENSUSXML))
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "invalid file extension, expected '" + FileTypes::typeToName(FileTypes::CONSENSUSXML) + "'");
    }

    if (!consensus_map.isMapConsistent(&OpenMS_Log_warn))
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
      OPENMS_LOG_INFO << String("ConsensusXMLFile::store():  found ") + invalid_unique_ids + " invalid unique ids" << std::endl;
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


    Internal::ConsensusXMLHandler handler(consensus_map, filename);
    handler.setOptions(options_);
    handler.setLogType(getLogType());
    save_(filename, &handler);
  }

  void ConsensusXMLFile::load(const String& filename, ConsensusMap& consensus_map)
  {
    consensus_map.clear(true); // clear map

    //set DocumentIdentifier
    consensus_map.setLoadedFileType(filename);
    consensus_map.setLoadedFilePath(filename);

    Internal::ConsensusXMLHandler handler(consensus_map, filename);
    handler.setOptions(options_);
    handler.setLogType(getLogType());
    parse_(filename, &handler);

    if (!consensus_map.isMapConsistent(&OpenMS_Log_warn)) // a warning is printed to LOG_WARN during isMapConsistent()
    {
      // don't throw exception for now, since this would prevent us from reading old files...
      // throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The ConsensusXML file contains invalid maps or references thereof. Please fix the file!");

    }

  }

} // namespace OpenMS
