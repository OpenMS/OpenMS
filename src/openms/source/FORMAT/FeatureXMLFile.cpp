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
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FeatureXMLFile.h>

#include <OpenMS/FORMAT/HANDLERS/FeatureXMLHandler.h>

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <fstream>

using namespace std;

namespace OpenMS
{
  FeatureXMLFile::FeatureXMLFile() :
    Internal::XMLFile("/SCHEMAS/FeatureXML_1_9.xsd", "1.9")
  {
  }

  FeatureXMLFile::~FeatureXMLFile()
  {
  }

  Size FeatureXMLFile::loadSize(const String& filename)
  {
    FeatureMap dummy;
    Internal::FeatureXMLHandler handler(dummy, filename);
    handler.setOptions(options_);
    handler.setSizeOnly(true);
    handler.setLogType(getLogType());
    parse_(filename, &handler);

    return handler.getSize();
  }

  void FeatureXMLFile::load(const String& filename, FeatureMap& feature_map)
  {
    feature_map.clear(true);
    //set DocumentIdentifier
    feature_map.setLoadedFileType(filename);
    feature_map.setLoadedFilePath(filename);

    Internal::FeatureXMLHandler handler(feature_map, filename);
    handler.setOptions(options_);
    handler.setLogType(getLogType());
    parse_(filename, &handler);

    // !!! Hack: set feature FWHM from meta info entries as
    // long as featureXML doesn't support a width entry.
    // See also hack in BaseFeature::setWidth().
    for (auto& feature : feature_map)
    {
      if (feature.metaValueExists("FWHM"))
      {
        feature.setWidth((double)feature.getMetaValue("FWHM"));
      }
    }

    // put ranges into defined state
    feature_map.updateRanges();
  }

  void FeatureXMLFile::store(const String& filename, const FeatureMap& feature_map)
  {

    if (!FileHandler::hasValidExtension(filename, FileTypes::FEATUREXML))
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "invalid file extension, expected '" + FileTypes::typeToName(FileTypes::FEATUREXML) + "'");
    }

    if (Size invalid_unique_ids = feature_map.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId))
    {

      // TODO Take care *outside* that this does not happen.
      // We can detect this here but it is too late to fix the problem;
      // there is no straightforward action to be taken in all cases.
      // Note also that we are given a const reference.
      OPENMS_LOG_INFO << String("FeatureXMLHandler::store():  found ") + invalid_unique_ids + " invalid unique ids" << std::endl;
    }

    // This will throw if the unique ids are not unique,
    // so we never create bad files in this respect.
    try
    {
      feature_map.updateUniqueIdToIndex();
    }
    catch (Exception::Postcondition& e)
    {
      OPENMS_LOG_FATAL_ERROR << e.getName() << ' ' << e.what() << std::endl;
      throw;
    }

    Internal::FeatureXMLHandler handler(feature_map, filename);
    handler.setOptions(options_);
    handler.setLogType(getLogType());
    save_(filename, &handler);
  }

  FeatureFileOptions& FeatureXMLFile::getOptions()
  {
    return options_;
  }

  const FeatureFileOptions& FeatureXMLFile::getOptions() const
  {
    return options_;
  }

  void FeatureXMLFile::setOptions(const FeatureFileOptions& options)
  {
    options_ = options;
  }


}
