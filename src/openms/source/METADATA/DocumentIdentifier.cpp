// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/DocumentIdentifier.h>

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <QDir>

namespace OpenMS
{

  DocumentIdentifier::DocumentIdentifier() :
    id_(),
    file_path_(),
    file_type_(FileTypes::UNKNOWN)
  {
  }

  DocumentIdentifier::~DocumentIdentifier()
  {
  }

  void DocumentIdentifier::setIdentifier(const String & id)
  {
    id_ = id;
  }

  const String & DocumentIdentifier::getIdentifier() const
  {
    return id_;
  }

  void DocumentIdentifier::setLoadedFilePath(const String & file_name)
  {
    // only change the path if we need to, otherwise low and upper case might be altered by Qt, making comparison in tests more tricky
    // i.e., a call to this will report unmatched strings
    //   FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"), e);
    //   TEST_STRING_EQUAL(e.getLoadedFilePath(), OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"));
    if (QDir::isRelativePath(file_name.toQString()))
    {
      file_path_ = File::absolutePath(file_name);
    }
    else
    {
      file_path_ = file_name;
    }
  }

  const String & DocumentIdentifier::getLoadedFilePath() const
  {
    return file_path_;
  }

  void DocumentIdentifier::setLoadedFileType(const String & file_name)
  {
    file_type_ = FileHandler::getTypeByContent(file_name);
  }

  const FileTypes::Type & DocumentIdentifier::getLoadedFileType() const
  {
    return file_type_;
  }

  void DocumentIdentifier::swap(DocumentIdentifier & from)
  {
    std::swap(id_, from.id_);
    std::swap(file_path_, from.file_path_);
    std::swap(file_type_, from.file_type_);
  }

  bool DocumentIdentifier::operator==(const DocumentIdentifier & rhs) const
  {
    return id_ == rhs.id_;
  }

}

