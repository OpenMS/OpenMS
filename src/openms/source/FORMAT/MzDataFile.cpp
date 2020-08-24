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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/VALIDATORS/MzDataValidator.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>

namespace OpenMS
{
  MzDataFile::MzDataFile() :
    XMLFile("/SCHEMAS/mzData_1_05.xsd", "1.05"),
    options_()
  {
  }

  MzDataFile::~MzDataFile()
  {
  }

  PeakFileOptions & MzDataFile::getOptions()
  {
    return options_;
  }

  const PeakFileOptions & MzDataFile::getOptions() const
  {
    return options_;
  }

  void MzDataFile::setOptions(const PeakFileOptions & options)
  {
      options_ = options;
  }

  bool MzDataFile::isSemanticallyValid(const String & filename, StringList & errors, StringList & warnings)
  {
    //load mapping
    CVMappings mapping;
    CVMappingFile().load(File::find("/MAPPING/mzdata-mapping.xml"), mapping);

    //load cvs
    ControlledVocabulary cv;
    cv.loadFromOBO("PSI", File::find("/CV/psi-mzdata.obo"));

    //validate
    Internal::MzDataValidator v(mapping, cv);
    bool result = v.validate(filename, errors, warnings);

    return result;
  }

  void MzDataFile::load(const String & filename, PeakMap & map)
  {
    map.reset();

    //set DocumentIdentifier
    map.setLoadedFileType(filename);
    map.setLoadedFilePath(filename);

    Internal::MzDataHandler handler(map, filename, schema_version_, *this);
    handler.setOptions(options_);
    parse_(filename, &handler);
  }

  void MzDataFile::store(const String & filename, const PeakMap & map) const
  {
    Internal::MzDataHandler handler(map, filename, schema_version_, *this);
    handler.setOptions(options_);
    save_(filename, &handler);
  }


} // namespace OpenMS
