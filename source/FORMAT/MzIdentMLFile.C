// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch, Mathias Walzer$
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/VALIDATORS/MzIdentMLValidator.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>
#include <OpenMS/FORMAT/HANDLERS/MzIdentMLHandler.h>
#include <OpenMS/SYSTEM/File.h>

namespace OpenMS
{

  MzIdentMLFile::MzIdentMLFile() :
    XMLFile("/SCHEMAS/mzIdentML1.1.0.xsd", "1.1.0")
  {
  }

  MzIdentMLFile::~MzIdentMLFile()
  {
  }

  void MzIdentMLFile::load(const String & filename, Identification & id)
  {
    Internal::MzIdentMLHandler handler(id, filename, schema_version_, *this);
    parse_(filename, &handler);
  }

  void MzIdentMLFile::store(const String & filename, const Identification & id) const
  {
    Internal::MzIdentMLHandler handler(id, filename, schema_version_, *this);
    save_(filename, &handler);
  }

  void MzIdentMLFile::store(const String & filename, const std::vector<ProteinIdentification> & poid, const std::vector<PeptideIdentification> & peid) const
  {
    Internal::MzIdentMLHandler handler(poid, peid, filename, schema_version_, *this);
    save_(filename, &handler);
  }

  bool MzIdentMLFile::isSemanticallyValid(const String & filename, StringList & errors, StringList & warnings)
  {
    //load mapping
    CVMappings mapping;
    CVMappingFile().load(File::find("/MAPPING/mzIdentML-mapping.xml"), mapping);

    //load cvs
    ControlledVocabulary cv;
    cv.loadFromOBO("MS", File::find("/CV/psi-ms.obo"));
    cv.loadFromOBO("PATO", File::find("/CV/quality.obo"));
    cv.loadFromOBO("UO", File::find("/CV/unit.obo"));
    cv.loadFromOBO("BTO", File::find("/CV/brenda.obo"));
    cv.loadFromOBO("GO", File::find("/CV/goslim_goa.obo"));

    //validate
    Internal::MzIdentMLValidator v(mapping, cv);
    bool result = v.validate(filename, errors, warnings);

    return result;
  }

} // namespace OpenMS
