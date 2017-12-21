// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Mathias Walzer $
// $Authors: Andreas Bertsch, Mathias Walzer $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MZIDENTMLFILE_H
#define OPENMS_FORMAT_MZIDENTMLFILE_H

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>

namespace OpenMS
{
  /**
      @brief File adapter for MzIdentML files

      This file adapter exposes the internal MzIdentML processing capabilities to the library. The file
      adapter interface is kept the same as idXML file adapter for downward capability reasons.
      For now, read-in will be performed with DOM write-out with STREAM

      @note due to the limited capabilities of idXML/PeptideIdentification/ProteinIdentification not all
        MzIdentML features can be supported. Development for these structures will be discontinued, a new
        interface with appropriate structures will be provided.
      @note If a critical error occurs due to the missing functionality, Exception::NotImplemented is thrown.

      @note All PSM will be read into PeptideIdentification, even the passThreshold=false, threshold will be
        read into ProteinIdentification (i.e. one id run), considered at writing also will only be the
        threshold set in ProteinIdentification
      @note All PSM will be read into PeptideIdentification, even the passThreshold=false

      @ingroup FileIO
  */
  class OPENMS_DLLAPI MzIdentMLFile :
    public Internal::XMLFile,
    public ProgressLogger
  {
public:
    ///Default constructor
    MzIdentMLFile();
    ///Destructor
    ~MzIdentMLFile() override;

    /**
        @brief Loads the identifications from a MzIdentML file.

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String& filename, Identification& id);

    /**
        @brief Loads the identifications from a MzIdentML file.

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String& filename, std::vector<ProteinIdentification>& poid, std::vector<PeptideIdentification>& peid);

    /**
        @brief Stores the identifications in a MzIdentML file.

        @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const String& filename, const std::vector<ProteinIdentification>& poid, const std::vector<PeptideIdentification>& peid) const;

    /**
        @brief Stores the identifications in a MzIdentML file.

        @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const String& filename, const Identification& id) const;

    /**
        @brief Checks if a file is valid with respect to the mapping file and the controlled vocabulary.

        @param filename File name of the file to be checked.
        @param errors Errors during the validation are returned in this output parameter.
        @param warnings Warnings during the validation are returned in this output parameter.

        @exception Exception::FileNotFound is thrown if the file could not be opened
    */
    bool isSemanticallyValid(const String& filename, StringList& errors, StringList& warnings);

private:

  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_MZIDENTMLFILE_H
