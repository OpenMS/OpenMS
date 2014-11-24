// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_XTANDEMXMLFILE_H
#define OPENMS_FORMAT_XTANDEMXMLFILE_H

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>

#include <vector>

namespace OpenMS
{
  class String;
  class ProteinIdentification;

  /**
    @brief Used to load XTandemXML files

    This class is used to load documents that implement
    the schema of XTandemXML files.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI XTandemXMLFile :
    protected Internal::XMLHandler,
    public Internal::XMLFile
  {
public:

    /// Default constructor
    XTandemXMLFile();

    /// Destructor
    virtual ~XTandemXMLFile();
    /**
      @brief loads data from a XTandemXML file

      @param filename the file to be loaded
      @param protein_identification protein identifications belonging to the whole experiment
      @param id_data the identifications with m/z and RT

      This class serves to read in a XTandemXML file. The information can be
      retrieved via the load function.

      @ingroup FileIO
    */
    void load(const String & filename, ProteinIdentification & protein_identification, std::vector<PeptideIdentification> & id_data);


    /// sets the valid modifications
    void setModificationDefinitionsSet(const ModificationDefinitionsSet & rhs);

protected:

    // Docu in base class
    void startElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname, const xercesc::Attributes & attributes);

    // Docu in base class
    void endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname);

    // Docu in base class
    void characters(const XMLCh * const chars, const XMLSize_t /*length*/);

    XTandemXMLFile(const XTandemXMLFile & rhs);

    XTandemXMLFile & operator=(const XTandemXMLFile & rhs);

private:

    ProteinIdentification * protein_identification_;

    // used to indicate that an protein tag is open
    bool protein_open_;

    // true if actual
    bool is_description_;

    // peptide hits of one spectrum
    Map<UInt, std::vector<PeptideHit> > peptide_hits_;

    // protein hits, sorted by the id
    Map<String, ProteinHit> protein_hits_;

    // id of the actual protein
    String actual_protein_id_;

    // charge of actual peptide
    Int actual_charge_;

    // id of actual peptide
    Int actual_id_;

    // tag
    String tag_;

    // actual start position of peptide in protein sequence
    UInt actual_start_;

    // actual stop position of peptide in protein sequence
    UInt actual_stop_;

    // modification definitions
    ModificationDefinitionsSet mod_def_set_;

  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_XTANDEMXMLFILE_H
