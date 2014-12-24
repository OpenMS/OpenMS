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

#ifndef OPENMS_FORMAT_CVMAPPINGFILE_H
#define OPENMS_FORMAT_CVMAPPINGFILE_H

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/DATASTRUCTURES/CVMappings.h>
#include <OpenMS/DATASTRUCTURES/CVMappingRule.h>

#include <vector>

namespace OpenMS
{
  class String;

  /**
    @brief Used to load CvMapping files

        This file contains the mapping of CV terms to the schema, which
        is used by PSI standard formats to semantically validate files.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI CVMappingFile :
    protected Internal::XMLHandler,
    public Internal::XMLFile
  {
public:

    /// Default constructor
    CVMappingFile();

    /// Destructor
    virtual ~CVMappingFile();

    /**
        @brief loads CvMappings from the given file

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
        @param strip_namespaces if enable, namespace definitions of the paths are eliminated, e.g. 'pf:cvParam' -> 'cvParam'
        @param cv_mappings  The CVMappings instance in which the rules, cvs and other content from the file should be stored
        @param filename  The filename to read from
    */
    void load(const String & filename, CVMappings & cv_mappings, bool strip_namespaces = false);

protected:

    // Docu in base class
    void startElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname, const xercesc::Attributes & attributes);

    // Docu in base class
    void endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname);

    // Docu in base class
    void characters(const XMLCh * const chars, const XMLSize_t /*length*/);

private:

    ///Not implemented
    CVMappingFile(const CVMappingFile & rhs);

    ///Not implemented
    CVMappingFile & operator=(const CVMappingFile & rhs);

    String tag_;

    bool strip_namespaces_;

    CVMappingRule actual_rule_;

    std::vector<CVMappingRule> rules_;

    std::vector<CVReference> cv_references_;

  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_CVMAPPINGFILE_H
