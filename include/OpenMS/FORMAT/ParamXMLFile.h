// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_PARAMXMLFILE_H
#define OPENMS_FORMAT_PARAMXMLFILE_H

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

namespace OpenMS
{
  /**
    @brief The file pendant of the Param class used to load and store the param
           datastructure as paramXML.
  */
  class OPENMS_DLLAPI ParamXMLFile :
    public Internal::XMLFile
  {
public:
    /// Constructor.
    ParamXMLFile();

    /**
      @brief Write XML file.

      @param filename The filename where the param data structure should be stored.
      @param param The Param class that should be stored in the file.

      @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const String& filename, const Param& param) const;

    /**
      @brief Write XML to output stream.

      @param os_ptr The stream where the param class should be written to.
      @param param The Param class that should be written to the stream.
    */
    void writeXMLToStream(std::ostream* os_ptr, const Param& param) const;

    /**
      @brief Read XML file.

      @param filename The file from where to read the Param object.
      @param param The param object where the read data should be stored.

      @exception Exception::FileNotFound is thrown if the file could not be found
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String& filename, Param& param);
  };

} // namespace OpenMS


#endif // OPENMS_FORMAT_PARAMXMLFILE_H
