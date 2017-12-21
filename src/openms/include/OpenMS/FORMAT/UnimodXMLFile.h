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
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_UNIMODXMLFILE_H
#define OPENMS_FORMAT_UNIMODXMLFILE_H

#include <OpenMS/FORMAT/XMLFile.h>

namespace OpenMS
{
  class ResidueModification;

  /**
    @brief Used to load XML files from unimod.org files

    @ingroup FileIO
  */
  class OPENMS_DLLAPI UnimodXMLFile :
    public Internal::XMLFile
  {
public:

    /// Default constructor
    UnimodXMLFile();

    /// Destructor
    ~UnimodXMLFile() override;
    /**
      @brief loads data from unimod.xml file

          @param filename the filename were the unimod xml file should be read from
          @param modifications the modifications which are read from the file
          @throw FileNotFound is thrown if the file could not be found
          @throw ParseError is thrown if the given file could not be parsed

      @ingroup FileIO
    */
    void load(const String & filename, std::vector<ResidueModification *> & modifications);

private:

    ///Not implemented
    UnimodXMLFile(const UnimodXMLFile & rhs);
    ///Not implemented
    UnimodXMLFile & operator=(const UnimodXMLFile & rhs);

  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_UNIMODXMLFILE_H
