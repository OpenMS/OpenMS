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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_PTMXMLFILE_H
#define OPENMS_FORMAT_PTMXMLFILE_H

#include <OpenMS/FORMAT/XMLFile.h>

#include <map>
#include <vector>

namespace OpenMS
{
  /**
      @brief Used to load and store PTMXML files

      This class is used to load and store documents that implement
      the schema of PTMXML files.

      @ingroup FileIO
  */
  class OPENMS_DLLAPI PTMXMLFile :
    public Internal::XMLFile
  {
public:
    /// Constructor
    PTMXMLFile();

    /**
        @brief Loads the informations of a PTMXML file

        @param filename The name of the file
        @param ptm_informations the PTM information from the file are stored herein
        @throw FileNotFound is thrown if the given file could not be found
        @throw ParseError is thrown if the given file could not be parsed
        The information is read in and stored in the corresponding variables
    */
    void load(const String & filename, std::map<String, std::pair<String, String> > & ptm_informations);

    /**
        @brief Stores the data in an PTMXML file

        @throw UnableToCreateFile is thrown if the given filename could not be created

        The data is read in and stored in the file 'filename'.
    */
    void store(String filename, std::map<String, std::pair<String, String> > & ptm_informations) const;
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_PTMXMLFILE_H
