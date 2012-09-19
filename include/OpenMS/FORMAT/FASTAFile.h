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
// $Maintainer: Sandro Andreotti $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_FASTAFILE_H
#define OPENMS_FORMAT_FASTAFILE_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief This class serves for reading in FASTA files
  */
  class OPENMS_DLLAPI FASTAFile
  {
public:

    /**
      @brief FASTA entry type (identifier, description and sequence)

      The first String corresponds to the identifier that is
      written after the > in the FASTA file. The part after the
      first whitespace is stored in description and the text
      from the next line until the next > (exclusive) is stored
      in sequence.
    */
    struct FASTAEntry
    {
      String identifier;
      String description;
      String sequence;

      FASTAEntry() :
        identifier(""),
        description(""),
        sequence("")
      {
      }

      FASTAEntry(String id, String desc, String seq) :
        identifier(id),
        description(desc),
        sequence(seq)
      {
      }

      bool operator==(const FASTAEntry & rhs) const
      {
        return identifier == rhs.identifier
               && description == rhs.description
               && sequence == rhs.sequence;
      }

    };

    /// Copy constructor
    FASTAFile();

    /// Destructor
    virtual ~FASTAFile();

    /**
      @brief loads a FASTA file given by 'filename' and stores the information in 'data'

      @exception Exception::FileNotFound is thrown if the file does not exists.
      @exception Exception::ParseError is thrown if the file does not suit to the standard.
    */
    void load(const String & filename, std::vector<FASTAEntry> & data);

    /**
      @brief stores the data given by 'data' at the file 'filename'

      @exception Exception::UnableToCreateFile is thrown if the process is not able to write the file.
    */
    void store(const String & filename, const std::vector<FASTAEntry> & data) const;

  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_FASTAFILE_H
