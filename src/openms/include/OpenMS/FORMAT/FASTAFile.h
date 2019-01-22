// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Authors: Nico Pfeifer, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/seq_io/guess_stream_format.h>
#include <seqan/seq_io/read_fasta_fastq.h>

#include <functional>
#include <fstream>
#include <memory>
//#include <utility>
#include <vector>

namespace OpenMS
{

  namespace Internal
  {
    template <bool SHRINK_SEQ> inline void doShrink(String& /*s*/);
    template <> inline void doShrink<true>(String& s)
    { // optimize memory usage (10-15% savings), at the cost of runtime;
      s.shrink_to_fit();
    }
    template <> inline void doShrink<false>(String& /*s*/) {} //do nothing
  }

  /**
    @brief This class serves for reading in and writing FASTA files

    If the protein/gene sequence contains unusual symbols (such as translation end (*)),
    they will be kept!

    You can use aggregate methods load() and store() to read/write a
    set of protein sequences at the cost of memory.
    
    Or use single read/write of protein sequences using readStart(), readNext()
    and writeStart(), writeNext(), writeEnd() for more memory efficiency.
    Reading from one and writing to another FASTA file can be handled by 
    one single FASTAFile instance.

  */

  class OPENMS_DLLAPI FASTAFile
  {
    typedef seqan::RecordReader<std::fstream, seqan::SinglePass<> > FASTARecordReader;

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
        identifier(),
        description(),
        sequence()
      {
      }

      FASTAEntry(String id, String desc, String seq) :
        identifier(id),
        description(desc),
        sequence(seq)
      {
      }
      
      FASTAEntry(const FASTAEntry& rhs)
        :
        identifier(rhs.identifier),
        description(rhs.description),
        sequence(rhs.sequence)
      {
      }

      FASTAEntry(FASTAEntry&& rhs) noexcept
       :
        identifier(::std::move(rhs.identifier)),
        description(::std::move(rhs.description)),
        sequence(::std::move(rhs.sequence)) 
      {
      }

      FASTAEntry& operator=(const FASTAEntry& rhs)
      {
        if (*this == rhs) return *this;
        identifier = rhs.identifier;
        description = rhs.description;
        sequence = rhs.sequence;
        return *this;
      }

      bool operator==(const FASTAEntry& rhs) const
      {
        return identifier == rhs.identifier
               && description == rhs.description
               && sequence == rhs.sequence;
      }
    
      bool headerMatches(const FASTAEntry& rhs) const
      {
        return identifier == rhs.identifier && 
  	     description == rhs.description;
      }
 
      bool sequenceMatches(const FASTAEntry& rhs) const
      {
        return sequence == rhs.sequence;
      }
    };

    /// Default constructor
    FASTAFile();

    /// Destructor
    virtual ~FASTAFile();

    /**
      @brief Prepares a FASTA file given by 'filename' for streamed reading using readNext().

      @exception Exception::FileNotFound is thrown if the file does not exists.
      @exception Exception::ParseError is thrown if the file does not suit to the standard.
    */
    void readStart(const String& filename);

    /**
    @brief Reads the next FASTA entry from file.

    If you want to read all entries in one go, use load().

    @param SHRINK_SEQ Use shrink_to_fit() on the sequence string. Can potentially save ~10-15% of memory, which is useful when loading
                      huge FASTA files into memory. However, if a single FASTAEntry is read over and over again, this will just slow down the program.

    @param protein The entry to fill.
    @return true if entry was read; false if Enf-Of-File was reached
    @exception Exception::FileNotFound is thrown if the file does not exists.
    @exception Exception::ParseError is thrown if the file does not suit to the standard.
    */
    template <bool SHRINK_SEQ>
    inline bool readNext(FASTAEntry& protein)
    {
      if (seqan::atEnd(*static_cast<FASTARecordReader*>(reader_.get())))
      {
        // do NOT close(), since we still might want to seek to certain positions
        return false;
      }
      String& id = protein.identifier;
      String& s = protein.sequence;
      if (readRecord(id, s, *static_cast<FASTARecordReader*>(reader_.get()), seqan::Fasta()) != 0)
      {
        if (entries_read_ == 0) s = "The first entry could not be read!";
        else s = "Only " + String(entries_read_) + " proteins could be read. The record after failed.";
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", "Error while parsing FASTA file! " + s + " Please check the file!");
      }
      ++entries_read_;
      s.removeWhitespaces();
      
      Internal::doShrink<SHRINK_SEQ>(s);

      // handle id
      id.trim();
      String::size_type position = id.find_first_of(" \v\t");
      if (position == String::npos)
      {
        // identifier is set implicitly (String& see above)
        protein.description.clear();
      }
      else
      {
        // do not swap lines (id contains the description)
        protein.description = id.suffix(id.size() - position - 1);
        id.resize(position);
      }

      return true;
    }

    /// current stream position
    std::streampos position() const;

    /// is stream at EOF?
    bool atEnd() const;

    /// seek stream to @p pos
    bool setPosition(const std::streampos& pos);

    /**
    @brief Prepares a FASTA file given by 'filename' for streamed writing using writeNext().

    @exception Exception::UnableToCreateFile is thrown if the process is not able to write to the file (disk full?).
    */
    void writeStart(const String& filename);

    /**
    @brief Stores the data given by @p protein. Call writeStart() once before calling writeNext().

    Call writeEnd() when done to close the file!

    @exception Exception::UnableToCreateFile is thrown if the process is not able to write the file.
    */
    void writeNext(const FASTAEntry& protein);

    /**
    @brief Closes the file (flush). Called implicitly when FASTAFile object does out of scope.

    */
    void writeEnd();
    

    /**
      @brief loads a FASTA file given by 'filename' and stores the information in 'data'

      This uses more RAM than readStart() and readNext().

      @exception Exception::FileNotFound is thrown if the file does not exists.
      @exception Exception::ParseError is thrown if the file does not suit to the standard.
    */
    void static load(const String& filename, std::vector<FASTAEntry>& data);

  /**
      @brief stores the data given by 'data' at the file 'filename'
      
      This uses more RAM than writeStart() and writeNext().

      @exception Exception::UnableToCreateFile is thrown if the process is not able to write the file.
    */
    void static store(const String& filename, const std::vector<FASTAEntry>& data);

protected:
    std::fstream infile_;   ///< filestream for reading; init using FastaFile::readStart()
    std::ofstream outfile_; ///< filestream for writing; init using FastaFile::writeStart()
    std::unique_ptr<void, std::function<void(void*) > > reader_; ///< filestream for reading; init using FastaFile::readStart(); needs to be a pointer, since its not copy-constructable; we use void* here, to avoid pulling in seqan includes
    Size entries_read_; ///< some internal book-keeping during reading
  };

} // namespace OpenMS

