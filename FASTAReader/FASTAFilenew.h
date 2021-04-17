#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <functional>
#include <fstream>
#include <memory>
#include <utility>
#include <vector>

namespace OpenMS
{
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
    @return true if entry was read; false if eof was reached
    @exception Exception::FileNotFound is thrown if the file does not exists.
    @exception Exception::ParseError is thrown if the file does not suit to the standard.
    */
    bool readNext(FASTAEntry& protein);

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

    //eigene Implementierung des readRecord
    bool FASTAFile::readRecordNew(std::string & id, std::string & seq); //die soll den vorhandenen infile_ benutzen

protected:
    std::fstream infile_;   ///< filestream for reading; init using FastaFile::readStart()
    std::ofstream outfile_; ///< filestream for writing; init using FastaFile::writeStart()
    Size entries_read_; ///< some internal book-keeping during reading
    unsigned fileSize_{};
  };

} // namespace OpenMS

