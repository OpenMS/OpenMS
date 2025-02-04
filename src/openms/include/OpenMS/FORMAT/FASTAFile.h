// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow, Nora Wild $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <fstream>
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

    class OPENMS_DLLAPI FASTAFile : public ProgressLogger
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

            FASTAEntry() = default;

            FASTAEntry(const String& id, const String& desc, const String& seq) :
                    identifier(id),
                    description(desc),
                    sequence(seq)
            {
            }

            FASTAEntry(const FASTAEntry& rhs) = default;

            FASTAEntry(FASTAEntry&& rhs) noexcept 
                    :
                    identifier(::std::move(rhs.identifier)),
                    description(::std::move(rhs.description)),
                    sequence(::std::move(rhs.sequence))
            {
            }


            FASTAEntry& operator=(const FASTAEntry& rhs) = default;

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
        FASTAFile() = default;

        /// Destructor
        ~FASTAFile() override = default;

        /**
          @brief Prepares a FASTA file given by 'filename' for streamed reading using readNext().
          @exception Exception::FileNotFound is thrown if the file does not exists.
          @exception Exception::ParseError is thrown if the file does not suit to the standard.
        */
        void readStart(const String& filename);

        /// same as readStart(), but does internal progress logging whenever readNextWithProgress() is called
        void readStartWithProgress(const String& filename, const String& progress_label);

        /**
        @brief Reads the next FASTA entry from file.
        If you want to read all entries in one go, use load().
        @return true if entry was read; false if EOF was reached
        @exception Exception::FileNotFound is thrown if the file does not exists.
        @exception Exception::ParseError is thrown if the file does not suit to the standard.
        */
        bool readNext(FASTAEntry& protein);

        /// same as readNext(), but does internal progress logging; use readStartWithProgress() to enable this
        /// Calls progressEnd() when EOF is reached (i.e. when returning false)
        bool readNextWithProgress(FASTAEntry& protein);

        /// current stream position when reading a file
        std::streampos position();

        /// is stream at EOF?
        bool atEnd();

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
        @brief Closes the file (flush). Called implicitly when FASTAFile object goes out of scope.
        */
        void writeEnd();


        /**
          @brief loads a FASTA file given by 'filename' and stores the information in 'data'
          This uses more RAM than readStart() and readNext().
          @exception Exception::FileNotFound is thrown if the file does not exists.
          @exception Exception::ParseError is thrown if the file does not suit to the standard.
        */
        void load(const String& filename, std::vector<FASTAEntry>& data) const;

        /**
            @brief stores the data given by 'data' at the file 'filename'

            This uses more RAM than writeStart() and writeNext().
            @exception Exception::UnableToCreateFile is thrown if the process is not able to write the file.
          */
        void store(const String& filename, const std::vector<FASTAEntry>& data) const;

    protected:
        /**
         @brief Reads a protein entry from the current file position and returns the ID and sequence
         @return Return true if the protein entry was read and saved successfully, false otherwise
         */
        bool readEntry_(std::string& id, std::string& description, std::string& seq);

        std::fstream infile_;       ///< filestream for reading; init using FastaFile::readStart()
        std::ofstream outfile_;     ///< filestream for writing; init using FastaFile::writeStart()
        Size entries_read_{0};      ///< some internal book-keeping during reading
        std::streampos fileSize_{}; ///< total number of characters of filestream
        std::string seq_;           ///< sequence of currently read protein
        std::string id_;            ///< identifier of currently read protein
        std::string description_;   ///< description of currently read protein
    };

} // namespace OpenMS
