// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow, Nora Wild $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FASTAFile.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/CONCEPT/LogStream.h>

#include <algorithm>
#include <fstream>
#include <iterator>
#ifdef _OPENMP
#include <omp.h>
#endif


namespace OpenMS
{
  using namespace std;

  namespace
  {
    /// Parses the FASTA records in text[begin, end) as FASTAFile::readEntry_() does with the same
    /// characters. @p end is either the end of the file or the '>' of a record whose preceding line
    /// does not start with '>' -- a record the reader also starts there, as that line is part of a
    /// sequence. Returns false for anything readEntry_() would not read (the caller then reads the
    /// file with the reader, which reports the problem).
    bool parseFastaRecords(const std::string& text, size_t begin, size_t end, bool end_of_file,
                           std::vector<FASTAFile::FASTAEntry>& entries)
    {
      size_t pos = begin;
      while (pos < end)
      {
        if (text[pos] != '>') return false; // records start right at a '>' after the first one
        ++pos;
        FASTAFile::FASTAEntry entry;
        // ID: up to a space or tab after it, or the end of the line
        bool description_exists = true;
        bool reading = true;
        while (reading)
        {
          if (pos == end) return false;
          const char c = text[pos++];
          switch (c)
          {
            case ' ':
            case '\t':
              if (!entry.identifier.empty()) reading = false;
              break;
            case '\n':
              reading = false;
              description_exists = false;
              break;
            case '\r':
              break;
            default:
              entry.identifier += c;
          }
        }
        if (entry.identifier.empty()) return false;
        // description: the rest of the line
        reading = description_exists;
        while (reading)
        {
          if (pos == end) return false;
          const char c = text[pos++];
          switch (c)
          {
            case '\n':
              reading = false;
              break;
            case '\r':
            case '\t':
              break;
            default:
              entry.description += c;
          }
        }
        // sequence: up to a line that starts with '>'
        reading = true;
        while (reading)
        {
          if (pos == end)
          {
            // the end of the file -- or, for a chunk, a '>' after a sequence line
            if (!end_of_file && (pos == begin || text[pos - 1] != '\n')) return false;
            break;
          }
          const char c = text[pos++];
          switch (c)
          {
            case '\n':
              if (pos < end ? text[pos] == '>' : !end_of_file) reading = false;
              break;
            case '\r':
            case ' ':
            case '\t':
              break;
            default:
              entry.sequence += c;
          }
        }
        if (entry.sequence.empty()) return false;
        entries.push_back(std::move(entry));
      }
      return true;
    }

    /// True if the '>' at text[pos] starts a chunk: it starts a line, and the line before it, which
    /// starts at or after @p first, does not start with '>' -- so the reader starts a record there.
    bool isChunkStart(const std::string& text, size_t pos, size_t first)
    {
      if (pos <= first || text[pos] != '>' || text[pos - 1] != '\n') return false;
      const size_t line_end = pos - 1; // the '\n' ending the preceding line
      const size_t previous_newline = line_end == 0 ? std::string::npos : text.rfind('\n', line_end - 1);
      const size_t line_start = previous_newline == std::string::npos ? 0 : previous_newline + 1;
      if (line_start < first) return false;
      return line_start == line_end || text[line_start] != '>';
    }

    /// The first chunk start in text[from, end), or npos.
    size_t nextChunkStart(const std::string& text, size_t from, size_t first, size_t end)
    {
      for (size_t pos = text.find("\n>", from == 0 ? 0 : from - 1); pos != std::string::npos && pos + 1 < end;
           pos = text.find("\n>", pos + 1))
      {
        if (isChunkStart(text, pos + 1, first)) return pos + 1;
      }
      return std::string::npos;
    }

    /// The last chunk start in text(first, end), or npos.
    size_t lastChunkStart(const std::string& text, size_t first, size_t end)
    {
      for (size_t pos = end; pos >= first + 3;)
      {
        const size_t newline = text.rfind("\n>", pos - 2);
        if (newline == std::string::npos || newline + 1 <= first) break;
        if (isChunkStart(text, newline + 1, first)) return newline + 1;
        pos = newline + 1;
      }
      return std::string::npos;
    }

    /// Parses the records in text[first, end) in pieces on @p threads threads and appends them to
    /// @p data. @p end is the end of the file (@p end_of_file) or a chunk start.
    bool parseInParallel(const std::string& text, size_t first, size_t end, bool end_of_file, int threads,
                         std::vector<FASTAFile::FASTAEntry>& data)
    {
      std::vector<size_t> starts{first};
      for (int t = 1; t < threads; ++t)
      {
        const size_t pos = nextChunkStart(text, std::max(starts.back() + 1, first + (end - first) / threads * t), first, end);
        if (pos == std::string::npos) break;
        starts.push_back(pos);
      }
      starts.push_back(end);

      const int chunks = static_cast<int>(starts.size()) - 1;
      std::vector<std::vector<FASTAFile::FASTAEntry>> parts(chunks);
      std::vector<char> ok(chunks, 1);
      #pragma omp parallel for schedule(static, 1) num_threads(threads)
      for (int c = 0; c < chunks; ++c)
      {
        // no exception may leave the parallel region: a failed chunk sends the file to the reader
        try
        {
          ok[c] = parseFastaRecords(text, starts[c], starts[c + 1], end_of_file && c + 1 == chunks, parts[c]);
        }
        catch (...)
        {
          ok[c] = 0;
          std::vector<FASTAFile::FASTAEntry>().swap(parts[c]);
        }
      }
      if (std::find(ok.begin(), ok.end(), 0) != ok.end()) return false;
      size_t total = data.size();
      for (const auto& part : parts) total += part.size();
      if (total > data.capacity()) data.reserve(std::max(total, 2 * data.capacity())); // grows over the segments
      for (auto& part : parts)
      {
        std::move(part.begin(), part.end(), std::back_inserter(data));
      }
      return true;
    }

    /// Reads all records of @p filename on several threads. Returns false if the file is small,
    /// cannot be read here or holds anything the parser above declines; the caller then reads it
    /// with FASTAFile::readNext().
    bool loadFastaInParallel(const std::string& filename, std::vector<FASTAFile::FASTAEntry>& data)
    {
#ifdef _OPENMP
      const int threads = omp_get_max_threads();
#else
      const int threads = 1;
#endif
      constexpr size_t min_bytes = size_t(1) << 22; // smaller files are read quickly anyway
      // The file is read and parsed in segments of this size, so that little text is held besides the
      // entries (only a record longer than a segment is held whole).
      constexpr size_t segment_bytes = size_t(1) << 22;
      if (threads < 2 || !File::exists(filename) || !File::readable(filename)) return false;
      std::ifstream in(filename, std::ios::binary);
      if (!in) return false;
      in.seekg(0, std::ios::end);
      const std::streamoff size = in.tellg();
      if (size < static_cast<std::streamoff>(min_bytes)) return false;
      in.seekg(0, std::ios::beg);

      std::string text; // the start of a record (or of the file) and what follows it
      size_t first = 0;
      bool start_of_file = true;
      while (true)
      {
        const size_t kept = text.size();
        text.resize(kept + segment_bytes);
        in.read(text.data() + kept, static_cast<std::streamsize>(segment_bytes));
        text.resize(kept + static_cast<size_t>(in.gcount()));
        if (in.bad()) return false;
        const bool end_of_file = in.eof();

        if (start_of_file)
        {
          // as readStart(): skip leading whitespace and '#' (PEFF header) lines
          while (first < text.size())
          {
            const char c = text[first];
            if (c == '#')
            {
              const size_t newline = text.find('\n', first);
              first = newline == std::string::npos ? text.size() : newline + 1;
            }
            else if (c == ' ' || c == '\t' || c == '\n' || c == '\r') ++first;
            else break;
          }
          if (first == text.size() || text[first] != '>') return false;
          start_of_file = false;
        }

        // parse up to the last record that starts in the text; the rest is kept for the next segment
        size_t end = text.size();
        if (!end_of_file)
        {
          end = lastChunkStart(text, first, text.size());
          if (end == std::string::npos) continue; // no record starts after the first one yet: read on
        }
        if (!parseInParallel(text, first, end, end_of_file, threads, data)) return false;
        if (end_of_file) return true;
        text.erase(0, end);
        first = 0;
      }
    }
  }

  bool FASTAFile::readEntry_(std::string& id, std::string& description, std::string& seq)
  {
    std::streambuf* sb = infile_.rdbuf();
    bool keep_reading = true;
    bool description_exists = true;
    
  // Skip leading whitespace (including newlines, tabs, spaces)
  int c;
  while ((c = sb->sgetc()) != std::streambuf::traits_type::eof() && 
         (c == ' ' || c == '\t' || c == '\n' || c == '\r')) 
  {
    sb->sbumpc(); 
  }

  if (sb->sbumpc() != '>') 
  {
    return false;     
  }

    while (keep_reading)// reading the ID
    {
      int c = sb->sbumpc();// get and advance to next char
      switch (c)
      {
        case ' ':
        case '\t':
          if (!id.empty())
          {
            keep_reading = false; // ID finished
          }
          break;
        case '\n':                // ID finished and no description available
          keep_reading = false;
          description_exists = false;
          break;
        case '\r':
          break;
        case std::streambuf::traits_type::eof():
          infile_.setstate(std::ios::eofbit);
          return false;
        default:
          id += (char) c;
      }
    }

    if (id.empty())
    {
      return false;
    }
      

    if (description_exists)
    {
      keep_reading = true;
    }

    // reading the description
    while (keep_reading)       
    {
      int c = sb->sbumpc();    // get and advance to next char
      switch (c)
      {
        case '\n':             // description finished
          keep_reading = false;
          break;
        case '\r': // .. or
        case '\t':
          break;
        case std::streambuf::traits_type::eof():
          infile_.setstate(std::ios::eofbit);
          return false;
        default:
          description += (char) c;
      }
    }
    
    // reading the sequence
    keep_reading = true;
    while (keep_reading)
    {
      int c = sb->sbumpc(); // get and advance to next char
      switch (c)
      {
        case '\n':
          if (sb->sgetc() == '>')// reaching the beginning of the next protein-entry
          {
            keep_reading = false;
          }
          break;
        case '\r': // not saving white spaces
        case ' ': 
        case '\t':
          break;
        case std::streambuf::traits_type::eof():
          infile_.setstate(std::ios::eofbit);
          return !seq.empty();
        default:
          seq += (char) c;
      }
    }
    return !seq.empty();
  }

  void FASTAFile::readStart(const std::string& filename)
  {

    if (!File::exists(filename))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    if (!File::readable(filename))
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    if (infile_.is_open()) infile_.close(); // precaution

    infile_.open(filename.c_str(), std::ios::binary | std::ios::in);
    infile_.seekg(0, infile_.end);
    fileSize_ = infile_.tellg();
    infile_.seekg(0, infile_.beg);

    std::streambuf *sb = infile_.rdbuf();
   // Skip leading whitespace and PEFF headers
    int c;
    while ((c = sb->sgetc()) != 
    std::streambuf::traits_type::eof())
  {
    if (c == '#') 
    {
      infile_.ignore(numeric_limits<streamsize>::max(), '\n');
    }
    else if (c == ' ' || c == '\t' || c == '\n' || c == '\r') 
    {
      sb->sbumpc(); 
    }
    else 
    {
      break;
    }
  }
  entries_read_ = 0;
  }

  void FASTAFile::readStartWithProgress(const std::string& filename, const std::string& progress_label)
  {
    readStart(filename);
    startProgress(0, fileSize_, progress_label);
  }

  bool FASTAFile::readNext(FASTAEntry &protein)
  {
    if (infile_.eof())
    {
      return false;
    }

    seq_.clear(); // Note: it is fine to clear() after std::move as it will turn the "valid but unspecified state" into a specified (the empty) one
    id_.clear();
    description_.clear();

    if (!readEntry_(id_, description_, seq_))
    {
      if (entries_read_ == 0)
      {
        seq_ = "The first entry could not be read!";
      }
      else
      {
        seq_ = "Only " + StringUtils::toStr(entries_read_) + " proteins could be read. Parsing next record failed.";
      }
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "",
                                  "Error while parsing FASTA file! " + seq_ + " Please check the file!");
    }
    ++entries_read_;

    protein.identifier = std::move(id_);
    protein.description = std::move(description_);
    protein.sequence = std::move(seq_);

    setProgress(infile_.tellg());

    return true;
  }

  bool FASTAFile::readNextWithProgress(FASTAEntry& protein)
  {
    if (readNext(protein))
    {
      setProgress(position());
      return true;
    }
    else
    {
      endProgress(); 
      return false;
    }
  }

  std::streampos FASTAFile::position()
  {
    return infile_.tellg();
  }

  bool FASTAFile::setPosition(const std::streampos &pos)
  {
    if (pos <= fileSize_)
    {
      infile_.clear(); // when end of file is reached, otherwise it gets -1
      infile_.seekg(pos);
      return true;
    }
    return false;
  }

  bool FASTAFile::atEnd()
  {
    return (infile_.peek() == std::streambuf::traits_type::eof());
  }

  void FASTAFile::load(const std::string &filename, vector<FASTAEntry> &data) const
  {
    startProgress(0, 1, "Loading FASTA file");
    data.clear();
    // Large files are parsed in pieces on several threads, with the same result.
    if (loadFastaInParallel(filename, data))
    {
      endProgress();
      return;
    }
    data.clear();
    FASTAEntry p;
    FASTAFile f;
    f.readStart(filename);
    while (f.readNext(p))
    {
      data.push_back(std::move(p));
    }
    endProgress();
  }

  void FASTAFile::writeStart(const std::string &filename)
  {
    if (!FileHandler::hasValidExtension(filename, FileTypes::FASTA))
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
                                          "invalid file extension; expected '" +
                                          FileTypes::typeToName(FileTypes::FASTA) + "'");
    }

    outfile_.open(filename.c_str(), ofstream::out);

    if (!outfile_.good())
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
  }

  void FASTAFile::writeNext(const FASTAEntry &protein)
  {
    outfile_ << '>' << protein.identifier << ' ' << protein.description << "\n";
    const std::string &tmp(protein.sequence);

    int chunks(tmp.size() / 80); // number of complete chunks
    Size chunk_pos(0);
    while (--chunks >= 0)
    {
      outfile_.write(&tmp[chunk_pos], 80);
      outfile_ << "\n";
      chunk_pos += 80;
    }

    if (tmp.size() > chunk_pos)
    {
      outfile_.write(&tmp[chunk_pos], tmp.size() - chunk_pos);
      outfile_ << "\n";
    }
  }

  void FASTAFile::writeEnd()
  {
    outfile_.close();
  }

  void FASTAFile::store(const std::string &filename, const vector<FASTAEntry> &data) const
  {
    startProgress(0, data.size(), "Writing FASTA file");
    FASTAFile f;
    f.writeStart(filename);
    for (const FASTAFile::FASTAEntry& it : data)
    {
      f.writeNext(it);
      nextProgress();
    }
    f.writeEnd(); // close file
    endProgress();
  }

}// namespace OpenMS