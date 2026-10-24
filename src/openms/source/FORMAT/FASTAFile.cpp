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


namespace OpenMS
{
  using namespace std;

  bool FASTAFile::readEntry_(std::string& id, std::string& description, std::string& seq)
  {
    std::streambuf* sb = infile_.rdbuf();
    constexpr int eof = std::streambuf::traits_type::eof();
    // the next character, counting the lines (if known, see setPosition())
    const auto next = [this, sb]()
    {
      const int c = sb->sbumpc();
      if (c == '\n' && line_ != 0) ++line_;
      return c;
    };

    // Skip leading whitespace (including newlines, tabs, spaces)
    int c;
    while ((c = sb->sgetc()) != eof && (c == ' ' || c == '\t' || c == '\n' || c == '\r'))
    {
      next();
    }

    entry_line_ = line_;
    entry_offset_ = infile_.tellg();
    c = next();
    if (c != '>')
    {
      if (c == eof)
      {
        infile_.setstate(std::ios::eofbit);
        error_ = "no FASTA entry found";
        return false;
      }
      std::string text(1, static_cast<char>(c));
      while (text.size() < 30 && (c = sb->sgetc()) != eof && c != '\n' && c != '\r') { text += static_cast<char>(next()); }
      error_ = "expected a header line starting with '>', found '" + text + "'";
      return false;
    }

    // reading the ID
    bool keep_reading = true;
    bool description_exists = true;
    while (keep_reading)
    {
      c = next();
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
        case eof:                 // a header without sequence at the end of the file
          infile_.setstate(std::ios::eofbit);
          keep_reading = false;
          description_exists = false;
          break;
        default:
          id += (char) c;
      }
    }

    if (id.empty())
    {
      error_ = "the header line has no identifier after '>'";
      return false;
    }

    // reading the description
    keep_reading = description_exists;
    while (keep_reading)
    {
      c = next();
      switch (c)
      {
        case '\n':             // description finished
          keep_reading = false;
          break;
        case '\r': // .. or
        case '\t':
          break;
        case eof:              // a header without sequence at the end of the file
          infile_.setstate(std::ios::eofbit);
          keep_reading = false;
          break;
        default:
          description += (char) c;
      }
    }

    // reading the sequence, up to a line that starts with '>' (the next entry): a header line directly followed by
    // the next header is an entry without sequence
    keep_reading = !infile_.eof() && sb->sgetc() != '>';
    while (keep_reading)
    {
      c = next();
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
        case eof:
          infile_.setstate(std::ios::eofbit);
          keep_reading = false;
          break;
        default:
          seq += (char) c;
      }
    }
    return true;
  }

  std::string FASTAFile::errorLocation_() const
  {
    return entry_line_ != 0 ? "line " + StringUtils::toStr(entry_line_)
                            : "at byte " + StringUtils::toStr(static_cast<Size>(static_cast<std::streamoff>(entry_offset_)));
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
    filename_ = filename;
    line_ = 1;

    std::streambuf *sb = infile_.rdbuf();
    // Skip a UTF-8 byte order mark (as written by some Windows editors)
    {
      char bom[3] = {0, 0, 0};
      const std::streamsize n = sb->sgetn(bom, 3);
      if (!(n == 3 && bom[0] == '\xEF' && bom[1] == '\xBB' && bom[2] == '\xBF'))
      {
        sb->pubseekpos(0, std::ios::in);
      }
    }
    // Skip leading whitespace and PEFF headers
    int c;
    while ((c = sb->sgetc()) != std::streambuf::traits_type::eof())
    {
      if (c == '#')
      {
        infile_.ignore(numeric_limits<streamsize>::max(), '\n');
        ++line_;
      }
      else if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
      {
        if (sb->sbumpc() == '\n') ++line_;
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
      std::string read = entries_read_ == 0 ? std::string("The first entry could not be read")
                                            : "Only " + StringUtils::toStr(entries_read_) + " proteins could be read";
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename_,
                                  "Error while parsing FASTA file, " + errorLocation_() + ": " + error_ + ". " + read
                                  + ". Please check the file!");
    }
    if (seq_.empty())
    {
      OPENMS_LOG_WARN << "Warning: the FASTA entry '" << id_ << "' (" << filename_ << ", " << errorLocation_()
                      << ") has no sequence." << std::endl;
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
      line_ = 0; // unknown: errors name the byte offset instead
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