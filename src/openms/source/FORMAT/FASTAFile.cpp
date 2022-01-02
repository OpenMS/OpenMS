// --------------------------------------------------------------------------
//           OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//  notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//  notice, this list of conditions and the following disclaimer in the
//  documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//  may be used to endorse or promote products derived from this software
//  without specific prior written permission.
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
    bool keep_reading = true;
    bool description_exists = true;

    if (sb->sbumpc() != '>')
    {
      return false;     // was in wrong position for reading ID
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
          if (seq.empty())
          {
            infile_.setstate(std::ios::badbit);
            return false;
          }
          return true;
        default:
          seq += (char) c;
      }
    }
    return !seq.empty();
  }

  void FASTAFile::readStart(const String& filename)
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
    while (sb->sgetc() == '#') // Skip the header of PEFF files (http://www.psidev.info/peff)
    {
      infile_.ignore(numeric_limits<streamsize>::max(), '\n');
    }
    entries_read_ = 0;
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
        seq_ = "Only " + String(entries_read_) + " proteins could be read. Parsing next record failed.";
      }
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "",
                                  "Error while parsing FASTA file! " + seq_ + " Please check the file!");
    }
    ++entries_read_;

    protein.identifier = std::move(id_);
    protein.description = std::move(description_);
    protein.sequence = std::move(seq_);

    return true;
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

  void FASTAFile::load(const String &filename, vector<FASTAEntry> &data) const
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

  void FASTAFile::writeStart(const String &filename)
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
    outfile_ << ">" << protein.identifier << " " << protein.description << "\n";
    const String &tmp(protein.sequence);

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

  void FASTAFile::store(const String &filename, const vector<FASTAEntry> &data) const
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

} // namespace OpenMS
