// --------------------------------------------------------------------------
//           OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Nico Pfeifer, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/CONCEPT/LogStream.h>


namespace OpenMS
{
    using namespace std;

    FASTAFile::FASTAFile()
            : entries_read_(0)
    {
    }

    FASTAFile::~FASTAFile()
    {
        // infile_ and outfile_ will close automatically when going out of scope. No need to do it explicitly here.
    }

    bool FASTAFile::readEntry_(std::string & id, std::string & seq)
    {
        std::streambuf* sb = infile_.rdbuf();
        bool keep_reading = true;


        if(sb->sbumpc() == '>')//not saving '>'
        {
            while(keep_reading)// reading the ID
            {
                int c = sb->sbumpc();// get and advance to next char
                switch (c)
                {
                    case '\n'://ID finished
                        keep_reading = false;
                        break;
                    case '\r':
                        if (sb->sgetc() == '\n') // peek current char
                        {
                            sb->sbumpc(); // consume it
                        }
                        break;
                    case std::streambuf::traits_type::eof():
                        infile_.setstate(std::ios::eofbit);
                        if (id.empty())
                        { // only if we just started a new line, we set the is.fail() == true, ie. is == false
                            infile_.setstate(std::ios::badbit);
                            return false;
                        }
                        return true;
                    default:
                        id += (char)c;
                }
            }
        }
        else return false;//was in wrong position for reading ID
        if(id.empty()==true) return false;

        keep_reading = true;
        while(keep_reading)//reading the sequence
        {
            int c = sb->sbumpc();// get and advance to next char
            switch (c)
            {
                case '\n':
                    if (sb->sgetc() == '>') //reaching the beginning of the next protein-entry
                    {
                        keep_reading = false;
                    }
                    break;
                case '\r':
                    if (sb->sgetc() == '\n') // peek current char
                    {
                        sb->sbumpc(); // consume it
                    }
                    break;
                case ' ': //not saving white spaces
                    break;
                case '\t':
                    break;
                case std::streambuf::traits_type::eof():
                    infile_.setstate(std::ios::eofbit);
                    if (seq.empty())
                    { // only if we just started a new line, we set the is.fail() == true, ie. is == false
                        infile_.setstate(std::ios::badbit);
                        return false;
                    }
                    return true;
                default:
                    seq += (char)c;
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

        std::streambuf* sb = infile_.rdbuf();
        while(sb->sgetc() == '#')// Skip the header of PEFF files (http://www.psidev.info/peff)
        {
            infile_.ignore(numeric_limits<streamsize>::max(),'\n');
        }
        entries_read_ = 0;
    }

    bool FASTAFile::readNext(FASTAEntry& protein)
    {
        if (infile_.eof())
        {
            return false;
        }

        String id, s;
        if (readEntry_(id, s) == false)
        {
            if (entries_read_ == 0) s = "The first entry could not be read!";
            else s = "Only " + String(entries_read_) + " proteins could be read. The record after failed.";
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", "Error while parsing FASTA file! " + s + " Please check the file!");
        }
        ++entries_read_;

        protein.sequence = s; // assign here, since 's' might have higher capacity, thus wasting memory (usually 10-15%)
        // handle id
        id.trim();
        String::size_type position = id.find_first_of(" \v\t");
        if (position == String::npos)
        {
            protein.identifier = std::move(id);
            protein.description = "";
        }
        else
        {
            protein.identifier = id.substr(0, position);
            protein.description = id.suffix(id.size() - position - 1);
        }
        return true;
    }

    std::streampos FASTAFile::position()
    {
        return infile_.tellg();
    }


    bool FASTAFile::setPosition(const std::streampos& pos)
    {
        if(pos <= fileSize_)
        {
            infile_.clear();//when end of file is reached, otherwise it gets -1
            infile_.seekg(pos);
            return true;
        }
        return false;
    }

    bool FASTAFile::atEnd() 
    {
        return (infile_.peek() == std::streambuf::traits_type::eof());
    }

    void FASTAFile::load(const String& filename, vector<FASTAEntry>& data) const
    {
        data.clear();
        FASTAEntry p;
        FASTAFile f;
        f.readStart(filename);
        while (f.readNext(p))
        {
            data.push_back(std::move(p));
        }
        return;
    }

    void FASTAFile::writeStart(const String& filename)
    {
        if (!FileHandler::hasValidExtension(filename, FileTypes::FASTA))
        {
            throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "invalid file extension; expected '" + FileTypes::typeToName(FileTypes::FASTA) + "'");
        }

        outfile_.open(filename.c_str(), ofstream::out);

        if (!outfile_.good())
        {
            throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
        }
    }

    void FASTAFile::writeNext(const FASTAEntry& protein)
    {
        outfile_ << ">" << protein.identifier << " " << protein.description << "\n";
        const String& tmp(protein.sequence);

        int chunks( tmp.size()/80 ); // number of complete chunks
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

    void FASTAFile::store(const String& filename, const vector<FASTAEntry>& data) const
    {
        FASTAFile f;
        f.writeStart(filename);
        for (vector<FASTAEntry>::const_iterator it = data.begin(); it != data.end(); ++it)
        {
            f.writeNext(*it);
        }
        f.writeEnd(); // close file
    }

} // namespace OpenMS
