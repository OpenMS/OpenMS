// --------------------------------------------------------------------------
//                       OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Authors: Marie Hoffmann $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FASTQFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/seq_io/guess_stream_format.h>
#include <seqan/seq_io/read_fasta_fastq.h>
#include <seqan/sequence.h>

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <ios>


namespace OpenMS
{

  FASTQFile::FASTQFile()
  {
  }

  FASTQFile::~FASTQFile()
  {
  }

  void FASTQFile::load(const String& filename, std::vector<FASTQEntry>& data)
  {
    data.clear();

    if (!File::exists(filename))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    if (!File::readable(filename))
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    std::fstream in(filename.c_str(), std::ios::binary | std::ios::in);
    seqan::RecordReader<std::fstream, seqan::SinglePass<> > reader(in);
    seqan::CharString id;
    seqan::CharString seq;
    seqan::CharString qual;
    String::size_type position = String::npos;
    Size size_read(0);

    while (!atEnd(reader))
    {
      if (readRecord(id, seq, qual, reader, seqan::Fastq()) != 0)
      {
        String msg;
        if (data.empty())
          msg = "The first entry could not be read!";
        else
          msg = "The last successful FASTQ record was: '@" + std::string(toCString(data.back().identifier)) +
                  "'. The record after failed.";
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", "Error while parsing FASTQ file '" +
             filename + "'! " + msg +  " Please check the file!");
      }

      FASTQEntry newEntry;
      newEntry.sequence = seq;

      // handle id
      String id_tmp = String(toCString(id)).trim();

      position = id_tmp.find_first_of(" \v\t");
      if (position == String::npos)
      {
        newEntry.identifier = seqan::CharString(id_tmp.c_str());
        newEntry.description = "";
      } else
      {
        newEntry.identifier = seqan::CharString(id_tmp.substr(0, position).c_str());
        newEntry.description = seqan::CharString(id_tmp.suffix(id_tmp.size() - position - 1).c_str());
      }

      // handle quality
      newEntry.quality = seqan::CharString(qual);
      id_tmp.clear();
      seqan::clear(id);
      seqan::clear(seq);
      seqan::clear(qual);
      data.push_back(newEntry);
      size_read += seqan::length(newEntry.sequence);
    }

    in.close();

    if (size_read > 0 && data.empty())
    {
      LOG_WARN << "No entries from FASTQ file read. Does the file have MacOS line endings? "
      << "Convert to Unix or Windows line endings to fix!" << std::endl;
    }
    return;
  }


  void FASTQFile::store(const String& filename, const std::vector<FASTQEntry>& data) const
  {
    std::ofstream outfile;
    outfile.open(filename.c_str(), std::ofstream::out);
    if (!outfile.good())
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    for (std::vector<FASTQEntry>::const_iterator it = data.begin(); it != data.end(); ++it)
    {
      outfile << "@" << it->identifier << " " << it->description << "\n";
      outfile << it->sequence << "\n+\n" << it->quality << "\n";
    }
    outfile.close();
  }

}  // namespace OpenMS
