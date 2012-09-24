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
// $Authors: Nico PFeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <fstream>

using namespace std;

namespace OpenMS
{
  FASTAFile::FASTAFile()
  {

  }

  FASTAFile::~FASTAFile()
  {

  }

  void FASTAFile::load(const String & filename, vector<FASTAEntry> & data)
  {
    String temp = "";
    string::size_type position = String::npos;

    data.clear();

    if (!File::exists(filename))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }

    if (!File::readable(filename))
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }

    String tag, seq;

    ifstream in(filename.c_str());
    String line;
    Size size_read(0);
    while (getline(in, line, '\n'))
    {
      size_read += line.size();
      if (line.size() > 0)
      {
        if (line[0] == '>')
        {
          if (tag != "" && seq != "")
          {
            FASTAEntry entry;
            position = tag.find_first_of(" \v\t");
            if (position == String::npos)
            {
              entry.identifier = tag;
              entry.description = "";
            }
            else
            {
              entry.identifier = tag.substr(0, position);
              entry.description = tag.suffix(tag.size() - position - 1);
            }
            entry.sequence = seq;
            entry.sequence.removeWhitespaces();
            data.push_back(entry);
            tag = "";
            seq = "";
          }

          line.erase(line.begin());
          tag = line.trim();
        }
        else
        {
          seq += line.trim();
        }
      }
    }

    if (tag != "" && seq != "")
    {
      FASTAEntry entry;
      // allow only " \t" and nothing else, especially not "|"=="\v",
      // since NCBI or Ensemble use >ENSG00000203832|ENST00000430395|ENSP00000412476|NBPF20
      //                         or >gi|129295|sp|P01013|OVAX_CHICK
      // as Identifier! (with description following)
      position = tag.find_first_of(" \t");
      if (position == String::npos)
      {
        entry.identifier = tag;
        entry.description = "";
      }
      else
      {
        entry.identifier = tag.substr(0, position);
        entry.description = tag.suffix(tag.size() - position - 1);
      }
      entry.sequence = seq;
      data.push_back(entry);
    }
    in.close();

    if (size_read > 0 && data.empty())
      LOG_WARN << "No entries from FASTA file read. Does the file have MacOS line endings? Convert to Unix or Windows line endings to fix!" << std::endl;

    return;
  }

  void FASTAFile::store(const String & filename, const vector<FASTAEntry> & data) const
  {
    ofstream outfile;
    outfile.open(filename.c_str(), ofstream::out);

    if (!outfile.good())
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }

    for (vector<FASTAEntry>::const_iterator it = data.begin(); it != data.end(); ++it)
    {

      outfile << ">" << it->identifier << " " << it->description << "\n";

      String tmp(it->sequence);
      while (tmp.size() > 80)       // surprisingly fast, even though its using erase(). For-loop with substr() is much SLOWER!
      {
        outfile << tmp.prefix(80) << "\n";
        tmp.erase(0, 80);
      }

      if (tmp.size() > 0)
      {
        outfile << tmp << "\n";
      }
    }
    outfile.close();
  }

} // namespace OpenMS
