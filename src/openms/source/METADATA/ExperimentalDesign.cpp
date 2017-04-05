// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
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
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/FORMAT/TextFile.h>

using namespace std;

namespace OpenMS
{
  void ExperimentalDesign::load(const string & tsv_file, ExperimentalDesign & design) const
  {
    TextFile tf(tsv_file, true);

    bool header_parsed(false);

    ExperimentalDesign ed;

    Size line_number = 0;

    for (TextFile::ConstIterator sit = tf.begin(); sit != tf.end(); ++sit)
    {
      String s = *sit;

      // skip empty lines
      if (s.trim().empty()) continue;

      // run-level header
      if (line_number == 0)
      {
        if (!s.hasPrefix("Run"))
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, tsv_file,
            "Error: Line does not contain the run header of the experimental design: " + String(s) + ".");
        }
        else
        {
          StringList cells;
          s.split("\t", cells);
          if (s.size() != 4)
          {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, tsv_file,
              "Error: Wrong number of columns in the experimental design header provided: " + String(s) + ".");
          }
          header_parsed = true;
          ++line_number;
          continue;
        }
      }

      // run-level rows
      StringList cells;
      s.split("\t", cells);

      if (s.size() < 4)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, tsv_file,
          "Error: Not all columns of experimental design provided: " + String(s) + ".");
      }

      MSRun r;
      r.file = cells[1];
      r.fraction = static_cast<unsigned>(cells[2].toInt());
      r.technical_replicate = static_cast<unsigned>(cells[3].toInt());
      
      ed.runs.push_back(r);

      ++line_number;
    }
  }

  map<unsigned, set<unsigned> > ExperimentalDesign::getFractionToRunsMapping() const
  {
    map<unsigned, set<unsigned> > ret;

    for (Size i = 0; i != runs.size(); ++i)
    {
      ret[runs[i].fraction].insert(i + 1);
    }

    return ret;
  }
}
