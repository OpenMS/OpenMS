// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <QtCore/QString>
#include <QtCore/QFileInfo>

#include <iostream>

using namespace std;

namespace OpenMS
{
  void ExperimentalDesign::load(const String & tsv_file, ExperimentalDesign & design) const
  {
    design.runs.clear();

    TextFile tf(tsv_file, true);

    int line_number = 0;

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
          if (cells.size() != 4)
          {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, tsv_file,
              "Error: Wrong number of columns (" + String(cells.size()) + ") in the experimental design header provided: " + String(s) + ".");
          }
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
      int run_number = cells[0].toInt();

      // read spectra file name
      String spec_file = cells[1];
      QFileInfo spectra_file_info(spec_file.toQString());

      if (spectra_file_info.isRelative())
      {
        // file name is relative so we need to figure out the correct folder

        // first check folder relative to folder of design file 
        // to allow, for example, a design in ./design.tsv and spectra in ./spectra/a.mzML
        // where ./ is the same folder
        QFileInfo design_file_info(tsv_file.toQString());
        QString design_file_relative(design_file_info.absolutePath());
        design_file_relative = design_file_relative + "/" + spec_file.toQString();

        if (File::exists(design_file_relative))
        {
          r.file = design_file_relative.toStdString();
        }
        else
        {
          // check current folder
          String f = File::absolutePath(spec_file);
          if (File::exists(f))
          {
            r.file = f;
          }
        }
      }     
      else
      {
        // set to absolute path
        r.file = spec_file;
      }

      if (!File::exists(r.file))
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, tsv_file,
          "Error: Spectra file does not exist: '" + String(r.file) + "'");
      }

      r.fraction = static_cast<unsigned>(cells[2].toInt());
      r.technical_replicate = static_cast<unsigned>(cells[3].toInt());
      
      design.runs.push_back(r);

      // validation: check if run number matches the line number in the design file
      if (line_number != run_number)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, tsv_file,
          "Error: Run index (" + String(run_number) + ") does not match row index (" + String(line_number) + ") in line: " + String(s) + ".");
      }

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

  bool ExperimentalDesign::sameNrOfRunsPerFraction() const
  {
    map<unsigned, set<unsigned> > frac2run = getFractionToRunsMapping();
    if (frac2run.size() > 1)
    {
      Size runs_per_fraction(0);
      for (map<unsigned, set<unsigned> >::const_iterator mit = frac2run.begin(); mit != frac2run.end(); ++mit)
      {
        if (mit == frac2run.begin()) // fraction 1
        {
          // initialize runs per fraction with the number of runs in the first fraction
          runs_per_fraction = mit->second.size();
        }
        else // fraction >= 2
        {
          if (mit->second.size() != runs_per_fraction)
          {
            return false;
          }
        }
      }
    }
    return true;
  }
}

