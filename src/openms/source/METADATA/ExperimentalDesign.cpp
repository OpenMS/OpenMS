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

using namespace std;

namespace OpenMS
{
  // static
  ExperimentalDesign ExperimentalDesign::load(const String & tsv_file)
  {
    ExperimentalDesign design;

    design.run_section_.clear();

    TextFile tf(tsv_file, true);

    int line_number(0);
    for (String s : tf)
    {
      ++line_number;      
      // skip empty lines
      if (s.trim().empty()) { continue; }

      // run-level header
      if (line_number == 1)
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
          continue;
        }
      }

      // run-level run_section_
      StringList cells;
      s.split("\t", cells);

      if (cells.size() < 3)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, tsv_file,
          "Error: Not all columns of experimental design provided: " + String(s) + ".");
      }

      RunRow row;

      row.run = cells[0].toInt();
      row.fraction = cells[1].toInt();

      // if sample and channel column missing, assign run index
      if (cells.size() <= 3)
      {
        row.channel = 1;  
        row.sample = row.run; 
      }
      else if (cells.size() <= 4) // at least channel given
      {
        row.channel = cells[3].toInt();

        if (cells.size() == 5)
        {
          row.sample = cells[4].toInt();;
        }
        else // assign sample index to channel number
        {
          row.sample = row.channel;
        }
      }

      String spec_file = cells[2];
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
          row.path = design_file_relative.toStdString();
        }
        else
        {
          // check current folder
          String f = File::absolutePath(spec_file);
          if (File::exists(f))
          {
            row.path = f;
          }
        }
      }     
      else
      {
        // set to absolute path
        row.path = spec_file;
      }

      if (!File::exists(row.path))
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, tsv_file,
          "Error: Spectra file does not exist: '" + String(row.path) + "'");
      }

      design.run_section_.push_back(row);
    }

    design.sort_();
    design.checkValidRunSection_();
  }
  
  ExperimentalDesign ExperimentalDesign::fromConsensusMap(const ConsensusMap& cm)
  {
    ExperimentalDesign ed;
    // path of the original MS run (mzML / raw file)
    StringList ms_run_paths;
    cm.getPrimaryMSRunPath(ms_run_paths);

    // no fractionation -> as many runs as samples
    // each consensus element corresponds to one sample abundance
    size_t sample(1);
    ExperimentalDesign::RunRows rows;
    for (auto const & f : ms_run_paths)
    {
      ExperimentalDesign::RunRow r;
      r.path = f;
      r.fraction = 1;
      r.sample = sample;
      r.run = sample;
      r.channel = 1; // TODO MULTIPLEXING: adapt for non-label-free
      rows.push_back(r);
      ++sample;
    }
    ed.setRunSection(rows);
    LOG_INFO << "Experimental design (ConsensusMap derived):\n" 
       << "  files: " << ed.getNumberOfMSFiles()
       << "  fractions: " << ed.getNumberOfFractions()
       << "  channels: " << ed.getNumberOfChannels()
       << "  samples: " << ed.getNumberOfSamples() << "\n"  
       << endl;
    return ed; 
  }

  ExperimentalDesign ExperimentalDesign::fromFeatureMap(const FeatureMap& fm)
  {
    ExperimentalDesign ed;
    // path of the original MS run (mzML / raw file)
    StringList ms_paths;
    fm.getPrimaryMSRunPath(ms_paths);
  
    if (ms_paths.size() != 1)
    {
        throw Exception::MissingInformation(
         __FILE__, 
         __LINE__, 
         OPENMS_PRETTY_FUNCTION, 
         "FeatureMap annotated with " + String(ms_paths.size()) + " MS files. Must be exactly one.");
    }
    
    // Feature map is simple. One file, one fraction, one sample, one run
    ExperimentalDesign::RunRow r;
    r.path = ms_paths[0];
    r.fraction = 1;
    r.sample = 1;
    r.run = 1;
    r.channel = 1;

    ExperimentalDesign::RunRows rows(1, r);
    ed.setRunSection(rows);
    LOG_INFO << "Experimental design (FeatureMap derived):\n" 
       << "  files: " << ed.getNumberOfMSFiles()
       << "  fractions: " << ed.getNumberOfFractions()
       << "  channels: " << ed.getNumberOfChannels()
       << "  samples: " << ed.getNumberOfSamples() << "\n"  
       << endl;
    return ed; 
  }

  ExperimentalDesign ExperimentalDesign::fromIdentifications(const vector<ProteinIdentification> & proteins)
  {
    ExperimentalDesign ed;
    // path of the original MS files (mzML / raw file)
    StringList ms_run_paths;
    for (auto const & p : proteins)
    {
      StringList tmp_ms_run_paths;
      p.getPrimaryMSRunPath(tmp_ms_run_paths);
      if (tmp_ms_run_paths.size() != 1)
      {
        throw Exception::MissingInformation(
         __FILE__, 
         __LINE__, 
         OPENMS_PRETTY_FUNCTION, 
         "ProteinIdentification annotated with " + String(tmp_ms_run_paths.size()) + " MS files. Must be exactly one.");
      }
      ms_run_paths.push_back(tmp_ms_run_paths[0]);
    }

    // no fractionation -> as many runs as samples
    // each identification run corresponds to one sample abundance
    size_t sample(1);
    ExperimentalDesign::RunRows rows;
    for (auto const & f : ms_run_paths)
    {
      ExperimentalDesign::RunRow r;
      r.path = f;
      r.fraction = 1;
      r.sample = sample;
      r.run = sample;
      r.channel = 1;

      rows.push_back(r);
      ++sample;
    }
    ed.setRunSection(rows);
    LOG_INFO << "Experimental design (Identification derived):\n" 
       << "  files: " << ed.getNumberOfMSFiles()
       << "  fractions: " << ed.getNumberOfFractions()
       << "  channels: " << ed.getNumberOfChannels()
       << "  samples: " << ed.getNumberOfSamples() << "\n"  
       << endl;
    return ed;
  }


  map<unsigned, vector<String> > ExperimentalDesign::getFractionToMSFilesMapping() const
  {
    map<unsigned, vector<String> > ret;

    for (RunRow const & r : run_section_)
    {
      ret[r.fraction].emplace_back(r.path);
    }

    return ret;
  }

  bool ExperimentalDesign::sameNrOfMSFilesPerFraction() const
  {
    map<unsigned, vector<String>> frac2files = getFractionToMSFilesMapping();
    if (frac2files.size() <= 1) { return true; }
 
    Size files_per_fraction(0);
    for (auto const & f : frac2files)
    {
      if (files_per_fraction == 0) // first fraction, initialize
      {
        files_per_fraction = f.second.size();
      }
      else // fraction >= 2
      {
        // different number of associated MS files?
        if (f.second.size() != files_per_fraction)
        {
          return false;
        }
      }
    }
    return true;
  }
}

