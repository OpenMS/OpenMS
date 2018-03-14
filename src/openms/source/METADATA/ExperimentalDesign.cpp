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

  // Parse Error of filename if test holds
  void parseErrorIf(const bool test, const String &filename, const String &message)
  {
    if (test)
    {
      throw Exception::ParseError(
          __FILE__,
          __LINE__,
           OPENMS_PRETTY_FUNCTION,
          filename,
          "Error: " + message);
    }
  }

  // Trims whitespace around all strings in vector
  void trimAll(std::vector<String> &vec)
  {
    for (Size i = 0; i < vec.size(); ++i)
    {
        vec[i] = vec[i].trim();
    }
  }

  inline bool contains(const StringList &vec, const String &item)
  {
    auto e = std::end(vec);
    return std::find(std::begin(vec), e, item) != e;
  }

  void parse_header( const StringList &header,
                     const String &filename,
                           std::map< String, Size> &column_map,
                     const std::set< String > &required,
                     const std::set< String > &optional,
                     const bool allow_other_header)
  {
    // Headers as set
    std::set< String > header_set(header.begin(), header.end());
    parseErrorIf(header_set.size() !=  header.size(), filename,
      "Some column headers of the table appear multiple times!");

    // Check that all required headers are there
    for (const String &req_header : required)
    {
      parseErrorIf( ! contains(header, req_header), filename,
        "Missing column header: " + req_header);
    }
    // Assign index in column map and check for weird headers
    for (Size i = 0; i < header.size(); ++i)
    {
      const String &h = header[i];
      const bool header_unexpected = required.find(h) == required.end() && optional.find(h) == optional.end();
      parseErrorIf(
        ! allow_other_header && header_unexpected,
        filename,
        "Header not allowed in this section of the Experimental Design: " + h
      );
      column_map[h] = i;
    }
  }

  String findSpectraFile(const String &spec_file, const String &tsv_file)
  {
    String result;
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
        result = design_file_relative.toStdString();
      }
      else
      {
        // check current folder
        String f = File::absolutePath(spec_file);
        if (File::exists(f))
        {
          result = f;
        }
      }

      // if result still empty, just use the provided value
      if (result.empty())
      {
        result = spec_file;
      }
    }
    else
    {
      // set to absolute path
      result = spec_file;
    }

    /* Currently, existence of file is not required
    if (!File::exists(result))
    {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, tsv_file,
              "Error: Spectra file does not exist: '" + result + "'");
    }
    */
    return result;
  }


  // static
  ExperimentalDesign ExperimentalDesign::load(const String & tsv_file)
  {
    ExperimentalDesign design;
    design.run_section_.clear();

    bool has_sample(false);
    bool has_channel(false);

    // Maps the column header string to the column index for
    // the run section
    std::map< String, Size> _run_column_header_to_index;

    unsigned line_number(0);
    int state(0);
    Size n_col = 0;

    const TextFile text_file(tsv_file, true);
    for (String s : text_file)
    {
      // skip empty lines (except in state 1, where the sample table is read)
      const String line(s.trim());

      if (line.empty() && state != 1) {
        continue;
      }

      // Now split the line into individual cells
      StringList cells;
      line.split("\t", cells);
      trimAll(cells);

      // State = 0 : Run Section header
      if (state == 0)
      {
        state = 1;
        parse_header(
          cells,
          tsv_file,
          _run_column_header_to_index,
          {"Run", "Fraction", "Path(Spectra File)"},
          {"Channel", "Sample"}, false
        );
        has_channel = _run_column_header_to_index.find("Channel") != _run_column_header_to_index.end();
        has_sample = _run_column_header_to_index.find("Sample") !=  _run_column_header_to_index.end();
        n_col = _run_column_header_to_index.size();
      }
      // End of run section lines, empty line separates run and sample table
      else if (state == 1 && line.empty())
      {
        // Next line is header of Sample table
        state = 2;
      }
      // Line is run line of run section
      else if (state == 1)
      {
        parseErrorIf(n_col != cells.size(), tsv_file,
            "Wrong number of records in line");

        RunRow row;

        // Assign run and fraction
        row.run = cells[_run_column_header_to_index["Run"]].toInt();
        row.fraction = cells[_run_column_header_to_index["Fraction"]].toInt();

        // Assign channel
        row.channel = has_channel ? cells[_run_column_header_to_index["Channel"]].toInt() : 1;

        // Assign sample number
        if (has_sample)
        {
          row.sample = cells[_run_column_header_to_index["Sample"]].toInt();
        }
        else if (has_channel)
        {
          row.sample = row.channel;
        }
        else
        {
          row.sample = row.run;
        }

        // Spectra files
        row.path = findSpectraFile(
            cells[_run_column_header_to_index["Path(Spectra File)"]],
            tsv_file);
        design.run_section_.push_back(row);
      }
      // Parse header of the Condition Table
      else if (state == 2)
      {
        state = 3;
        line_number = 0;
        parse_header(
          cells,
          tsv_file,
          design.sample_section_.columnname_to_columnindex_,
          {"Sample"}, {}, true
        );
        n_col = design.sample_section_.columnname_to_columnindex_.size();
      }
      // Parse Sample Row
      else if (state == 3)
      {
        // Parse Error if sample appears multiple times
        unsigned sample = cells[design.sample_section_.columnname_to_columnindex_["Sample"]].toInt();
        parseErrorIf(
          design.sample_section_.sample_to_rowindex_.find(sample) != design.sample_section_.sample_to_rowindex_.end(),
          tsv_file,
          "Sample: " + String(sample) + " appears multiple times in the sample table"
        );
        design.sample_section_.sample_to_rowindex_[sample] = line_number++;
        design.sample_section_.content_.push_back(cells);
      }
    }

    design.sort_();
    design.checkValidRunSection_();
    return design;
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
