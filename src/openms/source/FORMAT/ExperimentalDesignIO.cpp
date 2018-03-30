// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Authors: Timo Sachsenberg, Lukas Zimmermann $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/FORMAT/ExperimentalDesign.h>
#include <OpenMS/FORMAT/ExperimentalDesignIO.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <QtCore/QString>
#include <QtCore/QFileInfo>

#include <iostream>

using namespace std;

namespace OpenMS
{
    String findSpectraFile(const String &spec_file, const String &tsv_file, const bool require_spectra_files)
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

      // Fail if the existence of the spectra files is required but they do not exist
      if (require_spectra_files && !File::exists(result))
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, tsv_file,
                                    "Error: Spectra file does not exist: '" + result + "'");
      }
      return result;
    }

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

    void ExperimentalDesignIO::parseHeader_(
      const StringList &header,
      const String &filename,
      std::map <String, Size> &column_map,
      const std::set <String> &required,
      const std::set <String> &optional,
      const bool allow_other_header)
    {
      // Headers as set
      std::set <String> header_set(header.begin(), header.end());
      parseErrorIf(header_set.size() != header.size(), filename,
                   "Some column headers of the table appear multiple times!");

      // Check that all required headers are there
      for (const String &req_header : required)
      {
        parseErrorIf( ! ListUtils::contains(header, req_header), filename,
                      "Missing column header: " + req_header);
      }
      // Assign index in column map and check for weird headers
      for (Size i = 0; i < header.size(); ++i)
      {
        const String &h = header[i];

        // A header is unexpected if it is neither required nor optional and we do not allow other headers
        const bool header_unexpected = (required.find(h) == required.end()) && (optional.find(h) == optional.end());
        parseErrorIf(
          !allow_other_header && header_unexpected,
          filename,
          "Header not allowed in this section of the Experimental Design: " + h
        );
        column_map[h] = i;
      }
    }

    // static
    ExperimentalDesign ExperimentalDesignIO::load(const String &tsv_file, const bool require_spectra_file)
    {
      ExperimentalDesign design;
      design.run_section_.clear();

      bool has_sample(false);
      bool has_channel(false);

      // Maps the column header string to the column index for
      // the run section
      std::map <String, Size> run_column_header_to_index;

      unsigned line_number(0);

      enum ParseState { RUN_HEADER, RUN_CONTENT, SAMPLE_HEADER, SAMPLE_CONTENT };

      ParseState state(RUN_HEADER);
      Size n_col = 0;

      const TextFile text_file(tsv_file, true);
      for (String s : text_file)
      {
        // skip empty lines (except in state RUN_CONTENT, where the sample table is read)
        const String line(s.trim());

        if (line.empty() && state != RUN_CONTENT)
        {
          continue;
        }

        // Now split the line into individual cells
        StringList cells;
        line.split("\t", cells);

        // Trim whitespace from all cells (so , foo , and  ,foo, is the same)
        std::transform(cells.begin(), cells.end(), cells.begin(),
                       [](String cell) -> String { return cell.trim(); });

        if (state == RUN_HEADER)
        {
          state = RUN_CONTENT;
          parseHeader_(
            cells,
            tsv_file,
            run_column_header_to_index,
            {"Run", "Fraction", "Path(Spectra File)"},
            {"Channel", "Sample"}, false
          );
          has_channel = run_column_header_to_index.find("Channel") != run_column_header_to_index.end();
          has_sample = run_column_header_to_index.find("Sample") != run_column_header_to_index.end();
          n_col = run_column_header_to_index.size();
        }
          // End of run section lines, empty line separates run and sample table
        else if (state == RUN_CONTENT && line.empty())
        {
          // Next line is header of Sample table
          state = SAMPLE_HEADER;
        }
          // Line is run line of run section
        else if (state == RUN_CONTENT)
        {
          parseErrorIf(n_col != cells.size(), tsv_file,
                       "Wrong number of records in line");

          ExperimentalDesign::RunRow row;

          // Assign run and fraction
          row.run = cells[run_column_header_to_index["Run"]].toInt();
          row.fraction = cells[run_column_header_to_index["Fraction"]].toInt();

          // Assign channel
          row.channel = has_channel ? cells[run_column_header_to_index["Channel"]].toInt() : 1;

          // Assign sample number
          if (has_sample)
          {
            row.sample = cells[run_column_header_to_index["Sample"]].toInt();
          }
          else
          {
            row.sample = has_channel ? row.channel : row.run;
          }

          // Spectra files
          row.path = findSpectraFile(
            cells[run_column_header_to_index["Path(Spectra File)"]],
            tsv_file,
            require_spectra_file);
          design.run_section_.push_back(row);
        }
          // Parse header of the Condition Table
        else if (state == SAMPLE_HEADER)
        {
          state = SAMPLE_CONTENT;
          line_number = 0;
          parseHeader_(
            cells,
            tsv_file,
            design.sample_section_.columnname_to_columnindex_,
            {"Sample"}, {}, true
          );
          n_col = design.sample_section_.columnname_to_columnindex_.size();
        }
          // Parse Sample Row
        else if (state == SAMPLE_CONTENT)
        {
          // Parse Error if sample appears multiple times
          unsigned sample = cells[design.sample_section_.columnname_to_columnindex_["Sample"]].toInt();
          parseErrorIf(
            design.sample_section_.sample_to_rowindex_.find(sample) !=
            design.sample_section_.sample_to_rowindex_.end(),
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

    ExperimentalDesign ExperimentalDesignIO::fromConsensusMap(const ConsensusMap &cm)
    {
      ExperimentalDesign experimental_design;
      // path of the original MS run (mzML / raw file)
      StringList ms_run_paths;
      cm.getPrimaryMSRunPath(ms_run_paths);

      // no fractionation -> as many runs as samples
      // each consensus element corresponds to one sample abundance
      size_t sample(1);
      ExperimentalDesign::RunRows rows;
      for (const auto &f : ms_run_paths)
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
      experimental_design.setRunSection(rows);
      LOG_INFO << "Experimental design (ConsensusMap derived):\n"
               << "  files: " << experimental_design.getNumberOfMSFiles()
               << "  fractions: " << experimental_design.getNumberOfFractions()
               << "  channels: " << experimental_design.getNumberOfChannels()
               << "  samples: " << experimental_design.getNumberOfSamples() << "\n"
               << endl;
      return experimental_design;
    }

    ExperimentalDesign ExperimentalDesignIO::fromFeatureMap(const FeatureMap &fm)
    {
      ExperimentalDesign experimental_design;
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
      experimental_design.setRunSection(rows);
      LOG_INFO << "Experimental design (FeatureMap derived):\n"
               << "  files: " << experimental_design.getNumberOfMSFiles()
               << "  fractions: " << experimental_design.getNumberOfFractions()
               << "  channels: " << experimental_design.getNumberOfChannels()
               << "  samples: " << experimental_design.getNumberOfSamples() << "\n"
               << endl;
      return experimental_design;
    }

    ExperimentalDesign ExperimentalDesignIO::fromIdentifications(const vector <ProteinIdentification> &proteins)
    {
      ExperimentalDesign experimental_design;
      // path of the original MS files (mzML / raw file)
      StringList ms_run_paths;
      for (const auto &protein : proteins)
      {
        StringList tmp_ms_run_paths;
        protein.getPrimaryMSRunPath(tmp_ms_run_paths);
        if (tmp_ms_run_paths.size() != 1)
        {
          throw Exception::MissingInformation(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION,
            "ProteinIdentification annotated with " + String(tmp_ms_run_paths.size()) +
            " MS files. Must be exactly one.");
        }
        ms_run_paths.push_back(tmp_ms_run_paths[0]);
      }

      // no fractionation -> as many runs as samples
      // each identification run corresponds to one sample abundance
      unsigned sample(1);
      ExperimentalDesign::RunRows rows;
      for (const auto &f : ms_run_paths)
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
      experimental_design.setRunSection(rows);
      LOG_INFO << "Experimental design (Identification derived):\n"
               << "  files: " << experimental_design.getNumberOfMSFiles()
               << "  fractions: " << experimental_design.getNumberOfFractions()
               << "  channels: " << experimental_design.getNumberOfChannels()
               << "  samples: " << experimental_design.getNumberOfSamples() << "\n"
               << endl;
      return experimental_design;
    }
}
