// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/FORMAT/ExperimentalDesignFile.h>
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

    void ExperimentalDesignFile::parseHeader_(
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

    bool ExperimentalDesignFile::isOneTableFile_(const TextFile& text_file)
    {
      // To dermine if we have a *separate* sample table, we check if a row is present
      // that contains "Sample" but no "Fraction_Group".
      for (String s : text_file)
      {
        const String line(s.trim());
        if (line.empty()) { continue; }
        StringList cells;
        line.split("\t", cells);
        // check if we are outside of run section (no Fraction_Group column) 
        // but in sample section (Sample column present)
        if (std::count(cells.begin(), cells.end(), "Fraction_Group") == 0 
         && std::count(cells.begin(), cells.end(), "Sample") == 1)
        {
          return false;
        }
      }
      return true;
    }

    ExperimentalDesign ExperimentalDesignFile::parseOneTableFile_(const TextFile& text_file, const String& tsv_file, bool require_spectra_file)
    {
      ExperimentalDesign::MSFileSection msfile_section;
      bool has_sample(false);
      bool has_label(false);

      // Attributes of the sample section
      std::vector< std::vector < String > > sample_content_;
      std::map< unsigned, Size > sample_sample_to_rowindex_;
      std::map< String, Size > sample_columnname_to_columnindex_;

      // Maps the column header string to the column index for
      // the file section
      std::map <String, Size> fs_column_header_to_index;

      unsigned line_number(0);

      enum ParseState { RUN_HEADER, RUN_CONTENT };

      ParseState state(RUN_HEADER);
      Size n_col = 0;

      for (String s : text_file)
      {
        const String line(s.trim());
	
        if (line.empty()) { continue; }

        // Now split the line into individual cells
        StringList cells;
        line.split("\t", cells);

        // Trim whitespace from all cells (so , foo , and  ,foo, is the same)
        std::transform(cells.begin(), cells.end(), cells.begin(),
                       [](String cell) -> String { return cell.trim(); });

        if (state == RUN_HEADER)
        {
          state = RUN_CONTENT;
          line_number = 0;
          parseHeader_(
            cells,
            tsv_file,
            fs_column_header_to_index,
            {"Fraction_Group", "Fraction", "Spectra_Filepath"},
            {"Label", "Sample"}, true
          );
          has_label = fs_column_header_to_index.find("Label") != fs_column_header_to_index.end();
          has_sample = fs_column_header_to_index.find("Sample") != fs_column_header_to_index.end();
     
          // readd label column to end of header
          if (!has_label)
          {
            size_t hs = fs_column_header_to_index.size();
            fs_column_header_to_index["Label"] = hs;
            cells.push_back("Label");
          }

          // readd sample column to end of header
          if (!has_sample)
          {
            size_t hs = fs_column_header_to_index.size();
            fs_column_header_to_index["Sample"] = hs;
            cells.push_back("Sample");
          }
    
          n_col = fs_column_header_to_index.size();

          // determine columns with sample metainfo like condition or replication
          for (size_t i = 0; i != cells.size(); ++i)
          {
            const String& c = cells[i];
            if (c != "Fraction_Group" && c != "Fraction" 
             && c != "Spectra_Filepath" && c != "Label")
            { 
              sample_columnname_to_columnindex_[c] = i;
            }
          }
        }
        else if (state == RUN_CONTENT)
        {
          // readd label column as we already did in the header
          if (!has_label) { cells.push_back("1"); }

          // Assign label, fall back to 1 if column is missing
          int label = cells[fs_column_header_to_index["Label"]].toInt();
          int fraction = cells[fs_column_header_to_index["Fraction"]].toInt();
          int fraction_group = cells[fs_column_header_to_index["Fraction_Group"]].toInt();
                    // readd sample column
          if (!has_sample) 
          {
            int sample = fraction_group; // deducing the sample in the case of multiplexed could be done if label > 1 information is there (e.g., max(label) * (fraction_group - 1) + label 
            cells.push_back(String(sample)); 
          }

          int sample = cells[fs_column_header_to_index["Sample"]].toInt();
          parseErrorIf(sample < 1, tsv_file,
                       "Sample index may not be smaller than 1");


          parseErrorIf(n_col != cells.size(), tsv_file,
                       "Wrong number of records in line");

          ExperimentalDesign::MSFileSectionEntry e;

          // Assign fraction group and fraction
          e.fraction_group = fraction_group;
          e.fraction = fraction;
          e.label = label;
          e.sample = sample;

          sample_sample_to_rowindex_[e.sample] = line_number++;

          // get indices of sample metadata and store content in sample cells
          StringList sample_cells;
          for (const auto & c2i : sample_columnname_to_columnindex_)
          {
            sample_cells.push_back(cells[c2i.second]);
          }
          sample_content_.push_back(sample_cells);

          // Spectra files
          e.path = findSpectraFile(
            cells[fs_column_header_to_index["Spectra_Filepath"]],
            tsv_file,
            require_spectra_file);
       
          msfile_section.push_back(e);
        }
      }

      // Assign correct position in sample column (without "Fraction_Group", "Fraction", "Spectra_Filepath", "Label")
      int sample_index = 0;
      for (auto & c : sample_columnname_to_columnindex_)
      {
        c.second = sample_index++;
      }

      // Create Sample Section and set in design
      ExperimentalDesign::SampleSection sample_section(
        sample_content_,
        sample_sample_to_rowindex_,
        sample_columnname_to_columnindex_);

      // Create experimentalDesign
      ExperimentalDesign design(msfile_section, sample_section);

      return design;
    }

    ExperimentalDesign ExperimentalDesignFile::parseTwoTableFile_(const TextFile& text_file, const String& tsv_file, bool require_spectra_file)
    {
      ExperimentalDesign::MSFileSection msfile_section;
      bool has_sample(false);
      bool has_label(false);

      // Attributes of the sample section
      std::vector< std::vector < String > > sample_content_;
      std::map< unsigned, Size > sample_sample_to_rowindex_;
      std::map< String, Size > sample_columnname_to_columnindex_;

      // Maps the column header string to the column index for
      // the file section
      std::map <String, Size> fs_column_header_to_index;

      unsigned line_number(0);

      enum ParseState { RUN_HEADER, RUN_CONTENT, SAMPLE_HEADER, SAMPLE_CONTENT };

      ParseState state(RUN_HEADER);
      Size n_col = 0;

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
            fs_column_header_to_index,
            {"Fraction_Group", "Fraction", "Spectra_Filepath"},
            {"Label", "Sample"}, false
          );
          has_label = fs_column_header_to_index.find("Label") != fs_column_header_to_index.end();
          has_sample = fs_column_header_to_index.find("Sample") != fs_column_header_to_index.end();
          
          n_col = fs_column_header_to_index.size();
        }
          // End of file section lines, empty line separates file and sample section
        else if (state == RUN_CONTENT && line.empty())
        {
          // Next line is header of Sample table
          state = SAMPLE_HEADER;
        }
          // Line is file section line
        else if (state == RUN_CONTENT)
        {
          parseErrorIf(n_col != cells.size(), tsv_file,
                       "Wrong number of records in line");

          ExperimentalDesign::MSFileSectionEntry e;

          // Assign fraction group and fraction
          e.fraction_group = cells[fs_column_header_to_index["Fraction_Group"]].toInt();
          e.fraction = cells[fs_column_header_to_index["Fraction"]].toInt();

          // Assign label
          e.label = has_label ? cells[fs_column_header_to_index["Label"]].toInt() : 1;

          // Assign sample number
          if (has_sample)
          {
            e.sample = cells[fs_column_header_to_index["Sample"]].toInt();
          }
          else
          {
            e.sample = e.fraction_group; // deducing the sample in the case of multiplexed could be done if label > 1 information is there (e.g., max(label) * (fraction_group - 1) + label 
          }

          // Spectra files
          e.path = findSpectraFile(
            cells[fs_column_header_to_index["Spectra_Filepath"]],
            tsv_file,
            require_spectra_file);
          msfile_section.push_back(e);
        }
          // Parse header of the Condition Table
        else if (state == SAMPLE_HEADER)
        {
          state = SAMPLE_CONTENT;
          line_number = 0;
          parseHeader_(
            cells,
            tsv_file,
            sample_columnname_to_columnindex_,
            {"Sample"}, {}, true
          );
          n_col = sample_columnname_to_columnindex_.size();
        }
          // Parse Sample Row
        else if (state == SAMPLE_CONTENT)
        {
          // Parse Error if sample appears multiple times
          unsigned sample = cells[sample_columnname_to_columnindex_["Sample"]].toInt();
          parseErrorIf(
            sample_sample_to_rowindex_.find(sample) != sample_sample_to_rowindex_.end(),
            tsv_file,
            "Sample: " + String(sample) + " appears multiple times in the sample table"
          );
          sample_sample_to_rowindex_[sample] = line_number++;
          sample_content_.push_back(cells);
        }
      }

      // Create Sample Section and set in design
      ExperimentalDesign::SampleSection sample_section(
        sample_content_,
        sample_sample_to_rowindex_,
        sample_columnname_to_columnindex_);

      // Create experimentalDesign
      ExperimentalDesign design(msfile_section, sample_section);

      return design;
    }

    // static
    ExperimentalDesign ExperimentalDesignFile::load(const String &tsv_file, const bool require_spectra_file)
    {
      const TextFile text_file(tsv_file, true);

      // check if we have information stored in one or two files
      bool has_one_table = isOneTableFile_(text_file);

      if (has_one_table)
      {
        return parseOneTableFile_(text_file, tsv_file, require_spectra_file);
      }
      else // two tables
      {
        return parseTwoTableFile_(text_file, tsv_file, require_spectra_file);
      }
    }
}

