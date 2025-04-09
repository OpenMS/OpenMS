// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Lukas Zimmermann $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/ExperimentalDesignFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/File.h>
#include <QtCore/QFileInfo>
#include <QtCore/QString>
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
        // file name is relative, so we need to figure out the correct folder

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
    void ExperimentalDesignFile::parseErrorIf_(const bool test, const String &filename, const String &message)
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
      parseErrorIf_(header_set.size() != header.size(), filename, "Some column headers of the table appear multiple times!");

      // Check that all required headers are there
      for (const String &req_header : required)
      {
        parseErrorIf_(!ListUtils::contains(header, req_header), filename, "Missing column header: " + req_header);
      }
      // Assign index in column map and check for weird headers
      for (Size i = 0; i < header.size(); ++i)
      {
        const String &h = header[i];

        // A header is unexpected if it is neither required nor optional and we do not allow other headers
        const bool header_unexpected = (required.find(h) == required.end()) && (optional.find(h) == optional.end());
        parseErrorIf_(!allow_other_header && header_unexpected, filename, "Header not allowed in this section of the Experimental Design: " + h);
        column_map[h] = i;
      }
    }

    bool ExperimentalDesignFile::isOneTableFile_(const TextFile& text_file)
    {
      // To determine if we have a *separate* sample table, we check if a row is present
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

      /// Content of the sample section (vector of rows with column contents as strings)
      std::vector< std::vector < String > > sample_content_;
      /// Sample name to sample (row) index
      std::map< unsigned, Size > sample_sample_to_rowindex_;
      /// Inside each row, the value for a column can be accessed by name with this map
      std::map< String, Size > sample_columnname_to_columnindex_;

      /// Maps the column header string to the column index for
      /// the file section
      std::map <String, Size> fs_column_header_to_index;

      /// Maps the sample name to all values of sample-related columns
      std::map <String, std::vector<String>> sample_content_map;

      /// Maps the sample name to its index (i.e., the order how it gets added to the sample section)
      std::map<String, Size> samplename_to_index;
      enum ParseState { RUN_HEADER, RUN_CONTENT };

      ParseState state(RUN_HEADER);
      Size n_col = 0;

      for (String s : text_file)
      {
        const String line(s.trim());
	
        if (line.hasPrefix("#") || line.empty()) { continue; }

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
            {"Label", "Sample"}, true
          );
          has_label = fs_column_header_to_index.find("Label") != fs_column_header_to_index.end();
          has_sample = fs_column_header_to_index.find("Sample") != fs_column_header_to_index.end();

          if (!has_label) // add label column to end of header
          {
            size_t hs = fs_column_header_to_index.size();
            fs_column_header_to_index["Label"] = hs;
            cells.push_back("Label");
          }

          if (!has_sample) // add sample column to end of header
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
          // if no label column exists -> label free
          // -> add label column with label 1 at the end of every row
          if (!has_label) { cells.push_back("1"); }

          // Assign label
          int label = cells[fs_column_header_to_index["Label"]].toInt();
          int fraction = cells[fs_column_header_to_index["Fraction"]].toInt();
          int fraction_group = cells[fs_column_header_to_index["Fraction_Group"]].toInt();

          // read sample column
          Size sample = 1;
          String samplename = "";
          if (!has_sample) 
          {
            samplename = fraction_group; // deducing the sample in the case of multiplexed could be done if label > 1 information is there (e.g., max(label) * (fraction_group - 1) + label
            cells.push_back(samplename);
          }

          samplename = cells[fs_column_header_to_index["Sample"]];
          parseErrorIf_(n_col != cells.size(), tsv_file, "Wrong number of records in line");

          const auto& [it, inserted] = samplename_to_index.emplace(samplename, samplename_to_index.size());
          sample = it->second;

          // get indices of sample metadata and store content in sample cells
          StringList sample_cells;
          for (const auto & c2i : sample_columnname_to_columnindex_)
          {
            sample_cells.push_back(cells[c2i.second]);
          }

          if (inserted)
          {
            sample_content_.push_back(sample_cells);
          }
          else
          {
            if (sample_content_[it->second] != sample_cells)
            {
              OPENMS_LOG_WARN << "Warning: Factors for the same sample do not match." << std::endl;
            }
          }

          ExperimentalDesign::MSFileSectionEntry e;

          // Assign fraction group and fraction
          e.fraction_group = fraction_group;
          e.fraction = fraction;
          e.label = label;
          e.sample = sample;

          // Spectra files
          e.path = findSpectraFile(
            cells[fs_column_header_to_index["Spectra_Filepath"]],
            tsv_file,
            require_spectra_file);
       
          msfile_section.push_back(e);
        }
      }

      // Assign correct position for columns in sample section subset
      // (i.e., without "Fraction_Group", "Fraction", "Spectra_Filepath", "Label")
      int sample_index = 0;
      for (auto & c : sample_columnname_to_columnindex_)
      {
        c.second = sample_index++;
      }

      // Create Sample Section and set in design
      ExperimentalDesign::SampleSection sample_section(
        sample_content_,
        samplename_to_index,
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
      std::map< String, Size > sample_sample_to_rowindex_;
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
	      // also skip comment lines
        if (line.hasPrefix("#") || (line.empty() && state != RUN_CONTENT))
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
          parseErrorIf_(n_col != cells.size(), tsv_file, "Wrong number of records in line");

          ExperimentalDesign::MSFileSectionEntry e;

          // Assign fraction group and fraction
          e.fraction_group = cells[fs_column_header_to_index["Fraction_Group"]].toInt();
          e.fraction = cells[fs_column_header_to_index["Fraction"]].toInt();

          // Assign label
          e.label = has_label ? cells[fs_column_header_to_index["Label"]].toInt() : 1;

          // Assign sample number
          if (has_sample)
          {
            //e.sample has to be filled after the sample section was read
            e.sample_name = cells[fs_column_header_to_index["Sample"]];
          }
          else
          {
            e.sample = e.fraction_group; // TODO: deducing the sample in the case of multiplexed
                                         //  could be done if label > 1 information is there
                                         //  (e.g., max(label) * (fraction_group - 1) + label
            e.sample_name = "Fraction group " + String(e.fraction_group);
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
          const String& sample = cells[sample_columnname_to_columnindex_["Sample"]];
          parseErrorIf_(sample_sample_to_rowindex_.find(sample) != sample_sample_to_rowindex_.end(),
                        tsv_file,
                        "Sample: " + String(sample) + " appears multiple times in the sample table");
          sample_sample_to_rowindex_[sample] = line_number++;
          sample_content_.push_back(cells);
        }
      }

      for (auto& e : msfile_section)
      {
        e.sample = sample_sample_to_rowindex_.at(e.sample_name);
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
    ExperimentalDesign ExperimentalDesignFile::load(const TextFile &text_file, const bool require_spectra_file, String filename = "--no design file provided--")
    {
      // check if we have information stored in one or two files
      bool has_one_table = isOneTableFile_(text_file);

      if (has_one_table)
      {
        return parseOneTableFile_(text_file, filename, require_spectra_file);
      }
      else // two tables
      {
        return parseTwoTableFile_(text_file, filename, require_spectra_file);
      }
    }

    ExperimentalDesign ExperimentalDesignFile::load(const String &tsv_file, const bool require_spectra_file)
    {
      const TextFile text_file(tsv_file, true);
      return load(text_file, require_spectra_file, tsv_file);
    }

}

