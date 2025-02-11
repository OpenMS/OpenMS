// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Witold Wolski $
// --------------------------------------------------------------------------

#include <OpenMS/OPENSWATHALGO/DATAACCESS/DataFrameWriter.h>

#include <iostream>
#include <iomanip>

namespace OpenSwath
{

  IDataFrameWriter::~IDataFrameWriter()
  {
  }

  DataMatrix::DataMatrix() :
    colnames_(), rownames_(), store_()
  {
  }

  void DataMatrix::store(const std::string& rowname,
                         const std::vector<double>& values)
  {
    rownames_.push_back(rowname);
    store_.push_back(values);
  }

  void DataMatrix::colnames(const std::vector<std::string>& colnames)
  {
    colnames_ = colnames;
  }

  CSVWriter::CSVWriter(std::string filename) :
    sep_("\t"), eol_("\n")
  {
    file_stream_.open(filename.c_str());
  }

  void CSVWriter::store(const std::string& rowname,
                        const std::vector<double>& values)
  {
    file_stream_ << rowname;
    file_stream_ << sep_;
    std::size_t ncol = values.size();
    for (size_t i = 0; i < ncol; ++i)
    {
      file_stream_ << std::setprecision(5) << values[i];
      if (i < (ncol - 1))
      {
        file_stream_ << sep_;
      }
    }
    file_stream_ << eol_; //append line-end
  }

  CSVWriter::~CSVWriter()
  {
    file_stream_.flush();
    file_stream_.close();
    std::cout << "have flushed and closed the file stream" << std::endl;
  }

  void CSVWriter::colnames(const std::vector<std::string>& colnames)
  {
    std::size_t ncol = colnames.size();
    for (size_t i = 0; i < ncol; ++i)
    {
      file_stream_ << colnames[i];
      if (i < (ncol - 1))
      {
        file_stream_ << sep_;
      }
    }
    file_stream_ << eol_; //append line-end
  }

} // namespace OpenSwath
