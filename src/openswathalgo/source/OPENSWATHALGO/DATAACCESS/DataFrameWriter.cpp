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
        file_stream_ << sep_;
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
        file_stream_ << sep_;
    }
    file_stream_ << eol_; //append line-end
  }

} // namespace OpenSwath
