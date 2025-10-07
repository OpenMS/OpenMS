// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Witold Wolski $
// --------------------------------------------------------------------------

#pragma once

#include <fstream>
#include <string>
#include <vector>

#include <OpenMS/OPENSWATHALGO/OpenSwathAlgoConfig.h>

namespace OpenSwath
{
  struct OPENSWATHALGO_DLLAPI IDataFrameWriter
  {
    virtual ~IDataFrameWriter();
    virtual void colnames(const std::vector<std::string>& colnames) = 0;
    virtual void store(const std::string& rowname,
                       const std::vector<double>& values) = 0;
  };

  struct OPENSWATHALGO_DLLAPI DataMatrix :
    IDataFrameWriter
  {
private:
    std::vector<std::string> colnames_;
    std::vector<std::string> rownames_;
    std::vector<std::vector<double> > store_;

public:
    DataMatrix();

    void store(const std::string& rowname,
               const std::vector<double>& values) override;

    void colnames(const std::vector<std::string>& colnames) override;

  };

  struct OPENSWATHALGO_DLLAPI CSVWriter :
    IDataFrameWriter
  {
private:
    std::ofstream file_stream_;
    std::string sep_;
    std::string eol_;

public:
    explicit CSVWriter(std::string filename);

    void store(const std::string& rowname,
               const std::vector<double>& values) override;

    ~CSVWriter() override;

    void colnames(const std::vector<std::string>& colnames) override;

  };
}

