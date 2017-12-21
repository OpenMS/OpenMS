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
// $Authors: Witold Wolski $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_OPENSWATHALGO_DATAACCESS_DATAFRAMEWRITER_H
#define OPENMS_ANALYSIS_OPENSWATH_OPENSWATHALGO_DATAACCESS_DATAFRAMEWRITER_H

#include <fstream>
#include <string>
#include <vector>

#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/OpenSwathAlgoConfig.h>

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

#endif // OPENMS_ANALYSIS_OPENSWATH_OPENSWATHALGO_DATAACCESS_DATAFRAMEWRITER_H
