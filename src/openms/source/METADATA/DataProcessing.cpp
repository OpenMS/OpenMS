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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/DataProcessing.h>

using namespace std;

namespace OpenMS
{
  const std::string DataProcessing::NamesOfProcessingAction[] =
  {
    "Data processing action",
    "Charge deconvolution",
    "Deisotoping",
    "Smoothing",
    "Charge calculation",
    "Precursor recalculation",
    "Baseline reduction",
    "Peak picking",
    "Retention time alignment",
    "Calibration of m/z positions",
    "Intensity normalization",
    "Data filtering",
    "Quantitation",
    "Feature grouping",
    "Identification mapping",
    "File format conversion",
    "Conversion to mzData format",
    "Conversion to mzML format",
    "Conversion to mzXML format",
    "Conversion to DTA format"
  };

  DataProcessing::~DataProcessing()
  {

  }

  DataProcessing::DataProcessing(DataProcessing&& rhs) noexcept :
    MetaInfoInterface(std::move(rhs)),
    software_(std::move(rhs.software_)),
    processing_actions_(std::move(rhs.processing_actions_)),
    completion_time_(std::move(rhs.completion_time_))
  {
  }

  bool DataProcessing::operator==(const DataProcessing & rhs) const
  {
    return software_ == rhs.software_ &&
           processing_actions_ == rhs.processing_actions_ &&
           completion_time_ == rhs.completion_time_ &&
           MetaInfoInterface::operator==(rhs);
  }

  bool DataProcessing::operator!=(const DataProcessing & rhs) const
  {
    return !(operator==(rhs));
  }

  const Software & DataProcessing::getSoftware() const
  {
    return software_;
  }

  Software & DataProcessing::getSoftware()
  {
    return software_;
  }

  void DataProcessing::setSoftware(const Software & software)
  {
    software_ = software;
  }

  const DateTime & DataProcessing::getCompletionTime() const
  {
    return completion_time_;
  }

  void DataProcessing::setCompletionTime(const DateTime & completion_time)
  {
    completion_time_ = completion_time;
  }

  const set<DataProcessing::ProcessingAction> & DataProcessing::getProcessingActions() const
  {
    return processing_actions_;
  }

  set<DataProcessing::ProcessingAction> & DataProcessing::getProcessingActions()
  {
    return processing_actions_;
  }

  void DataProcessing::setProcessingActions(const set<DataProcessing::ProcessingAction> & processing_actions)
  {
    processing_actions_ = processing_actions;
  }

}

