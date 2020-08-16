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

#include <OpenMS/METADATA/Gradient.h>

#include <algorithm>

using namespace std;

namespace OpenMS
{

  Gradient::~Gradient()
  {

  }

  bool Gradient::operator==(const Gradient & rhs) const
  {
    return (eluents_ == rhs.eluents_) &&
           (times_ == rhs.times_) &&
           (percentages_ == rhs.percentages_);
  }

  bool Gradient::operator!=(const Gradient & rhs) const
  {
    return !(operator==(rhs));
  }

  void Gradient::addEluent(const String & eluent)
  {
    //check if the eluent is valid
    std::vector<String>::iterator elu_it = find(eluents_.begin(), eluents_.end(), eluent);
    if (elu_it != eluents_.end())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "A eluent with this name already exists!", eluent);
    }

    eluents_.push_back(eluent);
    // add zero values to percentages
    percentages_.push_back(vector<UInt>(times_.size(), 0));
  }

  void Gradient::clearEluents()
  {
    eluents_.clear();
  }

  const std::vector<String> & Gradient::getEluents() const
  {
    return eluents_;
  }

  void Gradient::addTimepoint(Int timepoint)
  {
    if ((times_.size() > 0) && (timepoint <= times_[times_.size() - 1]))
    {
      throw Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    times_.push_back(timepoint);

    // add zero values to percentages
    for (Size i = 0; i < eluents_.size(); ++i)
    {
      percentages_[i].push_back(0);
    }
  }

  void Gradient::clearTimepoints()
  {
    times_.clear();
  }

  const std::vector<Int> & Gradient::getTimepoints() const
  {
    return times_;
  }

  void Gradient::setPercentage(const String & eluent, Int timepoint, UInt percentage)
  {
    // (1) validity check

    //check if the eluent is valid
    std::vector<String>::iterator elu_it = find(eluents_.begin(), eluents_.end(), eluent);
    if (elu_it == eluents_.end())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The given eluent does not exist in the list of eluents!", eluent);
    }

    //check if the timepoint is valid
    std::vector<Int>::iterator time_it = find(times_.begin(), times_.end(), timepoint);
    if (time_it == times_.end())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The given timepoint does not exist in the list of timepoints!", String(timepoint));
    }

    // percentage is valid?
    if (percentage > 100)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The percentage should be between 0 and 100!", String(percentage));
    }

    // (2) Look up indices

    UInt elu_index(0), time_index(0);
    //look up eluents index
    for (std::vector<String>::iterator it = eluents_.begin(); it != eluents_.end(); ++it)
    {
      if (*it == eluent)
        break;
      ++elu_index;
    }
    //look up timepoint index
    for (std::vector<Int>::iterator it = times_.begin(); it != times_.end(); ++it)
    {
      if (*it == timepoint)
        break;
      ++time_index;
    }

    // (3) set percentage
    percentages_[elu_index][time_index] = percentage;
  }

  const std::vector<std::vector<UInt> > & Gradient::getPercentages() const
  {
    return percentages_;
  }

  UInt Gradient::getPercentage(const String & eluent, Int timepoint) const
  {
    // (1) validity check

    //check if the eluent is valid
    std::vector<String>::const_iterator elu_it = find(eluents_.begin(), eluents_.end(), eluent);
    if (elu_it == eluents_.end())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The given eluent does not exist in the list of eluents!", eluent);
    }

    //check if the timepoint is valid
    std::vector<Int>::const_iterator time_it = find(times_.begin(), times_.end(), timepoint);
    if (time_it == times_.end())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The given timepoint does not exist in the list of timepoints!", String(timepoint));
    }

    // (2) Look up indices

    UInt elu_index(0), time_index(0);
    //look up eluents index
    for (std::vector<String>::const_iterator it = eluents_.begin(); it != eluents_.end(); ++it)
    {
      if (*it == eluent)
        break;
      ++elu_index;
    }
    //look up timepoint index
    for (std::vector<Int>::const_iterator it = times_.begin(); it != times_.end(); ++it)
    {
      if (*it == timepoint)
        break;
      ++time_index;
    }

    // (3) return percentage
    return percentages_[elu_index][time_index];
  }

  void Gradient::clearPercentages()
  {
    percentages_.clear();
    // fill all percentages with 0
    percentages_.insert(percentages_.begin(), eluents_.size(), vector<UInt>(times_.size(), 0));
  }

  bool Gradient::isValid() const
  {
    for (Size j = 0; j < times_.size(); ++j)
    {
      Int sum = 0;
      for (Size i = 0; i < eluents_.size(); ++i)
      {
        sum += percentages_[i][j];
      }
      if (sum != 100)
      {
        return false;
      }
    }
    return true;
  }

} // namespace OpenMS
