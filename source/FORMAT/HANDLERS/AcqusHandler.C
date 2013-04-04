// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Guillaume Belz$
// $Authors: Guillaume Belz$
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/AcqusHandler.h>

#include <fstream>
#include <cmath>

using namespace std;

namespace OpenMS
{
  namespace Internal
  {

    AcqusHandler::AcqusHandler(const String & filename)
    {
      params_.clear();

      std::ifstream is(filename.c_str());
      if (!is)
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
      }

      //temporary variables
      String line;
      std::vector<String> strings(2);

      //read lines
      while (getline(is, line, '\n'))
      {
        if (line.size() < 5)
          continue;                    // minimal string = "##x=x"
        if (line.prefix(2) != String("##"))
          continue;

        if (line.split('=', strings))
        {
          if (strings.size() == 2)
          {
            params_[strings[0].substr(2)] = strings[1].trim();
          }
        }
      }

      // TOF calibration params
      dw_ = params_[String("$DW")].toDouble();
      delay_ = (Size)params_[String("$DELAY")].toInt();
      ml1_ = params_[String("$ML1")].toDouble();
      ml2_ = params_[String("$ML2")].toDouble();
      ml3_ = params_[String("$ML3")].toDouble();
      td_ = (Size) params_[String("$TD")].toInt();

      is.close();
    }

    AcqusHandler::~AcqusHandler()
    {
      params_.clear();
    }

    Size AcqusHandler::getSize()
    {
      return td_;
    }

    DoubleReal AcqusHandler::getPosition(const Size index)
    {
      DoubleReal sqrt_mz_;
      DoubleReal tof_ = dw_ * index + delay_;
      DoubleReal a_ = ml3_;
      DoubleReal b_ = sqrt(1000000000000.0 / ml1_);
      DoubleReal c_ = ml2_ - tof_;

      if (ml3_ == 0.0)
      {
        sqrt_mz_ = c_ / b_;
      }
      else
      {
        sqrt_mz_ = (sqrt(b_ * b_ - 4 * a_ * c_) - b_) / (2 * a_);
      }
      return sqrt_mz_ * sqrt_mz_;
    }

    String AcqusHandler::getParam(const String & param)
    {
      if (param == String("mzMax"))
      {
        return String(getPosition(td_ - 1));
      }
      else if (param == String("mzMin"))
      {
        return String(getPosition(0));
      }
      return params_[param];
    }

  }   // namespace Internal
} // namespace OpenMS
