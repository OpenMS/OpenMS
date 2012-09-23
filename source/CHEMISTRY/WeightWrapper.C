// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/WeightWrapper.h>

namespace OpenMS
{

  WeightWrapper::WeightWrapper() :
    weight_mode_(MONO)
  {
  }

  WeightWrapper::WeightWrapper(const WEIGHTMODE weight_mode) :
    weight_mode_(weight_mode)
  {
  }

  WeightWrapper::WeightWrapper(const WeightWrapper & source) :
    weight_mode_(source.weight_mode_)
  {
  }

  WeightWrapper::~WeightWrapper()
  {
  }

  void WeightWrapper::setWeightMode(const WEIGHTMODE mode)
  {
    if (mode >= WeightWrapper::SIZE_OF_WEIGHTMODE)
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "setWeightMode() received illegal 'mode' value!");
    weight_mode_ = mode;
  }

  WeightWrapper::WEIGHTMODE WeightWrapper::getWeightMode() const
  {
    return weight_mode_;
  }

  DoubleReal WeightWrapper::getWeight(const AASequence & aa) const
  {
    if (weight_mode_ == WeightWrapper::MONO)
      return aa.getMonoWeight();
    else
      return aa.getAverageWeight();
  }

  DoubleReal WeightWrapper::getWeight(const EmpiricalFormula & ef) const
  {
    if (weight_mode_ == WeightWrapper::MONO)
      return ef.getMonoWeight();
    else
      return ef.getAverageWeight();
  }

  DoubleReal WeightWrapper::getWeight(const Residue & r, Residue::ResidueType res_type) const
  {
    if (weight_mode_ == WeightWrapper::MONO)
      return r.getMonoWeight(res_type);
    else
      return r.getAverageWeight(res_type);
  }

}
