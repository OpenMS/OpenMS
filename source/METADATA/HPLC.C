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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/HPLC.h>

using namespace std;

namespace OpenMS
{
  HPLC::HPLC() :
//      PersistentObject(),
    instrument_(),
    column_(),
    temperature_(21),
    pressure_(0),
    flux_(0),
    comment_(),
    gradient_()
  {

  }

  HPLC::HPLC(const HPLC & source) :
//    PersistentObject(source),
    instrument_(source.instrument_),
    column_(source.column_),
    temperature_(source.temperature_),
    pressure_(source.pressure_),
    flux_(source.flux_),
    comment_(source.comment_),
    gradient_(source.gradient_)
  {

  }

  HPLC::~HPLC()
  {

  }

  HPLC & HPLC::operator=(const HPLC & source)
  {
    if (source == *this)
      return *this;

//    PersistentObject::operator = (source);
    instrument_ = source.instrument_;
    column_ = source.column_;
    temperature_ = source.temperature_;
    pressure_ = source.pressure_;
    flux_ = source.flux_;
    comment_ = source.comment_;
    gradient_ = source.gradient_;

    return *this;
  }

  bool HPLC::operator==(const HPLC & rhs) const
  {
    return (instrument_ == rhs.instrument_) &&
           (column_ == rhs.column_) &&
           (temperature_ == rhs.temperature_) &&
           (pressure_ == rhs.pressure_) &&
           (flux_ == rhs.flux_) &&
           (comment_ == rhs.comment_) &&
           (gradient_ == rhs.gradient_);
  }

  bool HPLC::operator!=(const HPLC & rhs) const
  {
    return !(operator==(rhs));
  }

  const String & HPLC::getInstrument() const
  {
    return instrument_;
  }

  void HPLC::setInstrument(const String & instrument)
  {
    instrument_ = instrument;
  }

  const String & HPLC::getColumn() const
  {
    return column_;
  }

  void HPLC::setColumn(const String & column)
  {
    column_ = column;
  }

  Int HPLC::getTemperature() const
  {
    return temperature_;
  }

  void HPLC::setTemperature(Int temperature)
  {
    temperature_ = temperature;
  }

  UInt HPLC::getPressure() const
  {
    return pressure_;
  }

  void HPLC::setPressure(UInt pressure)
  {
    pressure_ = pressure;
  }

  UInt HPLC::getFlux() const
  {
    return flux_;
  }

  void HPLC::setFlux(UInt flux)
  {
    flux_ = flux;
  }

  String HPLC::getComment() const
  {
    return comment_;
  }

  void HPLC::setComment(String comment)
  {
    comment_ = comment;
  }

  Gradient & HPLC::getGradient()
  {
    return gradient_;
  }

  const Gradient & HPLC::getGradient() const
  {
    return gradient_;
  }

  void HPLC::setGradient(const Gradient & gradient)
  {
    gradient_ = gradient;
  }

} // namespace OpenMS
