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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Digestion.h>

using namespace std;

namespace OpenMS
{


  Digestion::Digestion() :
    SampleTreatment("Digestion"),
    enzyme_(""),
    digestion_time_(0.0),
    temperature_(0.0),
    ph_(0.0)
  {

  }

  Digestion::Digestion(const Digestion & source) :
    SampleTreatment(source),
    enzyme_(source.enzyme_),
    digestion_time_(source.digestion_time_),
    temperature_(source.temperature_),
    ph_(source.ph_)
  {

  }

  Digestion::~Digestion()
  {

  }

  SampleTreatment * Digestion::clone() const
  {
    SampleTreatment * tmp = new Digestion(*this);
    return tmp;
  }

  Digestion & Digestion::operator=(const Digestion & source)
  {
    if (&source == this)
      return *this;

    SampleTreatment::operator=(source);
    enzyme_ = source.enzyme_;
    digestion_time_ = source.digestion_time_;
    temperature_ = source.temperature_;
    ph_ = source.ph_;

    return *this;
  }

  bool Digestion::operator==(const SampleTreatment & rhs) const
  {
    if (type_ != rhs.getType())
      return false;

    const Digestion * tmp = dynamic_cast<const Digestion *>(&rhs);
    return SampleTreatment::operator==(* tmp) &&
           enzyme_ == tmp->enzyme_ &&
           digestion_time_ == tmp->digestion_time_ &&
           temperature_ == tmp->temperature_ &&
           ph_ == tmp->ph_
    ;
  }

  const String & Digestion::getEnzyme() const
  {
    return enzyme_;
  }

  void Digestion::setEnzyme(const String & enzyme)
  {
    enzyme_ = enzyme;
  }

  DoubleReal Digestion::getDigestionTime() const
  {
    return digestion_time_;
  }

  void Digestion::setDigestionTime(DoubleReal digestion_time)
  {
    digestion_time_ = digestion_time;
  }

  DoubleReal Digestion::getTemperature() const
  {
    return temperature_;
  }

  void Digestion::setTemperature(DoubleReal temperature)
  {
    temperature_ = temperature;
  }

  DoubleReal Digestion::getPh() const
  {
    return ph_;
  }

  void Digestion::setPh(DoubleReal ph)
  {
    ph_ = ph;
  }

}
