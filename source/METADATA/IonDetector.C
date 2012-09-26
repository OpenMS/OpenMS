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

#include <OpenMS/METADATA/IonDetector.h>

using namespace std;

namespace OpenMS
{

  const std::string IonDetector::NamesOfType[] = {"Unknown", "Electron multiplier", "Photo multiplier", "Focal plane array", "Faraday cup", "Conversion dynode electron multiplier", "Conversion dynode photo multiplier", "Multi-collector", "Channel electron multiplier", "channeltron", "daly detector", "microchannel plate detector", "array detector", "conversion dynode", "dynode", "focal plane collector", "ion-to-photon detector", "point collector", "postacceleration detector", "photodiode array detector", "inductive detector", "electron multiplier tube"};

  const std::string IonDetector::NamesOfAcquisitionMode[] = {"Unknown", "Pulse counting", "Analog-digital converter", "Time-digital converter", "Transient recorder"};

  IonDetector::IonDetector() :
    MetaInfoInterface(),
    type_(TYPENULL),
    acquisition_mode_(ACQMODENULL),
    resolution_(0.0),
    ADC_sampling_frequency_(0.0),
    order_(0)
  {

  }

  IonDetector::IonDetector(const IonDetector & source) :
    MetaInfoInterface(source),
    type_(source.type_),
    acquisition_mode_(source.acquisition_mode_),
    resolution_(source.resolution_),
    ADC_sampling_frequency_(source.ADC_sampling_frequency_),
    order_(source.order_)
  {

  }

  IonDetector::~IonDetector()
  {

  }

  IonDetector & IonDetector::operator=(const IonDetector & source)
  {
    if (&source == this)
      return *this;

    order_ = source.order_;
    type_ = source.type_;
    acquisition_mode_ = source.acquisition_mode_;
    resolution_ = source.resolution_;
    ADC_sampling_frequency_ = source.ADC_sampling_frequency_;
    MetaInfoInterface::operator=(source);

    return *this;
  }

  bool IonDetector::operator==(const IonDetector & rhs) const
  {
    return order_ == rhs.order_ &&
           type_ == rhs.type_ &&
           acquisition_mode_ == rhs.acquisition_mode_ &&
           resolution_ == rhs.resolution_ &&
           ADC_sampling_frequency_ == rhs.ADC_sampling_frequency_ &&
           MetaInfoInterface::operator==(rhs);
  }

  bool IonDetector::operator!=(const IonDetector & rhs) const
  {
    return !(operator==(rhs));
  }

  IonDetector::Type IonDetector::getType() const
  {
    return type_;
  }

  void IonDetector::setType(IonDetector::Type type)
  {
    type_ = type;
  }

  IonDetector::AcquisitionMode IonDetector::getAcquisitionMode() const
  {
    return acquisition_mode_;
  }

  void IonDetector::setAcquisitionMode(IonDetector::AcquisitionMode acquisition_mode)
  {
    acquisition_mode_ = acquisition_mode;
  }

  DoubleReal IonDetector::getResolution() const
  {
    return resolution_;
  }

  void IonDetector::setResolution(DoubleReal resolution)
  {
    resolution_ = resolution;
  }

  DoubleReal IonDetector::getADCSamplingFrequency() const
  {
    return ADC_sampling_frequency_;
  }

  void IonDetector::setADCSamplingFrequency(DoubleReal ADC_sampling_frequency)
  {
    ADC_sampling_frequency_ = ADC_sampling_frequency;
  }

  Int IonDetector::getOrder() const
  {
    return order_;
  }

  void IonDetector::setOrder(Int order)
  {
    order_ = order;
  }

}
