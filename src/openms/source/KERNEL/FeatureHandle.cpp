// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Erhan Kenar $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/FeatureHandle.h>

#include <OpenMS/KERNEL/BaseFeature.h>

namespace OpenMS
{
  FeatureHandle::FeatureHandle() :
    Peak2D(),
    UniqueIdInterface(),
    map_index_(0),
    charge_(0),
    width_(0)
  {}

  FeatureHandle::FeatureHandle(UInt64 map_index, const Peak2D & point, UInt64 element_index) :
    Peak2D(point),
    map_index_(map_index),
    charge_(0),
    width_(0)
  {
    setUniqueId(element_index);
  }

  FeatureHandle::FeatureHandle(UInt64 map_index, const BaseFeature & feature) :
    Peak2D(feature),
    UniqueIdInterface(feature),
    map_index_(map_index),
    charge_(feature.getCharge()),
    width_(feature.getWidth())
  {}

  FeatureHandle::FeatureHandle(const FeatureHandle & rhs) :
    Peak2D(rhs),
    UniqueIdInterface(rhs),
    map_index_(rhs.map_index_),
    charge_(rhs.charge_),
    width_(rhs.width_)
  {}

  FeatureHandle & FeatureHandle::operator=(const FeatureHandle & rhs)
  {
    Peak2D::operator=(rhs);
    UniqueIdInterface::operator=(rhs);
    map_index_ = rhs.map_index_;
    charge_ = rhs.charge_;
    width_ = rhs.width_;

    return *this;
  }

  FeatureHandle::~FeatureHandle()
  {}

  UInt64 FeatureHandle::getMapIndex() const
  {
    return map_index_;
  }

  void FeatureHandle::setMapIndex(UInt64 i)
  {
    map_index_ = i;
  }

  void FeatureHandle::setCharge(FeatureHandle::ChargeType charge)
  {
    charge_ = charge;
  }

  FeatureHandle::ChargeType FeatureHandle::getCharge() const
  {
    return charge_;
  }

  void FeatureHandle::setWidth(FeatureHandle::WidthType width)
  {
    width_ = width;
  }

  FeatureHandle::WidthType FeatureHandle::getWidth() const
  {
    return width_;
  }

  bool FeatureHandle::operator==(const FeatureHandle & i) const
  {
    return (Peak2D::operator==(i))
           && (UniqueIdInterface::operator==(i))
           && (map_index_ == i.map_index_)
           && (charge_ == i.charge_)
           && (width_ == i.width_);
  }

  bool FeatureHandle::operator!=(const FeatureHandle & i) const
  {
    return !(operator==(i));
  }

  std::ostream & operator<<(std::ostream & os, const FeatureHandle & cons)
  {
    os  << "---------- FeatureHandle -----------------\n"
    << "RT: " << cons.getRT() << std::endl
    << "m/z: " << cons.getMZ() << std::endl
    << "Intensity: " << cons.getIntensity() << std::endl
    << "Map Index: " << cons.getMapIndex() << std::endl
    << "Element Id: " << cons.getUniqueId() << std::endl;
    return os;
  }
}
