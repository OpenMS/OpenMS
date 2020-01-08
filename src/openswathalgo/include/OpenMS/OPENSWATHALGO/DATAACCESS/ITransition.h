// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

#include <OpenMS/OPENSWATHALGO/OpenSwathAlgoConfig.h>

namespace OpenSwath
{
  // Datastructures for Scoring
  class OPENSWATHALGO_DLLAPI IFeature
  {
public:
    virtual ~IFeature(){}
    virtual void getRT(std::vector<double>& rt) const = 0;
    virtual void getIntensity(std::vector<double>& intens) const = 0;
    virtual float getIntensity() const = 0;
    virtual double getRT() const = 0;
  };

  class OPENSWATHALGO_DLLAPI IMRMFeature
  {
public:
    virtual ~IMRMFeature(){}
    virtual boost::shared_ptr<OpenSwath::IFeature> getFeature(std::string nativeID) = 0;
    virtual boost::shared_ptr<OpenSwath::IFeature> getPrecursorFeature(std::string nativeID) = 0;
    virtual std::vector<std::string> getNativeIDs() const = 0;
    virtual std::vector<std::string> getPrecursorIDs() const = 0;
    virtual float getIntensity() const = 0;
    virtual double getRT() const = 0;
    virtual size_t size() const = 0;
  };

  struct OPENSWATHALGO_DLLAPI ITransitionGroup
  {
    virtual ~ITransitionGroup() {}
    virtual std::size_t size() const = 0;
    virtual std::vector<std::string> getNativeIDs() const = 0;
    virtual void getLibraryIntensities(std::vector<double>& intensities) const = 0;
  };

  struct OPENSWATHALGO_DLLAPI ISignalToNoise
  {
    virtual ~ISignalToNoise() {}
    virtual double getValueAtRT(double RT) = 0; // cannot be const due to OpenMS implementation
  };
  typedef boost::shared_ptr<ISignalToNoise> ISignalToNoisePtr;


} //end Namespace OpenSwath

