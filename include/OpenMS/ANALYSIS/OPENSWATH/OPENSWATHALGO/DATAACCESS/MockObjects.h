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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENSWATH_DATAACCESS_MOCKOBJECTS_H_
#define OPENSWATH_DATAACCESS_MOCKOBJECTS_H_

#include "vector"
#include "map"
#include "boost/shared_ptr.hpp"
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ITransition.h>
namespace OpenSwath
{

  class MockFeature :
    public OpenSwath::IFeature
  {
public:

    MockFeature()
    {
    }

    ~MockFeature()
    {
    }

    void getRT(std::vector<double> & rt)
    {
      rt = m_rt_vec;
    }

    void getIntensity(std::vector<double> & intens)
    {
      intens = m_intensity_vec;
    }

    float getIntensity()
    {
      return m_intensity;
    }

    double getRT()
    {
      return m_rt;
    }

    std::vector<double> m_rt_vec;
    std::vector<double> m_intensity_vec;
    float m_intensity;
    double m_rt;
  };

  class MockMRMFeature :
    public OpenSwath::IMRMFeature
  {
public:

    MockMRMFeature() {}

    ~MockMRMFeature() {}

    boost::shared_ptr<OpenSwath::IFeature> getFeature(std::string nativeID)
    {
      return boost::static_pointer_cast<OpenSwath::IFeature>(m_features[nativeID]);
    }

    float getIntensity()
    {
      return m_intensity;
    }

    double getRT()
    {
      return m_rt;
    }

    std::map<std::string, boost::shared_ptr<MockFeature> > m_features;
    float m_intensity;
    double m_rt;
  };

  class MockTransitionGroup :
    public OpenSwath::ITransitionGroup
  {
public:

    MockTransitionGroup()
    {
    }

    ~MockTransitionGroup()
    {
    }

    std::size_t size()
    {
      return m_size;
    }

    std::vector<std::string> getNativeIDs()
    {
      return m_native_ids;
    }

    void getLibraryIntensities(std::vector<double> & intensities)
    {
      intensities = m_library_intensities;
    }

    std::size_t m_size;
    std::vector<std::string> m_native_ids;
    std::vector<double> m_library_intensities;
  };

  class MockSignalToNoise :
    public OpenSwath::ISignalToNoise
  {
public:
    MockSignalToNoise()
    {
    }

    double getValueAtRT(double RT)
    {
      // to prevent unused RT warning
      return m_sn_value * RT / RT;
    }

    double m_sn_value;
  };

} //end namespace OpenMS

#endif
