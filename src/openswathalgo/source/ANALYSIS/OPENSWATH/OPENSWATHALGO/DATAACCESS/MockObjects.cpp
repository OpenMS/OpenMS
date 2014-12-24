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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/MockObjects.h>

#include <string>

namespace OpenSwath
{

  MockFeature::MockFeature()
  {
  }

  MockFeature::~MockFeature()
  {
  }

  void MockFeature::getRT(std::vector<double>& rt)
  {
    rt = m_rt_vec;
  }

  void MockFeature::getIntensity(std::vector<double>& intens)
  {
    intens = m_intensity_vec;
  }

  float MockFeature::getIntensity()
  {
    return m_intensity;
  }

  double MockFeature::getRT()
  {
    return m_rt;
  }

  MockMRMFeature::MockMRMFeature()
  {
  }

  MockMRMFeature::~MockMRMFeature()
  {
  }

  boost::shared_ptr<OpenSwath::IFeature> MockMRMFeature::getFeature(std::string nativeID)
  {
    return boost::static_pointer_cast<OpenSwath::IFeature>(m_features[nativeID]);
  }

  boost::shared_ptr<OpenSwath::IFeature> MockMRMFeature::getPrecursorFeature(std::string nativeID)
  {
    return boost::static_pointer_cast<OpenSwath::IFeature>(m_precursor_features[nativeID]);
  }

  std::vector<std::string> MockMRMFeature::getPrecursorIDs() const
  {
    std::vector<std::string> v;
    for (std::map<std::string, boost::shared_ptr<MockFeature> >::const_iterator 
        it = m_precursor_features.begin(); it != m_precursor_features.end(); ++it) 
    {
      v.push_back(it->first);
    }
    return v;
  }

  float MockMRMFeature::getIntensity()
  {
    return m_intensity;
  }

  double MockMRMFeature::getRT()
  {
    return m_rt;
  }

  size_t MockMRMFeature::size()
  {
    return m_features.size();
  }

  MockTransitionGroup::MockTransitionGroup()
  {
  }

  MockTransitionGroup::~MockTransitionGroup()
  {
  }

  std::size_t MockTransitionGroup::size()
  {
    return m_size;
  }

  std::vector<std::string> MockTransitionGroup::getNativeIDs()
  {
    return m_native_ids;
  }

  void MockTransitionGroup::getLibraryIntensities(std::vector<double>& intensities)
  {
    intensities = m_library_intensities;
  }

  MockSignalToNoise::MockSignalToNoise()
  {
  }

  double MockSignalToNoise::getValueAtRT(double /* RT */)
  {
    return m_sn_value;
  }

}
