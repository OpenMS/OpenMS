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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessQuadMZTransforming.h>

namespace OpenMS
{

  SpectrumAccessQuadMZTransforming::SpectrumAccessQuadMZTransforming(
      OpenSwath::SpectrumAccessPtr sptr,
      double a, double b, double c, bool ppm) :
        SpectrumAccessTransforming(sptr), 
        a_(a), 
        b_(b), 
        c_(c), 
        ppm_(ppm)
    {}
        
    SpectrumAccessQuadMZTransforming::~SpectrumAccessQuadMZTransforming() {}

    boost::shared_ptr<OpenSwath::ISpectrumAccess> SpectrumAccessQuadMZTransforming::lightClone() const
    {
      // Create a light clone of *this by initializing a new
      // SpectrumAccessQuadMZTransforming with a light clone of the underlying
      // SpectrumAccess object and the parameters.
      return boost::shared_ptr<SpectrumAccessQuadMZTransforming>(
          new SpectrumAccessQuadMZTransforming(sptr_->lightClone(), a_, b_, c_, ppm_));
    }

    OpenSwath::SpectrumPtr SpectrumAccessQuadMZTransforming::getSpectrumById(int id)
    {
      OpenSwath::SpectrumPtr s = sptr_->getSpectrumById(id);
      for (size_t i = 0; i < s->getMZArray()->data.size(); i++)
      {
        // mz = a + b * mz + c * mz^2
        double predict = 
          a_ + 
          b_ * s->getMZArray()->data[i] +
          c_ * s->getMZArray()->data[i] * s->getMZArray()->data[i];

        // If ppm is true, we predicted the ppm deviation, not the actual new mass
        if (ppm_)
        {
          s->getMZArray()->data[i] = s->getMZArray()->data[i] - predict*s->getMZArray()->data[i]/1000000;
        }
        else
        {
          s->getMZArray()->data[i] = predict;
        }
      }
      return s;
    }

}
