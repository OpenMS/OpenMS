// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#ifndef OPENMS_ANALYSIS_OPENSWATH_DATAACCESS_SPECTRUMACCESSQUADMZTRANSFORMING_H
#define OPENMS_ANALYSIS_OPENSWATH_DATAACCESS_SPECTRUMACCESSQUADMZTRANSFORMING_H


#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessTransforming.h>

namespace OpenMS
{
  /**
   * @brief A transforming m/z wrapper around spectrum access using a quadratic equation. 
   *
   * For each spectrum access, each m/z value is transformed using the equation 
   *    mz = a + b * mz + c * mz^2
   * 
   * This can be used to implement an on-line mass correction for TOF
   * instruments (for example).
   *
   */
  class OPENMS_DLLAPI SpectrumAccessQuadMZTransforming :
    public SpectrumAccessTransforming
  {
public:

    /** @brief Constructor
     *
     * @param a Regression parameter 0
     * @param b Regression parameter 1
     * @param c Regression parameter 2
     * @param ppm Whether the transformation should be applied in ppm domain
     *            (if false, it is applied directly in m/z domain)
     *
    */
    explicit SpectrumAccessQuadMZTransforming(OpenSwath::SpectrumAccessPtr sptr,
        double a, double b, double c, bool ppm);
        
    ~SpectrumAccessQuadMZTransforming() override;

    boost::shared_ptr<OpenSwath::ISpectrumAccess> lightClone() const override;

    OpenSwath::SpectrumPtr getSpectrumById(int id) override;

private:

    double a_;
    double b_;
    double c_;
    bool ppm_;

  };
}

#endif // OPENMS_ANALYSIS_OPENSWATH_DATAACCESS_SPECTRUMACCESSQUADMZTRANSFORMING_H
