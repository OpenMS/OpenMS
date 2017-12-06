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
// $Maintainer: Mathias Walzer $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_NORMALIZER_H
#define OPENMS_FILTERING_TRANSFORMERS_NORMALIZER_H

#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Normalizes the peak intensities spectrum-wise.

    Either to a total intensity-sum of one (i.e. to total-ion-count; TIC) or to a maximum intensity of one.

    @htmlinclude OpenMS_Normalizer.parameters

    @ingroup SpectraPreprocessers
  */
  class OPENMS_DLLAPI Normalizer :
    public DefaultParamHandler
  {
public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    Normalizer();
    /// destructor
    ~Normalizer() override;

    /// assignment operator
    Normalizer & operator=(const Normalizer & source);
    /// copy constructor
    Normalizer(const Normalizer & source);

    // @}

    // @name Accessors
    // @{

    /**
      Workhorse of this class.

      @param spectrum Input/output spectrum containing peaks
      @throws Exception::InvalidValue if 'method_' has unknown value
    */
    template <typename SpectrumType>
    void filterSpectrum(SpectrumType& spectrum) const
    {
      if (spectrum.empty()) return;

      typedef typename SpectrumType::Iterator Iterator;
      typedef typename SpectrumType::ConstIterator ConstIterator;

      double divisor(0);
      // find divisor      
      if (method_ == "to_one")
      { // normalizes the max peak to 1 and the remaining peaks to values relative to max
        divisor = spectrum.begin()->getIntensity(); // safety measure: if all intensities are negative, divisor would stay 0 (as constructed)
        for (ConstIterator it = spectrum.begin(); it != spectrum.end(); ++it)
        {
          if (divisor < it->getIntensity()) divisor = it->getIntensity();
        }
      }
      else if (method_ == "to_TIC")
      { // normalizes the peak intensities to the TIC
        for (ConstIterator it = spectrum.begin(); it != spectrum.end(); ++it)
        {
          divisor += it->getIntensity();
        }
      }
      // method unknown
      else
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Method not known", method_);
      }

      // normalize
      for (Iterator it = spectrum.begin(); it != spectrum.end(); ++it)
      {
        it->setIntensity(it->getIntensity() / divisor);
      }

      return;

    }

    ///
    void filterPeakSpectrum(PeakSpectrum & spectrum) const;
    ///
    void filterPeakMap(PeakMap & exp) const;

    void updateMembers_() override;

    // @}

private:
    String method_;
  };


}
#endif //OPENMS_FILTERING_TRANSFORMERS_NORMALIZER_H
