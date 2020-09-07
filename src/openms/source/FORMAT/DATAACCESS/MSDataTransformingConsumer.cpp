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
// $Authors: Hannes Roest, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>

namespace OpenMS
{
  /**
        @brief Constructor
      */
      MSDataTransformingConsumer::MSDataTransformingConsumer()
       : lambda_spec_(nullptr),
         lambda_chrom_(nullptr),
         lambda_exp_settings_(nullptr)
      {
      }

      /// Default destructor
      MSDataTransformingConsumer::~MSDataTransformingConsumer()
      {
      }

      void MSDataTransformingConsumer::setExpectedSize(Size /* expectedSpectra */, Size /* expectedChromatograms */)
      {
      }

      void MSDataTransformingConsumer::consumeSpectrum(SpectrumType& s)
      {
        // apply the given function to it (unless nullptr)
        if (lambda_spec_) lambda_spec_(s);
      }

      void MSDataTransformingConsumer::setSpectraProcessingFunc( std::function<void (SpectrumType&)> f_spec )
      {
        lambda_spec_ = f_spec;
      }

      void MSDataTransformingConsumer::consumeChromatogram(ChromatogramType & c)
      {
        // apply the given function to it (unless nullptr)
        if (lambda_chrom_) lambda_chrom_(c);
      }

      void MSDataTransformingConsumer::setChromatogramProcessingFunc( std::function<void (ChromatogramType&)> f_chrom )
      {
        lambda_chrom_ = f_chrom;
      }
      
      void MSDataTransformingConsumer::setExperimentalSettingsFunc( std::function<void (const OpenMS::ExperimentalSettings&)> f_exp_settings )
      {
        lambda_exp_settings_ = f_exp_settings;
      }

      void MSDataTransformingConsumer::setExperimentalSettings(const OpenMS::ExperimentalSettings& es)
      { 
        // apply the given function to it (unless nullptr)
        if (lambda_exp_settings_) lambda_exp_settings_(es);        
      }
} // namespace OpenMS
