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

#pragma once

#include <OpenMS/INTERFACES/IMSDataConsumer.h>

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/ChromatogramPeak.h>

namespace OpenMS
{

    /**
      @brief Transforming consumer of MS data

      Is able to work on spectrum/chromatogram on the fly while it is read using a
      lambda function provided by the user.
      The lambda can transform the spectrum/chromatogram in-place or just extract some
      information and store it in locally captured variables. 
      Just make sure the captured variables are still in scope when this consumer is used.

    */
    class OPENMS_DLLAPI MSDataTransformingConsumer :
      public Interfaces::IMSDataConsumer
    {

    public:

      /**
        @brief Constructor
      */
      MSDataTransformingConsumer();

      /// Default destructor
      ~MSDataTransformingConsumer() override;

      /// ignored
      void setExpectedSize(Size /* expectedSpectra */, Size /* expectedChromatograms */) override;

      void consumeSpectrum(SpectrumType& s) override;

      /**
        @brief Sets the lambda function to be called for every spectrum which is passed to this interface

        The lambda can contain locally captured variables, which allow to change some state. Make sure
        that the captured variables are still in scope (i.e. not destroyed at the point of calling the lambda).
        Pass a nullptr if the spectrum should be left unchanged.
      */
      virtual void setSpectraProcessingFunc( std::function<void (SpectrumType&)> f_spec );

      void consumeChromatogram(ChromatogramType& c) override;

      /**
        @brief Sets the lambda function to be called for every chromatogram which is passed to this interface

        The lambda can contain locally captured variables, which allow to change some state. Make sure
        that the captured variables are still in scope (i.e. not destroyed at the point of calling the lambda).
        Pass a nullptr if the chromatogram should be left unchanged.
      */
      virtual void setChromatogramProcessingFunc( std::function<void (ChromatogramType&)> f_chrom );

      /**
        @brief Sets the lambda function to be called when setExperimentalSettings is called via this interface

        The lambda can contain locally captured variables, which allow to change some state. Make sure
        that the captured variables are still in scope (i.e. not destroyed at the point of calling the lambda).
        Pass a nullptr if nothing should happen (default).
      */
      virtual void setExperimentalSettingsFunc( std::function<void (const OpenMS::ExperimentalSettings&)> f_exp_settings );

      /// ignored
      void setExperimentalSettings(const OpenMS::ExperimentalSettings&) override;

    protected:
      std::function<void (SpectrumType&)> lambda_spec_;
      std::function<void (ChromatogramType&)> lambda_chrom_;
      std::function<void (const OpenMS::ExperimentalSettings&)> lambda_exp_settings_;
    };

} //end namespace OpenMS

