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

#ifndef OPENMS_FORMAT_DATAACCESS_MSDATATRANSFORMINGCONSUMER_H
#define OPENMS_FORMAT_DATAACCESS_MSDATATRANSFORMINGCONSUMER_H

#include <OpenMS/INTERFACES/IMSDataConsumer.h>

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/ChromatogramPeak.h>

namespace OpenMS
{

    /**
      @brief Empty (NOP) function
    */
    OPENMS_DLLAPI extern void FunctionSpectrumNOP (MSSpectrum & /* s */);

    /**
      @brief Empty (NOP) function
    */
    OPENMS_DLLAPI extern void FunctionChromatogramNOP (MSChromatogram & /* c */);

    /**
      @brief Transforming consumer of MS data

      Is able to transform a spectra on the fly while it is read using a
      function pointer that can be set on the object.

      Note that the spectrum gets transformed in-place.
    */
    class OPENMS_DLLAPI MSDataTransformingConsumer :
      public Interfaces::IMSDataConsumer
    {

    public:

      /**
        @brief Constructor
      */
      MSDataTransformingConsumer()
      {
        sprocessing_ptr_ = &FunctionSpectrumNOP; // setting default processing action to noop
        cprocessing_ptr_ = &FunctionChromatogramNOP; // setting default processing action to noop
      }

      /// Default destructor
      ~MSDataTransformingConsumer() override { }

      void setExpectedSize(Size /* expectedSpectra */, Size /* expectedChromatograms */) override
      {
        // do nothing
      }

      void consumeSpectrum(SpectrumType & s) override
      {
        // apply the given function to it
        (*sprocessing_ptr_)(s);
      }

      virtual void setSpectraProcessingPtr( void (*sproptr)(SpectrumType&) )
      {
        sprocessing_ptr_ = sproptr;
      }

      void consumeChromatogram(ChromatogramType & c) override
      {
        // apply the given function to it
        (*cprocessing_ptr_)(c);
      }

      virtual void setChromatogramProcessingPtr( void (*cproptr)(ChromatogramType&) )
      {
        cprocessing_ptr_ = cproptr;
      }

      void setExperimentalSettings(const OpenMS::ExperimentalSettings&) override {}

    protected:
      void (*sprocessing_ptr_)(SpectrumType&);
      void (*cprocessing_ptr_)(ChromatogramType&);
    };

} //end namespace OpenMS

#endif // OPENMS_FORMAT_DATAACCESS_MSDATATRANSFORMINGCONSUMER_H
