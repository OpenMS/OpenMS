// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

namespace OpenMS
{
    /**
      @brief Empty (NOP) function
    */
    void FunctionSpectrumNOP (MSSpectrum<Peak1D> & /* s */) {;}

    /**
      @brief Empty (NOP) function
    */
    void FunctionChromatogramNOP (MSChromatogram<ChromatogramPeak> & /* s */) {;}

    /**
      @brief Transforming consumer of MS data

      Is able to transform a spectra on the fly while it is read using a
      function pointer that can be set on the object.
    */
    class OPENMS_DLLAPI MSDataTransformingConsumer : 
      public Interfaces::IMSDataConsumer<>
    {

    public:
      typedef MSExperiment<> MapType;
      typedef MapType::SpectrumType SpectrumType;
      typedef MapType::ChromatogramType ChromatogramType;

      /**
        @brief Constructor
      */
      MSDataTransformingConsumer() 
      { 
        sprocessing_ptr_ = &FunctionSpectrumNOP; // setting default processing action to noop
        cprocessing_ptr_ = &FunctionChromatogramNOP; // setting default processing action to noop
      }

      /// Default destructor
      virtual ~MSDataTransformingConsumer() { }

      virtual void setExpectedSize(Size /* expectedSpectra */, Size /* expectedChromatograms */)
      {
        // do nothing
      }

      virtual void consumeSpectrum(SpectrumType & s)
      {
        // apply the given function to it
        (*sprocessing_ptr_)(s);
      }

      virtual void setSpectraProcessingPtr( void (*sproptr)(SpectrumType&) )
      {
        sprocessing_ptr_ = sproptr;
      }

      virtual void consumeChromatogram(ChromatogramType & c)
      {
        // apply the given function to it
        (*cprocessing_ptr_)(c);
      }

      virtual void setChromatogramProcessingPtr( void (*cproptr)(ChromatogramType&) )
      {
        cprocessing_ptr_ = cproptr;
      }

      virtual void setExperimentalSettings(OpenMS::ExperimentalSettings&) {};

    protected:
      void (*sprocessing_ptr_)(SpectrumType&);
      void (*cprocessing_ptr_)(ChromatogramType&);
    };

} //end namespace OpenMS

#endif // OPENMS_FORMAT_DATAACCESS_MSDATATRANSFORMINGCONSUMER_H
