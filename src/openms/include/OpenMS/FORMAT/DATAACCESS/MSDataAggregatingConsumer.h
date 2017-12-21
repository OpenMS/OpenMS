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

#ifndef OPENMS_FORMAT_DATAACCESS_MSDATAAGGREGATINGCONSUMER_H
#define OPENMS_FORMAT_DATAACCESS_MSDATAAGGREGATINGCONSUMER_H

#include <OpenMS/INTERFACES/IMSDataConsumer.h>

#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>

namespace OpenMS
{

    /**
      @brief Aggregates spectra by retention time

    */
    class OPENMS_DLLAPI MSDataAggregatingConsumer :
      public Interfaces::IMSDataConsumer
    {

      Interfaces::IMSDataConsumer* next_consumer_;
      double previous_rt_;
      bool rt_initialized_;
      SpectrumType s_tmp;
      std::vector<SpectrumType> s_list;

    public:

      /**
        @brief Constructor

        @note This does not transfer ownership of the consumer
      */
      MSDataAggregatingConsumer(Interfaces::IMSDataConsumer* next_consumer) :
        next_consumer_(next_consumer),
        previous_rt_(0.0),
        rt_initialized_(false)
      {}

      /**
        @brief Constructor

        Flushes data to next consumer

        @note It is essential to not delete the underlying next_consumer before
        deleting this object, otherwise we risk a memory error
      */
      ~MSDataAggregatingConsumer() override;

      void setExpectedSize(Size, Size) override {}

      void consumeSpectrum(SpectrumType & s) override;

      void consumeChromatogram(ChromatogramType & c) override;

      void setExperimentalSettings(const OpenMS::ExperimentalSettings&) override {}

    };

} //end namespace OpenMS

#endif // OPENMS_FORMAT_DATAACCESS_MSDATAAGGREGATINGCONSUMER_H

