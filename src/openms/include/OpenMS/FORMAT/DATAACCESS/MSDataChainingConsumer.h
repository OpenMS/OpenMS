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

#ifndef OPENMS_FORMAT_DATAACCESS_MSDATACHAININGCONSUMER_H
#define OPENMS_FORMAT_DATAACCESS_MSDATACHAININGCONSUMER_H

#include <OpenMS/INTERFACES/IMSDataConsumer.h>

namespace OpenMS
{

  /**
    @brief Consumer class that passes all operations on to a set of consumers

  */
  class OPENMS_DLLAPI MSDataChainingConsumer :
    public Interfaces::IMSDataConsumer< MSExperiment<> >
  {
    std::vector<IMSDataConsumer<> *> consumers_;

  public:

    /**
     * @brief Default Constructor
     *
     */
    MSDataChainingConsumer() {}

    /**
     * @brief Constructor
     *
     * Pass a list of consumers that should be called sequentially
     *
     */
    MSDataChainingConsumer(std::vector<IMSDataConsumer<> *> consumers) :
      consumers_(consumers)
    {}

    /**
     * @brief Destructor
     *
     */
    ~MSDataChainingConsumer() {}

    /**
     * @brief Append a consumer to the chain of consumers to be executed
     *
     * @note This does not transfers ownership - it is your responsibility to
     * delete the pointer to consumer afterwards.
     *
     */
    void appendConsumer(IMSDataConsumer<> * consumer)
    {
      consumers_.push_back(consumer);
    }

    void setExperimentalSettings(const ExperimentalSettings & settings)
    {
      for (Size i = 0; i < consumers_.size(); i++)
      {
        consumers_[i]->setExperimentalSettings(settings);
      }
    }

    void setExpectedSize(Size s_size, Size c_size) 
    {
      for (Size i = 0; i < consumers_.size(); i++)
      {
        consumers_[i]->setExpectedSize(s_size, c_size);
      }
    }

    void consumeSpectrum(SpectrumType & s)
    {
      for (Size i = 0; i < consumers_.size(); i++)
      {
        consumers_[i]->consumeSpectrum(s);
      }
    }

    void consumeChromatogram(ChromatogramType & c)
    {
      for (Size i = 0; i < consumers_.size(); i++)
      {
        consumers_[i]->consumeChromatogram(c);
      }
    }

  };
} //end namespace OpenMS

#endif
