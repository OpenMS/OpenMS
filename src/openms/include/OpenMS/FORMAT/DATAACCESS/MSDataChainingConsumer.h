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

#ifndef OPENMS_FORMAT_DATAACCESS_MSDATACHAININGCONSUMER_H
#define OPENMS_FORMAT_DATAACCESS_MSDATACHAININGCONSUMER_H

#include <OpenMS/INTERFACES/IMSDataConsumer.h>

#include <vector>

namespace OpenMS
{

  /**
    @brief Consumer class that passes all consumed data through a set of operations

    This consumer allows to chain multiple data consumers and applying them in
    a pre-specified order. This can be useful if a certain operation on a
    dataset needs to be performed but some pre-processing (data reduction etc.)
    or post-processing (writing to disk, caching on disk). The different
    processing steps can be added to the chaining consumer (in the correct
    order) without knowledge of the specific pre/post processing steps.

    Usage:

    @code
    MSDataTransformingConsumer * transforming_consumer_first = new MSDataTransformingConsumer(); // apply some transformation
    MSDataTransformingConsumer * transforming_consumer_second = new MSDataTransformingConsumer(); // apply second transformation
    MSDataWritingConsumer * writing_consumer = new MSDataWritingConsumer(outfile); // writing to disk

    std::vector<Interfaces::IMSDataConsumer *> consumer_list;
    consumer_list.push_back(transforming_consumer_first);
    consumer_list.push_back(transforming_consumer_second);
    consumer_list.push_back(writing_consumer);
    MSDataChainingConsumer * chaining_consumer = new MSDataChainingConsumer(consumer_list);

    // now chaining_consumer can be passed to a function expecting a IMSDataConsumer interface
    @endcode

  */
  class OPENMS_DLLAPI MSDataChainingConsumer :
    public Interfaces::IMSDataConsumer
  {
    std::vector<Interfaces::IMSDataConsumer *> consumers_;

  public:

    /**
     * @brief Default Constructor
     *
     */
    MSDataChainingConsumer();

    /**
     * @brief Constructor
     *
     * Pass a list of consumers that should be called sequentially
     *
     * @note By default, this does not transfers ownership - it is the callers
     * responsibility to delete the pointer to consumer afterwards.
     *
     */
    MSDataChainingConsumer(std::vector<Interfaces::IMSDataConsumer *> consumers);

    /**
     * @brief Destructor
     *
     * Does nothing. Does not destroy underlying consumers, therefore is the
     * responsibility of the caller to destroy all consumers.
     *
     */
    ~MSDataChainingConsumer() override;

    /**
     * @brief Append a consumer to the chain of consumers to be executed
     *
     * @note This does not transfer ownership - it is the callers
     * responsibility to delete the pointer to consumer afterwards.
     *
     */
    void appendConsumer(Interfaces::IMSDataConsumer * consumer);

    /**
     * @brief Set experimental settings for all consumers
     *
     * Will set the experimental settings for all chained consumers
     *
     */
    void setExperimentalSettings(const ExperimentalSettings & settings) override;

    /**
     * @brief Set expected size for all consumers
     *
     * Will set the expected size for all chained consumers
     *
     */
    void setExpectedSize(Size s_size, Size c_size) override;

    /**
     * @brief Call all consumers in the specified order for the given spectrum
     *
     */
    void consumeSpectrum(SpectrumType & s) override;

    /**
     * @brief Call all consumers in the specified order for the given chromatogram
     *
     */
    void consumeChromatogram(ChromatogramType & c) override;

  };

} //end namespace OpenMS

#endif // OPENMS_FORMAT_DATAACCESS_MSDATACHAININGCONSUMER_H

