// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/INTERFACES/IMSDataConsumer.h>

#include <vector>
#include <memory>

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


