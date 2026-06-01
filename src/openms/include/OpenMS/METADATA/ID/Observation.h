// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ID/InputFile.h>
#include <OpenMS/METADATA/ID/MetaData.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/member.hpp>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /*!
      @brief Representation of an observation, e.g. a spectrum or feature, in an input data file.
    */
    struct Observation: public MetaInfoInterface
    {
      /// Spectrum or feature ID (from the file referenced by @p input_file)
      String data_id;

      /// Reference to the input file
      InputFileRef input_file;

      double rt, mz; //< Position

      /// Constructor
      explicit Observation(
        const String& data_id,
        const InputFileRef& input_file,
        double rt = std::numeric_limits<double>::quiet_NaN(),
        double mz = std::numeric_limits<double>::quiet_NaN()):
        data_id(data_id), input_file(input_file), rt(rt), mz(mz)
      {
      }

      /// Merge in data from another object
      Observation& merge(const Observation& other)
      {
        // merge meta info - existing entries may be overwritten:
        addMetaValues(other);
        rt = other.rt;
        mz = other.mz;
        return *this;
      }
    };

    // combination of input file and data ID must be unique:
    typedef boost::multi_index_container<
      Observation,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<
          boost::multi_index::composite_key<
            Observation,
            boost::multi_index::member<Observation, InputFileRef,
                                       &Observation::input_file>,
            boost::multi_index::member<Observation, String,
                                       &Observation::data_id>>>>
      > Observations;
    typedef IteratorWrapper<Observations::iterator> ObservationRef;
  }
}
