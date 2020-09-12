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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ID/MetaData.h>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/optional.hpp>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /*!
      @brief Representation of a search query, e.g. spectrum or feature.
    */
    struct DataQuery: public MetaInfoInterface
    {
      /// Spectrum or feature ID (from the file referenced by @t input_file_opt)
      String data_id;

      /// (Optional) reference to the input file
      boost::optional<InputFileRef> input_file_opt;
      // @TODO: make this non-optional (i.e. required)?

      double rt, mz; //< Position

      /// Constructor
      explicit DataQuery(
        const String& data_id,
        boost::optional<InputFileRef> input_file_opt = boost::none,
        double rt = std::numeric_limits<double>::quiet_NaN(),
        double mz = std::numeric_limits<double>::quiet_NaN()):
        data_id(data_id), input_file_opt(input_file_opt), rt(rt), mz(mz)
      {
      }

      /// Merge in data from another object
      DataQuery& operator+=(const DataQuery& other)
      {
        // merge meta info - existing entries may be overwritten:
        std::vector<UInt> keys;
        other.getKeys(keys);
        for (const UInt key : keys)
        {
          setMetaValue(key, other.getMetaValue(key));
        }
        rt = other.rt;
        mz = other.mz;
        return *this;
      }

      // @TODO: do we need an "experiment label" (used e.g. in pepXML)?
      // if yes, should it be stored here or together with the input file?
    };

    // combination of input file and data ID must be unique:
    typedef boost::multi_index_container<
      DataQuery,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<
          boost::multi_index::composite_key<
            DataQuery,
            boost::multi_index::member<DataQuery, boost::optional<InputFileRef>,
                                       &DataQuery::input_file_opt>,
            boost::multi_index::member<DataQuery, String,
                                       &DataQuery::data_id>>>>
      > DataQueries;
    typedef IteratorWrapper<DataQueries::iterator> DataQueryRef;

  }
}
