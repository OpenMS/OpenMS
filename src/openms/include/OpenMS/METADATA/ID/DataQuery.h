// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

#include <boost/optional.hpp>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /** @brief Search query, e.g. spectrum or feature.
    */
    struct DataQuery: public MetaInfoInterface
    {
      /// spectrum or feature ID (from the file referenced by "input_file_ref"):
      String data_id;

      // @TODO: make this non-optional (i.e. required)?
      boost::optional<InputFileRef> input_file_opt;

      double rt, mz; // position

      explicit DataQuery(
        const String& data_id,
        boost::optional<InputFileRef> input_file_opt = boost::none,
        double rt = std::numeric_limits<double>::quiet_NaN(),
        double mz = std::numeric_limits<double>::quiet_NaN()):
        data_id(data_id), input_file_opt(input_file_opt), rt(rt), mz(mz)
      {
      }

      DataQuery(const DataQuery& other) = default;

      // ignore RT and m/z for comparisons to avoid issues with rounding:
      bool operator<(const DataQuery& other) const
      {
        // can't compare references directly, so compare addresses:
        const String* sp = input_file_opt ? &(**input_file_opt) : nullptr;
        const String* o_sp = other.input_file_opt ? &(**other.input_file_opt) :
          nullptr;
        return std::tie(sp, data_id) < std::tie(o_sp, other.data_id);
      }

      // ignore RT and m/z for comparisons to avoid issues with rounding:
      bool operator==(const DataQuery& other) const
      {
        return std::tie(input_file_opt, data_id) ==
          std::tie(other.input_file_opt, other.data_id);
      }

      // @TODO: do we need an "experiment label" (used e.g. in pepXML)?
      // if yes, should it be stored here or together with the input file?
    };

    typedef std::set<DataQuery> DataQueries;
    typedef IteratorWrapper<DataQueries::iterator> DataQueryRef;

  }
}
