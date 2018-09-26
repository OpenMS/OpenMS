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

#ifndef OPENMS_METADATA_ID_SCORETYPE_H
#define OPENMS_METADATA_ID_SCORETYPE_H

#include <OpenMS/METADATA/ID/MetaData.h>

#include <boost/optional.hpp>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /*!
      Information about a score type.
    */
    struct ScoreType: public MetaInfoInterface
    {
      CVTerm cv_term;

      String name;

      bool higher_better;

      // reference to the software that assigned the score:
      boost::optional<ProcessingSoftwareRef> software_opt;
      // @TODO: scores assigned by different software tools/versions are
      // considered as different scores (even if they have the same name) -
      // does that make sense?

      ScoreType():
        higher_better(true), software_opt()
      {
      }

      explicit ScoreType(const CVTerm& cv_term, bool higher_better,
                         boost::optional<ProcessingSoftwareRef> software_opt =
                         boost::none):
        cv_term(cv_term), name(cv_term.getName()), higher_better(higher_better),
        software_opt(software_opt)
      {
      }

      explicit ScoreType(const String& name, bool higher_better,
                         boost::optional<ProcessingSoftwareRef> software_opt =
                         boost::none):
        cv_term(), name(name), higher_better(higher_better),
        software_opt(software_opt)
      {
      }

      ScoreType(const ScoreType& other) = default;

      // don't include "higher_better" in the comparison:
      bool operator<(const ScoreType& other) const
      {
        return (std::tie(cv_term.getAccession(), name, software_opt) <
                std::tie(other.cv_term.getAccession(), other.name,
                         other.software_opt));
      }

      // don't include "higher_better" in the comparison:
      bool operator==(const ScoreType& other) const
      {
        return (std::tie(cv_term.getAccession(), name, software_opt) ==
                std::tie(other.cv_term.getAccession(), other.name,
                         other.software_opt));
      }
    };

    typedef std::set<ScoreType> ScoreTypes;
    typedef IteratorWrapper<ScoreTypes::iterator> ScoreTypeRef;

    // @TODO: use a "boost::multi_index_container" to allow efficient access in
    // sequence and by key?
    typedef std::vector<std::pair<ScoreTypeRef, double>> ScoreList;

  }
}

#endif
