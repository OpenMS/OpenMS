// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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

#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/ID/ProcessingSoftware.h>
#include <OpenMS/METADATA/ID/InputFile.h>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /** @brief Data processing step that is applied to the data (e.g. database search, PEP calculation, filtering, ConsensusID).
    */
    struct ProcessingStep: public MetaInfoInterface
    {
      ProcessingSoftwareRef software_ref;

      std::vector<InputFileRef> input_file_refs;

      DateTime date_time;

      // @TODO: add processing actions that are relevant for ID data
      std::set<DataProcessing::ProcessingAction> actions;

      explicit ProcessingStep(
        ProcessingSoftwareRef software_ref,
        const std::vector<InputFileRef>& input_file_refs = std::vector<InputFileRef>(),
        const DateTime& date_time = DateTime::now(),
        const std::set<DataProcessing::ProcessingAction>& actions = std::set<DataProcessing::ProcessingAction>())
        :
        software_ref(software_ref), input_file_refs(input_file_refs),
        date_time(date_time), actions(actions)
      {
      }

      ProcessingStep(const ProcessingStep& other) = default;

      // order by date/time first, don't compare meta data (?):
      bool operator<(const ProcessingStep& other) const
      {
        return (std::tie(date_time, software_ref, input_file_refs, actions) <
                std::tie(other.date_time, other.software_ref,
                         other.input_file_refs, other.actions));
      }

      // don't compare meta data (?):
      bool operator==(const ProcessingStep& other) const
      {
        return (std::tie(software_ref, input_file_refs, date_time, actions) ==
                std::tie(other.software_ref, other.input_file_refs,
                         other.date_time, other.actions));
      }
    };

    typedef std::set<ProcessingStep> ProcessingSteps;
    typedef IteratorWrapper<ProcessingSteps::iterator> ProcessingStepRef;

  }
}
