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
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <unordered_map>

using namespace std;
namespace OpenMS
{

  IDScoreSwitcherAlgorithm::IDScoreSwitcherAlgorithm() :
      IDScoreSwitcherAlgorithm::DefaultParamHandler("IDScoreSwitcherAlgorithm")
  {
    defaults_.setValue("new_score", "", "Name of the meta value to use as the new score");
    defaults_.setValue("new_score_orientation", "", "Orientation of the new score (are higher or lower values better?)");
    defaults_.setValidStrings("new_score_orientation", {"lower_better","higher_better"});
    defaults_.setValue("new_score_type", "", "Name to use as the type of the new score (default: same as 'new_score')");
    defaults_.setValue("old_score", "", "Name to use for the meta value storing the old score (default: old score type)");
    defaults_.setValue("proteins", "false", "Apply to protein scores instead of PSM scores");
    defaults_.setValidStrings("proteins", {"true","false"});
    defaultsToParam_();
    updateMembers_();
  }

  void IDScoreSwitcherAlgorithm::updateMembers_()
  {
    new_score_ = param_.getValue("new_score");
    new_score_type_ = param_.getValue("new_score_type");
    old_score_ = param_.getValue("old_score");
    higher_better_ = (param_.getValue("new_score_orientation").toString() ==
                      "higher_better");

    if (new_score_type_.empty()) new_score_type_ = new_score_;
  }

} // namespace OpenMS
