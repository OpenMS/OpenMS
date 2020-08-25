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
// $Maintainer: Johannes Veit $
// $Authors: Johannes Veit $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureNode.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>

namespace OpenMS
{

KDTreeFeatureNode::KDTreeFeatureNode(KDTreeFeatureMaps* data, Size idx) :
  data_(data),
  idx_(idx)
{
}

KDTreeFeatureNode::KDTreeFeatureNode(const KDTreeFeatureNode& rhs) :
  data_(rhs.data_),
  idx_(rhs.idx_)
{
}

KDTreeFeatureNode& KDTreeFeatureNode::operator=(KDTreeFeatureNode const& rhs)
{
  data_ = rhs.data_;
  idx_ = rhs.idx_;

  return *this;
}

KDTreeFeatureNode::~KDTreeFeatureNode()
{
}

Size KDTreeFeatureNode::getIndex() const
{
  return idx_;
}

KDTreeFeatureNode::value_type KDTreeFeatureNode::operator[](Size i) const
{
  if (i == 0)
  {
    return data_->rt(idx_);
  }
  else if (i == 1)
  {
    return data_->mz(idx_);
  }
  else
  {
    const String& err_msg = "Indices other than 0 (RT) and 1 (m/z) are not allowed!";
    throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, err_msg);
  }
}

}
