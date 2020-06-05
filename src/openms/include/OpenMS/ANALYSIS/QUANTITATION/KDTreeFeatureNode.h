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

#pragma once

#include <OpenMS/config.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Macros.h>

namespace OpenMS
{

class KDTreeFeatureMaps;

/// A node of the kD-tree with pointer to corresponding data and index
class OPENMS_DLLAPI KDTreeFeatureNode
{

public:

  /// Constructor
  KDTreeFeatureNode(KDTreeFeatureMaps* data, Size idx);

  /// Copy constructor - copy the pointer, use same data object
  KDTreeFeatureNode(const KDTreeFeatureNode& rhs);

  /// Assignment operator - copy the pointer, use same data object
  KDTreeFeatureNode& operator=(KDTreeFeatureNode const& rhs);

  /// Destructor
  virtual ~KDTreeFeatureNode();

  /// libkdtree++ needs this typedef
  typedef double value_type;

  /// Needed for 2D range queries using libkdtree++. [0] returns RT, [1] m/z.
  value_type operator[](Size i) const;

  /// Return index of corresponding feature in data_
  Size getIndex() const;

protected:

  /// Pointer to the actual data
  KDTreeFeatureMaps* data_;

  /// Index of this feature
  Size idx_;

private:

  /// Default constructor is not supposed to be called
  KDTreeFeatureNode();

};

}

