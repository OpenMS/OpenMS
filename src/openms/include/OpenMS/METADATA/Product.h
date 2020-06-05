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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/CVTermList.h>

namespace OpenMS
{
  /**
      @brief Product meta information.

      This class describes the product isolation window for special scan types, such as MRM.

      @ingroup Metadata
  */
  class OPENMS_DLLAPI Product :
    public CVTermList
  {

public:

    /// Constructor
    Product() = default;
    /// Copy constructor
    Product(const Product &) = default;
    /// Move constructor
    Product(Product&&) = default;
    /// Destructor
    ~Product() override = default;

    /// Assignment operator
    Product & operator=(const Product &) = default;
    /// Move assignment operator
    Product& operator=(Product&&) & = default;

    /// Equality operator
    bool operator==(const Product & rhs) const;
    /// Equality operator
    bool operator!=(const Product & rhs) const;

    /// returns the target m/z
    double getMZ() const;
    /// sets the target m/z
    void setMZ(double mz);

    /// returns the lower offset from the target m/z
    double getIsolationWindowLowerOffset() const;
    /// sets the lower offset from the target m/z
    void setIsolationWindowLowerOffset(double bound);

    /// returns the upper offset from the target m/z
    double getIsolationWindowUpperOffset() const;
    /// sets the upper offset from the target m/z
    void setIsolationWindowUpperOffset(double bound);

protected:

    double mz_ = 0.0;
    double window_low_ = 0.0;
    double window_up_ = 0.0;
  };
} // namespace OpenMS

