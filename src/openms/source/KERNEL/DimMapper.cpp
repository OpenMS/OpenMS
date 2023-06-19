// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/DimMapper.h>

#include <OpenMS/DATASTRUCTURES/String.h>

#include <QLocale>

using namespace std;


namespace OpenMS
{
  namespace
  {
    DimMapper<1> dims({DIM_UNIT::RT});
    DimMapper<1> d({DIM_UNIT::RT});
    DimMapper<1> d2(d);
    bool x = (d == dims);
    Area<2> area(nullptr);
  }

  String DimBase::formattedValue(const ValueType value) const
  {
    // hint: QLocale::c().toString adds group separators to better visualize large numbers (e.g. 23.009.646.54,3)
    return String(this->getDimNameShort()) + ": " + QLocale::c().toString(value, 'f', valuePrecision());
  }

  String DimBase::formattedValue(const ValueType value, const String& prefix) const
  {
    return prefix + formattedValue(value);
  }

  int DimBase::valuePrecision() const
  {
    // decide on precision depending on unit; add more units if you have some intuition
    constexpr auto precision_for_unit = [](DIM_UNIT u) {
      switch (u)
      {
        case DIM_UNIT::RT:
        case DIM_UNIT::INT:
          return 2;
        case DIM_UNIT::MZ:
          return 8;
        default:
          return 4;
      }
    };
    return precision_for_unit(this->getUnit());
  }
} // namespace OpenMS
