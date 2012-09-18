// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/Peak2D.h>

namespace OpenMS
{
  char const * const Peak2D::dimension_name_short_[] =
  {
    "RT",
    "MZ"
  };

  char const * const Peak2D::dimension_name_full_[] =
  {
    "retention time",
    "mass-to-charge"
  };

  char const * const Peak2D::dimension_unit_short_[] =
  {
    "sec",
    "Th"
  };

  char const * const Peak2D::dimension_unit_full_[] =
  {
    "Seconds",
    "Thomson"
  };

  char const * Peak2D::shortDimensionName(UInt const dim)
  {
    return dimension_name_short_[dim];
  }

  char const * Peak2D::shortDimensionNameRT()
  {
    return dimension_name_short_[RT];
  }

  char const * Peak2D::shortDimensionNameMZ()
  {
    return dimension_name_short_[MZ];
  }

  char const * Peak2D::fullDimensionName(UInt const dim)
  {
    return dimension_name_full_[dim];
  }

  char const * Peak2D::fullDimensionNameRT()
  {
    return dimension_name_full_[RT];
  }

  char const * Peak2D::fullDimensionNameMZ()
  {
    return dimension_name_full_[MZ];
  }

  char const * Peak2D::shortDimensionUnit(UInt const dim)
  {
    return dimension_unit_short_[dim];
  }

  char const * Peak2D::shortDimensionUnitRT()
  {
    return dimension_unit_short_[RT];
  }

  char const * Peak2D::shortDimensionUnitMZ()
  {
    return dimension_unit_short_[MZ];
  }

  char const * Peak2D::fullDimensionUnit(UInt const dim)
  {
    return dimension_unit_full_[dim];
  }

  char const * Peak2D::fullDimensionUnitRT()
  {
    return dimension_unit_full_[RT];
  }

  char const * Peak2D::fullDimensionUnitMZ()
  {
    return dimension_unit_full_[MZ];
  }

  std::ostream & operator<<(std::ostream & os, const Peak2D & point)
  {
    os << "RT: " << point.getRT() <<  " MZ: "  << point.getMZ() << " INT: " << point.getIntensity();
    return os;
  }
} // namespace OpenMS
