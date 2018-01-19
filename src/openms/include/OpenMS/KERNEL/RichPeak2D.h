// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_RICHPEAK2D_H
#define OPENMS_KERNEL_RICHPEAK2D_H

#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/CONCEPT/UniqueIdInterface.h>

namespace OpenMS
{

  /**
    @brief A 2-dimensional raw data point or peak with meta information.

    This data structure is intended for continuous data or peak data.
    If you do not need to annotated single peaks with meta data, use Peak2D instead.

    @ingroup Kernel
  */
  class OPENMS_DLLAPI RichPeak2D :
    public Peak2D,
    public MetaInfoInterface,
    public UniqueIdInterface
  {
public:

    /// Default constructor
    RichPeak2D() :
      Peak2D(),
      MetaInfoInterface(),
      UniqueIdInterface()
    {}

    /// Copy constructor
    RichPeak2D(const RichPeak2D& p) :
      Peak2D(p),
      MetaInfoInterface(p),
      UniqueIdInterface(p)
    {}

    /// Constructor from Peak2D
    explicit RichPeak2D(const Peak2D& p) :
      Peak2D(p),
      MetaInfoInterface()
    {
      UniqueIdInterface::clearUniqueId();
    }

    /// Member constructor
    explicit RichPeak2D(const PositionType& pos, const IntensityType in) :
      Peak2D(pos, in),
      MetaInfoInterface()
    {}

    /// Destructor
    ~RichPeak2D() override
    {}

    /// Assignment operator
    RichPeak2D & operator=(const RichPeak2D& rhs)
    {
      if (this == &rhs) return *this;

      Peak2D::operator=(rhs);
      MetaInfoInterface::operator=(rhs);
      UniqueIdInterface::operator=(rhs);

      return *this;
    }

    /// Assignment operator
    RichPeak2D & operator=(const Peak2D& rhs)
    {
      if (this == &rhs) return *this;

      Peak2D::operator=(rhs);
      clearMetaInfo();
      UniqueIdInterface::clearUniqueId();

      return *this;
    }

    /// Equality operator
    bool operator==(const RichPeak2D& rhs) const
    {
      return Peak2D::operator==(rhs) &&
             MetaInfoInterface::operator==(rhs) &&
             UniqueIdInterface::operator==(rhs);
    }

    /// Equality operator
    bool operator!=(const RichPeak2D& rhs) const
    {
      return !(operator==(rhs));
    }

  };

} // namespace OpenMS

#endif // OPENMS_KERNEL_RICHPEAK2D_H
