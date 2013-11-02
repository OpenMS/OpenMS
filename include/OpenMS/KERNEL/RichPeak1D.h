// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

#ifndef OPENMS_KERNEL_RICHPEAK1D_H
#define OPENMS_KERNEL_RICHPEAK1D_H

#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS
{

  /**
    @brief A 1-dimensional raw data point or peak with meta information.

    This data structure is intended for continuous data or peak data.
    If you do not need to annotated single peaks with meta data, use Peak1D instead.

    @ingroup Kernel
  */
  class OPENMS_DLLAPI RichPeak1D :
    public Peak1D,
    public MetaInfoInterface
  {
public:

    /// Default constructor
    inline RichPeak1D() :
      Peak1D(),
      MetaInfoInterface()
    {}

    /// Copy constructor
    inline RichPeak1D(const RichPeak1D & p) :
      Peak1D(p),
      MetaInfoInterface(p)
    {}

    /// Destructor
    ~RichPeak1D()
    {}

    /// Assignment operator
    inline RichPeak1D & operator=(const RichPeak1D & rhs)
    {
      if (this == &rhs) return *this;

      Peak1D::operator=(rhs);
      MetaInfoInterface::operator=(rhs);

      return *this;
    }

    /// Equality operator
    inline bool operator==(const RichPeak1D & rhs) const
    {
      return Peak1D::operator==(rhs) &&
             MetaInfoInterface::operator==(rhs);
    }

    /// Equality operator
    inline bool operator!=(const RichPeak1D & rhs) const
    {
      return !(operator==(rhs));
    }

  };

} // namespace OpenMS

#endif // OPENMS_KERNEL_RICHPEAK1D_H
