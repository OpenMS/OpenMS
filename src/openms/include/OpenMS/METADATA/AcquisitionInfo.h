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

#include <OpenMS/METADATA/Acquisition.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Description of the combination of raw data to a single spectrum

    Specification for combining raw scans ( Acquisition ) into a single spectrum.
    A list of acquisitions from the original raw file can be specified.

    @ingroup Metadata
  */
  class OPENMS_DLLAPI AcquisitionInfo :
    private std::vector<Acquisition>,
    public MetaInfoInterface
  {
private:
    typedef std::vector<Acquisition> ContainerType;

public:
    /// Constructor
    AcquisitionInfo() = default;
    /// Copy constructor
    AcquisitionInfo(const AcquisitionInfo&) = default;
    /// Move constructor
    AcquisitionInfo(AcquisitionInfo&&) = default;
    /// Destructor
    ~AcquisitionInfo() = default;

    /// Assignment operator
    AcquisitionInfo& operator=(const AcquisitionInfo&) = default;
    /// Move assignment operator
    AcquisitionInfo& operator=(AcquisitionInfo&&) & = default;

    /// Equality operator
    bool operator==(const AcquisitionInfo& rhs) const;
    /// Equality operator
    bool operator!=(const AcquisitionInfo& rhs) const;

    /// returns the method of combination
    const String& getMethodOfCombination() const;
    /// sets the method of combination
    void setMethodOfCombination(const String& method_of_combination);

    ///@name Export methods from private base std::vector<Acquisition>
    //@{
    using ContainerType::operator[];
    using ContainerType::begin;
    using ContainerType::end;
    using ContainerType::size;
    using ContainerType::push_back;
    using ContainerType::empty;
    using ContainerType::back;
    using ContainerType::insert;
    using ContainerType::resize;

    using ContainerType::iterator;
    using ContainerType::const_iterator;
    //@}

protected:
    String method_of_combination_;

  };
} // namespace OpenMS

