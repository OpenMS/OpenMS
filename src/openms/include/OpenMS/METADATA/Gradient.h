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

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Representation of a HPLC gradient

    It consists of several eluents and timepoints.
    Linear behaviour between timepoints is assumed.

    @ingroup Metadata
  */
  class OPENMS_DLLAPI Gradient
  {
public:
    /// Constructor
    Gradient() = default;
    /// Copy constructor
    Gradient(const Gradient &) = default;
    /// Move constructor
    Gradient(Gradient&&) = default;
    /// Destructor
    ~Gradient();

    /// Assignment operator
    Gradient & operator=(const Gradient &) = default;
    /// Move assignment operator
    Gradient& operator=(Gradient&&) & = default;

    /// Equality operator
    bool operator==(const Gradient & source) const;
    /// Equality operator
    bool operator!=(const Gradient & source) const;

    /**
        @brief Adds an eluent at the end of the eluent array

        @exception Exception::InvalidValue is thrown if the same eluent name is used twice.
    */
    void addEluent(const String & eluent);
    /// removes all eluents
    void clearEluents();
    /// returns a const reference to the list of eluents
    const std::vector<String> & getEluents() const;

    /**
        @brief Adds a timepoint at the end of the timepoint array

        @exception Exception::OutOfRange is thrown if the new timepoint is before the last timepoint.
    */
    void addTimepoint(Int timepoint);
    /// removes all timepoints
    void clearTimepoints();
    /// returns a const reference to the list of timepoints
    const std::vector<Int> & getTimepoints() const;

    /**
        @brief sets the percentage of eluent @p eluent at timepoint @p timepoint

        @exception Exception::InvalidValue is thrown if the eluent, timepoint or percentage is invalid.
    */
    void setPercentage(const String & eluent, Int timepoint, UInt percentage);

    /**
        @brief returns a const reference to the percentages

        First dimension of the vector is the eluents, second dimension is the timepoints.
    */
    const std::vector<std::vector<UInt> > & getPercentages() const;

    /**
        @brief returns the percentage of an @p eluent at a @p timepoint

        @exception Exception::InvalidValue is thrown if the eluent or timepoint is invalid.
    */
    UInt getPercentage(const String & eluent, Int timepoint) const;

    /// sets all percentage values to 0
    void clearPercentages();

    /// checks if the percentages of all timepoints add up to 100%
    bool isValid() const;

protected:
    std::vector<String> eluents_;
    std::vector<Int> times_;
    // first dimension is eluents, second is times
    std::vector<std::vector<UInt> > percentages_;
  };

} // namespace OpenMS

