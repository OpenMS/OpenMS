// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_CHARGEPAIR_H
#define OPENMS_DATASTRUCTURES_CHARGEPAIR_H

#include <OpenMS/DATASTRUCTURES/Compomer.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OpenMSConfig.h>

#include <iosfwd>
#include <vector>

namespace OpenMS
{

  /**
    @brief Representation of a (putative) link between two Features, which stem from the same compound
      but have different charge (including different adduct ions (H+, Na+, ..)

    A ChargePair represents an edge between two Features and specifies their respective charge and adducts,
    so that when decharged they can be explained as stemming from the same compound.


    @ingroup Datastructures
  */
  class OPENMS_DLLAPI ChargePair
  {

public:
    ///@name Constructors and destructor
    //@{
    /// Default constructor
    ChargePair();

    /// Constructor from map index, element index and Feature
    ChargePair(const Size & index0,
               const Size & index1,
               const Int & charge0,
               const Int & charge1,
               const Compomer & compomer,
               const double & mass_diff,
               const bool active);

    /// Copy constructor
    ChargePair(const ChargePair & rhs);

    /// Assignment operator
    ChargePair & operator=(const ChargePair & rhs);

    /// Destructor
    virtual ~ChargePair();

    //@}

    //@name Accessors
    //@{
    /// Returns the charge (for element 0 or 1)
    Int getCharge(UInt pairID) const;

    /// Set the charge (for element 0 or 1)
    void setCharge(UInt pairID, Int e);

    /// Returns the element index (for element 0 or 1)
    Size getElementIndex(UInt pairID) const;

    /// Set the element index (for element 0 or 1)
    void setElementIndex(UInt pairID, Size e);

    /// Returns the Id of the compomer that explains the mass difference
    const Compomer & getCompomer() const;

    /// Set the compomer id
    void setCompomer(const Compomer & compomer);

    /// Returns the mass difference
    double getMassDiff() const;

    /// Sets the mass difference
    void setMassDiff(double mass_diff);

    /// Returns the ILP edge score
    double getEdgeScore() const;

    /// Sets the ILP edge score
    void setEdgeScore(double score);

    /// is this pair realized?
    bool isActive() const;

    void setActive(const bool active);

    //@}

    /// Equality operator
    virtual bool operator==(const ChargePair & i) const;

    /// Equality operator
    virtual bool operator!=(const ChargePair & i) const;

protected:

    /// Int of the first element within the FeatureMap
    Size feature0_index_;
    /// Int of the second element within the FeatureMap
    Size feature1_index_;
    /// Assumed charge of the first feature
    Int feature0_charge_;
    /// Assumed charge of the second feature
    Int feature1_charge_;
    /// Compomer that explains the mass difference
    Compomer compomer_;
    /// mass difference (after explanation by compomer)
    double mass_diff_;
    /// Score of this edge used in ILP
    double score_;
    /// was this pair realized by ILP?
    bool is_active_;
  };

  ///Print the contents of a ChargePair to a stream.
  OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const ChargePair & cons);

} // namespace OpenMS

#endif // OPENMS_DATASTRUCTURES_CHARGEPAIR_H
