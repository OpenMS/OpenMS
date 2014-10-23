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

#ifndef OPENMS_DATASTRUCTURES_ADDUCT_H
#define OPENMS_DATASTRUCTURES_ADDUCT_H

#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{

  class OPENMS_DLLAPI Adduct
  {
public:

    typedef std::vector<Adduct> AdductsType;

    /// Default C'tor
    Adduct();

    /// C'tor with initial charge
    Adduct(Int charge);

    /// C'tor for all members
    Adduct(Int charge, Int amount, double singleMass, String formula, double log_prob, double rt_shift, const String label = "");

    /// Increase amount of this adduct by factor @param m
    Adduct operator*(const Int m) const;
    /// Add two adducts amount if they are equal (defined by equal formula)
    Adduct operator+(const Adduct & rhs);
    /// Add other adducts amount to *this (equal formula required!)
    void operator+=(const Adduct & rhs);


    /// Print the contents of an Adduct to a stream.
    friend OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const Adduct & a);

    /// Comparator
    friend OPENMS_DLLAPI bool operator==(const Adduct & a, const Adduct & b);

    //@{ Accessors
    const Int & getCharge() const;

    void setCharge(const Int & charge);

    const Int & getAmount() const;
    void setAmount(const Int & amount);

    const double & getSingleMass() const;
    void setSingleMass(const double & singleMass);

    const double & getLogProb() const;
    void setLogProb(const double & log_prob);

    const String & getFormula() const;
    void setFormula(const String & formula);

    const double & getRTShift() const;
    const String & getLabel() const;
    //}

private:
    Int charge_; //< usually +1
    Int amount_; //< number of entities
    double singleMass_; //< mass of a single entity
    double log_prob_;   //< log probability of observing a single entity of this adduct
    String formula_;   //< chemical formula (parsable by EmpiricalFormula)
    double rt_shift_;     //< RT shift induced by a single entity of this adduct (this is for adducts attached prior to ESI, e.g. labeling)
    String label_;     //< Label for this adduct (can be used to indicate heavy labels)

    String checkFormula_(const String & formula);

  };

} // namespace OpenMS


#endif
