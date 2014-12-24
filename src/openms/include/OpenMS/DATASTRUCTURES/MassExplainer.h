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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_MASSEXPLAINER_H
#define OPENMS_DATASTRUCTURES_MASSEXPLAINER_H

#include <OpenMS/DATASTRUCTURES/Adduct.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OpenMSConfig.h>

namespace OpenMS
{
  class Compomer;

  /**
    @brief computes empirical formulas for given mass differences using a set of allowed elements


    @ingroup Datastructures
  */
  class OPENMS_DLLAPI MassExplainer
  {

public:

    typedef Adduct::AdductsType AdductsType;     //vector<Adduct>
    typedef std::vector<Compomer>::const_iterator CompomerIterator;

    ///@name Constructors and destructor
    //@{
    /// Default constructor
    MassExplainer();

    /// Constructor
    MassExplainer(AdductsType adduct_base);

    /// Constructor
    MassExplainer(Int q_min, Int q_max, Int max_span, double thresh_logp);

    /// Constructor
    MassExplainer(AdductsType adduct_base, Int q_min, Int q_max, Int max_span, double thresh_logp, Size max_neutrals);


private:
    /// check consistency of input
    /// @param init_thresh_p set default threshold (set to "false" to keep current value)
    void init_(bool init_thresh_p);
public:
    /// Assignment operator
    MassExplainer & operator=(const MassExplainer & rhs);

    /// Destructor
    virtual ~MassExplainer();
    //@}


    /// fill map with possible mass-differences along with their explanation
    void compute();


    //@name Accessors
    //@{

    /// Sets the set of possible adducts
    void setAdductBase(AdductsType adduct_base);
    /// Returns the set of adducts
    AdductsType getAdductBase() const;

    /// return a compomer by its Id (useful after a query() ).
    const Compomer & getCompomerById(Size id) const;
    //@}


    /// search the mass database for explanations
    /// @param net_charge net charge of compomer seeked
    /// @param mass_to_explain mass in Da that needs explanation
    /// @param mass_delta allowed deviation from exact mass
    /// @param thresh_log_p  minimal log probability required
    /// @param firstExplanation begin range with candidates according to net_charge and mass
    /// @param lastExplanation  end range
    SignedSize query(const Int net_charge,
                     const float mass_to_explain,
                     const float mass_delta,
                     const float thresh_log_p,
                     std::vector<Compomer>::const_iterator & firstExplanation,
                     std::vector<Compomer>::const_iterator & lastExplanation) const;
protected:

    ///check if the generated compomer is valid judged by its probability, charges etc
    bool compomerValid_(const Compomer & cmp);

    /// create a proper adduct from formula and charge and probability
    Adduct createAdduct_(const String & formula, const Int charge, const double p) const;

    /// store possible explanations (as formula) for a certain ChargeDifference and MassDifference
    std::vector<Compomer> explanations_;
    /// all allowed adducts, whose combination explains the mass difference
    AdductsType adduct_base_;
    /// minimal expected charge
    Int q_min_;
    /// maximal expected charge
    Int q_max_;
    /// maximal span (in terms of charge) for co-features, e.g. a cluster with q={3,6} has span=4
    Int max_span_;
    /// minimum required probability of a compound (all other compounds are discarded)
    double thresh_p_;
    /// Maximum number of neutral(q=0) adducts
    Size max_neutrals_;

  };


} // namespace OpenMS

#endif // OPENMS_DATASTRUCTURES_MASSEXPLAINER_H
