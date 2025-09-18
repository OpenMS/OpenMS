// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/Adduct.h>

#include <OpenMS/CONCEPT/Types.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OpenMSConfig.h>

namespace OpenMS
{
  class Compomer;

  /**
    @brief Computes empirical formulas for given mass differences using a set of allowed elements
    
    MassExplainer is used to explain observed mass differences between features by
    determining the most likely combination of adducts that could cause such differences.
    
    The class works by:
    1. Taking a set of allowed adducts (elements, molecules, or ions)
    2. Computing all possible combinations of these adducts that could explain observed mass differences
    3. Providing a query interface to search for explanations for specific mass differences
    
    This is particularly useful in mass spectrometry data analysis for identifying
    related features that represent the same analyte but with different adducts or charge states.
    
    @ingroup Datastructures
  */
  class OPENMS_DLLAPI MassExplainer
  {

public:

    /// Type definition for a vector of Adduct objects
    typedef Adduct::AdductsType AdductsType;
    
    /// Type definition for an iterator over Compomer objects
    typedef std::vector<Compomer>::const_iterator CompomerIterator;

    ///@name Constructors and destructor
    //@{
    /**
      @brief Default constructor
      
      Initializes with default parameters:
      - No adducts
      - Charge range from -2 to +4
      - Maximum charge span of 4
      - Log probability threshold of -5.0
      - Maximum number of neutral adducts of 2
    */
    MassExplainer();

    /**
      @brief Constructor with custom adduct base
      
      @param adduct_base Set of allowed adducts to use for mass difference explanations
    */
    MassExplainer(AdductsType adduct_base);

    /**
      @brief Constructor with custom charge parameters
      
      @param q_min Minimum charge state to consider
      @param q_max Maximum charge state to consider
      @param max_span Maximum allowed charge span between related features
      @param thresh_logp Minimum log probability threshold for accepting explanations
    */
    MassExplainer(Int q_min, Int q_max, Int max_span, double thresh_logp);

    /**
      @brief Constructor with all custom parameters
      
      @param adduct_base Set of allowed adducts to use for mass difference explanations
      @param q_min Minimum charge state to consider
      @param q_max Maximum charge state to consider
      @param max_span Maximum allowed charge span between related features
      @param thresh_logp Minimum log probability threshold for accepting explanations
      @param max_neutrals Maximum number of neutral adducts allowed in an explanation
    */
    MassExplainer(AdductsType adduct_base, Int q_min, Int q_max, Int max_span, double thresh_logp, Size max_neutrals);


private:
    /**
      @brief Check consistency of input parameters and initialize internal data structures
      
      This method validates the input parameters and sets default values where needed.
      
      @param init_thresh_p Whether to initialize the probability threshold with default value
                          (set to "false" to keep current value)
    */
    void init_(bool init_thresh_p);
public:
    /**
      @brief Assignment operator
      
      @param rhs Source object to assign from
      @return Reference to this object
    */
    MassExplainer& operator=(const MassExplainer& rhs);

    /// Destructor
    virtual ~MassExplainer();
    //@}


    /**
      @brief Compute all possible mass differences and their explanations
      
      This method generates all possible combinations of adducts from the adduct base
      and stores them internally for later querying. This must be called after
      changing any parameters and before performing queries.
    */
    void compute();


    //@name Accessors
    //@{

    /**
      @brief Set the base set of allowed adducts
      
      @param adduct_base Vector of adducts to use for explanations
    */
    void setAdductBase(AdductsType adduct_base);
    
    /**
      @brief Get the current set of allowed adducts
      
      @return Vector of adducts currently used for explanations
    */
    AdductsType getAdductBase() const;

    /**
      @brief Get a specific compomer by its ID
      
      This is typically used after a query() to retrieve detailed information
      about a specific explanation.
      
      @param id ID of the compomer to retrieve
      @return Reference to the requested compomer
    */
    const Compomer& getCompomerById(Size id) const;
    //@}


    /**
      @brief Search for explanations of a given mass difference
      
      This method searches the precomputed explanations for those that match
      the given mass difference within the specified tolerance and have the
      required net charge.
      
      @param net_charge Net charge of the compomer being sought
      @param mass_to_explain Mass difference in Da that needs explanation
      @param mass_delta Allowed deviation from exact mass (tolerance)
      @param thresh_log_p Minimum log probability required for explanations
      @param firstExplanation Output iterator to the beginning of matching explanations
      @param lastExplanation Output iterator to the end of matching explanations
      @return Number of explanations found, or negative value if no explanations found
    */
    SignedSize query(const Int net_charge,
                     const float mass_to_explain,
                     const float mass_delta,
                     const float thresh_log_p,
                     std::vector<Compomer>::const_iterator& firstExplanation,
                     std::vector<Compomer>::const_iterator& lastExplanation) const;
protected:

    /**
      @brief Check if a generated compomer is valid based on its probability, charges, etc.
      
      @param cmp The compomer to validate
      @return True if the compomer is valid, false otherwise
    */
    bool compomerValid_(const Compomer& cmp) const;

    /**
      @brief Create a proper adduct from formula, charge, and probability
      
      @param formula Chemical formula of the adduct
      @param charge Charge of the adduct
      @param p Probability of the adduct
      @return Adduct object with the specified properties
    */
    Adduct createAdduct_(const String& formula, const Int charge, const double p) const;

    /// Vector storing all possible explanations for mass differences
    std::vector<Compomer> explanations_;
    
    /// Set of allowed adducts that can be combined to explain mass differences
    AdductsType adduct_base_;
    
    /// Minimum charge state to consider in explanations
    Int q_min_;
    
    /// Maximum charge state to consider in explanations
    Int q_max_;
    
    /// Maximum allowed charge span between related features (e.g., a cluster with q={3,6} has span=4)
    Int max_span_;
    
    /// Minimum required probability threshold for accepting explanations
    double thresh_p_;
    
    /// Maximum number of neutral (q=0) adducts allowed in an explanation
    Size max_neutrals_;

  };


} // namespace OpenMS

