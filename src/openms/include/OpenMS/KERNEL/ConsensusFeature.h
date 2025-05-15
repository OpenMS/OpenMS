// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/KERNEL/FeatureHandle.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OpenMSConfig.h>

#include <set>

namespace OpenMS
{
  class FeatureMap;
  class Peak2D;

  /**
    @brief A consensus feature spanning multiple LC-MS/MS experiments.

    A ConsensusFeature represents analytes that have been
    quantified across multiple LC-MS/MS experiments. Each analyte in a
    ConsensusFeature is linked to its original LC-MS/MS run through a
    unique identifier.

    A consensus feature represents corresponding features in multiple feature
    maps. The corresponding features are represented a set of @ref
    FeatureHandle instances.  Each ConsensusFeature "contains" zero or more
    FeatureHandles.

    @see ConsensusMap

    @ingroup Kernel
  */
  class OPENMS_DLLAPI ConsensusFeature :
    public BaseFeature
  {
public:
    ///Type definitions
    //@{
    typedef std::set<FeatureHandle, FeatureHandle::IndexLess> HandleSetType;
    typedef HandleSetType::const_iterator const_iterator;
    typedef HandleSetType::iterator iterator;
    typedef HandleSetType::const_reverse_iterator const_reverse_iterator;
    typedef HandleSetType::reverse_iterator reverse_iterator;
    //@}

    /// Compare by size(), the number of consensus elements
    struct SizeLess
    {
      inline bool operator()(ConsensusFeature const& left, ConsensusFeature const& right) const
      {
        return left.size() < right.size();
      }

      inline bool operator()(ConsensusFeature const& left, UInt64 const& right) const
      {
        return left.size() < right;
      }

      inline bool operator()(UInt64 const& left, ConsensusFeature const& right) const
      {
        return left < right.size();
      }

      inline bool operator()(const UInt64& left, const UInt64& right) const
      {
        return left < right;
      }

    };

    /// Compare by the sets of consensus elements (lexicographically)
    struct MapsLess
    {
      inline bool operator()(ConsensusFeature const& left, ConsensusFeature const& right) const
      {
        return std::lexicographical_compare(left.begin(), left.end(), right.begin(), right.end(), FeatureHandle::IndexLess());
      }

    };

    /// slim struct to feed the need for systematically storing of ratios.
    struct Ratio
    {
      Ratio()
      {
      }

      Ratio(const Ratio& rhs)
      {
        ratio_value_ = rhs.ratio_value_;
        denominator_ref_ = rhs.denominator_ref_;
        numerator_ref_ = rhs.numerator_ref_;
        description_ = rhs.description_;
      }

      virtual ~Ratio()
      {
      }

      Ratio& operator=(const Ratio& rhs)
      {
        if (&rhs != this)
        {
          ratio_value_ = rhs.ratio_value_;
          denominator_ref_ = rhs.denominator_ref_;
          numerator_ref_ = rhs.numerator_ref_;
          description_ = rhs.description_;
        }
        return *this;
      }

      // @TODO: members are public, names shouldn't end in underscores
      double ratio_value_;
      String denominator_ref_;
      String numerator_ref_;
      std::vector<String> description_;
      //TODO ratio cv info
    };

    ///@name Constructors and Destructor
    //@{
    /// Default constructor
    ConsensusFeature();

    /// Copy constructor
    ConsensusFeature(const ConsensusFeature& rhs) = default;

    /// Move constructor
    ConsensusFeature(ConsensusFeature&& rhs) = default;

    /// Constructor from basic feature
    explicit ConsensusFeature(const BaseFeature& feature);

    /**
      @brief Constructor with map and element index for a singleton consensus
      feature.

      Sets the consensus feature position and intensity to the values of @p element as well.
    */
    ConsensusFeature(UInt64 map_index, const Peak2D& element, UInt64 element_index);

    /**
      @brief Constructor with map index for a singleton consensus feature.

      Sets the consensus feature position, intensity, charge, quality, and peptide identifications
      to the values of @p element as well.
    */
    ConsensusFeature(UInt64 map_index, const BaseFeature& element);

    /// Assignment operator
    ConsensusFeature& operator=(const ConsensusFeature& rhs) = default;

    /// Move Assignment operator
    ConsensusFeature& operator=(ConsensusFeature&& rhs) = default;

    /// Destructor
    ~ConsensusFeature() override;
    //@}


    ///@name Management of feature handles
    //@{

    /**
      @brief Adds all feature handles (of the CF) into the consensus feature
    */
    void insert(const ConsensusFeature& cf);
    void insert(ConsensusFeature&& cf);

    /**
      @brief Adds an feature handle into the consensus feature

      @exception Exception::InvalidValue is thrown if a handle with the same map index and unique
      id already exists.
    */
    void insert(const FeatureHandle& handle);
    void insert(FeatureHandle&& handle);

    /// Adds all feature handles in @p handle_set to this consensus feature.
    void insert(const HandleSetType& handle_set);
    void insert(HandleSetType&& handle_set);

    /**
      @brief Creates a FeatureHandle and adds it

      @exception Exception::InvalidValue is thrown if a handle with the same map index and unique
      id already exists.
    */
    void insert(UInt64 map_index, const Peak2D& element, UInt64 element_index);

    /**
      @brief Creates a FeatureHandle and adds it

      @exception Exception::InvalidValue is thrown if a handle with the same map index and unique
      id already exists.
    */
    void insert(UInt64 map_index, const BaseFeature& element);

    /// Non-mutable access to the contained feature handles
    const HandleSetType& getFeatures() const;

    /// Mutable access to a copy of the contained feature handles
    std::vector<FeatureHandle> getFeatureList() const;

    /// Set the feature set to a new one
    void setFeatures(HandleSetType h);
    //@}

    ///@name Accessors
    //@{
    /// Returns the position range of the contained elements
    DRange<2> getPositionRange() const;
    /// Returns the intensity range of the contained elements
    DRange<1> getIntensityRange() const;

    //@}

    /**
      @brief Computes and updates the consensus position, intensity, and charge.

      The position and intensity of the contained feature handles is averaged.
      The most frequent charge state wins, while the tie breaking prefers
      smaller (absolute) charges.

      @note This method has to be called explicitly, <i>after</i> adding the feature handles.
    */
    void computeConsensus();

    /**
      @brief Computes and updates the consensus position, intensity, and charge.

      The m/z position is the lowest m/z value of the feature handles. The RT position and intensity
      of the contained feature handles is averaged. The most frequent charge state wins, while the
      tie breaking prefers smaller (absolute) charges.

      @note This method has to be called explicitly, <i>after</i> adding the feature handles.
    */
    void computeMonoisotopicConsensus();

    /**
      @brief Computes the uncharged parent RT & mass, assuming the handles are charge variants.

      The position of the feature handles (decharged) is averaged (using intensity as weights if
      @p intensity_weighted_averaging is true). Intensities are summed up. Charge is set to 0.
      Mass calculation: If the given features contain a metavalue "dc_charge_adduct_mass" then this
      will be used as adduct mass instead of weight(H+) * charge.

      @note This method has to be called explicitly, <i>after</i> adding the feature handles.

      @param fm Input feature map, which provides additional information on the features
      @param intensity_weighted_averaging Use unweighted averaging (default) or weighted by intensity
    */
    void computeDechargeConsensus(const FeatureMap& fm, bool intensity_weighted_averaging = false);

    /**
      @brief Add a ratio.

      Connects a ratio to the ConsensusFeature.

      @note still experimental. consensusfeaturehandler will ignore it.
    */
    void addRatio(const Ratio& r);

    /**
      @brief Add a ratio vector.

      Connects the ratios to the ConsensusFeature.

      @note still experimental. consensusfeaturehandler will ignore it.
    */
    void setRatios(std::vector<Ratio>& rs);

    /**
      @brief Get the ratio vector.
    */
    const std::vector<Ratio>& getRatios() const;

    /**
      @brief Get the ratio vector.
    */
    std::vector<Ratio>& getRatios();

    ///@name Accessors for set of FeatureHandles
    //@{
    Size size() const;

    const_iterator begin() const;

    iterator begin();

    const_iterator end() const;

    iterator end();

    const_reverse_iterator rbegin() const;

    reverse_iterator rbegin();

    const_reverse_iterator rend() const;

    reverse_iterator rend();

    void clear();

    bool empty() const;
    //@}

private:

    HandleSetType handles_;
    std::vector<Ratio> ratios_;


  };

  ///Print the contents of a ConsensusFeature to a stream
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const ConsensusFeature& cons);

} // namespace OpenMS
