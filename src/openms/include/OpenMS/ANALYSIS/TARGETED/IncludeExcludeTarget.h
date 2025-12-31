// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/CVTermList.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>
#include <OpenMS/CONCEPT/HashUtils.h>

#include <functional>
#include <vector>

namespace OpenMS
{
  /**
    @brief This class stores a SRM/MRM transition

    The default values for precursor and product m/z values are
    set to numeric_limits<double>::max(). Default values for
    precursor an product charge is set to numeric_limits<Int>::max().
  */
  class OPENMS_DLLAPI IncludeExcludeTarget :
    public CVTermList
  {

public:

    typedef TargetedExperimentHelper::Configuration Configuration;
    typedef TargetedExperimentHelper::RetentionTime RetentionTime;

    /** @name Constructors and destructors
    */
    //@{
    /// default constructor
    IncludeExcludeTarget();

    /// copy constructor
    IncludeExcludeTarget(const IncludeExcludeTarget & rhs);

    /// destructor
    ~IncludeExcludeTarget() override;
    //@}

    /// assignment operator
    IncludeExcludeTarget & operator=(const IncludeExcludeTarget & rhs);

    /** @name Accessors
    */
    //@{
    void setName(const String & name);

    const String & getName() const;

    void setPeptideRef(const String & peptide_ref);

    const String & getPeptideRef() const;

    void setCompoundRef(const String & compound_ref);

    const String & getCompoundRef() const;

    /// sets the precursor mz (Q1 value)
    void setPrecursorMZ(double mz);

    double getPrecursorMZ() const;

    void setPrecursorCVTermList(const CVTermList & list);

    void addPrecursorCVTerm(const CVTerm & cv_term);

    const CVTermList & getPrecursorCVTermList() const;

    void setProductMZ(double mz);

    double getProductMZ() const;

    void setProductCVTermList(const CVTermList & list);

    void addProductCVTerm(const CVTerm & cv_term);

    const CVTermList & getProductCVTermList() const;

    void setInterpretations(const std::vector<CVTermList> & interpretations);

    const std::vector<CVTermList> & getInterpretations() const;

    void addInterpretation(const CVTermList & interpretation);

    void setConfigurations(const std::vector<Configuration> & configuration);

    const std::vector<Configuration> & getConfigurations() const;

    void addConfiguration(const Configuration & configuration);

    void setPrediction(const CVTermList & prediction);

    void addPredictionTerm(const CVTerm & prediction);

    const CVTermList & getPrediction() const;

    void setRetentionTime(RetentionTime rt);

    const RetentionTime & getRetentionTime() const;
    //@}

    /** @name Predicates
    */
    //@{
    /// equality operator
    bool operator==(const IncludeExcludeTarget & rhs) const;

    /// inequality operator
    bool operator!=(const IncludeExcludeTarget & rhs) const;
    //@}

protected:

    void updateMembers_();

    String name_;

    double precursor_mz_;

    CVTermList precursor_cv_terms_;

    double product_mz_;

    CVTermList product_cv_terms_;

    std::vector<CVTermList> interpretation_list_;

    String peptide_ref_;

    String compound_ref_;

    std::vector<Configuration> configurations_;

    CVTermList prediction_;

    RetentionTime rts_;

  };
} // namespace OpenMS

namespace std
{
  /// std::hash specialization for OpenMS::IncludeExcludeTarget
  template <>
  struct hash<OpenMS::IncludeExcludeTarget>
  {
    std::size_t operator()(const OpenMS::IncludeExcludeTarget& t) const noexcept
    {
      std::size_t seed = 0;

      // Hash CVTermList base class (cv_terms_ map)
      // Hash by iterating through cv_terms and hashing accession + each CV term
      for (const auto& [accession, cv_terms] : t.getCVTerms())
      {
        OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(accession));
        for (const auto& cv_term : cv_terms)
        {
          OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cv_term.getAccession()));
          OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cv_term.getName()));
          OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cv_term.getValue().toString()));
        }
      }

      // Hash name_
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(t.getName()));

      // Hash precursor_mz_
      OpenMS::hash_combine(seed, OpenMS::hash_float(t.getPrecursorMZ()));

      // Hash precursor_cv_terms_
      for (const auto& [accession, cv_terms] : t.getPrecursorCVTermList().getCVTerms())
      {
        OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(accession));
        for (const auto& cv_term : cv_terms)
        {
          OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cv_term.getAccession()));
          OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cv_term.getName()));
          OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cv_term.getValue().toString()));
        }
      }

      // Hash product_mz_
      OpenMS::hash_combine(seed, OpenMS::hash_float(t.getProductMZ()));

      // Hash product_cv_terms_
      for (const auto& [accession, cv_terms] : t.getProductCVTermList().getCVTerms())
      {
        OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(accession));
        for (const auto& cv_term : cv_terms)
        {
          OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cv_term.getAccession()));
          OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cv_term.getName()));
          OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cv_term.getValue().toString()));
        }
      }

      // Hash interpretation_list_
      for (const auto& interp : t.getInterpretations())
      {
        for (const auto& [accession, cv_terms] : interp.getCVTerms())
        {
          OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(accession));
          for (const auto& cv_term : cv_terms)
          {
            OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cv_term.getAccession()));
            OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cv_term.getName()));
            OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cv_term.getValue().toString()));
          }
        }
      }

      // Hash peptide_ref_
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(t.getPeptideRef()));

      // Hash compound_ref_
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(t.getCompoundRef()));

      // Hash configurations_
      for (const auto& config : t.getConfigurations())
      {
        OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(config.contact_ref));
        OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(config.instrument_ref));
        // Hash CVTermList part of Configuration
        for (const auto& [accession, cv_terms] : config.getCVTerms())
        {
          OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(accession));
          for (const auto& cv_term : cv_terms)
          {
            OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cv_term.getAccession()));
            OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cv_term.getName()));
            OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cv_term.getValue().toString()));
          }
        }
        // Hash validations
        for (const auto& validation : config.validations)
        {
          for (const auto& [accession, cv_terms] : validation.getCVTerms())
          {
            OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(accession));
            for (const auto& cv_term : cv_terms)
            {
              OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cv_term.getAccession()));
              OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cv_term.getName()));
              OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cv_term.getValue().toString()));
            }
          }
        }
      }

      // Hash prediction_
      for (const auto& [accession, cv_terms] : t.getPrediction().getCVTerms())
      {
        OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(accession));
        for (const auto& cv_term : cv_terms)
        {
          OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cv_term.getAccession()));
          OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cv_term.getName()));
          OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cv_term.getValue().toString()));
        }
      }

      // Hash rts_ (RetentionTime)
      const auto& rt = t.getRetentionTime();
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(rt.software_ref));
      OpenMS::hash_combine(seed, OpenMS::hash_int(static_cast<int>(rt.retention_time_unit)));
      OpenMS::hash_combine(seed, OpenMS::hash_int(static_cast<int>(rt.retention_time_type)));
      OpenMS::hash_combine(seed, OpenMS::hash_int(rt.isRTset() ? 1 : 0));
      if (rt.isRTset())
      {
        OpenMS::hash_combine(seed, OpenMS::hash_float(rt.getRT()));
      }
      // Hash CVTermListInterface part of RetentionTime
      for (const auto& [accession, cv_terms] : rt.getCVTerms())
      {
        OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(accession));
        for (const auto& cv_term : cv_terms)
        {
          OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cv_term.getAccession()));
          OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cv_term.getName()));
          OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(cv_term.getValue().toString()));
        }
      }

      return seed;
    }
  };
} // namespace std

