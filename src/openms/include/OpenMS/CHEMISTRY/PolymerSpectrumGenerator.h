// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: OpenMS Team $
// $Authors: OpenMS Team $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/SequenceBase.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>

namespace OpenMS
{
  /**
   * @brief CRTP base class for spectrum generators working with polymer sequences
   * 
   * This template base class uses the Curiously Recurring Template Pattern (CRTP)
   * to provide a common interface for spectrum generators that work with different
   * polymer sequence types (AASequence, NASequence, etc.).
   * 
   * @tparam Derived The derived spectrum generator class
   * 
   * @ingroup Chemistry
   */
  template<typename Derived>
  class OPENMS_DLLAPI PolymerSpectrumGenerator :
    public DefaultParamHandler
  {
  public:
    
    /// Default constructor
    PolymerSpectrumGenerator(const String& name) : DefaultParamHandler(name) {}
    
    /// Destructor
    virtual ~PolymerSpectrumGenerator() = default;

    /// Get the derived instance (CRTP pattern)
    const Derived& derived() const noexcept
    {
      return static_cast<const Derived&>(*this);
    }

    /// Get the derived instance (CRTP pattern)
    Derived& derived() noexcept
    {
      return static_cast<Derived&>(*this);
    }

    /**
     * @brief Generic spectrum generation for any SequenceBase-derived type
     * 
     * @tparam SeqType Either AASequence or NASequence
     * @param spectrum Output spectrum
     * @param sequence Input polymer sequence  
     * @param min_charge Minimum charge state
     * @param max_charge Maximum charge state
     */
    template<typename SeqType>
    void generateSpectrum(MSSpectrum& spectrum, const SequenceBase<SeqType>& sequence, 
                         Int min_charge, Int max_charge) const
    {
      static_assert(is_sequence_type_v<SequenceBase<SeqType>>, "SeqType must be a sequence type");
      
      // Delegate to derived class implementation using CRTP
      derived().generateSpectrumImpl(spectrum, sequence.derived(), min_charge, max_charge);
    }

    /**
     * @brief Get the sequence type that this generator works with
     * 
     * This method should be implemented by derived classes to specify
     * which sequence type they work with.
     */
    virtual String getSequenceType() const = 0;

    /**
     * @brief Check if this generator supports the given sequence type
     * 
     * @tparam SeqType The sequence type to check
     * @return True if supported, false otherwise
     */
    template<typename SeqType>
    bool supportsSequenceType() const
    {
      if constexpr (std::is_same_v<SeqType, AASequence>)
      {
        return getSequenceType() == "AASequence" || getSequenceType() == "Peptide";
      }
      else if constexpr (std::is_same_v<SeqType, NASequence>)
      {
        return getSequenceType() == "NASequence" || getSequenceType() == "Oligonucleotide";
      }
      return false;
    }

  protected:
    /// Protected constructor to prevent direct instantiation
    PolymerSpectrumGenerator() = default;
  };

} // namespace OpenMS