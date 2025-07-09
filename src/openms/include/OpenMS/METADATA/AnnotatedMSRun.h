// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: David Voigt, Timo Sachsenberg $
// -------------------------------------------------------------------------------------------------------------------------------------

#pragma once

#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <boost/range/combine.hpp>

#include <vector>

namespace OpenMS
{
  class PeptideIdentification;

  class MSSpectrum;

  /**
   * @brief Class for storing MS run data with peptide and protein identifications
   *
   * This class stores an MSExperiment (containing spectra) along with peptide and protein
   * identifications. Each spectrum in the MSExperiment is associated with a single
   * PeptideIdentification object.
   *
   * The class provides methods to access and modify these identifications, as well as
   * iterators to traverse the spectra and their associated identifications together.
   */
  class OPENMS_DLLAPI AnnotatedMSRun
  {
  public:
    using SpectrumIdRef = std::pair<MSSpectrum&, PeptideIdentification&>;
    using ConstSpectrumIdRef = std::pair<const MSSpectrum&, const PeptideIdentification&>;
    using SpectrumType = MSExperiment::SpectrumType;
    using ChromatogramType = MSExperiment::ChromatogramType;


    /// Default constructor
    AnnotatedMSRun() = default;

    /**
     * @brief Move constructor for efficiently loading a MSExperiment without a deep copy
     * @param experiment The MSExperiment to move into this object
     */
    explicit AnnotatedMSRun(MSExperiment&& experiment) : data(std::move(experiment))
    {};

    /// Move constructor
    AnnotatedMSRun(AnnotatedMSRun&&) = default;

    /// Copy constructor
    AnnotatedMSRun(const AnnotatedMSRun&) = default;
    AnnotatedMSRun& operator=(const AnnotatedMSRun&) = default;
    AnnotatedMSRun& operator=(AnnotatedMSRun&&) = default;

    /// Destructor
    ~AnnotatedMSRun() = default;

    /**
     * @brief Get the protein identification
     * @return A reference to the protein identification
     */
    std::vector<ProteinIdentification>& getProteinIdentifications()
    {
      return protein_ids_;
    }

    /**
     * @brief Get the protein identification (const version)
     * @return A const reference to the protein identification
     */
    const std::vector<ProteinIdentification>& getProteinIdentifications() const
    {
      return protein_ids_;
    }

    /**
     * @brief set the protein identifications
     * @param ids Vector of protein identifications
     */
    void setProteinIdentifications(const std::vector<ProteinIdentification>& ids)
    {
      protein_ids_ = ids;
    }

    /**
     * @brief Set the protein identifications (move version)
     * @param ids Vector of protein identifications
     */
    void setProteinIdentifications(std::vector<ProteinIdentification>&& ids)
    {
      protein_ids_ = std::move(ids);
    }

    /**
     * @brief Get all peptide identifications for all spectra
     * @return A reference to the vector of peptide identifications
     */
    PeptideIdentificationList& getPeptideIdentifications();
    
    /**
     * @brief Get all peptide identifications for all spectra (const version)
     * @return A const reference to the vector of peptide identifications
     */
    const PeptideIdentificationList& getPeptideIdentifications() const;

    /**
     * @brief Set all peptide identifications for all spectra (move version)
     * @param ids Vector of peptide identifications
     */
    void setPeptideIdentifications(PeptideIdentificationList&& ids);

    /**
     * @brief Set all peptide identifications for all spectra
     * @param ids Vector of peptide identifications
     */
    void setPeptideIdentifications(const PeptideIdentificationList& ids);

    /**
     * @brief Get the MSExperiment
     * @return A reference to the MSExperiment
     */
    MSExperiment& getMSExperiment();
    
    /**
     * @brief Get the MSExperiment (const version)
     * @return A const reference to the MSExperiment
     */
    const MSExperiment& getMSExperiment() const;

    /**
     * @brief Set the MSExperiment
     * @param experiment The MSExperiment to set
     */
    void setMSExperiment(MSExperiment&& experiment);

    /**
     * @brief Set the MSExperiment
     * @param experiment The MSExperiment to set
     */
    void setMSExperiment(const MSExperiment& experiment);

    /**
     * @brief Get a const iterator to the beginning of the data
     * @return A const iterator to the beginning
     */
    inline auto cbegin() const
    {
      checkPeptideIdSize_(OPENMS_PRETTY_FUNCTION);
      return PairIterator(data.getSpectra().cbegin(), peptide_ids_.cbegin());
    }

    /**
     * @brief Get an iterator to the beginning of the data
     * @return An iterator to the beginning
     */
    inline auto begin()
    {
      checkPeptideIdSize_(OPENMS_PRETTY_FUNCTION);
      return PairIterator(data.getSpectra().begin(), peptide_ids_.begin());
    }

    /**
     * @brief Get a const iterator to the beginning of the data
     * @return A const iterator to the beginning
     */
    inline auto begin() const
    {
      checkPeptideIdSize_(OPENMS_PRETTY_FUNCTION);
      return PairIterator(data.getSpectra().cbegin(), peptide_ids_.cbegin());
    }

    /**
     * @brief Get an iterator to the end of the data
     * @return An iterator to the end
     */
    inline auto end()
    {
      return PairIterator(data.getSpectra().end(), peptide_ids_.end());
    }

    /**
     * @brief Get a const iterator to the end of the data
     * @return A const iterator to the end
     */
    inline auto end() const
    {
      return PairIterator(data.getSpectra().end(), peptide_ids_.end());
    }

    /**
     * @brief Get a const iterator to the end of the data
     * @return A const iterator to the end
     */
    inline auto cend() const
    {
      return PairIterator(data.getSpectra().cend(), peptide_ids_.cend());
    }

    /**
     * @brief Access a spectrum and its associated peptide identification
     * @param idx The index of the spectrum
     * @return A pair of references to the spectrum and its peptide identification
     */
    inline SpectrumIdRef operator[](size_t idx)
    {
      if (idx >= peptide_ids_.size())
      {
        throw Exception::IndexOverflow(__FILE__, __LINE__,
          OPENMS_PRETTY_FUNCTION,
          idx, peptide_ids_.size());
      }
      if (idx >= data.getSpectra().size())
      {
        throw Exception::IndexOverflow(__FILE__, __LINE__,
          OPENMS_PRETTY_FUNCTION,
          idx, data.getSpectra().size());
      }      
      return {data.getSpectra()[idx], peptide_ids_[idx]};
    }

    /**
     * @brief Access a spectrum and its associated peptide identification (const version)
     * @param idx The index of the spectrum
     * @return A pair of const references to the spectrum and its peptide identification
     */
    inline ConstSpectrumIdRef operator[](size_t idx) const
    {
      if (idx >= peptide_ids_.size())
      {
        throw Exception::IndexOverflow(__FILE__, __LINE__,
          OPENMS_PRETTY_FUNCTION,
          idx, peptide_ids_.size());
      }
      if (idx >= data.getSpectra().size())
      {
        throw Exception::IndexOverflow(__FILE__, __LINE__,
          OPENMS_PRETTY_FUNCTION,
          idx, data.getSpectra().size());
      }
      return {data.getSpectra()[idx], peptide_ids_[idx]};
    }

    /**
     * @brief Iterator for pairs of spectra and peptide identifications
     *
     * This iterator allows traversing the spectra and their associated peptide
     * identifications together.
     */
    template<typename T1, typename T2>
    struct PairIterator
    {     
      using iterator_category = std::forward_iterator_tag;
      using difference_type = std::ptrdiff_t;

      /**
       * @brief Constructor
       * @param ptr1 Iterator to the spectra
       * @param ptr2 Iterator to the peptide identifications
       */
      PairIterator(T1 ptr1, T2 ptr2) : m_ptr1(ptr1), m_ptr2(ptr2)
      {}

      /**
       * @brief Pre-increment operator
       * @return Reference to this iterator after incrementing
       */
      PairIterator& operator++()
      {
        ++m_ptr1;
        ++m_ptr2;
        return *this;
      }

      /**
       * @brief Post-increment operator
       * @return Copy of this iterator before incrementing
       */
      PairIterator operator++(int)
      {
        auto tmp(*this);
        ++(*this);
        return tmp;
      }

      /**
       * @brief Dereference operator
       * @return A pair of references to the current spectrum and peptide identification
       */
      auto operator*()
      {
        return std::make_pair(std::ref(*m_ptr1), std::ref(*m_ptr2));
      }

      /**
       * @brief Equality operator
       * @param a First iterator
       * @param b Second iterator
       * @return True if the iterators are equal
       */
      inline friend bool operator==(const PairIterator& a, const PairIterator& b)
      {
        return a.m_ptr1 == b.m_ptr1 && a.m_ptr2 == b.m_ptr2;
      }

      /**
       * @brief Inequality operator
       * @param a First iterator
       * @param b Second iterator
       * @return True if the iterators are not equal
       */
      inline friend bool operator!=(const PairIterator& a, const PairIterator& b)
      {
        return !(a == b);
      }

    private:
      T1 m_ptr1;
      T2 m_ptr2;
    };

    typedef AnnotatedMSRun::PairIterator<std::vector<MSSpectrum>::iterator, PeptideIdentificationList::iterator> Iterator;
    typedef AnnotatedMSRun::PairIterator<std::vector<MSSpectrum>::const_iterator, PeptideIdentificationList::const_iterator> ConstIterator;

  private:

    // Helper to enforce invariant
    void checkPeptideIdSize_(const char* function_name) const;

    PeptideIdentificationList peptide_ids_;
    std::vector<ProteinIdentification> protein_ids_;
    MSExperiment data;
  };
}
