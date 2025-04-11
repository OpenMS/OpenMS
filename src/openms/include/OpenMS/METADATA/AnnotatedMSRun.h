// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: David Voigt $
// -------------------------------------------------------------------------------------------------------------------------------------

#pragma once

#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/KERNEL/MSExperiment.h>
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
    typedef std::pair<MSSpectrum&, PeptideIdentification&> Mapping;
    typedef std::pair<const MSSpectrum&, const PeptideIdentification&> ConstMapping;

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
     * @brief Get all peptide identifications for all spectra
     * @return A reference to the vector of peptide identifications
     */
    std::vector<PeptideIdentification>& getPeptideIdentifications();
    
    /**
     * @brief Get all peptide identifications for all spectra (const version)
     * @return A const reference to the vector of peptide identifications
     */
    const std::vector<PeptideIdentification>& getPeptideIdentifications() const;

    /**
     * @brief Set all peptide identifications for all spectra
     * @param ids Vector of peptide identifications
     */
    void setPeptideIdentifications(std::vector<PeptideIdentification>&& ids);

    /**
     * @brief Set all peptide identifications for all spectra
     * @param ids Vector of peptide identifications
     */
    void setPeptideIdentifications(const std::vector<PeptideIdentification>& ids);

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
      return PairIterator(data.getSpectra().cbegin(), peptide_ids.cbegin());
    }

    /**
     * @brief Get an iterator to the beginning of the data
     * @return An iterator to the beginning
     */
    inline auto begin()
    {
      return PairIterator(data.getSpectra().begin(), peptide_ids.begin());
    }

    /**
     * @brief Get a const iterator to the beginning of the data
     * @return A const iterator to the beginning
     */
    inline auto begin() const
    {
      return PairIterator(data.getSpectra().cbegin(), peptide_ids.cbegin());
    }

    /**
     * @brief Get an iterator to the end of the data
     * @return An iterator to the end
     */
    inline auto end()
    {
      return PairIterator(data.getSpectra().end(), peptide_ids.end());
    }

    /**
     * @brief Get a const iterator to the end of the data
     * @return A const iterator to the end
     */
    inline auto end() const
    {
      return PairIterator(data.getSpectra().end(), peptide_ids.end());
    }

    /**
     * @brief Get a const iterator to the end of the data
     * @return A const iterator to the end
     */
    inline auto cend() const
    {
      return PairIterator(data.getSpectra().cend(), peptide_ids.cend());
    }

    /**
     * @brief Access a spectrum and its associated peptide identification
     * @param idx The index of the spectrum
     * @return A pair of references to the spectrum and its peptide identification
     */
    inline Mapping operator[](size_t idx)
    {
      return {data.getSpectra()[idx], peptide_ids[idx]};
    }

    /**
     * @brief Access a spectrum and its associated peptide identification (const version)
     * @param idx The index of the spectrum
     * @return A pair of const references to the spectrum and its peptide identification
     */
    inline ConstMapping operator[](size_t idx) const
    {
      return {data.getSpectra()[idx], peptide_ids[idx]};
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

    typedef AnnotatedMSRun::PairIterator<std::vector<MSSpectrum>::iterator, std::vector<PeptideIdentification>::iterator> Iterator;
    typedef AnnotatedMSRun::PairIterator<std::vector<MSSpectrum>::const_iterator, std::vector<PeptideIdentification>::const_iterator> ConstIterator;

  private:
    std::vector<PeptideIdentification> peptide_ids;
    std::vector<ProteinIdentification> protein_ids_;
    MSExperiment data;
  };
}