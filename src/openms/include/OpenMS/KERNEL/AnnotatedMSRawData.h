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

  class OPENMS_DLLAPI AnnotatedMSRawData
  {
  public:
    typedef std::pair<MSSpectrum&, std::vector<PeptideIdentification>&> Mapping;
    typedef std::pair<const MSSpectrum&, const std::vector<PeptideIdentification>&> ConstMapping;

    /// Default constructor
    AnnotatedMSRawData() = default;

    /// Move constructor for efficiently loading a MSExperiment without a deep copy.
    explicit AnnotatedMSRawData(MSExperiment&& experiment) : data(std::move(experiment))
    {};

    AnnotatedMSRawData(AnnotatedMSRawData&&) = default;

    ~AnnotatedMSRawData() = default;

    std::vector<ProteinIdentification>& getProteinIdentifications()
    {
      return protein_ids_;
    }

    const std::vector<ProteinIdentification>& getProteinIdentifications() const
    {
      return protein_ids_;
    }
    /// Get the peptide identifications for a single spectrum.
    std::vector<PeptideIdentification>& getPeptideIdentifications(size_t index);
    
    /// Get the peptide identifications for a single spectrum (const version).
    const std::vector<PeptideIdentification>& getPeptideIdentifications(size_t index) const;

    /// Get all peptide identifications for all spectra.
    std::vector<std::vector<PeptideIdentification>>& getAllPeptideIdentifications();
    
    /// Get all peptide identifications for all spectra (const version).
    const std::vector<std::vector<PeptideIdentification>>& getAllPeptideIdentifications() const;

    /// Set a single spectrum's peptide identification annotation
    void setPeptideIdentifications(std::vector<PeptideIdentification>&& ids, size_t index);

    /// Set all peptide identifications for all spectra
    void setAllPeptideIdentifications(std::vector<std::vector<PeptideIdentification>>&& ids);

    void clearAllPeptideIdentifications()
    {
      std::vector<std::vector<PeptideIdentification>> empty_ids;
      peptide_ids.swap(empty_ids);
    }

    MSExperiment& getMSExperiment();
    
    /// Get the MSExperiment (const version).
    const MSExperiment& getMSExperiment() const;

    inline auto cbegin() const
    {
      return PairIterator(data.getSpectra().cbegin(), peptide_ids.cbegin());
    }

    inline auto begin()
    {
      return PairIterator(data.getSpectra().begin(), peptide_ids.begin());
    }

    inline auto begin() const
    {
      return PairIterator(data.getSpectra().cbegin(), peptide_ids.cbegin());
    }

    inline auto end()
    {
      return PairIterator(data.getSpectra().end(), peptide_ids.end());
    }

    inline auto end() const
    {
      return PairIterator(data.getSpectra().end(), peptide_ids.end());
    }

    inline auto cend() const
    {
      return PairIterator(data.getSpectra().cend(), peptide_ids.cend());
    }

    inline Mapping operator[](size_t idx)
    {
      return {data.getSpectra()[idx], peptide_ids[idx]};
    }

    inline ConstMapping operator[](size_t idx) const
    {
      return {data.getSpectra()[idx], peptide_ids[idx]};
    }

    template<typename T1, typename T2>
    struct PairIterator
    {
      // TODO add check that both vectors are of the same length
      using iterator_category = std::forward_iterator_tag;
      using difference_type = std::ptrdiff_t;
      //using value_type = std::pair<T1, T2>;
      //using pointer = value_type*;
      //using reference = value_type&;

      PairIterator(T1 ptr1, T2 ptr2) : m_ptr1(ptr1), m_ptr2(ptr2)
      {}

      PairIterator& operator++()
      {
        ++m_ptr1;
        ++m_ptr2;
        return *this;
      }

      PairIterator operator++(int)
      {
        auto tmp(*this);
        ++(*this);
        return tmp;
      }

      auto operator*()
      {
        return std::make_pair(std::ref(*m_ptr1), std::ref(*m_ptr2));
      }

      inline friend bool operator==(const PairIterator& a, const PairIterator& b)
      {
        return a.m_ptr1 == b.m_ptr1 && a.m_ptr2 == b.m_ptr2;
      }

      inline friend bool operator!=(const PairIterator& a, const PairIterator& b)
      {
        return !(a == b);
      }

    private:
      T1 m_ptr1;
      T2 m_ptr2;
    };

    typedef AnnotatedMSRawData::PairIterator<std::vector<MSSpectrum>::iterator, std::vector<std::vector<PeptideIdentification>>::iterator> Iterator;
    typedef AnnotatedMSRawData::PairIterator<std::vector<MSSpectrum>::const_iterator, std::vector<std::vector<PeptideIdentification>>::const_iterator> ConstIterator;

  private:
    std::vector<std::vector<PeptideIdentification>> peptide_ids;
    std::vector<ProteinIdentification> protein_ids_;
    MSExperiment data;
  };
}
