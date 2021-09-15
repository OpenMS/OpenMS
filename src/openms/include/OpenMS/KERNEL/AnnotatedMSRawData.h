// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: David Voigt $
// $Authors: David Voigt $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <boost/range/combine.hpp>

#include <vector>

namespace OpenMS
{
  class PeptideIdentification;
  class MSSpectrum;

  class OPENMS_DLLAPI AnnotatedMSRawData : public MSExperiment
  {
  public:
    typedef std::pair<MSSpectrum&, std::vector<PeptideIdentification>&> Mapping;
    typedef std::pair<const MSSpectrum&, const std::vector<PeptideIdentification>&> ConstMapping;

    /// Default constructor
    AnnotatedMSRawData() = default;
    /// Move constructor for efficiently loading a MSExperiment without a deep copy.
    explicit AnnotatedMSRawData(MSExperiment&& experiment) : MSExperiment(std::move(experiment)) {};
    AnnotatedMSRawData(AnnotatedMSRawData&&) = default;
    ~AnnotatedMSRawData() = default;

    /// Get the peptide identifications for a single spectrum.
    std::vector<PeptideIdentification>& getPeptideIds();
    /// Get all peptide identifications for all spectra.
    std::vector<std::vector<PeptideIdentification>>& getAllPeptideIds();

    /// Set a single spectrum's peptide identification annotation
    void setPeptideIds(std::vector<PeptideIdentification>&& ids, size_t index);
    /// Set all peptide identifications for all spectra
    void setAllPeptideIds(std::vector<std::vector<PeptideIdentification>>&& ids);

    inline auto cbegin() const
    {
      return PairIterator(spectra_.cbegin(), peptide_ids.cbegin());
    }

    inline auto begin()
    {
      return PairIterator(spectra_.begin(), peptide_ids.begin());
    }

    inline auto end()
    {
      return PairIterator(spectra_.end(), peptide_ids.end());
    }

    inline auto cend() const
    {
      return PairIterator(spectra_.cend(), peptide_ids.cend());
    }

    inline Mapping operator[](size_t idx)
    {
      return { spectra_[idx], peptide_ids[idx] };
    }

    inline ConstMapping operator[](size_t idx) const
    {
      return { spectra_[idx], peptide_ids[idx] };
    }

    template<typename T1, typename T2>
    struct PairIterator
    {
      // TODO add check that both vectors are of the same length
      using iterator_category = std::forward_iterator_tag;
      using difference_type = std::ptrdiff_t;
      using value_type = std::pair<T1, T2>;
      //using pointer = value_type*;
      //using reference = value_type&;

      PairIterator(T1 ptr1, T2 ptr2) : m_ptr1(ptr1), m_ptr2(ptr2) {}

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
        return std::pair(*m_ptr1, *m_ptr2);
      }

      inline friend bool operator== (const PairIterator& a, const PairIterator& b)
      {
        return a.m_ptr1 == b.m_ptr1 && a.m_ptr2 == b.m_ptr2;
      }

      inline friend bool operator!= (const PairIterator& a, const PairIterator& b)
      {
        return !(a == b);
      }

    private:
      T1 m_ptr1;
      T2 m_ptr2;
  };

  typedef AnnotatedMSRawData::PairIterator<std::vector<MSSpectrum>::iterator , std::vector<std::vector<PeptideIdentification>>::iterator> Iterator;
  typedef AnnotatedMSRawData::PairIterator<std::vector<MSSpectrum>::const_iterator , std::vector<std::vector<PeptideIdentification>>::const_iterator> ConstIterator;

  private:
    std::vector<std::vector<PeptideIdentification>> peptide_ids;
  };
}