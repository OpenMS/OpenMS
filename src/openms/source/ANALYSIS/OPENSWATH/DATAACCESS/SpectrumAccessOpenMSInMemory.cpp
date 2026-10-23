// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSInMemory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessSqMass.h>

#include <algorithm>    // std::lower_bound, std::upper_bound, std::sort
#include <exception>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS
{

  SpectrumAccessOpenMSInMemory::SpectrumAccessOpenMSInMemory(OpenSwath::ISpectrumAccess & origin)
  {
    // special case: we can grab the data directly (and fast)
    if (dynamic_cast<SpectrumAccessSqMass*> (&origin))
    {
        SpectrumAccessSqMass* tmp = dynamic_cast<SpectrumAccessSqMass*> (&origin);
        tmp->getAllSpectra(spectra_, spectra_meta_);
    }
    else
    {
      const Size nr_spectra = origin.getNrSpectra();
      int nr_threads = 1;
#ifdef _OPENMP
      // Called from a single thread (e.g. a SWATH workflow that loads one SWATH at a time), read
      // ranges of the spectra in parallel, each through its own light clone of the access (reading
      // and transforming a spectrum is the expensive part, and a clone is safe to use in parallel).
      constexpr Size min_spectra_per_thread = 64;
      if (omp_in_parallel() == 0)
      {
        nr_threads = static_cast<int>(std::min<Size>(static_cast<Size>(omp_get_max_threads()), nr_spectra / min_spectra_per_thread));
      }
#endif
      if (nr_threads <= 1)
      {
        for (Size i = 0; i < nr_spectra; ++i)
        {
          spectra_.push_back( origin.getSpectrumById(i) );
          spectra_meta_.push_back( origin.getSpectrumMetaById(i) );
        }
      }
      else
      {
        spectra_.resize(nr_spectra);
        spectra_meta_.resize(nr_spectra);
        std::exception_ptr read_exception;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1) num_threads(nr_threads)
#endif
        for (int thread = 0; thread < nr_threads; ++thread)
        {
          try
          {
            OpenSwath::SpectrumAccessPtr clone = origin.lightClone();
            const Size begin = nr_spectra * thread / nr_threads;
            const Size end = nr_spectra * (thread + 1) / nr_threads;
            for (Size i = begin; i < end; ++i)
            {
              spectra_[i] = clone->getSpectrumById(i);
              spectra_meta_[i] = clone->getSpectrumMetaById(i);
            }
          }
          catch (...)
          {
#ifdef _OPENMP
#pragma omp critical (spectrum_access_in_memory_exception)
#endif
            {
              if (read_exception == nullptr) read_exception = std::current_exception();
            }
          }
        }
        if (read_exception != nullptr)
        {
          std::rethrow_exception(read_exception);
        }
      }
      for (Size i = 0; i < origin.getNrChromatograms(); ++i)
      {
        chromatograms_.push_back( origin.getChromatogramById(i) );
        chromatogram_ids_.push_back( origin.getChromatogramNativeID(i) );
      }
    }

    OPENMS_POSTCONDITION(spectra_.size() == spectra_meta_.size(), "Spectra and meta data needs to match")
    OPENMS_POSTCONDITION(chromatogram_ids_.size() == chromatograms_.size(), "Chromatograms and meta data needs to match")
  }

  SpectrumAccessOpenMSInMemory::~SpectrumAccessOpenMSInMemory() = default;

  SpectrumAccessOpenMSInMemory::SpectrumAccessOpenMSInMemory(const SpectrumAccessOpenMSInMemory & rhs) :
    spectra_(rhs.spectra_),
    spectra_meta_(rhs.spectra_meta_),
    chromatograms_(rhs.chromatograms_),
    chromatogram_ids_(rhs.chromatogram_ids_)
  {
    // this only copies the pointers and not the actual data ... 
  }

  std::shared_ptr<OpenSwath::ISpectrumAccess> SpectrumAccessOpenMSInMemory::lightClone() const
  {
    return std::shared_ptr<SpectrumAccessOpenMSInMemory>(new SpectrumAccessOpenMSInMemory(*this));
  }

  OpenSwath::SpectrumPtr SpectrumAccessOpenMSInMemory::getSpectrumById(int id)
  {
    OPENMS_PRECONDITION(id >= 0, "Id needs to be larger than zero");
    OPENMS_PRECONDITION(id < (int)getNrSpectra(), "Id cannot be larger than number of spectra");
    return spectra_[id];
  }

  OpenSwath::SpectrumMeta SpectrumAccessOpenMSInMemory::getSpectrumMetaById(int id) const
  {
    OPENMS_PRECONDITION(id >= 0, "Id needs to be larger than zero");
    OPENMS_PRECONDITION(id < (int)getNrSpectra(), "Id cannot be larger than number of spectra");
    return spectra_meta_[id];
  }

  std::vector<std::size_t> SpectrumAccessOpenMSInMemory::getSpectraByRT(double RT, double deltaRT) const
  {
    OPENMS_PRECONDITION(deltaRT >= 0, "Delta RT needs to be a positive number");

    // we first perform a search for the spectrum that is past the
    // beginning of the RT domain. Then we add this spectrum and try to add
    // further spectra as long as they are below RT + deltaRT.
    std::vector<std::size_t> result;
    OpenSwath::SpectrumMeta s;
    s.RT = RT - deltaRT;
    auto spectrum = std::lower_bound(spectra_meta_.begin(), spectra_meta_.end(), s, OpenSwath::SpectrumMeta::RTLess());
    if (spectrum == spectra_meta_.end()) return result;

    result.push_back(std::distance(spectra_meta_.begin(), spectrum));
    ++spectrum;
    while (spectrum->RT < RT + deltaRT && spectrum != spectra_meta_.end())
    {
      result.push_back(std::distance(spectra_meta_.begin(), spectrum));
      ++spectrum;
    }
    return result;
  }

  size_t SpectrumAccessOpenMSInMemory::getNrSpectra() const
  {
    OPENMS_PRECONDITION(spectra_.size() == spectra_meta_.size(), "Spectra and meta data needs to match")
    return spectra_.size();
  }

  OpenSwath::ChromatogramPtr SpectrumAccessOpenMSInMemory::getChromatogramById(int id)
  {
    OPENMS_PRECONDITION(id >= 0, "Id needs to be larger than zero");
    OPENMS_PRECONDITION(id < (int)getNrChromatograms(), "Id cannot be larger than number of chromatograms");

    return chromatograms_[id];
  }

  size_t SpectrumAccessOpenMSInMemory::getNrChromatograms() const
  {
    OPENMS_PRECONDITION(chromatogram_ids_.size() == chromatograms_.size(), "Chromatograms and meta data needs to match")
    return chromatograms_.size();
  }

  std::string SpectrumAccessOpenMSInMemory::getChromatogramNativeID(int id) const
  {
    OPENMS_PRECONDITION(id >= 0, "Id needs to be larger than zero");
    OPENMS_PRECONDITION(id < (int)getNrChromatograms(), "Id cannot be larger than number of spectra");
    return chromatogram_ids_[id];
  }

} //end namespace OpenMS

