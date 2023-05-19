// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest, Witold Wolski $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OPENSWATHALGO/OpenSwathAlgoConfig.h>

#include <OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <boost/shared_ptr.hpp>
#include <string>
#include <vector>

namespace OpenSwath
{

  using SpectrumSequence = std::vector<OpenSwath::SpectrumPtr>;
  /**
    @brief The interface of a mass spectrometry experiment.
  */
  class OPENSWATHALGO_DLLAPI ISpectrumAccess
  {
public:
    /// Destructor
    virtual ~ISpectrumAccess();

    /**
      @brief Light clone operator to produce a copy for concurrent read access.

      This function guarantees to produce a copy of the underlying object that
      provides thread-safe concurrent read access to the underlying data. It
      should be implemented with minimal copy-overhead to make this operation
      as fast as possible.

      To use this function, each thread should call this function to produce an
      individual copy on which it can operate.

    */
    virtual boost::shared_ptr<ISpectrumAccess> lightClone() const = 0;

    /// Return a pointer to a spectrum at the given id
    virtual SpectrumPtr getSpectrumById(int id) = 0;

    /// Return pointer to a spectrum at the given id, the spectrum will be filtered by drift time
    SpectrumPtr getSpectrumById(int id, double drift_start, double drift_end );

    /// Return a vector of ids of spectra that are within RT +/- deltaRT
    virtual std::vector<std::size_t> getSpectraByRT(double RT, double deltaRT) const = 0;
    /// Returns the number of spectra available
    virtual size_t getNrSpectra() const = 0;
    /// Returns the meta information for a spectrum
    virtual SpectrumMeta getSpectrumMetaById(int id) const = 0;

    /// Return a pointer to a chromatogram at the given id
    virtual ChromatogramPtr getChromatogramById(int id) = 0;
    /// Returns the number of chromatograms available
    virtual std::size_t getNrChromatograms() const = 0;
    /// Returns the native id of the chromatogram at the given id
    virtual std::string getChromatogramNativeID(int id) const = 0;

    /* @brief Fetches a spectrumSequence (multiple spectra pointers) closest to the given RT
     * @p RT = target RT
     * @p nr_spectra_to_fetch = # spectra around target RT to fetch (length of the spectrum sequence)
    */
    SpectrumSequence getMultipleSpectra(double RT, int nr_spectra_to_fetch);

    /* @brief Fetches a spectrumSequence (multiple spectra pointers) closest to the given RT. Filters all spectra by specified @p drift_start and @p drift_end
     * @p RT = target RT
     * @p nr_spectra_to_fetch = # spectra around target RT to fetch (length of the spectrum sequence)
    */
    SpectrumSequence getMultipleSpectra(double RT, int nr_spectra_to_fetch, double drift_start, double drift_end);

    /// filters a spectrum by drift time, spectrum pointer returned is a copy
    static SpectrumPtr filterByDrift(const SpectrumPtr& input, double drift_start, double drift_end)
    {
      // NOTE: this function is very inefficient because filtering unsorted array
      //OPENMS_PRECONDITION(drift_start <= 0, "Cannot filter by drift time if drift_start is not set");
      //OPENMS_PRECONDITION(drift_end - drift_start < 0, "Cannot filter by drift time if range is empty");
      //OPENMS_PRECONDITION(input->getDriftTimeArray() != nullptr, "Cannot filter by drift time if no drift time is available.");

      //if (input->getDriftTimeArray() == nullptr)
      //{
        //throw Exception::NullPointer(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      //}

      OpenSwath::SpectrumPtr output(new OpenSwath::Spectrum);

      OpenSwath::BinaryDataArrayPtr mz_arr = input->getMZArray();
      OpenSwath::BinaryDataArrayPtr int_arr = input->getIntensityArray();
      OpenSwath::BinaryDataArrayPtr im_arr = input->getDriftTimeArray();

      auto mz_it = mz_arr->data.cbegin();
      auto int_it = int_arr->data.cbegin();
      auto im_it = im_arr->data.cbegin();
      auto mz_end = mz_arr->data.cend();

      OpenSwath::BinaryDataArrayPtr mz_arr_out(new OpenSwath::BinaryDataArray);
      OpenSwath::BinaryDataArrayPtr intens_arr_out(new OpenSwath::BinaryDataArray);
      OpenSwath::BinaryDataArrayPtr im_arr_out(new OpenSwath::BinaryDataArray);
      im_arr_out->description = im_arr->description;

      while (mz_it != mz_end)
      {
        if ( (drift_start <= *im_it) & (drift_end >= *im_it) )
        {
          mz_arr_out->data.push_back( *mz_it );
          intens_arr_out->data.push_back( *int_it );
          im_arr_out->data.push_back( *im_it );
        }
        ++mz_it;
        ++int_it;
        ++im_it;
      }
      output->setMZArray(mz_arr_out);
      output->setIntensityArray(intens_arr_out);
      output->getDataArrays().push_back(im_arr_out);
      return output;
  }


   };

  typedef boost::shared_ptr<ISpectrumAccess> SpectrumAccessPtr;
}

