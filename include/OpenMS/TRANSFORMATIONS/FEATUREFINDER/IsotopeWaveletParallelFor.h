// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Rene Hussong$
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELETPARALLELFOR_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELETPARALLELFOR_H

#include <tbb/blocked_range.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmIsotopeWavelet.h>

#ifndef NULL
#define NULL 0
#endif


namespace OpenMS
{
  template <typename PeakType, typename FeatureType>
  class FeatureFinderAlgorithmIsotopeWavelet;

  /** @brief A class for distributing the data over several GPUs using Intel Threading Building Blocks. */
  template <typename PeakType, typename FeatureType>
  class IsotopeWaveletParallelFor
  {
public:

    /** @brief Constructor. */
    IsotopeWaveletParallelFor(std::vector<IsotopeWaveletTransform<PeakType> *> & iwts, FeatureFinderAlgorithmIsotopeWavelet<PeakType, FeatureType> * ff) :
      iwts_(iwts), ff_(ff)
    {
    }

    /** @brief The working horse of the class.
        * The operator initializes the computation on the individual GPU. */
    void operator()(const tbb::blocked_range<size_t> & r) const
    {
      for (size_t t = r.begin(); t != r.end(); ++t)       //this will be essentially one iteration
      {
        cudaSetDevice(ff_->gpu_ids_[t]);
        IsotopeWaveletTransform<PeakType> * c_iwt = iwts_[t];

        UInt num_gpus = ff_->gpu_ids_.size();
        UInt block_size = (int)(ff_->map_->size() / num_gpus); UInt additional = ff_->map_->size() - num_gpus * block_size;
        UInt from = t * block_size;
        UInt up_to = (t >= num_gpus - 1) ? from + block_size + additional : from + block_size;

        for (UInt i = from; i < up_to; ++i)
        {
          const MSSpectrum<PeakType> & c_ref((*ff_->map_)[i]);
          if (c_ref.size() <= 1)               //unable to transform anything
          {
            //need to do this atomic
            ff_->ff_->setProgress(ff_->progress_counter_ += 2);
            continue;
          }

          bool success = true;
          typename IsotopeWaveletTransform<PeakType>::TransSpectrum * c_trans(NULL);
          if (!ff_->hr_data_)               //LowRes data
          {
            std::cout << "Parallel for: here we are" << std::endl;

            c_trans = new typename IsotopeWaveletTransform<PeakType>::TransSpectrum(&(*ff_->map_)[i]);
            success = c_iwt->initializeScanCuda((*ff_->map_)[i]) == Constants::CUDA_INIT_SUCCESS;

            if (success)
            {
              for (UInt c = 0; c < ff_->max_charge_; ++c)
              {
                c_iwt->getTransformCuda(*c_trans, c);

#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
                std::stringstream stream;
                stream << "gpu_lowres_" << ((*ff_->map_)[i]).getRT() << "_" << c + 1 << ".trans\0";
                std::ofstream ofile(stream.str().c_str());
                for (UInt k = 0; k < c_trans->size(); ++k)
                {
                  ofile << c_trans->getMZ(k) << "\t" <<  c_trans->getTransIntensity(k) << "\t" << c_trans->getMZ(k) << "\t" << c_trans->getRefIntensity(k) << std::endl;
                }
                ofile.close();
#endif

#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
                std::cout << "cuda transform for charge " << c + 1 << "  O.K. ... "; std::cout.flush();
#endif
                ff_->ff_->setProgress(++ff_->progress_counter_);

                c_iwt->identifyChargeCuda(*c_trans, i, c, ff_->intensity_threshold_, ff_->check_PPMs_);

#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
                std::cout << "cuda charge recognition for charge " << c + 1 << " O.K." << std::endl;
#endif
                ff_->ff_->setProgress(++ff_->progress_counter_);
              }
              c_iwt->finalizeScanCuda();
            }
            else
            {
              std::cout << "Warning/Error generated at scan " << i << " (" << ((*ff_->map_)[i]).getRT() << ")." << std::endl;
            }
          }
          else               //HighRes data
          {
            c_trans = ff_->prepareHRDataCuda(i, c_iwt);
            for (UInt c = 0; c < ff_->max_charge_; ++c)
            {
              c_iwt->getTransformCuda(*c_trans, c);

#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
              std::stringstream stream;
              stream << "gpu_highres_" << ((*ff_->map_)[i]).getRT() << "_" << c + 1 << ".trans\0";
              std::ofstream ofile(stream.str().c_str());
              for (UInt k = 0; k < c_trans->size(); ++k)
              {
                ofile << c_trans->getMZ(k) << "\t" <<  c_trans->getTransIntensity(k) << "\t" << c_trans->getMZ(k) << "\t" << c_trans->getRefIntensity(k) << std::endl;
              }
              ofile.close();
#endif

#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
              std::cout << "cuda transform for charge " << c + 1 << "  O.K. ... "; std::cout.flush();
#endif
              ff_->ff_->setProgress(++ff_->progress_counter_);

              c_iwt->identifyChargeCuda(*c_trans, i, c, ff_->intensity_threshold_, ff_->check_PPMs_);

#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
              std::cout << "cuda charge recognition for charge " << c + 1 << " O.K." << std::endl;
#endif
              ff_->ff_->setProgress(++ff_->progress_counter_);
            }

            c_trans->destroy();
            c_iwt->finalizeScanCuda();
          }

          delete (c_trans); c_trans = NULL;

          c_iwt->updateBoxStates(*ff_->map_, i, ff_->RT_interleave_, ff_->real_RT_votes_cutoff_, from, up_to - 1);
#ifdef OPENMS_DEBUG_ISOTOPE_WAVELET
          std::cout << "updated box states." << std::endl;
#endif

          std::cout.flush();
        }

        c_iwt->updateBoxStates(*ff_->map_, INT_MAX, ff_->RT_interleave_, ff_->real_RT_votes_cutoff_);
      }
    }

private:

    std::vector<IsotopeWaveletTransform<PeakType> *> & iwts_;
    FeatureFinderAlgorithmIsotopeWavelet<PeakType, FeatureType> * ff_;

  };   //class


} //namespace
#endif
