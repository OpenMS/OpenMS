// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Sandro Andreotti $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_SVMTHEORETICALSPECTRUMGENERATORTRAINER_H
#define OPENMS_CHEMISTRY_SVMTHEORETICALSPECTRUMGENERATORTRAINER_H

#include <OpenMS/KERNEL/StandardDeclarations.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGenerator.h>
#include <OpenMS/FORMAT/TextFile.h>

namespace OpenMS
{
  /**
   @brief Train SVM models that are used by SvmTheoreticalSpectrumGenerator

   @htmlinclude OpenMS_SvmTheoreticalSpectrumGeneratorTrainer.parameters

   This class implements the algorithm used by the homonymous tool which can be
   used to train models for MS/MS spectrum simulation.\n
   For the primary ion types (y, b) a SVM is trained using the libSVM library.\n
   All important libSVM parameters are accessible as parameters.\n
   Please refer to the libSVM manuals for detailed description of the parameters.
   Default values are chosen as in the svm-training tool delivered with libSVM.\n

   For the secondary types (a, c, x, z, losses, b2, y2) a simple Bayesian model is used.

   @ingroup Chemistry
   */

  class OPENMS_DLLAPI SvmTheoreticalSpectrumGeneratorTrainer :
    public DefaultParamHandler
  {
    typedef SvmTheoreticalSpectrumGenerator::IonType IonType;
    typedef SvmTheoreticalSpectrumGenerator::DescriptorSet DescriptorSet;
    typedef std::map<std::pair<IonType, Size>, std::vector<double> > ObservedIntensMap;

    /// stores the observed intensities for each sector-type combination in a vector
    void countIntensities_(const PeakSpectrum & spectrum,
                           const AASequence & annotation,
                           IonType type,
                           std::map<std::pair<IonType, Size>, std::vector<double> > & observed_intensities,
                           double tolerance,
                           Size number_of_regions
                           );

    /// trains the Bayesian secondary peak types models
    void trainSecondaryTypes_(TextFile & info_outfile,
                              Size number_of_regions,
                              Size number_of_intensity_levels,
                              ObservedIntensMap & observed_intensities,
                              const std::vector<IonType> & ion_types,
                              const std::vector<bool> & is_primary
                              );
public:

    /** @name Constructors and Destructors
     */
    //@{
    /// Default constructor
    SvmTheoreticalSpectrumGeneratorTrainer();

    /// Copy constructor
    SvmTheoreticalSpectrumGeneratorTrainer(const SvmTheoreticalSpectrumGeneratorTrainer & source);

    /// Destructor
    ~SvmTheoreticalSpectrumGeneratorTrainer() override;
    //@}

    /// Assignment operator
    SvmTheoreticalSpectrumGeneratorTrainer & operator=(const SvmTheoreticalSpectrumGeneratorTrainer & tsg);

    /// trains an SVM for each ion_type and stores them in files \<filename\>_residue_loss_charge.svm
    void trainModel(const PeakMap & spectra, const std::vector<AASequence> & annotations, String filename, Int precursor_charge);

    /// Normalizes the intensity of the peaks in the input data
    void normalizeIntensity(PeakSpectrum & S) const;

protected:

    /// Write a training file that can be passed to libsvm command line tools
    void writeTrainingFile_(std::vector<DescriptorSet> & training_input, std::vector<double> & training_output, String filename);

  };
}

#endif // SvmTheoreticalSpectrumGeneratorTrainer_H
