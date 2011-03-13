// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Sandro Andreotti $
// $Authors: Sandro Andreotti $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_SVMTHEORETICALSPECTRUMGENERATORTRAINER_H
#define OPENMS_CHEMISTRY_SVMTHEORETICALSPECTRUMGENERATORTRAINER_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGenerator.h>
#include <OpenMS/SIMULATION/SimTypes.h>
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
   Default values are choses as in the svm-training tool delivered with libSVM.\n

   For the secondary types (a, c, x, z, losses, b2, y2) a simple Bayesian model is used.

   @ingroup Chemistry
   */

  class OPENMS_DLLAPI SvmTheoreticalSpectrumGeneratorTrainer : public DefaultParamHandler
  {
    typedef SvmTheoreticalSpectrumGenerator::IonType IonType;
    typedef SvmTheoreticalSpectrumGenerator::DescriptorSet DescriptorSet;
    typedef std::map<std::pair<IonType,Size>, std::vector<DoubleReal> > ObservedIntensMap;

    /// stores the observed intensities for each sector-type combination in a vector
    void countIntensities_(const PeakSpectrum &spectrum,
                          const AASequence &annotation,
                          IonType type,
                          std::map<std::pair<IonType, Size>, std::vector<DoubleReal> > & observed_intensities,
                          DoubleReal tolerance,
                          Size number_of_regions
                          );

    /// trains the Bayesian secondary peak types models
    void trainSecondaryTypes_(TextFile &info_outfile,
                              Size number_of_regions,
                              Size number_of_intensity_levels,
                              ObservedIntensMap &observed_intensities,
                              const std::vector<IonType> &ion_types,
                              const std::vector<bool> &is_primary
                              );
    public:

      /** @name Constructors and Destructors
       */
      //@{
      /// Default constructor
      SvmTheoreticalSpectrumGeneratorTrainer();

      /// Copy constructor
      SvmTheoreticalSpectrumGeneratorTrainer(const SvmTheoreticalSpectrumGeneratorTrainer& source);

      /// Destructor
      virtual ~SvmTheoreticalSpectrumGeneratorTrainer();
      //@}

      /// Assignment operator
      SvmTheoreticalSpectrumGeneratorTrainer& operator =(const SvmTheoreticalSpectrumGeneratorTrainer& tsg);

      /// trains an SVM for each ion_type and stores them in files <filename>_residue_loss_charge.svm
      void trainModel(const PeakMap &spectra, const std::vector<AASequence> & annotations, String filename, Size precursor_charge);

      /// Normalizes the intensity of the peaks in the input data
      void normalizeIntensity(PeakSpectrum &S) const;

    protected:

      /// Write a training file that can be passed to libsvm command line tools
      void write_training_file_(std::vector<DescriptorSet> &training_input, std::vector<DoubleReal> &training_output, String filename);

  };
}

#endif // SvmTheoreticalSpectrumGeneratorTrainer_H
