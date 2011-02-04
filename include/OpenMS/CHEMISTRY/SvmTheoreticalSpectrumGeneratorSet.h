// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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


#ifndef OPENMS_CHEMISTRY_SvmTheoreticalSpectrumGeneratorSet_H
#define OPENMS_CHEMISTRY_SvmTheoreticalSpectrumGeneratorSet_H

#include <OpenMS/SIMULATION/SimTypes.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGenerator.h>

namespace OpenMS
{
  /**
   @brief Loads SvmTheoreticalSpectrumGenerator instances for different charges

   @htmlinclude OpenMS_SvmTheoreticalSpectrumGeneratorSet.parameters

   @ingroup Chemistry
   */
  class OPENMS_DLLAPI SvmTheoreticalSpectrumGeneratorSet
  {
    public:

      /** @name Constructors and Destructors
       */
      //@{
      /// Default constructor
      SvmTheoreticalSpectrumGeneratorSet();

      /// Copy constructor
      SvmTheoreticalSpectrumGeneratorSet(const SvmTheoreticalSpectrumGeneratorSet& source);

      /// Destructor
      virtual ~SvmTheoreticalSpectrumGeneratorSet();
      //@}

      /// Assignment operator
      SvmTheoreticalSpectrumGeneratorSet& operator =(const SvmTheoreticalSpectrumGeneratorSet& tsg);

      /// Generate the MS/MS according to the model for the given precursor_charge
      void simulate(RichPeakSpectrum &spectrum, const AASequence &peptide, const gsl_rng *rng, Size precursor_charge);

      ///Load a trained Svm and Prob. models
      void load(String);

      ///Return precursor charges for which a model is contained in the set
      void getSupportedCharges(std::set<Size>&charges);

      ///return a modifiable reference to the SVM model with given charge. If charge is not supported throw exception
      SvmTheoreticalSpectrumGenerator & getSvmModel(Size);

    protected:
      //map containing the simulator for each charge variant
      std::map<Size, SvmTheoreticalSpectrumGenerator>simulators_;

  };


}



#endif
