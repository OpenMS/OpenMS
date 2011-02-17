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

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/SIMULATION/SimTypes.h>

#include<svm.h>


namespace OpenMS
{
  /**
   @brief Generates theoretical spectra according to a artificial neural network.

   The models are generated with the @ref UTILS_SpectrumGeneratorNetworkTrainer

   @htmlinclude OpenMS_SvmTheoreticalSpectrumGenerator.parameters

   @ingroup Chemistry
   */
  class OPENMS_DLLAPI SvmTheoreticalSpectrumGenerator: public TheoreticalSpectrumGenerator
  {
    friend class SvmTheoreticalSpectrumGeneratorTrainer;
    public:    

    /**
     @brief nested class
    */
    //@{
    ///IonType is defined by a ResidueType, a neutral loss and charge state
    struct IonType
    {
      Residue::ResidueType residue;
      EmpiricalFormula loss;
      Int charge;

      /** @name Constructors and Destructors*/
      //@{
      //Default constructor
      IonType():
        residue((Residue::ResidueType) 0),
        loss(),
        charge(0)
        {
        }

      //Custom construtor
      IonType(Residue::ResidueType residue, EmpiricalFormula loss= EmpiricalFormula(), Int charge=1):
        residue(residue),
        loss(loss),
        charge(charge)
        {
        }

      //Copy constructor
      IonType(const IonType & rhs):
        residue(rhs.residue),
        loss(rhs.loss),
        charge(rhs.charge)
        {
        }

      //Assignment operator
      IonType & operator = (const IonType & rhs)
      {
        if(this != &rhs)
        {
          residue = rhs.residue;
          loss = rhs.loss;
          charge = rhs.charge;
        }
        return *this;
      }

      bool operator < (const IonType &rhs) const
      {
        if(residue != rhs.residue)
          return residue<rhs.residue;
        else if(loss.getString()!=rhs.loss.getString())
          return loss.getString()<rhs.loss.getString();
        else
          return charge <rhs.charge;
      }
    };
    //@}

    struct DescriptorSet
    {
      typedef std::vector<svm_node> DescriptorSetContainerType;
      DescriptorSetContainerType descriptors;
    };

    struct SvmModelParameterSet
    {
      //pointers to the svm classification models (one per ion_type)
      std::vector<svm_model*>class_models;

      //pointers to the svm regression models (one per ion_type)
      std::vector<svm_model*>reg_models;

      //The intensity for each ion type for the SVC mode
      std::map<Residue::ResidueType, DoubleReal>static_intensities;

      ///The selected primary IonTypes
      std::vector<IonType>ion_types;

      ///The selected secondary IonTypes
      std::map<IonType, std::vector<IonType> >secondary_types;

      ///The number of intensity levels
      Size number_intensity_levels;

      ///The number of regions for every spectrum
      Size number_regions;

      ///upper limits (required for scaling)
      std::vector<DoubleReal>feature_max;

      ///lower limits (required for scaling)
      std::vector<DoubleReal>feature_min;

      ///lower bound for scaling
      double scaling_lower;

      ///upper bound for scaling
      double scaling_upper;

      ///border values for binning secondary types intensity
      std::vector<DoubleReal>intensity_bin_boarders;

      ///intensity values for binned secondary types intensity
      std::vector<DoubleReal>intensity_bin_values;

      ///conditional probabilities for secondary types
      std::map<std::pair<IonType, Size>, std::vector<std::vector<DoubleReal> > >conditional_prob;
    };


      /** @name Constructors and Destructors
       */
      //@{
      /// Default constructor
      SvmTheoreticalSpectrumGenerator();

      /// Copy constructor
      SvmTheoreticalSpectrumGenerator(const SvmTheoreticalSpectrumGenerator& source);

      /// Destructor
      virtual ~SvmTheoreticalSpectrumGenerator();
      //@}

      /// Assignment operator
      SvmTheoreticalSpectrumGenerator& operator =(const SvmTheoreticalSpectrumGenerator& tsg);

      /// Generate the MS/MS according to the given probabilistic model
      void simulate(RichPeakSpectrum &spectrum, const AASequence &peptide, const gsl_rng *rng, Size precursor_charge);

      ///Load a trained Svm and Prob. models
      void load();

      ///return the set of ion types that are modeled by the loaded SVMs
      const std::vector<IonType>& getIonTypes()
      {
        return mp_.ion_types;
      }


    protected:
      typedef std::map<IonType, DoubleReal>IntensityMap;      

      ///charge of the precursors used for training
      Size precursor_charge_;

      ///set of model parameters read from model file
      SvmModelParameterSet mp_;

      ///scale value to the intervall [lower,max] given the maximal and minimal entries for a feature
      inline void scaleSingleFeature_(double &value, double feature_min, double feature_max, double lower =-1.0, double upper=1.0);

      ///scale value to the intervall [lower,max] given the maximal and minimal entries for a feature
      void scaleDescriptorSet_(DescriptorSet &desc, double lower, double upper);      

      ///generate the desciptors for an input peptide and a given fragmentation position
      Size generateDescriptorSet(AASequence peptide, Size position, IonType type, Size precursor_charge, DescriptorSet &desc_set);

      ///Returns the ResidueType (e.g. AIon, BIon) as string for peak annotation
      String ResidueTypeToString_(Residue::ResidueType type);

      void updateMembers_();
  };

  void inline SvmTheoreticalSpectrumGenerator::scaleSingleFeature_(double &value, double lower, double upper, double feature_min, double feature_max)
  {    
    //std::cerr<<"vA: "<<value<< " l: "<<lower<<" u: "<<upper<<" fm: "<<feature_min<<" fma: "<<feature_max<<std::endl;
    /* skip single-valued attribute */
    double prev=value;
    if(feature_max == feature_min)
      return;

    if(value <= feature_min)
      value = lower;
    else if(value >= feature_max)
      value = upper;
    else
      value = lower + (upper-lower) *
              (value-feature_min)/
              (feature_max-feature_min);
    //std::cerr<<"vB: "<<value<< " l: "<<lower<<" u: "<<upper<<" fm: "<<feature_min<<" fma: "<<feature_max<<std::endl;
    if(value<0)
    {
      std::cerr<<"negative value!! "<<value<<"  l: "<<lower<<" u: "<<upper<<" fm: "<<feature_min<<" fma: "<<feature_max<<"  prev: "<<prev<<std::endl;
    }
  }

} // namespace OpenMS

#endif // #ifdef OPENMS_CHEMISTRY_SVMTHEORETICALSPECTRUMGENERATORTRAINER_H

