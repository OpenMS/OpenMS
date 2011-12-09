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


#ifndef OPENMS_CHEMISTRY_SVMTHEORETICALSPECTRUMGENERATOR_H
#define OPENMS_CHEMISTRY_SVMTHEORETICALSPECTRUMGENERATOR_H

#include <OpenMS/config.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/SIMULATION/SimTypes.h>
#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>
#include <boost/smart_ptr.hpp>





namespace OpenMS
{
  /**
   @brief Simulates ms2 spectra with support vector machines.

   The models are generated with the @ref UTILS_SvmTheoreticalSpectrumGeneratorTrainer. \n
   Two modes are supported:\n
   - Only a prediction of abundant/missing is performed and for abundant peaks are generated with user defined intensity. \n
   - The intensity is predicted using SVM-regression (only for the primary ion types b and y). For the secondary types a Bayesian model is used.

   <p>
   Currently, only a test model is shipped with OpenMS.<br>
   Please find trained models at: http://sourceforge.net/projects/open-ms/files/Supplementary/Simulation/.
   </p>


   @htmlinclude OpenMS_SvmTheoreticalSpectrumGenerator.parameters

   @ingroup Chemistry
   */
  class OPENMS_DLLAPI SvmTheoreticalSpectrumGenerator: public DefaultParamHandler
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
        if (this != &rhs)
        {
          residue = rhs.residue;
          loss = rhs.loss;
          charge = rhs.charge;
        }
        return *this;
      }

      bool operator < (const IonType &rhs) const
      {
        if (residue != rhs.residue)
          return residue<rhs.residue;
        else if (loss.getString()!=rhs.loss.getString())
          return loss.getString()<rhs.loss.getString();
        else
          return charge <rhs.charge;
      }
    };
    //@}

    /// A set of descriptors for a single training row
    struct DescriptorSet
    {
      typedef std::vector<svm_node> DescriptorSetType;
      DescriptorSetType descriptors;
    };


    /// Simple container storing the model parameters required for simulation
    struct SvmModelParameterSet
    {
      //pointers to the svm classification models (one per ion_type)
      std::vector<boost::shared_ptr<SVMWrapper> > class_models;

      //pointers to the svm regression models (one per ion_type)
      std::vector<boost::shared_ptr<SVMWrapper> > reg_models;

      //The intensity for each ion type for the SVC mode
      std::map<Residue::ResidueType, DoubleReal> static_intensities;

      //The selected primary IonTypes
      std::vector<IonType> ion_types;

      //The selected secondary IonTypes
      std::map<IonType, std::vector<IonType> > secondary_types;

      //The number of intensity levels
      Size number_intensity_levels;

      //The number of regions for every spectrum
      Size number_regions;

      //upper limits (required for scaling)
      std::vector<DoubleReal> feature_max;

      //lower limits (required for scaling)
      std::vector<DoubleReal> feature_min;

      //lower bound for scaling
      double scaling_lower;

      //upper bound for scaling
      double scaling_upper;

      //border values for binning secondary types intensity
      std::vector<DoubleReal> intensity_bin_boarders;

      //intensity values for binned secondary types intensity
      std::vector<DoubleReal> intensity_bin_values;

      //conditional probabilities for secondary types
      std::map<std::pair<IonType, Size>, std::vector<std::vector<DoubleReal> > >conditional_prob;
    };



      /** @name Constructors and Destructors
       */
      //@{
      /// Default constructor
      SvmTheoreticalSpectrumGenerator();

      /// Copy constructor
      SvmTheoreticalSpectrumGenerator(const SvmTheoreticalSpectrumGenerator& source);

      /// Assignment operator
      SvmTheoreticalSpectrumGenerator& operator =(const SvmTheoreticalSpectrumGenerator& tsg);


      /// Destructor
      virtual ~SvmTheoreticalSpectrumGenerator();
      //@}


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

      /// charge of the precursors used for training
      Size precursor_charge_;

      /// set of model parameters read from model file
      SvmModelParameterSet mp_;

      /// map AA to integers
      static std::map<String,Size> aa_to_index_;

      /// hydrophobicity values for each AA
      static std::map<String,DoubleReal> hydrophobicity_;

      /// helicity values for each AA
      static std::map<String,DoubleReal> helicity_;

      /// basicity values for each AA
      static std::map<String,DoubleReal> basicity_;

      /// whether ion types are hidden or not
      std::map<IonType, bool> hide_type_;

      /// scale value to the intervall [lower,max] given the maximal and minimal entries for a feature
      inline void scaleSingleFeature_(double &value, double feature_min, double feature_max, double lower =-1.0, double upper=1.0);

      /// scale value to the intervall [lower,max] given the maximal and minimal entries for a feature
      void scaleDescriptorSet_(DescriptorSet &desc, double lower, double upper);

      /// generate the desciptors for an input peptide and a given fragmentation position
      Size generateDescriptorSet_(AASequence peptide, Size position, IonType type, Size precursor_charge, DescriptorSet &desc_set);

      /// Returns the ResidueType (e.g. AIon, BIon) as string for peak annotation
      String ResidueTypeToString_(Residue::ResidueType type);

      /// initialized the maps
      static void initializeMaps_();

      /// flag to indicate if the hydrophobicity, helicity, and basicity maps were already initialized
      static bool initializedMaps_;

      void updateMembers_();
  };

  void inline SvmTheoreticalSpectrumGenerator::scaleSingleFeature_(double &value, double lower, double upper, double feature_min, double feature_max)
  {
    double prev = value;
    if (feature_max == feature_min)
    {
      return;
    }

    if (value <= feature_min)
    {
      value = lower;
    }
    else if (value >= feature_max)
    {
      value = upper;
    }
    else
    {
      value = lower + (upper-lower) *
              (value-feature_min)/
              (feature_max-feature_min);
    }

    if (value < 0)
    {
      std::cerr << "negative value!! "<<value<<"  l: "<<lower<<" u: "<<upper<<" fm: "<<feature_min<<" fma: "<<feature_max<<"  prev: "<<prev<<std::endl;
    }
  }

} // namespace OpenMS

#endif // #ifdef OPENMS_CHEMISTRY_SVMTHEORETICALSPECTRUMGENERATORTRAINER_H

