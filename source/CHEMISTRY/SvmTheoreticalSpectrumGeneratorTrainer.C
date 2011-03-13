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

#include <OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGeneratorTrainer.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/SqrtMower.h>
#include <stdio.h>
#include <gsl/gsl_statistics_double.h>
#include <sstream>

namespace OpenMS
{
  SvmTheoreticalSpectrumGeneratorTrainer::SvmTheoreticalSpectrumGeneratorTrainer() :
      DefaultParamHandler("SvmTheoreticalSpectrumGeneratorTrainer")
  {
    defaults_.setValue("write_training_files", "false", "If set to true no models are trained but files (<Filename>_<ion_type>_training.dat) are produced \
for the selected primary ion types. They can be used as input for LibSVM command line tools");

    defaults_.setValidStrings("write_training_files", StringList::create("true,false"));

    defaults_.setValue("number_intensity_levels", 7, "The number of intensity bins (for secondary type models)");
    defaults_.setValue("number_regions", 3, "The number of regions each spectrum is split to (for secondary type models)");
    defaults_.setValue("parent_tolerance", 2.5, "The maximum difference between theoretical and experimental parent mass to accept training spectrum");
    defaults_.setValue("peak_tolerance", 0.5, "The maximum mass error for a peak to the expected mass of some ion type");

    defaults_.setValue("add_b_ions", "true", "Train simulator for b-ions");
    defaults_.setValidStrings("add_b_ions", StringList::create("true,false"));
    defaults_.setValue("add_y_ions", "true", "Train simulator for y-ions");
    defaults_.setValidStrings("add_y_ions", StringList::create("true,false"));
    defaults_.setValue("add_a_ions", "false", "Train simulator for a-ions");
    defaults_.setValidStrings("add_a_ions", StringList::create("true,false"));
    defaults_.setValue("add_c_ions", "false", "Train simulator for c-ions");
    defaults_.setValidStrings("add_c_ions", StringList::create("true,false"));
    defaults_.setValue("add_x_ions", "false", "Train simulator for x-ions");
    defaults_.setValidStrings("add_x_ions", StringList::create("true,false"));
    defaults_.setValue("add_z_ions", "false", "Train simulator for z-ions");
    defaults_.setValidStrings("add_z_ions", StringList::create("true,false"));

    defaults_.setValue("add_losses", "false", "Train simulator for neutral losses of H2O and NH3 for b-ions and y-ions");
    defaults_.setValidStrings("add_losses", StringList::create("true,false"));
    defaults_.setValue("add_b2_ions", "false", "Train simulator for doubly charged b-ions");
    defaults_.setValidStrings("add_b2_ions", StringList::create("true,false"));
    defaults_.setValue("add_y2_ions", "false", "Train simulator for double charged y-ions");
    defaults_.setValidStrings("add_y2_ions", StringList::create("true,false"));

    defaults_.setValue("svm:svc_type", 0, "Type of the SVC: 0=C_SVC 1=NU_SVC");
    defaults_.setMinInt("svm:svc_type",0);
    defaults_.setMaxInt("svm:svc_type",1);

    defaults_.setValue("svm:svr_type", 1, "Type of the SVR: 0=EPSILON_SVR 1=NU_SVR");
    defaults_.setMinInt("svm:svr_type",0);
    defaults_.setMaxInt("svm:svr_type",1);

    defaults_.setValue("svm:svc:kernel_type", 2, "Type of the kernel:  0=LINEAR 1=POLY 2=RBF 3=SIGMOID");
    defaults_.setMinInt("svm:svc:kernel_type",0);
    defaults_.setMaxInt("svm:svc:kernel_type",3);

    defaults_.setValue("svm:svr:kernel_type", 2, "Type of the kernel:  0=LINEAR 1=POLY 2=RBF 3=SIGMOID");
    defaults_.setMinInt("svm:svr:kernel_type",0);
    defaults_.setMaxInt("svm:svr:kernel_type",3);


    defaults_.setValue("svm:svc:degree", 3, "For POLY");
    defaults_.setMinInt("svm:svc:degree",1);

    defaults_.setValue("svm:svr:degree", 3, "For POLY");
    defaults_.setMinInt("svm:svr:degree",1);

    defaults_.setValue("svm:svc:gamma", 0.0, "For POLY/RBF/SIGMOID");
    defaults_.setMinFloat("svm:svc:gamma",0.0);

    defaults_.setValue("svm:svr:gamma", 0.0, "For POLY/RBF/SIGMOID");
    defaults_.setMinFloat("svm:svr:gamma",0.0);

    defaults_.setValue("svm:svc:coef0", 0.0, "For POLY/SIGMOID");
    defaults_.setMinFloat("svm:svc:coef0",0.0);

    defaults_.setValue("svm:svr:coef0", 0.0, "For POLY/SIGMOID");
    defaults_.setMinFloat("svm:svr:coef0",0.0);

    defaults_.setValue("svm:svc:eps", 0.001, "Stopping criterion");
    defaults_.setValue("svm:svr:eps", 0.001, "Stopping criterion");

    defaults_.setValue("svm:svc:C", 1.0, "Cost of constraint violation");
    defaults_.setValue("svm:svr:C", 1.0, "Cost of constraint violation");

    defaults_.setValue("svm:svr:p", 0.1, "The epsilon for the loss function in epsilon-SVR");

    defaults_.setValue("svm:svc:nu", 0.5, "For NU_SVC, ONE_CLASS and NU_SVR");
    defaults_.setValue("svm:svr:nu", 0.5, "For NU_SVC, ONE_CLASS and NU_SVR");

    defaults_.setValue("svm:cache_size", 100, "Size of kernel cache in MB");
    defaults_.setMinInt("svm:cache_size",1);

    defaults_.setValue("svm:shrinking", "true", "Perform shrinking");
    defaults_.setValidStrings("svm:shrinking", StringList::create("true,false"));

    defaults_.setValue("svm:n_fold", 5, "n_fold cross validation is performed");
    defaults_.setMinInt("svm:n_fold",1);

    defaults_.setValue("svm:scaling", "true", "Apply scaling of feature values");
    defaults_.setValidStrings("svm:scaling", StringList::create("true,false"));

    defaults_.setValue("svm:scaling_lower", 0.0, "Lower bound for scaling");
    defaults_.setValue("svm:scaling_upper", 1.0, "Upper bound for scaling");

    defaults_.setValue("svm:svc:balancing", "true", "Use class balanced SVC training");
    defaults_.setValidStrings("svm:svc:balancing", StringList::create("true,false"));

    defaults_.setSectionDescription("svm", "Parameters controlling SVM trainig behaviour. All parameter names are chosen as in the libSVM library. Please refer to libSVM documentation for explanation");
    defaults_.setSectionDescription("svm:svc", "Parameters for svm - classification of missing/abundant");
    defaults_.setSectionDescription("svm:svr", "Parameters for svm - regression of peak intensities");


    defaultsToParam_();
  }

  SvmTheoreticalSpectrumGeneratorTrainer::SvmTheoreticalSpectrumGeneratorTrainer(const SvmTheoreticalSpectrumGeneratorTrainer& rhs) :
    DefaultParamHandler(rhs)
  {
    updateMembers_();
  }

  SvmTheoreticalSpectrumGeneratorTrainer& SvmTheoreticalSpectrumGeneratorTrainer::operator =(const SvmTheoreticalSpectrumGeneratorTrainer& rhs)
  {
    if (this != &rhs)
    {
      DefaultParamHandler::operator=(rhs);      
      updateMembers_();
    }
    return *this;
  }

  SvmTheoreticalSpectrumGeneratorTrainer::~SvmTheoreticalSpectrumGeneratorTrainer()
  {
  }


  void SvmTheoreticalSpectrumGeneratorTrainer::trainModel(const PeakMap &spectra, const std::vector<AASequence>&annotations, String filename, Size precursor_charge)
  {
    //----------- BEGIN OF PARAMETER READING-------------------------

    //upper and lower bounds for rescaling of features
    double upper = 0.0;
    double lower = 0.0;

    bool secondary_types = false;


    Size number_of_intensity_levels = (UInt) param_.getValue("number_intensity_levels");
    Size number_of_regions = (UInt) param_.getValue("number_regions");
    DoubleReal parent_tolerance = param_.getValue("parent_tolerance");
    DoubleReal peak_tolerance = param_.getValue("peak_tolerance");
    bool write_outfiles = param_.getValue("write_training_files").toBool();

    //build set of ion types
    std::vector<IonType> ion_types;
    std::vector<bool>is_primary;

    if (DataValue(param_.getValue("add_b_ions")).toBool())
    {
      ion_types.push_back(IonType(Residue::BIon, EmpiricalFormula(), 1));
      is_primary.push_back(true);
    }
    if (DataValue(param_.getValue("add_y_ions")).toBool())
    {
      ion_types.push_back(IonType(Residue::YIon, EmpiricalFormula(), 1));
      is_primary.push_back(true);
    }    
    if (DataValue(param_.getValue("add_x_ions")).toBool())
    {
      ion_types.push_back(IonType(Residue::XIon, EmpiricalFormula(), 1));
      is_primary.push_back(false);
      secondary_types=true;
    }
    if (DataValue(param_.getValue("add_a_ions")).toBool())
    {
      ion_types.push_back(IonType(Residue::AIon, EmpiricalFormula(), 1));
      is_primary.push_back(false);
      secondary_types=true;
    }
    if (DataValue(param_.getValue("add_z_ions")).toBool())
    {
      ion_types.push_back(IonType(Residue::ZIon, EmpiricalFormula(), 1));
      is_primary.push_back(false);
      secondary_types=true;
    }
    if (DataValue(param_.getValue("add_c_ions")).toBool())
    {
      ion_types.push_back(IonType(Residue::CIon, EmpiricalFormula(), 1));
      is_primary.push_back(false);
      secondary_types=true;
    }

    if (DataValue(param_.getValue("add_losses")).toBool())
    {
      EmpiricalFormula loss_ammonia("NH3");
      EmpiricalFormula loss_water("H2O");

      if (DataValue(param_.getValue("add_b_ions")).toBool())
      {
        ion_types.push_back(IonType(Residue::BIon, loss_ammonia, 1));
        ion_types.push_back(IonType(Residue::BIon, loss_water, 1));
        is_primary.push_back(false);
        is_primary.push_back(false);
        secondary_types=true;
      }
      if (DataValue(param_.getValue("add_y_ions")).toBool())
      {
        ion_types.push_back(IonType(Residue::YIon, loss_ammonia, 1));
        ion_types.push_back(IonType(Residue::YIon, loss_water, 1));
        is_primary.push_back(false);
        is_primary.push_back(false);
        secondary_types=true;
      }
    }

    if (DataValue(param_.getValue("add_y2_ions")).toBool())
    {
      ion_types.push_back(IonType(Residue::YIon, EmpiricalFormula(), 2));
      is_primary.push_back(false);
      secondary_types=true;
    }
    if (DataValue(param_.getValue("add_b2_ions")).toBool())
    {
      ion_types.push_back(IonType(Residue::BIon, EmpiricalFormula(), 2));
      is_primary.push_back(false);
      secondary_types=true;
    }

    //-----------------------LOADING THE SVM PARAMETERS------------------


    int n_fold = param_.getValue("svm:n_fold");

    bool scaling = param_.getValue("svm:scaling").toBool();
    if (scaling)
    {
      upper=param_.getValue("svm:scaling_upper");
      lower=param_.getValue("svm:scaling_lower");
    }


    //-----------------------loading svr parameters
    svm_parameter svm_param_reg;
    svm_param_reg.cache_size = param_.getValue("svm:cache_size");
    svm_param_reg.shrinking = param_.getValue("svm:shrinking").toBool();
    svm_param_reg.nr_weight = 0;

    svm_param_reg.svm_type = 3 + (Int)param_.getValue("svm:svr_type");
    svm_param_reg.kernel_type = param_.getValue("svm:svr:kernel_type");
    svm_param_reg.degree = param_.getValue("svm:svr:degree");
    svm_param_reg.coef0 = param_.getValue("svm:svr:coef0");

    svm_param_reg.eps = param_.getValue("svm:svr:eps");
    svm_param_reg.C = param_.getValue("svm:svr:C");
    svm_param_reg.gamma = param_.getValue("svm:svr:gamma");
    svm_param_reg.nu = param_.getValue("svm:svr:nu");
    svm_param_reg.p=param_.getValue("svm:svr:p");
    svm_param_reg.probability = 0;

    //-----------------------loading svc parameters
    svm_parameter svm_param_class;
    svm_param_class.cache_size = param_.getValue("svm:cache_size");
    svm_param_class.shrinking = param_.getValue("svm:shrinking").toBool();
    svm_param_class.nr_weight = 0;

    svm_param_class.svm_type = (Int)param_.getValue("svm:svc_type");
    svm_param_class.kernel_type = param_.getValue("svm:svc:kernel_type");
    svm_param_class.degree = param_.getValue("svm:svc:degree");
    svm_param_class.coef0 = param_.getValue("svm:svc:coef0");

    svm_param_class.eps = param_.getValue("svm:svc:eps");
    svm_param_class.C = param_.getValue("svm:svc:C");
    svm_param_class.gamma = param_.getValue("svm:svc:gamma");
    svm_param_class.nu = param_.getValue("svm:svc:nu");
    svm_param_class.probability = 0;

    bool balancing = param_.getValue("svm:svc:balancing").toBool();


    //----------- END OF PARAMETER READING-------------------------


    String info_outfile_name = filename + ".info";
    TextFile info_outfile;
    info_outfile.push_back("<PrecursorCharge>");
    info_outfile.push_back(precursor_charge);
    info_outfile.push_back("</PrecursorCharge>");

    info_outfile.push_back("<PrimaryTypes>");

    //count number of usable spectra
    Size usable_spectra=spectra.size();
    std::vector<bool>is_spectrum_usable(spectra.size(), true);

    for (Size index = 0; index < spectra.size(); ++index)
    {
      //test whether annotation for spectrum match. If theoretical and empirical mass differ too much skip this spectrum
      Int empirical_charge=spectra[index].getPrecursors()[0].getCharge();
      DoubleReal empirical_parent_mass=spectra[index].getPrecursors()[0].getMZ();
      DoubleReal theoretical_parent_mass=annotations[index].getMonoWeight(Residue::Full, empirical_charge)/empirical_charge;

      if (abs(empirical_parent_mass-theoretical_parent_mass) > parent_tolerance)
      {
        is_spectrum_usable[index] = false;
        --usable_spectra;
        std::cerr << "skipping spectrum " << index <<" due to parent mass missmatch"<<std::endl;
      }
      else if (precursor_charge != empirical_charge)
      {
        is_spectrum_usable[index] = false;
        --usable_spectra;
        std::cerr << "skipping spectrum " << index <<" due to wrong precursor charge"<<std::endl;
      }
    }

    //only required for generatrion of descriptors
    SvmTheoreticalSpectrumGenerator spec_gen;    

    //Amino acid sequences used for obtaining prefix and suffix mass    
    AASequence prefix, suffix;

    DescriptorSet tmp_desc;
    Size num_features = spec_gen.generateDescriptorSet_(annotations[0],0,ion_types[0], 1, tmp_desc);

    //use number of features to adapt gamma parameter if set to 0 (same is used in svm-train.c)
    if(svm_param_class.gamma == 0)
    {
      svm_param_class.gamma = 1.0/num_features;
    }
    if(svm_param_reg.gamma == 0)
    {
      svm_param_reg.gamma = 1.0/num_features;
    }

    //vectors to store the minimum and maximum value appearing in the training data for each feature (required for scaling)
    ObservedIntensMap observed_intensities;
    std::map<Size, std::vector<DescriptorSet> >training_input;
    std::map<Size, std::vector<double> >training_output;

    Size spec_index = 0;
    Size x = 1;
    //run over all input spectra
    for (PeakMap::const_iterator map_it = spectra.begin(); map_it < spectra.end(); map_it+=x, spec_index+=x)
    {      
      if (!is_spectrum_usable[spec_index])
      {
        continue;
      }      

      Size precursor_charge = map_it->getPrecursors()[0].getCharge();
      PeakSpectrum input_spec_norm(*map_it);

      normalizeIntensity(input_spec_norm);

      for (Size type_nr = 0; type_nr < ion_types.size(); ++type_nr)
      {        
        //store the intensities of the detected peaks for the given type in the actual spectrum (required for sec. type prob model)
        countIntensities_(input_spec_norm, annotations[spec_index], ion_types[type_nr], observed_intensities, peak_tolerance, number_of_regions);
      }

      for (Size type_nr = 0; type_nr < ion_types.size(); ++type_nr)
      {
        if(!is_primary[type_nr])
        {
          continue;
        }        

        DoubleReal true_offset_mass = 0.0;
        Residue::ResidueType residue = ion_types[type_nr].residue;
        Int charge = ion_types[type_nr].charge;
        EmpiricalFormula loss = ion_types[type_nr].loss;

        training_input[type_nr].reserve(usable_spectra);
        training_output[type_nr].reserve(usable_spectra);               

        for (Size frag_pos = 1; frag_pos < annotations[spec_index].size(); ++frag_pos)
        {
          prefix = annotations[spec_index].getPrefix(frag_pos);
          suffix = annotations[spec_index].getSuffix(annotations[spec_index].size() - frag_pos);

          //N-terminal fragments
          if (residue == Residue::AIon || residue == Residue::BIon || residue == Residue::CIon)
          {
            EmpiricalFormula loss_ion = prefix.getFormula(residue, charge) - loss;
            true_offset_mass = loss_ion.getMonoWeight() / charge;
          }
          //C-terminal fragments
          else if (residue == Residue::XIon || residue == Residue::YIon || residue == Residue::ZIon)
          {
            EmpiricalFormula loss_ion = suffix.getFormula(residue, charge) - loss;
            true_offset_mass = loss_ion.getMonoWeight() / charge;
          }

          //find the closest peak in the spectrum
          DoubleReal observed_peak_intensity = -1;
          Size true_nearest_peak_ind = input_spec_norm.findNearest(true_offset_mass);

          //check whether this peak is within the allowed mass range
          if (fabs(true_offset_mass - input_spec_norm[true_nearest_peak_ind].getMZ()) <= peak_tolerance)
          {
            observed_peak_intensity =  input_spec_norm[true_nearest_peak_ind].getIntensity();
          }

          DescriptorSet descriptors;
          spec_gen.generateDescriptorSet_(annotations[spec_index],frag_pos-1,ion_types[type_nr], precursor_charge, descriptors);

          training_input[type_nr].push_back(descriptors);
          training_output[type_nr].push_back(observed_peak_intensity);          
        }
      }//end of running over all types
    }//end of running over all spectra

    //scale the input data    
    std::vector<double>min_features;
    min_features.reserve(num_features);
    std::vector<double>max_features;
    max_features.reserve(num_features);

    min_features.assign(num_features, std::numeric_limits<double>::infinity());
    max_features.assign(num_features, -1 * std::numeric_limits<double>::infinity());

    //SCALING
    if(scaling)
    {
      for (Size type_nr = 0; type_nr < ion_types.size(); ++type_nr)
      {
        //first compute the minimal and maximal entries for each feature
        std::vector<svm_node>::iterator it;
        for(Size train_row=0; train_row<training_input[type_nr].size(); ++train_row)
        {
          Size index = 0;
          for(it = training_input[type_nr][train_row].descriptors.begin(); it != training_input[type_nr][train_row].descriptors.end(); ++it)
          {
            while( (Int) index < it->index-1)
            {
              min_features[index] = std::min(min_features[index], 0.0);
              max_features[index] = std::max(max_features[index],0.0);
              ++index;
            }
            min_features[index] = std::min(min_features[index], it->value);
            max_features[index] = std::max(max_features[index], it->value);
            ++index;
          }

          while(index<num_features)
          {
            min_features[index] = std::min(min_features[index], 0.0);
            max_features[index] = std::max(max_features[index],0.0);
            ++index;
          }
        }
      }

      //now use the min-max entries to scale the features
      //we have to set the private members feature_min and feature_max to call the scaleDescriptorSet method
      spec_gen.mp_.feature_min = min_features;
      spec_gen.mp_.feature_max = max_features;

      for (Size type_nr = 0; type_nr < ion_types.size(); ++type_nr)
      {
        if(!is_primary[type_nr])
        {
          continue;
        }

        //with the minimal and maximal now scale
        for(Size train_row=0; train_row<training_input[type_nr].size(); ++train_row)
        {
          spec_gen.scaleDescriptorSet_(training_input[type_nr][train_row],lower, upper);
        }        
      }
    }


    // Now 2 Steps - first we train a SVM classifier to predict abundance or miss of a peak.
    // In second step we train a SVM regression model to predict intensities. This svm is trained
    // only with rows for abundant peaks    

    //Within this loop the SVR and SVC are trained for the primary ion types
    for (Size type_nr = 0; type_nr < ion_types.size(); ++type_nr)
    {
      if(!is_primary[type_nr])
      {
        continue;
      }

      //------------------------------------------------------------------------------------------
      //----------------------------------Training of SVR-Model-----------------------------
      //------------------------------------------------------------------------------------------

      std::vector<DescriptorSet>training_input_reg;
      std::vector<double>training_output_reg;
      training_input_reg.reserve(training_input[type_nr].size());
      training_output_reg.reserve(training_output[type_nr].size());

      for(Size i = 0; i < training_output[type_nr].size(); ++i)
      {
        training_input_reg.push_back(training_input[type_nr][i]);
        training_output_reg.push_back(std::max(0.0,training_output[type_nr][i]));
      }

      String svm_model_file_reg;
      if(write_outfiles)
      {        
        write_training_file_(training_input_reg, training_output_reg, String("Training_") + Residue::getResidueTypeName(ion_types[type_nr].residue) + "_" +
                             ion_types[type_nr].loss.getString() + "_" + ion_types[type_nr].charge + "_reg.dat");
      }
      else
      {
        std::vector<double> predictions_reg(training_input_reg.size(), 0);
        svm_node ** input_training_reg = new svm_node * [training_input_reg.size()];
        double * output_training_reg = &training_output_reg[0];
        for(Size i = 0; i < training_input_reg.size(); ++i)
        {
          input_training_reg[i] = &(training_input_reg[i].descriptors[0]);
        }

        svm_problem svm_p_reg;
        svm_p_reg.l = (int) training_input_reg.size();
        svm_p_reg.y = output_training_reg;
        svm_p_reg.x = input_training_reg;        

        if(n_fold > 1)
        {
          //perform cross validation
          std::cerr<<"start cross validation for SVR model"<<std::endl;
          svm_cross_validation(&svm_p_reg, &svm_param_reg, n_fold, &predictions_reg[0]);

          //compute cross correlation coefficient of predicted and true intensities
          double coeff_reg =  gsl_stats_correlation(&training_output_reg[0], 1, &predictions_reg[0], 1, training_output_reg.size());
          std::cerr<<"Correlation Coefficient: "<<coeff_reg<<std::endl;
        }

        //train and store
        svm_model *model_reg =  svm_train(&svm_p_reg, &svm_param_reg);

        svm_model_file_reg = filename + "_" + Residue::getResidueTypeName(ion_types[type_nr].residue) + "_" +
                             ion_types[type_nr].loss.getString() + "_" + ion_types[type_nr].charge + "_reg.svm";

        if(svm_save_model(svm_model_file_reg.c_str(), model_reg)==-1)
        {
          throw Exception::FileNotWritable(__FILE__, __LINE__, __PRETTY_FUNCTION__,svm_model_file_reg);
        }

        svm_destroy_model(model_reg);
        delete [] input_training_reg;

      }//end of else


      //------------------------------------------------------------------------------------------
      //----------------------------------Training of SVC Model-----------------------------
      //------------------------------------------------------------------------------------------

      //make binary
      for(Size i=0; i<training_output[type_nr].size(); ++i)
      {
        training_output[type_nr][i] = training_output[type_nr][i] != -1 ? 1:0 ; //if peak was abundant class 1 else 0
      }

      //create balanced set      
      if(balancing)
      {
        std::vector<std::vector<Size> >training_data_by_class(2);
        std::vector<double>::const_iterator it;
        std::vector<DescriptorSet>tmp_training;
        std::vector<double>tmp_output;
        for(it = training_output[type_nr].begin(); it != training_output[type_nr].end(); ++it)
        {
          training_data_by_class[*it].push_back(it-training_output[type_nr].begin());
        }
        Size min_size=std::min(training_data_by_class[0].size(),training_data_by_class[1].size());
        std::random_shuffle(training_data_by_class[0].begin(),training_data_by_class[0].end() );
        std::random_shuffle(training_data_by_class[1].begin(),training_data_by_class[1].end() );

        for(Size num = 0; num < min_size; ++num)
        {
          for(Size intens = 0; intens < 2; ++intens)
          {
            tmp_training.push_back(training_input[type_nr][training_data_by_class[intens][num]]);
            tmp_output.push_back(intens);
          }
        }
        training_input[type_nr] = tmp_training;
        training_output[type_nr] = tmp_output;
      }

      if(write_outfiles)
      {
        write_training_file_(training_input[type_nr], training_output[type_nr], String("Training_") + Residue::getResidueTypeName(ion_types[type_nr].residue) + "_" +
                             ion_types[type_nr].loss.getString() + "_" + ion_types[type_nr].charge + "_class.dat");
      }

      else
      {        
        double * output_training_class = &training_output[type_nr][0];
        svm_node ** input_training_class = new svm_node * [training_input[type_nr].size()];
        for(Size i = 0; i < training_input[type_nr].size(); ++i)
        {
          input_training_class[i] = &(training_input[type_nr][i].descriptors[0]);
        }

        svm_problem svm_p_class;
        svm_p_class.l = (int) training_output[type_nr].size();
        svm_p_class.y = output_training_class;
        svm_p_class.x = input_training_class;

        //perform cross validation
        if(n_fold > 1)
        {
          std::cerr<<"start cross validation for SVC model"<<std::endl;
          std::vector<double> predictions_class(training_input[type_nr].size(), 0.);
          svm_cross_validation(&svm_p_class, &svm_param_class, n_fold, &predictions_class[0]);

          Size TP = 0, FP = 0, TN = 0, FN = 0;
          for(Size i=0; i<predictions_class.size(); ++i)
          {
            if(predictions_class[i] == 0 && output_training_class[i] == 0) ++TN;
            if(predictions_class[i] == 1 && output_training_class[i] == 1) ++TP;
            if(predictions_class[i] == 1 && output_training_class[i] == 0) ++FP;
            if(predictions_class[i] == 0 && output_training_class[i] == 1) ++FN;
          }
          std::cerr<<"TP: "<<TP<<"  TN: "<<TN<<"  FP: "<<FP<<" FN: "<<FN<<std::endl;
        }

        //train and store
        svm_model *model_class =  svm_train(&svm_p_class, &svm_param_class);

        String svm_model_file_class = filename + "_" + Residue::getResidueTypeName(ion_types[type_nr].residue) + "_" +
                                      ion_types[type_nr].loss.getString() + "_" + ion_types[type_nr].charge + "_class.svm";

        if(svm_save_model(svm_model_file_class.c_str(), model_class) == -1)
        {
          throw Exception::FileNotWritable(__FILE__, __LINE__, __PRETTY_FUNCTION__,svm_model_file_class);
        }

        svm_destroy_model(model_class);
        delete [] input_training_class;


        //add entries to info file
        info_outfile.push_back("<IonType>");
        info_outfile.push_back(ion_types[type_nr].residue);
        info_outfile.push_back(ion_types[type_nr].loss.getString());
        info_outfile.push_back(ion_types[type_nr].charge);
        info_outfile.push_back("<SvmModelFileClass>");
        info_outfile.push_back(svm_model_file_class);
        info_outfile.push_back("</SvmModelFileClass>");
        info_outfile.push_back("<SvmModelFileReg>");
        info_outfile.push_back(svm_model_file_reg);
        info_outfile.push_back("</SvmModelFileReg>");
        info_outfile.push_back("</IonType>");

      }//end of else

    }//End of training SVR and SVC for the primary types


    //If we are only generating the outfiles we can terminare here
    if(write_outfiles)
    {
      return;
    }

    info_outfile.push_back("</PrimaryTypes>");
    info_outfile.push_back("<ScalingUpper>");
    info_outfile.push_back(upper);
    info_outfile.push_back("</ScalingUpper>");
    info_outfile.push_back("<ScalingLower>");
    info_outfile.push_back(lower);
    info_outfile.push_back("</ScalingLower>");
    info_outfile.push_back("<MaxFeatures>");
    info_outfile.insert(info_outfile.end(), max_features.begin(), max_features.end());
    info_outfile.push_back("</MaxFeatures>");
    info_outfile.push_back("<MinFeatures>");
    info_outfile.insert(info_outfile.end(), min_features.begin(), min_features.end());
    info_outfile.push_back("</MinFeatures>");    

    //------------------------------------------------------------------------------------------
    //----------------------Training prob. model for secondary types------------------
    //------------------------------------------------------------------------------------------

    info_outfile.push_back("<SecondaryTypes>");
    if(secondary_types)
    {
      info_outfile.push_back("<IntensityLevels>");
      info_outfile.push_back(number_of_intensity_levels);
      info_outfile.push_back("</IntensityLevels>");

      trainSecondaryTypes_(info_outfile, number_of_regions, number_of_intensity_levels, observed_intensities, ion_types, is_primary);
    }

    info_outfile.push_back("</SecondaryTypes>");
    info_outfile.store(info_outfile_name);
  }


  //like in Cong Zhou paper
  void SvmTheoreticalSpectrumGeneratorTrainer::normalizeIntensity(PeakSpectrum &S) const
  {
    NLargest n_larg;
    Param larg_param=n_larg.getParameters();
    larg_param.setValue("n", (Int)(S.size()*0.8));
    n_larg.setParameters(larg_param);
    n_larg.filterPeakSpectrum(S);
    S.sortByPosition();

    Normalizer norm;
    Param norm_param=norm.getParameters();
    norm_param.setValue("method", "to_TIC");
    norm.setParameters(norm_param);
    norm.filterPeakSpectrum(S);

    DoubleReal min_intens=  std::numeric_limits<double>::infinity();
    DoubleReal max_intens = -1 * std::numeric_limits<double>::infinity();

    std::vector<DoubleReal>intensities(S.size());
    for(Size i=0; i<S.size(); ++i)
    {
      //std::cerr<<"before norm: "<<S[i].getIntensity()<<std::endl;
      if(S[i].getIntensity()>0)
      {
        intensities[i]=log(S[i].getIntensity() * 100);
        max_intens=std::max(max_intens, intensities[i]);
        min_intens=std::min(min_intens, intensities[i]);
      }
    }

    DoubleReal lower = 0;
    DoubleReal upper = 1;
    //normalize intensities to one
    for(Size i=0; i<S.size(); ++i)
    {
      if(S[i].getIntensity()>0)
      {
        DoubleReal intens = lower + (upper-lower) *
                             (intensities[i]-min_intens)/
                             (max_intens-min_intens);
        S[i].setIntensity(intens);
      }
      else
      {
        S[i].setIntensity(0);
      }
      //std::cerr<<"normed intens: "<<S[i].getIntensity()<<std::endl;
    }
  }


  void SvmTheoreticalSpectrumGeneratorTrainer::trainSecondaryTypes_(TextFile &info_outfile,
                                                                    Size number_of_regions,
                                                                    Size number_of_intensity_levels,
                                                                    ObservedIntensMap &observed_intensities,
                                                                    const std::vector<IonType> &ion_types,
                                                                    const std::vector<bool> &is_primary
                                                                    )
  {
    std::vector<DoubleReal>tmp;
    for(Size region = 0; region < number_of_regions; ++ region)
    {
    //we start by binning the intensities. We select the bin boarders such that
    //the intensities of the primary ions are equally split
      std::vector<DoubleReal> &observed_b=observed_intensities[std::make_pair(IonType(Residue::BIon), region)];
      tmp.insert(tmp.end(), observed_b.begin(), observed_b.end());
      std::vector<DoubleReal> &observed_y=observed_intensities[std::make_pair(IonType(Residue::YIon), region)];
      tmp.insert(tmp.end(),observed_y.begin(), observed_y.end());
    }

    std::vector<DoubleReal>bin_boarders(number_of_intensity_levels);
    std::vector<DoubleReal>bin_values(number_of_intensity_levels);
    std::sort(tmp.begin(), tmp.end());

    //find the first nonzero
    Size first_non_zero=0;
    while(first_non_zero<tmp.size()&& tmp[first_non_zero]==0)
    {
      ++first_non_zero;
    }

    Size non_zero_size = tmp.size()-first_non_zero;
    Size prev_index = 0;
    for(Size i=1; i<number_of_intensity_levels; ++i)
    {
      Size index = i* (DoubleReal) non_zero_size/(number_of_intensity_levels-1) + first_non_zero;
      DoubleReal count = 0;
      for(Size j = prev_index; j < index; ++j)
      {
        count += tmp[j];
      }
      bin_boarders[i-1] = tmp[index-1];
      bin_values[i-1] = count/(index-prev_index);
      prev_index = index;
    }

    info_outfile.push_back("<IntensityBinBoarders>");
    for(Size i = 0; i < number_of_intensity_levels-1; ++i)
    {
      info_outfile.push_back(bin_boarders[i]);
    }
    info_outfile.push_back("</IntensityBinBoarders>");
    info_outfile.push_back("<IntensityBinValues>");
    for(Size i = 0; i < number_of_intensity_levels-1; ++i)
    {
      info_outfile.push_back(bin_values[i]);
    }
    info_outfile.push_back("</IntensityBinValues>");

    //use the boarder values to bin the entries
    for(Size region = 0; region < number_of_regions; ++ region)
    {
      for(Size i=0; i<ion_types.size(); ++i)
      {
        std::vector<DoubleReal> & intensities = observed_intensities[std::make_pair(ion_types[i], region)];
        for(Size j=0; j<intensities.size(); ++j)
        {
          DoubleReal intens = intensities[j];
          if(intens == 0.0|| intens == -1.0)
            continue;

          Size k=1;
          while(k<number_of_intensity_levels-1 && intens>bin_boarders[k-1] )
            ++k;

          intensities[j]=k;
        }
      }
    }

    std::map<std::pair<IonType,Size>, std::vector<std::vector<DoubleReal> > >joint_counts;
    std::map<std::pair<IonType,Size>,  std::vector<DoubleReal> >background_counts;

    //count joint appearences of primary and secondary peaks
    for(Size i=0; i<ion_types.size(); ++i)
    {
      const IonType & type = ion_types[i];
      if(is_primary[i])
        continue;

      IonType primary_type=IonType(Residue::BIon);

      if(type.residue == Residue::YIon ||
         type.residue == Residue::XIon ||
         type.residue == Residue::ZIon)
      {
        primary_type=IonType(Residue::YIon);
      }
      for(Size region =0; region < number_of_regions; ++region)
      {
        joint_counts[std::make_pair(type, region)].assign(number_of_intensity_levels, std::vector<DoubleReal>(number_of_intensity_levels,1));
        background_counts[std::make_pair(type, region)].assign(number_of_intensity_levels, number_of_intensity_levels);
        const std::vector<DoubleReal> &secondary = observed_intensities[std::make_pair(type, region)];
        const std::vector<DoubleReal> &primary = observed_intensities[std::make_pair(primary_type, region)];

        for(Size j=0; j<primary.size(); ++j)
        {
          if(secondary[j]!=-1.0)
          {
            ++joint_counts[std::make_pair(type, region)][(Size)secondary[j]][(Size)primary[j]];
            ++background_counts[std::make_pair(type, region)][(Size)primary[j]];
          }
        }
      }
    }

    //compute conditional probabilities and store them in the outfile
    for(Size i=0; i<ion_types.size(); ++i)
    {
      const IonType & type = ion_types[i];
      if(is_primary[i])
        continue;

      info_outfile.push_back("<IonType>");
      info_outfile.push_back(type.residue);
      info_outfile.push_back(type.loss.getString());
      info_outfile.push_back(type.charge);
      info_outfile.push_back("<ConditionalProbabilities>");

      for(Size region = 0; region < number_of_regions; ++region)
      {
        info_outfile.push_back("<Region "+ String(region) + ">");
        std::vector<DoubleReal> & back_counts = background_counts[std::make_pair(type, region)];

        for(Size prim = 0; prim < number_of_intensity_levels; ++prim)
        {
          for(Size sec = 0; sec < number_of_intensity_levels; ++sec)
          {
            if(back_counts[prim]!=0)
            {
              joint_counts[std::make_pair(type, region)][sec][prim] = joint_counts[std::make_pair(type, region)][sec][prim]/back_counts[prim];
              //std::cerr<<"conditional prob  "<<type.residue<<" "<<sec<<"  "<<prim<<"  "<<joint_counts[std::make_pair(type, region)][sec][prim]<<std::endl;
            }
            info_outfile.push_back(joint_counts[std::make_pair(type, region)][sec][prim]);
          }
        }
        info_outfile.push_back("</Region " + String(region) + ">");
      }
      info_outfile.push_back("</ConditionalProbabilities>");
      info_outfile.push_back("</IonType>");
      }
  }


  void SvmTheoreticalSpectrumGeneratorTrainer::countIntensities_(const PeakSpectrum &spectrum,
                                           const AASequence &annotation,
                                           IonType type,
                                           std::map<std::pair<IonType, Size>, std::vector<DoubleReal> > & observed_intensities,
                                           DoubleReal tolerance,
                                           Size number_of_regions
                                           )
  {
    Residue::ResidueType residue = type.residue;
    EmpiricalFormula loss = type.loss;
    Int charge = type.charge;

    std::set<String> possible_n_term_losses;
    std::set<String> possible_c_term_losses;
    DoubleReal true_offset_mass = -1.;


    for (Size frag_pos = 1; frag_pos < annotation.size(); ++frag_pos)
    {

      AASequence prefix = annotation.getPrefix(frag_pos);
      AASequence suffix = annotation.getSuffix(annotation.size() - frag_pos);

      Size region = std::min(number_of_regions-1, (Size)floor(number_of_regions * prefix.getMonoWeight(Residue::Internal)/annotation.getMonoWeight()));

      if (annotation[frag_pos-1].hasNeutralLoss())
      {
        std::vector<EmpiricalFormula> loss_formulas = annotation[frag_pos-1].getLossFormulas();
        for (Size k = 0; k != loss_formulas.size(); ++k)
        {
          possible_n_term_losses.insert(loss_formulas[k].getString());
        }
      }
      //check for possible losses on the c-terminal ions
      possible_c_term_losses.clear();
      for(Size pos = frag_pos; pos < annotation.size(); ++pos)
      {
        if (annotation[pos].hasNeutralLoss())
        {
          std::vector<EmpiricalFormula> loss_formulas = annotation[pos].getLossFormulas();
          for (Size k = 0; k != loss_formulas.size(); ++k)
          {
            possible_c_term_losses.insert(loss_formulas[k].getString());
          }
        }
      }    

      //N-terminal fragments
      if (residue == Residue::AIon || residue == Residue::BIon || residue == Residue::CIon)
      {
        //if loss is not supported or no loss ions shall be generated -- continue
        if(!loss.isEmpty() && (!possible_n_term_losses.count(loss.getString())))
        {
          observed_intensities[std::make_pair(type, region)].push_back(-1);
          continue;
        }
        EmpiricalFormula loss_ion = prefix.getFormula(residue, charge) - loss;
        true_offset_mass = loss_ion.getMonoWeight() / charge;
      }
      //C-terminal fragments
      else if (residue == Residue::XIon || residue == Residue::YIon || residue == Residue::ZIon)
      {
        //if loss is not supported or no loss ions shall be generated -- continue
        if(!loss.isEmpty() && (!possible_c_term_losses.count(loss.getString())))
        {
          observed_intensities[std::make_pair(type, region)].push_back(-1);
          continue;
        }
        EmpiricalFormula loss_ion = suffix.getFormula(residue, charge) - loss;
        true_offset_mass = loss_ion.getMonoWeight() / charge;
      }

      //find the closest peak in the spectrum
      DoubleReal observed_peak_intensity = 0;
      Size true_nearest_peak_ind = spectrum.findNearest(true_offset_mass);

      //check whether this peak is within the allowed mass range
      if (fabs(true_offset_mass - spectrum[true_nearest_peak_ind].getMZ()) <= tolerance)
      {
        observed_peak_intensity = spectrum[true_nearest_peak_ind].getIntensity();
      }
      observed_intensities[std::make_pair(type, region)].push_back(observed_peak_intensity);
    }
  }


  void SvmTheoreticalSpectrumGeneratorTrainer::write_training_file_(std::vector<DescriptorSet> &training_input, std::vector<DoubleReal> &training_output, String filename)
  {
    std::cerr<<"Creating Training File.. "<< filename;
    TextFile file;
    for(Size i = 0; i < training_input.size(); ++i)
    {
      std::stringstream ss;
      ss<<training_output[i]<<" ";
      std::vector<svm_node>::iterator it_debug;
      for(it_debug=training_input[i].descriptors.begin(); it_debug < training_input[i].descriptors.end() - 1; ++it_debug)
      {
        ss<<" "<< it_debug->index << ":" << it_debug->value;
      }
      file.push_back (ss.str());
    }
    file.store(filename);
    std::cerr<<" Done" <<std::endl;
  }


}//Namespace




