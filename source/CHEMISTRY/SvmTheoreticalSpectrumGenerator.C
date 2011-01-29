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

#include <OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/CONCEPT/Constants.h>

#define DEBUG

namespace OpenMS
{
  String SvmTheoreticalSpectrumGenerator::ResidueTypeToString_(Residue::ResidueType type)
  {
    switch(type)
    {
      case(Residue::AIon):return("AIon");
      break;
      case(Residue::BIon):return("BIon");
      break;
      case(Residue::CIon):return("CIon");
      break;
      case(Residue::XIon):return("XIon");
      break;
      case(Residue::YIon):return("YIon");
      break;
      case(Residue::ZIon):return("ZIon");
      break;
      default: return("undefined ion type");
    }
  }


  SvmTheoreticalSpectrumGenerator::SvmTheoreticalSpectrumGenerator() :
    TheoreticalSpectrumGenerator()    
  {
    this->setName("SVMTheoreticalSpectrumGenerator");

//    defaults_.setValue("add_isotopes", "true", "If set to 1 isotope peaks of the product ion peaks are added");
//    defaults_.setValue("max_isotope", 2, "Defines the maximal isotopic peak which is added, add_isotopes must be set to 1");
//    defaults_.setValue("add_metainfo", "false", "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");
//    defaults_.setValue("add_losses", "false", "Adds common losses to those ion expect to have them, only water and ammonia loss is considered");
//    defaults_.setValue("add_precursor_peaks", "false", "Adds peaks of the precursor to the spectrum, which happen to occur sometimes");
    defaults_.setValue("model_file_name", "examples/simulation/SvmMSim.model", "Name of the probabilistic Model file");
    defaults_.setValue("add_first_nterminals", "false", "Whether first A, B and C ion shall be generated");
//    defaults_.setValue("precursor_intensity", 1.0, "Intensity of the precursor peak");
//    defaults_.setValue("precursor_H2O_intensity", 1.0, "Intensity of the H2O loss peak of the precursor");
//    defaults_.setValue("precursor_NH3_intensity", 1.0, "Intensity of the NH3 loss peak of the precursor");

    defaultsToParam_();

    // just in case someone wants the ion names;
    p_.metaRegistry().registerName("IonName", "Name of the ion");    
  }

  SvmTheoreticalSpectrumGenerator::SvmTheoreticalSpectrumGenerator(const SvmTheoreticalSpectrumGenerator& rhs) :
    TheoreticalSpectrumGenerator(rhs),
    class_models_(rhs.class_models_),
    reg_models_(rhs.reg_models_),
    ion_types_(rhs.ion_types_),
    number_intensity_levels_(rhs.number_intensity_levels_),
    feature_max_(rhs.feature_max_),
    feature_min_(rhs.feature_min_),
    scaling_lower_(rhs.scaling_lower_),
    scaling_upper_(rhs.scaling_upper_)
  {
    updateMembers_();
  }

  SvmTheoreticalSpectrumGenerator& SvmTheoreticalSpectrumGenerator::operator =(const SvmTheoreticalSpectrumGenerator& rhs)
  {
    if (this != &rhs)
    {
      TheoreticalSpectrumGenerator::operator=(rhs);
      class_models_ = rhs.class_models_;
      reg_models_ = rhs.reg_models_;
      number_intensity_levels_ = rhs.number_intensity_levels_;
      feature_max_ = rhs.feature_max_;
      feature_min_ = rhs.feature_min_;
      scaling_lower_ = rhs.scaling_lower_;
      scaling_upper_= rhs.scaling_upper_;
      updateMembers_();
    }
    return *this;
  }

  SvmTheoreticalSpectrumGenerator::~SvmTheoreticalSpectrumGenerator()
  {
    for(Size i=0; i<class_models_.size(); ++i)
    {
      svm_destroy_model(class_models_[i]);
    }
    for(Size i=0; i<reg_models_.size(); ++i)
    {
      svm_destroy_model(reg_models_[i]);
    }
  }


  Size SvmTheoreticalSpectrumGenerator::generateDescriptorSet2_
		(AASequence peptide, Size position, IonType type, Size /* precursor_charge */, DescriptorSet& desc_set)
  {
    ///map AA to integers
    static std::map<String,Size>aa_to_index;
    static std::map<String,DoubleReal>hydrophobicity, helicity, basicity;

    Int index = 0;
    //if the map was not yet generated (this is first call of this Function)
    if (aa_to_index.empty())
    {
      ResidueDB* res_db;
      res_db = ResidueDB::getInstance();
      std::set<const Residue*>all_aa = res_db->getResidues("Natural20");
      std::set<const Residue*>::const_iterator aa_it;
      for (aa_it=all_aa.begin(); aa_it!=all_aa.end(); ++aa_it)
      {
        aa_to_index[(*aa_it)->getOneLetterCode()] = index;
        ++index;
      }
    }

    if (hydrophobicity.empty())
    {
      hydrophobicity["A"] = 0.16;
      hydrophobicity["C"] = 2.5;
      hydrophobicity["D"] = -2.49;
      hydrophobicity["E"] = -1.5;
      hydrophobicity["F"] = 5;
      hydrophobicity["G"] = -3.31;
      hydrophobicity["H"] = -4.63;
      hydrophobicity["I"] = 4.76;
      hydrophobicity["K"] = -5;
      hydrophobicity["L"] = 4.76;
      hydrophobicity["M"] = 3.23;
      hydrophobicity["N"] = -3.79;
      hydrophobicity["P"] = -4.92;
      hydrophobicity["Q"] = -2.76;
      hydrophobicity["R"] = -2.77;
      hydrophobicity["S"] = -2.85;
      hydrophobicity["T"] = -1.08;
      hydrophobicity["V"] = 3.02;
      hydrophobicity["W"] = 4.88;
      hydrophobicity["Y"] = 2;      
    }

    if (helicity.empty())
    {
      helicity["A"] = 1.24;
      helicity["C"] = 0.79;
      helicity["D"] = 0.89;
      helicity["E"] = 0.85;
      helicity["F"] = 1.26;
      helicity["G"] = 1.15;
      helicity["H"] = 0.97;
      helicity["I"] = 1.28;
      helicity["K"] = 0.88;
      helicity["L"] = 1.28;
      helicity["M"] = 1.22;
      helicity["N"] = 0.94;
      helicity["P"] = 0.57;
      helicity["Q"] = 0.96;
      helicity["R"] = 0.95;
      helicity["S"] = 1;
      helicity["T"] = 1.09;
      helicity["V"] = 1.27;
      helicity["W"] = 1.07;
      helicity["Y"] = 1.11;      
    }


    if (basicity.empty())
    {
      basicity["A"] = 206.4;
      basicity["C"] = 206.2;
      basicity["D"] = 208.6;
      basicity["E"] = 215.5;
      basicity["F"] = 212.1;
      basicity["G"] = 202.7;
      basicity["H"] = 223.7;
      basicity["I"] = 209.6;
      basicity["K"] = 221.8;
      basicity["L"] = 209.6;
      basicity["M"] = 213.3;
      basicity["N"] = 212.8;
      basicity["P"] = 214.4;
      basicity["Q"] = 214.2;
      basicity["R"] = 237.0;
      basicity["S"] = 207.6;
      basicity["T"] = 211.7;
      basicity["V"] = 208.7;
      basicity["W"] = 216.1;
      basicity["Y"] = 213.1;      
    }



    std::vector<svm_node>descriptors_tmp;
    descriptors_tmp.reserve(50);

    Size charge=type.charge;
    Residue::ResidueType res_type=type.residue;
    EmpiricalFormula  loss =type.loss;

    AASequence fragment;
    if(res_type==Residue::AIon || res_type==Residue::BIon || res_type==Residue::CIon)
    {
      fragment=peptide.getPrefix(position+1);
    }

    if(res_type==Residue::XIon || res_type==Residue::YIon || res_type==Residue::ZIon)
    {
      fragment=peptide.getSuffix(peptide.size()-(position+1));
    }

    DoubleReal fragment_mass = (fragment.getMonoWeight(res_type, charge) - loss.getMonoWeight());

    Residue res_n = peptide.getResidue(position);
    Residue res_c = peptide.getResidue(position + 1);

    String res_n_string=res_n.getOneLetterCode();
    String res_c_string=res_c.getOneLetterCode();;

    Size num_aa=aa_to_index.size();
    index=1;
    svm_node node;

    //RB_C
    node.index = index + (Int)aa_to_index[peptide.getResidue(position+1).getOneLetterCode()];
    node.value=1;
    descriptors_tmp.push_back(node);
    index+=num_aa;

    //RB_N
    node.index=index + (Int)aa_to_index[peptide.getResidue(position).getOneLetterCode()];
    node.value = 1;
    descriptors_tmp.push_back(node);
    index += (Int)num_aa;

    //DB_N
    node.index=index++;
    node.value=position+1;
    descriptors_tmp.push_back(node);

    //DB_C
    node.index=index++;
    node.value=peptide.size()-position-1;
    descriptors_tmp.push_back(node);

    //DB_M
    node.index=index++;
    node.value=fabs(1+position-peptide.size()/2.0);
    descriptors_tmp.push_back(node);

    //B_E
    node.index=index++;
    node.value=(position==0 || position==peptide.size()-2);
    descriptors_tmp.push_back(node);

    //BaRB_N
    node.index=index++;
    node.value=basicity[res_n_string];
    descriptors_tmp.push_back(node);

    //BaRB_C
    node.index=index++;
    node.value=basicity[res_c_string];
    descriptors_tmp.push_back(node);

    //BaRB_A
    node.index=index++;
    node.value=(basicity[res_n_string]+basicity[res_c_string])/2.0;
    descriptors_tmp.push_back(node);

    //BaRB_D
    node.index=index++;
    node.value=basicity[res_n_string]-basicity[res_c_string];    
    descriptors_tmp.push_back(node);

    DoubleReal ba_p=0, hy_p=0, ba_yi=0, hy_yi=0, ba_bi=0, hy_bi=0;
    for(Size i=0; i<peptide.size();++i)
    {
      ba_p+=basicity[peptide.getResidue(i).getOneLetterCode()];
      hy_p+=hydrophobicity[peptide.getResidue(i).getOneLetterCode()];
    }
    for(Size i=0; i<position+1;++i)
    {
      ba_bi+=basicity[peptide.getResidue(i).getOneLetterCode()];
      hy_bi+=hydrophobicity[peptide.getResidue(i).getOneLetterCode()];
    }
    for(Size i=position+1; i<peptide.size();++i)
    {
      ba_yi+=basicity[peptide.getResidue(i).getOneLetterCode()];
      hy_yi+=hydrophobicity[peptide.getResidue(i).getOneLetterCode()];
    }

    //BaYI
    node.index=index++;
    node.value=ba_yi;
    descriptors_tmp.push_back(node);

    //BaBI
    node.index=index++;
    node.value=ba_bi;
    descriptors_tmp.push_back(node);

    //BaP
    node.index=index++;
    node.value=ba_p;
    descriptors_tmp.push_back(node);

    //HeRB_N
    node.index=index++;
    node.value=helicity[res_n_string];    
    descriptors_tmp.push_back(node);

    //HeRB_C
    node.index=index++;
    node.value=helicity[res_c_string];
    descriptors_tmp.push_back(node);

    //HeRB_A
    node.index=index++;
    node.value= (helicity[res_n_string]+helicity[res_c_string])/2;
    descriptors_tmp.push_back(node);

    //HeRB_D
    node.index=index++;
    node.value=helicity[res_n_string]-helicity[res_c_string];
    descriptors_tmp.push_back(node);

    //HyRB_N
    node.index=index++;
    node.value=hydrophobicity[res_n_string];    
    descriptors_tmp.push_back(node);

    //HyRB_C
    node.index=index++;
    node.value=hydrophobicity[res_c_string];
    descriptors_tmp.push_back(node);

    //HyRB_A
    node.index=index++;
    node.value= (hydrophobicity[res_n_string]+hydrophobicity[res_c_string])/2;
    descriptors_tmp.push_back(node);

    //HyRB_D
    node.index=index++;
    node.value=hydrophobicity[res_n_string]-hydrophobicity[res_c_string];
    descriptors_tmp.push_back(node);

    //Hy_YI
    node.index=index++;
    node.value=hy_yi;
    descriptors_tmp.push_back(node);

    //Hy_BI
    node.index=index++;
    node.value=hy_bi;
    descriptors_tmp.push_back(node);

    //Hy_P
    node.index=index++;
    node.value=hy_p;
    descriptors_tmp.push_back(node);

    //PIRB_N
    node.index=index++;
    node.value=res_n.getPiValue();
    descriptors_tmp.push_back(node);

    //PIRB_C
    node.index=index++;
    node.value=res_c.getPiValue();
    descriptors_tmp.push_back(node);

    //PIRB_A
    node.index=index++;
    node.value=(res_n.getPiValue()+res_c.getPiValue())/2;
    descriptors_tmp.push_back(node);

    //PIRB_D
    node.index=index++;
    node.value=res_n.getPiValue()-res_c.getPiValue();
    descriptors_tmp.push_back(node);

    //LP
    node.index=index++;
    node.value=peptide.size();
    descriptors_tmp.push_back(node);

    //LYI (works for B also)
    node.index=index++;
    node.value=fragment.size();
    descriptors_tmp.push_back(node);

    //RLIP (works for B also)
    node.index=index++;
    node.value=(DoubleReal)fragment.size()/peptide.size();
    descriptors_tmp.push_back(node);

    //NBaR_P
    node.index=index++;
    node.value=peptide.getNumberOf("H")+peptide.getNumberOf("K")+peptide.getNumberOf("R");
    descriptors_tmp.push_back(node);

    //NBaR_YI (works for B also)
    node.index=index++;
    node.value=fragment.getNumberOf("H")+fragment.getNumberOf("K")+fragment.getNumberOf("R");
    descriptors_tmp.push_back(node);

    //MP
    node.index=index++;
    node.value=peptide.getMonoWeight(Residue::Full);
    descriptors_tmp.push_back(node);

    //MYI
    node.index=index++;
    node.value=fragment_mass;
    descriptors_tmp.push_back(node);

    //RMIP
    node.index=index++;
    node.value=fragment_mass / peptide.getMonoWeight(Residue::Full);
    descriptors_tmp.push_back(node);

    //DBBa
    Size bas_dist_left=position, bas_dist_right=position+1;

    while(bas_dist_left>0)
    {
      if(peptide.getResidue(bas_dist_left).getOneLetterCode()=="H" ||
         peptide.getResidue(bas_dist_left).getOneLetterCode()=="R" ||
         peptide.getResidue(bas_dist_left).getOneLetterCode()=="K")
        break;
      --bas_dist_left;
    }
    while(bas_dist_right<peptide.size())
    {
      if(peptide.getResidue(bas_dist_right).getOneLetterCode()=="H" ||
         peptide.getResidue(bas_dist_right).getOneLetterCode()=="R" ||
         peptide.getResidue(bas_dist_right).getOneLetterCode()=="K")
        break;
      ++bas_dist_right;
    }

    //DBBa
    node.index=index++;
    node.value=std::min(position-bas_dist_left, bas_dist_right-position-1);
    descriptors_tmp.push_back(node);


    node.index=-1;
    descriptors_tmp.push_back(node);

    desc_set.descriptors=descriptors_tmp;

    return index;
  }


  void SvmTheoreticalSpectrumGenerator::load()
  {
    String svm_info_file = param_.getValue("model_file_name");
    if (! File::readable( svm_info_file ) )
    { // look in OPENMS_DATA_PATH
      svm_info_file = File::find( svm_info_file );
    }

    //delete all previously loaded models
    for(Size i=0; i<class_models_.size(); ++i)
    {
      svm_destroy_model(class_models_[i]);
    }
    //delete all previously loaded models
    for(Size i=0; i<reg_models_.size(); ++i)
    {
      svm_destroy_model(reg_models_[i]);
    }

    class_models_.clear();
    reg_models_.clear();
    ion_types_.clear();
    secondary_types_.clear();
    conditional_prob_.clear();

    TextFile info_file;
    String path_to_models = File().path(svm_info_file)+"/";
    info_file.load(svm_info_file);    

    TextFile::iterator left_marker = info_file.search("<PrimaryTypes>");
    TextFile::iterator right_marker = info_file.search("</PrimaryTypes>");
    if(left_marker==right_marker)
    {
      //Todo throw different exception (File Corrupt)
      throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__,svm_info_file);
    }

    //Now read the primary types and load the corresponding svm models
    left_marker = info_file.search(left_marker, "<IonType>");
    right_marker = info_file.search(left_marker, "<ScalingUpper>");
    while(left_marker < right_marker)
    {
      ++left_marker;
      Residue::ResidueType res = (Residue::ResidueType)(*left_marker++).toInt();
      EmpiricalFormula loss(*(left_marker++));
      UInt charge(left_marker++->toInt());
      ion_types_.push_back(IonType(res,loss,charge));      

      left_marker = info_file.search(left_marker, "<SvmModelFileClass>");      
      String svm_filename(*(++left_marker));
      LOG_INFO<<"SVM model file loaded: "<<svm_filename<<std::endl;
      class_models_.push_back(svm_load_model((path_to_models+svm_filename).c_str()));
      if(class_models_.back()==0)
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__,svm_filename);
      }

      left_marker = info_file.search(left_marker, "<SvmModelFileReg>");
      svm_filename= *(++left_marker);
      LOG_INFO<<"SVM model file loaded: "<<svm_filename<<std::endl;
      reg_models_.push_back(svm_load_model((path_to_models+svm_filename).c_str()));
      if(reg_models_.back()==0)
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__,svm_filename);
      }            

      left_marker = info_file.search(left_marker, "<IonType>");      
    }

    //Read the upper and lower bounds for the feature vectors (used for scaling)
    left_marker = right_marker;
    scaling_upper_=(++left_marker)->toDouble();    

    left_marker = info_file.search(left_marker, "<ScalingLower>");
    scaling_lower_=(++left_marker)->toDouble();    

    left_marker = info_file.search(left_marker, "<MaxFeatures>");
    right_marker = info_file.search(left_marker, "</MaxFeatures>");
    while(++left_marker!=right_marker)
    {
      feature_max_.push_back(left_marker->toDouble());
    }
    left_marker = info_file.search(right_marker, "<MinFeatures>");
    right_marker = info_file.search(right_marker, "</MinFeatures>");
    while(++left_marker!=right_marker)
    {
      feature_min_.push_back(left_marker->toDouble());
    }

    std::vector<IonType>sec_ion_types;
    //Load the secondary types
    left_marker = info_file.search(left_marker, "<SecondaryTypes>");    

    left_marker = info_file.search(left_marker, "<IntensityLevels>");
    right_marker = info_file.search(left_marker, "</IntensityLevels>");
    if(left_marker!=right_marker)
    {
      number_intensity_levels_= (++left_marker)->toInt();
    }

    left_marker = info_file.search(left_marker, "<IntensityBinBoarders>");
    right_marker = info_file.search(left_marker, "</IntensityBinBoarders>");
    while(++left_marker!=right_marker)
    {
      intensity_bin_boarders_.push_back(left_marker->toDouble());
    }

    left_marker = info_file.search(left_marker, "<IntensityBinValues>");
    right_marker = info_file.search(left_marker, "</IntensityBinValues>");
    while (++left_marker != right_marker)
    {
      intensity_bin_values_.push_back(left_marker->toDouble());
    }

    left_marker = info_file.search(left_marker, "<IonType>");
    right_marker = info_file.search(left_marker, "</IonType>");
    while (left_marker != right_marker)
    {
      //load type information
      ++left_marker;
      Residue::ResidueType res = (Residue::ResidueType)(*left_marker++).toInt();
      EmpiricalFormula loss(*(left_marker++));
      UInt charge(left_marker++->toInt());
      IonType actual_type(res,loss,charge);
      sec_ion_types.push_back(actual_type);      

      //load conditional probabilities      
      left_marker = info_file.search(left_marker, "<ConditionalProbabilities>");      

      Size region = 0;
      left_marker = info_file.search(left_marker, "<Region "+String(region)+">");
      right_marker = info_file.search(left_marker, "</IonType>");      
      while(left_marker<right_marker)
      {
        conditional_prob_[std::make_pair(actual_type, region)].assign(number_intensity_levels_, std::vector<DoubleReal>(number_intensity_levels_));
        for(Size prim=0; prim<number_intensity_levels_; ++prim)
        {
          for(Size sec=0; sec<number_intensity_levels_; ++sec)
          {            
            conditional_prob_[std::make_pair(actual_type,region)][prim][sec]=(++left_marker)->toDouble();
          }
        }
        ++region;
        left_marker = info_file.search(left_marker, "<Region "+String(region)+">");        
      }
      number_regions_=region;

      left_marker = info_file.search(right_marker, "<IonType>");
      right_marker = info_file.search(left_marker, "</IonType>");
    }


    //map the secondary types to their corresponding primary types
    for(Size i=0; i<sec_ion_types.size(); ++i)
    {
      const IonType &tmp= sec_ion_types[i];
      if(tmp.residue==Residue::BIon || tmp.residue==Residue::AIon || tmp.residue==Residue::CIon)
      {
        secondary_types_[IonType(Residue::BIon)].push_back(tmp);
      }
      else if(tmp.residue==Residue::YIon || tmp.residue==Residue::XIon || tmp.residue==Residue::ZIon)
      {
        secondary_types_[IonType(Residue::YIon)].push_back(tmp);
      }
    }

    /*
      //DEBUG STUFF
    //TEST parse successful
    std::cerr<<"Primary Types:"<<std::endl;
    for(Size i=0; i<ion_types_.size(); ++i)
    {
      std::cout<<ion_types_[i].residue<<" "<<ion_types_[i].loss<<std::endl;
    }
    std::cerr<<"scaling_lower "<<scaling_lower_<<std::endl;
    std::cerr<<"scaling_upper "<<scaling_upper_<<std::endl;
    std::cerr<<"MaxFeatures "<<std::endl;
    for(Size i=0; i<feature_max_.size(); ++i)
    {
      std::cerr<<feature_max_[i]<<std::endl;
    }
    std::cerr<<"MinFeatures "<<std::endl;
    for(Size i=0; i<feature_min_.size(); ++i)
    {
      std::cerr<<feature_min_[i]<<std::endl;
    }
    std::cerr<<"Binboarders "<<std::endl;
    for(Size i=0; i<intensity_bin_boarders_.size(); ++i)
    {
      std::cerr<<intensity_bin_boarders_[i]<<std::endl;
    }
    std::cerr<<"Binvalues "<<std::endl;
    for(Size i=0; i<intensity_bin_values_.size(); ++i)
    {
      std::cerr<<intensity_bin_values_[i]<<std::endl;
    }

    std::cerr<<"Secondary Types:"<<std::endl;
    for(Size i=0; i<ion_types_.size(); ++i)
    {
      for(Size j=0; j<secondary_types_[ion_types_[i]].size(); ++j)
      {
        std::cout<<secondary_types_[ion_types_[i]][j].residue<<" "<<secondary_types_[ion_types_[i]][j].loss<<std::endl;
        std::cerr<<"Conditional Probabilities: "<<std::endl;
        for(Size reg=0; reg<number_regions_; ++reg)
        {
          for(Size prim=0; prim<number_intensity_levels_; ++prim)
          {
            for(Size sec=0; sec<number_intensity_levels_; ++sec)
            {
              std::cerr<<reg<<"   "<<sec<<" | "<<prim<<"   "<<conditional_prob_[std::make_pair(secondary_types_[ion_types_[i]][j], reg)][prim][sec]<<std::endl;
            }
          }
        }
      }
    }
    */
  }


  void SvmTheoreticalSpectrumGenerator::simulate(RichPeakSpectrum &spectrum, const AASequence &peptide, const gsl_rng *rng, Size precursor_charge)
  {
    if(class_models_.empty() || reg_models_.empty() || ion_types_.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, __PRETTY_FUNCTION__, "no svm models loaded. Call load function before using simulate");
    }


    //load parameters
    bool add_isotopes = param_.getValue("add_isotopes").toBool();    
    Size max_isotope = (Size)param_.getValue("max_isotope");
    bool add_losses = param_.getValue("add_losses").toBool();    
    bool add_first_nterminals = param_.getValue("add_first_nterminals").toBool();
    bool add_metainfo = param_.getValue("add_metainfo").toBool();
    // bool add_precursor_peaks = param_.getValue("add_precursor_peaks").toBool();

    std::set<String>possible_n_term_losses;
    std::set<String>possible_c_term_losses;

    AASequence prefix, suffix;

    const AASequence * ion = 0;

    UInt ion_nr=0;
    gsl_ran_discrete_t * gsl_gen=0;


    for (Size i = 1; i < peptide.size(); ++i)
    {
      prefix=peptide.getPrefix(i);
      suffix=peptide.getSuffix(peptide.size()-i);

      Size region = std::min(number_regions_-1, (Size)floor(number_regions_ * prefix.getMonoWeight(Residue::Internal)/peptide.getMonoWeight()));

      //check for possible losses on the n-terminal ions
      if (peptide[i-1].hasNeutralLoss())
      {
        std::vector<EmpiricalFormula> loss_formulas = peptide[i-1].getLossFormulas();
        for (Size k = 0; k != loss_formulas.size(); ++k)
        {
          possible_n_term_losses.insert(loss_formulas[k].getString());
        }
      }

      //check for possible losses on the c-terminal ions
      possible_c_term_losses.clear();
      for(Size pos=i;pos<peptide.size();++pos)
      {
        if (peptide[pos].hasNeutralLoss())
        {
          std::vector<EmpiricalFormula> loss_formulas = peptide[pos].getLossFormulas();
          for (Size k = 0; k != loss_formulas.size(); ++k)
          {
            possible_c_term_losses.insert(loss_formulas[k].getString());
          }
        }
      }

      std::vector<std::pair<IonType, DoubleReal> >peaks_to_generate;

      for (Size type_nr=0; type_nr<ion_types_.size(); ++type_nr)
      {        
        Residue::ResidueType residue = ion_types_[type_nr].residue;
        EmpiricalFormula loss_formula = ion_types_[type_nr].loss;
        // Int charge = ion_types_[type_nr].charge;

        //std::cerr<<"prim type: "<<residue<<" "<<loss_formula<<"  "<<charge<<std::endl;

        //determine whether type is N- or C-terminal
        if (residue == Residue::AIon || residue == Residue::BIon || residue == Residue::CIon)
        {
          //usually no b1, a1, c1 ion added
          if((i < 2) && !add_first_nterminals)
          {
            continue;
          }          
          //if loss is not supported or no loss ions shall be generated -- continue
          if (!loss_formula.isEmpty() && (!possible_n_term_losses.count(loss_formula.getString()) || !add_losses))
          {
            continue;
          }          
          ion = &prefix;
          ion_nr = (UInt)i;
        }
        //same for c-terminal ions
        else if (residue == Residue::XIon || residue == Residue::YIon || residue == Residue::ZIon)
        {
          if (!loss_formula.isEmpty() && (!possible_c_term_losses.count(loss_formula.getString()) || !add_losses))
          {
            continue;
          }
          ion = &suffix;
          ion_nr = (UInt)(peptide.size()-i);
        }
        else
        {
          LOG_ERROR<< "requested unsupported ion type" << std::endl;
        }

        DescriptorSet descriptor;
        generateDescriptorSet2_(peptide,i-1,ion_types_[type_nr], precursor_charge, descriptor);
        if (scaling_lower_!=scaling_upper_)
        {
          scaleDescriptorSet_(descriptor, scaling_lower_, scaling_upper_);
        }

        DoubleReal predicted_intensity = 0.0;
        Size predicted_class = svm_predict(class_models_[type_nr], &descriptor.descriptors[0]);
        //predicted_class = 1;

        if (predicted_class==1)
        {
          predicted_intensity = std::max(0.0, svm_predict(reg_models_[type_nr], &descriptor.descriptors[0]));
          peaks_to_generate.push_back(std::make_pair(ion_types_[type_nr],predicted_intensity));          
        }

        //binning the predicted intensity
        Size bin=0;
        if (predicted_intensity > 0)
        {
          bin = 1;
          while ((bin < number_intensity_levels_ - 1) && predicted_intensity>intensity_bin_boarders_[bin-1])
					{
            ++bin;
					}
        }

        for (std::vector<IonType>::iterator it = secondary_types_[ion_types_[type_nr]].begin(); it!= secondary_types_[ion_types_[type_nr]].end(); ++it)
        {
          if (!add_losses && !it->loss.isEmpty())
          {
            continue;
          }

          //sample intensities for secondary types
          DoubleReal * condit_probs = &conditional_prob_[std::make_pair(*it, region)][bin][0];
          gsl_gen = gsl_ran_discrete_preproc(number_intensity_levels_, condit_probs);
          Size binned_int = gsl_ran_discrete(rng, gsl_gen);
          gsl_ran_discrete_free(gsl_gen);

          if(binned_int!=0)
          {
            peaks_to_generate.push_back(std::make_pair(*it,intensity_bin_values_[binned_int-1]));
          }
        }
      }


      for (Size i=0; i<peaks_to_generate.size(); ++i)
      {
        IonType type = peaks_to_generate[i].first;
        DoubleReal intensity = peaks_to_generate[i].second;
        Residue::ResidueType residue = type.residue;
        EmpiricalFormula loss_formula = type.loss;
        Int charge = type.charge;                

        //std::cerr<<"add: "<<residue<<"  "<<loss_formula<<"  "<<charge<<std::endl;
        EmpiricalFormula ion_formula = ion->getFormula(residue, charge) - loss_formula;
        DoubleReal mz_pos = ion_formula.getMonoWeight() / charge;

        String ion_name = ResidueTypeToString_(residue) +" "+ loss_formula.getString() +" "+ String(ion_nr)+String(charge, '+');

        if (add_isotopes)
        {          
          IsotopeDistribution dist = ion_formula.getIsotopeDistribution((Int)max_isotope);
          Size j = 0;
          for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it, ++j)
          {            
            p_.setMZ(mz_pos + (DoubleReal)j * Constants::NEUTRON_MASS_U / charge);
            p_.setIntensity(intensity * it->second);
            if (add_metainfo && j == 0)
            {
              p_.setMetaValue("IonName", ion_name);
            }
            spectrum.push_back(p_);            
          }
        }
        else
        {
          p_.setMZ(mz_pos);
          p_.setIntensity(intensity);
          if (add_metainfo)
          {
            p_.setMetaValue("IonName", ion_name);
          }
          spectrum.push_back(p_);
        }
      }
    }
    spectrum.sortByPosition();
  }

  void SvmTheoreticalSpectrumGenerator::scaleDescriptorSet_(DescriptorSet &desc, double lower, double upper)
  {
    std::vector<svm_node>tmp_desc;
    std::vector<svm_node>::iterator it;
    Size feature_index=1;
    Size num_features=feature_max_.size();    
    for(it=desc.descriptors.begin(); it!=(desc.descriptors.end()-1); ++it)
    {      
      while(feature_index<(Size)it->index)
      {
        double tmp_value = 0;
        scaleSingleFeature_(tmp_value, lower, upper, feature_min_[feature_index-1], feature_max_[feature_index-1]);
        if (tmp_value!=0)
        {
          svm_node tmp_node = {(Int)feature_index, tmp_value};
          tmp_desc.push_back(tmp_node);
        }
        ++feature_index;
      }      
      scaleSingleFeature_(it->value, lower, upper, feature_min_[feature_index-1], feature_max_[feature_index-1]);
      if(it->value!=0)
      {
        tmp_desc.push_back(*it);
      }
      ++feature_index;
    }

    while(feature_index<=num_features)
    {
      double tmp_value = 0;            
      scaleSingleFeature_(tmp_value, lower, upper, feature_min_[feature_index-1], feature_max_[feature_index-1]);
      if (tmp_value != 0)
      {
        svm_node tmp_node={(Int)feature_index, tmp_value};
        tmp_desc.push_back(tmp_node);
      }
      ++feature_index;
    }
    svm_node tmp_node={-1, -1};
    tmp_desc.push_back(tmp_node);

    desc.descriptors=tmp_desc;
  }

}//namespace




