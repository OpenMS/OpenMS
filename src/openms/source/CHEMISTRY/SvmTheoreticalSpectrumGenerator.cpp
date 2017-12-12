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

#include <OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <algorithm>
#include <iterator>

#include <boost/bind.hpp>
#include <boost/random/discrete_distribution.hpp>

#ifdef _OPENMP
#include <omp.h>
#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>
#include <boost/shared_ptr.hpp>
#endif

// #define DEBUG

namespace OpenMS
{

  std::map<String, Size> SvmTheoreticalSpectrumGenerator::aa_to_index_;
  std::map<String, double> SvmTheoreticalSpectrumGenerator::hydrophobicity_;
  std::map<String, double> SvmTheoreticalSpectrumGenerator::helicity_;
  std::map<String, double> SvmTheoreticalSpectrumGenerator::basicity_;

  bool SvmTheoreticalSpectrumGenerator::initializedMaps_ = false;

  String SvmTheoreticalSpectrumGenerator::ResidueTypeToString_(Residue::ResidueType type)
  {
    switch (type)
    {
    case (Residue::AIon): return "AIon";

    case (Residue::BIon): return "BIon";

    case (Residue::CIon): return "CIon";

    case (Residue::XIon): return "XIon";

    case (Residue::YIon): return "YIon";

    case (Residue::ZIon): return "ZIon";

    default: return "undefined ion type";
    }
  }

  void SvmTheoreticalSpectrumGenerator::initializeMaps_()
  {
    initializedMaps_ = true;

    if (aa_to_index_.empty())
    {
      ResidueDB* res_db;
      res_db = ResidueDB::getInstance();
      std::set<const Residue*> all_aa = res_db->getResidues("Natural20");
      std::set<String> residues;
      std::set<const Residue*>::const_iterator aa_it;
      for (aa_it = all_aa.begin(); aa_it != all_aa.end(); ++aa_it)
      {
        residues.insert((*aa_it)->getOneLetterCode());
      }
      std::set<String>::const_iterator aa2_it;
      Int index = 0;

      for (aa2_it = residues.begin(); aa2_it != residues.end(); ++aa2_it)
      {
        aa_to_index_[*aa2_it] = index;
        ++index;
      }
    }

    hydrophobicity_["A"] = 0.16;
    hydrophobicity_["C"] = 2.5;
    hydrophobicity_["D"] = -2.49;
    hydrophobicity_["E"] = -1.5;
    hydrophobicity_["F"] = 5;
    hydrophobicity_["G"] = -3.31;
    hydrophobicity_["H"] = -4.63;
    hydrophobicity_["I"] = 4.76;
    hydrophobicity_["K"] = -5;
    hydrophobicity_["L"] = 4.76;
    hydrophobicity_["M"] = 3.23;
    hydrophobicity_["N"] = -3.79;
    hydrophobicity_["P"] = -4.92;
    hydrophobicity_["Q"] = -2.76;
    hydrophobicity_["R"] = -2.77;
    hydrophobicity_["S"] = -2.85;
    hydrophobicity_["T"] = -1.08;
    hydrophobicity_["V"] = 3.02;
    hydrophobicity_["W"] = 4.88;
    hydrophobicity_["Y"] = 2;


    helicity_["A"] = 1.24;
    helicity_["C"] = 0.79;
    helicity_["D"] = 0.89;
    helicity_["E"] = 0.85;
    helicity_["F"] = 1.26;
    helicity_["G"] = 1.15;
    helicity_["H"] = 0.97;
    helicity_["I"] = 1.28;
    helicity_["K"] = 0.88;
    helicity_["L"] = 1.28;
    helicity_["M"] = 1.22;
    helicity_["N"] = 0.94;
    helicity_["P"] = 0.57;
    helicity_["Q"] = 0.96;
    helicity_["R"] = 0.95;
    helicity_["S"] = 1;
    helicity_["T"] = 1.09;
    helicity_["V"] = 1.27;
    helicity_["W"] = 1.07;
    helicity_["Y"] = 1.11;


    basicity_["A"] = 206.4;
    basicity_["C"] = 206.2;
    basicity_["D"] = 208.6;
    basicity_["E"] = 215.5;
    basicity_["F"] = 212.1;
    basicity_["G"] = 202.7;
    basicity_["H"] = 223.7;
    basicity_["I"] = 209.6;
    basicity_["K"] = 221.8;
    basicity_["L"] = 209.6;
    basicity_["M"] = 213.3;
    basicity_["N"] = 212.8;
    basicity_["P"] = 214.4;
    basicity_["Q"] = 214.2;
    basicity_["R"] = 237.0;
    basicity_["S"] = 207.6;
    basicity_["T"] = 211.7;
    basicity_["V"] = 208.7;
    basicity_["W"] = 216.1;
    basicity_["Y"] = 213.1;
  }

  SvmTheoreticalSpectrumGenerator::SvmTheoreticalSpectrumGenerator() :
    DefaultParamHandler("SvmTheoreticalSpectrumGenerator")
  {
    if (!initializedMaps_)
      initializeMaps_();

    defaults_.setValue("svm_mode", 1, "whether to predict abundant/missing using SVC (0) or predict intensities using SVR (1)");
    defaults_.setValue("model_file_name", "examples/simulation/SvmMSim.model", "Name of the probabilistic Model file");

    defaults_.setValue("add_isotopes", "false", "If set to 1 isotope peaks of the product ion peaks are added");
    defaults_.setValidStrings("add_isotopes", ListUtils::create<String>("true,false"));

    defaults_.setValue("max_isotope", 2, "Defines the maximal isotopic peak which is added, add_isotopes must be set to 1");

    defaults_.setValue("add_metainfo", "false", "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");
    defaults_.setValidStrings("add_metainfo", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_first_prefix_ion", "false", "If set to true e.g. b1 ions are added");
    defaults_.setValidStrings("add_first_prefix_ion", ListUtils::create<String>("true,false"));

    defaults_.setValue("hide_y_ions", "false", "Add peaks of y-ions to the spectrum");
    defaults_.setValidStrings("hide_y_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("hide_y2_ions", "false", "Add peaks of y-ions to the spectrum");
    defaults_.setValidStrings("hide_y2_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("hide_b_ions", "false", "Add peaks of b-ions to the spectrum");
    defaults_.setValidStrings("hide_b_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("hide_b2_ions", "false", "Add peaks of b-ions to the spectrum");
    defaults_.setValidStrings("hide_b2_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("hide_a_ions", "false", "Add peaks of a-ions to the spectrum");
    defaults_.setValidStrings("hide_a_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("hide_c_ions", "false", "Add peaks of c-ions to the spectrum");
    defaults_.setValidStrings("hide_c_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("hide_x_ions", "false", "Add peaks of  x-ions to the spectrum");
    defaults_.setValidStrings("hide_x_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("hide_z_ions", "false", "Add peaks of z-ions to the spectrum");
    defaults_.setValidStrings("hide_z_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("hide_losses", "false", "Adds common losses to those ion expect to have them, only water and ammonia loss is considered");
    defaults_.setValidStrings("hide_losses", ListUtils::create<String>("true,false"));

    // intensity options of the ions
    defaults_.setValue("y_intensity", 1.0, "Intensity of the y-ions");
    defaults_.setValue("b_intensity", 1.0, "Intensity of the b-ions");
    defaults_.setValue("a_intensity", 1.0, "Intensity of the a-ions");
    defaults_.setValue("c_intensity", 1.0, "Intensity of the c-ions");
    defaults_.setValue("x_intensity", 1.0, "Intensity of the x-ions");
    defaults_.setValue("z_intensity", 1.0, "Intensity of the z-ions");

    defaults_.setValue("relative_loss_intensity", 0.1, "Intensity of loss ions, in relation to the intact ion intensity");

    defaultsToParam_();
  }

  SvmTheoreticalSpectrumGenerator::SvmTheoreticalSpectrumGenerator(const SvmTheoreticalSpectrumGenerator& rhs) :
    DefaultParamHandler(rhs),
    mp_(rhs.mp_)
  {
    updateMembers_();
  }

  SvmTheoreticalSpectrumGenerator& SvmTheoreticalSpectrumGenerator::operator=(const SvmTheoreticalSpectrumGenerator& rhs)
  {
    if (this != &rhs)
    {
      DefaultParamHandler::operator=(rhs);
      mp_ = rhs.mp_;
      updateMembers_();
    }
    return *this;
  }

  SvmTheoreticalSpectrumGenerator::~SvmTheoreticalSpectrumGenerator()
  {
  }

  Size SvmTheoreticalSpectrumGenerator::generateDescriptorSet_(AASequence peptide, Size position, IonType type, Size /* precursor_charge */, DescriptorSet& desc_set)
  {

    std::vector<svm_node> descriptors_tmp;
    descriptors_tmp.reserve(50);

    Int charge = type.charge;
    Residue::ResidueType res_type = type.residue;
    EmpiricalFormula  loss = type.loss;

    AASequence fragment;
    if (res_type == Residue::AIon || res_type == Residue::BIon || res_type == Residue::CIon)
    {
      fragment = peptide.getPrefix(position + 1);
    }

    if (res_type == Residue::XIon || res_type == Residue::YIon || res_type == Residue::ZIon)
    {
      fragment = peptide.getSuffix(peptide.size() - (position + 1));
    }

    double fragment_mass = (fragment.getMonoWeight(res_type, charge) - loss.getMonoWeight());

    Residue res_n = peptide.getResidue(position);
    Residue res_c = peptide.getResidue(position + 1);

    String res_n_string = res_n.getOneLetterCode();
    String res_c_string = res_c.getOneLetterCode();

    Size num_aa = aa_to_index_.size();
    Int index = 1;
    svm_node node;

    //RB_C
    node.index = index + (Int)aa_to_index_[peptide.getResidue(position + 1).getOneLetterCode()];
    node.value = 1;
    descriptors_tmp.push_back(node);
    index += (Int)num_aa;

    //RB_N
    node.index = index + (Int)aa_to_index_[peptide.getResidue(position).getOneLetterCode()];
    node.value = 1;
    descriptors_tmp.push_back(node);
    index += (Int)num_aa;

    //DB_N
    node.index = index++;
    node.value = position + 1;
    descriptors_tmp.push_back(node);

    //DB_C
    node.index = index++;
    node.value = peptide.size() - position - 1;
    descriptors_tmp.push_back(node);

    //DB_M
    node.index = index++;
    node.value = fabs(1 + position - peptide.size() / 2.0);
    descriptors_tmp.push_back(node);

    //B_E
    node.index = index++;
    node.value = (position == 0 || position == peptide.size() - 2);
    descriptors_tmp.push_back(node);

    //BaRB_N
    node.index = index++;
    node.value = basicity_[res_n_string];
    descriptors_tmp.push_back(node);

    //BaRB_C
    node.index = index++;
    node.value = basicity_[res_c_string];
    descriptors_tmp.push_back(node);

    //BaRB_A
    node.index = index++;
    node.value = (basicity_[res_n_string] + basicity_[res_c_string]) / 2.0;
    descriptors_tmp.push_back(node);

    //BaRB_D
    node.index = index++;
    node.value = basicity_[res_n_string] - basicity_[res_c_string];
    descriptors_tmp.push_back(node);

    double ba_p = 0, hy_p = 0, ba_yi = 0, hy_yi = 0, ba_bi = 0, hy_bi = 0;
    for (Size i = 0; i < peptide.size(); ++i)
    {
      ba_p += basicity_[peptide.getResidue(i).getOneLetterCode()];
      hy_p += hydrophobicity_[peptide.getResidue(i).getOneLetterCode()];
    }
    for (Size i = 0; i < position + 1; ++i)
    {
      ba_bi += basicity_[peptide.getResidue(i).getOneLetterCode()];
      hy_bi += hydrophobicity_[peptide.getResidue(i).getOneLetterCode()];
    }
    for (Size i = position + 1; i < peptide.size(); ++i)
    {
      ba_yi += basicity_[peptide.getResidue(i).getOneLetterCode()];
      hy_yi += hydrophobicity_[peptide.getResidue(i).getOneLetterCode()];
    }

    //BaYI
    node.index = index++;
    node.value = ba_yi;
    descriptors_tmp.push_back(node);

    //BaBI
    node.index = index++;
    node.value = ba_bi;
    descriptors_tmp.push_back(node);

    //BaP
    node.index = index++;
    node.value = ba_p;
    descriptors_tmp.push_back(node);

    //HeRB_N
    node.index = index++;
    node.value = helicity_[res_n_string];
    descriptors_tmp.push_back(node);

    //HeRB_C
    node.index = index++;
    node.value = helicity_[res_c_string];
    descriptors_tmp.push_back(node);

    //HeRB_A
    node.index = index++;
    node.value = (helicity_[res_n_string] + helicity_[res_c_string]) / 2;
    descriptors_tmp.push_back(node);

    //HeRB_D
    node.index = index++;
    node.value = helicity_[res_n_string] - helicity_[res_c_string];
    descriptors_tmp.push_back(node);

    //HyRB_N
    node.index = index++;
    node.value = hydrophobicity_[res_n_string];
    descriptors_tmp.push_back(node);

    //HyRB_C
    node.index = index++;
    node.value = hydrophobicity_[res_c_string];
    descriptors_tmp.push_back(node);

    //HyRB_A
    node.index = index++;
    node.value = (hydrophobicity_[res_n_string] + hydrophobicity_[res_c_string]) / 2;
    descriptors_tmp.push_back(node);

    //HyRB_D
    node.index = index++;
    node.value = hydrophobicity_[res_n_string] - hydrophobicity_[res_c_string];
    descriptors_tmp.push_back(node);

    //Hy_YI
    node.index = index++;
    node.value = hy_yi;
    descriptors_tmp.push_back(node);

    //Hy_BI
    node.index = index++;
    node.value = hy_bi;
    descriptors_tmp.push_back(node);

    //Hy_P
    node.index = index++;
    node.value = hy_p;
    descriptors_tmp.push_back(node);

    //PIRB_N
    node.index = index++;
    node.value = res_n.getPiValue();
    descriptors_tmp.push_back(node);

    //PIRB_C
    node.index = index++;
    node.value = res_c.getPiValue();
    descriptors_tmp.push_back(node);

    //PIRB_A
    node.index = index++;
    node.value = (res_n.getPiValue() + res_c.getPiValue()) / 2;
    descriptors_tmp.push_back(node);

    //PIRB_D
    node.index = index++;
    node.value = res_n.getPiValue() - res_c.getPiValue();
    descriptors_tmp.push_back(node);

    //LP
    node.index = index++;
    node.value = peptide.size();
    descriptors_tmp.push_back(node);

    //LYI (works for B also)
    node.index = index++;
    node.value = fragment.size();
    descriptors_tmp.push_back(node);

    //RLIP (works for B also)
    node.index = index++;
    node.value = (double)fragment.size() / peptide.size();
    descriptors_tmp.push_back(node);

    //NBaR_P
    String umps = peptide.toUnmodifiedString(); // unmodified peptide string
    node.index = index++;
    node.value = std::count(umps.begin(), umps.end(), 'H') + std::count(umps.begin(), umps.end(), 'K') + std::count(umps.begin(), umps.end(), 'R');
    descriptors_tmp.push_back(node);

    //NBaR_YI (works for B also)
    node.index = index++;
    String umfs = fragment.toUnmodifiedString(); // unmodified fragment string
    node.value = std::count(umfs.begin(), umfs.end(), 'H') + std::count(umfs.begin(), umfs.end(), 'K') + std::count(umfs.begin(), umfs.end(), 'R');
    descriptors_tmp.push_back(node);

    //MP
    node.index = index++;
    node.value = peptide.getMonoWeight(Residue::Full);
    descriptors_tmp.push_back(node);

    //MYI
    node.index = index++;
    node.value = fragment_mass;
    descriptors_tmp.push_back(node);

    //RMIP
    node.index = index++;
    node.value = fragment_mass / peptide.getMonoWeight(Residue::Full);
    descriptors_tmp.push_back(node);

    //DBBa
    Size bas_dist_left = position, bas_dist_right = position + 1;

    while (bas_dist_left > 0)
    {
      if (peptide.getResidue(bas_dist_left).getOneLetterCode() == "H" ||
          peptide.getResidue(bas_dist_left).getOneLetterCode() == "R" ||
          peptide.getResidue(bas_dist_left).getOneLetterCode() == "K")
        break;
      --bas_dist_left;
    }
    while (bas_dist_right < peptide.size())
    {
      if (peptide.getResidue(bas_dist_right).getOneLetterCode() == "H" ||
          peptide.getResidue(bas_dist_right).getOneLetterCode() == "R" ||
          peptide.getResidue(bas_dist_right).getOneLetterCode() == "K")
        break;
      ++bas_dist_right;
    }

    //DBBa
    node.index = index++;
    node.value = std::min(position - bas_dist_left, bas_dist_right - position - 1);
    descriptors_tmp.push_back(node);


    node.index = -1;
    descriptors_tmp.push_back(node);

    desc_set.descriptors = descriptors_tmp;

    return index;
  }

  void SvmTheoreticalSpectrumGenerator::load()
  {
    String svm_info_file = param_.getValue("model_file_name");
    if (!File::readable(svm_info_file)) // look in OPENMS_DATA_PATH
    {
      svm_info_file = File::find(svm_info_file);
    }

    //delete all previously loaded models
    //for(Size i = 0; i < mp_.class_models.size(); ++i)
    //{
    //  mp_.class_models[i];
    //svm_destroy_model(mp_.class_models[i]);
    //}
    //delete all previously loaded models
    //for(Size i = 0; i < mp_.reg_models.size(); ++i)
    //{
    //  delete mp_.class_models[i];
    //svm_destroy_model(mp_.reg_models[i]);
    //}

    mp_.class_models.clear();
    mp_.reg_models.clear();
    mp_.ion_types.clear();
    mp_.secondary_types.clear();
    mp_.conditional_prob.clear();
    mp_.feature_max.clear();
    mp_.feature_min.clear();
    mp_.intensity_bin_boarders.clear();
    mp_.intensity_bin_values.clear();
    mp_.static_intensities.clear();

    TextFile info_file;
    String path_to_models = File().path(svm_info_file) + "/";
    info_file.load(svm_info_file);

    TextFile::ConstIterator left_marker = StringListUtils::searchPrefix(info_file.begin(), info_file.end(), "<PrecursorCharge>");
    TextFile::ConstIterator right_marker = StringListUtils::searchPrefix(info_file.begin(), info_file.end(), "</PrecursorCharge>");
    if (left_marker == right_marker)
    {
      //Todo throw different exception (File Corrupt)
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, svm_info_file);
    }
    ++left_marker;
    precursor_charge_ = left_marker->toInt();

    left_marker = StringListUtils::searchPrefix(info_file.begin(), info_file.end(), "<PrimaryTypes>");
    right_marker = StringListUtils::searchPrefix(info_file.begin(), info_file.end(), "</PrimaryTypes>");
    if (left_marker == right_marker)
    {
      //Todo throw different exception (File Corrupt)
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, svm_info_file);
    }

    //Now read the primary types and load the corresponding svm models
    left_marker = StringListUtils::searchPrefix(left_marker, info_file.end(), "<IonType>");
    right_marker = StringListUtils::searchPrefix(left_marker, info_file.end(), "<ScalingUpper>");
    while (left_marker < right_marker)
    {
      ++left_marker;
      Residue::ResidueType res = (Residue::ResidueType)(*left_marker++).toInt();
      EmpiricalFormula loss(*(left_marker++));
      UInt charge = (left_marker++)->toInt();
      mp_.ion_types.push_back(IonType(res, loss, charge));

      left_marker = StringListUtils::searchPrefix(left_marker, info_file.end(), "<SvmModelFileClass>");
      String svm_filename(*(++left_marker));

      boost::shared_ptr<SVMWrapper> sh_ptr_c(new SVMWrapper);
      sh_ptr_c.get()->loadModel((path_to_models + svm_filename).c_str());

      mp_.class_models.push_back(sh_ptr_c);
      LOG_INFO << "SVM model file loaded: " << svm_filename << std::endl;


      left_marker = StringListUtils::searchPrefix(left_marker, info_file.end(), "<SvmModelFileReg>");
      svm_filename = *(++left_marker);

      boost::shared_ptr<SVMWrapper> sh_ptr_r(new SVMWrapper);
      sh_ptr_r.get()->loadModel((path_to_models + svm_filename).c_str());

      mp_.reg_models.push_back(sh_ptr_r);
      LOG_INFO << "SVM model file loaded: " << svm_filename << std::endl;

      left_marker = StringListUtils::searchPrefix(left_marker, info_file.end(), "<IonType>");
    }

    //Read the upper and lower bounds for the feature vectors (used for scaling)
    left_marker = right_marker;
    mp_.scaling_upper = (++left_marker)->toDouble();

    left_marker = StringListUtils::searchPrefix(left_marker, info_file.end(), "<ScalingLower>");
    mp_.scaling_lower = (++left_marker)->toDouble();

    left_marker = StringListUtils::searchPrefix(left_marker, info_file.end(), "<MaxFeatures>");
    right_marker = StringListUtils::searchPrefix(left_marker, info_file.end(), "</MaxFeatures>");
    while (++left_marker != right_marker)
    {
      mp_.feature_max.push_back(left_marker->toDouble());
    }
    left_marker = StringListUtils::searchPrefix(right_marker, info_file.end(), "<MinFeatures>");
    right_marker = StringListUtils::searchPrefix(right_marker, info_file.end(), "</MinFeatures>");
    while (++left_marker != right_marker)
    {
      mp_.feature_min.push_back(left_marker->toDouble());
    }

    std::vector<IonType> sec_ion_types;

    //Load the secondary types
    left_marker = StringListUtils::searchPrefix(left_marker, info_file.end(), "<SecondaryTypes>");
    right_marker = StringListUtils::searchPrefix(left_marker, info_file.end(), "</SecondaryTypes>");

    if (++left_marker == right_marker)
    {
      mp_.number_intensity_levels = 0;
      return; //if no secondary types are set we quit
    }

    left_marker = StringListUtils::searchPrefix(left_marker, info_file.end(), "<IntensityLevels>");
    right_marker = StringListUtils::searchPrefix(left_marker, info_file.end(), "</IntensityLevels>");
    if (left_marker != right_marker)
    {
      mp_.number_intensity_levels = (++left_marker)->toInt();
    }

    left_marker = StringListUtils::searchPrefix(left_marker, info_file.end(), "<IntensityBinBoarders>");
    right_marker = StringListUtils::searchPrefix(left_marker, info_file.end(), "</IntensityBinBoarders>");
    while (++left_marker != right_marker)
    {
      mp_.intensity_bin_boarders.push_back(left_marker->toDouble());
    }

    left_marker = StringListUtils::searchPrefix(left_marker, info_file.end(), "<IntensityBinValues>");
    right_marker = StringListUtils::searchPrefix(left_marker, info_file.end(), "</IntensityBinValues>");
    while (++left_marker != right_marker)
    {
      mp_.intensity_bin_values.push_back(left_marker->toDouble());
    }

    left_marker = StringListUtils::searchPrefix(left_marker, info_file.end(), "<IonType>");
    right_marker = StringListUtils::searchPrefix(left_marker, info_file.end(), "</IonType>");
    while (left_marker != right_marker)
    {
      //load type information
      ++left_marker;
      Residue::ResidueType res = (Residue::ResidueType)(*left_marker++).toInt();
      EmpiricalFormula loss(*(left_marker++));
      UInt charge = (left_marker++)->toInt();
      IonType actual_type(res, loss, charge);
      sec_ion_types.push_back(actual_type);

      //load conditional probabilities
      left_marker = StringListUtils::searchPrefix(left_marker, info_file.end(), "<ConditionalProbabilities>");

      Size region = 0;
      left_marker = StringListUtils::searchPrefix(left_marker, info_file.end(), "<Region " + String(region) + ">");
      right_marker = StringListUtils::searchPrefix(left_marker, info_file.end(), "</IonType>");
      while (left_marker < right_marker)
      {
        mp_.conditional_prob[std::make_pair(actual_type, region)].assign(mp_.number_intensity_levels, std::vector<double>(mp_.number_intensity_levels));
        for (Size prim = 0; prim < mp_.number_intensity_levels; ++prim)
        {
          for (Size sec = 0; sec < mp_.number_intensity_levels; ++sec)
          {
            mp_.conditional_prob[std::make_pair(actual_type, region)][prim][sec] = (++left_marker)->toDouble();
          }
        }
        ++region;
        left_marker = StringListUtils::searchPrefix(left_marker, info_file.end(), "<Region " + String(region) + ">");
      }
      mp_.number_regions = region;

      left_marker = StringListUtils::searchPrefix(right_marker, info_file.end(), "<IonType>");
      right_marker = StringListUtils::searchPrefix(left_marker, info_file.end(), "</IonType>");
    }


    //map the secondary types to their corresponding primary types
    for (Size i = 0; i < sec_ion_types.size(); ++i)
    {
      const IonType& tmp = sec_ion_types[i];
      if (tmp.residue == Residue::BIon || tmp.residue == Residue::AIon || tmp.residue == Residue::CIon)
      {
        mp_.secondary_types[IonType(Residue::BIon)].push_back(tmp);
      }
      else if (tmp.residue == Residue::YIon || tmp.residue == Residue::XIon || tmp.residue == Residue::ZIon)
      {
        mp_.secondary_types[IonType(Residue::YIon)].push_back(tmp);
      }
    }
  }

  void SvmTheoreticalSpectrumGenerator::simulate(PeakSpectrum& spectrum, const AASequence& peptide, boost::random::mt19937_64& rng, Size precursor_charge)
  {
    if (mp_.class_models.empty() || mp_.reg_models.empty() || mp_.ion_types.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "no svm models loaded. Call load function before using simulate");
    }

    if (peptide.empty())
    {
      return;
    }

    //load parameters
    bool add_isotopes = param_.getValue("add_isotopes").toBool();
    Size max_isotope = (Size)param_.getValue("max_isotope");
    bool add_losses = !(param_.getValue("hide_losses").toBool());
    double relative_loss_intens = param_.getValue("relative_loss_intensity");
    bool add_first_nterminals = param_.getValue("add_first_prefix_ion").toBool();
    bool add_metainfo = param_.getValue("add_metainfo").toBool();

    Int simulation_type = (Int)param_.getValue("svm_mode");

    if (add_metainfo)
    {
      if (spectrum.getIntegerDataArrays().size() == 0)
      {
        spectrum.getIntegerDataArrays().resize(1);
        spectrum.getIntegerDataArrays()[0].setName("Charges");
      }
      if (spectrum.getStringDataArrays().size() == 0)
      {
        spectrum.getStringDataArrays().resize(1);
        spectrum.getStringDataArrays()[0].setName("IonNames");
      }
    }

    std::vector<std::set<String> > possible_n_term_losses(peptide.size());
    std::vector<std::set<String> > possible_c_term_losses(peptide.size());

    UInt ion_nr = 0;

    for (Size i = 1; i < peptide.size(); ++i)
    {
      possible_n_term_losses[i] = possible_n_term_losses[i - 1];

      //check for possible losses on the n-terminal ions
      if (peptide[i - 1].hasNeutralLoss())
      {
        std::vector<EmpiricalFormula> loss_formulas = peptide[i - 1].getLossFormulas();
        for (Size k = 0; k != loss_formulas.size(); ++k)
        {
          possible_n_term_losses[i].insert(loss_formulas[k].toString());
        }
      }

      //check for possible losses on the c-terminal ions
      for (Size pos = i; pos < peptide.size(); ++pos)
      {
        if (peptide[pos].hasNeutralLoss())
        {
          std::vector<EmpiricalFormula> loss_formulas = peptide[pos].getLossFormulas();
          for (Size k = 0; k != loss_formulas.size(); ++k)
          {
            possible_c_term_losses[i].insert(loss_formulas[k].toString());
          }
        }
      }
    }

    std::vector<std::pair<std::pair<IonType, double>, Size> > peaks_to_generate;

    for (Size type_nr = 0; type_nr < mp_.ion_types.size(); ++type_nr)
    {
      if (hide_type_[mp_.ion_types[type_nr]])
      {
        continue;
      }

      Residue::ResidueType residue = mp_.ion_types[type_nr].residue;
      EmpiricalFormula loss_formula = mp_.ion_types[type_nr].loss;

      std::vector<double> predicted_intensity(peptide.size(), 0.);
      std::vector<bool> predicted_class(peptide.size(), false);

#pragma omp parallel for
      for (SignedSize i = 1; i < (SignedSize)peptide.size(); ++i)
      {
        //determine whether type is N- or C-terminal
        if (residue == Residue::AIon || residue == Residue::BIon || residue == Residue::CIon)
        {
          //usually no b1, a1, c1 ion added
          if ((i < 2) && !add_first_nterminals)
          {
            continue;
          }
          //if loss is not supported or no loss ions shall be generated -- continue
          if (!loss_formula.isEmpty() && (!possible_n_term_losses[i].count(loss_formula.toString()) || !add_losses))
          {
            continue;
          }
        }
        //same for c-terminal ions
        else if (residue == Residue::XIon || residue == Residue::YIon || residue == Residue::ZIon)
        {
          if (!loss_formula.isEmpty() && (!possible_c_term_losses[i].count(loss_formula.toString()) || !add_losses))
          {
            continue;
          }
        }
        else
        {
          LOG_ERROR << "Requested unsupported ion type" << std::endl;
        }

        DescriptorSet descriptor;
        generateDescriptorSet_(peptide, i - 1, mp_.ion_types[type_nr], precursor_charge, descriptor);
        if (mp_.scaling_lower != mp_.scaling_upper)
        {
          scaleDescriptorSet_(descriptor, mp_.scaling_lower, mp_.scaling_upper);
        }

        if (simulation_type == 0)
        {
          std::vector<double> tmp_out;
          std::vector<svm_node*> tmp_in(1, &descriptor.descriptors[0]);

          mp_.class_models[type_nr].get()->predict(tmp_in, tmp_out);
          predicted_class[i] = tmp_out[0];
        }

        if (simulation_type == 1)
        {
          std::vector<double> tmp_out;
          std::vector<svm_node*> tmp_in(1, &descriptor.descriptors[0]);

          mp_.reg_models[type_nr].get()->predict(tmp_in, tmp_out);
          predicted_intensity[i] = std::min(std::max(0., tmp_out[0]), 1.0);
        }
      }
      //end of parallel execution



      for (Size i = 1; i < peptide.size(); ++i)
      {
        AASequence prefix = peptide.getPrefix(i);

        Size bin = 0;

        if (simulation_type == 0)
        {
          if (predicted_class[i] == 1 && mp_.static_intensities[residue] != 0)
          {
            double intens = mp_.static_intensities[residue];
            if (!loss_formula.isEmpty())
            {
              intens = intens * relative_loss_intens;
            }
            if (intens > 0)
            {
              peaks_to_generate.push_back(std::make_pair(std::make_pair(mp_.ion_types[type_nr], intens), i));
            }
          }
        }
        else if (simulation_type == 1)
        {
          peaks_to_generate.push_back(std::make_pair(std::make_pair(mp_.ion_types[type_nr], predicted_intensity[i]), i));

          //binning the predicted intensity
          if (predicted_intensity[i] > 0)
          {
            bin = 1;
            while ((1 + bin < mp_.number_intensity_levels) && predicted_intensity[i] > mp_.intensity_bin_boarders[bin - 1])
            {
              ++bin;
            }
          }
        }

        //now generate all secondary types for the primary one
        for (std::vector<IonType>::iterator it = mp_.secondary_types[mp_.ion_types[type_nr]].begin(); it != mp_.secondary_types[mp_.ion_types[type_nr]].end(); ++it)
        {
          if (hide_type_[*it] || (!add_losses && !it->loss.isEmpty()))
          {
            continue;
          }

          if (simulation_type == 0)
          {
            if (predicted_class[i] != 0 && mp_.static_intensities[it->residue] != 0  && hide_type_[*it] == false) //no secondary ion if primary is classified as missing
            {
              double intens = mp_.static_intensities[it->residue];
              if (!it->loss.isEmpty())
              {
                intens = intens * relative_loss_intens;
              }
              if (intens > 0)
              {
                peaks_to_generate.push_back(std::make_pair(std::make_pair(*it, intens), i));
              }
            }
          }
          else if (simulation_type == 1)
          {
            //sample intensities for secondary types
            Size region = std::min(mp_.number_regions - 1, (Size)floor(mp_.number_regions * prefix.getMonoWeight(Residue::Internal) / peptide.getMonoWeight()));
            std::vector<double>& props = mp_.conditional_prob[std::make_pair(*it, region)][bin];
            std::vector<double> weights;
            std::transform(props.begin(), props.end(), std::back_inserter(weights), boost::bind(std::multiplies<double>(), _1, 10));
            boost::random::discrete_distribution<Size> ddist(weights.begin(), weights.end());
            Size binned_int = ddist(rng);

            if (binned_int != 0)
            {
              peaks_to_generate.push_back(std::make_pair(std::make_pair(*it, mp_.intensity_bin_values[binned_int - 1]), i));
            }
          }
        }
      }
    }

    //Now create the actual spectrum containing also isotopic peaks if selected and precursor peaks etc.
    for (Size i = 0; i < peaks_to_generate.size(); ++i)
    {
      IonType type = peaks_to_generate[i].first.first;
      double intensity = peaks_to_generate[i].first.second;
      Residue::ResidueType residue = type.residue;
      EmpiricalFormula loss_formula = type.loss;
      Int charge = type.charge;
      Size pos = peaks_to_generate[i].second;

      AASequence ion;

      if (residue == Residue::AIon || residue == Residue::BIon || residue == Residue::CIon)
      {
        ion = peptide.getPrefix(pos);
      }
      else if (residue == Residue::XIon || residue == Residue::YIon || residue == Residue::ZIon)
      {
        ion = peptide.getSuffix(peptide.size() - pos);
      }

      EmpiricalFormula ion_formula = ion.getFormula(residue, charge) - loss_formula;
      double mz_pos = ion_formula.getMonoWeight() / charge;

      String ion_name = ResidueTypeToString_(residue) + " " + loss_formula.toString() + " " + String(ion_nr) + String(charge, '+');

      if (add_isotopes)
      {
        IsotopeDistribution dist = ion_formula.getIsotopeDistribution((Int)max_isotope);
        Size j = 0;
        for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it, ++j)
        {
          spectrum.push_back(Peak1D(mz_pos + (double)j * Constants::C13C12_MASSDIFF_U / charge, intensity * it->second));
          if (add_metainfo)
          {
            spectrum.getStringDataArrays()[0].push_back(ion_name);
            spectrum.getIntegerDataArrays()[0].push_back(charge);
          }
        }
      }
      else
      {
        spectrum.push_back(Peak1D(mz_pos, intensity));
        if (add_metainfo)
        {
          spectrum.getStringDataArrays()[0].push_back(ion_name);
          spectrum.getIntegerDataArrays()[0].push_back(charge);
        }
      }
    }

    spectrum.sortByPosition();
  }

  void SvmTheoreticalSpectrumGenerator::scaleDescriptorSet_(DescriptorSet& desc, double lower, double upper)
  {
    std::vector<svm_node> tmp_desc;
    std::vector<svm_node>::iterator it;
    Size feature_index = 1;
    Size num_features = mp_.feature_max.size();
    for (it = desc.descriptors.begin(); it != (desc.descriptors.end() - 1); ++it)
    {
      while (feature_index < (Size)it->index)
      {
        double tmp_value = 0;
        scaleSingleFeature_(tmp_value, lower, upper, mp_.feature_min[feature_index - 1], mp_.feature_max[feature_index - 1]);
        if (tmp_value != 0)
        {
          svm_node tmp_node = {(Int)feature_index, tmp_value};
          tmp_desc.push_back(tmp_node);
        }
        ++feature_index;
      }
      scaleSingleFeature_(it->value, lower, upper, mp_.feature_min[feature_index - 1], mp_.feature_max[feature_index - 1]);
      if (it->value != 0)
      {
        tmp_desc.push_back(*it);
      }
      ++feature_index;
    }

    while (feature_index <= num_features)
    {
      double tmp_value = 0;
      scaleSingleFeature_(tmp_value, lower, upper, mp_.feature_min[feature_index - 1], mp_.feature_max[feature_index - 1]);
      if (tmp_value != 0)
      {
        svm_node tmp_node = {(Int)feature_index, tmp_value};
        tmp_desc.push_back(tmp_node);
      }
      ++feature_index;
    }
    svm_node tmp_node = {-1, -1};
    tmp_desc.push_back(tmp_node);

    desc.descriptors = tmp_desc;
  }

  void SvmTheoreticalSpectrumGenerator::updateMembers_()
  {
    hide_type_.clear();

    hide_type_[IonType(Residue::BIon, EmpiricalFormula(""), 1)] = param_.getValue("hide_b_ions").toBool();
    hide_type_[IonType(Residue::YIon, EmpiricalFormula(""), 1)] = param_.getValue("hide_y_ions").toBool();

    hide_type_[IonType(Residue::BIon, EmpiricalFormula(""), 2)] = param_.getValue("hide_b2_ions").toBool();
    hide_type_[IonType(Residue::YIon, EmpiricalFormula(""), 2)] = param_.getValue("hide_y2_ions").toBool();

    hide_type_[IonType(Residue::AIon, EmpiricalFormula(""), 1)] = param_.getValue("hide_a_ions").toBool();
    hide_type_[IonType(Residue::CIon, EmpiricalFormula(""), 1)] = param_.getValue("hide_c_ions").toBool();
    hide_type_[IonType(Residue::XIon, EmpiricalFormula(""), 1)] = param_.getValue("hide_x_ions").toBool();
    hide_type_[IonType(Residue::ZIon, EmpiricalFormula(""), 1)] = param_.getValue("hide_z_ions").toBool();

    mp_.static_intensities[Residue::BIon] = !hide_type_[IonType(Residue::BIon)] ? (double)param_.getValue("b_intensity") : 0;
    mp_.static_intensities[Residue::YIon] = !hide_type_[IonType(Residue::YIon)] ? (double)param_.getValue("y_intensity") : 0;
    mp_.static_intensities[Residue::AIon] = !hide_type_[IonType(Residue::AIon)] ? (double)param_.getValue("a_intensity") : 0;
    mp_.static_intensities[Residue::CIon] = !hide_type_[IonType(Residue::CIon)] ? (double)param_.getValue("c_intensity") : 0;
    mp_.static_intensities[Residue::XIon] = !hide_type_[IonType(Residue::XIon)] ? (double)param_.getValue("x_intensity") : 0;
    mp_.static_intensities[Residue::ZIon] = !hide_type_[IonType(Residue::ZIon)] ? (double)param_.getValue("z_intensity") : 0;
  }

} //namespace
