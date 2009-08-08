// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/AdvancedTheoreticalSpectrumGenerator.h>
#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/FORMAT/TextFile.h>

namespace OpenMS
{

  //Function to compute the structure of the probabilistic network by performing a minimal
  //spanning tree computation
  void AdvancedTheoreticalSpectrumGenerator::TreeAugmentedNetwork::generateTree(std::vector<Int> &has_parent)
  {
    typedef std::vector<TanEdge>::const_iterator EdgeConstIter;
    typedef std::map<UInt, UInt>::iterator LabelIter;

    std::map<UInt, UInt>tmp_node_labels;
    for(EdgeConstIter it = edges_.begin(); it!=edges_.end(); ++it)
    {
      tmp_node_labels[it->left_node]=it->left_node;
      tmp_node_labels[it->right_node] = it->right_node;
    }

    UInt min_id=tmp_node_labels.begin()->first;
    UInt max_id=tmp_node_labels.rbegin()->first;

    has_parent.clear();
    has_parent.assign(max_id+1, -1);


    //sort edges in increasing order
    std::sort(edges_.begin(), edges_.end());

    EdgeConstIter edges_iter = edges_.begin();
    EdgeConstIter edges_end = edges_.end();

    //vector holding the edges selected for the spanning tree
    std::vector<TanEdge> selected_edges;

    //selecting the edges of the minimum spanning tree
    while (edges_iter != edges_end)
    {
      //if both ends of this edge do belong to different components
      //add the edge to the tree and adapt labels.
      if (tmp_node_labels[edges_iter->left_node] != tmp_node_labels[edges_iter->right_node])
      {
        UInt min_label = std::min(tmp_node_labels[edges_iter->left_node], tmp_node_labels[edges_iter->right_node]);
        UInt max_label = std::max(tmp_node_labels[edges_iter->left_node], tmp_node_labels[edges_iter->right_node]);

        //relabel all nodes in the connected component with id max_label to min_label
        LabelIter label_it = tmp_node_labels.begin();
        LabelIter labels_end = tmp_node_labels.end();

        while (label_it != labels_end)
        {
          if (label_it->second == max_label)
          {
            label_it->second = min_label;
          }
          ++label_it;
        }

        //add the edge to the set of selected edges
        selected_edges.push_back(*edges_iter);
      }
      ++edges_iter;
    }

    LabelIter label_it = tmp_node_labels.begin();
    UInt first_elem=label_it->second;
    LabelIter labels_end = tmp_node_labels.end();
    while(label_it !=labels_end && label_it->second == first_elem)
    {
      ++label_it;
    }

    if(label_it!=labels_end)
    {
      throw(Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Input graph is no connected component"));
    }

    //with the selected edges now build the tree
    //vector of adjacency lists for easy reconstruction
    std::vector<std::list<UInt> > adjacency_lists(max_id+1);
    typedef std::list<UInt>::const_iterator AdjaListConstIter;

    EdgeConstIter selected_iter = selected_edges.begin();
    EdgeConstIter selected_end = selected_edges.end();

    while (selected_iter != selected_end)
    {
      adjacency_lists[selected_iter->left_node].push_back(selected_iter->right_node);
      adjacency_lists[selected_iter->right_node].push_back(selected_iter->left_node);
      ++selected_iter;
    }
    selected_edges.clear();

    //creating the rooted tree in a DFS manner using a stack
    std::vector<UInt> node_stack;
    node_stack.push_back(min_id); //use the node with lowest id as root
    nodes_in_dfs_order_.clear();
    nodes_in_dfs_order_.push_back(min_id);

    std::vector<bool> discovered(max_id, false);

    while (!node_stack.empty())
    {
      UInt actual_node = node_stack.back();
      node_stack.pop_back();

      AdjaListConstIter adja_iter = adjacency_lists[actual_node].begin();
      AdjaListConstIter adja_end = adjacency_lists[actual_node].end();

      while (adja_iter != adja_end)
      {
        if (!discovered[*adja_iter])
        {
          node_stack.push_back(*adja_iter);
          nodes_in_dfs_order_.push_back(*adja_iter);
        }
        else
        {
          has_parent[actual_node] = *adja_iter;
        }
        ++adja_iter;
      }
      discovered[actual_node] = true;
    }
  }

  AdvancedTheoreticalSpectrumGenerator::AdvancedTheoreticalSpectrumGenerator() :
    TheoreticalSpectrumGenerator(),
    conditional_probabilities_(0),
    tan_(0),
    number_of_intensity_levels_(0),
    number_of_sectors_(0),
    ion_types_(0)
  {
    this->setName("AdvancedTheoreticalSpectrumGenerator");
    defaults_.setValue("add_isotopes", 0, "If set to 1 isotope peaks of the product ion peaks are added");
    defaults_.setValue("max_isotope", 2, "Defines the maximal isotopic peak which is added, add_isotopes must be set to 1");
    defaults_.setValue("add_metainfo", 0, "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");
    defaults_.setValue("add_losses", 0, "Adds common losses to those ion expect to have them, only water and ammonia loss is considered");
    defaults_.setValue("add_precursor_peaks", 0, "Adds peaks of the precursor to the spectrum, which happen to occur sometimes");
    defaults_.setValue("model_file_name", "examples/simulation/MSMSSim.model", "Name of the probabilistic Model file");

    //defaults_.setValue("relative_loss_intensity", 0.1, "Intensity of loss ions, in relation to the intact ion intensity");

    // precursor intensity
    defaults_.setValue("precursor_intensity", 1.0, "Intensity of the precursor peak");
    defaults_.setValue("precursor_H2O_intensity", 1.0, "Intensity of the H2O loss peak of the precursor");
    defaults_.setValue("precursor_NH3_intensity", 1.0, "Intensity of the NH3 loss peak of the precursor");

    defaultsToParam_();

    // just in case someone wants the ion names;
    p_.metaRegistry().registerName("IonName", "Name of the ion");
  }

  AdvancedTheoreticalSpectrumGenerator::AdvancedTheoreticalSpectrumGenerator(const AdvancedTheoreticalSpectrumGenerator& rhs) :
    TheoreticalSpectrumGenerator(rhs),
    conditional_probabilities_(rhs.conditional_probabilities_),
    tan_(rhs.tan_),
    number_of_intensity_levels_(rhs.number_of_intensity_levels_),
    number_of_sectors_(rhs.number_of_sectors_),
    ion_types_(rhs.ion_types_)
  {
    updateMembers_();
  }

  AdvancedTheoreticalSpectrumGenerator& AdvancedTheoreticalSpectrumGenerator::operator =(const AdvancedTheoreticalSpectrumGenerator& rhs)
  {
    if (this != &rhs)
    {
      TheoreticalSpectrumGenerator::operator=(rhs);
      conditional_probabilities_ = rhs.conditional_probabilities_;
      tan_ = rhs.tan_;
      number_of_intensity_levels_ = rhs.number_of_intensity_levels_;
      number_of_sectors_ = rhs.number_of_sectors_;
      ion_types_ = rhs.ion_types_;

      updateMembers_();
    }
    return *this;
  }

  AdvancedTheoreticalSpectrumGenerator::~AdvancedTheoreticalSpectrumGenerator()
  {
  }

  void AdvancedTheoreticalSpectrumGenerator::simulate(RichPeakSpectrum &spectrum, const AASequence &peptide, const gsl_rng *rng, Int charge)
  {
    //load parameters
    bool add_isotopes = (UInt)param_.getValue("add_isotopes");
    UInt max_isotope = (UInt)param_.getValue("max_isotope");
    bool add_metainfo = (UInt)param_.getValue("add_metainfo");
    bool add_precursor_peaks = (UInt)param_.getValue("add_precursor_peaks");

    //get an instance of the IndexConverter functor
    IndexConverter index_converter;

    //store the structure of the tree (for each sector)
    std::vector<std::vector<Int> > is_child_of(number_of_sectors_);
    //store a safe ordering of the nodes (for each sector)
    std::vector<std::vector<UInt> >ordered_nodes(number_of_sectors_);

    gsl_ran_discrete_t * gsl_gen=0;

    //generate the structures for each sector
    for(Size sector=0; sector<number_of_sectors_; ++sector)
    {
      tan_[sector].generateTree(is_child_of[sector]);
      tan_[sector].getDFSOrder(ordered_nodes[sector]);
    }

    std::set<String>possible_n_term_losses;
    std::set<String>possible_c_term_losses;

    AASequence prefix, suffix;

    const AASequence * ion = 0;

    DoubleReal parent_mass = peptide.getMonoWeight(Residue::Full);
    for (Size i = 1; i < peptide.size(); ++i)
    {
      prefix=peptide.getPrefix(i);
      suffix=peptide.getSuffix(peptide.size()-i);

      DoubleReal prefix_mass = prefix.getMonoWeight(Residue::NTerminal);

      if (peptide[i-1].hasNeutralLoss())
      {
        std::vector<EmpiricalFormula> loss_formulas = peptide[i-1].getLossFormulas();
        for (Size k = 0; k != loss_formulas.size(); ++k)
        {
          possible_n_term_losses.insert(loss_formulas[k].getString());
        }
      }

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

      //compute the sector of the prefix mass
      UInt sector = 1;
      while (prefix_mass / parent_mass > (DoubleReal) sector / number_of_sectors_)
      {
        ++sector;
      }
      --sector;

      std::vector<UInt>generated_intensity(ordered_nodes[sector].size(),0);

      for(Size node=0; node<ordered_nodes[sector].size(); ++node)
      {
        UInt type_id = ordered_nodes[sector][node];
        IonType ion_type = ion_types_[type_id];
        Residue::ResidueType residue = ion_type.residue;
        EmpiricalFormula loss_formula = ion_type.loss;
        Int charge = ion_type.charge;

        //determine whether type is N- or C-terminal
        if (residue == Residue::AIon || residue == Residue::BIon || residue == Residue::CIon)
        {
          //no b,a,c ion for position 1
          if(i<2 || (!loss_formula.isEmpty() && possible_n_term_losses.find(loss_formula.getString())==possible_n_term_losses.end()))
          {
            continue;
          }
          ion = &prefix;
        }
        else if (residue == Residue::XIon || residue == Residue::YIon || residue == Residue::ZIon)
        {
          if(!loss_formula.isEmpty() && possible_c_term_losses.find(loss_formula.getString())==possible_c_term_losses.end())
          {
            continue;
          }
          ion = &suffix;
        }
        else
        {
          std::cerr << "Cannot create peaks of that ion type" << std::endl;
        }

        EmpiricalFormula ion_formula = ion->getFormula(residue, charge) - loss_formula;
        DoubleReal mz_pos = ion_formula.getMonoWeight() / charge;

        UInt conditional_intensity=0;

        if(is_child_of[sector][type_id]!=-1)
          conditional_intensity = generated_intensity[is_child_of[sector][type_id]];

        DoubleReal * tmp_pointer = &(conditional_probabilities_[sector][index_converter(type_id, 0, conditional_intensity, number_of_intensity_levels_)]);
        //std::cerr<<tmp_pointer[0]<<"  "<<tmp_pointer[1]<<" "<<tmp_pointer[2]<<"  "<<tmp_pointer[3]<<"  "<<tmp_pointer[4]<<std::endl; //DEBUG

        gsl_gen = gsl_ran_discrete_preproc(number_of_intensity_levels_, tmp_pointer);
        size_t intensity = gsl_ran_discrete(rng, gsl_gen);

        generated_intensity[type_id]=intensity;

        String ion_name = String(residue) + loss_formula.getString() + String(i) + String(charge, '+');

        if (add_isotopes)
        {
          IsotopeDistribution dist = ion_formula.getIsotopeDistribution(max_isotope);
          UInt j(0);
          for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it, ++j)
          {
            p_.setMZ(mz_pos + j / charge);
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
          //std::cout<<"added ion type: "<<ion_name<<"  with intensity: "<<intensity<<"  at position: "<<mz_pos<<std::endl;
        }
      }
    }

    //add precursor peaks if selected
    if(add_precursor_peaks)
    {
      addPrecursorPeaks(spectrum, peptide, charge);
    }
    spectrum.sortByPosition();
  }


  void AdvancedTheoreticalSpectrumGenerator::loadProbabilisticModel()
  {
    //std::cout << "check point readwrite" << std::endl;

    //store the model in text_file
    String model_file_name = (String)param_.getValue("model_file_name");

    if (! File::readable( model_file_name ) )
    { // look in OPENMS_DATA_PATH
      model_file_name = File::find( model_file_name );
    }

    TextFile text_file;
    text_file.load(model_file_name);

    TextFile::iterator left_marker = text_file.search("<IonTypes>");
    TextFile::iterator right_marker = text_file.search("<IntensityLevels>");

    ++left_marker;
    //read the set of IonTypes
    while(left_marker<right_marker-2)
    {
      Residue::ResidueType res = (Residue::ResidueType)(*left_marker++).toInt();
      EmpiricalFormula loss(*(left_marker++));
      UInt charge(left_marker++->toInt());
      ion_types_.push_back(IonType(res,loss,charge));
    }

    //read the number of intensity levels
    number_of_intensity_levels_=  (++right_marker)->toInt();

    //read the number of sectors
    left_marker=text_file.search("<Sectors>");
    number_of_sectors_=(++left_marker)->toInt();

    //read the tree structures for each sector as is_parent_of relations
    left_marker=text_file.search("<TreeStructures>");
    std::vector<TreeAugmentedNetwork::TanEdge>edges;
    for(Size sector = 0; sector < number_of_sectors_; ++sector)
    {
      left_marker=text_file.search("<Sector "+ String(sector)+">");
      ++left_marker;
      for(Size i=0; i<ion_types_.size();++i)
      {
        if(left_marker->toInt()>-1)
        {
          TreeAugmentedNetwork::TanEdge edge = {i, left_marker->toInt(), -1};
          edges.push_back(edge);
        }
        ++left_marker;
      }
      tan_.push_back(TreeAugmentedNetwork(edges));
      edges.clear();
    }

    //read the conditional probabilities for each sector
    IndexConverter index_converter;
    Size number_of_prob_entries=index_converter(ion_types_.size()-1, number_of_intensity_levels_-1, number_of_intensity_levels_-1, number_of_intensity_levels_)+1;
    conditional_probabilities_.reserve(number_of_sectors_);

    //temporary storage
    std::vector<DoubleReal>tmp_condit;
    tmp_condit.reserve(number_of_prob_entries);

    left_marker=text_file.search("<ConditionalProbabilities>");
    for(Size sector = 0; sector < number_of_sectors_; ++sector)
    {
      //std::cout<<"Out Sector: "<<sector<<std::endl;
      left_marker=text_file.search(left_marker, "<Sector");
      for(Size i=0; i<number_of_prob_entries;++i)
      {
        tmp_condit.push_back((++left_marker)->toDouble());
      }
      conditional_probabilities_.push_back(tmp_condit);
      tmp_condit.clear();
    }
  }



  void AdvancedTheoreticalSpectrumGenerator::writeProbabilisticModel(const String &file_name)
  {
    //store the model in text_file
    TextFile text_file;

    //store the set of selected IonTypes
    text_file.push_back("<IonTypes>\n");
    for(std::vector<IonType>::iterator m_it=ion_types_.begin(); m_it!=ion_types_.end(); ++m_it)
    {
      text_file.push_back(m_it->residue);
      text_file.push_back(m_it->loss.getString());
      text_file.push_back(m_it->charge);
    }

    //store the number of intensity levels
    text_file.push_back("<IntensityLevels>");
    text_file.push_back(number_of_intensity_levels_);

    //store the number of sectors
    text_file.push_back("<Sectors>");
    text_file.push_back(number_of_sectors_);

    std::vector<Int>has_parent;

    //store the tree structures for each sector as is_parent_of relations
    text_file.push_back("<TreeStructures>");
    for(Size sector = 0; sector < number_of_sectors_; ++sector)
    {
      tan_[sector].generateTree(has_parent);
      text_file.push_back("<Sector "+ String(sector)+">");
      for(std::vector<Int>::iterator it=has_parent.begin(); it!=has_parent.end(); ++it)
      {
        text_file.push_back(*it);
      }
    }

    //store the conditional probabilities for each sector
    text_file.push_back("<ConditionalProbabilities>");
    for(Size sector = 0; sector < number_of_sectors_; ++sector)
    {
      text_file.push_back("<Sector " + String(sector) + ">");
      for(Size i=0; i<conditional_probabilities_[sector].size();++i)
      {
        text_file.push_back(conditional_probabilities_[sector][i]);
      }
    }

    text_file.store(file_name);
  }

}//Namespace
