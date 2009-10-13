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
// $Authors: Sandro Andreotti $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/AdvancedTheoreticalSpectrumGenerator.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>

//#define SGNT_DEBUG
#undef SGNT_DEBUG

using namespace std;
using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page UTILS_SpectrumGeneratorNetworkTrainer SpectrumGeneratorNetworkTrainer
	
	@brief Trainer for probabilistic network as input for AdvancedSpectrumGenerator.
	
  This application requires a list of annotated spectra and generates a bayesian network
  with tree structure. For each pair of ion types (i.e. a,b,c,x,y,z + losses) the mutual
  information is computed. Finally the application computes a spanning tree that maximizes
  the total mutual information content. In the resulting bayesian network the probability
  for each ion type to occur with a certain intensity depends only on his parent ion type
  in the tree.
	
	@note This tool is experimental!
		
	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_SpectrumGeneratorNetworkTrainer.cli 
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


  class SpectrumGeneratorNetworkTrainer
  	: public TOPPBase
  {
    typedef std::vector<Int> IntVec;
    typedef std::vector<IntVec> IntMatrix;
    typedef std::vector<UInt> UIntVec;
    typedef std::vector<UIntVec> UIntMatrix;
    typedef std::vector<DoubleReal> DRealVec;
    typedef std::vector<DRealVec> DRealMatrix;
    typedef AdvancedTheoreticalSpectrumGenerator::TreeAugmentedNetwork TreeAugmentedNetwork;
    typedef AdvancedTheoreticalSpectrumGenerator::IonType IonType;
    typedef AdvancedTheoreticalSpectrumGenerator::IndexConverter IndexConverter;

    public:
      SpectrumGeneratorNetworkTrainer() :
        TOPPBase("SpectrumGeneratorNetworkTrainer", "Trainer for Probabilistic network as input for AdvancedSpectrumGenerator", false)
      {
      }

    protected:
      void registerOptionsAndFlags_()
      {
        // I/O settings
        registerInputFile_("in_spectra", "<file>", "", "Input Training Spectra in mzData", true);
        registerInputFile_("in_identifications", "<file>", "", "Input file with corresponding sequences in IdXML", "", true);
        registerOutputFile_("out_network_model", "<file>", "", "Output model of probabilistic network as textfile", true);

        //considered ion types
        registerIntOption_("add_y_ions", "<Int>", 1, "If set to 1 y-ion peaks will be considered", false, true);
        registerIntOption_("add_b_ions", "<Int>", 1, "If set to 1 b-ion peaks will be considered", false, true);
        registerIntOption_("add_y2_ions", "<Int>", 1, "If set to 1 doubly charged y-ion peaks will be considered", false, true);
        registerIntOption_("add_b2_ions", "<Int>", 1, "If set to 1 doubly charged b-ion peaks will be considered", false, true);
        registerIntOption_("add_a_ions", "<Int>", 1, "If set to 1 a-ion peaks will be considered", false, true);
        registerIntOption_("add_c_ions", "<Int>", 1, "If set to 1 c-ion peaks will be considered", false, true);
        registerIntOption_("add_x_ions", "<Int>", 1, "If set to 1 x-ion peaks will be considered", false, true);
        registerIntOption_("add_z_ions", "<Int>", 1, "If set to 1 z-ion peaks will be considered", false, true);

        //losses
        registerIntOption_("add_losses", "<Int>", 1, "Considers common losses to those ion expect to have them, only water and ammonia loss is considered", false, true);

        //model parameters
        //registerDoubleOption_("frequency_cutoff", "<Double>", 0.05, "Only ion types with a frequency > frequency_cutoff will be added to the model", false);
        registerIntOption_("number_of_sectors", "<Int>", 3, "Each spectrum is split into sectors and probabilities are separately learned for each sector", false, true);
        //registerIntOption_("number_of_intensity_levels", "<Int>", 5, "Peak intensities are discretized into the given number of bins (at least 2)", true, true);
        registerDoubleList_("intensity_level_bins", "<Int>", DoubleList::create("0.05,2.0,10.0"),
            "Borders of normalized intensities for intensity discretization. n border values result in n+1 intensity levels", false, true);
        registerDoubleOption_("delta", "<Double>", 0.5, "Error intervall for each peak", false, true);
      }

      void trainModel()
      {
        //read the options
        UInt number_of_sectors = getIntOption_("number_of_sectors");
        DoubleReal delta = getDoubleOption_("delta");
        DRealVec intensity_limits = getDoubleList_("intensity_level_bins");
        UInt number_of_intensity_levels = (UInt) intensity_limits.size() + 1;

        //file options
        String mzdata_file = getStringOption_("in_spectra");
        String idxml_file = getStringOption_("in_identifications");

        std::vector<IonType> ion_types;

        //build set of ion typed
        if (getIntOption_("add_y_ions"))
          ion_types.push_back(IonType(Residue::YIon, EmpiricalFormula(), 1));
        if (getIntOption_("add_b_ions"))
          ion_types.push_back(IonType(Residue::BIon, EmpiricalFormula(), 1));
        if (getIntOption_("add_x_ions"))
          ion_types.push_back(IonType(Residue::XIon, EmpiricalFormula(), 1));
        if (getIntOption_("add_a_ions"))
          ion_types.push_back(IonType(Residue::AIon, EmpiricalFormula(), 1));
        if (getIntOption_("add_z_ions"))
          ion_types.push_back(IonType(Residue::ZIon, EmpiricalFormula(), 1));
        if (getIntOption_("add_c_ions"))
          ion_types.push_back(IonType(Residue::CIon, EmpiricalFormula(), 1));

        if (getIntOption_("add_losses"))
        {
          EmpiricalFormula loss_ammonia("NH3");
          EmpiricalFormula loss_water("H2O");

          ion_types.push_back(IonType(Residue::BIon, loss_ammonia, 1));
          ion_types.push_back(IonType(Residue::BIon, loss_water, 1));

          ion_types.push_back(IonType(Residue::YIon, loss_ammonia, 1));
          ion_types.push_back(IonType(Residue::YIon, loss_water, 1));
        }

        if (getIntOption_("add_y2_ions"))
          ion_types.push_back(IonType(Residue::YIon, EmpiricalFormula(), 2));
        if (getIntOption_("add_b2_ions"))
          ion_types.push_back(IonType(Residue::BIon, EmpiricalFormula(), 2));

        //Amino acid sequences used for obtaining prefix and suffix mass
        AASequence prefix, suffix;

        //loading data
        RichPeakMap spectra_map;
        std::vector<PeptideIdentification> pep_id_vec;
        std::vector<ProteinIdentification> prot_id_vec;

        String tmp_str;

        MzDataFile().load(mzdata_file, spectra_map);
        IdXMLFile().load(idxml_file, prot_id_vec, pep_id_vec, tmp_str);
        IDMapper().annotate(spectra_map, pep_id_vec, prot_id_vec);

        UInt number_of_ion_types = (UInt) ion_types.size();

        //stores the tans for the each sector
        IntMatrix has_parent_all_sectors;
        has_parent_all_sectors.reserve(number_of_sectors);

        //stores the conditional probs for each sector
        DRealMatrix condit_prob_all_sectors;
        condit_prob_all_sectors.reserve(number_of_sectors);

        for (UInt sector = 0; sector < number_of_sectors; ++sector)
        {
          //counter for the number of spectra that were actually used for training
          UInt train_set_size = 0;

          //only prefix masses with a ratio to parent mass between the lower and upper bounds are considered
          DoubleReal sector_lower_bound = (DoubleReal) sector / number_of_sectors;
          DoubleReal sector_upper_bound = (DoubleReal) (sector + 1) / number_of_sectors;

          UInt matrix_dim_lim = calc_index(number_of_ion_types - 1, number_of_intensity_levels - 1, number_of_intensity_levels) + 1;

          //generate the two_dimensional vector used for storing the conditional mutual information
          DRealMatrix mi(number_of_ion_types, DRealVec(number_of_ion_types, 0.0));

          //matrices for counting the pairwise abundances of ion-with-intensity-level Pairs
          DRealMatrix pairwise_count_true(matrix_dim_lim, DRealVec(matrix_dim_lim, 1.0)); //set each entry to one --> pseudo-count

          //Vectors to store for each ion type the observed intensity level (0=missing)
          UIntVec peak_list_true(matrix_dim_lim, 0);

          //run over all input training spectra
          for (RichPeakMap::iterator map_it = spectra_map.begin(); map_it != spectra_map.end(); ++map_it)
          {
            //TODO optimize this by internalizing the sector loop and get rid of repeated normalizations
            MSSpectrum<RichPeak1D> input_spec_norm(*map_it);

            normalizeIntensity(input_spec_norm, number_of_intensity_levels, intensity_limits);

            //test whether annotation for spectrum is available
            if (!map_it->getPeptideIdentifications()[0].getHits()[0].getSequence().isValid())
            {
              std::cerr << "no annotation available spectrum number " << map_it - spectra_map.begin() << std::endl;
            }
            else
            {
              //this spectrum can be used for training
              ++train_set_size;

              //get the annotation for the given input spectrum file
              AASequence annot(map_it->getPeptideIdentifications()[0].getHits()[0].getSequence());

              //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              //----------------------generate the true positive counts for true prefix positions-------------------------
              //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

              DoubleReal parent_mass = (annot.getMonoWeight(Residue::Full) + 1);

              for (Size frag_pos = 1; frag_pos < annot.size(); ++frag_pos)
              {
                //reset the observations
                std::fill(peak_list_true.begin(), peak_list_true.end(), 0);

                prefix = annot.getPrefix(frag_pos);
                suffix = annot.getSuffix(annot.size() - frag_pos);
                DoubleReal true_prefix_mass = (annot.getPrefix(frag_pos).getMonoWeight(Residue::Internal));

                if (!(true_prefix_mass / parent_mass > sector_upper_bound || true_prefix_mass / parent_mass <= sector_lower_bound))
                {
                  //now check for each ion type whether a peak is abundant at corresponding mass offset
                  for (Size type_nr = 0; type_nr < ion_types.size(); ++type_nr)
                  {
                    DoubleReal true_offset_mass = 0.0;
                    Residue::ResidueType residue = ion_types[type_nr].residue;
                    Int charge = ion_types[type_nr].charge;
                    EmpiricalFormula loss = ion_types[type_nr].loss;

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
                    Size true_nearest_peak_ind = input_spec_norm.findNearest(true_offset_mass);

                    //check whether this peak is within the allowed mass range
                    if (fabs(true_offset_mass - input_spec_norm[true_nearest_peak_ind].getMZ()) <= delta)
                    {
                      UInt intensity_index = (UInt) input_spec_norm[true_nearest_peak_ind].getIntensity();
                      peak_list_true[type_nr] = intensity_index;
                    }
                  }

                  for (Size left_type_nr = 0; left_type_nr < ion_types.size(); ++left_type_nr)
                  {
                    for (Size right_type_nr = 0; right_type_nr < ion_types.size(); ++right_type_nr)
                    {
                      ++pairwise_count_true[calc_index((UInt) left_type_nr, peak_list_true[left_type_nr], number_of_intensity_levels)][calc_index((UInt) right_type_nr,
                          peak_list_true[right_type_nr], number_of_intensity_levels)];
                    }
                  }
                }
              }
            }
          }//end of running over all spectra

#ifdef SGNT_DEBUG
          //DEBUG OUTPUT::
          for (Size left_type_nr = 0; left_type_nr < ion_types.size(); ++left_type_nr)
          {
            for (Size right_type_nr = 0; right_type_nr < ion_types.size(); ++right_type_nr)
            {
              for (UInt level_left = 0; level_left < number_of_intensity_levels; ++level_left)
              {
                for (UInt level_right = 0; level_right < number_of_intensity_levels; ++level_right)
                {
                  std::cout << left_type_nr << "  " << level_left << " " << right_type_nr << "   " << level_right << " " << pairwise_count_true[calc_index(left_type_nr,
                      level_left, number_of_intensity_levels)][calc_index(right_type_nr, level_right, number_of_intensity_levels)] << std::endl;
                }
              }
            }
          }
#endif

          //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          //----------------------generate the statistics of the observed data to build a TAN Network--------------------
          //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          DRealVec background_probs_true(matrix_dim_lim, 0.0);

          //matrices for storing the pairwise probabilities of ion-with-intensity-level Pairs
          DRealMatrix pairwise_prob_true(matrix_dim_lim, DRealVec(matrix_dim_lim, 0.0));

          //compute the pairwise probabilities
          for (Size left_type_nr = 0; left_type_nr < ion_types.size(); ++left_type_nr)
          {
            for (Size right_type_nr = 0; right_type_nr < ion_types.size(); ++right_type_nr)
            {
              UInt true_sum = 0;

              //first compute the total counts
              for (UInt level_left = 0; level_left < number_of_intensity_levels; ++level_left)
              {
                for (UInt level_right = 0; level_right < number_of_intensity_levels; ++level_right)
                {
                  true_sum += (UInt) pairwise_count_true[calc_index((UInt) left_type_nr, level_left, number_of_intensity_levels)][calc_index((UInt) right_type_nr, level_right,
                      number_of_intensity_levels)];
                }
              }

              //now use the total counts to generate the relative amounts
              for (UInt level_left = 0; level_left < number_of_intensity_levels; ++level_left)
              {
                for (UInt level_right = 0; level_right < number_of_intensity_levels; ++level_right)
                {
                  UInt left_index = calc_index((UInt) left_type_nr, level_left, number_of_intensity_levels);
                  UInt right_index = calc_index((UInt) right_type_nr, level_right, number_of_intensity_levels);
                  if (true_sum == 0)
                    pairwise_prob_true[left_index][right_index] = 0;
                  else
                    pairwise_prob_true[left_index][right_index] = (DoubleReal) pairwise_count_true[left_index][right_index] / true_sum;
                }
              }
            }
          }

#ifdef SGNT_DEBUG
          //DEBUG OUTPUT::
          for (Size left_type_nr = 0; left_type_nr < number_of_ion_types; ++left_type_nr)
          {
            for (Size right_type_nr = 0; right_type_nr < number_of_ion_types; ++right_type_nr)
            {
              for (UInt level_left = 0; level_left < number_of_intensity_levels; ++level_left)
              {
                for (UInt level_right = 0; level_right < number_of_intensity_levels; ++level_right)
                {
                  UInt left_index = calc_index(left_type_nr, level_left, number_of_intensity_levels);
                  UInt right_index = calc_index(right_type_nr, level_right, number_of_intensity_levels);
                  std::cout << left_type_nr << "  " << level_left << " " << right_type_nr << "   " << level_right << "  " << pairwise_prob_true[left_index][right_index]
                  << std::endl;
                }
              }
            }
          }
#endif

          //compute background probabilities
          UInt reference_type = 1;
          for (Size left_type_nr = 0; left_type_nr < ion_types.size(); ++left_type_nr)
          {
            for (UInt level_left = 0; level_left < number_of_intensity_levels; ++level_left)
            {
              DoubleReal back_prob_true = 0;

              for (UInt level_ref = 0; level_ref < number_of_intensity_levels; ++level_ref)
              {
                back_prob_true += pairwise_prob_true[calc_index((UInt) left_type_nr, level_left, number_of_intensity_levels)][calc_index((UInt) reference_type, level_ref,
                    number_of_intensity_levels)];
              }
              background_probs_true[calc_index((UInt) left_type_nr, level_left, number_of_intensity_levels)] = back_prob_true;
            }
            reference_type = 0;
          }

#ifdef SGNT_DEBUG
          //DEBUG OUTPUT::
          for (Size type_nr = 0; type_nr < ion_types.size(); ++type_nr)
          {
            for (UInt level = 0; level < number_of_intensity_levels; ++level)
            {
              std::cout << "background:  " << type_nr << "  " << level << "  true: " << background_probs_true[calc_index(type_nr, level, number_of_intensity_levels)] << std::endl;
            }
          }
#endif
          //now compute mutual information between types to as weight for the tan training
          //DoubleReal mi_old;

          for (Size left_type_nr = 0; left_type_nr < ion_types.size(); ++left_type_nr)
          {
            for (Size right_type_nr = left_type_nr + 1; right_type_nr < ion_types.size(); ++right_type_nr)
            {
              for (UInt level_left = 0; level_left < number_of_intensity_levels; ++level_left)
              {
                for (UInt level_right = 0; level_right < number_of_intensity_levels; ++level_right)
                {
                  UInt left_index = calc_index((UInt) left_type_nr, level_left, number_of_intensity_levels);
                  UInt right_index = calc_index((UInt) right_type_nr, level_right, number_of_intensity_levels);

                  DoubleReal background_true_factor = background_probs_true[left_index] * background_probs_true[right_index];
                  DoubleReal pairwise_prob = pairwise_prob_true[left_index][right_index];

                  if (background_true_factor != 0 && pairwise_prob != 0)
                  {
                    //mi_old = mi[left_type_nr][right_type_nr];
                    mi[left_type_nr][right_type_nr] += pairwise_prob * log(pairwise_prob / background_true_factor);
                    //std::cout << mi[left_type_nr][right_type_nr] - mi_old << std::endl;
                  }
                }
              }
            }
          }

#ifdef SGNT_DEBUG
          for (Size left_type_nr = 0; left_type_nr < ion_types.size(); ++left_type_nr)
          {
            for (Size right_type_nr = left_type_nr + 1; right_type_nr < ion_types.size(); ++right_type_nr)
            {
              std::cout << "mi: " << left_type_nr << "  " << right_type_nr << " " << mi[left_type_nr][right_type_nr] << std::endl;
            }
          }
#endif

          //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          //----------------------with the pairwise mutual information build the TAN Network--------------------------
          //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          //generate a set of edges for computation of the TAN
          std::vector<TreeAugmentedNetwork::TanEdge> tan_input_edges;

          for (Size left_type_nr = 0; left_type_nr < ion_types.size(); ++left_type_nr)
          {
            for (Size right_type_nr = left_type_nr + 1; right_type_nr < ion_types.size(); ++right_type_nr)
            {
              TreeAugmentedNetwork::TanEdge edge =
              { (UInt) left_type_nr, (UInt) right_type_nr, -mi[left_type_nr][right_type_nr] };
              tan_input_edges.push_back(edge);
            }
          }

          TreeAugmentedNetwork t_net(tan_input_edges);
          IntVec has_parent;
          t_net.generateTree(has_parent);

#ifdef SGNT_DEBUG
          for (Size type_nr = 0; type_nr < ion_types.size(); ++type_nr)
          {
            std::cout << "has_parent " << type_nr << "  :" << has_parent[type_nr] << std::endl;
          }
#endif

          //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          //----------------------calculate and store the conditional probability network-----------------------------
          //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          IndexConverter indexConverter;

          //vector to store the conditional probabilities
          std::vector<DoubleReal> condit_prob(indexConverter(number_of_ion_types - 1, number_of_intensity_levels - 1, number_of_intensity_levels - 1, number_of_intensity_levels)
              + 1, 0.0);

          //parameter for smoothing the probabilites as described in Bern Paper
          const UInt smoother = 50; //This is *still* Magic

          for (Size child_type_nr = 0; child_type_nr < number_of_ion_types; ++child_type_nr)
          {
            Size parent_type_nr = 0;

            //for the root node which depends on no other type (by default we always choose Type 0 as root)
            if (has_parent[child_type_nr] == -1)
            {
              if (child_type_nr == 0)
              {
                parent_type_nr = 1;
              }
            }
            else
            {
              parent_type_nr = has_parent[child_type_nr];
            }

            for (UInt level_child = 0; level_child < number_of_intensity_levels; ++level_child)
            {
              for (UInt level_parent = 0; level_parent < number_of_intensity_levels; ++level_parent)
              {
                UInt index_condit = indexConverter((UInt) child_type_nr, level_child, level_parent, number_of_intensity_levels);
                UInt index_child = calc_index((UInt) child_type_nr, level_child, number_of_intensity_levels);
                UInt index_parent = calc_index((UInt) parent_type_nr, level_parent, number_of_intensity_levels);

                if (has_parent[child_type_nr] == -1)
                {
                  condit_prob[index_condit] = background_probs_true[index_child];
                }
                else
                {
                  condit_prob[index_condit] = pairwise_prob_true[index_child][index_parent] / background_probs_true[index_parent];
                }

                DoubleReal background = background_probs_true[index_child];
                condit_prob[index_condit] = ((DoubleReal) train_set_size / (train_set_size + smoother)) * condit_prob[index_condit] + (((DoubleReal) smoother / (train_set_size
                    + smoother)) * background);
              }
            }
          }

          //store the tan for the actual sector
          has_parent_all_sectors.push_back(has_parent);
          //store the conditional probs for the actual sector
          condit_prob_all_sectors.push_back(condit_prob);

#ifdef SGNT_DEBUG
          //DEBUG OUTPUT OF CONDITIONAL PROBS + COMPARISON TO BACKGROUND PROBS
          std::vector<Int> tmp_has_parent = has_parent;
          std::cout << "SECTOR: " << sector << std::endl;
          for (Size i = 0; i < tmp_has_parent.size(); ++i)
          {
            for (UInt lev_i = 0; lev_i < number_of_intensity_levels; ++lev_i)
            {
              std::cout << "Background for level  " << lev_i << "   :" << background_probs_true[calc_index(i, lev_i, number_of_intensity_levels)] << std::endl;
              for (UInt lev_parent = 0; lev_parent < number_of_intensity_levels; ++lev_parent)
              {
                UInt index_true = indexConverter(i, lev_i, lev_parent, number_of_intensity_levels);
                std::cout << "condit_prob for parent level " << lev_parent << "  : " << condit_prob[index_true] << std::endl;
              }
            }
          }
#endif

        }//end loop over sectors

#ifdef SGNT_DEBUG
        std::cout << "Generating Parameter File" << std::endl;
#endif

        //store the model in text_file
        String file_name = getStringOption_("out_network_model");

        TextFile text_file;

        //store the set of selected IonTypes
        text_file.push_back("<IonTypes>\n");
        for (std::vector<IonType>::iterator m_it = ion_types.begin(); m_it != ion_types.end(); ++m_it)
        {
          text_file.push_back(m_it->residue);
          text_file.push_back(m_it->loss.getString());
          text_file.push_back(m_it->charge);
        }

        //store the number of intensity levels
        text_file.push_back("<IntensityLevels>");
        text_file.push_back(number_of_intensity_levels);

        //store the number of sectors
        text_file.push_back("<Sectors>");
        text_file.push_back(number_of_sectors);

        //store the tree structures for each sector as is_parent_of relations
        text_file.push_back("<TreeStructures>");
        for (Size sector = 0; sector < number_of_sectors; ++sector)
        {
          text_file.push_back("<Sector " + String(sector) + ">");
          for (std::vector<Int>::iterator it = has_parent_all_sectors[sector].begin(); it != has_parent_all_sectors[sector].end(); ++it)
          {
            text_file.push_back(*it);
          }
        }

        //store the conditional probabilities for each sector
        text_file.push_back("<ConditionalProbabilities>");
        for (Size sector = 0; sector < number_of_sectors; ++sector)
        {
          text_file.push_back("<Sector " + String(sector) + ">");
          for (Size i = 0; i < condit_prob_all_sectors[sector].size(); ++i)
          {
            text_file.push_back(condit_prob_all_sectors[sector][i]);
          }
        }

        text_file.store(file_name);
      }

      UInt calc_index(const UInt &type_id, const UInt &intensity_level, const UInt &number_intensity_levels) const
      {
        return type_id * number_intensity_levels + intensity_level;
      }

      void normalizeIntensity(MSSpectrum<RichPeak1D> &S, const UInt &number_of_intensity_levels, const DRealVec &intensity_limits) const
      {
        //see PepNovo Paper
        DoubleReal baseline_grass_intens(0.0), total_intens(0.0);
        const Size weak_third = S.size() / 3;

        //compute baseline_grass_intensity
        S.sortByIntensity();

        MSSpectrum<RichPeak1D>::iterator fwit = S.begin();
        for (UInt idx = 0; idx <= weak_third; ++idx)
        {
          total_intens += fwit->getIntensity();
          fwit++;
        }

        //the average intensity of the weakest 33% peaks
        baseline_grass_intens = total_intens / weak_third;
        //normalization and discretization of intensity
        MSSpectrum<RichPeak1D>::iterator fwit_disc = S.begin();

        for (fwit = S.begin(); fwit != S.end(); ++fwit)
        {
          DoubleReal orig_intens = fwit->getIntensity();
          //normalize
          fwit->setIntensity(orig_intens / baseline_grass_intens);
          //DoubleReal norm_intens = fwit->getIntensity();
          //discretize
          UInt level = 0;
          //if an indexing error occurs here then intensity is +Infinity
          while (level < number_of_intensity_levels - 1 && fwit->getIntensity() >= intensity_limits[level])
          {
            ++level;
          }
          fwit->setIntensity(level);
#ifdef SGNT_DEBUG
          std::cout << "original: " << orig_intens << "   normalized: " << norm_intens << "  discretized:  " << fwit->getIntensity() << std::endl;
#endif
        }
        S.sortByPosition();
      }

      ExitCodes main_(int, const char**)
      {
        trainModel();
        return EXECUTION_OK;
      }
  };


int main(int argc, const char** argv)
{
  SpectrumGeneratorNetworkTrainer tool;
  return tool.main(argc, argv);
}

/// @endcond
