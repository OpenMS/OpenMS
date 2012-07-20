// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Alexandra Zerck $
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/TARGETED/OfflinePrecursorIonSelection.h>

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/ANALYSIS/TARGETED/PrecursorIonSelectionPreprocessing.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>

namespace OpenMS
{

  OfflinePrecursorIonSelection::OfflinePrecursorIonSelection() :
    DefaultParamHandler("OfflinePrecursorIonSelection")
  {
    defaults_.setValue("ms2_spectra_per_rt_bin", 5, "Number of allowed MS/MS spectra in a retention time bin.");
    defaults_.setMinInt("ms2_spectra_per_rt_bin", 1);

    defaults_.setValue("min_peak_distance", 3., "The minimal distance (in Da) of two peaks in one spectrum so that they can be selected.");
    defaults_.setMinFloat("min_peak_distance", 0.);

    defaults_.setValue("selection_window", 2., "All peaks within a mass window (in Da) of a selected peak are also selected for fragmentation.");
    defaults_.setMinFloat("selection_window", 0.);

    defaults_.setValue("exclude_overlapping_peaks", "false", "If true overlapping or nearby peaks (within min_peak_distance) are excluded for selection.");
    defaults_.setValidStrings("exclude_overlapping_peaks", StringList::create("true,false"));

    defaults_.setValue("Exclusion:use_dynamic_exclusion", "false", "If true dynamic exclusion is applied.");
    defaults_.setValidStrings("Exclusion:use_dynamic_exclusion", StringList::create("true,false"));

    defaults_.setValue("Exclusion:exclusion_time", 100., "The time (in seconds) a feature is excluded.");
    defaults_.setMinFloat("Exclusion:exclusion_time", 0.);

    defaults_.insert("ProteinBasedInclusion:",PSLPFormulation().getDefaults());
    defaults_.remove("ProteinBasedInclusion:mz_tolerance");
    defaults_.remove("ProteinBasedInclusion:combined_ilp:");
    defaults_.remove("ProteinBasedInclusion:thresholds:min_protein_probability");
    defaults_.remove("ProteinBasedInclusion:thresholds:min_pred_pep_prob");
    defaults_.remove("ProteinBasedInclusion:thresholds:min_rt_weight");
    defaults_.setValue("ProteinBasedInclusion:max_list_size",1000,"The maximal number of precursors in the inclusion list.");
    defaults_.setMinInt("ProteinBasedInclusion:max_list_size",1);

    defaultsToParam_();
  }

  OfflinePrecursorIonSelection::~OfflinePrecursorIonSelection()
  {

  }

  void OfflinePrecursorIonSelection::createProteinSequenceBasedLPInclusionList(String include, String rt_model_file, String pt_model_file,
                                                                               FeatureMap<> & precursors)
  {
    PrecursorIonSelectionPreprocessing pisp;
    Param pisp_param = pisp.getParameters();
    pisp_param.setValue("store_peptide_sequences", "true");
    pisp.setParameters(pisp_param);
    pisp.dbPreprocessing(include, rt_model_file, pt_model_file, false);
    //  std::cout << "now learn rt probabilities"<<std::endl;
    //pisp.learnRTProbabilities(f_map,rt_model,0.5);
    //  pisp.setGaussianParameters(3,-1);
    PSLPFormulation ilp_wrapper;
    Param opis_param = param_.copy("ProteinBasedInclusion:", true);
    opis_param.remove("max_list_size");
    ilp_wrapper.setParameters(opis_param);
    ilp_wrapper.setLPSolver(solver_);
    // std::cout << "nun die inclusion liste erstellen"<<std::endl;
    // std::cout << param_.getValue("ms2_spectra_per_rt_bin") <<std::endl;
    // std::cout << param_.getValue("ProteinBasedInclusion:max_list_size") <<std::endl;
    ilp_wrapper.createAndSolveILPForInclusionListCreation(pisp, param_.getValue("ms2_spectra_per_rt_bin"),
                                                          param_.getValue("ProteinBasedInclusion:max_list_size"), precursors, true); //,960.,3840.,30.);

  }

}
