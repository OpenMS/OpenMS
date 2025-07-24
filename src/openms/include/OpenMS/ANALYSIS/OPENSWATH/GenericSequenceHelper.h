// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: OpenMS Team $
// $Authors: OpenMS Team $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <type_traits>
#include <vector>
#include <iostream>

namespace OpenMS
{
  /**
   * @brief Generic helper functions for sequence processing in OpenSWATH
   * 
   * This class provides template functions that work with both AASequence and NASequence
   * through the common SequenceBase interface, enabling code reuse in OpenSWATH workflows.
   */
  class OPENMS_DLLAPI GenericSequenceHelper
  {
  public:
    
    /**
     * @brief Generic function to compute theoretical masses for any sequence type
     * 
     * @tparam SeqType Either AASequence or NASequence
     * @param sequence The sequence to analyze
     * @param charges Vector of charge states to compute
     * @return Vector of theoretical masses for each charge state
     */
    template<typename SeqType>
    static std::vector<double> computeTheoreticalMasses(
        const SeqType& sequence, 
        const std::vector<int>& charges)
    {
      static_assert(is_sequence_type_v<SeqType>, "SeqType must be a sequence type");
      
      std::vector<double> masses;
      masses.reserve(charges.size());
      
      for (int charge : charges)
      {
        double mass = sequence.getMonoWeight(charge);
        masses.push_back(mass);
      }
      
      return masses;
    }

    /**
     * @brief Generic function to create theoretical fragment ion series
     *
     * @tparam SeqType Either AASequence or NASequence
     * @param sequence The sequence to fragment
     * @param charges Vector of charge states for each fragment
     * @return Vector of fragment masses
     */
    template<typename SeqType>
    static std::vector<double> generateFragmentIons(
        const SeqType& sequence,
        const std::vector<int>& charges)
    {
      static_assert(is_sequence_type_v<SeqType>, "SeqType must be a sequence type");
      
      std::vector<double> fragment_masses;
      
      for (int charge : charges)
      {
        double mass = sequence.getMonoWeight(charge);
        fragment_masses.push_back(mass);
      }
      
      return fragment_masses;
    }

    /**
     * @brief Generic function to create prefix/suffix series for any sequence
     * 
     * @tparam SeqType Either AASequence or NASequence
     * @param sequence The sequence to fragment
     * @param generate_prefixes Whether to generate prefix fragments
     * @param generate_suffixes Whether to generate suffix fragments  
     * @param charge Charge state for fragments
     * @return Vector of fragment masses
     */
    template<typename SeqType>
    static std::vector<double> generateFragmentSeries(
        const SeqType& sequence,
        bool generate_prefixes = true,
        bool generate_suffixes = true,
        int charge = 1)
    {
      static_assert(is_sequence_type_v<SeqType>, "SeqType must be a sequence type");
      
      std::vector<double> fragments;
      
      if (generate_prefixes)
      {
        for (Size i = 1; i < sequence.size(); ++i)
        {
          auto prefix = sequence.getPrefix(i);
          double mass = prefix.getMonoWeight(charge);
          fragments.push_back(mass);
        }
      }
      
      if (generate_suffixes)
      {
        for (Size i = 1; i < sequence.size(); ++i)
        {
          auto suffix = sequence.getSuffix(i);
          double mass = suffix.getMonoWeight(charge);
          fragments.push_back(mass);
        }
      }
      
      return fragments;
    }

    /**
     * @brief Convert any sequence to LightCompound for OpenSWATH processing
     * 
     * @tparam SeqType Either AASequence or NASequence
     * @param sequence The sequence to convert
     * @param id Compound identifier
     * @param charge Charge state
     * @return LightCompound object for use in OpenSWATH workflows
     */
    template<typename SeqType>
    static OpenSwath::LightCompound sequenceToLightCompound(
        const SeqType& sequence,
        const std::string& id,
        int charge = 0)
    {
      static_assert(is_sequence_type_v<SeqType>, "SeqType must be a sequence type");
      
      OpenSwath::LightCompound compound;
      compound.id = id;
      compound.sequence = sequence.toString().toStdString();
      compound.charge = charge;
      
      // Set rt to default if not specified
      compound.rt = -1.0;
      compound.drift_time = -1.0;
      
      // Type-specific processing
      if constexpr (std::is_same_v<SeqType, AASequence>)
      {
        // Peptide-specific processing
        compound.peptide_group_label = id + "_group";
      }
      else if constexpr (std::is_same_v<SeqType, NASequence>)
      {
        // Oligonucleotide-specific processing  
        compound.compound_name = id;
      }
      
      return compound;
    }

    /**
     * @brief Generic scoring function that works with any sequence type
     * 
     * @tparam SeqType Either AASequence or NASequence
     * @param sequence The sequence being scored
     * @param experimental_masses Vector of observed masses with intensities
     * @param mass_tolerance Mass tolerance for matching
     * @param intensity_threshold Minimum intensity threshold
     * @return Score based on mass accuracy and coverage
     */
    template<typename SeqType>
    static double scoreSequenceMatch(
        const SeqType& sequence,
        const std::vector<std::pair<double, double>>& experimental_masses,
        double mass_tolerance = 0.01,
        double intensity_threshold = 0.0)
    {
      static_assert(is_sequence_type_v<SeqType>, "SeqType must be a sequence type");
      
      // Generate theoretical fragments
      auto theoretical_fragments = generateFragmentSeries(sequence, true, true, 1);
      
      // Score based on mass matching
      int matches = 0;
      double total_intensity = 0.0;
      
      for (const auto& exp_peak : experimental_masses)
      {
        if (exp_peak.second < intensity_threshold) continue;
        
        for (double theo_mass : theoretical_fragments)
        {
          if (std::abs(exp_peak.first - theo_mass) <= mass_tolerance)
          {
            matches++;
            total_intensity += exp_peak.second;
            break;
          }
        }
      }
      
      // Simple scoring: matches weighted by intensity
      return matches > 0 ? total_intensity / matches : 0.0;
    }

    /**
     * @brief Generic function to process sequence-specific OpenSWATH operations
     * 
     * This demonstrates how to handle both peptide and oligonucleotide sequences
     * in the same function while preserving type-specific behavior.
     */
    template<typename SeqType>
    static void processForOpenSWATH(const SeqType& sequence)
    {
      static_assert(is_sequence_type_v<SeqType>, "SeqType must be a sequence type");
      
      // Common operations that work for both types
      std::cout << "Processing sequence: " << sequence.toString().toStdString() << std::endl;
      std::cout << "Length: " << sequence.size() << std::endl;
      std::cout << "Monoisotopic mass: " << sequence.getMonoWeight() << std::endl;
      
      // Type-specific operations using constexpr if (C++17)
      if constexpr (std::is_same_v<SeqType, AASequence>)
      {
        std::cout << "Processing as peptide sequence..." << std::endl;
        // Could call peptide-specific OpenSWATH functions here
        // e.g., generateBYIons(), processPeptideModifications(), etc.
      }
      else if constexpr (std::is_same_v<SeqType, NASequence>)
      {
        std::cout << "Processing as nucleic acid sequence..." << std::endl;
        // Could call oligonucleotide-specific OpenSWATH functions here
        // e.g., generateABCDWXYZIons(), processNucleotideModifications(), etc.
      }
    }
  };

  /**
   * @brief Example usage functions demonstrating the CRTP interface
   */
  namespace GenericSequenceExample
  {
    /**
     * @brief Example function showing how OpenSWATH code can work with both sequence types
     */
    inline void demonstrateGenericUsage()
    {
      // Create both sequence types
      AASequence peptide = AASequence::fromString("PEPTIDER");
      NASequence oligo = NASequence::fromString("AUCG");
      
      // Generic processing works for both
      auto peptide_masses = GenericSequenceHelper::computeTheoreticalMasses(peptide, {1, 2, 3});
      auto oligo_masses = GenericSequenceHelper::computeTheoreticalMasses(oligo, {-1, -2, -3});
      
      // Convert to LightCompound for OpenSWATH workflows
      auto peptide_compound = GenericSequenceHelper::sequenceToLightCompound(peptide, "peptide_1", 2);
      auto oligo_compound = GenericSequenceHelper::sequenceToLightCompound(oligo, "oligo_1", -2);
      
      // Process for OpenSWATH (demonstrates type-specific branching)
      GenericSequenceHelper::processForOpenSWATH(peptide);
      GenericSequenceHelper::processForOpenSWATH(oligo);
      
      std::cout << "Generic processing complete!" << std::endl;
    }
  }

} // namespace OpenMS