// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: OpenMS Team $
// $Authors: OpenMS Team $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/OPENSWATH/GenericSequenceHelper.h>
#include <OpenMS/CHEMISTRY/SequenceBase.h>

namespace OpenMS
{
  /**
   * @brief Example demonstrating the CRTP sequence interface in OpenSWATH context
   * 
   * This example shows how the same OpenSWATH workflow functions can process
   * both AASequence (peptides) and NASequence (oligonucleotides) using the
   * common SequenceBase interface.
   */
  class OPENMS_DLLAPI SequenceUsageExample
  {
  public:
    
    /**
     * @brief Example OpenSWATH workflow function that works with any sequence type
     * 
     * This template function demonstrates how OpenSWATH algorithms can be written
     * generically to handle both peptides and oligonucleotides.
     */
    template<typename SeqType>
    static void openSwathWorkflow(const SeqType& sequence, const std::string& sequence_id)
    {
      static_assert(is_sequence_type_v<SeqType>, "Must be a sequence type");
      
      std::cout << "\n=== OpenSWATH Workflow for " << sequence_id << " ===" << std::endl;
      
      // Step 1: Basic sequence validation
      if (sequence.empty())
      {
        std::cout << "ERROR: Empty sequence!" << std::endl;
        return;
      }
      
      // Step 2: Compute basic properties (works for both types)
      std::cout << "Sequence: " << sequence.toString().toStdString() << std::endl;
      std::cout << "Length: " << sequence.size() << " residues" << std::endl;
      
      // Step 3: Mass calculations
      auto default_mass = sequence.getMonoWeight();
      std::cout << "Monoisotopic mass: " << default_mass << " Da" << std::endl;
      
      // Step 4: Generate theoretical fragment masses
      std::vector<int> charges;
      if constexpr (std::is_same_v<SeqType, AASequence>)
      {
        charges = {1, 2, 3}; // Positive charges for peptides
        std::cout << "Processing as PEPTIDE sequence" << std::endl;
      }
      else
      {
        charges = {-1, -2, -3}; // Negative charges for oligonucleotides
        std::cout << "Processing as OLIGONUCLEOTIDE sequence" << std::endl;
      }
      
      auto theoretical_masses = GenericSequenceHelper::computeTheoreticalMasses(sequence, charges);
      
      // Step 5: Fragment series generation
      auto fragments = GenericSequenceHelper::generateFragmentSeries(sequence, true, true, charges[0]);
      std::cout << "Generated " << fragments.size() << " fragment ions" << std::endl;
      
      // Step 6: Convert to LightCompound for OpenSWATH processing
      auto light_compound = GenericSequenceHelper::sequenceToLightCompound(sequence, sequence_id, charges[0]);
      std::cout << "Created LightCompound with ID: " << light_compound.id << std::endl;
      
      // Step 7: Type-specific processing
      if constexpr (std::is_same_v<SeqType, AASequence>)
      {
        processPeptideSpecific(sequence, light_compound);
      }
      else if constexpr (std::is_same_v<SeqType, NASequence>)
      {
        processOligonucleotideSpecific(sequence, light_compound);
      }
      
      std::cout << "=== Workflow complete ===" << std::endl;
    }
    
    /**
     * @brief Example of running both peptide and oligonucleotide workflows
     */
    static void runExamples()
    {
      std::cout << "CRTP Sequence Interface Example" << std::endl;
      std::cout << "===============================" << std::endl;
      
      try
      {
        // Create example sequences
        AASequence peptide = AASequence::fromString("PEPTIDER");
        NASequence oligo = NASequence::fromString("AUCG");
        
        // Run the same workflow for both types
        openSwathWorkflow(peptide, "example_peptide");
        openSwathWorkflow(oligo, "example_oligonucleotide");
        
        // Demonstrate scoring
        std::cout << "\n=== Scoring Example ===" << std::endl;
        
        // Simulate experimental data
        std::vector<std::pair<double, double>> experimental_peaks = {
          {100.0, 1000.0}, {200.0, 500.0}, {300.0, 2000.0}, {400.0, 800.0}
        };
        
        auto peptide_score = GenericSequenceHelper::scoreSequenceMatch(peptide, experimental_peaks, 0.1);
        auto oligo_score = GenericSequenceHelper::scoreSequenceMatch(oligo, experimental_peaks, 0.1);
        
        std::cout << "Peptide score: " << peptide_score << std::endl;
        std::cout << "Oligonucleotide score: " << oligo_score << std::endl;
      }
      catch (const std::exception& e)
      {
        std::cout << "Error in example: " << e.what() << std::endl;
      }
    }
    
  private:
    
    /**
     * @brief Peptide-specific processing
     */
    static void processPeptideSpecific(const AASequence& peptide, const OpenSwath::LightCompound& compound)
    {
      std::cout << "  -> Peptide-specific processing:" << std::endl;
      std::cout << "     - Could generate b/y ion series" << std::endl;
      std::cout << "     - Could check for modifications" << std::endl;
      std::cout << "     - Could validate amino acid sequence" << std::endl;
      std::cout << "     - Peptide group label: " << compound.peptide_group_label << std::endl;
    }
    
    /**
     * @brief Oligonucleotide-specific processing
     */
    static void processOligonucleotideSpecific(const NASequence& oligo, const OpenSwath::LightCompound& compound)
    {
      std::cout << "  -> Oligonucleotide-specific processing:" << std::endl;
      std::cout << "     - Could generate a/b/c/d/w/x/y/z ion series" << std::endl;
      std::cout << "     - Could check for nucleotide modifications" << std::endl;
      std::cout << "     - Could analyze secondary structure" << std::endl;
      std::cout << "     - Compound name: " << compound.compound_name << std::endl;
    }
  };

  /**
   * @brief Template function showing how OpenSWATH functions can be generalized
   * 
   * This is an example of how existing OpenSWATH functions that currently work
   * only with AASequence can be templated to work with both sequence types.
   */
  template<typename SeqType>
  void processSequenceForDIA(const SeqType& sequence)
  {
    static_assert(is_sequence_type_v<SeqType>, "Must be a sequence type");
    
    // Common DIA processing
    if (sequence.empty()) return;
    
    auto mass = sequence.getMonoWeight();
    auto formula = sequence.getFormula();
    
    // Generate theoretical spectrum
    auto fragments = GenericSequenceHelper::generateFragmentSeries(sequence);
    
    // Type-specific DIA scoring
    if constexpr (std::is_same_v<SeqType, AASequence>)
    {
      // Call existing peptide-specific DIA functions
      // e.g., dia_by_ion_score(), etc.
    }
    else if constexpr (std::is_same_v<SeqType, NASequence>)
    {
      // Call new oligonucleotide-specific DIA functions
      // e.g., dia_nucleotide_ion_score(), etc.
    }
  }

} // namespace OpenMS