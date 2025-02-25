// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/NUXL/NuXLPresets.h>
#include <OpenMS/SYSTEM/File.h>
#include <nlohmann/json.hpp>
#include <fstream>

using json = nlohmann::json;

namespace OpenMS
{
  namespace NuXLPresets
  {
    void getAllPresetsNames(StringList& built_in, StringList& custom)
    {
      // Start with built-in presets
      built_in = StringList(presets_names.begin(), presets_names.end());
      
      // Try to load custom presets from JSON file
      String exec_path = File::getExecutablePath();
      String json_path = File::path(exec_path) + "/nuxl_presets.json";
      
      if (File::exists(json_path))
      {
        OPENMS_LOG_INFO << "Found custom presets file: " << json_path << std::endl;
      }      

      if (File::exists(json_path))
      {
        try
        {
          std::ifstream file(json_path.c_str());
          json j;
          file >> j;
          
          // Add custom presets to the list
          for (auto it = j.begin(); it != j.end(); ++it)
          {
            custom.push_back(it.key());
            OPENMS_LOG_INFO << "Found custom preset: " << it.key() << std::endl;
          }
        }
        catch (const std::exception& e)
        {
          OPENMS_LOG_WARN << "Error reading custom presets from " << json_path << ": " << e.what() << std::endl;
        }
      }      
      return;
    }
    
    void getPresets(const String& p, 
      StringList& nucleotides, 
      StringList& mapping, 
      StringList& modifications, 
      StringList& fragment_adducts, 
      String& can_cross_link)
    {

      StringList built_in;
      StringList custom;
      getAllPresetsNames(built_in, custom);

      // sanity check: preset name needs to be in the list of supported presets
      bool found_in_built_in = find(built_in.begin(), built_in.end(), p) != built_in.end();
      bool found_in_custom = find(custom.begin(), custom.end(), p) != custom.end();
      if (!found_in_built_in && !found_in_custom)
      {
        throw std::runtime_error("Error: unknown preset.");
      }

      if (found_in_built_in)
      {
        // set NTs for RNA / DNA
        if (p.hasPrefix("RNA"))
        {
          nucleotides = StringList(RNA_nucleotides.begin(), RNA_nucleotides.end());
          mapping = StringList(RNA_mapping.begin(), RNA_mapping.end());
        }
        else if (p.hasPrefix("DNA"))
        {
          nucleotides = StringList(DNA_nucleotides.begin(), DNA_nucleotides.end());
          mapping = StringList(DNA_mapping.begin(), DNA_mapping.end());
        }

        // initialize all StringLists from constexpr arrays
        // note: we do this here as this raises a logic error if e.g., size of the array doesn't match the reserved size.
        //       This can easily happen if a comma is omitted and two string literals on two lines joined
        StringList RNA_UV_modifications(modifications_RNA_UV.begin(), modifications_RNA_UV.end());
        StringList RNA_UV_EXTENDED_modifications(modifications_RNA_UV_EXTENDED.begin(), modifications_RNA_UV_EXTENDED.end());
        StringList RNA_UV_fragments(fragments_RNA_UV.begin(), fragments_RNA_UV.end());

        StringList RNA_UV_4SU_modifications(modifications_RNA_UV_4SU.begin(), modifications_RNA_UV_4SU.end());
        StringList RNA_UV_4SU_EXTENDED_modifications(modifications_RNA_UV_4SU_EXTENDED.begin(), modifications_RNA_UV_4SU_EXTENDED.end());
        StringList RNA_UV_4SU_fragments(fragments_RNA_UV_4SU.begin(), fragments_RNA_UV_4SU.end());

        StringList RNA_UV_6SG_modifications(modifications_RNA_UV_6SG.begin(), modifications_RNA_UV_6SG.end());
        StringList RNA_UV_6SG_EXTENDED_modifications(modifications_RNA_UV_6SG_EXTENDED.begin(), modifications_RNA_UV_6SG_EXTENDED.end());
        StringList RNA_UV_6SG_fragments(fragments_RNA_UV_6SG.begin(), fragments_RNA_UV_6SG.end());

        StringList DNA_UV_modifications(modifications_DNA_UV.begin(), modifications_DNA_UV.end());
        StringList DNA_UV_EXTENDED_modifications(modifications_DNA_UV_EXTENDED.begin(), modifications_DNA_UV_EXTENDED.end());
        StringList DNA_UV_fragments(fragments_DNA_UV.begin(), fragments_DNA_UV.end());

        StringList RNA_DEB_modifications(modifications_RNA_DEB.begin(), modifications_RNA_DEB.end());
        StringList RNA_DEB_EXTENDED_modifications(modifications_RNA_DEB_EXTENDED.begin(), modifications_RNA_DEB_EXTENDED.end());
        StringList RNA_DEB_fragments(fragments_RNA_DEB.begin(), fragments_RNA_DEB.end());

        StringList RNA_NM_modifications(modifications_RNA_NM.begin(), modifications_RNA_NM.end());
        StringList RNA_NM_EXTENDED_modifications(modifications_RNA_NM_EXTENDED.begin(), modifications_RNA_NM_EXTENDED.end());
        StringList RNA_NM_EXTENDED_H2O_modifications(modifications_RNA_NM_EXTENDED_H2O.begin(), modifications_RNA_NM_EXTENDED_H2O.end());

        StringList RNA_NM_fragments(fragments_RNA_NM.begin(), fragments_RNA_NM.end()); 
        StringList RNA_NM_fragments_H2O(fragments_RNA_NM_H2O.begin(), fragments_RNA_NM_H2O.end()); 

        StringList DNA_DEB_modifications(modifications_DNA_DEB.begin(), modifications_DNA_DEB.end());
        StringList DNA_DEB_EXTENDED_modifications(modifications_DNA_DEB_EXTENDED.begin(), modifications_DNA_DEB_EXTENDED.end());
        StringList DNA_DEB_fragments(fragments_DNA_DEB.begin(), fragments_DNA_DEB.end());

        StringList DNA_NM_modifications(modifications_DNA_NM.begin(), modifications_DNA_NM.end());
        StringList DNA_NM_EXTENDED_modifications(modifications_DNA_NM_EXTENDED.begin(), modifications_DNA_NM_EXTENDED.end());
        StringList DNA_NM_fragments(fragments_DNA_NM.begin(), fragments_DNA_NM.end());
        
        StringList RNA_FA_modifications(modifications_RNA_FA.begin(), modifications_RNA_FA.end());
        StringList RNA_FA_fragments(fragments_RNA_FA.begin(), fragments_RNA_FA.end());
        StringList RNA_FA_EXTENDED_modifications(modifications_RNA_FA_EXTENDED.begin(), modifications_RNA_FA_EXTENDED.end());

        StringList DNA_FA_modifications(modifications_DNA_FA.begin(), modifications_DNA_FA.end());
        StringList DNA_FA_fragments(fragments_DNA_FA.begin(), fragments_DNA_FA.end());
        StringList DNA_FA_EXTENDED_modifications(modifications_DNA_FA_EXTENDED.begin(), modifications_DNA_FA_EXTENDED.end());
      
        StringList DNA_UV_BrU_modifications(modifications_DNA_BrU_UV.begin(), modifications_DNA_BrU_UV.end());
        StringList DNA_UV_BrU_fragments(fragments_DNA_BrU_UV.begin(), fragments_DNA_BrU_UV.end());

        const String RNA_U = "U";
        const String RNA_UCGA = "UCGA";
        const String DNA_TCGAd = "TCGAd";
        const String RNA_CGA = "CGA";
        const String DNA_CGAd = "CGAd";

        // set precursor + fragment adducts and cross-linked nucleotide
        if (p == "RNA-UV (U)" || p  == "RNA-UV (UCGA)")
        {
          modifications = RNA_UV_modifications;
          fragment_adducts = RNA_UV_fragments;
          can_cross_link = (p == "RNA-UV (U)") ? RNA_U : RNA_UCGA;
          return;
        }
        else if (p == "RNA-UV Extended (U)" || p  == "RNA-UV Extended (UCGA)")
        {
          modifications = RNA_UV_EXTENDED_modifications; 
          fragment_adducts = RNA_UV_fragments;
          can_cross_link = (p == "RNA-UV Extended (U)") ? RNA_U : RNA_UCGA ;
          return;
        }
        else if (p == "RNA-UV (4SU)")
        {
          nucleotides.push_back("S=C9H13N2O8PS"); // include 4-Thio-UMP
          mapping.push_back("S->S");
          modifications = RNA_UV_4SU_modifications;
          fragment_adducts = RNA_UV_4SU_fragments;
          can_cross_link = "S";
          return;
        }
        else if (p == "RNA-UV Extended (4SU)")
        {
          nucleotides.push_back("S=C9H13N2O8PS"); // include 4-Thio-UMP
          mapping.push_back("S->S");
          modifications = RNA_UV_4SU_EXTENDED_modifications;
          fragment_adducts = RNA_UV_4SU_fragments;
          can_cross_link = "S";
          return;
        }
        else if (p == "RNA-UV (6SG)")
        {
          nucleotides.push_back("X=C10H14N5O7PS"); // include 6-Thio-GMP
          mapping.push_back("X->X");
          modifications = RNA_UV_6SG_modifications;
          fragment_adducts = RNA_UV_6SG_fragments;
          can_cross_link = "X";
          return;
        }
        else if (p == "RNA-UV Extended (6SG)")
        {
          nucleotides.push_back("X=C10H14N5O7PS"); // include 6-Thio-GMP
          mapping.push_back("X->X");
          modifications = RNA_UV_6SG_EXTENDED_modifications;
          fragment_adducts = RNA_UV_6SG_fragments;
          can_cross_link = "X";
          return;
        }    
        else if (p == "DNA-UV")
        {
          modifications = DNA_UV_modifications;
          fragment_adducts = DNA_UV_fragments;
          can_cross_link = DNA_TCGAd;
          return;
        }
        else if (p == "DNA-UV Extended")
        {
          modifications = DNA_UV_EXTENDED_modifications;
          fragment_adducts = DNA_UV_fragments;
          can_cross_link = DNA_TCGAd;
          return;
        }
        else if (p == "DNA-UV (BrU)")
        {
          nucleotides.push_back("B=C9H11N2O8P"); // include BrU
          mapping.push_back("B->B");
          modifications = DNA_UV_BrU_modifications;
          fragment_adducts = DNA_UV_BrU_fragments;
          can_cross_link = "B";
          return;
        }    
        else if (p == "RNA-FA")
        {
          modifications = RNA_FA_modifications;     
          fragment_adducts = RNA_FA_fragments;
          can_cross_link = RNA_CGA;
          return;
        }
        else if (p == "RNA-FA Extended")
        {
          modifications = RNA_FA_EXTENDED_modifications;     
          fragment_adducts = RNA_FA_fragments;
          can_cross_link = RNA_CGA;
          return;
        }
        else if (p == "DNA-FA")
        {
          modifications = DNA_FA_modifications;     
          fragment_adducts = DNA_FA_fragments;
          can_cross_link = DNA_CGAd;
          return;
        }
        else if (p == "DNA-FA Extended")
        {
          modifications = DNA_FA_EXTENDED_modifications;     
          fragment_adducts = DNA_FA_fragments;
          can_cross_link = DNA_CGAd;
          return;
        }
        else if (p == "RNA-DEB")
        {
          // add special methionine loss
          auto r_ptr = const_cast<Residue*>(ResidueDB::getInstance()->getResidue('M'));
          r_ptr->addLossFormula(EmpiricalFormula("CH4S1"));

          modifications = RNA_DEB_modifications;     
          fragment_adducts = RNA_DEB_fragments;
          can_cross_link = RNA_UCGA;
          return;
        }
        else if (p == "RNA-DEB Extended")
        {
          // add special methionine loss
          auto r_ptr = const_cast<Residue*>(ResidueDB::getInstance()->getResidue('M'));
          r_ptr->addLossFormula(EmpiricalFormula("CH4S1"));

          modifications = RNA_DEB_EXTENDED_modifications;
          fragment_adducts = RNA_DEB_fragments;
          can_cross_link = RNA_UCGA;
          return;
        }
        else if (p == "RNA-NM")
        {
          // add special methionine loss
          auto r_ptr = const_cast<Residue*>(ResidueDB::getInstance()->getResidue('M'));
          r_ptr->addLossFormula(EmpiricalFormula("CH4S1"));

          modifications = RNA_NM_modifications;
          fragment_adducts = RNA_NM_fragments; 
          can_cross_link = RNA_UCGA;
          return;
        }
        else if (p == "RNA-NM Extended")
        {
          // add special methionine loss
          auto r_ptr = const_cast<Residue*>(ResidueDB::getInstance()->getResidue('M'));
          r_ptr->addLossFormula(EmpiricalFormula("CH4S1"));

          modifications = RNA_NM_EXTENDED_modifications;
          fragment_adducts = RNA_NM_fragments; 
          can_cross_link = RNA_UCGA;
          return;
        }
        else if (p == "RNA-NM Extended (+H2O)")
        {
          // add special methionine loss
          auto r_ptr = const_cast<Residue*>(ResidueDB::getInstance()->getResidue('M'));
          r_ptr->addLossFormula(EmpiricalFormula("CH4S1"));

          modifications = RNA_NM_EXTENDED_H2O_modifications;
          fragment_adducts = RNA_NM_fragments_H2O; 
          can_cross_link = RNA_UCGA;
          return;
        }
        else if (p == "DNA-DEB")
        {
          modifications = DNA_DEB_modifications;
          fragment_adducts = DNA_DEB_fragments;
          can_cross_link = DNA_TCGAd;
          return;
        }
        else if (p == "DNA-DEB Extended")
        {
          modifications = DNA_DEB_EXTENDED_modifications;
          fragment_adducts = DNA_DEB_fragments;
          can_cross_link = DNA_TCGAd;
          return;
        }    
        else if (p == "DNA-NM")
        {
          modifications = DNA_NM_modifications;
          fragment_adducts = DNA_NM_fragments;
          can_cross_link = DNA_TCGAd;
          return;
        }
        else if (p == "DNA-NM Extended")
        {
          modifications = DNA_NM_EXTENDED_modifications;
          fragment_adducts = DNA_NM_fragments;
          can_cross_link = DNA_TCGAd;
          return;
        }      
      }
      else // found in custom (json)
      {
        // Try to load custom presets from JSON file
        String exec_path = File::getExecutablePath();
        String json_path = File::path(exec_path) + "/nuxl_presets.json";
        
        if (File::exists(json_path))
        {
          try
          {
            std::ifstream file(json_path.c_str());
            json j;
            file >> j;
            
            // Check if the requested preset exists in the JSON file
            if (j.contains(p.c_str()))
            {
              const auto& preset = j[p.c_str()];
              
              // Load nucleotides
              if (preset.contains("target_nucleotides"))
              {
                nucleotides.clear();
                for (const auto& nuc : preset["target_nucleotides"])
                {
                  nucleotides.push_back(nuc.get<std::string>());
                }
              }
              
              // Load mapping
              if (preset.contains("mapping"))
              {
                mapping.clear();
                for (const auto& map : preset["mapping"])
                {
                  mapping.push_back(map.get<std::string>());
                }
              }
              
              // Load modifications
              if (preset.contains("modifications"))
              {
                modifications.clear();
                for (const auto& mod : preset["modifications"])
                {
                  modifications.push_back(mod.get<std::string>());
                }
              }
              
              // Load fragment adducts
              if (preset.contains("fragment_adducts"))
              {
                fragment_adducts.clear();
                for (const auto& frag : preset["fragment_adducts"])
                {
                  fragment_adducts.push_back(frag.get<std::string>());
                }
              }
              
              // Load can_cross_link
              if (preset.contains("can_cross_link"))
              {
                can_cross_link = preset["can_cross_link"].get<std::string>();
              }
              
              // Custom preset loaded successfully, return
              OPENMS_LOG_INFO << "Using custom preset '" << p << "' from " << json_path << std::endl;
              return;
            }
          }
          catch (const std::exception& e)
          {
            // If there's an error reading the JSON file, fall back to built-in presets
            OPENMS_LOG_WARN << "Error reading custom presets from " << json_path << ": " << e.what() << std::endl;
          }
        }      
      } // end of: found in custom

    }
  }
}