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
    StringList getAllPresetsNames()
    {
      StringList presets;
      
      // Try to load presets from JSON file
      String json_path = File::getOpenMSDataPath() + "/NUXL/nuxl_presets.json";
      
      if (File::exists(json_path))
      {
        OPENMS_LOG_INFO << "Found presets file: " << json_path << std::endl;
      }      

      if (File::exists(json_path))
      {
        try
        {
          std::ifstream file(json_path.c_str());
          json j;
          file >> j;
          
          // Add all presets to the list
          for (auto it = j.begin(); it != j.end(); ++it)
          {
            presets.push_back(it.key());
            OPENMS_LOG_INFO << "Found preset: " << it.key() << std::endl;
          }
        }
        catch (const std::exception& e)
        {
          OPENMS_LOG_WARN << "Error reading presets from " << json_path << ": " << e.what() << std::endl;
        }
      }
      else
      {
        OPENMS_LOG_WARN << "Presets file not found: " << json_path << std::endl;
      }
      return presets;
    }
    
    void getPresets(const String& p, 
      StringList& nucleotides, 
      StringList& mapping, 
      StringList& modifications, 
      StringList& fragment_adducts, 
      String& can_cross_link)
    {
      StringList presets = getAllPresetsNames();

      // Check if preset exists
      bool found = find(presets.begin(), presets.end(), p) != presets.end();
      if (!found)
      {
        throw std::runtime_error("Error: unknown preset '" + p + "'.");
      }

      // Try to load presets from JSON file
      String share_path = File::getOpenMSDataPath();
      String json_path = share_path + "/NUXL/nuxl_presets.json";
      
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
            
            // Special handling for DEB and NM presets that need methionine loss
            if (p.hasSubstring("DEB") || p.hasSubstring("NM"))
            {
              // add special methionine loss
              auto r_ptr = const_cast<Residue*>(ResidueDB::getInstance()->getResidue('M'));
              r_ptr->addLossFormula(EmpiricalFormula("CH4S1"));
            }
            
            // Preset loaded successfully, return
            OPENMS_LOG_INFO << "Using preset '" << p << "' from " << json_path << std::endl;
            return;
          }
        }
        catch (const std::exception& e)
        {
          // If there's an error reading the JSON file, throw an error
          OPENMS_LOG_WARN << "Error reading presets from " << json_path << ": " << e.what() << std::endl;
          throw std::runtime_error("Error reading presets.");
        }
      }
      else
      {
        throw std::runtime_error("Error: presets file not found.");
      }
    }
  }
}