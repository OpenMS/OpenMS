// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/DECHARGING/FeatureDecharger.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
   @page Decharger Decharger

   @brief Decharges a feature map by clustering charge variants of a peptide to zero-charge entities.

   The Decharger uses an hierarchical clustering (complete linkage) to group charge variants of the same peptide, which
   usually occur in ESI ionization mode. The resulting zero-charge peptides, which are defined by RT and mass
   are written to a featureXML file. Intensities of charge variants are summed up. The position of the zero charge
   variant is the average of all clustered peptides in each dimension.
   If several peptides with the same charge variant are grouped (which is clearly not allowed), a heuristic is used:
   <ul>
   <li>cluster consists of only one charge variant (but several peptides) -> split cluster into single elements</li>
   <li>cluster consists of several charge variants -> dispose cluster</li>
   </ul>
   

   @ingroup TOPP
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPDecharger
  : public TOPPBase
{
 public:
  TOPPDecharger()
    : TOPPBase("Decharger","decharging of feature maps")
  {
  }

 protected:
  void registerOptionsAndFlags_()
  {
    registerStringOption_("in","<file>","","input feat file");
    registerStringOption_("out","<file>","","output feat file");

    addEmptyLine_();
    addText_("All other options of the Decharger depend on the FeatureDecharger and HierarchicalClustering used.\n"
             "They can be given only in the 'algorithm' section  of the INI file.\n"
             "For a detailed description, please have a look at the doxygen documentation.\n"
             "How the docu can be built is explained in OpenMS/doc/index.html."); 
    
    registerSubsection_("algorithm");
  }


  ExitCodes main_(int , char**)
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");

    FeatureDecharger fdc;
    Param const& dc_param = getParam_().copy("algorithm:FeatureDecharger:",true);

    writeDebug_("Parameters passed to Decharger", dc_param, 3);
    
    if (dc_param.empty())
    {
      writeLog_("No parameters for Decharger module given. Aborting!");
      return ILLEGAL_PARAMETERS;
    }
    fdc.setParameters(dc_param); 

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    writeDebug_("Loading input file", 1);
    
    typedef FeatureMap<> FeatureMapType;
    FeatureMapType map;
    FeatureXMLFile().load(in,map);
    //map.sortByPosition();

    FeatureMapType map_sm;
    
    // removing items with charge 1
    for (uint i = 0; i<map.size(); ++i)
    {
      if (map[i].getCharge() > 1)
      {
        map_sm.push_back (map[i]);              
      }
    }
  
    std::cout << "removing features with charge 1...  #Features before: " << map.size() << "#Features after: " << map_sm.size() << "\n";

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    fdc.compute(map_sm);
    
    FeatureMapType feature_map = fdc.getFeatureMap();

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    
    writeDebug_("Saving output file", 1);

    FeatureXMLFile().store(out,feature_map);

    return EXECUTION_OK;
  }
};


Param TOPPBase::getSubsectionDefaults_(const String& /*section*/) const
{
  // there is only one subsection: 'algorithm' (s.a) .. and in it belongs the FeatureDecharger param
  FeatureDecharger fdc;
  Param tmp;
  tmp.insert("FeatureDecharger:",fdc.getParameters());
  return tmp;
}

int main( int argc, char ** argv )
{
    TOPPDecharger tool;
    return tool.main(argc,argv);
}

/// @endcond
