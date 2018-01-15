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
// $Maintainer: Douglas McCloskey $
// $Authors: Douglas McCloskey $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_ABSOLUTEQUANTITATIONSTANDARDS_H
#define OPENMS_METADATA_ABSOLUTEQUANTITATIONSTANDARDS_H

#include <OpenMS/config.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <cstddef> // for size_t & ptrdiff_t
#include <vector>
#include <string>

namespace OpenMS
{

  /**
    @brief AbsoluteQuantitationStandards is a class to handle the relationship between
    runs, components, and their actual concentrations.

    A mapping between a run, the components in the run, and the actual concentration
    of the components in the run are required to build a calibration curve that is
    required for absolute quantitation.
  */
  class OPENMS_DLLAPI AbsoluteQuantitationStandards
  {

public:    
    //@{
    /// Constructor
    AbsoluteQuantitationStandards();

    /// Destructor
    ~AbsoluteQuantitationStandards();
    //@}
   
    /**
      @brief Structure to map runs to components to known concentrations

    */ 
    struct runConcentration
    {
      String sample_name;
      String component_name;
      String IS_component_name;
      double actual_concentration;
      double IS_actual_concentration;
      String concentration_units;
      double dilution_factor;
    };

    /**
      @brief Structure to hold a single component and its corresponding known concentration.

    */ 
    struct featureConcentration
    {
      Feature feature;
      Feature IS_feature;
      double actual_concentration;
      double IS_actual_concentration;
      String concentration_units;
      double dilution_factor;
    };
    
     /**
       @brief Method to map runs to components to known concentrations

       Note that for the method to work, the features must be annotated with
         a metaValue for "sample_name"

      @param run_concentrations a list of runConcentration structs (e.g., from file upload).
      @param features a list of corresponding features for each of the unique runs in run_concentrations
      @param components_to_concentrations A map that links run data to feature data
 
     */ 
     void mapConcentrationsToComponents(const std::vector<runConcentration> & run_concentrations,
      const std::vector<FeatureMap> & features,
      std::map<String,std::vector<featureConcentration>> components_to_concentrations);

    // members
    std::map<String, std::vector<featureConcentration>> components_to_concentrations;

  };
}
#endif // OPENMS_METADATA_ABSOLUTEQUANTITATIONSTANDARDS_H

