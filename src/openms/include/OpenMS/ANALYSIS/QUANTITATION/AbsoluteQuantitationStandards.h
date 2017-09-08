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

#ifndef OPENMS_ANALYSIS_OPENSWATH_ABSOLUTEQUANTIFICATIONSTANDARDS_H
#define OPENMS_ANALYSIS_OPENSWATH_ABSOLUTEQUANTIFICATIONSTANDARDS_H

#include <OpenMS/config.h>

#include <cstddef> // for size_t & ptrdiff_t
#include <vector>
#include <string>

namespace OpenMS
{

  /**
    @brief AbsoluteQuantitationStandards is a class to handle the relationship between
    runs, components, and actual concentrations.
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
      @brief Structure to map runs to features to known concentrations

    */ 
    struct runConcentrations
    {
      str::string run_id;
      str::string feature_id;
      double actual_concentration;
      str::string concentration_units;
    }

    /**
      @brief Structure to hold all features for a single component
        with their corresponding known concentrations.

    */ 
    struct featureConcentrations
    {
      std::vector<Feature> features;
      std::vector<double> actual_concentrations;
      std::vector<str::string> concentration_units;
    }
                                      
    // members
    std::map<std::string,featureConcentrations> features_to_oncentrations;

  };

}
#endif // OPENMS_ANALYSIS_OPENSWATH_ABSOLUTEQUANTIFICATIONSTANDARDS_H

