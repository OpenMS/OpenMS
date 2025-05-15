// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <vector>
#include <set>
#include <algorithm>
#include <iostream>

namespace OpenMS
{
  /**
   * @brief data structure for mass shift pattern
   * 
   * Groups of labelled peptides appear with characteristic mass shifts.
   * 
   * For example, for an Arg6 labeled SILAC peptide pair we expect to see
   * mass shifts of 0 and 6 Da. Or as second example, for a 
   * peptide pair of a dimethyl labelled sample with a single lysine
   * we will see mass shifts of 56 Da and 64 Da.
   * 28 Da (N-term) + 28 Da (K) and 34 Da (N-term) + 34 Da (K)
   * for light and heavy partners respectively.
   * 
   * The data structure stores the mass shifts and corresponding labels
   * for a group of matching peptide features. 
   */
  class OPENMS_DLLAPI MultiplexDeltaMasses
  {
    public:
    
    /**
     * @brief set of labels associated with a mass shift
     * 
     * For example, a set of SILAC labels [Lys8, Lys8, Arg10] would
     * result in a +26 Da mass shift.
     */
    typedef std::multiset<String> LabelSet;

    /**
     * @brief mass shift with corresponding label set
     */
    struct OPENMS_DLLAPI DeltaMass
    {
      double delta_mass;
      LabelSet label_set;
      
      DeltaMass(double dm, LabelSet ls);
      
      // delta mass with a label set containing a single label
      DeltaMass(double dm, const String& l);
    };

    /**
     * @brief constructor
     */
    MultiplexDeltaMasses();
    
    /**
     * @brief constructor
     */
    MultiplexDeltaMasses(const std::vector<DeltaMass>& dm);
        
    /**
     * @brief returns delta masses
     */
    std::vector<DeltaMass>& getDeltaMasses();
    
    /**
     * @brief returns delta masses
     */
    const std::vector<DeltaMass>& getDeltaMasses() const;
    
    /**
     * @brief converts a label set to a string
     */
    static String labelSetToString(const LabelSet& ls);
    
    private:
   
    /**
     * @brief mass shifts between peptides
     * (including zero mass shift for first peptide)
     */
    std::vector<DeltaMass> delta_masses_;
    
 };
 
 bool operator<(const MultiplexDeltaMasses &dm1, const MultiplexDeltaMasses &dm2);
  
}


