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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_MULTIPLEXDELTAMASSES_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_MULTIPLEXDELTAMASSES_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/LogStream.h>

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
      DeltaMass(double dm, String l);
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

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_MULTIPLEXDELTAMASSES_H

