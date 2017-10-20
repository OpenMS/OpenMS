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
// $Maintainer:	Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_EXPERIMENTALDESIGN_H
#define OPENMS_KERNEL_EXPERIMENTALDESIGN_H

#include <OpenMS/DATASTRUCTURES/String.h>

#include <vector>
#include <string>
#include <map>
#include <set>

namespace OpenMS
{
  /**
  @brief A TSV and user friendly representation of the experimental design.
  Used for loading and storing the experimental design in OpenMS.

  @ingroup Metadata
  */
  class OPENMS_DLLAPI ExperimentalDesign
  {
  public:
    /// 1) Mandatory section with run-level information of the experimental design.
    ///    Required to process fractionated data and technical replicates.
    ///
    /// Run	Spectra File	Fraction	Technical Replicate   
    ///   1	humanA.mzML	      1	              1          
    ///   2	humanB.mzML	      1	              1          
    ///   3	humanC.mzML	      2	              1              
    ///   4	humanD.mzML	      2	              1
    ///   5	humanE.mzML	      3	              1
    ///   6	humanF.mzML	      3	              1          
    ///   7	humanG.mzML	      3	              2       (<- example how a 2nd technical replicate is stored for fraction 3)
    /// TODO: add possibility to provide optional columns that map to additional run-level specific meta data (see, for example, the mzTab specification)
    class OPENMS_DLLAPI MSRun
    {
    public:
      MSRun():
        file("NA"), 
        fraction(1), 
        technical_replicate(1) 
      {
      }
      std::string file; //< file name, mandatory
      unsigned fraction; //< fraction 1..m, mandatory, 1 if not set
      unsigned technical_replicate; //< technical replicate 1..k of a fraction, 1 if not set
    };
    std::vector<MSRun> runs;  //< run 1..n (index + 1 determines run id of the first column)

    /// return fraction index to MSRuns (e.g., Fraction 1 maps to MSRun 1 and 2 in the example above)
    std::map<unsigned int, std::set<unsigned int> > getFractionToRunsMapping() const;

    /// return if each fraction number is associated with the same number of runs 
    bool sameNrOfRunsPerFraction() const;

    /// TODO:
    /// 2) Optional section with assay-level information of the experimental design.
    ///    Required to perform statistical down-stream processing.
    ///
    /// Assay	Spectra File	Assay Quantification	Numeric Factor "concentration [mmol/l]"
    ///  1	  humanA.mzML	    iTRAQ reagent 114	                 0.434
    ///  2   	humanA.mzML	    iTRAQ reagent 115	               223.34
    ///  3  	humanA.mzML	    iTRAQ reagent 116	               984.52

    /// Loads an experimental design from a tabular separated file
    void load(const String & tsv_file, ExperimentalDesign & design) const;

  };
}


#endif // header guard
