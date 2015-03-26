// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Andreas Bertsch, Marc Sturm, Sven Nahnsen, Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_ID_CONSENSUSIDALGORITHM_H
#define OPENMS_ANALYSIS_ID_CONSENSUSIDALGORITHM_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <map>
#include <vector>

namespace OpenMS
{
  /**
    @brief Base class for all ConsensusID algorithms (that calculate a consensus from multiple ID runs)

    @htmlinclude OpenMS_ConsensusIDAlgorithm.parameters
    
    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI ConsensusIDAlgorithm :
    public DefaultParamHandler
  {
  public:
    /**
        @brief Calculates the consensus ID for a set of PeptideIdentification instances of the same spectrum

        @note Make sure that the score orientation (PeptideIdentification::isHigherScoreBetter()) is set properly!
    */
    void apply(std::vector<PeptideIdentification>& ids, 
               Size number_of_runs = 0);

    /// virtual destructor
    virtual ~ConsensusIDAlgorithm();

  protected:
    /// mapping: peptide sequence -> (charge, scores)
    typedef std::map<AASequence, std::pair<Int, std::vector<double> > > 
      SequenceGrouping;

   /// Number of peptide hits considered per ID run (input parameter)
    Size considered_hits_;

    /// Number of ID runs
    Size number_of_runs_;

    /// fraction of required support by other ID runs (input parameter)
    double min_support_;

    /// count empty runs in "min_support" calculation? (input parameter)
    bool count_empty_;
   
    /// Default constructor
    ConsensusIDAlgorithm();

    /// consensus computation (to be implemented by subclasses)
    virtual void apply_(std::vector<PeptideIdentification>& ids,
                        SequenceGrouping& results) = 0;

    /// Docu in base class
    virtual void updateMembers_();

    /// compare (and possibly update) charge state information
    void compareChargeStates_(Int& recorded_charge, Int new_charge, 
                              const AASequence& peptide);

  private:
   /// Not implemented
    ConsensusIDAlgorithm(const ConsensusIDAlgorithm&);

    /// Not implemented
    ConsensusIDAlgorithm& operator=(const ConsensusIDAlgorithm&);

  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_CONSENSUSIDALGORITHM_H
