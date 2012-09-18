// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Sandro Andreotti $
// $Authors: Sandro Andreotti $
// --------------------------------------------------------------------------


#ifndef OPENMS_CHEMISTRY_SVMTHEORETICALSPECTRUMGENERATORSET_H
#define OPENMS_CHEMISTRY_SVMTHEORETICALSPECTRUMGENERATORSET_H

#include <OpenMS/SIMULATION/SimTypes.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGenerator.h>

namespace OpenMS
{
  /**
   @brief Loads SvmTheoreticalSpectrumGenerator instances for different charges

   The input file contains pairs of charge and svm models separated by a ":"
   (see share/OpenMS/examples/simulation/SvmModelSet.model)

   <p>
   Actually only a test model is shipped with OpenMS.<br>
   Please find trained models at: http://sourceforge.net/projects/open-ms/files/Supplementary/Simulation/.
   </p>

   @ingroup Chemistry
   */
  class OPENMS_DLLAPI SvmTheoreticalSpectrumGeneratorSet
  {
    public:

      /** @name Constructors and Destructors
       */
      //@{
      /// Default constructor
      SvmTheoreticalSpectrumGeneratorSet();

      /// Copy constructor
      SvmTheoreticalSpectrumGeneratorSet(const SvmTheoreticalSpectrumGeneratorSet& source);

      /// Destructor
      virtual ~SvmTheoreticalSpectrumGeneratorSet();
      //@}

      /// Assignment operator
      SvmTheoreticalSpectrumGeneratorSet& operator =(const SvmTheoreticalSpectrumGeneratorSet& tsg);

      /// Generate the MS/MS according to the model for the given precursor_charge
      void simulate(RichPeakSpectrum &spectrum, const AASequence &peptide, const gsl_rng *rng, Size precursor_charge);

      ///Load a trained Svm and Prob. models
      void load(String);

      ///Return precursor charges for which a model is contained in the set
      void getSupportedCharges(std::set<Size>&charges);

      ///return a modifiable reference to the SVM model with given charge. If charge is not supported throw exception
      SvmTheoreticalSpectrumGenerator & getSvmModel(Size);

    protected:
      //map containing the simulator for each charge variant
      std::map<Size, SvmTheoreticalSpectrumGenerator>simulators_;

  };


}



#endif
