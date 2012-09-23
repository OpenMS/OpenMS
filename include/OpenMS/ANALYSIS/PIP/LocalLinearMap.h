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
// $Maintainer: Alexandra Scherbart $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_PIP_LOCALLINEARMAP_H
#define OPENMS_ANALYSIS_PIP_LOCALLINEARMAP_H

#include <OpenMS/DATASTRUCTURES/Matrix.h>

namespace OpenMS
{

  /**
      @brief Trained Local Linear %Map (LLM) model for peak intensity prediction

      This class offers a model for predictions of peptide peak heights
      (referred to as intensities) by a Local Linear %Map (LLM) model and
      is the basis of PeakIntensityPredictor.

      A general introduction to the Peak Intensity Predictor (PIP)
      can be found in the <A HREF="tutorial_pip.html">PIP Tutorial</A>.

      The model trained needs two files for storing the position of the
      codebook vectors and the linear mappings (codebooks.data, linearMapping.data)
      This is the default model used by PeakIntensityPredictor.
*/
  class OPENMS_DLLAPI LocalLinearMap
  {

public:

    /**
    @brief Define parameters needed by the Local Linear %Map (LLM) model

    Parameters xdim and ydim define the size of the two dimensional
    grid structure. Parameter radius gives the width of the Gaussian
    neighborhood function.
*/
    struct OPENMS_DLLAPI LLMParam
    {
      UInt xdim;           /**< size of first coordinate */
      UInt ydim;           /**< size of second coordinate */
      DoubleReal radius;           /**< width of Gaussian neighborhood function */
    };

    /// default constructor
    LocalLinearMap();
    /// destructor
    virtual ~LocalLinearMap();

    ///return parameters of the LocalLinearMap model
    const LLMParam & getLLMParam() const;
    ///return position of the codebook vectors (18-dim)
    const Matrix<DoubleReal> & getCodebooks() const;
    ///return linear mappings of the codebooks
    const Matrix<DoubleReal> & getMatrixA() const;
    ///return linear bias
    const std::vector<DoubleReal> & getVectorWout() const;
    ///return coordinates of codebook vectors on the 2-d grid
    const Matrix<UInt> & getCord() const;
    ///calculate and return normalized amino acid index variables from string representation of peptide
    void normalizeVector(std::vector<DoubleReal> & aaIndexVariables);
    ///calculate neighborhood function based on distance of prototypes to winner prototype on two-dimensional grid structure and neighborhood width.
    std::vector<DoubleReal> neigh(const Matrix<UInt> & cord, Size win, DoubleReal radius);

private:

    LLMParam param_;                                            ///<parameters of the model
    Matrix<DoubleReal> code_;                           ///<codebook vectors
    Matrix<DoubleReal> A_;                              ///<linear mappings
    std::vector<DoubleReal> wout_;              ///<linear bias
    Matrix<UInt> cord_;                                     ///<coordinates of codebooks on grid

    /// needed to store prototype coordinates
    Matrix<UInt> genCord_(Size xdim, Size ydim);
    ///calculate distance between two prototypes
    DoubleReal dist_(const Matrix<UInt> & u, const Matrix<UInt> & v, Size a, Size b);

    ///Copy constructor not implemented => private
    LocalLinearMap(LocalLinearMap & rhs);
    ///Assignment operator not implemented => private
    LocalLinearMap & operator=(const LocalLinearMap & llm);

  };


}
#endif
