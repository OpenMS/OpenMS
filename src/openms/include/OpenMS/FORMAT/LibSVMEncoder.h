// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>
#include <svm.h>

#include <vector>
#include <utility>

namespace OpenMS
{
  /**
    @brief Serves for encoding sequences into feature vectors

    The class can be used to construct composition vectors for
    sequences. Additionally the vectors can be encoded into
    the libsvm format.

  */
  class OPENMS_DLLAPI LibSVMEncoder
  {
public:
    /// Constructor
    LibSVMEncoder();
    /// Destructor
    ~LibSVMEncoder();

    /**
              @brief stores a composition vector of 'sequence' in 'encoded_vector'

              The allowed characters given by 'allowed_characters' are counted in the sequence 'sequence'
              and the relative frequency of the letters are stored in the composition vector.
              The first entry of the vector (<UInt, double>) corresponds to the first letter of
              'allowed_characters' that has a non zero frequency in 'sequence' and its corresponding
              relative frequency...
          */
    void encodeCompositionVector(const String & sequence, std::vector<std::pair<Int, double> > & encoded_vector, const String & allowed_characters = "ACDEFGHIKLMNPQRSTVWY");

    /**
      @brief stores composition vectors of the sequences given by 'sequence' in 'composition_vectors'

      The allowed characters given by 'allowed_characters' are counted in the sequences 'sequences'
      and the relative frequency of the letters are stored in the composition vectors.
      The first entry of the first vector (<UInt, double>) corresponds to the first letter of
      'allowed_characters' that has a non zero frequency in the first 'sequence' and its corresponding
      relative frequency...
    */
    void encodeCompositionVectors(const std::vector<String> & sequences, const String & allowed_characters, std::vector<std::vector<std::pair<Int, double> > > & composition_vectors);
    /// encodes the feature vector in LibSVM compliant format
    svm_node * encodeLibSVMVector(const std::vector<std::pair<Int, double> > & feature_vector);

    /// encodes the feature vectors in LibSVM compliant format
    void encodeLibSVMVectors(const std::vector<std::vector<std::pair<Int, double> > > & feature_vectors, std::vector<svm_node *> & libsvm_vectors);

    /// encodes the LibSVM compliant vectors into a LibSVM compliant structure
    svm_problem * encodeLibSVMProblem(const std::vector<svm_node *> & vectors,
                                      std::vector<double> & labels);

    /// creates composition vectors for 'sequences' and stores them in LibSVM compliant format
    svm_problem * encodeLibSVMProblemWithCompositionVectors(const std::vector<String> & sequences,
                                                            std::vector<double> & labels,
                                                            const String & allowed_characters);

    /**
      @brief creates composition vectors with additional length information for 'sequences' and stores them in LibSVM compliant format

      @note The resulting svm_problem needs to be destroyed by calling LibSVMEncoder::destroyProblem
    */
    svm_problem * encodeLibSVMProblemWithCompositionAndLengthVectors(const std::vector<String> & sequences,
                                                                     std::vector<double> & labels,
                                                                     const String & allowed_characters,
                                                                     UInt maximum_sequence_length);

    /**
      @brief creates composition vectors with additional length and average weight information for 'sequences' and stores them in LibSVM compliant format

      @note The resulting svm_problem needs to be destroyed by calling LibSVMEncoder::destroyProblem
    */
    svm_problem * encodeLibSVMProblemWithCompositionLengthAndWeightVectors(const std::vector<String> & sequences,
                                                                           std::vector<double> & labels,
                                                                           const String & allowed_characters);

    /// stores the LibSVM-encoded data in a text file that can be used by the LibSVM applications (svm-scale, svm-train,...)
    bool storeLibSVMProblem(const String & filename, const svm_problem * problem) const;

    /// loads the LibSVM-encoded data stored in 'filename'
    svm_problem * loadLibSVMProblem(const String & filename);

    /// encodes the borders of the sequence as k_mer oligos and stores them in 'libsvm_vector'
    void encodeOligoBorders(String                                                                                       sequence,
                            UInt                                                                             k_mer_length,
                            const String & allowed_characters,
                            UInt                                        border_length,
                            std::vector<std::pair<Int, double> > & libsvm_vector,
                            bool                                                                                             strict = false,
                            bool                                                                                             unpaired = false,
                            bool                                                                                             length_encoding = false);

    /// creates oligo border vectors vectors for 'sequences' and stores them in LibSVM compliant format
    svm_problem * encodeLibSVMProblemWithOligoBorderVectors(const std::vector<String> & sequences,
                                                            std::vector<double> & labels,
                                                            UInt                           k_mer_length,
                                                            const String & allowed_characters,
                                                            UInt                           border_length,
                                                            bool                                                   strict = false,
                                                            bool                                                   unpaired = false,
                                                            bool                                                   length_encoding = false);

    /// creates oligo border vectors vectors for 'sequences' and stores them in 'vectors'
    void encodeProblemWithOligoBorderVectors(const std::vector<AASequence> & sequences,
                                             UInt                                                       k_mer_length,
                                             const String & allowed_characters,
                                             UInt                                                       border_length,
                                             std::vector<std::vector<std::pair<Int, double> > > & vectors);

    /**
              @brief stores a string representation of the encoded sequence 'vector' in 'output'

              This function can be used if one wants to print one feature vector that is used in
              the libsvm.
          */
    void libSVMVectorToString(svm_node * vector, String & output);

    /**
              @brief stores a string representation of the encoded sequences in 'vectors' in 'output'

              This function can be used if one wants to print the feature vectors that are used in
              the libsvm.
          */
    void libSVMVectorsToString(svm_problem * vector, String & output);

    /**
              @brief encodes an AASequence instance in oligo encoding

              This function is used to get the oligo encoding for AASequence
      'sequence'. If a residue is modified, it gets an extra oligo function.
          */
    void encodeOligo(const AASequence & sequence,
                     UInt k_mer_length,
                     const String & allowed_characters,
                     std::vector<std::pair<Int, double> > & values,
                     bool is_right_border = false);

    /**
              @brief frees all the memory of the svm_problem instance

              This function is used to free all the memory used by 'problem'
          */
    static void destroyProblem(svm_problem * problem);

    static std::vector<double> predictPeptideRT(const std::vector<String> & sequences,
                                                SVMWrapper& svm,
                                                const String & allowed_characters = "ACDEFGHIKLMNPQRSTVWY",
                                                UInt maximum_sequence_length = 50)
    {
      std::vector<double> predicted_retention_times;

      LibSVMEncoder encoder;
      std::vector<double> temp_rts;
      temp_rts.resize(sequences.size(), 0);
      svm_problem * prediction_data =
          encoder.encodeLibSVMProblemWithCompositionAndLengthVectors(sequences,
                                                                     temp_rts,
                                                                     allowed_characters,
                                                                     maximum_sequence_length);
      svm.predict(prediction_data, predicted_retention_times);
      LibSVMEncoder::destroyProblem(prediction_data);
      return predicted_retention_times;
    }

private:
    /// comparator for oligos encoded by encodeOligo
    static bool cmpOligos_(std::pair<Int, double> a,
                           std::pair<Int, double> b);

  };

} // namespace OpenMS

