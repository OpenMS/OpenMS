// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

#include <iostream>
#include <fstream>

using namespace std;

namespace OpenMS
{
  LibSVMEncoder::LibSVMEncoder()
  {

  }

  LibSVMEncoder::~LibSVMEncoder()
  {

  }

  void LibSVMEncoder::encodeCompositionVector(const String& sequence,
                                              vector<pair<Int, double> >& composition_vector,
                                              const String& allowed_characters)
  {

    Size number_of_different_letters = allowed_characters.size();
    Size* counts = new Size[number_of_different_letters];
    Size total_count = 0;

    composition_vector.clear();

    for (Size i = 0; i < number_of_different_letters; i++)
    {
      counts[i] = 0;
    }

    for (Size i = 0; i < sequence.size(); i++)
    {
      if (allowed_characters.find(sequence[i]) != String::npos)
      {
        counts[allowed_characters.find(sequence[i])]++;
        total_count++;
      }
    }

    for (Size i = 0; i < number_of_different_letters; i++)
    {
      if (counts[i] > 0)
      {
        composition_vector.push_back(make_pair(Int(i + 1), (((double) counts[i]) / total_count)));
      }
    }
    delete[] counts;

  }

  void LibSVMEncoder::encodeCompositionVectors(const vector<String>& sequences,
                                               const String& allowed_characters,
                                               vector<vector<pair<Int, double> > >& composition_vectors)
  {
    vector<pair<Int, double> > composition_vector;

    composition_vectors.clear();

    for (Size i = 0; i < sequences.size(); i++)
    {
      encodeCompositionVector(sequences[i],
                              composition_vector,
                              allowed_characters);
      composition_vectors.push_back(composition_vector);
    }

  }

  svm_node* LibSVMEncoder::encodeLibSVMVector(const vector<pair<Int, double> >& feature_vector)
  {

    vector<pair<Int, double> >::const_iterator vector_iterator;
    svm_node* nodes;
    UInt i = 0;

    nodes = new svm_node[feature_vector.size() + 1];

    for (vector_iterator = feature_vector.begin();
         vector_iterator != feature_vector.end();
         ++vector_iterator)
    {
      nodes[i].index = vector_iterator->first;
      nodes[i].value = vector_iterator->second;
      i++;
    }
    nodes[feature_vector.size()].index = -1;
    nodes[feature_vector.size()].value = 0;

    return nodes;
  }

  void LibSVMEncoder::encodeLibSVMVectors(const vector<vector<pair<Int, double> > >& feature_vectors,
                                          vector<svm_node*>& libsvm_vectors)
  {
    libsvm_vectors.clear();

    for (Size i = 0; i < feature_vectors.size(); i++)
    {
      libsvm_vectors.push_back(encodeLibSVMVector(feature_vectors[i]));
    }
  }

  svm_problem* LibSVMEncoder::encodeLibSVMProblem(const vector<svm_node*>& vectors,
                                                  vector<double>& labels)
  {
    svm_problem* problem;
    svm_node** node_vectors;

    if (labels.size() != vectors.size())
    {
      return nullptr;
    }

    problem = new svm_problem;
    problem->l = (int) vectors.size();
    if (problem->l < 0) // dubious. Just makes sense if vectors.size() is larger than int and overflows
    {
      return nullptr;
    }

    problem->y = new double[problem->l];
    for (Size i = 0; i < vectors.size(); ++i)
    {
      problem->y[i] = labels[i];
    }

    node_vectors = new svm_node*[problem->l];

    for (Size i = 0; i < vectors.size(); i++)
    {
      node_vectors[i] = vectors[i];
    }
    problem->x = node_vectors;

// cout << "Problem encoded" << "\n";

    return problem;
  }

  svm_problem* LibSVMEncoder::encodeLibSVMProblemWithCompositionVectors(const vector<String>& sequences,
                                                                        std::vector<double>& labels,
                                                                        const String& allowed_characters)
  {
    vector<svm_node*> vectors;
    vector<pair<Int, double> > encoded_vector;

    for (Size i = 0; i < sequences.size(); i++)
    {
      encodeCompositionVector(sequences[i], encoded_vector, allowed_characters);
      svm_node* libsvm_vector = encodeLibSVMVector(encoded_vector);
      vectors.push_back(libsvm_vector);
    }

    return encodeLibSVMProblem(vectors, labels);
  }

  svm_problem* LibSVMEncoder::encodeLibSVMProblemWithCompositionAndLengthVectors(const vector<String>& sequences,
                                                                                 std::vector<double>& labels,
                                                                                 const String& allowed_characters,
                                                                                 UInt                           maximum_sequence_length)
  {
    vector<svm_node*> vectors;
    vector<pair<Int, double> > encoded_vector;

    for (Size i = 0; i < sequences.size(); i++)
    {

      encodeCompositionVector(sequences[i], encoded_vector, allowed_characters);
      encoded_vector.emplace_back(Int(allowed_characters.size() + 1), ((double) sequences[i].length()) / maximum_sequence_length);
      svm_node* libsvm_vector = encodeLibSVMVector(encoded_vector);
      vectors.push_back(libsvm_vector);
    }

    return encodeLibSVMProblem(vectors, labels);
  }

  svm_problem* LibSVMEncoder::encodeLibSVMProblemWithCompositionLengthAndWeightVectors(const vector<String>& sequences,
                                                                                       std::vector<double>& labels,
                                                                                       const String& allowed_characters)
  {
    vector<svm_node*> vectors;
    vector<pair<Int, double> > encoded_vector;

    for (Size i = 0; i < sequences.size(); i++)
    {

      encodeCompositionVector(sequences[i], encoded_vector, allowed_characters);
      encoded_vector.push_back(make_pair(Int(allowed_characters.size() + 1), (double) sequences[i].length()));
      encoded_vector.push_back(make_pair(Int(allowed_characters.size() + 2), AASequence::fromString(sequences[i]).getAverageWeight()));
      svm_node* libsvm_vector = encodeLibSVMVector(encoded_vector);
      vectors.push_back(libsvm_vector);
    }

    return encodeLibSVMProblem(vectors, labels);
  }

  bool LibSVMEncoder::storeLibSVMProblem(const String& filename, const svm_problem* problem) const
  {
    if (problem == nullptr)
    {
      return false;
    }

    ofstream output_file(filename.c_str());
    Int j = 0;

    // checking if file is writable
    if (!File::writable(filename))
    {
      return false;
    }

    // writing feature vectors
    for (Int i = 0; i < problem->l; i++)
    {
      j = 0;
      output_file << problem->y[i] << " ";
      while (problem->x[i][j].index != -1)
      {
        output_file << problem->x[i][j].index << ":" << problem->x[i][j].value << " ";
        ++j;
      }
      output_file << "\n";
    }
    output_file.flush();
    output_file.close();
    cout.flush();
    return true;
  }

  svm_problem* LibSVMEncoder::loadLibSVMProblem(const String& filename)
  {
    svm_problem* data = nullptr;
    UInt counter = 0;
    vector<String> parts;
    vector<String> temp_parts;

    if (!File::exists(filename))
    {
      return nullptr;
    }
    if (!File::readable(filename))
    {
      return nullptr;
    }
    if (File::empty(filename))
    {
      return nullptr;
    }

    TextFile text_file(filename.c_str(), true);
    TextFile::ConstIterator it;

    it = text_file.begin();

    data = new svm_problem;
    data->l = (Int)(text_file.end() - text_file.begin());
    data->x = new svm_node*[(text_file.end() - text_file.begin())];
    data->y = new double[(text_file.end() - text_file.begin())];
    while (it != text_file.end())
    {
      it->split(' ', parts);
      data->y[counter] = parts[0].trim().toFloat();
      data->x[counter] = new svm_node[parts.size()];
      for (Size j = 1; j < parts.size(); ++j)
      {
        parts[j].split(':', temp_parts);
        if (temp_parts.size() < 2)
        {
          delete data;
          return nullptr;
        }
        data->x[counter][j - 1].index = temp_parts[0].trim().toInt();
        data->x[counter][j - 1].value = temp_parts[1].trim().toFloat();
      }
      data->x[counter][parts.size() - 1].index = -1;
      data->x[counter][parts.size() - 1].value = 0;
      ++counter;
      ++it;
    }
    return data;
  }

  void LibSVMEncoder::encodeOligoBorders(String sequence,
                                         UInt k_mer_length,
                                         const String& allowed_characters,
                                         UInt border_length,
                                         vector<pair<Int, double> >& libsvm_vector,
                                         bool strict,
                                         bool unpaired,
                                         bool length_encoding)
  {
    multimap<Int, Int>                        ordered_tree;
    multimap<Int, Int>::const_iterator    elements;
    multimap<Int, Int>::const_iterator    elements_start;
    multimap<Int, Int>::const_iterator    elements_end;
    pair<Int, Int>                            values;
    Size                                                  oligo_value = 0;
    Size                                              factor      = 1;
    map<String::value_type, Size>             residue_values;
    Size                                                              counter         = 0;
    Size                                                              number_of_residues = 0;
    Size                                                              left_border = 0;
    Size                                                              right_border = 0;
    Size                                                              sequence_length = 0;
    bool                                                              wrong_characters = false;

    number_of_residues = allowed_characters.size();

    if (number_of_residues == 0)
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, number_of_residues);
    }

    libsvm_vector.clear();
    sequence_length = sequence.size();

    for (Size i = 0; i < sequence.size(); ++i)
    {
      if (!allowed_characters.has(sequence.at(i)))
      {
        wrong_characters = true;
      }
    }
    if (k_mer_length <= sequence_length && !wrong_characters)
    {
      // if a border must not be longer than half of the peptide
      if (strict)
      {
        if (border_length > (sequence_length - k_mer_length + 1) / 2)
        {
          left_border = (UInt) (floor((sequence_length - k_mer_length + 1.0) / 2));
          right_border = (UInt) (ceil((sequence_length - k_mer_length + 1.0) / 2));
        }
        else
        {
          left_border = border_length;
          right_border = sequence_length - k_mer_length + 1 - border_length;
        }
      }
      else
      {
        if (border_length >= sequence_length - k_mer_length + 1)
        {
          left_border = sequence_length - k_mer_length + 1;
          right_border = 0;
        }
        else
        {
          left_border = border_length;
          right_border = sequence_length - k_mer_length + 1 - border_length;
        }
      }

      for (Size i = 0; i < number_of_residues; ++i)
      {
        residue_values.insert(make_pair(allowed_characters[i], counter));
        ++counter;
      }
      for (Int k = k_mer_length - 1; k >= 0; k--)
      {
        oligo_value += factor * residue_values[sequence[k]];
        factor *= number_of_residues;
      }
      factor /= number_of_residues;
      values.first = ((Int) (oligo_value + 2));
      values.second = 1;
      ordered_tree.insert(values);

      for (Size j = 1; j < left_border; j++)
      {
        oligo_value -= factor * residue_values[sequence[j - 1]];
        oligo_value = oligo_value * number_of_residues + residue_values[sequence[j + k_mer_length - 1]];

        values.first = ((Int) (oligo_value + 2));
        values.second = (Int)j + 1;

        ordered_tree.insert(values);
      }
      oligo_value = 0;
      factor = 1;

      if (k_mer_length > 1)
      {
        for (Int k = k_mer_length; k > 0; k--)
        {
          oligo_value += factor * residue_values[sequence[sequence_length - k]];
          factor *= number_of_residues;
        }
        factor /= number_of_residues;
        if (unpaired)
        {
          values.first = ((Int) (oligo_value + 2)) * -1;
        }
        else
        {
          values.first = ((Int) (oligo_value + 2));
        }
        values.second = 1;
        ordered_tree.insert(values);

        for (Size j = 1; j < left_border; j++)
        {
          oligo_value -= factor * residue_values[sequence[sequence_length - j]];
          oligo_value = oligo_value * number_of_residues + residue_values[sequence[sequence_length - k_mer_length - j]];

          if (unpaired)
          {
            values.first = ((Int) (oligo_value + 2)) * -1;
          }
          else
          {
            values.first = ((Int) (oligo_value + 2));
          }
          values.second = (Int)j + 1;

          ordered_tree.insert(values);
        }
      }
      else
      {
        for (Size k = right_border + k_mer_length; k > right_border; k--)
        {
          oligo_value += factor * residue_values[sequence[k - 1]];
          factor *= number_of_residues;
        }
        factor /= number_of_residues;
        if (unpaired)
        {
          values.first = ((Int) (oligo_value + 2)) * -1;
        }
        else
        {
          values.first = ((Int) (oligo_value + 2));
        }
        values.second = (Int)(right_border - sequence_length) * -1;
        ordered_tree.insert(values);
        for (Size j = right_border + 1; j < sequence_length - k_mer_length + 1; j++)
        {
          oligo_value -= factor * residue_values[sequence[j - 1]];
          oligo_value = oligo_value * number_of_residues + residue_values[sequence[j + k_mer_length - 1]];

          if (unpaired)
          {
            values.first = ((Int) (oligo_value + 2)) * -1;
          }
          else
          {
            values.first = ((Int) (oligo_value + 2));
          }
          values.second = (Int)(j - sequence_length) * -1;

          ordered_tree.insert(values);
        }
      }
      pair<multimap<Int, Int>::const_iterator, multimap<Int, Int>::const_iterator> range_iterators;
      vector<Int> temp_positions;
      elements = ordered_tree.begin();
      while (elements != ordered_tree.end())
      {
        temp_positions.clear();
        range_iterators = ordered_tree.equal_range(elements->first);
        elements_start = range_iterators.first;
        elements_end = range_iterators.second;
        while (elements_start != elements_end)
        {
          temp_positions.push_back(elements_start->second);
          ++elements_start;
        }
        sort(temp_positions.begin(), temp_positions.end());
        for (Size i = 0; i < temp_positions.size(); ++i)
        {
          libsvm_vector.push_back(make_pair(elements->first, temp_positions[i]));
        }
        elements = elements_end;
      }
      if (length_encoding)
      {
        libsvm_vector.push_back(make_pair((Int) sequence.size(), (double)pow((double)k_mer_length, (double) number_of_residues) + 1));
      }
    }

  }

  svm_problem* LibSVMEncoder::encodeLibSVMProblemWithOligoBorderVectors(const vector<String>& sequences,
                                                                        vector<double>& labels,
                                                                        UInt                                            k_mer_length,
                                                                        const String& allowed_characters,
                                                                        UInt                                            border_length,
                                                                        bool                                            strict,
                                                                        bool                                            unpaired,
                                                                        bool                                            length_encoding)
  {
    vector<svm_node*> vectors;
    vector<pair<Int, double> > encoded_vector;

    for (Size i = 0; i < sequences.size(); i++)
    {
      encodeOligoBorders(sequences[i], k_mer_length, allowed_characters, border_length, encoded_vector, strict, unpaired, length_encoding);
      svm_node* libsvm_vector = encodeLibSVMVector(encoded_vector);
      vectors.push_back(libsvm_vector);
    }

    return encodeLibSVMProblem(vectors, labels);
  }

  void LibSVMEncoder::encodeProblemWithOligoBorderVectors(const vector<AASequence>& sequences,
                                                          UInt                                                                                k_mer_length,
                                                          const String& allowed_characters,
                                                          UInt                                                                                border_length,
                                                          vector<vector<pair<Int, double> > >& vectors)
  {
    vector<pair<Int, double> > temp_encoded_vector_left;
    vector<pair<Int, double> > temp_encoded_vector_right;

    vectors.clear();
    for (Size i = 0; i < sequences.size(); i++)
    {
      if (sequences[i].size() > border_length)
      {
        encodeOligo(sequences[i].getPrefix(border_length), k_mer_length, allowed_characters, temp_encoded_vector_left, false);
        encodeOligo(sequences[i].getSuffix(border_length), k_mer_length, allowed_characters, temp_encoded_vector_right, true);
      }
      else
      {
        encodeOligo(sequences[i], k_mer_length, allowed_characters, temp_encoded_vector_left, false);
        encodeOligo(sequences[i], k_mer_length, allowed_characters, temp_encoded_vector_right, true);
      }
      temp_encoded_vector_left.insert(temp_encoded_vector_left.end(), temp_encoded_vector_right.begin(), temp_encoded_vector_right.end());
      stable_sort(temp_encoded_vector_left.begin(), temp_encoded_vector_left.end(), cmpOligos_);
      vectors.push_back(temp_encoded_vector_left);
    }
  }

  void LibSVMEncoder::libSVMVectorToString(svm_node* vector, String& output)
  {
    UInt i = 0;

    output.clear();
    while (vector[i].index != -1)
    {
      output = output + "(" + String(vector[i].index) + ", " + String(vector[i].value) + ") ";
      ++i;
    }
  }

  void LibSVMEncoder::libSVMVectorsToString(svm_problem* vector, String& output)
  {
    String temp_string = "";

    output.clear();
    if (vector != nullptr)
    {
      for (Int i = 0; i < vector->l; ++i)
      {
        libSVMVectorToString(vector->x[i], temp_string);
        output = output + temp_string + "\n";
        temp_string = "";
      }
    }
  }

  void LibSVMEncoder::encodeOligo(const AASequence& sequence,
                                  UInt k_mer_length,
                                  const String& allowed_characters,
                                  vector<pair<Int, double> >& values,
                                  bool is_right_border)
  {
    map<String, UInt> residue_values;
    UInt counter = 0;
    Size number_of_residues = allowed_characters.size();
    Size sequence_length = sequence.size();
    bool sequence_ok = true;
    ModificationsDB* modifications = ModificationsDB::getInstance();
    Size number_of_modifications = modifications->getNumberOfModifications();

    // checking if sequence contains illegal characters
    for (Size i = 0; i < sequence.size(); ++i)
    {
      if (allowed_characters.find(sequence[i].getOneLetterCode()) == String::npos)
      {
        sequence_ok = false;
      }
    }

    if (sequence_ok && k_mer_length <= sequence_length)
    {
      double oligo_value = 0.;
      double factor_simple = double(number_of_residues * (number_of_modifications + 1));
      double factor = 1.;

      values.resize(sequence_length - k_mer_length + 1, pair<Int, double>());
      for (Size i = 0; i < number_of_residues; ++i)
      {
        residue_values.insert(make_pair(String(allowed_characters[i]), counter));
        ++counter;
      }

      if (!is_right_border || k_mer_length == 1)
      {

        for (SignedSize k = (Int) k_mer_length - 1; k >= 0; k--)
        {
          if (sequence[k].isModified())
          {
            oligo_value += factor * (residue_values[(sequence.getResidue(k)).getOneLetterCode()]
                                     + (modifications->findModificationIndex(sequence.getResidue(k).getModification()->getFullId()) + 1) * number_of_residues);
          }
          else
          {
            oligo_value += factor * residue_values[sequence[k].getOneLetterCode()];
          }
          factor *= factor_simple;
        }
        factor /= factor_simple;
        counter = 0;
        if (is_right_border)
        {
          values[counter].first = Int(sequence_length - k_mer_length + 1);
        }
        else
        {
          values[counter].first = 1;
        }
        values[counter].second = oligo_value;
        ++counter;

        for (Size j = 1; j < sequence_length - k_mer_length + 1; j++)
        {
          if (sequence[j - 1].isModified())
          {
            oligo_value -= factor * (residue_values[(sequence.getResidue(j - 1)).getOneLetterCode()]
                                     + (modifications->findModificationIndex(sequence.getResidue(j - 1).getModification()->getFullId()) + 1) * number_of_residues);
          }
          else
          {
            oligo_value -= factor * residue_values[sequence[j - 1].getOneLetterCode()];
          }
          if (sequence[j + k_mer_length - 1].isModified())
          {
            oligo_value = oligo_value * factor_simple + (residue_values[sequence[j + k_mer_length - 1].getOneLetterCode()]
                                                         + (modifications->findModificationIndex(sequence[j + k_mer_length - 1].getModification()->getFullId()) + 1) * number_of_residues);
          }
          else
          {
            oligo_value = oligo_value * factor_simple + residue_values[sequence[j + k_mer_length - 1].getOneLetterCode()];
          }
          if (is_right_border)
          {
            values[counter].first = Int(sequence_length - k_mer_length - j + 1);
          }
          else
          {
            values[counter].first = (Int)j + 1;
          }
          values[counter].second = oligo_value;
          ++counter;
        }
        stable_sort(values.begin(), values.end(), cmpOligos_);
      }
      else
      {
        for (SignedSize k = sequence_length - k_mer_length; k < (Int) sequence_length; k++)
        {
          if (sequence[k].isModified())
          {
            oligo_value += factor * (residue_values[(sequence.getResidue(k)).getOneLetterCode()]
                                     + (modifications->findModificationIndex(sequence.getResidue(k).getModification()->getFullId()) + 1) * number_of_residues);
          }
          else
          {

            oligo_value += factor * residue_values[sequence[k].getOneLetterCode()];
          }
          factor *= factor_simple;
        }
        factor /= factor_simple;
        counter = 0;
        values[counter].first = 1;
        values[counter].second = oligo_value;
        ++counter;

        for (SignedSize j = sequence_length - k_mer_length - 1; j >= 0; j--)
        {
          if (sequence[(Size) j + k_mer_length].isModified())
          {
            oligo_value -= factor * (residue_values[(sequence.getResidue((Size) j + k_mer_length)).getOneLetterCode()]
                                     + (modifications->findModificationIndex(sequence.getResidue((Size) j + k_mer_length).getModification()->getFullId()) + 1) * number_of_residues);
          }
          else
          {
            oligo_value -= factor * residue_values[sequence[(Size) j + k_mer_length].getOneLetterCode()];
          }
          if (sequence[j].isModified())
          {
            oligo_value = oligo_value * factor_simple + (residue_values[sequence[j].getOneLetterCode()]
                                                         + (modifications->findModificationIndex(sequence[j].getModification()->getFullId()) + 1) * number_of_residues);
          }
          else
          {
            oligo_value = oligo_value * factor_simple + residue_values[sequence[j].getOneLetterCode()];
          }

          values[counter].first = (Int)(sequence_length - k_mer_length - j + 1);

          values[counter].second = oligo_value;
          ++counter;
        }
        stable_sort(values.begin(), values.end(), cmpOligos_);

      }
    }
    else
    {
      values.clear();
    }
  }

  bool LibSVMEncoder::cmpOligos_(pair<Int, double> a, pair<Int, double> b)
  {
    return (a.second == b.second) ? (a.first < b.first) : (a.second < b.second);
  }

  void LibSVMEncoder::destroyProblem(svm_problem* problem)
  {
    if (problem != nullptr)
    {
      for (Int  i = 0; i < problem->l; i++)
      {
        delete[] problem->x[i];
      }
      delete[] problem->y;
      delete[] problem->x;
      delete problem;
    }
    problem = nullptr;
  }

} // namespace OpenMS
