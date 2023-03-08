// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <iterator>
#include <utility>

// helper class that can be used to create iterators that iterate over single members of a struct
template <typename InputIterator, typename UnaryFunction>
class TransformIterator {
 public:
  using iterator_category = typename std::iterator_traits<InputIterator>::iterator_category;
  using value_type = typename std::result_of<UnaryFunction(typename std::iterator_traits<InputIterator>::value_type)>::type;
  using difference_type = typename std::iterator_traits<InputIterator>::difference_type;
  using pointer = value_type*;
  using reference = value_type&;

  TransformIterator(InputIterator it, UnaryFunction f)
      : it_(it), f_(f) {}

  reference operator*() const {
    return f_(*it_);
  }

  TransformIterator& operator++() {
    ++it_;
    return *this;
  }

  TransformIterator operator++(int) {
    TransformIterator tmp(*this);
    ++(*this);
    return tmp;
  }

  bool operator==(const TransformIterator& other) const {
    return it_ == other.it_;
  }

  bool operator!=(const TransformIterator& other) const {
    return !(*this == other);
  }

 private:
  InputIterator it_;
  UnaryFunction f_;
};

template <typename InputIterator, typename UnaryFunction>
TransformIterator<InputIterator, UnaryFunction> make_transform_iterator(InputIterator it, UnaryFunction f) {
  return TransformIterator<InputIterator, UnaryFunction>(it, f);
}
