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
//

#include <OpenMS/ANALYSIS/PIP/LocalLinearMap.h>
#include <OpenMS/SYSTEM/File.h>

#include <fstream>

using namespace std;


namespace OpenMS
{

  LocalLinearMap::LocalLinearMap()
  {
    String codefile = "/PIP/codebooks.data";
    String a_file = "/PIP/linearMapping.data";
    UInt xdim = 1;
    UInt ydim = 2;
    DoubleReal radius = 0.4;

    param_.xdim = xdim;
    param_.ydim = ydim;
    param_.radius = radius;

    code_ = Matrix<DoubleReal>(param_.xdim * param_.ydim, 18);
    A_ = Matrix<DoubleReal>(param_.xdim * param_.ydim, 18);
    wout_ = vector<DoubleReal>(param_.xdim * param_.ydim);

    // path to codefile + a_file
    codefile = File::find(codefile);
    a_file = File::find(a_file);

    // read in file containing codebook vectors
    ifstream inputstream_c(codefile.c_str());
    string line;
    UInt i = 0, j = 0, k = 0;
    //open stream
    if (inputstream_c.good())
    {
      DoubleReal pos = 0.0;
      while (getline(inputstream_c, line, '\n'))
      {
        istringstream linestream(line);
        string proto;
        while (getline(linestream, proto, ' '))
        {
          stringstream(proto) >> pos;
          i = (UInt)k / 18;
          j = (UInt)k % 18;
          code_.setValue(i, j, pos);
          k++;
        }
      }
      inputstream_c.close();
      //reading codefile is done with success
    }
    else
    {
      //Throw Exception when file not found for codebooks.data
      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("LocalLinearMap could not open 'codebooks.data' at: ") + codefile);
    }

    // read in file containing Matrix< DoubleReal > A
    ifstream inputstream_a(a_file.c_str());
    //string line;
    i = 0;
    k = 0;
    //open stream
    if (inputstream_a.good())
    {
      DoubleReal pos = 0.0;
      while (getline(inputstream_a, line, '\n'))
      {
        istringstream linestream(line);
        string map;
        while (getline(linestream, map, ' '))
        {
          stringstream(map) >> pos;
          i = (UInt)k / (19);
          j = (UInt)k % (19);
          if (j > 0)       //<18)
          {
            A_.setValue((UInt)(k - 1) / (19), (UInt)(k - 1) % (19), pos);
          }
          else if (j == 0)
          {
            wout_[i] = pos;
          }
          k++;
        }
      }
      inputstream_a.close();
      //reading a_file is done with success
    }
    else
    {
      //Throw Exception when file not found for codebooks.data
      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("LocalLinearMap could not open 'linearMapping.data' at: ") + a_file);
    }
    //init 2 dimensional grid of size xdim x ydim and store in coordinates
    cord_ = genCord_(param_.xdim, param_.ydim);

  }

  LocalLinearMap::~LocalLinearMap()
  {
  }

  Matrix<UInt> LocalLinearMap::genCord_(Size xdim, Size ydim)
  {
    //store 2 dim coordinates
    Matrix<UInt> foo(xdim * ydim, 2);
    for (Size i = 0; i < xdim * ydim; i++)
    {
      foo.setValue(i, 0, UInt(i / ydim));
      foo.setValue(i, 1, UInt(i % ydim));
    }
    return foo;
  }

  vector<DoubleReal> LocalLinearMap::neigh(const Matrix<UInt> & cord, Size win, DoubleReal radius)
  {
    vector<DoubleReal> neighborhood(cord.rows());

    for (Size i = 0; i < cord.rows(); ++i)
    {
      // get dist for code i to winner code on grid structure
      DoubleReal dd = dist_(cord, cord, i, win);
      //Gaussian neighborhood function
      dd = exp(-dd / 2.0 / radius / radius);
      neighborhood[i] = dd;
    }

    return neighborhood;
  }

  DoubleReal LocalLinearMap::dist_(const Matrix<UInt> & u, const Matrix<UInt> & v, Size a, Size b)
  {
    DoubleReal dd = 0.0;
    //get euclidean distance of instances a of u and b of v
    for (Size i = 0; i < u.cols(); ++i)
    {
      dd += (u.getValue(a, i) - v.getValue(b, i)) * (u.getValue(a, i) - v.getValue(b, i));
    }
    return dd;
  }

  const Real normMeanFactors[18] = // hint: remove 'f' IFF this should ever be DoubleReal
  {
    0.5967742f, 11.5440323f, 0.4193548f, 1.2177419f, 11.9581452f,
    1399.2211022f, 0.1935484f, 412.0838710f, 0.1209677f, 1358.0966317f,
    160.5080645f, 475.8736559f, -14.4842204f, 0.4892473f, 1.6975806f,
    3.0309624f, 14.0243817f, 0.3118280f
  };


  const Real normStdFactors[18] = // hint: remove 'f' IFF this should ever be DoubleReal
  {
    0.5179165f, 5.7367444f, 0.6780753f, 0.4962471f, 5.1953755f,
    51.6311526f, 0.4527976f, 205.0635677f, 0.3727817f, 571.4667323f,
    208.2837647f, 389.9339603f, 18.0231208f, 0.7647155f, 10.0989402f,
    1.4787198f, 10.3548547f, 0.5635562f
  };

  //center and scale by variance
  void LocalLinearMap::normalizeVector(vector<DoubleReal> & aaIndexVariables)
  {
    for (Size i = 0; i < aaIndexVariables.size(); ++i)
    {
      //subtract precalculated mean of instances the model was trained on
      aaIndexVariables[i] = aaIndexVariables[i] - normMeanFactors[i];
      //and divide by sd
      aaIndexVariables[i] = aaIndexVariables[i] / normStdFactors[i];
    }
  }

  const LocalLinearMap::LLMParam & LocalLinearMap::getLLMParam() const
  {
    return param_;
  }

  const Matrix<DoubleReal> & LocalLinearMap::getCodebooks() const
  {
    return code_;
  }

  const Matrix<DoubleReal> & LocalLinearMap::getMatrixA() const
  {
    return A_;
  }

  const vector<DoubleReal> & LocalLinearMap::getVectorWout() const
  {
    return wout_;
  }

  const Matrix<UInt> & LocalLinearMap::getCord() const
  {
    return cord_;
  }

} //namespace OpenMS
