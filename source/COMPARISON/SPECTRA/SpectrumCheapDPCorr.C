// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/COMPARISON/SPECTRA/SpectrumCheapDPCorr.h>

#include <cmath>

#define SPECTRUMCHEAPDPCORR_DEBUG
#undef  SPECTRUMCHEAPDPCORR_DEBUG

#ifdef SPECTRUMCHEAPDPCORR_DEBUG
#include <iostream>
#endif

#include <boost/math/distributions.hpp>

using namespace std;

namespace OpenMS
{
  SpectrumCheapDPCorr::SpectrumCheapDPCorr() :
    PeakSpectrumCompareFunctor(),
    lastconsensus_()
  {
    setName(SpectrumCheapDPCorr::getProductName());
    defaults_.setValue("variation", 0.001, "Maximum difference in position (in percent of the current m/z).\nNote that big values of variation ( 1 being the maximum ) result in consideration of all possible pairings which has a running time of O(n*n)");
    defaults_.setValue("int_cnt", 0, "How the peak heights are used in the score.\n0 = product\n1 = sqrt(product)\n2 = sum\n3 = agreeing intensity\n");
    defaults_.setValue("keeppeaks", 0, "Flag that states if peaks without alignment partner are kept in the consensus spectrum.");
    factor_ = 0.5;
    defaultsToParam_();
  }

  SpectrumCheapDPCorr::SpectrumCheapDPCorr(const SpectrumCheapDPCorr & source) :
    PeakSpectrumCompareFunctor(source),
    lastconsensus_(source.lastconsensus_),
    factor_(source.factor_)
  {
  }

  SpectrumCheapDPCorr::~SpectrumCheapDPCorr()
  {
  }

  SpectrumCheapDPCorr & SpectrumCheapDPCorr::operator=(const SpectrumCheapDPCorr & source)
  {
    if (this != &source)
    {
      PeakSpectrumCompareFunctor::operator=(source);
      lastconsensus_ = source.lastconsensus_;
      factor_ = source.factor_;
    }
    return *this;
  }

  void SpectrumCheapDPCorr::setFactor(double f)
  {
    if (f < 1 && f > 0)
    {
      factor_ = f;
    }
    else
    {
      throw Exception::OutOfRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    }
  }

  double SpectrumCheapDPCorr::operator()(const PeakSpectrum & csa) const
  {
    return operator()(csa, csa);
  }

  /**
        looks for peak pairs where there is just one or none possibility for alignment
        and aligns them (if possible). The rest is aligned using dynprog_
  */
  double SpectrumCheapDPCorr::operator()(const PeakSpectrum & x, const PeakSpectrum & y) const
  {
    double var = (double)param_.getValue("variation");
    double score(0);
    bool keeppeaks_ = (int)param_.getValue("keeppeaks");

    lastconsensus_ = PeakSpectrum();
    Precursor p1, p2;
    if (!x.getPrecursors().empty())
      p1 = x.getPrecursors()[0];
    if (!y.getPrecursors().empty())
      p2 = y.getPrecursors()[0];
    lastconsensus_.getPrecursors().resize(1);
    lastconsensus_.getPrecursors()[0].setMZ((p1.getMZ() + p2.getMZ()) / 2);
    lastconsensus_.getPrecursors()[0].setCharge(p1.getCharge());
    peak_map_.clear();

    int xpos = 0;
    int ypos = 0;
    PeakSpectrum::ConstIterator xit = x.begin();
    PeakSpectrum::ConstIterator yit = y.begin();
    while (xit != x.end() && yit != y.end())
    {
      double variation = (xit->getMZ() + yit->getMZ()) / 2 * var;

      //ignore pairs that cannot be paired
      if (fabs(xit->getMZ() - yit->getMZ()) > variation)
      {
        if (xit->getMZ() < yit->getMZ()) // while effizienter?
        {
          Peak1D consensuspeak;
          consensuspeak.setMZ(xit->getMZ());
          consensuspeak.setIntensity((xit->getIntensity()) * (1 - factor_));
          if (keeppeaks_)
            lastconsensus_.push_back(consensuspeak);
          ++xit;
          ++xpos;
        }
        else
        {
          Peak1D consensuspeak;
          consensuspeak.setMZ(yit->getMZ());
          consensuspeak.setIntensity((yit->getIntensity()) * (factor_));
          if (keeppeaks_)
            lastconsensus_.push_back(consensuspeak);
          ++yit;
          ++ypos;
        }
      }
      else
      {
        //x/yrun represents the number of peaks in both spectra that could be paired
        int xrun = 1;
        int yrun = 1;
        while (xit + xrun != x.end() && yit + yrun != y.end() &&
               (!((xit + xrun - 1)->getMZ() + variation < (yit + yrun)->getMZ()) ||
                !((yit + yrun - 1)->getMZ() + variation < (xit + xrun)->getMZ())))
        {
          if ((yit + yrun - 1)->getMZ() + variation > (xit + xrun)->getMZ())
          {
            xrun++;
          }
          else if ((xit + xrun - 1)->getMZ() + variation > (yit + yrun)->getMZ())
          {
            yrun++;
          }
          else
          {
            xrun++;
            yrun++;
          }
          if (xit + xrun == x.end())
            break;
          if (yit + yrun == y.end())
            break;
        }

        //dynamic programming necessary to calculate optimal pairing
        if (xrun > 1 && yrun > 1)
        {
          score += dynprog_(x, y, xpos, xpos + xrun - 1, ypos, ypos + yrun - 1);
          xit = xit + xrun;
          yit = yit + yrun;
          xpos += xrun;
          ypos += yrun;
        }
        //the optimal pairing of 2 peaks is easy...
        else
        {
          // calculate consensus peak
          Peak1D consensuspeak;
          consensuspeak.setMZ((xit->getMZ() * (1 - factor_) + yit->getMZ() * (factor_)));
          consensuspeak.setIntensity((xit->getIntensity() * (1 - factor_) + yit->getIntensity() * factor_));
          lastconsensus_.push_back(consensuspeak);

          if (!peak_map_.has(xit - x.begin()))
          {
            peak_map_[xit - x.begin()] = yit - y.begin();
          }
          else
          {
            peak_map_[xit - x.begin()] = yit - y.begin() > xit - x.begin() ? xit - x.begin() : yit - y.begin();
          }

          variation = (xit->getMZ() + yit->getMZ()) / 2 * var;
          score += comparepeaks_(xit->getMZ(), yit->getMZ(), xit->getIntensity(), yit->getIntensity());
          ++xit;
          ++yit;
          ++xpos;
          ++ypos;
        }
      }
    }
    factor_ = 0.5;
    return score;
  }

  double SpectrumCheapDPCorr::dynprog_(const PeakSpectrum & x, const PeakSpectrum & y, int xstart, int xend, int ystart, int yend) const
  {
#ifdef SPECTRUMCHEAPDPCORR_DEBUG
    cerr << "SpectrumCheapDPCorr::dynprog_(const DDiscreteSpectrum<1>& x, const DDiscreteSpectrum<1>& y, " << xstart << ", " << xend << ", " <<  ystart << ", " << yend << ")" <<  endl;
#endif
    double var = (double)param_.getValue("variation");
    vector<vector<double> > dparray(xend - xstart + 2, vector<double>(yend - ystart + 2));
    vector<vector<int> > trace(xend - xstart + 2, vector<int>(yend - ystart + 2));
    double align;
    for (int i = 1; i < xend - xstart + 2; ++i)
    {
      for (int j = 1; j < yend - ystart + 2; ++j)
      {
        double variation = (y[ystart + j - 1].getMZ() + x[xstart + i - 1].getMZ()) / 2 * var;
        //positions too different
        if (fabs(x[xstart + i - 1].getMZ() - y[ystart + j - 1].getMZ()) > variation)
          align = 0;

        //calculate score of alignment
        else
          align = comparepeaks_(x[xstart + i - 1].getMZ(),
                                y[ystart + j - 1].getMZ(),
                                x[xstart + i - 1].getIntensity(),
                                y[ystart + j - 1].getIntensity());
        //dynaminc programming step
        if ((((dparray[i][j - 1]) > (dparray[i - 1][j - 1] + align)) ? (dparray[i][j - 1]) : (dparray[i - 1][j - 1] + align)) /*== max*/ > dparray[i - 1][j])
        {
          if (dparray[i - 1][j - 1] + align > dparray[i][j - 1])
          {
            dparray[i][j] = dparray[i - 1][j - 1] + align;
            trace[i][j] = 5;
          }
          else
          {
            dparray[i][j] = dparray[i][j - 1];
            trace[i][j] = -1;
          }
        }
        else
        {
          dparray[i][j] = dparray[i - 1][j];
          trace[i][j] = 1;
        }
      }
    }

    unsigned int i = xend - xstart + 1;
    unsigned int j = yend - ystart + 1;
    for (;; )
    {
      if (trace[i][j] == 5)
      {
        Peak1D consensuspeak;
        consensuspeak.setMZ((y[ystart + j - 1].getMZ() * (1 - factor_) + x[xstart + i - 1].getMZ() * factor_));
        consensuspeak.setIntensity((y[ystart + j - 1].getIntensity() * (1 - factor_) + x[xstart + i - 1].getIntensity() * factor_));
        lastconsensus_.push_back(consensuspeak);
        if (!peak_map_.has(xstart + i - 1))
        {
          peak_map_[xstart + i - 1] = ystart + j - 1;
        }
        else
        {
          peak_map_[xstart + i - 1] = ystart + j - 1 > peak_map_[xstart + i - 1] ? peak_map_[xstart + i - 1] : ystart + j - 1;
        }
        i--;
        j--;
      }
      else if (trace[i][j] == 1)
      {
        Peak1D consensuspeak;
        consensuspeak.setMZ(x[xstart + i - 1].getMZ());
        consensuspeak.setIntensity((x[xstart + i - 1].getIntensity()) * (1 - factor_));
        if (keeppeaks_)
          lastconsensus_.push_back(consensuspeak);
        i--;
      }
      else if (trace[i][j] == -1)
      {
        Peak1D consensuspeak;
        consensuspeak.setMZ(y[ystart + j - 1].getMZ());
        consensuspeak.setIntensity((y[ystart + j - 1].getIntensity()) * factor_);
        if (keeppeaks_)
          lastconsensus_.push_back(consensuspeak);
        j--;
      }
      if (!i || !j)
        break;
    }

    return dparray[xend - xstart + 1][yend - ystart + 1];
  }

  const PeakSpectrum & SpectrumCheapDPCorr::lastconsensus() const
  {
    return lastconsensus_;
  }

  Map<UInt, UInt> SpectrumCheapDPCorr::getPeakMap() const
  {
    return peak_map_;
  }

  /**
   @param posa position of peak a
   @param posb position of peak b
   @param inta intensity of peak a
   @param intb intensity of peak b
   @return score
   */
  double SpectrumCheapDPCorr::comparepeaks_(double posa, double posb, double inta, double intb) const
  {
    double variation = (posa + posb) / 2 * (double)param_.getValue("variation");
    boost::math::normal_distribution<double> normal(0., variation);


    unsigned int int_cnt = (unsigned int)param_.getValue("int_cnt");
    if (int_cnt == 0)
    {
      double p = boost::math::pdf(normal, posa - posb);
      return p * inta * intb;
    }
    else
    {
      if (int_cnt == 1)
      {
        return boost::math::pdf(normal, posa - posb) * sqrt(inta * intb);
      }
      else
      {
        if (int_cnt == 2)
        {
          return boost::math::pdf(normal, posa - posb) * (inta + intb);
        }
        else
        {
          if (int_cnt == 3)
          {
            return max(0.0, boost::math::pdf(normal, posa - posb) * ((inta + intb) / 2 - fabs(inta - intb)));
          }
          else
          {
            // TODO exception
            // cerr << "int_cnt is not in [0,1,2,3]\n";
            return -1;
          }
        }
      }
    }
  }

}
