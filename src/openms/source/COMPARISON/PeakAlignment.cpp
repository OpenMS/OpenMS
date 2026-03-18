// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/PeakAlignment.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace std;

namespace OpenMS
{
  PeakAlignment::PeakAlignment() :
    PeakSpectrumCompareFunctor()
  {
    defaults_.setValue("epsilon", 0.2, "defines the absolute error of the mass spectrometer");
    defaults_.setValue("normalized", 1, "is set 1 if the similarity-measurement is normalized to the range [0,1]");
    defaults_.setValue("heuristic_level", 0, "set 0 means no heuristic is applied otherwise the given value is interpreted as unsigned integer, the number of strongest peaks considered for heurisitcs - in those sets of peaks has to be at least one match to conduct comparison");
    defaults_.setValue("precursor_mass_tolerance", 3.0, "Mass tolerance of the precursor peak, defines the distance of two PrecursorPeaks for which they are supposed to be from different peptides");
    defaultsToParam_();
  }

  PeakAlignment::PeakAlignment(const PeakAlignment& source) = default;

  PeakAlignment::~PeakAlignment() = default;

  PeakAlignment& PeakAlignment::operator=(const PeakAlignment& source)
  {
    if (this != &source)
    {
      PeakSpectrumCompareFunctor::operator=(source);
    }
    return *this;
  }

  double PeakAlignment::operator()(const PeakSpectrum& spec) const
  {
    return operator()(spec, spec);
  }

  double PeakAlignment::operator()(const PeakSpectrum& spec1, const PeakSpectrum& spec2) const
  {

    PeakSpectrum s1(spec1), s2(spec2);

    // shortcut similarity calculation by comparing PrecursorPeaks (PrecursorPeaks more than delta away from each other are supposed to be from another peptide)
    double pre_mz1 = 0.0;
    if (!spec1.getPrecursors().empty())
    {
      pre_mz1 = spec1.getPrecursors()[0].getMZ();
    }
    double pre_mz2 = 0.0;
    if (!spec2.getPrecursors().empty())
    {
      pre_mz2 = spec2.getPrecursors()[0].getMZ();
    }
    if (fabs(pre_mz1 - pre_mz2) > (double)param_.getValue("precursor_mass_tolerance"))
    {
      return 0;
    }

    // heuristic shortcut
    const double epsilon = (double)param_.getValue("epsilon");
    const UInt heuristic_level = (UInt)param_.getValue("heuristic_level");
    bool heuristic_filters(true);
    if (heuristic_level)
    {
      s1.sortByIntensity(true);
      s2.sortByIntensity(true);

      //heuristic filters (and shortcuts) if spec1 and spec2 have NOT at least one peak in the sets of |heuristic_level|-many highest peaks in common
      for (PeakSpectrum::ConstIterator it_s1 = s1.begin(); Size(it_s1 - s1.begin()) < heuristic_level && it_s1 != s1.end(); ++it_s1)
      {
        for (PeakSpectrum::ConstIterator it_s2 = s2.begin(); Size(it_s2 - s2.begin()) < heuristic_level && it_s2 != s2.end(); ++it_s2)
        {
          // determine if it is a match, i.e. mutual peak at certain m/z with epsilon tolerance
          if (fabs((*it_s2).getMZ() - (*it_s1).getMZ()) < epsilon)
          {
            heuristic_filters = false;
            break;
          }
        }
      }
    }
    if (heuristic_filters && heuristic_level)
    {
      return 0;
    }

    //TODO gapcost dependence on distance ?
    const double gap = (double)param_.getValue("epsilon");

    //initialize alignment matrix with 0 in (0,0) and a multiple of gapcost in the first row/col matrix(row,col,values)
    Matrix<double> matrix(spec1.size() + 1, spec2.size() + 1, 0);
    for (long int i = 1; i < matrix.rows(); i++)
    {
      matrix(i, 0) = -gap * i;
    }
    for (long int i = 1; i < matrix.cols(); i++)
    {
      matrix(0, i) = -gap * i;
    }

    //get sigma - the standard deviation (sqrt of variance)
    double mid(0);
    for (Size i = 0; i < spec1.size(); ++i)
    {
      for (Size j = 0; j < spec2.size(); ++j)
      {
        double pos1(spec1[i].getMZ()), pos2(spec2[j].getMZ());
        mid += fabs(pos1 - pos2);
      }
    }
    // average peak distance
    mid /= (spec1.size() * spec2.size());

    /* to manually retrace
    cout << "average peak distance " << mid << endl;
    */


    double var(0);
    for (Size i = 0; i < spec1.size(); ++i)
    {
      for (Size j = 0; j < spec2.size(); ++j)
      {
        double pos1(spec1[i].getMZ()), pos2(spec2[j].getMZ());
        var += (fabs(pos1 - pos2) - mid) * (fabs(pos1 - pos2) - mid);
      }
    }
    // peak distance variance
    var /= (spec1.size() * spec2.size());

    /* to manually retrace
    cout << "peak distance variance " << var << endl;
    */

    //only in case of only two equal peaks in the spectra sigma is 0


    const double sigma((var == 0) ? numeric_limits<double>::min() : sqrt(var));

    /* to manually retrace
    cout << "peak standard deviation " << sigma << endl;
    */

    //fill alignment matrix
    for (Size i = 1; i < spec1.size() + 1; ++i)
    {
      for (Size j = 1; j < spec2.size() + 1; ++j)
      {
        double pos1(spec1[i - 1].getMZ()), pos2(spec2[j - 1].getMZ());
        //only if peaks are in reasonable proximity alignment is considered else only gaps
        if (fabs(pos1 - pos2) <= epsilon)
        {
          // actual cell = max(upper left cell+score, left cell-gap, upper cell-gap)
          double from_left(matrix(i, j - 1) - gap);
          double from_above(matrix(i - 1, j) - gap);
          double int1(spec1[i - 1].getIntensity()), int2(spec2[j - 1].getIntensity());
          double from_diagonal(matrix(i - 1, j - 1) + peakPairScore_(pos1, int1, pos2, int2, sigma));
          matrix(i, j) = max(from_left, max(from_above, from_diagonal));
        }
        else
        {
          // actual cell = max(left cell-gap, upper cell-gap)
          double from_left(matrix(i, j - 1) - gap);
          double from_above(matrix(i - 1, j) - gap);
          matrix(i, j) = max(from_left, from_above);
        }
      }
    }

    /* to manually retrace
    cout << endl << matrix << endl;
    */

    //get best overall score and return
    double best_score(numeric_limits<double>::min());
    for (long int i = 0; i < matrix.cols(); i++)
    {
      best_score = max(best_score, matrix(matrix.rows() - 1, i));
    }
    for (long int i = 0; i < matrix.rows(); i++)
    {
      best_score = max(best_score, matrix(i, matrix.cols() - 1));
    }

    //calculate self-alignment scores for both input spectra
    double score_spec1(0), score_spec2(0);
    for (Size i = 0; i < spec1.size(); ++i)
    {
      double int_i(spec1[i].getIntensity());
      double pos_i(spec1[i].getMZ());
      score_spec1 += peakPairScore_(pos_i, int_i, pos_i, int_i, sigma);
    }
    for (Size i = 0; i < spec2.size(); ++i)
    {
      double int_i(spec2[i].getIntensity());
      double pos_i(spec2[i].getMZ());
      score_spec2 += peakPairScore_(pos_i, int_i, pos_i, int_i, sigma);
    }


    /* to manually retrace
    cout << "score_spec1: " << score_spec1 << "score_spec2: " << score_spec2 << endl;
    */

    //normalize score to interval [0,1] with geometric mean
    double best_score_normalized(best_score / sqrt(score_spec1 * score_spec2));

    /*
    cout << "score_spec1: " << score_spec1 << " score_spec2: " << score_spec2 <<  " best_score: " << best_score << endl;

    //normalize score to interval [0,1] with arithmetic mean
    double best_score_normalized( (best_score*2) / (score_spec1 + score_spec2) );
    */

    return best_score_normalized;
  }

  vector<pair<Size, Size> > PeakAlignment::getAlignmentTraceback(const PeakSpectrum& spec1, const PeakSpectrum& spec2) const
  {
    const double epsilon = (double)param_.getValue("epsilon");

    //TODO gapcost dependence on distance ?
    const double gap = (double)param_.getValue("epsilon");

    //initialize alignment matrix with 0 in (0,0) and a multiple of gapcost in the first row/col matrix(row,col,values)
    Matrix<double> matrix(spec1.size() + 1, spec2.size() + 1, 0);
    for (long int i = 1; i < matrix.rows(); i++)
    {
      matrix(i, 0) = -gap * i;
    }
    for (long int i = 1; i < matrix.cols(); i++)
    {
      matrix(0, i) = -gap * i;
    }

    // gives the direction of the matrix cell that originated the respective cell
    // e.g. matrix(i+1,j+1) could have originated from matrix(i,j), matrix(i+1,j) or matrix(i,j+1)
    // so traceback(i,j) represents matrix(i+1,j+1) and contains a "1"-from diagonal, a "0"-from left or a "2"-from above
    Matrix<Size> traceback(spec1.size(), spec2.size());

    //get sigma - the standard deviation (sqrt of variance)
    double mid(0);
    for (Size i = 0; i < spec1.size(); ++i)
    {
      for (Size j = 0; j < spec2.size(); ++j)
      {
        double pos1(spec1[i].getMZ()), pos2(spec2[j].getMZ());
        mid += fabs(pos1 - pos2);
      }
    }
    mid /= (spec1.size() * spec2.size());

    /* to manually retrace
        cout << mid << endl;
    */

    double var(0);
    for (Size i = 0; i < spec1.size(); ++i)
    {
      for (Size j = 0; j < spec2.size(); ++j)
      {
        double pos1(spec1[i].getMZ()), pos2(spec2[j].getMZ());
        var += (fabs(pos1 - pos2) - mid) * (fabs(pos1 - pos2) - mid);
      }
    }
    var /= (spec1.size() * spec2.size());

    /* to manually retrace
        cout << var << endl;
    */

    const double sigma(sqrt(var));

    /* to manually retrace
        cout << sigma << endl;
    */


    //fill alignment matrix
    for (Size i = 1; i < spec1.size() + 1; ++i)
    {
      for (Size j = 1; j < spec2.size() + 1; ++j)
      {
        double pos1(spec1[i - 1].getMZ()), pos2(spec2[j - 1].getMZ());
        //only if peaks are in reasonable proximity alignment is considered else only gaps
        if (fabs(pos1 - pos2) <= epsilon)
        {
          // actual cell = max(upper left cell+score, left cell-gap, upper cell-gap)
          double from_left(matrix(i, j - 1) - gap);
          double from_above(matrix(i - 1, j) - gap);
          double int1(spec1[i - 1].getIntensity()), int2(spec2[j - 1].getIntensity());
          double from_diagonal(matrix(i - 1, j - 1) + peakPairScore_(pos1, int1, pos2, int2, sigma));
          matrix(i, j) = max(from_left, max(from_above, from_diagonal));

          // TODO the cases where all or two values are equal
          if (from_diagonal > from_left && from_diagonal > from_above)
          {
            traceback(i - 1, j - 1) = 1;
          }
          else
          {
            if (from_left > from_diagonal && from_left > from_above)
            {
              traceback(i - 1, j - 1) = 0;
            }
            else
            {
              if (from_above > from_diagonal && from_above > from_left)
              {
                traceback(i - 1, j - 1) = 2;
              }
            }
          }
        }
        else
        {
          // actual cell = max(left cell-gap, upper cell-gap)
          double from_left(matrix(i, j - 1) - gap);
          double from_above(matrix(i - 1, j) - gap);
          matrix(i, j) = std::max(from_left, from_above);
          if (from_left > from_above)
          {
            traceback(i - 1, j - 1) = 0;
          }
          else           //from_left <= from_above
          {
            traceback(i - 1, j - 1) = 2;
          }
        }
      }
    }
    //return track from best alloverscore to 0,0
    vector<pair<Size, Size> > ret_val;

    //get matrix coordinates from best alloverscore
    Size row_index(0), col_index(0);
    double best_score(numeric_limits<double>::min());
    for (long int i = 0; i < matrix.cols(); i++)
    {
      if (best_score < matrix(matrix.rows() - 1, i))
      {
        best_score = matrix(matrix.rows() - 1, i);
        row_index = matrix.rows() - 1;
        col_index = i;
      }
    }
    for (long int i = 0; i < matrix.rows(); i++)
    {
      if (best_score < matrix(i, matrix.cols() - 1))
      {
        best_score = matrix(i, matrix.cols() - 1);
        row_index = i;
        col_index = matrix.cols() - 1;
      }
    }

    // TODO check the invariant!
    while (row_index > 0 && col_index > 0)
    {
      //from diagonal - peaks aligned
      if (traceback(row_index - 1, col_index - 1) == 1)
      {
        //register aligned peaks only
        ret_val.insert(ret_val.begin(), pair<Size, Size>(row_index - 1, col_index - 1));
        row_index = row_index - 1;
        col_index = col_index - 1;
      }
      // gap alignment
      else if (traceback(row_index - 1, col_index - 1) == 0)
      {
        col_index = col_index - 1;
      }
      else
      {
        row_index = row_index - 1;
      }
    }

    /* to manually retrace
    cout << endl << matrix << endl << traceback << endl;
    */

    return ret_val;
  }

  double PeakAlignment::peakPairScore_(double& pos1, double& intens1, double& pos2, double& intens2, const double& sigma) const
  {
    //scoring formula : peak intensity score * peak position score
    double pi(sqrt(intens1 * intens2));
    double pp((1 / (sigma * sqrt(2 * Constants::PI))) * exp(-(fabs(pos1 - pos2)) / 2 * sigma * sigma));

    /* to manually retrace
    cout << fabs(pos1-pos2) << " - "<< pi*pp << endl;
    */

    return pi * pp;
  }

}
