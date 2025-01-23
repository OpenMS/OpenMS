// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{
  /**
    @brief Computes the Savitzky-Golay filter coefficients using QR decomposition.

    This class represents a Savitzky-Golay lowpass-filter. The idea of the Savitzky-Golay filter
    is to find filtercoefficients that preserve higher moments, which means to approximate the underlying
    function within the moving window by a polynomial of higher order (typically quadratic or quartic).
    Therefore we least-squares fit for each data point a polynomial to all points \f$ f_i \f$ in the window
    and set \f$g_i\f$ to be the value of that polynomial at position \f$ i \f$. This method is superior
    to adjacent averaging because it tends to preserve features of the data such as peak height and width, which
    are usually 'washed out' by adjacent averaging.

    Because of the linearity of the problem, we can reduce the work computing by fitting in advance, for fictious data
    consisting of all zeros except for a singe 1 and then do the fits on the real data just by taking linear
    combinations. There are a particular sets of filter coefficients \f$ c_n \f$ which accomplish the process of
    polynomial least-squares fit inside a moving window. To get the symmetric coefficient-matrix
    \f$C \in R^{frameSize \times frameSize}\f$ with
    \f[ C= \left[ \begin{array}{cccc} c_{0,0} & c_{0,1} & \cdots & c_{0,frameSize-1} \\
                               \vdots  &         &         & \vdots            \\
             c_{frameSize-1,0} & c_{frameSize-1,2} & \ldots  & c_{frameSize-1,frameSize-1} \end{array} \right]\f]
    The first (last) \f$ \frac{frameSize}{2} \f$ rows of \f$ C \f$ we need to smooth the first (last)
    \f$ frameSize \f$ data points of the signal. So we use for the smoothing of the first data point the data
    point itself and the next \f$ frameSize-1 \f$ future points. For the second point we take the first datapoint,
    the data point itself and \f$ frameSize-2 \f$ of rightward data points... .
    We compute the Matrix \f$ C \f$ by solving the underlying least-squares problems with the singular value decomposition.
    Here we demonstrate the computation of the first row of a coefficient-matrix \f$ C \f$ for a Savitzky-Golay Filter
    of order=3 and frameSize=5:
    The design-matrix for the least-squares fit of a linear combination of 3 basis functions to 5 data points is:
    \f[ A=\left[ \begin{array}{ccc} x_0^0 & x_0^1 & x_0^2 \\ \\
                             x_1^0 & x_1^1 & x_1^2 \\ \\
           x_2^0 & x_2^1 & x_2^2 \\ \\
           x_3^0 & x_3^1 & x_3^2 \\ \\
           x_4^0 & x_4^1 & x_4^2 \end{array} \right]. \f]
    To smooth the first data point we have to create a design-matrix with \f$ x=[0,\ldots, frameSize-1] \f$.
    Now we have to solve the over-determined set of \f$ frameSize \f$ linear equations
    \f[ Ac=b \f]
    where \f$ b=[1,0,\ldots,0] \f$ represents the fictious data.
    Therefore we solve the normal equations of the least-squares problem
    \f[ A^TAc=A^Tb. \f]
    Now, it is possible to get
    \f[ c_n=\sum_{m=0}^8 \{(A^TA)^{-1}\}_{0,m} n^m, \f]
    with \f$ 0\le n \le 8\f$. Because we only need one row of the inverse matrix, it is possible to use LU decomposition with
    only a single backsubstitution.
    The vector \f$c=[c_0,\ldots,c_8] \f$ represents the wanted coefficients.
    Note that the solution of a least-squares problem directly from the normal equations is faster than the singular value
    decomposition but rather susceptible to roundoff error!

    @note This filter works only for uniform profile data!
          A polynomial order of 4 is recommended.
          The bigger the frame size the smoother the signal (the more detail information get lost!). The frame size corresponds to the number
          of filter coefficients, so the width of the smoothing interval is given by frame_size*spacing of the profile data.

        @note The data must be sorted according to ascending m/z!

        @htmlinclude OpenMS_SavitzkyGolayFilter.parameters

    @ingroup SignalProcessing
  */
  class OPENMS_DLLAPI SavitzkyGolayFilter :
    public ProgressLogger,
    public DefaultParamHandler
  {
public:
    /// Constructor
    SavitzkyGolayFilter();

    /// Destructor
    ~SavitzkyGolayFilter() override;

    // low level template to filters spectra and chromatograms
    // raw data and meta data needs to be copied to the output container before calling this function
    template<class InputIt, class OutputIt>
    void filter(InputIt first, InputIt last, OutputIt d_first)
    {
      size_t n = std::distance(first, last);

      if (frame_size_ > n) { return; }

      int i;
      UInt j;
      int mid = (frame_size_ / 2);
      double help;

      // compute the transient on
      OutputIt out_it = d_first;

      for (i = 0; i <= mid; ++i)
      {
        InputIt it_forward = (first - i);
        help = 0;
        for (j = 0; j < frame_size_; ++j)
        {
          help += it_forward->getIntensity() * coeffs_[(i + 1) * frame_size_ - 1 - j];
          ++it_forward;
        }

        out_it->setPosition(first->getPosition());
        out_it->setIntensity(std::max(0.0, help));
        ++out_it;
        ++first;
      }

      // compute the steady state output
      InputIt it_help = (last - mid);

      while (first != it_help)
      {
        InputIt it_forward = (first - mid);
        help = 0;

        for (j = 0; j < frame_size_; ++j)
        {
          help += it_forward->getIntensity() * coeffs_[mid * frame_size_ + j];
          ++it_forward;
        }

        out_it->setPosition(first->getPosition());
        out_it->setIntensity(std::max(0.0, help));
        ++out_it;
        ++first;
      }

      // compute the transient off
      for (i = (mid - 1); i >= 0; --i)
      {
        InputIt it_forward = (first - (frame_size_ - i - 1));
        help = 0;

        for (j = 0; j < frame_size_; ++j)
        {
          help += it_forward->getIntensity() * coeffs_[i * frame_size_ + j];
          ++it_forward;
        }

        out_it->setPosition(first->getPosition());
        out_it->setIntensity(std::max(0.0, help));
        ++out_it;
        ++first;
      }

    }

    /**
      @brief Removed the noise from an MSSpectrum containing profile data.
    */
    void filter(MSSpectrum & spectrum)
    {
      // copy the data AND META DATA to the output container
      MSSpectrum output = spectrum;
      // filter
      filter(spectrum.begin(), spectrum.end(), output.begin());
      // swap back
      std::swap(spectrum, output);
    }

    /**
      @brief Removed the noise from an MSChromatogram
    */
    void filter(MSChromatogram & chromatogram)
    {
      // copy the data AND META DATA to the output container
      MSChromatogram output = chromatogram;
      // filter
      filter(chromatogram.begin(), chromatogram.end(), output.begin());
      // swap back
      std::swap(chromatogram, output);
    }

    /**
      @brief Removed the noise from an MSExperiment containing profile data.
    */
    void filterExperiment(PeakMap & map)
    {
      Size progress = 0;
      startProgress(0, map.size() + map.getChromatograms().size(), "smoothing data");
      for (Size i = 0; i < map.size(); ++i)
      {
        filter(map[i]);
        setProgress(++progress);
      }
      for (Size i = 0; i < map.getChromatograms().size(); ++i)
      {
        filter(map.getChromatogram(i));
        setProgress(++progress);
      }
      endProgress();
    }

protected:
    /// Coefficients
    std::vector<double> coeffs_;

    /// UInt of the filter kernel (number of pre-tabulated coefficients)
    UInt frame_size_;

    /// The order of the smoothing polynomial.
    UInt order_;

    // Docu in base class
    void updateMembers_() override;
  };

} // namespace OpenMS
