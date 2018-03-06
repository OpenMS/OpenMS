// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors: Vipul Patel $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMSPECTRUMALIGNMENT_H
#define OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMSPECTRUMALIGNMENT_H

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{
  /**
      @brief A map alignment algorithm based on spectrum similarity (dynamic programming).

      @htmlinclude OpenMS_MapAlignmentAlgorithmSpectrumAlignment.parameters

      @experimental This algorithm is work in progress and might change.

      @ingroup MapAlignment
  */
  class OPENMS_DLLAPI MapAlignmentAlgorithmSpectrumAlignment :
    public DefaultParamHandler,
    public ProgressLogger
  {
public:
    /// Default constructor
    MapAlignmentAlgorithmSpectrumAlignment();

    /// Destructor
    ~MapAlignmentAlgorithmSpectrumAlignment() override;

    /// Align peak maps
    virtual void align(std::vector<PeakMap >&, std::vector<TransformationDescription>&);

private:
    /// Copy constructor is not implemented -> private
    MapAlignmentAlgorithmSpectrumAlignment(const MapAlignmentAlgorithmSpectrumAlignment&);
    /// Assignment operator is not implemented -> private
    MapAlignmentAlgorithmSpectrumAlignment & operator=(const MapAlignmentAlgorithmSpectrumAlignment &);

    /**
      @brief inner class necessary for using the sort algorithm.

      Defines several overloaded operator() for usage with the std::sort algorithm

      @experimental This algorithm is work in progress and might change.

    */
    class OPENMS_DLLAPI Compare
    {
    protected:
      bool flag;

    public:

      /// Default constructor with an order flag
      explicit Compare(bool b = false) :
        flag(b)
      {
      }

      /**
       * @brief overloaded operator() for comparing maps of maps std::pair<std::pair<Int,float>,float>.
       *
       * If the order flag is false, the second argument of the outer map is
       * selected. The output is an ascending order. If the order flag is true,
       * the first argument of the inner class is selected to get a descending
       * order.
       *
       */
      inline bool operator()(const std::pair<std::pair<Int, float>, float>& c1, const std::pair<std::pair<Int, float>, float>& c2)
      {
        if (!flag)
        {
          return c1.second > c2.second;
        }
        else
        {
          return (c1.first).first < (c2.first).first;
        }
      }

      /** 
       * @brief overloaded operator() for comparing pairs of float, float std::pair<float,float>.
       *
       * If the order flag is false, an ascending order are returned else a
       * descending. The comparison is done by the first argument of the map.
       *
       */
      inline bool operator()(const std::pair<float, float>& c1, const std::pair<float, float>& c2)
      {
        if (!flag)
        {
          return c1.first > c2.first;
        }
        else
        {
          return c1.first < c2.first;
        }
      }

    };

    /**
        @brief A function to prepare the sequence for the alignment. It calls intern the main function for the alignment.

        This function takes two arguments. These argument types are two MSExperiments.
        The first argument should have been filtered, so that only the type of MSLevel 1 exists in the Sequence.
        The second argument doesn't have to fulfill this restriction. It's going to be filtered automatically.
        With these two arguments a pre-calculation is done to find some corresponding data points(maximum 4) for building alignment blocks.
        After the alignment a re-transformation is done, the new Retention Times appear in the original data.

        The parameters are MSExperiments.

        @param pattern template map.
        @param aligned map which has to be aligned.
        @param transformation container for rebuilding the alignment only by specific data-points
    */
    void prepareAlign_(const std::vector<MSSpectrum*>& pattern, PeakMap& aligned, std::vector<TransformationDescription>& transformation);

    /**
        @brief filtered the MSLevel to gain only MSLevel 1
      
        The alignment works only on MSLevel 1 data, so a filter has to be run.

        @param peakmap map which has to be filtered
        @param spectrum_pointer_container output container, where pointers of the MSSpectrum are saved (only with MS level 1)

        @exception Exception::IllegalArgument is thrown if no spectra are contained in @p peakmap
    */
    void msFilter_(PeakMap& peakmap, std::vector<MSSpectrum*>& spectrum_pointer_container);

    /**
        @brief function for the test if cell i,j of the grid is inside the band

        The function returns true if the cell underlie these conditions:
        -k<=i-j<=k+n-m
        else return false.
        @param i coordinate i
        @param j coordinate j
        @param n size of column
        @param m size of row
        @param k_ size of k_
    */
    bool insideBand_(Size i, Size j, Size n, Size m, Int k_);

    /**
        @brief calculate the size of the band for the alignment for two given Sequence

        This function calculates the size of the band for the alignment. It
        takes three samples from the aligned sequence and tries to find the
        highscore pairs (matching against the template sequence). The highscore
        pair with the worst distance is to be chosen as the size of k.

        @param pattern vector of pointers of the template sequence
        @param aligned vector of pointers of the aligned sequence
        @param buffer holds the calculated score of index i,j.
        @param column_row_orientation indicate the order of the matrix
        @param xbegin indicate the beginning of the template sequence
        @param xend indicate the end of the template sequence
        @param ybegin indicate the beginning of the aligned sequence
        @param yend indicate the end of the aligned sequence
    */
    Int bestk_(const std::vector<MSSpectrum*>& pattern,
              std::vector<MSSpectrum*>& aligned, std::map<Size, std::map<Size, float> >& buffer,
              bool column_row_orientation, Size xbegin, Size xend, Size ybegin, Size yend);

    /**
        @brief calculate the score of two given MSSpectra calls intern scoring_

        This function calculates the score from two MSSpectra. These two
        MSSpectra are chosen by the coordinates i,j.  The two coordinates i,j
        indicate the index in the matrix. To find the right index on the
        sequence, each beginning is also given to the function.  A flag
        indicates the labeling of the axes. The buffermatrix stores the result
        of the scoring. If the band expands only a lookup of known scores is
        done.

        @param i is a index from the matrix.
        @param j is a index from the matrix.
        @param patternbegin indicate the beginning of the template sequence
        @param alignbegin  indicate the beginning of the aligned sequence
        @param pattern vector of pointers of the template sequence
        @param aligned vector of pointers of the aligned sequence
        @param buffer  holds the calculated score of index i,j.
        @param column_row_orientation indicate the order of the matrix
    */
    float scoreCalculation_(Size i, Size j, Size patternbegin, Size alignbegin,
                            const std::vector<MSSpectrum*>& pattern, std::vector<MSSpectrum*>& aligned,
                            std::map<Size, std::map<Size, float> >& buffer, bool column_row_orientation);

    /**
        @brief return the score of two given MSSpectra by calling the scorefunction
    */
    float scoring_(const MSSpectrum& a, MSSpectrum& b);

    /**
        @brief affine gap cost Alignment

        This Alignment is based on the Needleman Wunsch Algorithm.
        To improve the time complexity a banded version was implemented, known as k - alignment.
        To save some space, the alignment is going to be calculated by position xbegin to xend of one sequence and ybegin
        and yend by another given sequence. The result of the alignment is stored in the second argument.
        The first sequence is used as a template for the alignment.

        @param xbegin coordinate for the beginning of the template sequence.
        @param ybegin coordinate for the beginning of the aligned sequence .
        @param xend coordinate for the end of the template sequence.
        @param yend coordinate for the end of the aligned sequence.
        @param pattern template map.
        @param aligned map to be aligned.
        @param xcoordinate save the position of anchor points
        @param ycoordinate save the retentiontimes of an anchor points
        @param xcoordinatepattern save the reference position of the anchor points from the pattern

        @exception Exception::OutOfRange if a out of bound appear @p pattern or @p aligned
    */
    void affineGapalign_(Size xbegin, Size ybegin, Size xend, Size yend,
                        const std::vector<MSSpectrum*>& pattern,
                        std::vector<MSSpectrum*>& aligned,
                        std::vector<int>& xcoordinate, std::vector<float>& ycoordinate, 
                        std::vector<int>& xcoordinatepattern);

    /**
        @brief  preparation function of data points to construct later the spline function.

        This function reduced the amount of data values for the next step. The
        reduction is done by using a number of buckets, where the data points a
        selected.  Within the buckets, only defined number a selected, to be
        written back as a data point.  The selection within the buckets is done
        by scoring.

        @param pattern template map.
        @param aligned map to be aligned.
        @param xcoordinate save the position of anchor points
        @param ycoordinate  save the retention times of an anchor points
        @param xcoordinatepattern save the reference position of the anchor points from the pattern
    */
    void bucketFilter_(const std::vector<MSSpectrum*>& pattern,
                       std::vector<MSSpectrum*>& aligned, std::vector<Int>& xcoordinate,
                       std::vector<float>& ycoordinate, std::vector<Int>& xcoordinatepattern);

    /**
        @brief Creates files for the debugging

        This function is only active if the debug_ flag is true. The
        debugfileCreator creates following files:

        - debugtraceback.txt(gnuplotScript),
        - debugscoreheatmap.r and
        - debugRscript. 
        
        Debugscoreheatmap.r contains the scores of the Spectra to each other
        from the alignment and also the traceback. DebugRscript is the R script
        which reads those data. So both files are only working under R. Start R
        and type main(location of debugscoreheatmap.r). The output will be a
        heatmap of each sub-alignment. Debugtraceback.txt shows the way of the
        Traceback by using gnuplot.

        @param pattern template map.
        @param aligned map to be aligned.
    */
    void debugFileCreator_(const std::vector<MSSpectrum*>& pattern, std::vector<MSSpectrum*>& aligned);

    /**
        @brief Rounding the score of two spectra, only necessary for debugging

        This function rounded the score of two spectra. This is necessary for some function in the Debug-Mode
    */
    void debugscoreDistributionCalculation_(float score);

    ///Represent the gap cost for opening or closing a gap in the alignment
    float gap_;
    ///Extension cost after a gap is open
    float e_;
    ///Pointer holds the scoring function, which can be selected
    PeakSpectrumCompareFunctor* c1_;
    ///This is the minimal score to be count as a mismatch(range 0.0 - 1.0)
    float cutoffScore_;
    ///Defines the size of one bucket
    Size bucketsize_;
    ///Defines the amount of anchor points which are selected within one bucket.
    Size anchorPoints_;
    ///Debug mode flag default: False
    bool debug_;
    ///Represent the cost of a mismatch in the alignment
    float mismatchscore_;
    ///This is the minimum score for counting as a match(1-cutoffScore_)
    float threshold_;
    ///Container holding the score of the matchmatrix and also the insertmatrix
    std::vector<std::vector<float> > debugmatrix_;
    ///Container holding the only the score of Spectra
    std::vector<std::vector<float> > debugscorematrix_;
    ///Container holding the path of the traceback
    std::vector<std::pair<float, float> > debugtraceback_;
    ///Container holding the score of each cell(matchmatrix,insertmatrix, traceback)
    std::vector<float> scoredistribution_; //save the cell i, j , matchscore, insertscore, traceback
    //docu in base class
    void updateMembers_() override;
  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMSPECTRUMALIGNMENT_H
