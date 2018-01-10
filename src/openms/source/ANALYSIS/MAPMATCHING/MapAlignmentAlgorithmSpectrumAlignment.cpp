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

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmSpectrumAlignment.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/CONCEPT/Factory.h>

#include <fstream>

namespace OpenMS
{

  MapAlignmentAlgorithmSpectrumAlignment::MapAlignmentAlgorithmSpectrumAlignment() :
    DefaultParamHandler("MapAlignmentAlgorithmSpectrumAlignment"),
    ProgressLogger(), c1_(nullptr)
  {
    defaults_.setValue("gapcost", 1.0, "This Parameter stands for the cost of opining a gap in the Alignment. A gap means that one spectrum can not be aligned directly to another spectrum in the Map. This happens, when the similarity of both spectra a too low or even not present. Imagine it as a insert or delete of the spectrum in the map (similar to sequence alignment). The gap is necessary for aligning, if we open a gap there is a possibility that an another spectrum can be correct aligned with a higher score as before without gap. But to open a gap is a negative event and needs to carry a punishment, so a gap should only be opened if the benefits outweigh the downsides. The Parameter is to giving as a positive number, the implementation convert it to a negative number.");
    defaults_.setMinFloat("gapcost", 0.0);
    defaults_.setValue("affinegapcost", 0.5, "This Parameter controls the cost of extension a already open gap. The idea behind the affine gapcost lies under the assumption, that it is better to get a long distance of connected gaps than to have a structure of gaps interspersed with matches (gap match gap match etc.).  Therefor the punishment for the extension of a gap generally should be lower than the normal gapcost. If the result of the alignment shows high compression, it is a good idea to lower either the affine gapcost or gap opening cost.");
    defaults_.setMinFloat("affinegapcost", 0.0);
    defaults_.setValue("cutoff_score", 0.70, "The Parameter defines the threshold which filtered spectra, these spectra are high potential candidate for deciding the interval of a sub-alignment.  Only those pair of spectra are selected, which has a score higher or same of the threshold.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("cutoff_score", 0.0);
    defaults_.setMaxFloat("cutoff_score", 1.0);
    defaults_.setValue("bucketsize", 100, "Defines the numbers of buckets. It is a quantize of the interval of those points, which defines the main alignment (match points). These points have to filtered, to reduce the amount of points for the calculating a smoother spline curve.", ListUtils::create<String>("advanced"));
    defaults_.setMinInt("bucketsize", 1);
    defaults_.setValue("anchorpoints", 100, "Defines the percent of numbers of match points which a selected from one bucket. The high score pairs are previously selected. The reduction of match points helps to get a smoother spline curve.", ListUtils::create<String>("advanced"));
    defaults_.setValue("debug", "false", "Activate the debug mode, there a files written starting with debug prefix.", ListUtils::create<String>("advanced"));
    defaults_.setMinInt("anchorpoints", 1);
    defaults_.setMaxInt("anchorpoints", 100);
    defaults_.setValidStrings("debug", ListUtils::create<String>("true,false"));
    defaults_.setValue("mismatchscore", -5.0, "Defines the score of two spectra if they have no similarity to each other. ", ListUtils::create<String>("advanced"));
    defaults_.setMaxFloat("mismatchscore", 0.0);
    defaults_.setValue("scorefunction", "SteinScottImproveScore", "The score function is the core of an alignment. The success of an alignment depends mostly of the elected score function. The score function return the similarity of two spectra. The score influence defines later the way of possible traceback. There are multiple spectra similarity scores available..");
    defaults_.setValidStrings("scorefunction", ListUtils::create<String>("SteinScottImproveScore,ZhangSimilarityScore")); //Factory<PeakSpectrumCompareFunctor>::registeredProducts());
    defaultsToParam_();
    setLogType(CMD);
  }

  MapAlignmentAlgorithmSpectrumAlignment::~MapAlignmentAlgorithmSpectrumAlignment()
  {
    delete c1_;
  }

  void MapAlignmentAlgorithmSpectrumAlignment::align(std::vector<PeakMap >& peakmaps, std::vector<TransformationDescription>& transformation)
  {
    transformation.clear();
    TransformationDescription trafo;
    trafo.fitModel("identity");
    transformation.push_back(trafo); // transformation of reference map
    try
    {
      std::vector<MSSpectrum*> spectrum_pointers;
      msFilter_(peakmaps[0], spectrum_pointers);
      startProgress(0, (peakmaps.size() - 1), "Alignment");
      for (Size i = 1; i < peakmaps.size(); ++i)
      {
        prepareAlign_(spectrum_pointers, peakmaps[i], transformation);
        setProgress(i);
      }
      endProgress();
    }
    catch (Exception::OutOfRange& /*e*/)
    {
      throw Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
  }

  void MapAlignmentAlgorithmSpectrumAlignment::prepareAlign_(const std::vector<MSSpectrum*>& pattern, PeakMap& aligned, std::vector<TransformationDescription>& transformation)
  {
    //tempalign ->container for holding only MSSpectrums with MS-Level 1
    std::vector<MSSpectrum*> tempalign;
    msFilter_(aligned, tempalign);

    //if it's possible, built 4 blocks. These can be individually be aligned.
    std::vector<Size> alignpoint;
    //saving the first coordinates
    alignpoint.push_back(0);
    alignpoint.push_back(0);
    //4 blocks : 0-0.25 ,0.25-50,0.50-0.75,1 The data points must have a high similarity score
    for (float i = 0.25; i <= 0.75; i += 0.25)
    {
      Size y = (Size)(tempalign.size() * i);
      Size x = 0;
      float maxi = -999.0;

      for (Size k = 0; k < pattern.size(); ++k)
      {
        float s =    scoring_(*pattern[k], *(tempalign[y]));
        if (s > maxi && s > cutoffScore_)
        {
          x = k;
          maxi = s;
        }
      }
      if (x >= alignpoint[alignpoint.size() - 2] + 3 && y >= alignpoint[alignpoint.size() - 1] + 3)
      {
        alignpoint.push_back(x);
        alignpoint.push_back(y);
      }

      Size xn = (Size)(pattern.size() * i);
      Size yn = 0;
      for (Size k = 0; k < tempalign.size(); ++k)
      {
        float s =    scoring_(*pattern[xn], *(tempalign[k]));
        if (s > maxi && s > cutoffScore_)
        {
          yn = k;
          maxi = s;
        }
      }
      if (xn >= alignpoint[alignpoint.size() - 2] + 3 && yn >= alignpoint[alignpoint.size() - 1] + 3)
      {
        alignpoint.push_back(xn);
        alignpoint.push_back(yn);
      }
      //only save possible data points, if they are not already contained
    }
    //save also the endpoint as a data point
    alignpoint.push_back(pattern.size() - 1);
    alignpoint.push_back(tempalign.size() - 1);
    //the distance of two data points have to be greater than 3, if not the spline would thrown an Exception
    //do a affine gap alignment of the block of the data points x1,y1,x2,y2

    std::vector<Int> xcoordinate;
    std::vector<Int> xcoordinatepattern;
    std::vector<float> ycoordinate;
    debugmatrix_.clear();

    for (Size i = 0; i < alignpoint.size() - 2; i += 2)
    {
      affineGapalign_(alignpoint[i], alignpoint[i + 1], alignpoint[i + 2], alignpoint[i + 3], pattern, tempalign, xcoordinate, ycoordinate, xcoordinatepattern);
    }

    //affineGapalign_(0,0,(pattern.size()-1),(tempalign.size()-1), pattern, tempalign,xcoordinate,ycoordinate,xcoordinatepattern);
    if (debug_)
    {
      debugFileCreator_(pattern, tempalign);
    }

    /*
        for (Size i = 0; i< xcoordinate.size(); ++i)
        {
            std::cout<< xcoordinate[i] << " " << ycoordinate[i] << " x  y  anchorpunkte " << std::endl;
        }
        std::cout << std::endl;
        */
    bucketFilter_(pattern, tempalign, xcoordinate, ycoordinate, xcoordinatepattern);
    /*std::cout << xcoordinate.size()<< std::endl;
            for (Size i = 0; i< xcoordinate.size(); ++i)
            {
                std::cout<< xcoordinate[i] << " " << ycoordinate[i] << " x  y  anchorpunkte " << std::endl;
            }*/

    // store the data points defining the transformation:
    TransformationDescription::DataPoints data;
    for (Size i = 0; i < xcoordinate.size(); ++i)
    {
      double rt = tempalign[xcoordinate[i]]->getRT();
      data.push_back(std::make_pair(rt, double(ycoordinate[i])));
    }
    transformation.push_back(TransformationDescription(data));
  }

  void MapAlignmentAlgorithmSpectrumAlignment::affineGapalign_(Size xbegin, Size ybegin, Size xend, Size yend, const std::vector<MSSpectrum*>& pattern, std::vector<MSSpectrum*>& aligned, std::vector<int>& xcoordinate, std::vector<float>& ycoordinate, std::vector<int>& xcoordinatepattern)
  {
    //affine gap alignment needs two matrices
    std::map<Size, std::map<Size, float> > firstcolummatchmatrix;
    std::map<Size, std::map<Size, float> > secondcolummatchmatrix;
    Size n = std::max((xend - xbegin), (yend - ybegin)) + 1; //column
    Size m = std::min((xend - xbegin), (yend - ybegin)) + 1; //row
    // std::cout<< n << " n " << m << " m " <<  xbegin << " " <<xend<< " " << ybegin << " " << yend <<std::endl;
    //log the Progress of the subaligmnet
    String temp = "sub-alignment of interval: template sequence " + String(xbegin) + " " + String(xend) + " interval: alignsequence " + String(ybegin) + " " + String(yend);
    startProgress(0, n, temp);

    bool column_row_orientation = false;
    if (n != (xend - xbegin) + 1)
    {
      column_row_orientation = true;
    }
    //matrix for holding calculated scores
    std::map<Size, std::map<Size, float> > buffermatrix;
    std::map<Size, std::map<Size, Size> > traceback;
    //calculate the value of k
    Int k_ = bestk_(pattern, aligned, buffermatrix, column_row_orientation, xbegin, xend, ybegin, yend) + 2;
    float score_ = -99999999.0f;
    //flag if we have to calculate again the alignment in step k+1
    bool finish = false;
    while (!finish)
    {
      traceback.clear();
      for (Size i = 0; i <= n; ++i)
      {
        setProgress(i);
        for (Size j = 0; j <= m; ++j) //if( j >=1 && (Size)j<=m)
        {
          if (insideBand_(i, j, n, m, k_))
          {
            if (i == 0 || j == 0)
            {
              if (i == 0)
              {
                firstcolummatchmatrix[i][j] = (-gap_) * j;
                //std::cout << firstcolummatchmatrix[i][j] << " i j firstcolum "<< std::endl;
              }
              else if (j == 0)
              {
                secondcolummatchmatrix[i][j] = (-gap_) * i;
                //std::cout << secondcolummatchmatrix[i][j] << " i j secondcolum "<< std::endl;
              }
            }
            else
            {
              try
              {
                double s = scoreCalculation_(i, j, xbegin, ybegin, pattern, aligned, buffermatrix, column_row_orientation);
                if (debug_)
                {
                  std::vector<float> ltemp;
                  if (!column_row_orientation)
                  {
                    ltemp.push_back((float)i + xbegin - 1);
                    ltemp.push_back((float)j + ybegin - 1);
                    ltemp.push_back(s);
                    ltemp.push_back(0);
                    debugscorematrix_.push_back(ltemp);
                  }
                  else
                  {
                    ltemp.push_back((float)j + xbegin - 1);
                    ltemp.push_back((float)i + ybegin - 1);
                    ltemp.push_back(s);
                    ltemp.push_back(0);
                    debugscorematrix_.push_back(ltemp);
                  }
                }
                float mv = -999.0;
                float mh = -999.0;
                i = i - 1;
                if (insideBand_(i, j, n, m, k_))
                {
                  mh = firstcolummatchmatrix[i][j] - gap_;
                  //std::cout <<firstcolummatchmatrix[i][j] << " " << i << " " << j<<" firstcolumnnn "<<std::endl;
                  i += 1;
                }
                else
                {
                  i += 1;
                }
                j = j - 1;
                if (insideBand_(i, j, n, m, k_))
                {
                  mv = secondcolummatchmatrix[i][j] - gap_;
                  j += 1;
                }
                else
                {
                  j += 1;
                }
                //mv=matchmatrix[i][j-1]-gap_;
                //mh=matchmatrix[i-1][j]-gap_;
                float md = firstcolummatchmatrix[i - 1][j - 1] + s;
                //std::cout << i << " " << j << " i j " << mv << " mv " << mh << " mh " << md << " md " << std::endl;
                secondcolummatchmatrix[i][j] = std::max((float)md, std::max((float)(mv), (float)(mh)));
                //std::cout << secondcolummatchmatrix[i][j] << " " << i << " "<< j << " i j zweiter colum" << std::endl;
                if (secondcolummatchmatrix[i][j] == mh)
                {
                  traceback[i][j] = 1;
                }
                else if (secondcolummatchmatrix[i][j] == mv)
                {
                  traceback[i][j] = 2;
                }
                else
                  traceback[i][j] = 0;
              }
              catch (Exception::OutOfRange /*&e*/)
              {
                throw Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
              }
            }
          }
        }
        if (i != 0)
        {
          firstcolummatchmatrix = secondcolummatchmatrix;
          secondcolummatchmatrix.clear();
        }
      }
      if (score_ >= firstcolummatchmatrix[n][m] || k_ == (Int)n + 2) // || matchmatrix[n][m] >= (2*(k_+1)+n-m)*(-gap_)+(n-(k_+1))*3/*(m-(k_+1))*cutoffScore_ -2* gap_ - (2*(k_) +((Int)m-(Int)n))*e_*/)
      {
        finish = true;
        firstcolummatchmatrix.clear();
        secondcolummatchmatrix.clear();
      }
      else
      {
        score_ = firstcolummatchmatrix[n][m];
        k_ *= 2;
        if (k_ > (Int)n + 2)
          k_ = (Int)n + 2;
      }
    }
    // matchmatrix.clear();
    /*for(Size i=0; i <=n;++i)
        {
            for(Size j=0; j<=m;++j)
            {
                if(insideBand_(i,j,n,m,k_))
                std::cout << i << " "<< j << " i j " << traceback[i][j] << " Match  Trace" << std::endl;
            }
        std::cout<< std::endl;
        }*/
    //traceback
    bool endtraceback = false;
    int i = (int) n;
    int j = (int) m;
    //float maximum = -999.0;
    //container necessary for collecting the positions of both sequence to gain later the correct datapoints for the spline
    std::vector<int> xvar;
    std::vector<int> xxvar;
    std::vector<double> yvar;

    while (!endtraceback)
    {
      //std::cout << i << " " << j << " i j " << std::endl;
      if (i <= 0 || j <= 0)
        endtraceback = true;
      else
      {
        if (traceback[i][j] == 0)
        {
          if (!column_row_orientation)
          {
            if (debug_)
            {
              debugtraceback_.push_back(std::make_pair(float(i + xbegin - 1), float(j + ybegin - 1)));
            }
            xvar.push_back(j + (int)ybegin - 1);
            yvar.push_back((*pattern[i + xbegin - 1]).getRT());
            xxvar.push_back(i + (int)xbegin - 1);
          }
          else
          {
            if (debug_)
            {
              debugtraceback_.push_back(std::make_pair(float(j + xbegin - 1), float(i + ybegin - 1)));
            }
            xvar.push_back(i + (int)ybegin - 1);
            yvar.push_back((*pattern[j + xbegin - 1]).getRT());
            xxvar.push_back(j + (int)xbegin - 1);
          }
          i = i - 1;
          j = j - 1;
        }
        else if (traceback[i][j] == 1)
          i = i - 1;
        else if (traceback[i][j] == 2)
          j = j - 1;
      }
    }

    for (Size k = 0; k < xvar.size(); ++k)
    {
      if (xcoordinate.size() > 0)
      {
        if (xvar[xvar.size() - 1 - k] != xcoordinate[xcoordinate.size() - 1])
        {
          xcoordinate.push_back(xvar[xvar.size() - 1 - k]);
          ycoordinate.push_back(yvar[yvar.size() - 1 - k]);
          xcoordinatepattern.push_back(xxvar[xxvar.size() - 1 - k]);
        }
      }
      else
      {
        xcoordinate.push_back(xvar[xvar.size() - 1 - k]);
        ycoordinate.push_back(yvar[yvar.size() - 1 - k]);
        xcoordinatepattern.push_back(xxvar[xxvar.size() - 1 - k]);
      }
    }
    //std::cout<< xcoordinate.size()<< std::endl;
    endProgress();
  }

  void MapAlignmentAlgorithmSpectrumAlignment::msFilter_(PeakMap& peakmap, std::vector<MSSpectrum*>& spectrum_pointer_container)
  {
    std::vector<UInt> pattern;
    peakmap.updateRanges(-1);
    pattern = peakmap.getMSLevels();

    if (!pattern.empty())
    {
      for (Size i = 0; i < peakmap.size(); ++i)
      {
        if (peakmap[i].getMSLevel() == 1)
        {
          spectrum_pointer_container.push_back(&(peakmap[i]));
        }
      }
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No spectra contained");
    }
  }

  inline bool MapAlignmentAlgorithmSpectrumAlignment::insideBand_(Size i, Size j, Size n, Size m, Int k_)
  {
    if ((Int)(m - n - k_) <= (Int)(i - j) && (Int) (i - j) <= k_) //  if((Int)(-k_)<=(Int)(i-j) &&(Int) (i-j) <=(Int)(k_+n-m))
    {
      //std::cout << i << " i " << j << " j " << " innerhalb der Bande " << std::endl;
      return true;
    }
    else
    {
      //std::cout << i << " i " << j << " j " << "NICHT innerhalb der Bande " << std::endl;
      return false;
    }
  }

  inline Int MapAlignmentAlgorithmSpectrumAlignment::bestk_(const std::vector<MSSpectrum*>& pattern, std::vector<MSSpectrum*>& aligned, std::map<Size, std::map<Size, float> >& buffer, bool column_row_orientation, Size xbegin, Size xend, Size ybegin, Size yend)
  {
    Int ktemp = 2;
    for (float i = 0.25; i <= 0.75; i += 0.25)
    {
      Size temp = (Size)((yend - ybegin) * i);
      float    maxi = -999.0;
      float    s = -999.0;
      for (Size k = 0; k <= (xend - xbegin); ++k)
      {
        Size x;
        Int y;
        if (column_row_orientation)
        {
          x = temp + 1;
          y = (Int)k + 1;
          s = scoreCalculation_(x, y, xbegin, ybegin, pattern, aligned, buffer, column_row_orientation);
        }
        else
        {
          x = k + 1;
          y = (Int)(temp + 1);
          s = scoreCalculation_(x, y, xbegin, ybegin, pattern, aligned, buffer, column_row_orientation);
        }
        if (s > maxi && s > cutoffScore_)
        {
          maxi = s;
          if (ktemp < std::abs((Int)x - y) + 1)
          {
            ktemp = std::abs((Int)x - y) + 1;
          }
        }
      }
    }
    return ktemp;
  }

  inline float MapAlignmentAlgorithmSpectrumAlignment::scoreCalculation_(Size i, Size j, Size patternbegin, Size alignbegin, const std::vector<MSSpectrum*>& pattern, std::vector<MSSpectrum*>& aligned, std::map<Size, std::map<Size, float> >& buffer, bool column_row_orientation)
  {
    if (!column_row_orientation)
    {
      if (buffer[i][j] == 0)
      {
        float score = scoring_(*pattern[i + patternbegin - 1], *aligned[j + alignbegin - 1]);
        if (score > 1)
          score = 1;
        if (debug_)
        {
          debugscoreDistributionCalculation_(score);
        }
        if (score < threshold_)
          score = mismatchscore_;
        else
          score = 2 + score;
        buffer[i][j] = score;
      }
      return buffer[i][j];
    }
    else
    {
      if (buffer[j][i] == 0)
      {
        float score = scoring_(*pattern[j + patternbegin - 1], *aligned[i + alignbegin - 1]);
        if (score > 1)
          score = 1;
        if (debug_)
        {
          debugscoreDistributionCalculation_(score);
        }
        if (score < threshold_)
          score = mismatchscore_;
        else
          score = 2 + score;
        buffer[j][i] = score;
      }
      return buffer[j][i];
    }
  }

  inline float MapAlignmentAlgorithmSpectrumAlignment::scoring_(const MSSpectrum& a, MSSpectrum& b)
  {
    return c1_->operator()(a, b);
  }

  inline void MapAlignmentAlgorithmSpectrumAlignment::bucketFilter_(const std::vector<MSSpectrum*>& pattern, std::vector<MSSpectrum*>& aligned, std::vector<int>& xcoordinate, std::vector<float>& ycoordinate, std::vector<int>& xcoordinatepattern)
  {
    std::vector<std::pair<std::pair<Int, float>, float> > tempxy;
    Size size = 0;
    //std::cout <<bucketsize_  << " bucketsize " <<xcoordinate.size() << " xsize()" << std::endl;
    if (bucketsize_ >= xcoordinate.size())
    {
      bucketsize_ = xcoordinate.size() - 1;
      size = 1;
    }
    else
      size = xcoordinate.size() / bucketsize_;

    if (size == 1)
      bucketsize_ = xcoordinate.size() - 1;
    //std::cout << size << " size "<< xcoordinate.size() << " xcoordinate.size() " << std::endl;
    for (Size i = 0; i < size; ++i)
    {
      std::vector<std::pair<std::pair<Int, float>, float> > temp;
      for (Size j = 0; j < bucketsize_; ++j)
      {
        //std::cout<< j << " j " << std::endl;
        float score = scoring_(*pattern[xcoordinatepattern[(i * bucketsize_) + j]], *aligned[xcoordinate[(i * bucketsize_) + j]]);
        //modification only view as a possible data point if the score is higher than 0
        if (score >= threshold_)
        {
          temp.push_back(std::make_pair(std::make_pair(xcoordinate[(i * bucketsize_) + j], ycoordinate[(i * bucketsize_) + j]), score));
        }
      }
      /*for(Size i=0; i < temp.size();++i)
          {
          std::cout<< (temp[i].first).first << " " << (temp[i].first).second << " in temp"<< std::endl;
          }
          std::cout << std::endl;
      */
      std::sort(temp.begin(), temp.end(), Compare(false));
      //Int anchor=(Int)(size*anchorPoints_/100);
      float anchor = (temp.size() * anchorPoints_ / 100);
      if (anchor <= 0 && !temp.empty())
      {
        anchor = 1;
      }

      //std::cout << anchor << " anchorpoints "  << anchorPoints_<< std::endl;
      /*for(UInt i=0; i< temp.size();++i)
          {
          std::cout<< (temp[i].first).first << "first" << (temp[i].first).second << " second" <<std::endl;
          }
      */
      for (Size k = 0; k < (Size)anchor; ++k)
      {
        tempxy.push_back(temp[k]);
      }
    }
    std::sort(tempxy.begin(), tempxy.end(), Compare(true));
    xcoordinate.clear();
    ycoordinate.clear();
    for (Size i = 0; i < tempxy.size(); ++i)
    {
      if (i != 0)
      {
        if (xcoordinate[xcoordinate.size() - 1] != (tempxy[i].first).first)
        {
          xcoordinate.push_back((tempxy[i].first).first);
          ycoordinate.push_back((tempxy[i].first).second);
        }
      }
      else
      {
        xcoordinate.push_back((tempxy[i].first).first);
        ycoordinate.push_back((tempxy[i].first).second);
      }
    }
    /*
    for(Size i=0; i < xcoordinate.size();++i)
    {
        std::cout << xcoordinate[i] << " xcoordinate "<< std::endl;
    }
    */
  }

  inline void MapAlignmentAlgorithmSpectrumAlignment::debugFileCreator_(const std::vector<MSSpectrum*>& pattern, std::vector<MSSpectrum*>& aligned)
  {
    //plotting scores of the alignment
    /*std::ofstream tempfile3;
    float maximum=2.0;
    tempfile3.open("debugscore.txt",std::ios::trunc);
    tempfile3 << "set xrange[0:"<< pattern.size()-1<<  "]" << "\n set yrange[0:"<< aligned.size()-1 << "]" << "\n set zrange[0:"
    << maximum << "] \n set view 45,20,1.0,2.5 \n"<< "splot \'-\'" <<std::endl;
    for (Size i =0; i< debugscorematrix_.size();++i)
    {
        tempfile3<< debugscorematrix_[i][0] << " " << debugscorematrix_[i][1] << " " << debugscorematrix_[i][2] << std::endl;
    }
    tempfile3 << "e" << std::endl;
    tempfile3.close();*/

    //gnuplot of the traceback
    std::ofstream myfile;
    myfile.open("debugtraceback.txt", std::ios::trunc);
    myfile << "set xrange[0:" << pattern.size() - 1 <<  "]" << "\n set yrange[0:" << aligned.size() - 1 << "] \n plot \'-\' with lines " << std::endl;
    std::sort(debugtraceback_.begin(), debugtraceback_.end(), Compare(false));

    for (Size i = 0; i < debugtraceback_.size(); ++i)
    {
      myfile << debugtraceback_[i].first << " " << debugtraceback_[i].second << std::endl;
      for (Size p = 0; p < debugscorematrix_.size(); ++p)
      {
        if (debugscorematrix_[p][0] == debugtraceback_[i].first && debugscorematrix_[p][1] == debugtraceback_[i].second)
        {
          debugscorematrix_[p][3] = 1;
          break;
        }
      }
    }
    myfile << "e" << std::endl;
    myfile.close();
    //R heatplot score of both sequence
    // std::map<Size, std::map<Size, float> > debugbuffermatrix;

    float scoremaximum = -2;
    //precalculation for the heatmap
    //getting the maximum score
    for (Size i = 0; i < debugscorematrix_.size(); ++i)
    {
      if (scoremaximum < debugscorematrix_[i][2] + 2)
        scoremaximum = debugscorematrix_[i][2] + 2;
      //shift all score about 2 (to get 0, the default score in the debugbuffermatrix is -2 )
      debugscorematrix_[i][2] += 2;
    }
    //to get the intvall [0,1] divide all score to the global maximum
    for (Size i = 0; i < debugscorematrix_.size(); ++i)
    {
      if (debugscorematrix_[i][2] != 0)
        debugscorematrix_[i][2] /= scoremaximum;
    }
    //write the score in a file
    /*
    for (Size i=0; i < debugscorematrix_.size();++i)
    {
        debugbuffermatrix[(UInt)debugscorematrix_[i][0]][(UInt)debugscorematrix_[i][1]]=debugscorematrix_[i][2];
    }
    */
    std::ofstream scorefile;
    scorefile.open("debugscoreheatmap.r", std::ios::trunc);
    /*
    for (Size i=0; i < debugbuffermatrix.size();++i)
    {
        for (Size j=0; j< debugbuffermatrix[i].size(); ++j)
        {
        scorefile<< i << " "<< j << " "<<  debugbuffermatrix[i][j] << std::endl;
        }
    }
    */
    for (Size i = 0; i < debugscorematrix_.size(); ++i)
    {
      scorefile << debugscorematrix_[i][0] << " " << debugscorematrix_[i][1] << " " << debugscorematrix_[i][2] << " " << debugscorematrix_[i][3] << std::endl;
    }
    scorefile.close();

    std::ofstream rscript;
    rscript.open("debugRscript.r", std::ios::trunc);

    rscript << "#Name: LoadFile \n #transfer data from file into a matrix \n #Input: Filename \n #Output Matrix \n LoadFile<-function(fname){\n temp<-read.table(fname); \n temp<-as.matrix(temp); \n return(temp); \n } " << std::endl;
    rscript << "#Name: ScoreHeatmapPlot \n #plot the score in a way of a heatmap \n #Input: Scorematrix \n #Output Heatmap \n ScoreHeatmapPlot<-function(matrix) { \n xcord<-as.vector(matrix[,1]); \n ycord<-as.vector(matrix[,2]); \n color<-rgb(as.vector(matrix[,4]),as.vector(matrix[,3]),0);\n  plot(xcord,ycord,col=color, main =\"Heatplot of scores included the traceback\" , xlab= \" Template-sequence \", ylab=\" Aligned-sequence \", type=\"p\" ,phc=22)\n } \n main<-function(filenamea) { \n a<-Loadfile(filenamea) \n X11() \n ScoreHeatmapPlot(a) \n  " << std::endl;
    rscript.close();
    /*
    float matchmaximum=-999.0;
    float insertmaximum=-999.0;
    for (Size i =0; i< debugmatrix_.size();++i)
    {
        debugbuffermatrix[debugmatrix_[i][0]][debugmatrix_[i][1]]=(float)debugmatrix_[i][2];
        if(matchmaximum <debugmatrix_[i][2]) matchmaximum = (float)debugmatrix_[i][2];
            //debuginsertmatrix[debugmatrix_[i][0]][debugmatrix_[i][1]]=debugmatrix_[i][3];
        if(insertmaximum <debugmatrix_[i][3]) insertmaximum = (float)debugmatrix_[i][3];
            //debugtracebackmatrix[debugmatrix_[i][0]][debugmatrix_[i][1]]=debugmatrix_[i][4];
    }

    //todo option for gnuplot
    myfile<< "#plotting matchmatrix" <<std::endl;
    myfile<< "set multiplot layout 1,1 columnsfirst" << std::endl;
    myfile<< "set title \"Heat Map of the matchmatrix\"  \n "
    "unset key \n"
    "set tic scale 0 \n set palette rgbformula -7,2,-7 \n set cbrange [-999.0:" << matchmaximum <<"] \n set cblabel \"Score\" \n unset cbtics"<< std::endl;
    myfile<< "p \'-\' using 1:2:3 with image" <<std::endl;
    //output matchmatrix
    for (Size i=0; i< debugbuffermatrix.size();++i)
    {
        if(i!=0) myfile << std::endl;
        for (Size j=0; j < debugbuffermatrix[0].size();++j)
        {
        myfile<< i << " " << j << " " << debugbuffermatrix[i][j] << std::endl;
        }
    }
    myfile <<"e"<< std::endl;
    myfile << std::endl;
    std::cout << " complete matchmatrix " << std::endl;
    for (Size i =0; i< debugmatrix_.size();++i)
    {
        debugbuffermatrix[debugmatrix_[i][0]][debugmatrix_[i][1]]=debugmatrix_[i][3];
    }
    myfile<< "set title \"Heat Map of the insertmatrix\"  \n "
    "unset key \n"
    "set tic scale 0 \n set palette rgbformula -7,2,-7 \n set cbrange [-999.0:"<< insertmaximum<<"] \n set cblabel \"Score\" \n unset cbtics"<< std::endl;
    myfile<< "p \'-\' using 1:2:3 with image" <<std::endl;
    for (Size i=0; i< debugbuffermatrix.size();++i)
    {
        if(i!=0) myfile << std::endl;
        for (Size j=0; j < debugbuffermatrix[0].size();++j)
        {
        myfile<< i << " " << j << " " << debugbuffermatrix[i][j] << std::endl;
        }
    }
    myfile <<"e"<< std::endl;
    std::cout << "complete insertmatrix " << std::endl;
    myfile << std::endl;

    myfile<< "set multiplot layout 1,1 columnsfirst \n set title \"Heat Map of the tracebackmatrix\" \n unset key \n set tic scale 0 \n set palette rgbformula -7,2,-7 \n set cbrange [0:5] \n set cblabel \"Score\" \n unset cbtics \n set xrange [0: " << pattern.size()-1<< "] \n set yrange[0:" << tempalign.size()<< "] \n set view map \n"<< std::endl;
    myfile<< "p \'-\' using 1:2:3 with image" <<std::endl;
    for (Size i=0;i <pattern.size();++i)
    {
        for (Size j=0; j<tempalign.size();++j)
        {
        debugbuffermatrix[i][j]=0;
        }
    }
    for (Size i =0; i< debugmatrix_.size();++i)
    {
        float score = 0;
        if(debugmatrix_[i][4] == 0)
        {
        score = 5;
        }
        else if(debugmatrix_[i][4] == 1)
        {
        score = 3;
        }
        else if(debugmatrix_[i][4] == 2)
        {
        score = 1;
        }
        debugbuffermatrix[debugmatrix_[i][0]][debugmatrix_[i][1]]=score;
    }
    for (Size i=0; i< debugbuffermatrix.size();++i)
    {
        if(i!=0) myfile << std::endl;
        for (Size j=0; j < debugbuffermatrix[0].size();++j)
        {
        myfile<< i << " " << j << " " << debugbuffermatrix[i][j] << std::endl;
        }
    }
    myfile <<"e"<< std::endl;
    std::cout << " complete tracematrix " << std::endl;
    myfile << std::endl;

    myfile <<"unset multiplot"<<std::endl;

    for (Size i=0; i<debugmatrix_.size();++i)
    {
        myfile<<debugmatrix_[i][0] << " " << debugmatrix_[i][1] << " " << debugmatrix_[i][4] << std::endl;
    }
    myfile.close();
    */
    debugmatrix_.clear();
    debugtraceback_.clear();
    debugscorematrix_.clear();
  }

  void MapAlignmentAlgorithmSpectrumAlignment::debugscoreDistributionCalculation_(float score)
  {
    Int index = (Int)(score + 0.5);
    scoredistribution_.push_back(index);
  }

  void MapAlignmentAlgorithmSpectrumAlignment::updateMembers_()
  {
    gap_    = (float)param_.getValue("gapcost");
    e_      = (float)param_.getValue("affinegapcost");

    // create spectrum compare functor if it does not yet exist
    if (c1_ == nullptr || c1_->getName() != (String)param_.getValue("scorefunction"))
    {
      c1_ = Factory<PeakSpectrumCompareFunctor>::create((String)param_.getValue("scorefunction"));
    }

    cutoffScore_ = (float)param_.getValue("cutoff_score");
    bucketsize_ = (Int)param_.getValue("bucketsize");
    mismatchscore_ = (float)param_.getValue("mismatchscore");
    anchorPoints_ = (Int)param_.getValue("anchorpoints");

    // this is a percentage, should always be between 0 and 100
    if (anchorPoints_ > 100)
    {
      anchorPoints_ = 100;
    }

    String tmp = (String)param_.getValue("debug");
    if (tmp == "true")
    {
      debug_ = true;
    }
    else
    {
      debug_ = false;
    }
    threshold_ = 1 - cutoffScore_;
  }

}
