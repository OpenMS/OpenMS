// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------
//

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>

namespace OpenMS
{
  PeakPickerCWT::PeakPickerCWT()
      : PeakPicker(),
      radius_(3),
      scale_(0.15),
      peak_corr_bound_(0),
      noise_level_(10),
      optimization_(false)
  {

    wt_= new ContinuousWaveletTransformNumIntegration<1>();
    double precision;
    DataValue dv = param_.getValue("Thresholds:Precision");
    if (dv.isEmpty() || dv.toString() == "") precision = 1e-5;
    else precision = (double)dv;

    /** Initialize the Wavelet Transform **/
    double wavelet_spacing;
    dv = param_.getValue("WaveletTransform:Spacing");
    if (dv.isEmpty() || dv.toString() == "") wavelet_spacing= 0.001;
    else wavelet_spacing = (double)dv;

    // estimate the peak bound in the wavelet transform concerning the peak bound in the original signal
    calculatePeakBoundCWT_();

  }

  PeakPickerCWT::PeakPickerCWT(const String& filename) : PeakPicker(filename)
  {
    init();
  }

  PeakPickerCWT::PeakPickerCWT(const Param& parameters) : PeakPicker(parameters)
  {
    init();
  }

  void PeakPickerCWT::init()
  {
    //std::cout << param_ << std::endl;
    // if a peak picking parameter is missed in the param object the value should be substituted by a default value
    DataValue dv =  param_.getValue("Thresholds:Correlation");
    if (dv.isEmpty() || dv.toString() == "") peak_corr_bound_ = 0;
    else peak_corr_bound_ = (float)dv;

    dv = (param_.getValue("Optimization:SkipOptimization"));
    if (dv.isEmpty() || dv.toString() == "") optimization_ = false;
    else optimization_ = (dv.toString() == "no");

    dv = param_.getValue("WaveletTransform:Scale");
    if (dv.isEmpty() || dv.toString() == "") scale_ = 0.15;
    else scale_ = (float)dv;

    dv = param_.getValue("Thresholds:NoiseLevel");
    if (dv.isEmpty() || dv.toString() == "") noise_level_ = 10.;
    else noise_level_ = (float)dv;

    //    std::cout << "Noise Level " << noise_level_ << " scale " << scale_ << std::endl;

    dv =param_.getValue("Thresholds:SearchRadius");
    if (dv.isEmpty() || dv.toString() == "") radius_ = 3;
    else radius_ = (int)dv;

    wt_= new ContinuousWaveletTransformNumIntegration<1>();

    // estimate the peak bound in the wavelet transform concerning the peak bound in the original signal
    calculatePeakBoundCWT_();
  }

  PeakPickerCWT::~PeakPickerCWT()
  {
    if (wt_)
    {
      delete wt_;
#ifdef DEBUG_PEAK_PICKING
      std::cout << "delete wt" << std::endl;
#endif
      wt_=0;
    }
  }

  bool PeakPickerCWT::getMaxPosition_
  ( RawDataPointIterator first,
    RawDataPointIterator last,
    ContinuousWaveletTransform<1>* wt,
    PeakArea_& area,
    int distance_from_scan_border,
    int ms_level,
    int direction)
  {
    // ATTENTION! It is assumed that the resolution==1 (no resolution higher than 1).
    double noise_level=0.;
    double noise_level_cwt=0.;
    if (ms_level==1)
    {
      noise_level = peak_bound_;
      noise_level_cwt = peak_bound_cwt_;
    }
    else
    {
      noise_level = peak_bound_ms2_level_;
      noise_level_cwt = peak_bound_ms2_level_cwt_;
    }

    int zeros_left_index  = wt->getLeftPaddingIndex();
    int zeros_right_index = wt->getRightPaddingIndex();

    // Points to most intensive data point in the signal
    RawDataPointIterator it_max_pos;
    double max_value;

    // Given direction, start the search from left or right
    int start = (direction > 0) ? ((zeros_left_index + 2) + distance_from_scan_border) : ((zeros_right_index - 2) - distance_from_scan_border) ;
    int end   = (direction > 0) ? (zeros_right_index - 1)  : zeros_left_index+1;

    int i=0, j=0, k, max_pos;
    for(i=start, k=0; i!=end; i+=direction, ++k)
    {
      // Check for maximum in cwt at position i
      if((((*wt)[i-1] - (*wt)[i]  ) < 0)
          && (((*wt)[i]   - (*wt)[i+1]) > 0)
          && ( (*wt)[i]   >  noise_level_cwt))
      {
        max_pos = (direction > 0) ? (i - distance_from_scan_border)  : i;
#ifdef DEBUG_PEAK_PICKING
        std::cout << "MAX in CWT at " << (first + max_pos)->getPos()<< " with " << (*wt)[i]
        << std::endl;
#endif
        max_value=(first + max_pos)->getIntensity();


        // search for the corresponding maximum in the signal (consider the radius left and right adjacent points)
        int start_intervall = ((max_pos - (int)radius_) < 0 ) ? 0 : (max_pos - (int)radius_);
        int end_intervall= ((max_pos + (int)radius_) >= distance(first,last)) ? 0 : (max_pos + (int)radius_);

        for(j = start_intervall; j <= end_intervall; ++j)
        {
          if((first + j)->getIntensity() > max_value)
          {
            max_pos = j;
            max_value = (first + j)->getIntensity();
          }

        }

        // if the maximum position is high enough and isn't one of the border points, we return it
        if(((first + max_pos)->getIntensity() >= noise_level)
            && (((first + max_pos) != first)
                && (first + max_pos)   != (last-1)))
        {
          area.max = first + max_pos;
#ifdef DEBUG_PEAK_PICKING
          std::cout << "_________Max in data at__________ " << area.max->getPos()
          << std::endl;
#endif
          return true;
        }

#ifdef DEBUG_PEAK_PICKING
        std::cout << "NO REAL MAX!!!!!!!" << std::endl;
#endif

      }
    }

    // No relevant peak was found
    return false;
  }



  bool PeakPickerCWT::getPeakEndPoints_(RawDataPointIterator first,
                                        RawDataPointIterator last,
                                        PeakArea_& area,
                                        int& peak_left_index,
                                        int& peak_right_index)
  {
    // the Maximum may neither be the first or last point in the signal
    if ((area.max <= first) || (area.max >= last-1))
    {
      return false;
    }

    /*  The algorithm does the following:
     *    - let x_m be the position of the maximum in the data and let (x_l, x_r) be
     *      the left and right neighbours
     *
     *
     *  (1) starting from x_l', walk left until one of the following happens
     *         - the new point is lower than the original bound
     *              => we found our left endpoint
     *
     *         - the new point is larger than the last, but the point left from
     *           the new point is smaller. In that case, we either ran into another
     *           peak, or we encounter some noise. Therefore we now look in the cwt
     *           at the position corresponding to this value. If the cwt here is
     *           monotonous, we consider the point as noise and continue further to the
     *           left. Otherwise, we probably found the beginning of a new peak and
     *           therefore stop here.
     *
     *  (2) analogous procedure to the right of x_r
     */

    RawDataPointIterator it_help=area.max-1;
    PositionType vec_pos;
    int cwt_pos;
    int ep_radius=2;
    int start;
    int stop;
    bool monoton;

    int zeros_left_index  = wt_->getLeftPaddingIndex();

    // search for the left endpoint
    while (((it_help-1) > first) && (it_help->getIntensity() > noise_level_))
    {
      // if the values are still falling to the left, everything is ok.
      if ((it_help-1)->getIntensity() < it_help->getIntensity())
      {
        --it_help;
      }
      // if the values are _rising_, we have to check the cwt
      else
      {
        if ((it_help-2) <= first)
        {
          break;
        }
        // now check the value to the left of the problematic value
        if ((it_help-2)->getIntensity() > (it_help-1)->getIntensity()) // we probably ran into another peak
        {
          break;
        }


        // to the left, the values are falling again => let the cwt decide if we
        // are seeing a new peak or just noise

        // compute the position of the corresponding point in the cwt
        cwt_pos = distance(first, it_help);
        vec_pos = it_help->getPos();

        // since the cwt is pretty smooth usually, we consider the point as noise
        // if the cwt is monotonous in this region
        // TODO: better monotonicity test... say two or three points more
        monoton=true;
        start   =   ((cwt_pos-ep_radius) < 0)
                    ? zeros_left_index+1
                    : cwt_pos-ep_radius+zeros_left_index+1;
        stop    =   ((cwt_pos+ep_radius) > wt_->getSignalLength())
                    ?  (wt_->getSignalLength() + zeros_left_index)
                    : (cwt_pos+ep_radius+zeros_left_index);

        for (; start < stop; ++start)
        {
          if (   ((*wt_)[start-1] - (*wt_)[start]  )
                 * ((*wt_)[start]   - (*wt_)[start+1]) < 0 )
          {
            // different slopes at the sides => stop here
            monoton=false;
            break;
          }
        }

        if (!monoton)
        {
          break;
        }
        --it_help;
      }
    }
    area.left=it_help;

    it_help=area.max+1;
    // search for the right endpoint ???
    while (((it_help+1) < last) && (it_help->getIntensity() > noise_level_))
    {
      // if the values are still falling to the right, everything is ok.
      if (it_help->getIntensity() > (it_help+1)->getIntensity())
      {
        ++it_help;
      }
      // if the values are _rising_, we have to check the cwt
      else
      {
        if ((it_help+2) >= last)
        {
          break;
        }
        // now check the value to the right of the problematic value
        if ((it_help+2)->getIntensity() > (it_help+1)->getIntensity()) // we probably ran into another peak
        {
          break;
        }

        // to the left, the values are falling again => let the cwt decide if we
        // are seeing a new peak or just noise
        // compute the position of the corresponding point in the cwt
        cwt_pos = distance(first, it_help);
        //cwt_pos = distance(first, it_help);
        vec_pos=it_help->getPos();

        // since the cwt is pretty smooth usually, we consider the point as noise
        // if the cwt is monotonous in this region
        // TODO: better monotonicity test... say two or three points more
        monoton = true;

        start   =   ((cwt_pos-ep_radius) < 0)
                    ? zeros_left_index+1
                    : cwt_pos-ep_radius+zeros_left_index+1;
        stop    =   ((cwt_pos+ep_radius) > wt_->getSignalLength())
                    ?  (wt_->getSignalLength() + zeros_left_index)
                    : (cwt_pos+ep_radius+zeros_left_index);

        for (; start < stop; ++start)
        {
          if (   ((*wt_)[start-1] - (*wt_)[start])
                 * ((*wt_)[start]   - (*wt_)[start+1]) < 0 )
          {
            // different slopes at the sides => stop here
            monoton=false;
            break;
          }
        }

        if (!monoton)
        {
          break;
        }
        ++it_help;
      }
    }
    area.right=it_help;

    peak_left_index=distance(first,area.left);
    peak_right_index=distance(first,area.right);

    // The minimal raw data points per peak should be 2
    if ((distance(area.left,area.max) > 0) && (distance(area.max,area.right) > 0))
    {
      return true;
    }

    return false;
  }

  void PeakPickerCWT::getPeakCentroid_(PeakArea_& area)
  {
    RawDataPointIterator left_it=area.max-1, right_it=area.max;
    double max_intensity=area.max->getIntensity();
    double rel_peak_height=max_intensity*0.6;
    double sum=0., w=0.;
    area.centroid_position=area.max->getPos();

    // compute the centroid position (use weighted mean)
    while ((left_it >= area.left) && (left_it->getIntensity() >=rel_peak_height) )
    {
      if (left_it->getIntensity() >=rel_peak_height)
      {
        w+=left_it->getIntensity()*left_it->getPos();
        sum+=left_it->getIntensity();
        --left_it;
      }
    }

    while ((right_it < area.right) && (right_it->getIntensity() >=rel_peak_height) )
    {
      if (right_it->getIntensity() >=rel_peak_height)
      {
        w+=right_it->getIntensity()*right_it->getPos();
        sum+=right_it->getIntensity();
        ++right_it;
      }
    }

    area.centroid_position = w / sum;

#ifdef DEBUG_PEAK_PICKING
    std::cout << "________Centroid is___________ " << area.centroid_position << std::endl;
#endif

  }


  double PeakPickerCWT::lorentz_(double height, double lambda, double pos, double x)
  {
    return height/(1+pow(lambda*(x-pos),2));
  }

  void PeakPickerCWT::calculatePeakBoundCWT_()
  {
#ifdef DEBUG_PEAK_PICKING
    std::cout << "PeakPickerCWT<D>::calculatePeakBoundCWT_ peak_bound_" << peak_bound_ <<  std::endl;
#endif

    // build a lorentz peak of height peak_bound_
    // compute its cwt, and compute the resulting height
    // of the transformed peak

    //compute the peak in the intervall [-2*scale,2*scale]
    double spacing=0.001;
    int n = (int)((4*scale_)/spacing)+1;

    // compute the width parameter using height=peak_bound_ and the peak endpoints should be -scale and +scale, so at
    // positions -scale and +scale the peak value should correspond to the noise_level_
    double lambda = sqrt((-noise_level_*(-peak_bound_+noise_level_)))/(noise_level_*scale_);

    RawDataArrayType lorentz_peak(n);
    RawDataArrayType lorentz_peak2(n);

    /** TODO: switch the type of the transform **/

    ContinuousWaveletTransform<1>* lorentz_cwt;
    ContinuousWaveletTransform<1>* lorentz_ms2_cwt;

    lorentz_cwt = new ContinuousWaveletTransformNumIntegration<1>();
    lorentz_ms2_cwt = new ContinuousWaveletTransformNumIntegration<1>();

    lorentz_cwt->init(scale_, spacing, 0);
    lorentz_ms2_cwt->init(scale_, spacing, 0);
    double start = -2*scale_;
    for (int i=0; i < n; ++i)
    {
      DPosition<1> p;
      p = i*spacing + start;
      lorentz_peak[i].setPosition(p);
      lorentz_peak[i].setIntensity(lorentz_(peak_bound_,lambda,0,i*spacing + start));
      lorentz_peak2[i].setPosition(p);
      lorentz_peak2[i].setIntensity(lorentz_(peak_bound_ms2_level_,lambda,0,i*spacing + start));
    }

    float resolution = 1.;
    lorentz_cwt->transform(lorentz_peak.begin(), lorentz_peak.end(),resolution);
    lorentz_ms2_cwt->transform(lorentz_peak2.begin(), lorentz_peak2.end(),resolution);

    float peak_max=0;
    float peak_max2=0;

    for (int i=0; i<lorentz_cwt->getSignalLength(); i++)
    {
      if ((*lorentz_cwt)[i] > peak_max)
      {
        peak_max = (*lorentz_cwt)[i];
      }
      if ((*lorentz_ms2_cwt)[i] > peak_max2)
      {
        peak_max2 = (*lorentz_ms2_cwt)[i];
      }
    }

    peak_bound_cwt_ = peak_max;
    peak_bound_ms2_level_cwt_ = peak_max2;
#ifdef DEBUG_PEAK_PICKING
    std::cout << "PEAK BOUND IN CWT " << peak_bound_cwt_ << std::endl;
    std::cout << "PEAK BOUND IN CWT (MS 2 Level)" << peak_bound_ms2_level_cwt_ << std::endl;
#endif

    delete lorentz_cwt;
    delete lorentz_ms2_cwt;
  }

  PeakPickerCWT::RawDataPointIterator
  PeakPickerCWT::getIteratorLeftDataPoint_(RawDataPointIterator first,
      RawDataPointIterator last,
      double value)
  {
    int length = distance(first,last);
    double origin = first->getPos();

    OPENMS_PRECONDITION(((origin < value) && (value < (last-1)->getPos())),
                        "The position can't be found in this peak array.");

    double spacing = ((last-1)->getPos()-origin)/(length-1);
    double distance=value-origin;

    int value_index=(int)(distance/spacing);

    RawDataPointIterator it_pos=first+value_index;
    while (true)
    {
      if (it_pos->getPos() < value)
      {
        if ((it_pos+1)->getPos() < value)
        {
          ++it_pos;
        }
        else
        {
          return it_pos;
        }
      }
      else
      {
        --it_pos;
      }
    }
  }

  void PeakPickerCWT::getPeakArea_(const PeakPickerCWT::PeakArea_& area, double& area_left, double& area_right)
  {
    area_left += area.left->getIntensity() * ((area.left+1)->getPos() - area.left->getPos()) * 0.5;
    area_left += area.max->getIntensity() *  (area.max->getPos() - (area.max-1)->getPos()) * 0.5;

    for (RawDataPointIterator pi=area.left+1; pi<area.max; pi++)
    {
      double step = ((pi)->getPos() - (pi-1)->getPos());
      area_left += step * pi->getIntensity();
    }

    area_right += area.right->getIntensity() * ((area.right)->getPos() - (area.right-1)->getPos()) * 0.5;
    area_right += (area.max+1)->getIntensity() *  ((area.max+2)->getPos() - (area.max+1)->getPos()) * 0.5;

    for (RawDataPointIterator pi=area.max+2; pi<area.right; pi++)
    {
      double step = ((pi)->getPos() - (pi-1)->getPos());
      area_right += step * pi->getIntensity();
    }
  }

  PeakShape PeakPickerCWT::fitPeakShape_
  (const PeakPickerCWT::PeakArea_& area,
   bool enable_centroid_fit)
  {

#ifdef DEBUG_PEAK_PICKING
    std::cout << "Left end point: " << area.left->getPos() << std::endl;
#endif
    double max_intensity   =   area.max->getIntensity();
    double left_intensity  =  area.left->getIntensity();
    double right_intensity = area.right->getIntensity();

    //avoid zero width
    float minimal_endpoint_centroid_distance=0.01;
    if (  (fabs( area.left->getPos()-area.centroid_position[0]) < minimal_endpoint_centroid_distance)
          ||(fabs(area.right->getPos()-area.centroid_position[0]) < minimal_endpoint_centroid_distance) )
    {
#ifdef DEBUG_PEAK_PICKING
      std::cout << "The distance between centroid and the endpoints is too small!" << std::endl;
#endif
      return PeakShape();
    }

    if (enable_centroid_fit)
    {
#ifdef DEBUG_PEAK_PICKING
      std::cout << "Fit at the peak centroid" << std::endl;
#endif
      // the maximal position was taken directly from the cwt.
      // first we do a "regular" fit of the left half
      // TODO: avoid zero width!

#ifdef DEBUG_PEAK_PICKING
      std::cout << "Left end point: "         << area.left->getPos()
      << " centroid: "              << area.centroid_position
      << " right end point: "       << area.right->getPos()
      << std::endl;
      std::cout << " point left of centroid:" << (area.left_behind_centroid)->getPos()
      << std::endl;
#endif

      // lorentzian fit

      // estimate the width parameter of the left peak side
      RawDataPointIterator left_it=area.left_behind_centroid;
      double x0=area.centroid_position[0];
      double l_sqrd=0.;
      int n=0;
      while(left_it-1 >= area.left)
      {
        double x1=left_it->getPos();
        double x2=(left_it-1)->getPos();
        double c=(left_it-1)->getIntensity()/left_it->getIntensity();
        l_sqrd+=(1-c)/(c*(pow((x2-x0),2))-pow((x1-x0),2));
        --left_it;
        ++n;
      }
      double left_heigth=area.left_behind_centroid->getIntensity()/(1+l_sqrd*pow(area.left_behind_centroid->getPos()-area.centroid_position[0],2));

      // estimate the width parameter of the right peak side
      RawDataPointIterator right_it=area.left_behind_centroid+1;
      l_sqrd=0.;
      n=0;
      while(right_it+1 <= area.right)
      {
        double x1=right_it->getPos();
        double x2=(right_it+1)->getPos();
        double c=(right_it+1)->getIntensity()/right_it->getIntensity();
        l_sqrd+=(1-c)/(c*(pow((x1-x0),2))-pow((x2-x0),2));
        ++right_it;
        ++n;
      }

      //estimate the heigth
      double right_heigth=(area.left_behind_centroid+1)->getIntensity()/(1+l_sqrd*pow((area.left_behind_centroid+1)->getPos()-area.centroid_position[0],2));

      double height=std::min(left_heigth,right_heigth);

      // compute the left and right area
      double peak_area_left = 0.;
      peak_area_left += area.left->getIntensity() * (  (area.left+1)->getPos()
                        -    area.left->getPos()  ) * 0.5;
      peak_area_left += height * (area.centroid_position[0]-area.left_behind_centroid->getPos()) * 0.5;

      for (RawDataPointIterator pi=area.left+1; pi <= area.left_behind_centroid; pi++)
      {
        double step = ((pi)->getPos() - (pi-1)->getPos());
        peak_area_left += step * pi->getIntensity();
      }

      double peak_area_right = 0.;
      peak_area_right += area.right->getIntensity() * ((area.right)->getPos()
                         - (area.right-1)->getPos()  ) * 0.5;
      peak_area_right += height * ( (area.left_behind_centroid+1)->getPos()-area.centroid_position[0]) * 0.5;

      for (RawDataPointIterator pi=area.left_behind_centroid+1; pi < area.right; pi++)
      {
        double step = ((pi)->getPos() - (pi-1)->getPos());
        peak_area_right += step * pi->getIntensity();
      }

      double left_width =    height/peak_area_left
                             * atan( sqrt( height/area.left->getIntensity() - 1. ) );
      double right_width =  height/peak_area_right
                            * atan( sqrt( height/area.right->getIntensity() - 1. ) );


      // TODO: test different heights; recompute widths; compute area
      PeakShape lorentz(height, area.centroid_position[0], left_width, right_width,
                        peak_area_left + peak_area_right, PeakShapeType::LORENTZ_PEAK);

      lorentz.r_value = correlate_(lorentz, area);

#ifdef DEBUG_PEAK_PICKING
      std::cout << "lorentz r: " << lorentz.r_value << " " << "pos: " << lorentz.mz_position << std::endl;
      std::cout << "w1, w2: " << lorentz.left_width << " " << lorentz.right_width << std::endl;
      std::cout << "h: " << lorentz.height << std::endl;
#endif
      return lorentz;
    }

    else // no fitting on centroids
    {
#ifdef DEBUG_PEAK_PICKING
      std::cout << "fit at the peak maximum " << std::endl;
#endif
      // determine the left half of the area of the PeakArea_...
      double peak_area_left = 0.;
      peak_area_left += area.left->getIntensity() * ((area.left+1)->getPos() - area.left->getPos()) * 0.5;
      peak_area_left += area.max->getIntensity() *  (area.max->getPos() - (area.max-1)->getPos()) * 0.5;

      for (RawDataPointIterator pi=area.left+1; pi<area.max; pi++)
      {
        double step = ((pi)->getPos() - (pi-1)->getPos());
        peak_area_left += step * pi->getIntensity();
      }

      double peak_area_right = 0.;
      peak_area_right += area.right->getIntensity() * ((area.right)->getPos() - (area.right-1)->getPos()) * 0.5;
      peak_area_right += area.max->getIntensity() *  ((area.max+1)->getPos() - (area.max)->getPos()) * 0.5;

      for (RawDataPointIterator pi=area.max+1; pi<area.right; pi++)
      {
        double step = ((pi)->getPos() - (pi-1)->getPos());
        peak_area_right += step * pi->getIntensity();
      }

      // first the lorentz-peak...

      double left_width = max_intensity / peak_area_left * atan(sqrt(max_intensity / left_intensity - 1.));
      double right_width = max_intensity / peak_area_right * atan(sqrt(max_intensity / right_intensity - 1.));



      PeakShape lorentz(max_intensity, area.max->getPos(),
                        left_width, right_width, peak_area_left + peak_area_right,
                        PeakShapeType::LORENTZ_PEAK);

      lorentz.r_value = correlate_(lorentz, area);

      // now the sech-peak...
      left_width  = max_intensity /peak_area_left * sqrt(1. - left_intensity / max_intensity);
      right_width  = max_intensity /peak_area_right * sqrt(1. - right_intensity / max_intensity);


      PeakShape sech(max_intensity, area.max->getPos(),
                     left_width, right_width,
                     peak_area_left + peak_area_right,
                     PeakShapeType::SECH_PEAK);

      sech.r_value = correlate_(sech, area);

#ifdef DEBUG_PEAK_PICKING
      std::cout << "r: " << lorentz.r_value << " " << sech.r_value << std::endl;
      std::cout << "pos: " << lorentz.mz_position <<  " " << sech.mz_position << std::endl;
      std::cout << "w1, w2: " << lorentz.left_width << " " << lorentz.right_width << " "
      << sech.left_width << " " << sech.right_width << std::endl;
      std::cout << "h: " << lorentz.height << std::endl;
#endif

      if ((lorentz.r_value > sech.r_value) && isnan(sech.r_value))
      {
        return lorentz;
      }
      else
      {
        return sech;
      }
    }
  }

  double PeakPickerCWT::correlate_(const PeakShape& peak,
                                   const PeakPickerCWT::PeakArea_& area,
                                   int direction) const
  {
    double SSxx = 0., SSyy = 0., SSxy = 0.;

    // compute the averages
    double data_average=0., fit_average=0.;
    double data_sqr=0., fit_sqr=0.;
    double cross=0.;

    int number_of_points = 0;
    RawDataPointIterator corr_begin=area.left;
    RawDataPointIterator corr_end=area.right;

    // for separate overlapping peak correlate until the max position...
    if (direction > 0)
      corr_end=area.max;
    else
      if (direction < 0)
        corr_begin=area.max;

    for (RawDataPointIterator pi = corr_begin; pi<=corr_end; pi++)
    {
      double data_val = pi->getIntensity();
      double peak_val = peak(pi->getPos());

      data_average += data_val;
      fit_average  += peak_val;

      data_sqr += data_val * data_val;
      fit_sqr  += peak_val * peak_val;

      cross += data_val * peak_val;

      number_of_points++;
    }

    if (number_of_points == 0)
      return 0.;

    data_average /= number_of_points;
    fit_average  /= number_of_points;

    SSxx = data_sqr - number_of_points * (data_average * data_average);
    SSyy = fit_sqr - number_of_points * (fit_average * fit_average);
    SSxy = cross - number_of_points * (data_average * fit_average);

    return (SSxy * SSxy) / (SSxx * SSyy);
  }
  
  template <>
  void PeakPickerCWT::fillPeak_< DPickedPeak<1> >(const PeakShape& peak_shape, DPickedPeak<1>& picked_peak)
  {
    picked_peak.getRValue() = peak_shape.r_value;
    picked_peak.getArea() = peak_shape.area;
    picked_peak.getFWHM() = peak_shape.getFWHM();
    picked_peak.getLeftWidthParameter() = peak_shape.left_width;
    picked_peak.getRightWidthParameter() = peak_shape.right_width;
    picked_peak.getPeakShape() = peak_shape.type;
    picked_peak.getSN() = peak_shape.signal_to_noise;
  }
}
