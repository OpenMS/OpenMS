// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
#include <cmath>

using namespace std;

namespace OpenMS
{
  PeakPickerCWT::PeakPickerCWT()
		: PeakPicker(),
      	
      radius_(0),
      scale_(0.0),
      peak_bound_cwt_(0.0),
      peak_bound_ms2_level_cwt_(0.0),
      peak_corr_bound_(0.0),
      noise_level_(0.0),
      optimization_(false)
  {
    // if a peak picking parameter is missed in the param object the value should be substituted by a default value
  	defaults_.setValue("thresholds:correlation",0.5);
  	defaults_.setValue("wavelet_transform:scale",0.15);
  	defaults_.setValue("wavelet_transform:spacing",0.001);
  	defaults_.setValue("thresholds:noise_level",0.1);
   	defaults_.setValue("thresholds:search_radius",3); 	
		defaults_.setValue("thresholds:signal_to_noise",3);
		
		//Optimization parameters
  	defaults_.setValue("Optimization:optimization","no");
		defaults_.setValue("Optimization:penalties:position",0.0);
		defaults_.setValue("Optimization:penalties:left_width",1.0); 	
		defaults_.setValue("Optimization:penalties:right_width",1.0); 	
		defaults_.setValue("Optimization:iterations",15); 	
		defaults_.setValue("Optimization:delta_abs_error",1e-04f); 	
		defaults_.setValue("Optimization:delta_rel_error",1e-04f);

		// deconvolution parameters
		defaults_.setValue("deconvolution:penalties:position",0.0);
		defaults_.setValue("deconvolution:penalties:height",1.0);
		defaults_.setValue("deconvolution:penalties:left_width",0.0);
		defaults_.setValue("deconvolution:penalties:right_width",0.0);
		defaults_.setValue("deconvolution:skip_deconvolution","yes");
		defaults_.setValue("deconvolution:width_threshold",1.1);
		defaults_.setValue("deconvolution:asym_threshold",0.3);
		defaults_.setValue("deconvolution:fwhm_threshold",0.7);
		defaults_.setValue("deconvolution:delta_abs_error",1e-05f);
		defaults_.setValue("deconvolution:delta_rel_error",1e-05f);
		defaults_.setValue("deconvolution:left_width",2);
		defaults_.setValue("deconvolution:right_width",2);
		defaults_.setValue("deconvolution:scaling",0.12);
		defaults_.setValue("deconvolution:iterations",10);

		// 2D optimization parameters
		defaults_.setValue("2D_optimization:penalties:position",0.0);
		defaults_.setValue("2D_optimization:penalties:height",1.0);
		defaults_.setValue("2D_optimization:penalties:left_width",0.0);
		defaults_.setValue("2D_optimization:penalties:right_width",0.0);
		defaults_.setValue("2D_optimization:thresholds:tolerance_mz",0.2);
		defaults_.setValue("2D_optimization:thresholds:max_peak_distance",1.0);
		defaults_.setValue("2D_optimization:skip_optimization","yes");
		defaults_.setValue("2D_optimization:delta_abs_error",1e-05f);
		defaults_.setValue("2D_optimization:delta_rel_error",1e-05f);
		defaults_.setValue("2D_optimization:iterations",10);

		
		
  	
		setParam(Param());
  }

  PeakPickerCWT::~PeakPickerCWT()
  {}

  void PeakPickerCWT::setParam(Param param)
  {
    PeakPicker::setParam(param);

    peak_corr_bound_ = (float)param_.getValue("thresholds:correlation");
    String opt = param_.getValue("Optimization:optimization").toString();
    if (opt=="one_dimensional")
			{
				optimization_ = true;
				two_d_optimization_ = false;
			}
    else if (opt=="two_dimensional")
			{
				two_d_optimization_ = true;
				optimization_ = false;
			}
    else 
			{
				optimization_ = false;
				two_d_optimization_ = false;
    	
			}


    scale_ = (float)param_.getValue("wavelet_transform:scale");
    noise_level_ = (float)param_.getValue("thresholds:noise_level");
    radius_ = (int)param_.getValue("thresholds:search_radius");
		signal_to_noise_ = (float)param_.getValue("thresholds:signal_to_noise");

		opt = param_.getValue("deconvolution:skip_deconvolution").toString();
		if (opt=="yes")
			{
				deconvolution_ = false;
			}
		else if (opt=="no")
			{
				deconvolution_ = true;
			}
		else
			{
				cerr << "Warning: PeakPickerCWT option 'deconvolution:skip_deconvolution' should be 'yes' or 'no'!"
						 << " It is set to '" << opt << "'" << endl;
			}
		
	}
	

  bool PeakPickerCWT::getMaxPosition_
  ( RawDataPointIterator first,
    RawDataPointIterator last,
    const ContinuousWaveletTransform& wt,
    PeakArea_& area,
    int distance_from_scan_border,
    int ms_level,
    int direction)
  {
    // ATTENTION! It is assumed that the resolution==1 (no resolution higher than 1).
    // Comment: Who cares ??
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

    int zeros_left_index  = wt.getLeftPaddingIndex();
    int zeros_right_index = wt.getRightPaddingIndex();

    // Points to most intensive data point in the signal
    RawDataPointIterator it_max_pos;
    double max_value;

    // Given direction, start the search from left or right
    int start = (direction > 0) ? ((zeros_left_index + 2) + distance_from_scan_border) : ((zeros_right_index - 2) - distance_from_scan_border) ;
    int end   = (direction > 0) ? (zeros_right_index - 1)  : zeros_left_index+1;

    int i=0, j=0, k, max_pos;
    for(i=start, k=0; i!=end; i+=direction, ++k)
			{
#ifdef DEBUG_PEAK_PICKING
				std::cout << "Search for max pos in cwt " <<  noise_level_cwt <<std::endl;
#endif
				// Check for maximum in cwt at position i
				if(((wt[i-1] - wt[i]  ) < 0)
					 && ((wt[i] - wt[i+1]) > 0)
					 && ( wt[i]  >  noise_level_cwt))
					{
						max_pos = (direction > 0) ? (i - distance_from_scan_border)  : i;
#ifdef DEBUG_PEAK_PICKING

						std::cout << "MAX in CWT at " << (first + max_pos)->getPos()<< " with " << wt[i]
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
                                        int distance_from_scan_border,
                                        int& peak_left_index,
                                        int& peak_right_index)
  {
    // the Maximum may neither be the first or last point in the signal
    if ((area.max <= first) || (area.max >= last-1))
			{
				return false;
			}

    RawDataPointIterator it_help=area.max-1;
    PositionType vec_pos;
    int cwt_pos;
    int ep_radius=2;
    int start;
    int stop;
    bool monoton;

    int zeros_left_index  = wt_.getLeftPaddingIndex();

    // search for the left endpoint
    while (((it_help-1) > first) && (it_help->getIntensity() > noise_level_))
			{
	
#ifdef DEBUG_PEAK_PICKING
				std::cout << "while left endpoint " << std::endl;
#endif
				// if the values are still falling to the left, everything is ok.
				if ((it_help-1)->getIntensity() < it_help->getIntensity())
					{
						--it_help;

#ifdef DEBUG_PEAK_PICKING
						std::cout << "it_help " << it_help->getPos() << std::endl;
#endif

					}
				// if the values are _rising_, we have to check the cwt
				else
					{
						if ((it_help-2) <= first)
							{
#ifdef DEBUG_PEAK_PICKING        	
								std::cout << "it_help-2) <= first"  << std::endl;
#endif      
								break;
							}
						// now check the value to the left of the problematic value
						if ((it_help-2)->getIntensity() > (it_help-1)->getIntensity()) // we probably ran into another peak
							{
#ifdef DEBUG_PEAK_PICKING        	
								std::cout << "((it_help-2)->getIntensity() > (it_help-1)->getIntensity()"  << std::endl;
#endif        		
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

						start   =   (cwt_pos < ep_radius)
							? (distance_from_scan_border + zeros_left_index) + 2
							: cwt_pos - ep_radius + (distance_from_scan_border + zeros_left_index + 2);
						stop    =   ((cwt_pos + ep_radius) > distance(it_help,last))
							?  (wt_.getSize() - 2)
							: cwt_pos + ep_radius + (distance_from_scan_border + zeros_left_index + 2);

#ifdef DEBUG_PEAK_PICKING
						std::cout << "start " << start << " stop " << stop << std::endl;
#endif					
						for (; start < stop; ++start)
							{
								if (   (wt_[start-1] - wt_[start]  )
											 * (wt_[start]   - wt_[start+1]) < 0 )
									{
										// different slopes at the sides => stop here
#ifdef DEBUG_PEAK_PICKING            
										std::cout << "monoton test " << wt_.getSignal()[start-1].getPos() 
															<< " " <<  wt_.getSignal()[start].getPos()
															<< " " <<  wt_.getSignal()[start+1].getPos() << std::endl;
#endif            						
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
#ifdef DEBUG_PEAK_PICKING    	
				std::cout << "while right endpoint " << std::endl;
#endif    		
				//      if the values are still falling to the right, everything is ok.
				if (it_help->getIntensity() > (it_help+1)->getIntensity())
					{
						++it_help;
#ifdef DEBUG_PEAK_PICKING        
						std::cout << "it_help " << it_help->getPos() << std::endl;
#endif         	
					}
				// if the values are _rising_, we have to check the cwt
				else
					{
						if ((it_help+2) >= last)
							{
#ifdef DEBUG_PEAK_PICKING        	
								std::cout << "it_help+2) <= first"  << std::endl;
#endif        		
								break;
							}
						// now check the value to the right of the problematic value
						if ((it_help+2)->getIntensity() > (it_help+1)->getIntensity()) // we probably ran into another peak
							{
#ifdef DEBUG_PEAK_PICKING        	
								std::cout << "(it_help+2)->getIntensity() > (it_help+1)->getIntensity())"  << std::endl;
#endif        		
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

						start   =   (cwt_pos < ep_radius)
							? (distance_from_scan_border + zeros_left_index) + 2
							: cwt_pos - ep_radius + (distance_from_scan_border + zeros_left_index + 2);
						stop    =   ((cwt_pos + ep_radius) > distance(it_help,last))
							?  (wt_.getSize() - 2)
							: cwt_pos + ep_radius + (distance_from_scan_border + zeros_left_index + 2);

#ifdef DEBUG_PEAK_PICKING
						std::cout << "start " << start << " stop " << stop << std::endl;
#endif					
						for (; start < stop; ++start)
							{
								if (   (wt_[start-1] - wt_[start])
											 * (wt_[start]  - wt_[start+1]) < 0 )
									{
										// different slopes at the sides => stop here
#ifdef DEBUG_PEAK_PICKING            
										std::cout << "monoton test " << wt_.getSignal()[start-1].getPos() 
															<< " " <<  wt_.getSignal()[start].getPos()
															<< " " <<  wt_.getSignal()[start+1].getPos() << std::endl;
#endif            						
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
#ifdef DEBUG_PEAK_PICKING
				std::cout << "OK !!!! The endpoints are "
									<< area.left->getPosition()
									<< " and "
									<< area.right->getPosition()
									<< std::endl;
#endif

				return true;
			}
#ifdef DEBUG_PEAK_PICKING
    std::cout << "NO!!! The endpoints are "
							<< area.left->getPosition()
							<< " and "
							<< area.right->getPosition()
							<< std::endl;
#endif

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

  void PeakPickerCWT::initializeWT_()
  {
#ifdef DEBUG_PEAK_PICKING
    std::cout << "PeakPickerCWT<D>::initialize_ peak_bound_" << peak_bound_ <<  std::endl;
#endif
    //initialize wavelet transformer
    wt_.init(scale_, (double)param_.getValue("wavelet_transform:spacing"));

    //calculate peak bound in CWT

    // build a lorentz peak of height peak_bound_
    // compute its cwt, and compute the resulting height
    // of the transformed peak

    //compute the peak in the intervall [-2*scale,2*scale]
    double spacing=0.001;
    int n = (int)((4*scale_)/spacing)+1;

    double lambda = 2. / scale_;
    // compute the width parameter using height=peak_bound_ and the peak endpoints should be -scale and +scale, so at
    // positions -scale and +scale the peak value should correspond to the noise_level_
    //double lambda = sqrt((-noise_level_*(-peak_bound_+noise_level_)))/(noise_level_*scale_);

    RawDataArrayType lorentz_peak(n);
    RawDataArrayType lorentz_peak2(n);

    // TODO: switch the type of the transform

    ContinuousWaveletTransformNumIntegration lorentz_cwt;
    ContinuousWaveletTransformNumIntegration lorentz_ms2_cwt;

    lorentz_cwt.init(scale_, spacing);
    lorentz_ms2_cwt.init(scale_, spacing);
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
    lorentz_cwt.transform(lorentz_peak.begin(), lorentz_peak.end(), resolution);
    lorentz_ms2_cwt.transform(lorentz_peak2.begin(), lorentz_peak2.end(), resolution);

    float peak_max=0;
    float peak_max2=0;

    for (int i=0; i < lorentz_cwt.getSignalLength(); i++)
			{
				if (lorentz_cwt[i] > peak_max)
					{
						peak_max = lorentz_cwt[i];
					}
				if (lorentz_ms2_cwt[i] > peak_max2)
					{
						peak_max2 = lorentz_ms2_cwt[i];
					}
			}

    peak_bound_cwt_ = peak_max;
    peak_bound_ms2_level_cwt_ = peak_max2;
#ifdef DEBUG_PEAK_PICKING

    std::cout << "PEAK BOUND IN CWT " << peak_bound_cwt_ << std::endl;
    std::cout << "PEAK BOUND IN CWT (MS 2 Level)" << peak_bound_ms2_level_cwt_ << std::endl;
#endif

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

    if (enable_centroid_fit)
			{
#ifdef DEBUG_PEAK_PICKING
				std::cout << "Fit at the peak centroid" << std::endl;
#endif

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

				if ((lorentz.r_value > sech.r_value) && std::isnan(sech.r_value))
					{
						return lorentz;
					}
				else
					{
						return sech;
					}
			}
  } 

	void PeakPickerCWT::deconvolutePeak_(PeakShape& shape,
																			 PeakPickerCWT::PeakArea_& area,
																			 std::vector<double> peak_endpoints)
	{
		// scaling for charge one
		float scaling_DC = (float) param_.getValue("deconvolution:scaling");

		
		// init and calculate the transform of the signal in the convoluted region
		// first take the scaling for charge 2
		wtDC_.init(scaling_DC/2,  wt_.getSpacing());
		wtDC_.transform(area.left,area.right,2);


#ifdef DEBUG_DECONV
		std::cout << "------------------\n---------------------\nconvoluted area begin "<<area.left->getPos()<<"\tend "<<area.right->getPos()<<std::endl;
#endif
		
		
		int charge=2,old_peaks;
		std::vector<double> peak_values,old_peak_values;
		std::vector<PeakShape> peaks_DC;
		int peaks = getNumberOfPeaks_(area.left,area.right,peak_values,/*charge*/1,2,wtDC_);
		
																											
		bool correct_scale = false;

 		//	if(peaks > 1) correct_scale = true;
		std::vector<double> distances;
		while(!correct_scale)
			{
				old_peaks = peaks;
				old_peak_values = peak_values;
#ifdef DEBUG_DECONV
				std::cout << "current charge: "<< charge << "\t #peaks: " << peaks << std::endl;
				std::cout << "current scale: "<< scaling_DC/charge<<std::endl;
#endif
				correct_scale = true;
				distances.clear();
				//check distances between the peaks, if there is a larger one than try a higher charge
				for(unsigned int i=1; i<peak_values.size()/2;++i)
					{
						distances.push_back(peak_values[2*i+1] - peak_values[2*i-1]);
#ifdef DEBUG_DECONV
						std::cout << "distances["<< distances.size()-1<<"] = "
											<< peak_values[2*i+1] - peak_values[2*i-1]
											<< " ---- "<<peak_values[2*i+1]<< " - " <<peak_values[2*i-1]
											<< " ---- "<<distances[distances.size()-1]  <<std::endl;
#endif
					}
				if(peaks <=1)
					{
						//std::cout << "peaks "<< peaks << "\told_peaks "<< old_peaks << std::endl;
				
						++charge;
						wtDC_.init(scaling_DC/charge,  wt_.getSpacing());
						wtDC_.transform(area.left,area.right,2);
						peak_values.clear();
						peaks = getNumberOfPeaks_(area.left,area.right,peak_values,/*charge,*/1,2,wtDC_);
						//	std::cout << "peaks "<< peaks << "\told_peaks "<< old_peaks << std::endl;
						if(peaks <= old_peaks)
							{
								//	std::cout << "peaks "<< peaks << "\told_peaks "<< old_peaks << std::endl;
								peaks = old_peaks;
								peak_values = old_peak_values;
								correct_scale = true;
							}
						continue;
					}
				sort(distances.begin(),distances.end());
				unsigned int index = static_cast<unsigned int>(floor((distances.size()-1)/2));
				double median = distances[index];
#ifdef DEBUG_DECONV
				std::cout << "median: distance[" << index<<"]= "<<median << std::endl;
#endif
				for(;index<distances.size();++index)
					{
						// if there is one differing distance try next scale
						if(distances[index] - median > 0.2 || distances[index] > 1.003/charge+0.15)
							{
#ifdef DEBUG_DECONV
								std::cout << "distance too large: " << distances[index] - median
													<< " > " <<  0.2 << " or "<<distances[index]<<" > "
													<<1.003/charge+0.15<<std::endl;
#endif

								correct_scale = false;
								break;
							}
					}
				if(!correct_scale)
					{
						++charge;
						wtDC_.init(scaling_DC/charge,  wt_.getSpacing());
						wtDC_.transform(area.left,area.right,2);
						peak_values.clear();
						peaks = getNumberOfPeaks_(area.left,area.right,peak_values,/*charge,*/1,2,wtDC_);
						if(peaks <= old_peaks)
							{
								peaks = old_peaks;
								peak_values = old_peak_values;
								correct_scale = true;
							}
						
					}
			}
#ifdef DEBUG_PEAK_PICKING
		std::cout << "Number of peaks: "<<peaks << std::endl;
#endif
		// determine the probable charge state
		// the best result from [charge-1,charge,charge+1] will be taken
		// if charge equals zero, something bad happened
		charge = determineChargeState_(peak_values);
#ifdef DEBUG_PEAK_PICKING
		std::cout << "charge "<<charge<<std::endl;
#endif
		// one peak needn't be deconvoluted
		if (peaks > 1 && charge >0)
			{
				
				OptimizationFunctions::positions_DC_.clear();
				OptimizationFunctions::signal_DC_.clear();
				OptimizationFunctions::peaks_DC_.clear();

				// enter zero-intensity at the left margin
				OptimizationFunctions::positions_DC_.push_back((area.left)->getPos()-0.2);
				OptimizationFunctions::signal_DC_.push_back(0);	
				
				for (unsigned int i = 0; area.left+i != area.right ;++i)
					{
						OptimizationFunctions::positions_DC_.push_back((area.left+i)->getPos());
						OptimizationFunctions::signal_DC_.push_back((area.left+i)->getIntensity());	
					}
				OptimizationFunctions::positions_DC_.push_back((area.right)->getPos());
				OptimizationFunctions::signal_DC_.push_back((area.right)->getIntensity());	

				OptimizationFunctions::positions_DC_.push_back((area.right)->getPos()+0.2);
				OptimizationFunctions::signal_DC_.push_back(0);	
				
				
			
				// initial parameters for the optimization
				double leftwidth = (float)param_.getValue("deconvolution:left_width");
				double rightwidth = (float)param_.getValue("deconvolution:right_width");
			
				double dist = 1.003 / charge;

				
				peaks_DC.resize(peaks);
				PeakShape peak(peak_values[0],peak_values[1],leftwidth,rightwidth,0,PeakShapeType::SECH_PEAK);
				peaks_DC[0]=peak;
				//	std::cout<<"peak.mz_position "<<peak.mz_position<<std::endl;
				for(int i=1;i<peaks;++i)
					{
						PeakShape peak(peak_values[2*i],peak_values[1]+i*dist,leftwidth,rightwidth,0,PeakShapeType::SECH_PEAK);
						peaks_DC[i]=peak;
						
#ifdef DEBUG_DECONV
						std::cout<<"peak.mz_position "<<peak.mz_position<<std::endl;
#endif
					}
				OpenMS::OptimizationFunctions::PenaltyFactorsInt penalties;
				penalties.height = (float)param_.getValue("deconvolution:penalties:height");
 				penalties.pos = (float)param_.getValue("deconvolution:penalties:position");
				penalties.lWidth = param_.getValue("deconvolution:penalties:left_width");
				penalties.rWidth = param_.getValue("deconvolution:penalties:right_width");

				unsigned int max_iteration = (int)param_.getValue("deconvolution:iterations");

				double eps_abs = (float)param_.getValue("deconvolution:delta_abs_error");
				double eps_rel = (float)param_.getValue("deconvolution:delta_rel_error");

				OptimizePeakDeconvolution opt(penalties,max_iteration,eps_abs,eps_rel,charge);
				
				int runs=0;
				
#ifdef DEBUG_DECONV
				std::cout<<"OptimizationType: Levenberg-Marquardt mit "<<peaks_DC.size()
								 <<"peaks\n";
#endif
						
				// if the optimization fails (peaks are too broad) try entering an additional peak
				while(!opt.optimize(peaks_DC,param_,runs))
					{
						++runs;
#ifdef DEBUG_DECONV
						std::cout << "\n\n\nat least one peak too broad -- adding one peak\n\n"<<std::endl;
#endif
						addPeak_(peaks_DC,area,leftwidth,rightwidth);
					
					}
						
      
						
				
#ifdef DEBUG_DECONV
				
				std::cout<<" \nLM results:\n";
				
				//print gsl results:
				for(int i=0;i<(int)peaks_DC.size();++i)
					{
						std::cout<<"\nposLM("<<i+1<<")="<<peaks_DC[i].mz_position;
						std::cout<<"\nheightLM("<<i+1<<")="<<peaks_DC[i].height;
						if(peaks_DC[i].type == PeakShapeType::LORENTZ_PEAK)
							std::cout<<"\npeak("<<i+1<<")=0\n";
						else
							std::cout<<"\npeak("<<i+1<<")=1\n";
					}
				std::cout<<"\nlw="<<peaks_DC[0].left_width;			
				std::cout<<"\nrw="<<peaks_DC[0].right_width<<"\n\n";			
				
#endif
		    double left_endpoint,right_endpoint;
				for(unsigned int curr_peak=0;curr_peak<peaks_DC.size();++curr_peak)
					{
						peak_shapes_.push_back(peaks_DC[curr_peak]);
		
						PeakShape p = peaks_DC[curr_peak];
						if (peaks_DC[curr_peak].type == PeakShapeType::LORENTZ_PEAK)
							{
								left_endpoint=p.mz_position+1/p.left_width*sqrt(p.height/1-1);
								right_endpoint=p.mz_position+1/p.right_width*sqrt(p.height/1-1);
							}
						else
							{
								left_endpoint=p.mz_position+1/p.left_width*acosh(sqrt(p.height/0.001));
								right_endpoint=p.mz_position+1/p.right_width*acosh(sqrt(p.height/0.001));
							}
						peak_endpoints.push_back(left_endpoint);
						peak_endpoints.push_back(right_endpoint);
					}
				
				OptimizationFunctions::peaks_DC_.clear();
				OptimizationFunctions::signal_DC_.clear();
				OptimizationFunctions::positions_DC_.clear();
				peaks_DC.clear();
			}
		else
			{
				peak_shapes_.push_back(shape);
				peak_endpoints.push_back(area.left->getPos());
				peak_endpoints.push_back(area.right->getPos());
			}

		
	}


	void PeakPickerCWT::addPeak_(std::vector<PeakShape>& peaks_DC,PeakArea_& area,double left_width,double right_width)
	{
		// just enter a peak using equally spaced peak positions

		double peak_width = area.right->getPos() - area.left->getPos();
		int num_peaks = peaks_DC.size()+1;

		double dist = peak_width / (num_peaks+1);

		// put peak into peak vector using default values for the widths and peak type
		peaks_DC.push_back(PeakShape(0,0,left_width,right_width,0,PeakShapeType::SECH_PEAK));

		// adjust the positions and get their initial intensities from the raw data
		for(int i=0; i < num_peaks; ++i)
			{
				peaks_DC[i].mz_position = area.left->getPos() + dist/2 + i*dist;

				std::vector<double>::iterator it_help = lower_bound(OptimizationFunctions::positions_DC_.begin(),
																														OptimizationFunctions::positions_DC_.end(),
																														peaks_DC[i].mz_position);
				if(it_help != OptimizationFunctions::positions_DC_.end())
					{
						peaks_DC[i].height =
							OptimizationFunctions::signal_DC_[distance(OptimizationFunctions::positions_DC_.begin(),it_help)]/10;
#ifdef DEBUG_DECONV
						std::cout << "height "<<i<<"   "<<	peaks_DC[i].height<<"\t"
											<<distance(OptimizationFunctions::positions_DC_.begin(),it_help)<<"\n";
#endif
					}
				else
					{
						peaks_DC[i].height =
							OptimizationFunctions::signal_DC_[OptimizationFunctions::positions_DC_.size()-1];
#ifdef DEBUG_DECONV
						std::cout << "else height "<<i<<"   "<<	peaks_DC[i].height<<"\n";
#endif
					}
			}
		
		
	}
	
	int PeakPickerCWT::getNumberOfPeaks_(RawDataPointIterator& first,
																			 RawDataPointIterator& last,
																			 std::vector<double>& peak_values,
																			 int direction,
																			 int resolution,
																			 ContinuousWaveletTransformNumIntegration& wt)
  {
    double noise_level=0.;
    double noise_level_cwt=0.;
    
		noise_level = peak_bound_;
		noise_level_cwt = peak_bound_cwt_;
    
#ifdef DEBUG_DECONV
    std::cout<<"noise_level = "<<noise_level<<";\tnoise_level_cwt = "<<noise_level_cwt<<";\n";
		std::cout << "resoltuion "<<resolution <<"\n";
		std::cout <<"first und last "<< first->getPos() << "\t"<<last->getPos() << std::endl;
				
#endif

    int found = 0;
    
    int zeros_left_index  = wt.getLeftPaddingIndex();
    int zeros_right_index = wt.getRightPaddingIndex();

    // The maximum intensity in the signal
    RawDataPointIterator it_max_pos;
    //double max_value;T
    int start = (direction>0) ? zeros_left_index+2 : zeros_right_index-2;
    int end   = (direction>0) ? zeros_right_index-1  : zeros_left_index+1;
	
    int i=0, max_pos;
    int k=0;

    std::vector<double>::iterator checker;
		while(wt.getSignal()[start+1].getPos() <= first->getPos())     ++start;
		//k=i;
		int offset = start;
		while(wt.getSignal()[end].getPos() > last->getPos())     --end;
		for(i=start; i!=end; i+=direction,k+=direction) 
			{

				// Does a maximum in the wavelet transform at position i exists?
				if(((wt[i-1] - wt[i]  ) < 0) 
					 && ((wt[i]   - wt[i+1]) > 0) 
					 && ( wt[i]   >  noise_level_cwt))					
					{
#ifdef DEBUG_DECONV
						std::cout << "MAX in CWT at " << (first+(i-offset)/resolution)->getPosition()[0]<< " with " << wt[i]
											<< std::endl;

#endif
						max_pos=(i-offset)/resolution;

						// if the maximum position is high enough and isn't one of the border points, we return it
						if(((first+max_pos)->getIntensity() >= noise_level) 
							 && (((first+max_pos) != first) 
									 && (first+max_pos)   != (last-1))
							 )
							{
#ifdef DEBUG_DECONV
								std::cout << "_________REAL Max in CWT at__________ " << (first+max_pos)->getPosition()[0] << "\n";
#endif
								peak_values.push_back((first+max_pos)->getIntensity());
								peak_values.push_back((first+max_pos)->getPos());
								
								++found;
							}
					
					}
      }
#ifdef DEBUG_DECONV
		std::cout<<"number of peaks found: "<<found<<"\n";
#endif
    return found;
  }

  int PeakPickerCWT::determineChargeState_(std::vector<double>& peak_values)
  {
    int charge;
    int peaks = (int)peak_values.size() / 2;
    if(peaks>1)
      {
				double dif = 0;
				int i=peaks-1;
				while(i>0)
					{
						dif += fabs(peak_values[2*i+1] - peak_values[2*(i-1)+1]);
#ifdef DEBUG_DECONV
						std::cout << "peak_values["<<i <<"] "<< peak_values[2*i+1]
											<< " peak_values[" <<i-1 <<"] "<< peak_values[2*(i-1)+1]
											<<" dif["<<i<<"]=" << dif<<std::endl;
#endif
						--i;
					}
				dif /= peaks-1;
				charge = (int) round(1/dif);
				if(isnan(charge) || isinf(charge)) charge = 0;
#ifdef DEBUG_DECONV
				std::cout<<"1/dif = "<<1/dif<<";\tcharge = "<<charge<<std::endl;
#endif
      }
    else charge = 1;
    
    return charge;
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
