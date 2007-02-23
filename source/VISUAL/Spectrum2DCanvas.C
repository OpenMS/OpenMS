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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/KERNEL/DFeature.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum2DCanvasPDP.h>
#include <OpenMS/CONCEPT/TimeStamp.h>

//STL
#include <algorithm>	

//QT
#include <QtGui/QWheelEvent>
#include <QtGui/QMouseEvent>
#include <QtGui/QPainter>
#include <QtGui/QBitmap>
#include <QtGui/QPolygon>
#include <QtCore/QTime>

using namespace std;

namespace OpenMS
{
	using namespace Internal;

	Spectrum2DCanvas::Spectrum2DCanvas(QWidget* parent)
		: SpectrumCanvas(parent),
		marching_squares_matrices_(),
		max_values_(),
		show_contours_(),
		show_surface_(),
		show_dots_(),
		selected_peak_(0),
		measurement_start_(0),
		measurement_stop_(0),
		tmp_peak_(),
		dot_gradient_(),
		surface_gradient_()
	{
		projection_mz_.resize(1);
		projection_rt_.resize(1);
	}
	
	Spectrum2DCanvas::~Spectrum2DCanvas()
	{
	}
	
	void Spectrum2DCanvas::showContours(bool on)
	{
		if (on != show_contours_[current_layer_])
		{
			show_contours_[current_layer_] = on;
			update_buffer_ = true;
			update();
		}
	}
	
	void Spectrum2DCanvas::showSurface(bool on)
	{
		if (on != show_surface_[current_layer_])
		{
			show_surface_[current_layer_] = on;
			update_buffer_ = true;
			update();
		}
	}
	
	void Spectrum2DCanvas::showPoints(bool on)
	{
		if (on != show_dots_[current_layer_])
		{
			show_dots_[current_layer_] = on;
			update_buffer_ = true;
			update();
		}
	}
	
	bool Spectrum2DCanvas::contoursAreShown()
	{
		return show_contours_[current_layer_];
	}
	
	bool Spectrum2DCanvas::surfaceIsShown()
	{
		return show_surface_[current_layer_];
	}
	
	bool Spectrum2DCanvas::dotsAreShown()
	{
		return show_dots_[current_layer_];
	}
	
	void Spectrum2DCanvas::highlightPeak_(QPainter& painter, DFeature<2>* peak)
	{
		if (!peak) return;
		painter.save();
		painter.setPen(QPen(Qt::red, 2));
		QPoint pos;
		dataToWidget_(peak->getPosition(),pos);
		painter.drawEllipse(pos.x() - 5, pos.y() - 5, 10, 10);
		painter.restore();
	}
	
	DFeature<2>* Spectrum2DCanvas::findNearestPeak_(const QPoint& pos)
	{
		//Constructing the area corrects swapped mapping of RT and m/z
		AreaType area (widgetToData_(pos - QPoint(5,5)),widgetToData_(pos + QPoint(5,5)));

		DFeature<2>* max_peak = 0;
		float max_int = -1 * numeric_limits<float>::max();
		
		//cout << "findNearestPeak_: Int range -- " << getCurrentLayer().min_int << " "  << getCurrentLayer().max_int << endl;
		
		if (getCurrentLayer().type==LayerData::DT_PEAK)
		{
			for (ExperimentType::ConstAreaIterator i = getCurrentPeakData().areaBeginConst(area.min()[1],area.max()[1],area.min()[0],area.max()[0]); 
					 i != getCurrentPeakData().areaEndConst(); 
					 ++i)
			{
				if (i->getIntensity() > max_int && i->getIntensity()>=getCurrentLayer().min_int && i->getIntensity()<=getCurrentLayer().max_int)
				{
					//cout << "new max: " << i.getRetentionTime() << " " << i->getPos() << endl;
					max_int = i->getIntensity();
					
					tmp_peak_.setIntensity(i->getIntensity());
					tmp_peak_.getPosition()[0] = i->getPos();
					tmp_peak_.getPosition()[1] = i.getRetentionTime();
					
					max_peak = &tmp_peak_;
				}
			}
	 	}
	 	else
	 	{
			for (FeatureMapType::ConstIterator i = getCurrentLayer().features.begin();
				   i != getCurrentLayer().features.end();
				   ++i)
			{
				if ( i->getPosition()[RT] >= area.min()[1] &&
						 i->getPosition()[RT] <= area.max()[1] &&
						 i->getPosition()[MZ] >= area.min()[0] &&
						 i->getPosition()[MZ] <= area.max()[0] &&
						 i->getIntensity()    >= getCurrentLayer().min_int &&
						 i->getIntensity()    <= getCurrentLayer().max_int )
				{
					if (i->getIntensity() > max_int)
					{
						max_int = i->getIntensity();
						
						tmp_peak_.setIntensity(i->getIntensity());
						tmp_peak_.getPosition()[0] = i->getPosition()[1];
						tmp_peak_.getPosition()[1] = i->getPosition()[0];
						tmp_peak_.getConvexHulls() = i->getConvexHulls();
						
						max_peak = &tmp_peak_;
					}				
				}
			}	 	
		}
//		if (max_peak!=0)
//		{
//			cout << "MAX PEAK: " << max_peak->getPosition() << endl;
//		}
		return max_peak;
	}
	
	float Spectrum2DCanvas::betweenFactor_(float v1, float v2, float val)
	{
		float lo = min(v1, v2);
		float hi = max(v1, v2);
		return (hi - lo == 0) ? 1 : (val - lo) / (hi - lo);
	}
	
	const QColor& Spectrum2DCanvas::heightColor_(float val, const MultiGradient& gradient)
	{
		switch (intensity_mode_)
		{
			case IM_NONE:
				return gradient.precalculatedColorAt(val);
				break;
			case IM_LOG:
				return gradient.precalculatedColorAt(log(val+1)); //prevent log of numbers samller than 1
				break;
			case IM_PERCENTAGE:
				return gradient.precalculatedColorAt(val*percentage_factor_);
				break;
			case IM_SNAP:
				return gradient.precalculatedColorAt(val*snap_factor_);
				break;
			default:
				throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
	}
		
	void Spectrum2DCanvas::calculateMarchingSquareMatrix_(UnsignedInt layer_index)
	{
		SignedInt steps = getPrefAsInt("Preferences:2D:MarchingSquaresSteps");
		const double cell_width = visible_area_.width() / steps;
		const double cell_height = visible_area_.height() / steps;
		const double half_width = cell_width / 2.0f;
		const double half_height = cell_height / 2.0f;
	
		float y = visible_area_.minY();
		int i, j;
		for (i = 0; i <= steps; i++)
		{
			float x = visible_area_.minX();
			vector<float> line;
			for (j = 0; j <= steps; j++)
			{
				// build sum of all peak heights in the current cell
				float sum = 0.0f;

				for (ExperimentType::ConstAreaIterator i = getPeakData(layer_index).areaBeginConst(y - half_height, y + half_height, x - half_width, x + half_width); 
						 i != getPeakData(layer_index).areaEndConst(); 
						 ++i)
				{
					sum += i->getIntensity();
				}
				// log mode
				if (intensity_mode_ == IM_LOG)
				{
					sum = log(sum+1); // prevent log of numbers smaller than one
				}
				// store max
				if (sum > max_values_[layer_index])
				{
					max_values_[layer_index] = sum;
				}
				
				line.push_back(sum);
				//cout << sum << " ";
				x += cell_width;
			}
			//cout << endl;
			marching_squares_matrices_[layer_index].push_back(line);
			y += cell_height;
		}
		//cout << "rows: " << marching_squares_matrices_[layer_index].size() << " cols: " << marching_squares_matrices_[layer_index][0].size() << endl;
	}
	

	
	void Spectrum2DCanvas::paintDots_(UnsignedInt layer_index, QPainter& painter)
	{
#ifdef TIMING_TOPPVIEW
		QTime timer;
		timer.start();
#endif	
//		// Create this dot shape: 
//		// .##.
//		// ####
//		// ####
//		// .##.
//		// This is much faster than always drawing circles
//		QBitmap dotshape(4,4);
//		dotshape.fill(Qt::color1);
//		QPainter dotshape_painter(&dotshape);
//		dotshape_painter.setPen(Qt::color0);
//		dotshape_painter.drawPoint(0,0);
//		dotshape_painter.drawPoint(3,3);
//		dotshape_painter.drawPoint(0,3);
//		dotshape_painter.drawPoint(3,0);
//		
//		//shift of the bitmap
//		QPoint shift(2,2);
		
		if (intensity_mode_ == IM_PERCENTAGE)
		{
			if (getLayer(layer_index).type == LayerData::DT_PEAK)
			{
				percentage_factor_ = overall_data_range_.max()[2]/getPeakData(layer_index).getMaxInt();
			}
			else
			{
				percentage_factor_ = overall_data_range_.max()[2]/getLayer(layer_index).features.getMaxInt();
			}
		}
		else 
		{
			percentage_factor_ = 1.0;
		}
		
		//tmporary variable
		QPoint pos;
		
		double min_int = getLayer(layer_index).min_int;
		double max_int = getLayer(layer_index).max_int;
		SignedInt mode = getDotMode();
		
		painter.setPen(Qt::black);

		if (getLayer(layer_index).type==LayerData::DT_PEAK) //peaks
		{
			for (ExperimentType::ConstAreaIterator i = getPeakData(layer_index).areaBeginConst(visible_area_.min()[1],visible_area_.max()[1],visible_area_.min()[0],visible_area_.max()[0]); 
					 i != getPeakData(layer_index).areaEndConst(); 
					 ++i)
			{
				if (i->getIntensity()>=min_int && i->getIntensity()<=max_int)
				{
					if (mode==DOT_GRADIENT)
					{
						painter.setPen(heightColor_(i->getIntensity(), dot_gradient_));
					}
					dataToWidget_(i->getPos(), i.getRetentionTime(),pos);
					painter.drawLine(pos.x(),pos.y()-1,pos.x(),pos.y()+1);
					painter.drawLine(pos.x()-1,pos.y(),pos.x()+1,pos.y());
				}
			}
		}
		else //features
		{
			for (FeatureMapType::ConstIterator i = getLayer(layer_index).features.begin();
				   i != getLayer(layer_index).features.end();
				   ++i)
			{
				if ( i->getPosition()[RT] >= visible_area_.min()[1] &&
						 i->getPosition()[RT] <= visible_area_.max()[1] &&
						 i->getPosition()[MZ] >= visible_area_.min()[0] &&
						 i->getPosition()[MZ] <= visible_area_.max()[0] &&
						 i->getIntensity()>=min_int &&
						 i->getIntensity()<=max_int)
				{
					if (mode==DOT_GRADIENT)
					{
						painter.setPen(heightColor_(i->getIntensity(), dot_gradient_));
					}
					dataToWidget_(i->getPosition()[MZ],i->getPosition()[RT],pos);
					painter.drawLine(pos.x(),pos.y()-1,pos.x(),pos.y()+1);
					painter.drawLine(pos.x()-1,pos.y(),pos.x()+1,pos.y());
				}
			}
		}

#ifdef TIMING_TOPPVIEW
		cout << "paintDots_ took " << timer.elapsed() << " ms" << endl;
#endif	
	}

	void Spectrum2DCanvas::paintConvexHulls_(UnsignedInt layer_index, QPainter& painter)
	{
		painter.setPen(Qt::black);

		double min_int = getLayer(layer_index).min_int;
		double max_int = getLayer(layer_index).max_int;

		for (FeatureMapType::ConstIterator i = getLayer(layer_index).features.begin();
			   i != getLayer(layer_index).features.end();
			   ++i)
		{
			if ( i->getPosition()[RT] >= visible_area_.min()[1] &&
					 i->getPosition()[RT] <= visible_area_.max()[1] &&
					 i->getPosition()[MZ] >= visible_area_.min()[0] &&
					 i->getPosition()[MZ] <= visible_area_.max()[0] &&
					 i->getIntensity()>=min_int &&
					 i->getIntensity()<=max_int)
			{
				paintConvexHulls_(i->getConvexHulls(),painter);
			}
		}
	}            

	void Spectrum2DCanvas::paintFeaturePairConnections_(UnsignedInt layer_index, QPainter& painter)
	{
		painter.setPen(Qt::black);

		double min_int = getLayer(layer_index).min_int;
		double max_int = getLayer(layer_index).max_int;

		QPoint line_begin, line_end;
		FeatureMapType::ConstIterator i2;
	
		for (FeatureMapType::ConstIterator i1 = getLayer(layer_index).features.begin();
			   i1 != getLayer(layer_index).features.end();
			   i1+=2)
		{
			//get second feature
			i2 = i1 + 1;
			
			if ( i1->getPosition()[RT] >= visible_area_.min()[1] &&
					 i1->getPosition()[RT] <= visible_area_.max()[1] &&
					 i1->getPosition()[MZ] >= visible_area_.min()[0] &&
					 i1->getPosition()[MZ] <= visible_area_.max()[0] &&
					 i1->getIntensity()>=min_int &&
					 i1->getIntensity()<=max_int &&
					 i2->getPosition()[RT] >= visible_area_.min()[1] &&
					 i2->getPosition()[RT] <= visible_area_.max()[1] &&
					 i2->getPosition()[MZ] >= visible_area_.min()[0] &&
					 i2->getPosition()[MZ] <= visible_area_.max()[0] &&
					 i2->getIntensity()>=min_int &&
					 i2->getIntensity()<=max_int
					 )
			{
				dataToWidget_(i1->getPosition()[MZ],i1->getPosition()[RT], line_begin);
				dataToWidget_(i2->getPosition()[MZ],i2->getPosition()[RT], line_end);
				painter.drawLine(line_begin, line_end);
			}
		}
	}      

  void Spectrum2DCanvas::paintConvexHulls_(const DFeature<2>::ConvexHullVector& hulls, QPainter& painter)
  {
		QPolygon points;
		
		//iterate over all convex hulls
		for (UnsignedInt hull=0; hull<hulls.size(); ++hull)
		{
			points.resize(hulls[hull].getPoints().size());
			UnsignedInt index=0;
			QPoint pos;
			//iterate over hull points
			for(DFeature<2>::ConvexHullType::PointArrayType::const_iterator it=hulls[hull].getPoints().begin(); it!=hulls[hull].getPoints().end(); ++it, ++index)
			{
				dataToWidget_(it->Y(), it->X(),pos);
				points.setPoint(index, pos);
			}	
			//cout << "Hull: " << hull << " Points: " << points.size()<<endl;
			painter.drawPolygon(points);
		}
  }

	void Spectrum2DCanvas::paintContours_(UnsignedInt layer_index, QPainter& painter)
	{
		if (max_values_[layer_index] == 0) return;
		
		painter.setPen(Qt::black);
		
		//intensity steps where lines are drawn (valid for all the layer)
		float intensity_step = max_values_[layer_index] / getPrefAsInt("Preferences:2D:Contour:Lines");
		
		//calculate data/pixel width and height or a cell
		SignedInt steps = getPrefAsInt("Preferences:2D:MarchingSquaresSteps");
		SignedInt pixel_width = width() / steps;
		SignedInt pixel_height = height() / steps;
		float data_width = visible_area_.width() / steps;
		float data_height = visible_area_.height() / steps;
		
		QPoint cell_pos;
		
		// draw the lines
		float y = visible_area_.minY();
		for (int i = 0; i < steps; i++)
		{
			float x = visible_area_.minX();
			for (int j = 0; j < steps; j++)
			{
				float left_bottom = marching_squares_matrices_[layer_index][i][j];
				float right_bottom = 0;
				float left_top = 0;
				float right_top = marching_squares_matrices_[layer_index][i + 1][j + 1];

				if ( isMzToXAxis() )
				{
				  right_bottom = marching_squares_matrices_[layer_index][i][j + 1];
				  left_top = marching_squares_matrices_[layer_index][i + 1][j];
				}
				else
				{
				  left_top = marching_squares_matrices_[layer_index][i][j + 1];
				  right_bottom = marching_squares_matrices_[layer_index][i + 1][j];		
				}
	
				dataToWidget_(x, y, cell_pos);
				cell_pos.setY(cell_pos.y() - (pixel_height+1));

				const float minimum = min(left_top, min(right_top, min(left_bottom, right_bottom)));
				const float maximum = max(left_top, max(right_top, max(left_bottom, right_bottom)));
				for (float height = ceil(minimum / intensity_step) * intensity_step; height <= maximum; height += intensity_step)
				{
					// this bitset indicates which points are above the height threshold
					int state = (left_top > height) << 3 |
											(right_top > height) << 2 |
											(left_bottom > height) << 1 |
											(right_bottom > height);
	
					// this is the ugly marching squares case differentiation.
					switch (state)
					{
						default:
						case 0:
						case 15:
						{
							// no line to draw
							break;
						}
	
						case 1:
						{
							painter.drawLine(cell_pos.x() + int(betweenFactor_(left_bottom, right_bottom, height) * pixel_width),
													cell_pos.y() + pixel_height,
													cell_pos.x() + pixel_width,
													cell_pos.y() + int(betweenFactor_(right_top, right_bottom, height) * pixel_height));
							break;
						}
						case 14:
						{
							painter.drawLine(cell_pos.x() + pixel_width - int(betweenFactor_(left_bottom, right_bottom, height) * pixel_width),
													cell_pos.y() + pixel_height,
													cell_pos.x() + pixel_width,
													cell_pos.y() + pixel_height - int(betweenFactor_(right_top, right_bottom, height) * pixel_height));
							break;
						}
						case 2:
						{
							painter.drawLine(cell_pos.x(),
													cell_pos.y() + int(betweenFactor_(left_top, left_bottom, height) * pixel_height),
													cell_pos.x() + pixel_width - int(betweenFactor_(left_bottom, right_bottom, height) * pixel_width),
													cell_pos.y() + pixel_height);
							break;
						}
						case 13:
						{
							painter.drawLine(cell_pos.x(),
													cell_pos.y() + pixel_height - int(betweenFactor_(left_top, left_bottom, height) * pixel_height),
													cell_pos.x() + int(betweenFactor_(left_bottom, right_bottom, height) * pixel_width),
	
													cell_pos.y() + pixel_height);
							break;
						}
						case 3:
						{
							painter.drawLine(cell_pos.x(),
													cell_pos.y() + int(betweenFactor_(left_top, left_bottom, height) * pixel_height),
													cell_pos.x() + pixel_width,
													cell_pos.y() + int(betweenFactor_(right_top, right_bottom, height) * pixel_height));
	
							break;
						}
						case 12:
						{
							painter.drawLine(cell_pos.x(),
	
													cell_pos.y() + pixel_height - int(betweenFactor_(left_top, left_bottom, height) * pixel_height),
													cell_pos.x() + pixel_width,
													cell_pos.y() + pixel_height - int(betweenFactor_(right_top, right_bottom, height) * pixel_height));
							break;
						}
						case 4:
						{
							painter.drawLine(cell_pos.x() + int(betweenFactor_(left_top, right_top, height) * pixel_width),
													cell_pos.y(),
													cell_pos.x() + pixel_width,
													cell_pos.y() + pixel_height - int(betweenFactor_(right_top, right_bottom, height) * pixel_height));
							break;
						}
						case 11:
						{
							painter.drawLine(cell_pos.x() + pixel_width - int(betweenFactor_(left_top, right_top, height) * pixel_width),
													cell_pos.y(),
													cell_pos.x() + pixel_width,
													cell_pos.y() + int(betweenFactor_(right_top, right_bottom, height) * pixel_height));
							break;
						}
						case 5:
						{
							painter.drawLine(cell_pos.x() + int(betweenFactor_(left_top, right_top, height) * pixel_width),
													cell_pos.y(),
													cell_pos.x() + int(betweenFactor_(left_bottom, right_bottom, height) * pixel_width),
													cell_pos.y() + pixel_height);
							break;
						}
						case 10:
						{
							painter.drawLine(cell_pos.x() + pixel_width - int(betweenFactor_(left_top, right_top, height) * pixel_width),
													cell_pos.y(),
													cell_pos.x() + pixel_width - int(betweenFactor_(left_bottom, right_bottom, height) * pixel_width),
													cell_pos.y() + pixel_height);
							break;
						}
						case 6:
						{
							painter.drawLine(cell_pos.x(),
													cell_pos.y() + int(betweenFactor_(left_top, left_bottom, height) * pixel_height),
													cell_pos.x() + pixel_width - int(betweenFactor_(left_bottom, right_bottom, height) * pixel_width),
													cell_pos.y() + pixel_height);
							painter.drawLine(cell_pos.x() + int(betweenFactor_(left_top, right_top, height) * pixel_width),
													cell_pos.y(),
													cell_pos.x() + pixel_width,
													cell_pos.y() + pixel_height - int(betweenFactor_(right_top, right_bottom, height) * pixel_height));
							break;
						}
						case 9:
						{
							painter.drawLine(cell_pos.x() + int(betweenFactor_(left_bottom, right_bottom, height) * pixel_width),
													cell_pos.y() + pixel_height,
													cell_pos.x() + pixel_width,
													cell_pos.y() + int(betweenFactor_(right_top, right_bottom, height) * pixel_height));
							painter.drawLine(cell_pos.x(),
													cell_pos.y() + pixel_height - int(betweenFactor_(left_top, left_bottom, height) * pixel_height),
													cell_pos.x() + pixel_width - int(betweenFactor_(left_top, right_top, height) * pixel_width),
													cell_pos.y());
							break;
						}
	
						case 7:
						{
							painter.drawLine(cell_pos.x(),
													cell_pos.y() + int(betweenFactor_(left_top, left_bottom, height) * pixel_height),
													cell_pos.x() + int(betweenFactor_(left_top, right_top, height) * pixel_width),
													cell_pos.y());
							break;
						}
						case 8:
						{
							painter.drawLine(cell_pos.x(),
													cell_pos.y() + pixel_height - int(betweenFactor_(left_top, left_bottom, height) * pixel_height),
													cell_pos.x() + pixel_width - int(betweenFactor_(left_top, right_top, height) * pixel_width),
													cell_pos.y());
	
							break;
						}
					}
				}
				x += data_width;
			}
			y += data_height;
		}
	}
	
	void Spectrum2DCanvas::paintSurface_(UnsignedInt layer_index, QPainter& painter)
	{
		if (max_values_[layer_index] == 0) return;
		
		QImage image(buffer_.width(), buffer_.height(),QImage::Format_RGB32);
		QRgb* image_start = reinterpret_cast<QRgb*>(image.scanLine(0));
		const uint image_line_diff = reinterpret_cast<QRgb*>(image.scanLine(1)) - image_start;

		//calculate data/pixel width and height or a cell
		SignedInt steps = getPrefAsInt("Preferences:2D:MarchingSquaresSteps");
		SignedInt pixel_width = width() / steps;
		SignedInt pixel_height = height() / steps;
		float data_width = visible_area_.width() / steps;
		float data_height = visible_area_.height() / steps;
		
		//cout << "Pixel size: " << pixel_width << " " << pixel_height <<endl;
		//cout << "Data size: " << data_width << " " << data_height <<endl;
		
		//construct color matrix
		vector<vector<const QColor*> > color_matrix;
		//cout << "color matrix for layer " << layer_index << ": " << endl;
		for (SignedInt i = 0; i <= steps; i++)
		{
			color_matrix.insert(color_matrix.end(), vector<const QColor*>() );
			for (SignedInt j = 0; j <= steps; j++)
			{
				color_matrix.back().push_back(&heightColor_(marching_squares_matrices_[layer_index][i][j] / max_values_[layer_index] * overall_data_range_.max()[2], surface_gradient_));
				//cout << 255.0 - color_matrix.back().back()->red() << " ";
			}
			//cout << endl;
		}
		
		//draw
		QRgb* pixel;
		float y = visible_area_.minY();
		QPoint cell_pos;
		for (int i = 0; i < steps; i++)
		{
			float x = visible_area_.minX();
			for (int j = 0; j < steps; j++)
			{
				const QColor* left_top = 0;
				const QColor* right_top = color_matrix[i + 1][j + 1];
				const QColor* left_bottom = color_matrix[i][j];
				const QColor* right_bottom = 0;
	
				if ( isMzToXAxis() )
				{
					right_bottom = color_matrix[i][j + 1];
					left_top = color_matrix[i + 1][j];
				}
				else
				{
					left_top = color_matrix[i][j + 1];
					right_bottom = color_matrix[i + 1][j];					
				}
				dataToWidget_(x, y, cell_pos);
				cell_pos.setY(cell_pos.y() - (pixel_height+1));
				
				//cout << cell_pos.x() << " " << cell_pos.y() << endl;
				
				int left_red = left_top->red() << 8;
				int left_green = left_top->green() << 8;
				int left_blue = left_top->blue() << 8;
				int right_red = right_top->red() << 8;
				int right_green = right_top->green() << 8;
				int right_blue = right_top->blue() << 8;
	
				const int left_d_red = ((left_bottom->red() << 8) - left_red) / pixel_height;
				const int left_d_green = ((left_bottom->green() << 8) - left_green) / pixel_height;
				const int left_d_blue = ((left_bottom->blue() << 8) - left_blue) / pixel_height;
				const int right_d_red = ((right_bottom->red() << 8) - right_red) / pixel_height;
				const int right_d_green = ((right_bottom->green() << 8) - right_green) / pixel_height;
				const int right_d_blue = ((right_bottom->blue() << 8) - right_blue) / pixel_height;
	
				pixel = image_start;
				pixel += cell_pos.y() * buffer_.width() + cell_pos.x();
	
				for (int py = 0; py !=pixel_height + 1; py++)
				{
					QRgb* start_pixel = pixel;
	
					// vertical clipping
					if (cell_pos.y() + py >= 0 && cell_pos.y() + py < buffer_.height())
					{
						const int d_red = (right_red - left_red) / pixel_width;
						const int d_green = (right_green - left_green) / pixel_width;
	
						const int d_blue = (right_blue - left_blue) / pixel_width;
	
						int c_red = left_red;
						int c_green = left_green;
						int c_blue = left_blue;
	
						for (int px = 0; px != pixel_width + 1; px++)
						{
							// horizontal clipping
							if (cell_pos.x() + px >= 0 && cell_pos.x() + px < buffer_.width())
							{
								*pixel = qRgb(c_red >> 8, c_green >> 8, c_blue >> 8);
							}
	
							pixel++;
	
							c_red += d_red;
							c_green += d_green;
							c_blue += d_blue;
						}
					}
	
					left_red += left_d_red;
					left_green += left_d_green;
					left_blue += left_d_blue;
	
					right_red += right_d_red;
					right_green += right_d_green;
					right_blue += right_d_blue;
	
					// next line
					pixel = start_pixel + image_line_diff;
				}
				x += data_width;
			}
			y += data_height;
		}
		painter.drawImage(0, 0, image);
	}

	void Spectrum2DCanvas::intensityDistributionChange_()
	{
		update_buffer_ = true;
		update();
	}
	
	void Spectrum2DCanvas::intensityModeChange_()
	{
		recalculateDotGradient_();
		recalculateSurfaceGradient_();
		SpectrumCanvas::intensityModeChange_();
	}
	
	void Spectrum2DCanvas::recalculateDotGradient_()
	{
		//cout << "recalculateDotGradient_" << endl;
		if (intensity_mode_ == IM_LOG)
		{
			//cout << "LOG:" <<" "<< log(overall_data_range_.min()[2]) <<" "<< log(overall_data_range_.max()[2])<<" "<<getPrefAsInt("Preferences:2D:InterpolationSteps")<<endl;
			dot_gradient_.activatePrecalculationMode(0, log(overall_data_range_.max()[2]+1), getPrefAsInt("Preferences:2D:InterpolationSteps"));
		}
		else
		{
			//cout << "NORMAL:" << overall_data_range_.min()[2] <<" "<< overall_data_range_.max()[2]<<" "<<getPrefAsInt("Preferences:2D:InterpolationSteps")<<endl;
			dot_gradient_.activatePrecalculationMode(0, overall_data_range_.max()[2], getPrefAsInt("Preferences:2D:InterpolationSteps"));
		}	
	}
	
	void Spectrum2DCanvas::recalculateSurfaceGradient_()
	{
		if (intensity_mode_ == IM_LOG)
		{
			surface_gradient_.activatePrecalculationMode(0, log(overall_data_range_.max()[2]+1), getPrefAsInt("Preferences:2D:InterpolationSteps"));
		}
		else
		{
			surface_gradient_.activatePrecalculationMode(0, overall_data_range_.max()[2], getPrefAsInt("Preferences:2D:InterpolationSteps"));		
		}	
	}
	
	void Spectrum2DCanvas::createProjections_(const AreaType& area, bool shift_pressed, bool ctrl_pressed)
	{
		ExperimentType::CoordinateType rt_l = area.minY();
		ExperimentType::CoordinateType rt_h = area.maxY();
		ExperimentType::CoordinateType mz_l = area.minX();
		ExperimentType::CoordinateType mz_h = area.maxX();

		//cout << "Projection: "<< rt_l << " " << rt_h << " " << mz_l << " " << mz_h << endl;
		
		if (shift_pressed)
		{
			if (isMzToXAxis())
			{
				mz_l = visible_area_.min()[0];
				mz_h = visible_area_.max()[0];	
			}
			else
			{
				rt_l = visible_area_.min()[1];
				rt_h = visible_area_.max()[1];
			}
		}
		else if (ctrl_pressed)
		{
			if (isMzToXAxis())
			{
				rt_l = visible_area_.min()[1];
				rt_h = visible_area_.max()[1];	
			}
			else
			{
				mz_l = visible_area_.min()[0];
				mz_h = visible_area_.max()[0];
			}
		}
		
		//cout << "Projection (buttons): "<< rt_l << " " << rt_h << " " << mz_l << " " << mz_h << endl;
		
		//create projection data
		map<float, float> mz, rt;
		for (ExperimentType::ConstAreaIterator i = getCurrentPeakData().areaBeginConst(rt_l,rt_h,mz_l,mz_h); 
				 i != getCurrentPeakData().areaEndConst();
				 ++i)
		{
			if (i->getIntensity()>=getCurrentLayer().min_int && i->getIntensity()<=getCurrentLayer().max_int)
			{
				mz[i->getPos()] += i->getIntensity();
				rt[i.getRetentionTime()] += i->getIntensity();
			}
		}
		
		// write to spectra
		MSExperiment<>::SpectrumType::ContainerType& cont_mz = projection_mz_[0].getContainer();
		MSExperiment<>::SpectrumType::ContainerType& cont_rt = projection_rt_[0].getContainer();
		
		//resize and add boundary peaks		
		cont_mz.resize(mz.size()+2);
		cont_mz[0].setPos(mz_l);
		cont_mz[0].setIntensity(0.0);
		cont_mz[1].setPos(mz_h);
		cont_mz[1].setIntensity(0.0);
		cont_rt.resize(rt.size()+2);
		cont_rt[0].setPos(rt_l);
		cont_rt[0].setIntensity(0.0);
		cont_rt[1].setPos(rt_h);
		cont_rt[1].setIntensity(0.0);
		
		UnsignedInt i = 2;
		for (map<float, float>::iterator it = mz.begin(); it != mz.end(); ++it)
		{
			cont_mz[i].setPos(it->first);
			cont_mz[i].setIntensity(it->second);
			++i;
		}

		i = 2;
		for (map<float, float>::iterator it = rt.begin(); it != rt.end(); ++it)
		{
			cont_rt[i].setPos(it->first);
			cont_rt[i].setIntensity(it->second);
			++i;
		}
		
		if (isMzToXAxis())
		{
			emit showProjectionHorizontal(projection_mz_);
			emit showProjectionVertical(projection_rt_);	
		}
		else
		{
			emit showProjectionHorizontal(projection_rt_);
			emit showProjectionVertical(projection_mz_);
		}
	}

	
	PreferencesDialogPage* Spectrum2DCanvas::createPreferences(QWidget* parent)
	{
		return new Spectrum2DCanvasPDP(this, parent);
	}
	
	void Spectrum2DCanvas::setDotMode(SignedInt mode)
	{
		prefs_.setValue("Preferences:2D:Dot:Mode", mode);
	}
	
	SignedInt Spectrum2DCanvas::getDotMode()
	{
		if (prefs_.getValue("Preferences:2D:Dot:Mode").isEmpty())
		{
			return 0;
		}
		
		return SignedInt(prefs_.getValue("Preferences:2D:Dot:Mode"));
	}
	
	void Spectrum2DCanvas::setDotGradient(const string& gradient)
	
	{
		prefs_.setValue("Preferences:2D:Dot:Gradient",gradient);
		dot_gradient_.fromString(gradient);
		recalculateDotGradient_();
	}
	
	void Spectrum2DCanvas::setSurfaceGradient(const string& gradient)
	{
		prefs_.setValue("Preferences:2D:Surface:Gradient",gradient);
		surface_gradient_.fromString(gradient);
		recalculateSurfaceGradient_();
	}
	
	void Spectrum2DCanvas::setMainPreferences(const Param& prefs)
	{
		SpectrumCanvas::setMainPreferences(prefs);
		surface_gradient_.fromString(getPrefAsString("Preferences:2D:Surface:Gradient"));
		recalculateSurfaceGradient_();
		dot_gradient_.fromString(getPrefAsString("Preferences:2D:Dot:Gradient"));
		recalculateDotGradient_();
		if (getPrefAsString("Preferences:2D:Mapping:MappingOfMzTo") != "X-Axis")
		{
			mzToXAxis(false);
		}
	}
	
	SignedInt Spectrum2DCanvas::finishAdding(float low_intensity_cutoff)
	{
		current_layer_ = getLayerCount()-1;

		//set visibility to true
		show_contours_.push_back(false);
		show_surface_.push_back(false);
		show_dots_.push_back(true);
		
		if (layers_.back().type==LayerData::DT_PEAK) //Feature data
		{
			currentPeakData_().sortSpectra(true);
			currentPeakData_().updateRanges(1);
			getCurrentLayer_().min_int = low_intensity_cutoff;
			getCurrentLayer_().max_int = getCurrentPeakData().getMaxInt();
			update_buffer_ = true;
		}
		else //peak data
		{

			getCurrentLayer_().features.updateRanges();
			getCurrentLayer_().max_int = getCurrentLayer().features.getMaxInt();
		}
		
		//overall values update
		updateRanges_(current_layer_,0,1,2);

		if (getLayerCount()==1)
		{
			AreaType tmp_area;
			tmp_area.assign(overall_data_range_);
			visible_area_ = tmp_area;
			emit visibleAreaChanged(tmp_area);
		}
		else
		{
			resetZoom();
		}

		intensityModeChange_();
		
		emit sendStatusMessage("",0);
		
		emit layerActivated(this);

		return current_layer_;
	}
	
	void Spectrum2DCanvas::removeLayer(int layer_index )
	{
		if (layer_index<0 || layer_index >= int(getLayerCount()))
		{
			return;
		}

		//unselect all peaks
		delete(selected_peak_);
		selected_peak_ = 0;
		delete(measurement_start_);
		measurement_start_ = 0;
		delete(measurement_stop_);
		measurement_stop_ = 0;
	
		//remove the data
		layers_.erase(layers_.begin()+layer_index);
		
		//remove settings
		show_contours_.erase(show_contours_.begin()+layer_index);
		show_surface_.erase(show_surface_.begin()+layer_index);
		show_dots_.erase(show_dots_.begin()+layer_index);
		
		//update visible area and boundaries
		recalculateRanges_(0,1,2);

		AreaType tmp;
		tmp.assign(overall_data_range_);
		if (tmp != visible_area_)
		{
			visible_area_.assign(overall_data_range_);
		}

		//update current layer
		if (current_layer_!=0 && current_layer_ >= getLayerCount())
		{
			current_layer_ = getLayerCount()-1;
		}

		if (layers_.empty())
		{
			return;
		}

		intensityModeChange_();
		emit layerActivated(this);
	}
	
	//change the current layer
	void Spectrum2DCanvas::activateLayer(int layer_index)
	{
		if (layer_index<0 || layer_index >= int(getLayerCount()) || layer_index==int(current_layer_))
		{
			return ;
		}
		current_layer_ = layer_index;
		emit layerActivated(this);
		
		//unselect all peaks
		delete(selected_peak_);
		selected_peak_ = 0;
		delete(measurement_start_);
		measurement_start_ = 0;
		delete(measurement_stop_);
		measurement_stop_ = 0;
		
		update();
	}
	
	void Spectrum2DCanvas::repaintAll()
	{				
		recalculateDotGradient_();
		recalculateSurfaceGradient_();
		SpectrumCanvas::repaintAll();
	}
	
	void Spectrum2DCanvas::recalculateSnapFactor_()
	{
		if (intensity_mode_ == IM_SNAP) 
		{
			double local_max  = -numeric_limits<double>::max();
			for (UnsignedInt i=0; i<getLayerCount(); i++)
			{
				if (getLayer(i).visible)
				{
					if (getLayer(i).type==LayerData::DT_PEAK)
					{
						for (ExperimentType::ConstAreaIterator it = getPeakData(i).areaBeginConst(visible_area_.min()[1],visible_area_.max()[1],visible_area_.min()[0],visible_area_.max()[0]); 
								 it != getPeakData(i).areaEndConst(); 
								 ++it)
						{
							if (it->getIntensity() > local_max && it->getIntensity()>=getLayer(i).min_int && it->getIntensity()<=getLayer(i).max_int)
							{
								local_max = it->getIntensity();
							}
						}
					}
					else //features
					{
						for (FeatureMapType::ConstIterator it = getLayer(i).features.begin();
							   it != getLayer(i).features.end();
							   ++it)
						{
							if ( it->getPosition()[RT] >= visible_area_.min()[1] &&
									 it->getPosition()[RT] <= visible_area_.max()[1] &&
									 it->getPosition()[MZ] >= visible_area_.min()[0] &&
									 it->getPosition()[MZ] <= visible_area_.max()[0] &&
									 it->getIntensity()>=getLayer(i).min_int && 
									 it->getIntensity()<=getLayer(i).max_int &&
									 it->getIntensity() > local_max)
							{
								local_max = it->getIntensity();
							}
						}
					}
				}
			}
			snap_factor_ = overall_data_range_.max()[2]/local_max;			
		}
		else
		{ 
			snap_factor_ = 1.0;
		}
	}

	void Spectrum2DCanvas::updateScrollbars_()
	{
		if (isMzToXAxis())
		{
			emit updateHScrollbar(overall_data_range_.min()[0],visible_area_.min()[0],visible_area_.max()[0],overall_data_range_.max()[0]);
			emit updateVScrollbar(overall_data_range_.min()[1],visible_area_.min()[1],visible_area_.max()[1],overall_data_range_.max()[1]);
		}
		else
		{
			emit updateVScrollbar(overall_data_range_.min()[0],visible_area_.min()[0],visible_area_.max()[0],overall_data_range_.max()[0]);
			emit updateHScrollbar(overall_data_range_.min()[1],visible_area_.min()[1],visible_area_.max()[1],overall_data_range_.max()[1]);
		}
	}

	void Spectrum2DCanvas::horizontalScrollBarChange(int value)
	{
		AreaType new_area = visible_area_;
		if (isMzToXAxis())
		{
			new_area.setMinX(value);
			new_area.setMaxX(value + (visible_area_.maxX() - visible_area_.minX()));
			changeVisibleArea_(new_area);
		}
		else
		{
			new_area.setMinY(value);
			new_area.setMaxY(value + (visible_area_.maxY() - visible_area_.minY()));
			changeVisibleArea_(new_area);
		}
	}

	void Spectrum2DCanvas::verticalScrollBarChange(int value)
	{
		AreaType new_area = visible_area_;
		if (!isMzToXAxis())
		{
			double range = (overall_data_range_.maxX() - overall_data_range_.minX())- (visible_area_.maxX() - visible_area_.minX());
			double newval = (1.0 - (double(value) - overall_data_range_.minX()) / range )* range + overall_data_range_.minX();
			//cout << value << " " << newval << " " << newval + (visible_area_.maxX() - visible_area_.minX()) << endl;
			//cout << "Min: " <<  overall_data_range_.minX() << " Range: " << range << endl << endl;
			new_area.setMinX(newval);
			new_area.setMaxX(newval + (visible_area_.maxX() - visible_area_.minX()));
			changeVisibleArea_(new_area);
		}
		else
		{
			double range = (overall_data_range_.maxY() - overall_data_range_.minY())- (visible_area_.maxY() - visible_area_.minY());
			double newval = (1.0 - (double(value) - overall_data_range_.minY()) / range )* range + overall_data_range_.minY();
			//cout << value << " " << newval << " " << newval + (visible_area_.maxY() - visible_area_.minY()) << endl;
			//cout << "Min: " <<  overall_data_range_.minY() << " Range: " << range << endl << endl;
			new_area.setMinY(newval);
			new_area.setMaxY(newval + (visible_area_.maxY() - visible_area_.minY()));
			changeVisibleArea_(new_area);
		}
	}

	void Spectrum2DCanvas::paintEvent(QPaintEvent* e)
	{		
#ifdef DEBUG_TOPPVIEW
		cout << "BEGIN " << __PRETTY_FUNCTION__ << endl;
	  cout << "  Visible area -- m/z: " << visible_area_.minX() << " - " << visible_area_.maxX() << " int: " << visible_area_.minY() << " - " << visible_area_.maxY() << endl;
	  cout << "  Overall area -- m/z: " << overall_data_range_.min()[0] << " - " << overall_data_range_.max()[0] << " int: " << overall_data_range_.min()[1] << " - " << overall_data_range_.max()[1] << endl; 
#endif
#ifdef TIMING_TOPPVIEW
		QTime timer;
 		timer.start();
#endif
		
		QPainter painter;
		
		if (update_buffer_)
		{
			update_buffer_ = false;
			
			//recalculate marching squares matices for layers with visible surfaces and contour lines
			marching_squares_matrices_.clear();
			marching_squares_matrices_.resize(getLayerCount());
			max_values_.clear();
			max_values_.resize(getLayerCount());
			for (UnsignedInt i=0; i<getLayerCount(); i++)
			{
				if ( getLayer(i).type == LayerData::DT_PEAK && getLayer(i).visible && (show_surface_[i] || show_contours_[i]))
				{
					calculateMarchingSquareMatrix_(i);
				}
			}
			
			//recalculate snap factor
			recalculateSnapFactor_();
			
			//draw buffer in image in order to have composition
			//QImage image(buffer_.width(),buffer_.height(),QImage::Format_ARGB32_Premultiplied);
			buffer_.fill(QColor(getPrefAsString("Preferences:2D:BackgroundColor").c_str()).rgb());
			painter.begin(&buffer_);
			for (UnsignedInt i=0; i<getLayerCount(); i++)
			{
				if (getLayer(i).visible)
				{
					if (getLayer(i).type==LayerData::DT_PEAK)
					{
						if (show_surface_[i])
						{
							//cout << "surface peak layer: " << i << endl;
							paintSurface_(i, painter);
						}
						if (show_contours_[i])
						{
							//cout << "countour peak layer: " << i << endl;
							paintContours_(i, painter);
						}
						if (show_dots_[i])
						{
							//cout << "dot peak layer: " << i << endl;
							paintDots_(i, painter);
						}
					}
					else if (getLayer(i).type==LayerData::DT_FEATURE)
					{
						//cout << "dot feature layer: " << i << endl;
						paintDots_(i, painter);
						paintConvexHulls_(i, painter);
					}
					else if (getLayer(i).type==LayerData::DT_FEATURE_PAIR)
					{
						//cout << "dot feature pair layer: " << i << endl;
						paintDots_(i, painter);
						paintConvexHulls_(i, painter);
						paintFeaturePairConnections_(i, painter);
					}
				}
			}
			paintGridLines_(painter);
			painter.end();
		}
		
		painter.begin(this);
		
		//copy peak data from buffer
		QVector<QRect> rects = e->region().rects();
		for (int i = 0; i < (int)rects.size(); ++i)
		{
			painter.drawPixmap(rects[i].topLeft(), buffer_, rects[i]);
		}
		
		//draw mesaurement peak
		if (measurement_start_)
		{
			painter.setPen(Qt::black);
			
			QPoint line_begin, line_end;
			
			if (measurement_stop_)
			{
				 dataToWidget_(measurement_stop_->getPosition(), line_end);
				//cout << "Line end: " << line_end << endl;
			}
			else
			{
				line_end = last_mouse_pos_;
				//cout << "Ende: " << line_end.x() << " " << line_end.y() << endl;
			}
			dataToWidget_(measurement_start_->getPosition(), line_begin);
			painter.drawLine(line_begin, line_end);
		}
		highlightPeak_(painter, measurement_start_);
		highlightPeak_(painter, measurement_stop_);
		
		//draw selected peak
		highlightPeak_(painter, selected_peak_);
		
		//draw convex hull of selected peak
		if (selected_peak_ && getCurrentLayer().type!=LayerData::DT_PEAK)
		{
			painter.setPen(QPen(Qt::red, 2));
			paintConvexHulls_(selected_peak_->getConvexHulls(),painter);
		}
		
		painter.end();
#ifdef DEBUG_TOPPVIEW
		cout << "END   " << __PRETTY_FUNCTION__ << endl;
#endif
#ifdef TIMING_TOPPVIEW
		cout << "2D PaintEvent took " << timer.elapsed() << " ms" << endl << endl;
#endif	
	}

	void Spectrum2DCanvas::mousePressEvent(QMouseEvent* e)
	{
		last_mouse_pos_ = e->pos();
		
		switch (action_mode_)
		{
			case AM_SELECT:
				if (e->button() == Qt::LeftButton)
				{
					rubber_band_.setGeometry(e->pos().x(),e->pos().y(),0,0);
					rubber_band_.show();
				}
				break;
			case AM_MEASURE:
				if (e->button() == Qt::LeftButton)
				{
					if (selected_peak_)
					{
						delete(measurement_start_);
						measurement_start_ = new DFeature<2>(*selected_peak_);
					}
					else
					{
						delete(measurement_start_);
						measurement_start_ = 0;
					}
					delete(measurement_stop_);
					measurement_stop_ = 0;
				}
				break;
			case AM_ZOOM:
				if (e->button() == Qt::LeftButton)
				{
					rubber_band_.setGeometry(e->pos().x(),e->pos().y(),0,0);
					rubber_band_.show();
				}
				break;
			case AM_TRANSLATE:
				if (e->button() == Qt::LeftButton)
				{
					setCursor(cursor_translate_in_progress_);
				}
				break;
			default:
				break;
		}

		e->accept();
	}

	void Spectrum2DCanvas::mouseMoveEvent(QMouseEvent* e)
	{
		QPoint pos = e->pos();
		
		switch (action_mode_)
		{
			case AM_SELECT:
			{
				// highlight nearest peak
				if (e->buttons() == Qt::NoButton)
				{
					
					DFeature<2>* max_peak = findNearestPeak_(pos);
					
					if (max_peak)
					{
						// show Peak Coordinates (with intensity)
						emit sendCursorStatus(max_peak->getPosition()[0], max_peak->getIntensity(), max_peak->getPosition()[1]);
						//show lable
						string meta = max_peak->getMetaValue(3).toString();
						if (meta!="") sendStatusMessage(meta, 0);
					}
					else
					{
						//show Peak Coordinates (without intensity)
						PointType pnt = widgetToData_(pos);
						emit sendCursorStatus( pnt[0], -1.0, pnt[1]);				
					}
					
					selected_peak_ = max_peak;
					update();
				}
				else if (e->buttons() & Qt::LeftButton) //projection
				{
					if (e->modifiers() & Qt::ShiftModifier)
					{
						rubber_band_.setGeometry(QRect(QPoint(0, last_mouse_pos_.y()), QPoint(width(), pos.y())));
					}
					else if (e->modifiers() & Qt::ControlModifier)
					{
						rubber_band_.setGeometry(QRect(QPoint(last_mouse_pos_.x(), 0), QPoint(pos.x(), height())));
					}
					else
					{
						rubber_band_.setGeometry(QRect(last_mouse_pos_, pos));
					}
					update();
				}
				break;
			}
			case AM_MEASURE:
			{
				// highlight nearest peak
				if (e->buttons() == Qt::NoButton)
				{
					DFeature<2>* max_peak = findNearestPeak_(pos);
					
					if (max_peak && max_peak != selected_peak_ && !measurement_start_)
					{
						//show Peak Coordinates
						emit sendCursorStatus(max_peak->getPosition()[RT], max_peak->getIntensity(), max_peak->getPosition()[MZ]);
						string meta = max_peak->getMetaValue(3).toString();
						if (meta!="")
							sendStatusMessage(meta, 0);
					}
					
					selected_peak_ = max_peak;
					update();
				}
				else if (e->buttons() & Qt::LeftButton && measurement_start_)
				{
					measurement_stop_ = findNearestPeak_(pos);
					last_mouse_pos_ = pos;
					update();
				
					if (measurement_stop_)
					{
						emit sendCursorStatus(measurement_stop_->getPosition()[RT] - measurement_start_->getPosition()[RT],
						                      measurement_stop_->getIntensity() / measurement_start_->getIntensity(),
						                      measurement_stop_->getPosition()[MZ] - measurement_start_->getPosition()[MZ]);
					}
					else
					{
						emit sendCursorStatus(measurement_start_->getPosition()[RT], measurement_start_->getIntensity(), measurement_start_->getPosition()[MZ]);
					}
				}
				break;
			}
			case AM_ZOOM:
			{
				//show Peak Coordinates
				PointType pnt = widgetToData_(pos);
				emit sendCursorStatus( pnt[RT], -1.0, pnt[MZ]);
				
				if (e->buttons() & Qt::LeftButton)
				{
					rubber_band_.setGeometry(last_mouse_pos_.x(), last_mouse_pos_.y(), pos.x() - last_mouse_pos_.x(), pos.y() - last_mouse_pos_.y());
					update();
				}
				break;
			}
			case AM_TRANSLATE:
			{
				if (e->buttons() & Qt::LeftButton)
				{
					//caldulate data coordinates of shift
					PointType old_data = widgetToData_(last_mouse_pos_);
					PointType new_data = widgetToData_(pos);
					//calculate x shift
					double shift = old_data.X() - new_data.X();
					double newLoX = visible_area_.minX() + shift;
					double newHiX = visible_area_.maxX() + shift;
					// check if we are falling out of bounds
					if (newLoX < overall_data_range_.minX())
					{
						newLoX = overall_data_range_.minX();
						newHiX = newLoX + visible_area_.width();
					}
					if (newHiX > overall_data_range_.maxX())
					{
						newHiX = overall_data_range_.maxX();
						newLoX = newHiX - visible_area_.width();
					}
					//calculate y shift
					shift = old_data.Y() - new_data.Y();
					double newLoY = visible_area_.minY() + shift;
					double newHiY = visible_area_.maxY() + shift;
					// check if we are falling out of bounds
					if (newLoY < overall_data_range_.minY())
					{
						newLoY = overall_data_range_.minY();
						newHiY = newLoY + visible_area_.height();
					}
					if (newHiY > overall_data_range_.maxY())
					{
						newHiY = overall_data_range_.maxY();
						newLoY = newHiY - visible_area_.height();
					}
	     		
	     		//change area
					changeVisibleArea_(AreaType(newLoX,newLoY,newHiX,newHiY));
	
					last_mouse_pos_ = pos;
				}
				break;
			}
			default:
				break;
		}
		e->accept();
	}
	
	void Spectrum2DCanvas::mouseReleaseEvent(QMouseEvent* e)
	{
		QPoint pos = e->pos();

		switch (action_mode_)
		{
			case AM_SELECT:
			{
				rubber_band_.hide();
				if (e->button() == Qt::LeftButton && getCurrentLayer().type==LayerData::DT_PEAK)
				{
					//determine data coordiantes
					QRect rect = rubber_band_.geometry();	
					AreaType area(widgetToData_(rect.topLeft()), widgetToData_(rect.bottomRight()));
					createProjections_(area, e->buttons() & Qt::ShiftModifier, e->buttons() & Qt::ControlModifier);
				}
				break;
			}
			case AM_MEASURE:
			{
				if (e->button() == Qt::LeftButton)
				{
					if (!measurement_stop_)
					{
						delete(measurement_start_);
						measurement_start_ = 0;
					}
					else
					{
						measurement_stop_ = new DFeature<2>(*measurement_stop_);
					}
					
					update();
					
					if (measurement_start_)
					{
						emit sendStatusMessage(QString("Measured: dRT = %1, dMZ = %3, Intensity ratio = %2")
																	.arg(measurement_stop_->getPosition()[MZ] - measurement_start_->getPosition()[MZ])
																	.arg(measurement_stop_->getIntensity() / measurement_start_->getIntensity())
																	.arg(measurement_stop_->getPosition()[RT] - measurement_start_->getPosition()[RT]).toAscii().data(), 0);
					}
				}
				break;
			}
			case AM_ZOOM:
			{
				rubber_band_.hide();
				if (e->button() == Qt::LeftButton) //zoom to rubber band area
				{
					QRect rect = rubber_band_.geometry();
					if (rect.width()!=0 && rect.height()!=0) //probably double click -> mouseDoubleClickEvent
					{
						AreaType area(widgetToData_(rect.topLeft()), widgetToData_(rect.bottomRight()));
						changeVisibleArea_(area, true);
					}
				}
				break;
			}
			case AM_TRANSLATE:
			{
	      setCursor(cursor_translate_);
	      break;
			}
			default:
				break;
		}
		e->accept();
	}

	void Spectrum2DCanvas::mouseDoubleClickEvent(QMouseEvent* e)
	{
		// left-doubleclick shows the whole spectrum
		if (e->button() == Qt::LeftButton && action_mode_ == AM_ZOOM)
		{
			resetZoom();
		}
	}

	void Spectrum2DCanvas::wheelEvent(QWheelEvent* e)
	{
		switch (action_mode_)
		{
			case AM_ZOOM:
				if (e->delta() > 0) // forward rotation -> zoom in
				{
					PointType new_pos = visible_area_.center();

					// adjust new width (we don't want it bigger than overall_data_range_)
					float new_width = visible_area_.width() * 0.9;
					float new_height = visible_area_.height() * 0.9;	 
					if (new_width >= overall_data_range_.width()) new_width = overall_data_range_.width();
					if (new_height >= overall_data_range_.height()) new_height = overall_data_range_.height();
			
					float half_width = new_width / 2.0f;
					float half_height = new_height / 2.0f;	
					if (new_pos.X() < overall_data_range_.minX() + half_width)   new_pos.setX(overall_data_range_.minX() + half_width);
					if (new_pos.Y() < overall_data_range_.minY() + half_height)  new_pos.setY(overall_data_range_.minY() + half_height);
					if (new_pos.X() > overall_data_range_.maxX() - half_width)   new_pos.setX(overall_data_range_.maxX() - half_width);
					if (new_pos.Y() > overall_data_range_.maxY() - half_height)  new_pos.setY(overall_data_range_.maxY() - half_height);
			
					// set visible area accordingly and redraw
					changeVisibleArea_(AreaType(new_pos.X() - half_width, new_pos.Y() - half_height, new_pos.X() + half_width, new_pos.Y() + half_height), true);
				}
				else // backward rotation -> zoom out
				{
					zoomBack_();
				}
				break;
			default:
				e->ignore();
		}
		e->accept();
	}

} //namespace

