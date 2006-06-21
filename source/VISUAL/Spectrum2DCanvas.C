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
// $Id: Spectrum2DCanvas.C,v 1.55 2006/06/09 22:00:08 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS
#include <OpenMS/VISUAL/Spectrum2DWidget.h>
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/KERNEL/DFeature.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum2DCanvasPDP.h>

//STL
#include <algorithm>	

//QT
#include <qimage.h>
#include <qpainter.h>


using namespace std;

namespace OpenMS
{
	using namespace Internal;

	Spectrum2DCanvas::Spectrum2DCanvas(QWidget* parent, const char* name)
		: SpectrumCanvas(parent, name),
		trees_(),
		marching_squares_matrices_(),
		max_values_(),
		show_contours_(),
		show_colors_(),
		show_points_(),
	  intensity_scaled_dots_(false),
		nearest_peak_(0),
		measurement_start_(0),
		measurement_stop_(0),
		tmp_peak_(),
		dot_gradient_()
	{
		// prevents errors caused by too small width,height values
		setMinimumSize(200,200);
		
		setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
		viewport()->setMouseTracking(true);
		action_mode_ = AM_SELECT;
//		min_x_ = numeric_limits<float>::max();
//		max_x_ = -1 * numeric_limits<float>::max();
//		min_y_ = numeric_limits<float>::max();
//		max_y_ = -1 * numeric_limits<float>::max();
	}
	
	Spectrum2DCanvas::~Spectrum2DCanvas()
	{
		for (UnsignedInt i=0; i<trees_.size(); i++)
		{
			if (trees_[i]) delete trees_[i];
		}
	}
	
	void Spectrum2DCanvas::print(QPainter* p, int width, int height)
	{
		for (UnsignedInt i=0; i<trees_.size(); i++)
		{
			if (layer_visible_[i])
			{
				paintContent_(i, p, width, height);
			}
		}
	}
	
	void Spectrum2DCanvas::showContours(bool on)
	{
		recalculate_ = show_contours_[current_data_] != on;
		show_contours_[current_data_] = on;
		invalidate_();
	}
	
	void Spectrum2DCanvas::showColors(bool on)
	{
		recalculate_ = show_colors_[current_data_] != on;
		show_colors_[current_data_] = on;
		invalidate_();
	}
	
	void Spectrum2DCanvas::showPoints(bool on)
	{
		recalculate_ = show_points_[current_data_] != on;
		show_points_[current_data_] = on;
		invalidate_();
	}
	
	void Spectrum2DCanvas::changeShowContours()
	{
		recalculate_ = true;
		show_contours_[current_data_] = !(show_contours_[current_data_]);
		invalidate_();
	}
	
	void Spectrum2DCanvas::changeShowColors()
	{
		recalculate_ = true;
		show_colors_[current_data_] = !(show_colors_[current_data_]);
		invalidate_();
	}
	
	void Spectrum2DCanvas::changeShowPoints()
	{
		recalculate_ = true;
		show_points_[current_data_] = !(show_points_[current_data_]);
		invalidate_();
	}
	
	bool Spectrum2DCanvas::getShowContours()
	{
		return show_contours_[current_data_];
	}
	
	bool Spectrum2DCanvas::getShowColors()
	{
		return show_colors_[current_data_];
	}
	
	bool Spectrum2DCanvas::getShowPoints()
	{
		return show_points_[current_data_];
	}
	
	void Spectrum2DCanvas::setIntensityScaledDots(bool on)
	{
	  intensity_scaled_dots_ = on;
	  invalidate_();
	}
	
	void Spectrum2DCanvas::contentsMousePressEvent(QMouseEvent* e)
	{
		// when pressing down the mouse, only the position
		// has to be stored, since with no tool any action
		// is performed when the mouse button is pressed
		// but the mouse is not (yet) moved
		mouse_pos_ = contentsToViewport(e->pos());
	  if (e->button() == LeftButton)
		{
			if (action_mode_ == AM_TRANSLATE)
			{
				viewport()->setCursor(cursor_translate_in_progress_);
			}
			else if (action_mode_ == AM_MEASURE)
			{
				DPeak<2>* tmp = nearest_peak_;
				if (tmp)
				{
					delete(measurement_start_);
					measurement_start_ = new DPeak<2>(*tmp);
				}
				else
				{
					delete(measurement_start_);
					measurement_start_ = 0;
				}
				delete(measurement_stop_);
				measurement_stop_ = 0;
			}
		}
		e->accept();
	}
	
	void Spectrum2DCanvas::contentsMouseReleaseEvent(QMouseEvent* e)
	{
		QPoint pos = contentsToViewport(e->pos());
		
		if (e->button() == Qt::RightButton)
		{
			// context menu
			emit contextMenu(contentsToViewport(e->globalPos()));
			return;
		}
		
		switch (action_mode_)
		{
			default:
			case AM_SELECT:
			{
				if (e->button() == Qt::LeftButton)
				{
					// clear rubber band
					AreaType area(widgetToChart_(mouse_pos_), widgetToChart_(pos));
					//area.normalize();
					
					if (e->state() & Qt::ShiftButton)
					{
						createHorzScan_(area.minY(), area.maxY());
					}
					else if (e->state() & Qt::ControlButton)
					{
						createVertScan_(area.minX(), area.maxX());
					}
					else
					{
						createHorzScan_(area.minY(), area.maxY());
						createVertScan_(area.minX(), area.maxX());
					}
					refresh_();
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
						measurement_stop_ = new DPeak<2>(*measurement_stop_);
					}
					
					refresh_();
					
					if (measurement_start_)
					{
						emit sendStatusMessage(QString("Measured: dRT = %3, dMZ = %1, dInt = %2")
																	.arg(measurement_stop_->getPosition()[MZ] - measurement_start_->getPosition()[MZ])
																	.arg(measurement_stop_->getIntensity() - measurement_start_->getIntensity())
																	.arg(measurement_stop_->getPosition()[RT] - measurement_start_->getPosition()[RT]).ascii(), 0);
					}
				}
				break;
			}
			case AM_ZOOM:
			{
				if (mouse_pos_ == pos)
				{
					if (e->button() == Qt::LeftButton)
					{
						// left button means zoom in
						zoomIn_(widgetToChart_(pos));
					}
					else if (e->button() == Qt::MidButton)
					{
						// middle button means zoom out
						zoomBack_();
					}
				}
				else
				{
					// we get here, if the user has been dragging a
					// rectangular area in the diagram. This will
					// the be the whole visible area, and this always
					// means that we zoom in.
	
					if (e->button() == Qt::LeftButton)
					{
						// sort coordinates ascending
						int min_x = mouse_pos_.x();
						int max_x = pos.x();
						int min_y = mouse_pos_.y();
						int max_y = pos.y();
						if (min_x > max_x) swap(min_x, max_x);
						if (min_y > max_y) swap(min_y, max_y);
	
						// transform widget coordinates to chart coordinates
						PointType left_top = widgetToChart_(QPoint(min_x, min_y));
						PointType right_bottom = widgetToChart_(QPoint(max_x, max_y));
						
						// eventually adjust values to the QuadTree area's
						// borders
						if (left_top.X() < trees_[current_data_]->getArea().minX()) 
							left_top.setX(trees_[current_data_]->getArea().minX());
						if (left_top.Y() < trees_[current_data_]->getArea().minY()) 
							left_top.setY(trees_[current_data_]->getArea().minY());
						if (right_bottom.X() > trees_[current_data_]->getArea().maxX()) 
							right_bottom.setX(trees_[current_data_]->getArea().maxX());
						if (right_bottom.Y() > trees_[current_data_]->getArea().maxY()) 
							right_bottom.setY(trees_[current_data_]->getArea().maxY());
					 
						//emit visibleAreaChanged(visible_area_);
						changeVisibleArea_(AreaType(left_top, right_bottom));
					}
				}
				break;
			}
			case AM_TRANSLATE:
			{
	      viewport()->setCursor(cursor_translate_);  // open-hand cursor
				// do nothing, because releasing the mouse while
				// moving only means that we're done moving, and
				// so we don't have to do anything here.
				break;
			}
		}
	
		e->accept();
	}
	
	void Spectrum2DCanvas::highlightPeaks_()
	{
		QPainter p(viewport());
		
		if (measurement_start_)
		{
			p.setPen(Qt::black);
			
			QPoint line_end;
			
			if (measurement_stop_)
			{
				line_end = chartToWidget_(measurement_stop_->getPosition());
				//cout << "Line end: " << line_end << endl;
			}
			else
			{
				line_end = mouse_pos_;
				//cout << "Ende: " << line_end.x() << " " << line_end.y() << endl;
			}
			
			p.drawLine(chartToWidget_(measurement_start_->getPosition()), line_end);
		}
		
		highlightPeak_(&p, nearest_peak_);
		highlightPeak_(&p, measurement_start_);
		highlightPeak_(&p, measurement_stop_);
	}
	
	void Spectrum2DCanvas::highlightPeak_(QPainter* p, DPeak<2>* peak)
	{
		if (!peak)
			return;
		
		const QPoint diff(5, 5);
		p->setPen(QPen(Qt::red, 2));
		
		QPoint peak_pos(chartToWidget_(peak->getPosition()));
		QRect peak_rect(peak_pos - diff, peak_pos + diff);
		
		//cout << "Highlight: " << peak_rect.x() << " " << peak_rect.y() << endl;
		
		p->drawEllipse(peak_rect);
	}
	
	DPeak<2>* Spectrum2DCanvas::findNearestPeak_(QPoint pos)
	{
		const QPoint diff(5, 5);
		QRect rect(pos - diff, pos + diff);
		AreaType area(widgetToChart_(rect.topLeft()), widgetToChart_(rect.bottomRight()));
				
		DPeak<2>* max_peak = 0;
		float max_int = -1 * numeric_limits<float>::max();
		
		for (QuadTreeType_::Iterator i = trees_[current_data_]->begin(area); i != trees_[current_data_]->end(); ++i)
		{
			//cout << "second: " << i->second->getPosition()[RT] << " " << i->second->getPosition()[MZ] << endl;
			if (i->second->getIntensity() > max_int)
			{
				max_int = i->second->getIntensity();
				
				tmp_peak_.setIntensity(i->second->getIntensity());
				tmp_peak_.getPosition()[0] = i->second->getPosition()[0];
				tmp_peak_.getPosition()[1] = i->first[1];
				
				//cout << "Peak MZ: " << tmp_peak_.getPosition()[MZ] << " RT: " << tmp_peak_.getPosition()[RT] << endl;
				
				max_peak = &tmp_peak_;
			}
		}
		
		return max_peak;
	}
	
	void Spectrum2DCanvas::contentsMouseMoveEvent(QMouseEvent* e)
	{
		QPoint pos = contentsToViewport(e->pos());
		
		switch (action_mode_)
		{
			default:
			case AM_SELECT:
			{
	      viewport()->setCursor(Qt::ArrowCursor);
				// highlight nearest peak
				if (e->state() == Qt::NoButton)
				{
					DPeak<2>* max_peak = findNearestPeak_(pos);
					
					if (max_peak)
					{
						// show Peak Coordinates (with intensity)
						emit sendCursorStatus(max_peak->getPosition()[0], max_peak->getIntensity(), max_peak->getPosition()[1]);
						//show lable
						string meta = max_peak->getMetaValue(3).toString();
						if (meta!="") sendStatusMessage(meta, 0);
					}
					else if (!max_peak)
					{
						//show Peak Coordinates (without intensity)
						PointType pnt = widgetToChart_(pos);
						emit sendCursorStatus( pnt[0], -1.0, pnt[1]);				
					}
					
					nearest_peak_ = max_peak;
					refresh_();
				}
				else if (e->state() & Qt::LeftButton)
				{
					// select 1D spectrum
					QRect rect_horz(QPoint(0, mouse_pos_.y()), QPoint(viewport()->width(), pos.y()));
					QRect rect_vert(QPoint(mouse_pos_.x(), 0), QPoint(pos.x(), viewport()->height()));
					QRect rect_mid(mouse_pos_, pos);
	
					// draw rubber band(s)
					refresh_();
					QPainter p(viewport());
					p.setPen(Qt::NoPen);
					p.setBrush(Qt::red);
					p.setRasterOp(Qt::XorROP);
	
					if (e->state() & Qt::ShiftButton)
					{
						p.drawRect(rect_horz);
					}
					else if (e->state() & Qt::ControlButton)
					{
						p.drawRect(rect_vert);
					}
					else
					{
						p.drawRect(rect_horz);
						p.drawRect(rect_vert);
						p.drawRect(rect_mid);
					}
				}
				break;
			}
			case AM_MEASURE:
			{
	      viewport()->setCursor(Qt::ArrowCursor);
				// highlight nearest peak
				if (e->state() == Qt::NoButton)
				{
					DPeak<2>* max_peak = findNearestPeak_(pos);
					
					if (max_peak && max_peak != nearest_peak_ && !measurement_start_)
					{
						//show Peak Coordinates (with intensity)
						emit sendCursorStatus(max_peak->getPosition()[MZ], max_peak->getIntensity(), max_peak->getPosition()[RT]);
						string meta = max_peak->getMetaValue(3).toString();
						if (meta!="")
							sendStatusMessage(meta, 0);
					}
					
					nearest_peak_ = max_peak;
					refresh_();
				}
				else if (e->state() & Qt::LeftButton && measurement_start_)
				{
					measurement_stop_ = findNearestPeak_(pos);
					mouse_pos_ = pos;
					refresh_();
				
					if (measurement_stop_)
					{
						emit sendCursorStatus(measurement_stop_->getPosition()[MZ] - measurement_start_->getPosition()[MZ],
						                      measurement_stop_->getIntensity() - measurement_start_->getIntensity(),
						                      measurement_stop_->getPosition()[RT] - measurement_start_->getPosition()[RT]);
					}
					else
					{
						emit sendCursorStatus(measurement_start_->getPosition()[MZ], measurement_start_->getIntensity(), measurement_start_->getPosition()[RT]);
					}
				}
				break;
			}
			case AM_ZOOM:
			{
				//show Peak Coordinates
				PointType pnt = widgetToChart_(pos);
				emit sendCursorStatus( pnt[MZ], -1.0, pnt[RT]);
	      //set Cursor
	      viewport()->setCursor(Qt::CrossCursor);
				
				if (e->state() & Qt::LeftButton)
				{
					// draw zoom rect
					refresh_();
					QPainter p(viewport());
					p.setBrush(Qt::red);
					p.setRasterOp(Qt::XorROP);
					p.drawRect(mouse_pos_.x(), mouse_pos_.y(), pos.x() - mouse_pos_.x(), pos.y() - mouse_pos_.y());
				}
				break;
			}
			case AM_TRANSLATE:
			{
				//show Peak Coordinates
				PointType pnt = widgetToChart_(pos);
				emit sendCursorStatus( pnt[MZ], -1.0, pnt[RT]);
				//set Cursor
	      viewport()->setCursor(cursor_translate_);
				
				if (e->state() & Qt::LeftButton)
				{
	        viewport()->setCursor(cursor_translate_in_progress_);
					// move the visible area
					
					QPoint pmove = mouse_pos_ - pos;
	
					PointType move = widgetToChart_(pmove);
	
					PointType left_top = PointType(move.X(), move.Y());
	
					// eventually adjust values to the QuadTree area's borders
					if (left_top.X() < trees_[current_data_]->getArea().minX()) 
						left_top.setX(trees_[current_data_]->getArea().minX());
					if (left_top.Y() < trees_[current_data_]->getArea().minY()) 
						left_top.setY(trees_[current_data_]->getArea().minY());
					if (left_top.X() + visible_area_.width() > trees_[current_data_]->getArea().maxX()) 
						left_top.setX(trees_[current_data_]->getArea().maxX() - visible_area_.width());
					if (left_top.Y() + visible_area_.height() > trees_[current_data_]->getArea().maxY()) 
						left_top.setY(trees_[current_data_]->getArea().maxY() - visible_area_.height());
	
					PointType right_bottom = PointType(left_top.X() + visible_area_.width(), left_top.Y() + visible_area_.height());
	
					changeVisibleArea_(AreaType(left_top, right_bottom));
	
					mouse_pos_ = pos;
				}
				if (e->button() == Qt::RightButton)
				{
					// context menu
					emit contextMenu(contentsToViewport(e->globalPos()));
				}
				break;
			}
		}
	
		e->accept();
	}
	
	void Spectrum2DCanvas::contentsWheelEvent(QWheelEvent* e)
	{
	  QPoint pos = contentsToViewport(e->pos());
	
		switch (action_mode_)
		{
		case AM_ZOOM:
			if (e->delta() > 0)
			{
				// forward rotation -> zoom in
				zoomIn_(widgetToChart_(pos));
			}
	
			else
			{
				// backward rotation -> zoom out
				zoomOut_(widgetToChart_(pos));
			}
	
			e->accept();
			break;
	
		default:
			e->ignore();
		}
	}
	
	float Spectrum2DCanvas::betweenFactor_(float v1, float v2, float val)
	{
		float lo = min(v1, v2);
		float hi = max(v1, v2);
		return (hi - lo == 0) ? 1 : (val - lo) / (hi - lo);
	}
	
	const QColor& Spectrum2DCanvas::heightColor_(float val)
	{
		if (intensity_modification_ == IM_LOG)
		{
			return surface_gradient_.precalculatedColorAt(log(val+1)); //prevent log of numbers samller than 1
		}
		else
		{
			return surface_gradient_.precalculatedColorAt(val);
		}
	}
	
	Spectrum2DCanvas::AreaType Spectrum2DCanvas::getLeftTopCell_(UnsignedInt data_set)
	{
		PointType cell_size(trees_[data_set]->getArea().width() / getPrefAsInt("Preferences:2D:MarchingSquaresSteps"), trees_[data_set]->getArea().height() / getPrefAsInt("Preferences:2D:MarchingSquaresSteps"));
		PointType pos(float(int(visible_area_.minX() / cell_size.X())) * cell_size.X(), float(int(visible_area_.minY() / cell_size.Y())) * cell_size.Y());
		
		return AreaType(pos.X(), pos.Y(), pos.X() + cell_size.X(), pos.Y() + cell_size.Y() /*/ 10.0f*/);
	}
	
	void Spectrum2DCanvas::getMarchingSquareMatrix_(UnsignedInt data_set)
	{
		//cout << "Marching squares matrix for dataset " << data_set <<":" <<endl;
		AreaType cell = getLeftTopCell_(data_set);
	
		const float cell_width = cell.width();
		const float cell_height = cell.height();
		const float half_width = cell_width / 2;
		const float half_height = cell_height / 2;
	
		float y = cell.minY();
	
		int i, j;
		for (i = 0; y <= visible_area_.maxY() + cell_height; i++)
		{
			float x = cell.minX();
			vector<float> line;
			for (j = 0; x <= visible_area_.maxX() + cell_width; j++)
			{
				// build sum of all peak heights in the current cell
				float sum = 0.0f;
				AreaType area(x - half_width, y - half_height, x + half_width, y + half_height);
				for (QuadTreeType_::Iterator i = trees_[data_set]->begin(area); i != trees_[data_set]->end(); ++i)
				{
					sum += i->second->getIntensity();
				}
				// log mode
				if (intensity_modification_ == IM_LOG)
				{
					sum += log(sum+1); // prevent log of numbers smaller than one
				}
				// store max
				if (sum > max_values_[data_set])
				{
					max_values_[data_set] = sum;
				}
				
				line.push_back(sum);
				//cout << sum << " ";
				x += cell_width;
			}
			//cout << endl;
			marching_squares_matrices_[data_set].push_back(line);
			y += cell_height;
		}
	}
	
	void Spectrum2DCanvas::paintContent_(UnsignedInt data_set, QPainter* p, int width, int height)
	{
		emit sendStatusMessage("repainting", 0);
		if (trees_[data_set] && !(trees_[data_set]->begin(trees_[data_set]->getArea()) == trees_[data_set]->end()))
		{
			if (show_colors_[data_set])
			{
				paintColorMap_(data_set, p, width, height);
			}
	
			if (show_contours_[data_set])
			{
				paintContourLines_(data_set, p, width, height);
			}
	
			if (show_points_[data_set])
			{
				paintPoints_(data_set, p,width,height);
			}
		}
		emit sendStatusMessage("", 0);
	}
	
	void Spectrum2DCanvas::paintPoints_(UnsignedInt data_set, QPainter* p, int width, int height)
	{
		p->setPen(Qt::black);
		p->setBrush(Qt::black);
		QColor inter = Qt::black;
	
		// TODO: ConstIterator
		bool isFeature = getDataSet(data_set).metaValueExists("FeatureDrawMode");
		for (QuadTreeType_::SortedIterator i = trees_[data_set]->sortedBegin(visible_area_);
				 i != trees_[data_set]->sortedEnd(); ++i)
		{
			if (!isFeature)
			{
				QPoint pos = chartToContext_(i->first, width, height);
				if (getDotMode()==DOT_GRADIENT)
				{
					if (intensity_modification_ == IM_LOG)
					{
						//prevent log of numbers samller than 1
						inter = dot_gradient_.precalculatedColorAt(log(i->second->getIntensity()+1)); 
					}
					else
					{
						inter = dot_gradient_.precalculatedColorAt(i->second->getIntensity());
						//cout << "Int: " << i->second->getIntensity() << " -> " << inter.name() << endl;
					}
					p->setPen(inter);
					p->setBrush(inter);
				}
				
				if (intensity_scaled_dots_)
				{
					// points get scaled relative to the minimum displayed intensity
					int radius = static_cast<int>(log10(i->second->getIntensity() - overall_data_range_.min()[2])/2);
					p->drawEllipse(pos.x()- radius, pos.y() - radius, 2*radius, 2*radius);  
				}  else
				{
					p->drawEllipse(pos.x() - 2, pos.y() - 2, 4, 4);
				}
			}
			else //Draw special feature attributes like convex hull
			{
				const DFeature<2>* feature = dynamic_cast<const DFeature<2>* >(&(*(i->second)));
				if (feature)
				{
					p->setBrush(Qt::NoBrush);
					String mode = getDataSet(data_set).getMetaValue("FeatureDrawMode").toString();
	
					if (mode=="ConvexHulls")
						for (UnsignedInt hull=0; hull<feature->getConvexHulls().size(); ++hull) // Draw the convex hulls
						{
							QPointArray points(feature->getConvexHulls()[hull].size());
							UnsignedInt index=0;
							for(DFeature<2>::ConvexHullType::const_iterator it=feature->getConvexHulls()[hull].begin(); 
									it!=feature->getConvexHulls()[hull].end(); ++it, ++index)
							{
								QPoint p = chartToContext_(*it, width, height);
								points.setPoint(index, p.x(), p.y());
							}	
							p->drawPolygon(points);
						}
				}
			}
		}                
	}
	
	void Spectrum2DCanvas::paintContourLines_(UnsignedInt data_set, QPainter* p, int width, int height)
	{
		const float height_difference = max_values_[data_set] / 10.0f; //TODO
	
		if (height_difference == 0) return;
	
	
		p->setPen(Qt::black);
	
		AreaType cell = getLeftTopCell_(data_set);
	
		const float cell_width = cell.width();
		const float cell_height = cell.height();
	
		QPoint cell_size = chartToContext_(PointType(visible_area_.minX() + cell_width, visible_area_.minY() + cell_height), width, height);
	
		// draw the lines
		float y = cell.minY();
		for (int i = 0; y <= visible_area_.maxY(); i++)
		{
			float x = cell.minX();
			for (int j = 0; x <= visible_area_.maxX(); j++)
			{

				float left_top = marching_squares_matrices_[data_set][i][j];
				float right_top = 0;
				float left_bottom = 0;
				float right_bottom = marching_squares_matrices_[data_set][i + 1][j + 1];

				if ( getMappingInfo().isMzToXAxis() )
				{
				  right_top = marching_squares_matrices_[data_set][i][j + 1];
				  left_bottom = marching_squares_matrices_[data_set][i + 1][j];
				}
				else
				{
				  left_bottom = marching_squares_matrices_[data_set][i][j + 1];
				  right_top = marching_squares_matrices_[data_set][i + 1][j];		
				}
	
				const float minimum = min(left_top, min(right_top, min(left_bottom, right_bottom)));
				const float maximum = max(left_top, max(right_top, max(left_bottom, right_bottom)));
	
				const float first = static_cast<float>(static_cast<int>(minimum / height_difference)) * height_difference;
	
				QPoint cell_pos = chartToContext_(PointType(x, y), width, height);
	
				for (float height = first; height <= maximum; height += height_difference)
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
							p->drawLine(cell_pos.x() + int(betweenFactor_(left_bottom, right_bottom, height) * cell_size.x()),
													cell_pos.y() + cell_size.y(),
													cell_pos.x() + cell_size.x(),
													cell_pos.y() + int(betweenFactor_(right_top, right_bottom, height) * cell_size.y()));
							break;
						}
						case 14:
						{
							p->drawLine(cell_pos.x() + cell_size.x() - int(betweenFactor_(left_bottom, right_bottom, height) * cell_size.x()),
													cell_pos.y() + cell_size.y(),
													cell_pos.x() + cell_size.x(),
													cell_pos.y() + cell_size.y() - int(betweenFactor_(right_top, right_bottom, height) * cell_size.y()));
							break;
						}
						case 2:
						{
							p->drawLine(cell_pos.x(),
													cell_pos.y() + int(betweenFactor_(left_top, left_bottom, height) * cell_size.y()),
													cell_pos.x() + cell_size.x() - int(betweenFactor_(left_bottom, right_bottom, height) * cell_size.x()),
													cell_pos.y() + cell_size.y());
							break;
						}
						case 13:
						{
							p->drawLine(cell_pos.x(),
													cell_pos.y() + cell_size.y() - int(betweenFactor_(left_top, left_bottom, height) * cell_size.y()),
													cell_pos.x() + int(betweenFactor_(left_bottom, right_bottom, height) * cell_size.x()),
	
													cell_pos.y() + cell_size.y());
							break;
						}
						case 3:
						{
							p->drawLine(cell_pos.x(),
													cell_pos.y() + int(betweenFactor_(left_top, left_bottom, height) * cell_size.y()),
													cell_pos.x() + cell_size.x(),
													cell_pos.y() + int(betweenFactor_(right_top, right_bottom, height) * cell_size.y()));
	
							break;
						}
						case 12:
						{
							p->drawLine(cell_pos.x(),
	
													cell_pos.y() + cell_size.y() - int(betweenFactor_(left_top, left_bottom, height) * cell_size.y()),
													cell_pos.x() + cell_size.x(),
													cell_pos.y() + cell_size.y() - int(betweenFactor_(right_top, right_bottom, height) * cell_size.y()));
							break;
						}
						case 4:
						{
							p->drawLine(cell_pos.x() + int(betweenFactor_(left_top, right_top, height) * cell_size.x()),
													cell_pos.y(),
													cell_pos.x() + cell_size.x(),
													cell_pos.y() + cell_size.y() - int(betweenFactor_(right_top, right_bottom, height) * cell_size.y()));
							break;
						}
						case 11:
						{
							p->drawLine(cell_pos.x() + cell_size.x() - int(betweenFactor_(left_top, right_top, height) * cell_size.x()),
													cell_pos.y(),
													cell_pos.x() + cell_size.x(),
													cell_pos.y() + int(betweenFactor_(right_top, right_bottom, height) * cell_size.y()));
							break;
						}
						case 5:
						{
							p->drawLine(cell_pos.x() + int(betweenFactor_(left_top, right_top, height) * cell_size.x()),
													cell_pos.y(),
													cell_pos.x() + int(betweenFactor_(left_bottom, right_bottom, height) * cell_size.x()),
													cell_pos.y() + cell_size.y());
							break;
						}
						case 10:
						{
							p->drawLine(cell_pos.x() + cell_size.x() - int(betweenFactor_(left_top, right_top, height) * cell_size.x()),
													cell_pos.y(),
													cell_pos.x() + cell_size.x() - int(betweenFactor_(left_bottom, right_bottom, height) * cell_size.x()),
													cell_pos.y() + cell_size.y());
							break;
						}
						case 6:
						{
							p->drawLine(cell_pos.x(),
													cell_pos.y() + int(betweenFactor_(left_top, left_bottom, height) * cell_size.y()),
													cell_pos.x() + cell_size.x() - int(betweenFactor_(left_bottom, right_bottom, height) * cell_size.x()),
													cell_pos.y() + cell_size.y());
							p->drawLine(cell_pos.x() + int(betweenFactor_(left_top, right_top, height) * cell_size.x()),
													cell_pos.y(),
													cell_pos.x() + cell_size.x(),
													cell_pos.y() + cell_size.y() - int(betweenFactor_(right_top, right_bottom, height) * cell_size.y()));
							break;
						}
						case 9:
						{
							p->drawLine(cell_pos.x() + int(betweenFactor_(left_bottom, right_bottom, height) * cell_size.x()),
													cell_pos.y() + cell_size.y(),
													cell_pos.x() + cell_size.x(),
													cell_pos.y() + int(betweenFactor_(right_top, right_bottom, height) * cell_size.y()));
							p->drawLine(cell_pos.x(),
													cell_pos.y() + cell_size.y() - int(betweenFactor_(left_top, left_bottom, height) * cell_size.y()),
													cell_pos.x() + cell_size.x() - int(betweenFactor_(left_top, right_top, height) * cell_size.x()),
													cell_pos.y());
							break;
						}
	
						case 7:
						{
							p->drawLine(cell_pos.x(),
													cell_pos.y() + int(betweenFactor_(left_top, left_bottom, height) * cell_size.y()),
													cell_pos.x() + int(betweenFactor_(left_top, right_top, height) * cell_size.x()),
													cell_pos.y());
							break;
						}
						case 8:
						{
							p->drawLine(cell_pos.x(),
													cell_pos.y() + cell_size.y() - int(betweenFactor_(left_top, left_bottom, height) * cell_size.y()),
													cell_pos.x() + cell_size.x() - int(betweenFactor_(left_top, right_top, height) * cell_size.x()),
													cell_pos.y());
	
							break;
						}
					}
				}
				x += cell_width;
			}
			y += cell_height;
		}
	}
	
	void Spectrum2DCanvas::paintColorMap_(UnsignedInt data_set, QPainter* p, int width, int height)
	{
		QImage image(buffer_->width(), buffer_->height(), 32);
		QRgb* image_start = reinterpret_cast<QRgb*>(image.scanLine(0));
		const uint image_line_diff = reinterpret_cast<QRgb*>(image.scanLine(1)) - image_start;
	
		AreaType cell = getLeftTopCell_(data_set);
		
		const float cell_width = cell.width();
		const float cell_height = cell.height();
	
		QPoint cell_size = chartToContext_(PointType(visible_area_.minX() + cell_width, visible_area_.minY() + cell_height), width, height);
	
		float y = cell.minY();
		vector<vector<const QColor*> > color_matrix;
		for (int i = 0; y <= visible_area_.maxY() + cell_height; i++)
		{
			color_matrix.insert(color_matrix.end(), vector<const QColor*>() );
			float x = cell.minX();
			for (int j = 0; x <= visible_area_.maxX() + cell_width; j++)
			{
				color_matrix.back().push_back(&heightColor_(marching_squares_matrices_[data_set][i][j]));
				x += cell_width;
			}
			y += cell_height;
		}
	
		QRgb* pixel;
	
		y = cell.minY();
		for (int i = 0; y <= visible_area_.maxY(); i++)
		{
			float x = cell.minX();
			for (int j = 0; x <= visible_area_.maxX(); j++)
			{
				const QColor* left_top = color_matrix[i][j];
				const QColor* right_top = 0;
				const QColor* left_bottom = 0;
				const QColor* right_bottom = color_matrix[i + 1][j + 1];
	
				if ( getMappingInfo().isMzToXAxis() )
				{
					right_top = color_matrix[i][j + 1];
					left_bottom = color_matrix[i + 1][j];
				}
				else
				{
					left_bottom = color_matrix[i][j + 1];
					right_top = color_matrix[i + 1][j];					
				}
				QPoint cell_pos = chartToContext_(PointType(x, y), width, height);
		
				int left_red = left_top->red() << 8;
				int left_green = left_top->green() << 8;
				int left_blue = left_top->blue() << 8;
				int right_red = right_top->red() << 8;
				int right_green = right_top->green() << 8;
				int right_blue = right_top->blue() << 8;
	
				const int left_d_red = ((left_bottom->red() << 8) - left_red) / cell_size.y();
				const int left_d_green = ((left_bottom->green() << 8) - left_green) / cell_size.y();
				const int left_d_blue = ((left_bottom->blue() << 8) - left_blue) / cell_size.y();
				const int right_d_red = ((right_bottom->red() << 8) - right_red) / cell_size.y();
				const int right_d_green = ((right_bottom->green() << 8) - right_green) / cell_size.y();
				const int right_d_blue = ((right_bottom->blue() << 8) - right_blue) / cell_size.y();
	
				pixel = image_start;
				pixel += cell_pos.y() * buffer_->width() + cell_pos.x();
	
				for (int py = 0; py != cell_size.y() + 1; py++)
				{
					QRgb* start_pixel = pixel;
	
					// vertical clipping
					if (cell_pos.y() + py >= 0 && cell_pos.y() + py < buffer_->height())
					{
						const int d_red = (right_red - left_red) / cell_size.x();
						const int d_green = (right_green - left_green) / cell_size.x();
	
						const int d_blue = (right_blue - left_blue) / cell_size.x();
	
						int c_red = left_red;
						int c_green = left_green;
						int c_blue = left_blue;
	
						for (int px = 0; px != cell_size.x() + 1; px++)
						{
							// horizontal clipping
							if (cell_pos.x() + px >= 0 && cell_pos.x() + px < buffer_->width())
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
	
				x += cell_width;
			}
			y += cell_height;
		}
	
		p->drawImage(0, 0, image);
	}
	
	void Spectrum2DCanvas::refresh_()
	{
		bitBlt(viewport(), 0, 0, buffer_, 0, 0, viewport()->width(), viewport()->height(), Qt::CopyROP, true);
		
		highlightPeaks_();
	}
	
	void Spectrum2DCanvas::invalidate_()
	{
		if (recalculate_)
		{
			marching_squares_matrices_.clear();
			max_values_.clear();
	
			for (UnsignedInt i=0; i<getDataSetCount(); i++)
			{
				marching_squares_matrices_.push_back(std::vector< std::vector<float> >());
				max_values_.push_back(0.0f);
				if (layer_visible_[i] && (show_colors_[i] || show_contours_[i]))
				{
					getMarchingSquareMatrix_(i);
				}
			}
			recalculate_ = false;
		}
		buffer_->fill(QColor(getPrefAsString("Preferences:2D:BackgroundColor").c_str()));
		QPainter p(buffer_);
		p.setRasterOp(Qt::AndROP);
		print(&p, viewport()->width(), viewport()->height());
	  
		p.setRasterOp(Qt::CopyROP);
		paintGridLines_(&p);
		refresh_();
	}
	
	void Spectrum2DCanvas::zoom_(const PointType& pos, float factor)
	{
		// calculate new width
		float new_width = visible_area_.width() * factor;
		float new_height = visible_area_.height() * factor;
		PointType new_pos = pos;
	 
		// adjust new width (we don't want it bigger than the complete area covered by the QuadTree)
		if (new_width >= trees_[current_data_]->getArea().width()) new_width = trees_[current_data_]->getArea().width();
		if (new_height >= trees_[current_data_]->getArea().height()) new_height = trees_[current_data_]->getArea().height();
	
		float half_width = new_width / 2;
		float half_height = new_height / 2;
	
		// calculate new position
		if (new_pos.X() < trees_[current_data_]->getArea().minX() + half_width) 
			new_pos.setX(trees_[current_data_]->getArea().minX() + half_width);
		if (new_pos.Y() < trees_[current_data_]->getArea().minY() + half_height) 
			new_pos.setY(trees_[current_data_]->getArea().minY() + half_height);
		if (new_pos.X() > trees_[current_data_]->getArea().maxX() - half_width) 
			new_pos.setX(trees_[current_data_]->getArea().maxX() - half_width);
		if (new_pos.Y() > trees_[current_data_]->getArea().maxY() - half_height) 
			new_pos.setY(trees_[current_data_]->getArea().maxY() - half_height);
	
		// set visible area accordingly and redraw
		changeVisibleArea_(AreaType(new_pos.X() - half_width, new_pos.Y() - half_height, new_pos.X() + half_width, new_pos.Y() + half_height));
	}
	
	void Spectrum2DCanvas::zoomIn_(const PointType& /*pos*/)
	{
		float zoom_in_factor = 0.95;
		zoom_(PointType(visible_area_.center()), zoom_in_factor);
	}
	
	void Spectrum2DCanvas::zoomOut_(const PointType& /*pos*/)
	{
		float zoom_out_factor = 1.05;
		zoom_(PointType(visible_area_.center()), zoom_out_factor);
	}
	
	void Spectrum2DCanvas::intensityDistributionChange_()
	{
		emit sendStatusMessage("reconstructing quad trees", 0);
	
		QuadTreeType_* new_tree = new QuadTreeType_(visible_area_);
		for (ExperimentType::Iterator exp_it = currentDataSet().begin(); exp_it != currentDataSet().end(); ++exp_it)
		{
			if (exp_it->getMSLevel()!=1)
			{
				continue;
			}
			for (SpectrumIteratorType i = exp_it->begin(); i != exp_it->end(); ++i)
			{
				if (i->getIntensity() >= getMinDispInt() && i->getIntensity() <= getMaxDispInt())
				{
					try
					{
						new_tree->insert(PointType(i->getPosition()[0],exp_it->getRetentionTime()), &(*i));  
					}
					catch (Exception::IllegalTreeOperation& e)
					{
						// removes problems with identical peak positions
					}   
				}
			}
		}
	  
		delete trees_[current_data_];
		trees_[current_data_] = new_tree;
	
		recalculate_ = true;
		invalidate_();
		emit sendStatusMessage("", 0);
		emit visibleAreaChanged(visible_area_);
	}
	
	void Spectrum2DCanvas::intensityModificationChange_()
	{
		recalculateDotGradient_();
		recalculateSurfaceGradient_();
		SpectrumCanvas::intensityModificationChange_();
	}
	
	void Spectrum2DCanvas::recalculateDotGradient_()
	{
		//cout << "recalculateDotGradient_" << endl;
		if (intensity_modification_ == IM_LOG)
		{
			//cout << "LOG:" <<" "<< log(overall_data_range_.min()[2]) <<" "<< log(overall_data_range_.max()[2])<<" "<<getPrefAsInt("Preferences:2D:Dot:InterpolationSteps")<<endl;
			dot_gradient_.activatePrecalculationMode(log(overall_data_range_.min()[2]+1), log(overall_data_range_.max()[2]+1), getPrefAsInt("Preferences:2D:Dot:InterpolationSteps"));
		}
		else
		{
			//cout << "NORMAL:" << overall_data_range_.min()[2] <<" "<< overall_data_range_.max()[2]<<" "<<getPrefAsInt("Preferences:2D:Dot:InterpolationSteps")<<endl;
			dot_gradient_.activatePrecalculationMode(overall_data_range_.min()[2], overall_data_range_.max()[2], getPrefAsInt("Preferences:2D:Dot:InterpolationSteps"));
		}	
	}
	
	void Spectrum2DCanvas::recalculateSurfaceGradient_()
	{
		if (intensity_modification_ == IM_LOG)
		{
			surface_gradient_.activatePrecalculationMode(log(overall_data_range_.min()[2]+1), log(overall_data_range_.max()[2]+1), getPrefAsInt("Preferences:2D:Surface:InterpolationSteps"));
		}
		else
		{
			surface_gradient_.activatePrecalculationMode(overall_data_range_.min()[2], overall_data_range_.max()[2], getPrefAsInt("Preferences:2D:Surface:InterpolationSteps"));		
		}	
	}
	
	void Spectrum2DCanvas::createHorzScan_(float min, float max)
	{
		DSpectrum<1> spectrum;
		DPeakArray<1>& array = spectrum.getContainer();
		
		for (int x = 0; x != viewport()->width(); x++)
		{
			// get chart coordinates for this pixel column
			PointType px1 = widgetToChart_(QPoint(x, 0));
			PointType px2 = widgetToChart_(QPoint(x + 1, 0));
			
			AreaType area(px1.X(), min, px2.X(), max);
			
			float sum = 0.0f;
			for (QuadTreeType_::Iterator i = trees_[current_data_]->begin(area); i != trees_[current_data_]->end(); ++i)
			{
				sum += i->second->getIntensity();
			}
			
			if (sum > 0.0f)
	
			{
				DPeak<1> peak;
				peak.getPosition()[0] = px1.X();
				peak.setIntensity(sum);
				
				array.push_back(peak);
			}
		}
		
		emit selectedHorz(spectrum);
	}
	
	void Spectrum2DCanvas::createVertScan_(float min, float max)
	{
		DSpectrum<1> spectrum;
		DPeakArray<1>& array = spectrum.getContainer();
		
		for (int y = 0; y != viewport()->height(); y++)
		{
			// get chart coordinates for this pixel column
			PointType py1 = widgetToChart_(QPoint(0, y));
			PointType py2 = widgetToChart_(QPoint(0, y + 1));
			
			AreaType area(min, py1.Y(), max, py2.Y());
			
			float sum = 0.0f;
			for (QuadTreeType_::Iterator i = trees_[current_data_]->begin(area); i != trees_[current_data_]->end(); ++i)
			{
				sum += i->second->getIntensity();
			}
			
			if (sum > 0.0f)
			{
				DPeak<1> peak;
				peak.getPosition()[0] = py1.Y();
				peak.setIntensity(sum);
				
				array.push_back(peak);
			}
		}
		
		emit selectedVert(spectrum);
	}
	
	PreferencesDialogPage* Spectrum2DCanvas::createPreferences(QWidget* parent)
	{
		return new Spectrum2DCanvasPDP(this, parent);
	}
	
	void Spectrum2DCanvas::setDotMode(SignedInt mode)
	{
		prefs_.setValue("Preferences:Dot:Mode", mode);
	}
	
	SignedInt Spectrum2DCanvas::getDotMode()
	{
		if (prefs_.getValue("Preferences:Dot:Mode").isEmpty())
		{
			return 0;
		}
		
		return SignedInt(prefs_.getValue("Preferences:Dot:Mode"));
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
		mapping_info_.setParam(prefs.copy("Preferences:2D:Mapping:",true));
	}
	
	SignedInt Spectrum2DCanvas::finishAdding()
	{
		//cout << "2DCanvas::finishAdding()" << endl;

		current_data_ = getDataSetCount()-1;
		currentDataSet().updateRanges(1);

		//set visibility to true
		layer_visible_.push_back(true);
		show_contours_.push_back(false);
		show_colors_.push_back(false);
		show_points_.push_back(true);
		
		trees_.push_back(0);
		
		
		//if there are spectra with MS-level 1
		if (std::find(currentDataSet().getMSLevels().begin(), currentDataSet().getMSLevels().end(),UnsignedInt(1))!=currentDataSet().getMSLevels().end())
		{
			recalculate_ = true;
			emit sendStatusMessage("constructing quad tree",0);
			
			// find lower left and upper right bound (position and intensity)		
			//values for the current dataset
			disp_ints_.push_back(pair<float,float>(currentDataSet().getMinInt(),currentDataSet().getMaxInt()));
			
			//overall values
			updateRanges_(current_data_,0,1,2);
			
			//cout<<"dataset boudaries MZ: "<< currentDataSet().getMinMZ() << " " << currentDataSet().getMaxMZ() << " RT: " << currentDataSet().getMinRT() << " " << currentDataSet().getMaxRT() << endl;
			//cout<<"Overall new boudaries MZ: "<< overall_data_range_ << endl;
			
			AreaType tmp;
			tmp.assign(overall_data_range_);
			
			if (tmp != visible_area_)
			{ 
				visible_area_.assign(overall_data_range_);
				
				bool insertion_error = false;
				
				for (UnsignedInt data_set=0; data_set<getDataSetCount(); data_set++)
				{
					QuadTreeType_* new_tree = new QuadTreeType_(visible_area_);
	
					for (ExperimentType::Iterator exp_it = getDataSet(data_set).begin(); exp_it != getDataSet(data_set).end(); ++exp_it)
					{
						if (exp_it->getMSLevel()!=1)
						{
							continue;
						}
						
						//cout << endl << endl << "new Spectrum Peaks: " << exp_it->size() << endl;
						for (SpectrumIteratorType i = exp_it->begin(); i != exp_it->end(); ++i)
						{
							//cout << endl << "Peak (RT: " << exp_it->getRetentionTime() << " MZ: " << i->getPosition()[0] << ")" << endl;
							if (i->getIntensity() >= disp_ints_[data_set].first && i->getIntensity() <= disp_ints_[data_set].second)
							{
								try
								{
									new_tree->insert(PointType(i->getPosition()[0],exp_it->getRetentionTime()), &(*i));      
								}
								catch (Exception::IllegalTreeOperation& e)
								{
									// removes problems with identical peak positions
									insertion_error = true;
								}
							}
						}
					}
					
					delete trees_[data_set];
					trees_[data_set] = new_tree;
				}
				if (insertion_error)
				{
					cout << "Warning: Multiple similar peak positions in one data set!" << endl;
				}
				//cout << "Peaks added!" << endl;
			}
			else
			{
				intensityDistributionChange_();
			}
			
			intensityModificationChange_();
			emit sendStatusMessage("",0);
		}
		else
		{
			// create empty tree for empty data sets
			trees_[trees_.size()-1] = new QuadTreeType_(AreaType(0, 0, 0, 0));
		}
	  viewport()->setCursor(Qt::ArrowCursor);
		
		emit layerActivated(this);

		if (parentWidget()!=0)
		{
			dynamic_cast<Spectrum2DWidget*>(parentWidget())->recalculateAxes();
		}
		invalidate_();
		
		return current_data_;
	}
	
	void Spectrum2DCanvas::removeDataSet(int data_set )
	{
		if (data_set >= int(getDataSetCount()))
		{
			return;
		}
	
		//remove the data
		datasets_.erase(datasets_.begin()+data_set);
		delete trees_[data_set];
		trees_.erase(trees_.begin()+data_set);
	
		//remove settings
		layer_visible_.erase(layer_visible_.begin()+data_set);
		show_contours_.erase(show_contours_.begin()+data_set);
		show_colors_.erase(show_colors_.begin()+data_set);
		show_points_.erase(show_points_.begin()+data_set);
		disp_ints_.erase(disp_ints_.begin()+data_set);
		
		//update visible area and boundaries
		recalculateRanges_(0,1,2);

		//cout<<"Overall new boudaries: "<< overall_data_range_<< endl;
		
		AreaType tmp;
		tmp.assign(overall_data_range_);
		if (tmp != visible_area_)
		{ 
			visible_area_.assign(overall_data_range_);
			
			for (UnsignedInt data_set=0; data_set<getDataSetCount(); data_set++)
			{
				QuadTreeType_* new_tree = new QuadTreeType_(visible_area_);		
				//cout << endl << endl << "new Experiment (Spectra: " << getDataSet(data_set).size() << ")" << endl;
				for (ExperimentType::Iterator exp_it = getDataSet(data_set).begin(); exp_it != getDataSet(data_set).end(); ++exp_it)
				{
					if (exp_it->getMSLevel()!=1)
					{
						continue;
					}
					//cout << endl << "new Spectrum (Peaks: " << exp_it->size() << ")" << endl;
					for (SpectrumIteratorType i = exp_it->begin(); i != exp_it->end(); ++i)
					{
						//cout << "Peak (RT: " << exp_it->getRetentionTime() << " MZ: " << i->getPosition()[0] << ")" << endl;
						if (i->getIntensity() >= disp_ints_[data_set].first && i->getIntensity() <= disp_ints_[data_set].second)
						{
							try
							{
								new_tree->insert(PointType(i->getPosition()[0],exp_it->getRetentionTime()), &(*i));        
							}
							catch (Exception::IllegalTreeOperation& e)
							{
								// removes problems with identical peak positions
							}   
						}
					}
				}
				delete trees_[data_set];
				trees_[data_set] = new_tree;
			}
		}
		intensityModificationChange_();
		emit sendStatusMessage("",0);
	
		//update current data set
		if (current_data_ >= getDataSetCount())
		{
			current_data_ = getDataSetCount()-1;
		}
	
		if (datasets_.empty())
		{
			return;
		}
	
		emit layerActivated(this);
		invalidate_();
	}
	
	//change the current data set
	void Spectrum2DCanvas::activateDataSet(int data_set)
	{
		if ((data_set >= int(getDataSetCount())) || data_set==int(current_data_))
		{
			return ;
		}
		current_data_ = data_set;
		emit layerActivated(this);
		
		// no peak is selected
		nearest_peak_ = 0;
	
		invalidate_();
	}

} //namespace

