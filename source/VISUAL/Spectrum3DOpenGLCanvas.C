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
// $Id: Spectrum3DOpenGLCanvas.C,v 1.29 2006/06/09 11:47:50 cfriedle Exp $
// $Author: cfriedle $
// $Maintainer: Cornelia Friedle $
// --------------------------------------------------------------------------

//Qt
#include<qcolor.h>

//Open_GL
#include<qgl.h>
#include<GL/glu.h>

//STL
#include<iostream.h>
//#include <stdio.h> // Header File For Standard Input/Output
//#include <stdarg.h> // Header File For Variable Argument Routines

#include<vector>
#include<iterator>
//OpenMS
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/VISUAL/Spectrum3DOpenGLCanvas.h>
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/VISUAL/AxisTickCalculator.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

namespace OpenMS
{
	
Spectrum3DOpenGLCanvas::Spectrum3DOpenGLCanvas(QWidget *parent, const char* name,Spectrum3DCanvas &canvas_3d)
  : QGLWidget(parent,name),
		canvas_3d_(canvas_3d)
{
	//	setFocusPolicy(QWidget::TabFocus);
	//new and old values for the angels
	xRot_ = 0;
	yRot_ = 0;
	zRot_ = 0;
	xRot_old_ = 0;
	yRot_old_ = 0;
	zRot_old_ = 0;

	//Boundingbox-values
	corner_ = 100.0;
	near_=0.0;
	far_ = 6*corner_;
	zoom_ = 1.25;
	x_1_ = 0.0;
	y_1_ = 0.0;
	x_2_ = 0.0;
	y_2_ = 0.0;	
	//prefernces-values
	showbackview_=false;
	view_mode_ = VIEW_SELECT;
	grid_exists_=false;
	intensity_scale_=false;
}
	
Spectrum3DOpenGLCanvas::  ~Spectrum3DOpenGLCanvas()
{
 	grid_rt_.erase(	grid_rt_.begin(),	grid_rt_.end());
	grid_mz_.erase(grid_mz_.begin(),grid_mz_.end());
	grid_intensity_.erase(grid_intensity_.begin(),grid_intensity_.end());
}

void Spectrum3DOpenGLCanvas::calculateGridLines_()
{
	if(intensity_scale_)
	{
		grid_intensity_log_=  AxisTickCalculator::calcLogGridLines_(log10(int_scale_.min_[0]),log10(int_scale_.max_[0])); 
		grid_intensity_=  AxisTickCalculator::calcGridLines_(int_scale_.min_[0],int_scale_.max_[0],3); 

	}
	else
	{
		grid_intensity_log_=  AxisTickCalculator::calcLogGridLines_(log10(canvas_3d_.overall_data_range_.min_[2]),log10(canvas_3d_.overall_data_range_.max_[2])); 
		grid_intensity_=  AxisTickCalculator::calcGridLines_(canvas_3d_.overall_data_range_.min_[2],canvas_3d_.overall_data_range_.max_[2],3); 
	}
	grid_rt_=  AxisTickCalculator::calcGridLines_(canvas_3d_.overall_data_range_.min_[0],canvas_3d_.overall_data_range_.max_[0],3);
	grid_mz_=  AxisTickCalculator::calcGridLines_(canvas_3d_.overall_data_range_.min_[1],canvas_3d_.overall_data_range_.max_[1],3);
	grid_exists_ = true;
}

void Spectrum3DOpenGLCanvas::resizeGL(int w,int h)
{
	width_ = (float)w;
	heigth_ = (float)h;
	glViewport(0,0,(GLsizei) w, (GLsizei) h); 
}

void Spectrum3DOpenGLCanvas::initializeGL()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-corner_*zoom_,
					corner_*zoom_, 
					-corner_*zoom_,
					corner_*zoom_ ,
					near_,
					far_);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0.0, 0.0, -20.0, 
						0.0, 0.0, 0.0, 
						0.0, 1.0, 0.0);

	if(canvas_3d_.getDataSetCount()!=0)
	{
		
		calculateGridLines_();
		coord_ = makeCoordinates();
		ground_ = makeGround();
		switch (view_mode_)
		{
		case VIEW_SELECT:
			calculateGridLines_();
			if(canvas_3d_.getPrefAsInt("Preferences:3D:IntScale:Mode"))
			{
				stickdata_ = makeDataAsStickLog();
			}
 			else
 			{ 
					stickdata_ =  makeDataAsStick();
			}
			axeslabel_ = makeAxesLabel();	
	 		break;
		case VIEW_TOP:
			stickdata_ = makeDataAsTopView();
			break;
		case VIEW_ZOOM:
			stickdata_ = makeDataAsTopView();
			if(show_zoom_selection_)
			{
				zoomselection_ = makeZoomSelection();
			}
			break;
		}
	}
}

void Spectrum3DOpenGLCanvas::paintGL()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	glLoadIdentity();
	
	glTranslated(0.0, 0.0,-3.0*corner_);
	glRotated(xRot_ / 16.0, 1.0, 0.0, 0.0);
	glRotated(yRot_ / 16.0, 0.0, 1.0, 0.0);
	glRotated(zRot_/16.0, 0.0, 0.0, 1.0);
	glTranslated(0.0, 0.0,3.0*corner_);
	
	QColor color(canvas_3d_.getPrefAsString("Preferences:3D:BackgroundColor").c_str());
 	qglClearColor(color);
	glEnable(GL_DEPTH_TEST);
	
	if(canvas_3d_.getDataSetCount()!=0)
	{
		glCallList(ground_);
		glCallList(stickdata_);	
		
		if( view_mode_ !=VIEW_ZOOM && view_mode_ != VIEW_TOP)
			{
				glCallList(axeslabel_);
				if(grid_exists_)
					{
						makeFont();
						paintAxesScale();
					}
			}
		if(show_zoom_selection_)
			
			{
				glCallList(zoomselection_);
			}
		glCallList(coord_);
		
	}
	
}

void Spectrum3DOpenGLCanvas::paintAxesScale()
{
	GLfloat black[3] = { 0.0, 0.0, 0.0 };
	glColor3fv(black);	//rt-axe
	if(yRot_> 280*16 ^ yRot_<80*16)
		{
			if(zoom_<3)
				{
					for(UnsignedInt i = 0;i<grid_rt_[0].size();i++)
						{
							String result(grid_rt_[0][i]);
							glPushAttrib (GL_LIST_BIT);
							glListBase(fontOffset);
							
							glRasterPos3d(-corner_-result.length()+scaledRT(grid_rt_[0][i]), 
														-corner_-9.0,
										-near_-2*corner_+5.0);
							glCallLists(result.length(), GL_UNSIGNED_BYTE,(GLubyte*)result.c_str());
							
							glPopAttrib ();
						}
				}
			if(zoom_<2.0)
				{
					for(UnsignedInt i = 0;i<grid_rt_[1].size();i++)
						{
							String result(grid_rt_[1][i]);
							glPushAttrib (GL_LIST_BIT);
							glListBase(fontOffset);
							glRasterPos3d(-corner_-result.length()+scaledRT(grid_rt_[1][i]), 
														-corner_-9.0,
														-near_-2*corner_+5.0);
							glCallLists(result.length(), GL_UNSIGNED_BYTE,(GLubyte*)result.c_str());
							glPopAttrib ();
							
						}
				}
			if(width_>1000 && heigth_>700&& zoom_<1.5)
				{
					for(UnsignedInt i = 0;i<grid_rt_[2].size();i++)
						{
							String result(grid_rt_[2][i]);
							glPushAttrib (GL_LIST_BIT);
							glListBase(fontOffset);
							glRasterPos3d(-corner_-result.length()+scaledRT(grid_rt_[2][i]), 
														-corner_-9.0,
														-near_-2*corner_+5.0);
							glCallLists(result.length(), GL_UNSIGNED_BYTE,(GLubyte*)result.c_str());
							glPopAttrib ();
							
						}
				}
		}

	if(yRot_>10*16 && yRot_<190*16 )
		{
			if(zoom_<3)
			{
				for(UnsignedInt i = 0;i<grid_mz_[0].size();i++)
				{
					String result(grid_mz_[0][i]);
					glPushAttrib (GL_LIST_BIT);
					glListBase(fontOffset);
						glRasterPos3d(-corner_-result.length()-10.0, 
													-corner_-9.0,
													-near_-2*corner_-scaledMZ(grid_mz_[0][i]));
						glCallLists(result.length(), GL_UNSIGNED_BYTE,(GLubyte*)result.c_str());
						glPopAttrib ();
				}
			}
			if(zoom_<2.0)
				{
				for(UnsignedInt i = 0;i<grid_mz_[1].size();i++)
					{
						String result(grid_mz_[1][i]);
						glPushAttrib (GL_LIST_BIT);
						glListBase(fontOffset);
						glRasterPos3d(-corner_-result.length()-10.0, 
													-corner_-9.0,
													-near_-2*corner_-scaledMZ(grid_mz_[1][i]));
						glCallLists(result.length(), GL_UNSIGNED_BYTE,(GLubyte*)result.c_str());
						glPopAttrib ();
						
					}
			}
			if(width_>1000 && heigth_>700&& zoom_<1.5)
				{
					for(UnsignedInt i = 0;i<grid_mz_[2].size();i++)
						{
							String result(grid_mz_[2][i]);
							glPushAttrib (GL_LIST_BIT);
							glListBase(fontOffset);
							glRasterPos3d(-corner_-result.length()-10.0, 
														-corner_-9.0,
														-near_-2*corner_-scaledMZ(grid_mz_[2][i]));
							glCallLists(result.length(), GL_UNSIGNED_BYTE,(GLubyte*)result.c_str());
							glPopAttrib ();
							
						}
				}
	}
	if(canvas_3d_.getPrefAsInt("Preferences:3D:IntScale:Mode"))
	{
			if(zoom_<3)
			{
				for(UnsignedInt i = 0;i<grid_intensity_log_[0].size();i++)
				{
					String result(grid_intensity_log_[0][i]);
					glPushAttrib (GL_LIST_BIT);
					glListBase(fontOffset);
					glRasterPos3d(-corner_-result.length()-3.0, 
												-corner_+scaledLogIntensity(grid_intensity_log_[0][i]),
												-near_-2*corner_);
					glCallLists(result.length(), GL_UNSIGNED_BYTE,(GLubyte*)result.c_str());
					glPopAttrib ();
				}
			}
		}
	else
	{
		if(zoom_<3)
			{
			for(UnsignedInt i = 0;i<grid_intensity_[0].size();i++)
				{
					String result(grid_intensity_[0][i]);
					glPushAttrib (GL_LIST_BIT);
					glListBase(fontOffset);
					glRasterPos3d(-corner_-result.length()-8.0, 
												-corner_+scaledIntensity(grid_intensity_[0][i]),
												-near_-2*corner_);
					glCallLists(result.length(), GL_UNSIGNED_BYTE,(GLubyte*)result.c_str());
					glPopAttrib ();
				}
			}
		if(zoom_<2.0)
			{
				for(UnsignedInt i = 0;i<grid_intensity_[1].size();i++)
					{
						String result(grid_intensity_[1][i]);
						glPushAttrib (GL_LIST_BIT);
						glListBase(fontOffset);
						glRasterPos3d(-corner_-result.length()-8.0, 
													-corner_+scaledIntensity(grid_intensity_[1][i]),
													-near_-2*corner_);
						glCallLists(result.length(), GL_UNSIGNED_BYTE,(GLubyte*)result.c_str());
						glPopAttrib ();
						
					}
			}
		if(width_>1000 && heigth_>700&& zoom_<1.5)
			{
				for(UnsignedInt i = 0;i<grid_intensity_[2].size();i++)
					{
						String result(grid_intensity_[2][i]);
						glPushAttrib (GL_LIST_BIT);
						glListBase(fontOffset);
						glRasterPos3d(-corner_-result.length()-8.0, 
												-corner_+scaledIntensity(grid_intensity_[2][i]),
													-near_-2*corner_);
						glCallLists(result.length(), GL_UNSIGNED_BYTE,(GLubyte*)result.c_str());
						glPopAttrib ();
				
					}
			}
	}
}

void Spectrum3DOpenGLCanvas::makeFont()
{	
	GLubyte rasters[][13] = {
		{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}, 
		{0x00, 0x00, 0x18, 0x18, 0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18}, 
		{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x36, 0x36, 0x36, 0x36}, 
		{0x00, 0x00, 0x00, 0x66, 0x66, 0xff, 0x66, 0x66, 0xff, 0x66, 0x66, 0x00, 0x00}, 
		{0x00, 0x00, 0x18, 0x7e, 0xff, 0x1b, 0x1f, 0x7e, 0xf8, 0xd8, 0xff, 0x7e, 0x18},		
		{0x00, 0x00, 0x0e, 0x1b, 0xdb, 0x6e, 0x30, 0x18, 0x0c, 0x76, 0xdb, 0xd8, 0x70}, 
		{0x00, 0x00, 0x7f, 0xc6, 0xcf, 0xd8, 0x70, 0x70, 0xd8, 0xcc, 0xcc, 0x6c, 0x38}, 
		{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x18, 0x1c, 0x0c, 0x0e}, 
		{0x00, 0x00, 0x0c, 0x18, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x18, 0x0c}, 
		{0x00, 0x00, 0x30, 0x18, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x18, 0x30},
		{0x00, 0x00, 0x00, 0x00, 0x99, 0x5a, 0x3c, 0xff, 0x3c, 0x5a, 0x99, 0x00, 0x00}, 
		{0x00, 0x00, 0x00, 0x18, 0x18, 0x18, 0xff, 0xff, 0x18, 0x18, 0x18, 0x00, 0x00}, 
		{0x00, 0x00, 0x30, 0x18, 0x1c, 0x1c, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}, 
		{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00, 0x00}, 
		{0x00, 0x00, 0x00, 0x38, 0x38, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},		
		{0x00, 0x60, 0x60, 0x30, 0x30, 0x18, 0x18, 0x0c, 0x0c, 0x06, 0x06, 0x03, 0x03}, 
		{0x00, 0x00, 0x3c, 0x66, 0xc3, 0xe3, 0xf3, 0xdb, 0xcf, 0xc7, 0xc3, 0x66, 0x3c}, 
		{0x00, 0x00, 0x7e, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x78, 0x38, 0x18}, 
		{0x00, 0x00, 0xff, 0xc0, 0xc0, 0x60, 0x30, 0x18, 0x0c, 0x06, 0x03, 0xe7, 0x7e}, 
		{0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x07, 0x7e, 0x07, 0x03, 0x03, 0xe7, 0x7e},
		{0x00, 0x00, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0xff, 0xcc, 0x6c, 0x3c, 0x1c, 0x0c}, 
		{0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x07, 0xfe, 0xc0, 0xc0, 0xc0, 0xc0, 0xff}, 
		{0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc7, 0xfe, 0xc0, 0xc0, 0xc0, 0xe7, 0x7e}, 
		{0x00, 0x00, 0x30, 0x30, 0x30, 0x30, 0x18, 0x0c, 0x06, 0x03, 0x03, 0x03, 0xff}, 
		{0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xe7, 0x7e, 0xe7, 0xc3, 0xc3, 0xe7, 0x7e},
		{0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x03, 0x7f, 0xe7, 0xc3, 0xc3, 0xe7, 0x7e}, 
		{0x00, 0x00, 0x00, 0x38, 0x38, 0x00, 0x00, 0x38, 0x38, 0x00, 0x00, 0x00, 0x00}, 
		{0x00, 0x00, 0x30, 0x18, 0x1c, 0x1c, 0x00, 0x00, 0x1c, 0x1c, 0x00, 0x00, 0x00}, 
		{0x00, 0x00, 0x06, 0x0c, 0x18, 0x30, 0x60, 0xc0, 0x60, 0x30, 0x18, 0x0c, 0x06}, 
		{0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0x00, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00}, 
		{0x00, 0x00, 0x60, 0x30, 0x18, 0x0c, 0x06, 0x03, 0x06, 0x0c, 0x18, 0x30, 0x60}, 
		{0x00, 0x00, 0x18, 0x00, 0x00, 0x18, 0x18, 0x0c, 0x06, 0x03, 0xc3, 0xc3, 0x7e}, 
		{0x00, 0x00, 0x3f, 0x60, 0xcf, 0xdb, 0xd3, 0xdd, 0xc3, 0x7e, 0x00, 0x00, 0x00}, 
		{0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xff, 0xc3, 0xc3, 0xc3, 0x66, 0x3c, 0x18}, 
		{0x00, 0x00, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe}, 
		{0x00, 0x00, 0x7e, 0xe7, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xe7, 0x7e}, 
		{0x00, 0x00, 0xfc, 0xce, 0xc7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc7, 0xce, 0xfc}, 
		{0x00, 0x00, 0xff, 0xc0, 0xc0, 0xc0, 0xc0, 0xfc, 0xc0, 0xc0, 0xc0, 0xc0, 0xff}, 
		{0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xfc, 0xc0, 0xc0, 0xc0, 0xff}, 
		{0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xcf, 0xc0, 0xc0, 0xc0, 0xc0, 0xe7, 0x7e}, 
		{0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xff, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3}, 
		{0x00, 0x00, 0x7e, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x7e}, 
		{0x00, 0x00, 0x7c, 0xee, 0xc6, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06}, 
		{0x00, 0x00, 0xc3, 0xc6, 0xcc, 0xd8, 0xf0, 0xe0, 0xf0, 0xd8, 0xcc, 0xc6, 0xc3}, 
		{0x00, 0x00, 0xff, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0}, 
		{0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xdb, 0xff, 0xff, 0xe7, 0xc3}, 
		{0x00, 0x00, 0xc7, 0xc7, 0xcf, 0xcf, 0xdf, 0xdb, 0xfb, 0xf3, 0xf3, 0xe3, 0xe3}, 
		{0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xe7, 0x7e}, 
		{0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe}, 
		{0x00, 0x00, 0x3f, 0x6e, 0xdf, 0xdb, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0x66, 0x3c}, 
		{0x00, 0x00, 0xc3, 0xc6, 0xcc, 0xd8, 0xf0, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe}, 
		{0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x07, 0x7e, 0xe0, 0xc0, 0xc0, 0xe7, 0x7e}, 
		{0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0xff}, 
		{0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3}, 
		{0x00, 0x00, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3}, 
		{0x00, 0x00, 0xc3, 0xe7, 0xff, 0xff, 0xdb, 0xdb, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3}, 
		{0x00, 0x00, 0xc3, 0x66, 0x66, 0x3c, 0x3c, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3}, 
		{0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3}, 
		{0x00, 0x00, 0xff, 0xc0, 0xc0, 0x60, 0x30, 0x7e, 0x0c, 0x06, 0x03, 0x03, 0xff}, 
		{0x00, 0x00, 0x3c, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x3c}, 
		{0x00, 0x03, 0x03, 0x06, 0x06, 0x0c, 0x0c, 0x18, 0x18, 0x30, 0x30, 0x60, 0x60}, 
		{0x00, 0x00, 0x3c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x3c}, 
		{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xc3, 0x66, 0x3c, 0x18}, 
		{0xff, 0xff, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}, 
		{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x18, 0x38, 0x30, 0x70}, 
		{0x00, 0x00, 0x7f, 0xc3, 0xc3, 0x7f, 0x03, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00}, 
		{0x00, 0x00, 0xfe, 0xc3, 0xc3, 0xc3, 0xc3, 0xfe, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0}, 
		{0x00, 0x00, 0x7e, 0xc3, 0xc0, 0xc0, 0xc0, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00}, 
		{0x00, 0x00, 0x7f, 0xc3, 0xc3, 0xc3, 0xc3, 0x7f, 0x03, 0x03, 0x03, 0x03, 0x03}, 
		{0x00, 0x00, 0x7f, 0xc0, 0xc0, 0xfe, 0xc3, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00}, 
		{0x00, 0x00, 0x30, 0x30, 0x30, 0x30, 0x30, 0xfc, 0x30, 0x30, 0x30, 0x33, 0x1e}, 
		{0x7e, 0xc3, 0x03, 0x03, 0x7f, 0xc3, 0xc3, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00}, 
		{0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xfe, 0xc0, 0xc0, 0xc0, 0xc0}, 
		{0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x00, 0x00, 0x18, 0x00}, 
		{0x38, 0x6c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x00, 0x00, 0x0c, 0x00}, 
		{0x00, 0x00, 0xc6, 0xcc, 0xf8, 0xf0, 0xd8, 0xcc, 0xc6, 0xc0, 0xc0, 0xc0, 0xc0}, 
		{0x00, 0x00, 0x7e, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x78}, 
		{0x00, 0x00, 0xdb, 0xdb, 0xdb, 0xdb, 0xdb, 0xdb, 0xfe, 0x00, 0x00, 0x00, 0x00}, 
		{0x00, 0x00, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0xfc, 0x00, 0x00, 0x00, 0x00}, 
		{0x00, 0x00, 0x7c, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0x7c, 0x00, 0x00, 0x00, 0x00}, 
		{0xc0, 0xc0, 0xc0, 0xfe, 0xc3, 0xc3, 0xc3, 0xc3, 0xfe, 0x00, 0x00, 0x00, 0x00}, 
		{0x03, 0x03, 0x03, 0x7f, 0xc3, 0xc3, 0xc3, 0xc3, 0x7f, 0x00, 0x00, 0x00, 0x00}, 
		{0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xe0, 0xfe, 0x00, 0x00, 0x00, 0x00}, 
		{0x00, 0x00, 0xfe, 0x03, 0x03, 0x7e, 0xc0, 0xc0, 0x7f, 0x00, 0x00, 0x00, 0x00}, 
		{0x00, 0x00, 0x1c, 0x36, 0x30, 0x30, 0x30, 0x30, 0xfc, 0x30, 0x30, 0x30, 0x00}, 
		{0x00, 0x00, 0x7e, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0x00, 0x00, 0x00, 0x00}, 
		{0x00, 0x00, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3, 0xc3, 0x00, 0x00, 0x00, 0x00}, 
		{0x00, 0x00, 0xc3, 0xe7, 0xff, 0xdb, 0xc3, 0xc3, 0xc3, 0x00, 0x00, 0x00, 0x00}, 
		{0x00, 0x00, 0xc3, 0x66, 0x3c, 0x18, 0x3c, 0x66, 0xc3, 0x00, 0x00, 0x00, 0x00}, 
		{0xc0, 0x60, 0x60, 0x30, 0x18, 0x3c, 0x66, 0x66, 0xc3, 0x00, 0x00, 0x00, 0x00}, 
		{0x00, 0x00, 0xff, 0x60, 0x30, 0x18, 0x0c, 0x06, 0xff, 0x00, 0x00, 0x00, 0x00}, 
		{0x00, 0x00, 0x0f, 0x18, 0x18, 0x18, 0x38, 0xf0, 0x38, 0x18, 0x18, 0x18, 0x0f}, 
		{0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18}, 
		{0x00, 0x00, 0xf0, 0x18, 0x18, 0x18, 0x1c, 0x0f, 0x1c, 0x18, 0x18, 0x18, 0xf0}, 
		{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x06, 0x8f, 0xf1, 0x60, 0x00, 0x00, 0x00} 
	};
	GLuint i;
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	
	fontOffset = glGenLists (128);
	
	for ( i = 32; i < 127; i++) {
		glNewList(i+fontOffset, GL_COMPILE);
		glBitmap(8, 13, 0.0, 0.0, 8.0, 0.0, rasters[i-32]);
		glEndList();
	}
}

GLuint Spectrum3DOpenGLCanvas::makeGround()
{
	GLuint list = glGenLists(1);
	glNewList(list,GL_COMPILE);
	QColor color(canvas_3d_.getPrefAsString("Preferences:3D:BackgroundColor").c_str());
 	qglColor(color);
	glBegin(GL_QUADS);
	glVertex3d(-corner_-1.0, 
						 -corner_-1.0,
						 -near_-2*corner_+1.0);
	glVertex3d( -corner_-1.0, 
							-corner_-1.0,
							-far_+2*corner_);
	glVertex3d( corner_-1.0, 
							-corner_-1.0,
							-far_+2*corner_);
	glVertex3d(corner_-1.0, 
						 -corner_-1.0,
						 -near_-2*corner_+1.0);
	glEnd();
	glEndList();
	return list;
}


GLuint Spectrum3DOpenGLCanvas::makeZoomSelection()
{
	GLuint list = glGenLists(1);
	glNewList(list,GL_COMPILE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glColor4ub(100, 0, 0, 80);	
	glBegin(GL_QUADS);
	glVertex3d((GLfloat)x_1_,
 						 corner_,
 						 (GLfloat)y_1_);
	glVertex3d((GLfloat)x_1_,
						 corner_,
						 (GLfloat)y_2_);
	glVertex3d((GLfloat)x_2_,
						 corner_,
						 (GLfloat)y_2_);
	glVertex3d((GLfloat)x_2_,
						 corner_,
						 (GLfloat)y_1_);
	glEnd();
	glEndList();

	return list;
}

GLuint Spectrum3DOpenGLCanvas::makeCoordinates()
{
	GLuint list = glGenLists(1);
	glNewList(list,GL_COMPILE);
	glLineWidth(3.0);
	glShadeModel(GL_FLAT);
	glBegin(GL_LINES);
	QColor color(canvas_3d_.getPrefAsString("Preferences:3D:AxesColor").c_str());
	qglColor(color);
	//x_achse
	glVertex3d(-corner_-1.0, 
						 -corner_-1.0,
						 -near_-2*corner_+1.0);
	glVertex3d( corner_+1.0, 
							-corner_-1.0,
							-near_-2*corner_+1.0);
	//z-achse
	glVertex3d(-corner_-1.0, 
						 -corner_-1.0,
						 -near_-2*corner_+1.0);
	glVertex3d( -corner_-1.0, 
							-corner_-1.0,
							-far_+2*corner_);
	//y-achse
	glVertex3d(-corner_-1.0, 
						 -corner_-1.0,
						 -near_-2*corner_+1.0);
	glVertex3d( -corner_-1.0, 
							corner_+1.0,
							-near_-2*corner_+1.0);
	glEnd();
	glEndList();
	return list;
}

GLuint Spectrum3DOpenGLCanvas::makeDataAsTopView()
{
	GLuint list = glGenLists(1);
	glNewList(list,GL_COMPILE);
	glPointSize(3.0);
	if(canvas_3d_.getPrefAsInt("Preferences:3D:Shade:Mode"))
		{
			glShadeModel(GL_SMOOTH); 
		}
	else
		{
			glShadeModel(GL_FLAT); 
		}
		for(UnsignedInt i =0;i<canvas_3d_.getDataSetCount();i++)
			{
			recalculateDotGradient_(i);		 
			if(canvas_3d_.layer_visible_[i]==true)
				{	
					canvas_3d_.getDataSet(i).sortSpectra(false);
					
					for (Spectrum3DCanvas::ExperimentType::Iterator spec_it = canvas_3d_.getDataSet(i).RTBegin(canvas_3d_.overall_data_range_.min_[0]); 
						 spec_it != canvas_3d_.getDataSet(i).RTEnd(canvas_3d_.overall_data_range_.max_[0]); 
						 ++spec_it)
					{
						//	canvas_3d_.getDataSet(i).sortSpectra(true);
					
						for (BaseSpectrum::Iterator it = spec_it->MZBegin(canvas_3d_.overall_data_range_.min_[1]); it!=spec_it->MZEnd(canvas_3d_.overall_data_range_.max_[1]); ++it)
						{	
						
							if(it->getIntensity()>= canvas_3d_.disp_ints_[i].first && it->getIntensity()<= canvas_3d_.disp_ints_[i].second)
							{
								
								if(int(canvas_3d_.getPref("Preferences:3D:Dot:Mode")))
								{
									qglColor(QColor( gradient_.precalculatedColorAt(it->getIntensity())));
									glBegin(GL_POINTS);
									qglColor(QColor( gradient_.precalculatedColorAt(it->getIntensity())));
									glVertex3d(-corner_+(GLfloat)scaledRT(spec_it->getRetentionTime()), 
														 -corner_,
														 -near_-2*corner_-(GLfloat)scaledMZ(it->getPosition()[0]));
									glEnd();					
								}
								else
								{
									glBegin(GL_POINTS);
									qglColor(black);			
									glVertex3d(-corner_+(GLfloat)scaledRT(spec_it->getRetentionTime()), 
														 -corner_,
														 -near_-2*corner_-(GLfloat)scaledMZ(it->getPosition()[0]));
									glEnd();					
								}
							}
						}
					}
				}
			}
		glEndList();
		return list; 
}

GLuint Spectrum3DOpenGLCanvas::makeDataAsStick()
{
	GLuint list = glGenLists(1);
	glNewList(list,GL_COMPILE);
	glLineWidth(canvas_3d_.getPref("Preferences:3D:Dot:LineWidth"));
	if(canvas_3d_.getPrefAsInt("Preferences:3D:Shade:Mode"))
		{
			glShadeModel(GL_SMOOTH); 
		}
	else
		{
			glShadeModel(GL_FLAT); 
		}
	
	for(UnsignedInt i =0;i<canvas_3d_.getDataSetCount();i++)
		{
			recalculateDotGradient_(i);		 
			if(canvas_3d_.isDataSetVisible(i))
			{	
				canvas_3d_.getDataSet(i).sortSpectra(false);
		
				for (Spectrum3DCanvas::ExperimentType::Iterator spec_it = canvas_3d_.getDataSet(i).RTBegin(canvas_3d_.overall_data_range_.min_[0]); 
						 spec_it != canvas_3d_.getDataSet(i).RTEnd(canvas_3d_.overall_data_range_.max_[0]); 
						 ++spec_it)
				{
					//	canvas_3d_.getDataSet(i).sortSpectra(true);
		
					for (BaseSpectrum::Iterator it = spec_it->MZBegin(canvas_3d_.overall_data_range_.min_[1]); it!=spec_it->MZEnd(canvas_3d_.overall_data_range_.max_[1]); ++it)
					{			
						if(it->getIntensity()>= canvas_3d_.disp_ints_[i].first && it->getIntensity()<= canvas_3d_.disp_ints_[i].second)
						{
							glBegin(GL_LINES);
							if(int(canvas_3d_.getPref("Preferences:3D:Dot:Mode")))
							{
								qglColor(QColor( gradient_.precalculatedColorAt(canvas_3d_.overall_data_range_.min_[2])));
								glVertex3d(-corner_+(GLfloat)scaledRT(spec_it->getRetentionTime()), 
													 -corner_,
													 -near_-2*corner_-(GLfloat)scaledMZ(it->getPosition()[0]));
								qglColor(QColor( gradient_.precalculatedColorAt(it->getIntensity())));
							
								glVertex3d(-corner_+(GLfloat)scaledRT(spec_it->getRetentionTime()),
													 -corner_+(GLfloat)scaledIntensity(it->getIntensity()),
													 -near_-2*corner_-(GLfloat)scaledMZ(it->getPosition()[0]));
							}
							else
							{
								qglColor(black);			
								glVertex3d(-corner_+(GLfloat)scaledRT(spec_it->getRetentionTime()), 
													 -corner_,
												 -near_-2*corner_-(GLfloat)scaledMZ(it->getPosition()[0]));
								glVertex3d(-corner_+(GLfloat)scaledRT(spec_it->getRetentionTime()),
													 -corner_+(GLfloat)scaledIntensity(it->getIntensity()),
													 -near_-2*corner_-(GLfloat)scaledMZ(it->getPosition()[0]));
							}
							glEnd();
						}
					}
					
				}
			}
			
		}
	glEndList();
	return list; 
}

GLuint Spectrum3DOpenGLCanvas::makeDataAsStickLog()
{	
	GLuint list = glGenLists(1);
	glNewList(list,GL_COMPILE);
	glLineWidth(canvas_3d_.getPref("Preferences:3D:Dot:LineWidth"));
	if(canvas_3d_.getPrefAsInt("Preferences:3D:Shade:Mode"))
	{
		glShadeModel(GL_SMOOTH); 
	}
	else
	{
		glShadeModel(GL_FLAT); 
	}
	for(UnsignedInt i =0;i<canvas_3d_.getDataSetCount();i++)
		{
			canvas_3d_.getDataSet(i).sortSpectra(false);
			if(canvas_3d_.isDataSetVisible(i))
			{		
				recalculateDotGradientLog_(i);	
				
				for (Spectrum3DCanvas::ExperimentType::Iterator spec_it = canvas_3d_.getDataSet(i).RTBegin(canvas_3d_.overall_data_range_.min_[0]); 
						 spec_it != canvas_3d_.getDataSet(i).RTEnd(canvas_3d_.overall_data_range_.max_[0]); 
						 ++spec_it)
				{
					//	canvas_3d_.getDataSet(i).sortSpectra(true);
					for (BaseSpectrum::Iterator it = spec_it->MZBegin(canvas_3d_.overall_data_range_.min_[1]); it!=spec_it->MZEnd(canvas_3d_.overall_data_range_.max_[1]); ++it)
					{			
						if(it->getIntensity()>= canvas_3d_.disp_ints_[i].first && it->getIntensity()<= canvas_3d_.disp_ints_[i].second)
						{
							glBegin(GL_LINES);
							if(int(canvas_3d_.getPref("Preferences:3D:Dot:Mode")))
							{
									qglColor(QColor( gradient_.precalculatedColorAt(log10(canvas_3d_.overall_data_range_.min_[2]))));
								
								glVertex3d(-corner_+(GLfloat)scaledRT(spec_it->getRetentionTime()), 
													 -corner_,
													 -near_-2*corner_-(GLfloat)scaledMZ(it->getPosition()[0]));
								qglColor(QColor( gradient_.precalculatedColorAt(log10(it->getIntensity()))));
								glVertex3d(-corner_+(GLfloat)scaledRT(spec_it->getRetentionTime()),
													 -corner_+(GLfloat)scaledLogIntensity(log10(it->getIntensity())),
													 -near_-2*corner_-(GLfloat)scaledMZ(it->getPosition()[0]));
							}
							else
							{
								qglColor(black);			
															glVertex3d(-corner_+(GLfloat)scaledRT(spec_it->getRetentionTime()), 
																				 -corner_,
																				 -near_-2*corner_-(GLfloat)scaledMZ(it->getPosition()[0]));
															glVertex3d(-corner_+(GLfloat)scaledRT(spec_it->getRetentionTime()),
																				 -corner_+(GLfloat)scaledLogIntensity(log10(it->getIntensity())),
																				 -near_-2*corner_-(GLfloat)scaledMZ(it->getPosition()[0]));
							}
							glEnd();
						}
					}
					
				}
			}
		}
	glEndList();
	return list; 
}

GLuint Spectrum3DOpenGLCanvas::makeAxesLabel()
{
	GLuint list = glGenLists(1);
	glNewList(list,GL_COMPILE);
	glShadeModel(GL_FLAT);
	glLineWidth(2.0);
	glBegin(GL_LINES);
	QColor color(canvas_3d_.getPrefAsString("Preferences:3D:AxesColor").c_str());
	qglColor(color);//x_achse
	//RT
	for(UnsignedInt i = 0;i<grid_rt_[0].size();i++)
	{
		glVertex3d(-corner_+scaledRT(grid_rt_[0][i]), 
							 -corner_-1.0,
							 -near_-2*corner_+1.0);
		glVertex3d( -corner_+scaledRT(grid_rt_[0][i]), 
								-corner_+4.0,
								-near_-2*corner_+1.0);
	}
	for(UnsignedInt i = 0;i<grid_rt_[1].size();i++)
	{
		glVertex3d(-corner_+scaledRT(grid_rt_[1][i]), 
							 -corner_-1.0,
							 -near_-2*corner_+1.0);
		glVertex3d( -corner_+scaledRT(grid_rt_[1][i]), 
								-corner_+3.0,
								-near_-2*corner_+1.0);
	}
	for(UnsignedInt i = 0;i<grid_rt_[2].size();i++)
	{
		glVertex3d(-corner_+scaledRT(grid_rt_[2][i]), 
							 -corner_-1.0,
							 -near_-2*corner_+1.0);
		glVertex3d( -corner_+scaledRT(grid_rt_[2][i]), 
								-corner_+2.0,
								-near_-2*corner_+1.0);
	}
	//MZ
	for(UnsignedInt i = 0;i<grid_mz_[0].size();i++)
	{
		glVertex3d(-corner_-1.0, 
							 -corner_-1.0,
							 -near_-2*corner_-scaledMZ(grid_mz_[0][i]));
		glVertex3d( -corner_-1.0,
								-corner_+4.0,
								-near_-2*corner_-scaledMZ(grid_mz_[0][i]));
	}
	for(UnsignedInt i = 0;i<grid_mz_[1].size();i++)
	{
		glVertex3d(-corner_-1.0, 
							 -corner_-1.0,
							 -near_-2*corner_-scaledMZ(grid_mz_[1][i]));
		glVertex3d( -corner_-1.0,
								-corner_+3.0,
								-near_-2*corner_-scaledMZ(grid_mz_[1][i]));
	}
	for(UnsignedInt i = 0;i<grid_mz_[2].size();i++)
	{
		glVertex3d(-corner_-1.0, 
							 -corner_-1.0,
							 -near_-2*corner_-scaledMZ(grid_mz_[2][i]));
		glVertex3d( -corner_-1.0,
								-corner_+2.0,
								-near_-2*corner_-scaledMZ(grid_mz_[2][i]));
	}
	
	if(canvas_3d_.getPrefAsInt("Preferences:3D:IntScale:Mode"))
		{
			for(UnsignedInt i = 0;i<grid_intensity_log_[0].size();i++)
				{	
					glVertex3d(-corner_-1.0, 
										 -corner_+scaledLogIntensity(grid_intensity_log_[0][i]),
										 -near_-2*corner_+1.0);
					glVertex3d( -corner_+3.0, 
											-corner_+scaledLogIntensity(grid_intensity_log_[0][i]),
											-near_-2*corner_-3.0);
				}
			for(UnsignedInt i = 0;i<grid_intensity_log_[1].size();i++)
				{
					glVertex3d(-corner_-1.0, 
										 -corner_+scaledLogIntensity(grid_intensity_log_[1][i]),
										 -near_-2*corner_+1.0);
					glVertex3d( -corner_+2.0, 
											-corner_+scaledLogIntensity(grid_intensity_log_[1][i]),
											-near_-2*corner_-2.0);
				}
		}
	else
		{
			for(UnsignedInt i = 0;i<grid_intensity_[0].size();i++)
			{	
				glVertex3d(-corner_-1.0, 
									 -corner_+scaledIntensity(grid_intensity_[0][i]),
									 -near_-2*corner_+1.0);
				glVertex3d( -corner_+4.0, 
										-corner_+scaledIntensity(grid_intensity_[0][i]),
										-near_-2*corner_-4.0);
			}
			for(UnsignedInt i = 0;i<grid_intensity_[1].size();i++)
				{
					glVertex3d(-corner_-1.0, 
										 -corner_+scaledIntensity(grid_intensity_[1][i]),
										 -near_-2*corner_+1.0);
					glVertex3d( -corner_+3.0, 
											-corner_+scaledIntensity(grid_intensity_[1][i]),
											-near_-2*corner_-3.0);
			}
			for(UnsignedInt i = 0;i<grid_intensity_[2].size();i++)
				{ 
					glVertex3d(-corner_-1.0, 
										 -corner_+scaledIntensity(grid_intensity_[2][i]),
										 -near_-2*corner_+1.0);
				glVertex3d( -corner_+2.0, 
										-corner_+scaledIntensity(grid_intensity_[2][i]),
										-near_-2*corner_-2.0);
				}
		} 
	
	glEnd();
	glEndList();
	return list; 
}

double Spectrum3DOpenGLCanvas::scaledRT(double rt)
{
	double scaledrt = rt - canvas_3d_.overall_data_range_.min_[0];
	scaledrt = scaledrt * 2.0 * corner_ /(canvas_3d_.overall_data_range_.max_[0]-canvas_3d_.overall_data_range_.min_[0]);
	return scaledrt;
}

double Spectrum3DOpenGLCanvas::scaledInversRT(double rt)
{
	double i_rt =(rt* canvas_3d_.overall_data_range_.max_[0] -canvas_3d_.overall_data_range_.min_[0]*rt);
	i_rt = i_rt/200.0;
	i_rt = i_rt + canvas_3d_.overall_data_range_.min_[0]; 
	return i_rt;
}

double Spectrum3DOpenGLCanvas::scaledMZ(double mz)
{
	double scaledmz = mz - canvas_3d_.overall_data_range_.min_[1];
	scaledmz = scaledmz * 2.0 * corner_/(canvas_3d_.overall_data_range_.max_[1]-canvas_3d_.overall_data_range_.min_[1])/*dis_mz_*/;
	return scaledmz;
}


double Spectrum3DOpenGLCanvas::scaledInversMZ(double mz)
{
	double i_mz =(mz *canvas_3d_.overall_data_range_.max_[1] - mz *canvas_3d_.overall_data_range_.min_[1]);
	i_mz = i_mz/200;
	i_mz = i_mz + canvas_3d_.overall_data_range_.min_[1]; 
	return i_mz;
}

double Spectrum3DOpenGLCanvas::scaledIntensity(double intensity)
{
	double scaledintensity;
	if(intensity_scale_)
	{
		scaledintensity = intensity -int_scale_.min_[0];
		scaledintensity = ( scaledintensity * 2.0 * corner_)/(int_scale_.max_[0]-int_scale_.min_[0]);
	}
	else
	{
		scaledintensity = intensity -canvas_3d_.overall_data_range_.min_[2];
		scaledintensity = ( scaledintensity * 2.0 * corner_)/(canvas_3d_.overall_data_range_.max_[2]-canvas_3d_.overall_data_range_.min_[2]);
	}
	return scaledintensity;
}

double Spectrum3DOpenGLCanvas::scaledLogIntensity(double intensity)
{
	double scaledint = intensity -log10(canvas_3d_.overall_data_range_.min_[2]);	
	scaledint =(scaledint * 2.0 * corner_)/(log10(canvas_3d_.overall_data_range_.max_[2])-log10(canvas_3d_.overall_data_range_.min_[2]));
	if(scaledint<0)
	{
		return 0.0;
	}
	return scaledint ;
}

void Spectrum3DOpenGLCanvas::setRotationX(int angle)
{
	normalizeAngle(&angle);
	if (angle != xRot_) 
	{
		xRot_ = angle;
		updateGL();
	}
}

void Spectrum3DOpenGLCanvas::setRotationY(int angle)
{
	normalizeAngle(&angle);
	if (angle != yRot_) 
	{
		yRot_ = angle;
		updateGL();
	}
}

void Spectrum3DOpenGLCanvas::setRotationZ(int angle)
{
	normalizeAngle(&angle);
	if (angle != zRot_)
	{
		zRot_ = angle;
		updateGL();
	}
}	

void Spectrum3DOpenGLCanvas::normalizeAngle(int *angle)
{
	while (*angle < 0)
	{
		*angle += 360 * 16;
	}
	while (*angle > 360 * 16)
	{
		*angle -= 360 * 16;
	}
}

void Spectrum3DOpenGLCanvas::setZoomFactor(double zoom)
{
	zoom_ = zoom;
	initializeGL();
	updateGL();
}

///////////////wheel- and MouseEvents//////////////////
void Spectrum3DOpenGLCanvas::wheelEvent ( QWheelEvent * e )
{
	if(view_mode_!= VIEW_ZOOM)
		{
			int wheel = e->delta();
			double distance = double(wheel/480.0);
			double zoom = zoom_+distance;
			if(zoom>0.0)
				{	
					setZoomFactor( zoom);
				}
			else
				{
					zoom = 0.25;
					setZoomFactor( zoom);
				}
		}
}

void Spectrum3DOpenGLCanvas::mouseMoveEvent ( QMouseEvent * e)
{
	if(view_mode_ == VIEW_TOP)
	{
		view_mode_ = VIEW_SELECT;
		initializeGL();
	}
	if(view_mode_!= VIEW_ZOOM)
		{
			int d_x = e->x() - lastMousePos_.x();
			int d_y = e->y() - lastMousePos_.y();
			setRotationX(xRot_ + 8 * d_y);
			setRotationY(yRot_ + 8 * d_x);
			lastMousePos_ = e->pos();
		}
  
	if(view_mode_== VIEW_ZOOM)
		{
 			lastMousePos_ = e->pos();
 			x_1_ = ((firstMousePos_.x()- width_/2) * corner_ * 2) / width_;
 			y_1_ = -300 + (((firstMousePos_.y()-heigth_/2) * corner_* 2) / heigth_);
 			x_2_ = ((lastMousePos_.x()- width_/2) * corner_ * 2) / width_;
 			y_2_ = -300 + (((lastMousePos_.y()-heigth_/2) * corner_* 2) / heigth_);
			show_zoom_selection_ = true;
 			initializeGL();
 			updateGL();
		}
	
}
void Spectrum3DOpenGLCanvas::mousePressEvent ( QMouseEvent * e)
{
	firstMousePos_ = e->pos();
	lastMousePos_ = e->pos();
	if(view_mode_!= VIEW_ZOOM)
	{
		lastMousePos_ = e->pos();
	}
}

void Spectrum3DOpenGLCanvas::mouseReleaseEvent ( QMouseEvent * e)
{
	if(e->button()==Qt::RightButton)
		{
			emit rightButton(e->globalPos());
		}

	if(view_mode_== VIEW_ZOOM )
		{			
			dataToZoomArray(x_1_, y_1_, x_2_, y_2_);
			show_zoom_selection_ = false;
			initializeGL();
			updateGL();
			 
		}
}

void Spectrum3DOpenGLCanvas::dataToZoomArray(double x_1, double y_1, double x_2, double y_2)
{
 double scale_x1 = scaledInversRT(x_1+100.0);
 double scale_x2 = scaledInversRT(x_2+100.0);
 double scale_y1 = scaledInversMZ(-200-y_1);
 double scale_y2 = scaledInversMZ(-200-y_2);
 if(scale_x1<=scale_x2)
	{
		canvas_3d_.overall_data_range_.min_[0]= scale_x1;
		canvas_3d_.overall_data_range_.max_[0]= scale_x2;
	} 
 else
	 {
		 canvas_3d_.overall_data_range_.min_[0]= scale_x2;
		 canvas_3d_.overall_data_range_.max_[0]= scale_x1;
	 }
 if(scale_y1<=scale_y2)
	{
		canvas_3d_.overall_data_range_.min_[1]= scale_y1;
		canvas_3d_.overall_data_range_.max_[1]= scale_y2;
	} 
 else
	 {
		 canvas_3d_.overall_data_range_.min_[1]= scale_y2;
		 canvas_3d_.overall_data_range_.max_[1]= scale_y1;
	 } 
}
void Spectrum3DOpenGLCanvas::updateIntensityScale()
{
	int_scale_.min_[0]= canvas_3d_.overall_data_range_.max_[2];
	int_scale_.max_[0]= canvas_3d_.overall_data_range_.min_[2];
	
		for(UnsignedInt i =0;i<canvas_3d_.getDataSetCount();i++)
		{
				for (Spectrum3DCanvas::ExperimentType::Iterator spec_it = canvas_3d_.getDataSet(i).RTBegin(canvas_3d_.overall_data_range_.min_[0]); 
						 spec_it != canvas_3d_.getDataSet(i).RTEnd(canvas_3d_.overall_data_range_.max_[0]); 
						 ++spec_it)
				{
					for (BaseSpectrum::Iterator it = spec_it->MZBegin(canvas_3d_.overall_data_range_.min_[1]); it!=spec_it->MZEnd(canvas_3d_.overall_data_range_.max_[1]); ++it)
					{
						if(	int_scale_.min_[0]>= it->getIntensity())
						{
							int_scale_.min_[0]= it->getIntensity(); 
						}
						if(	int_scale_.max_[0]<= it->getIntensity())
							{
							int_scale_.max_[0]= it->getIntensity(); 
						} 
					}
				}
		}
}

// sets the Multigradient gradient_
void Spectrum3DOpenGLCanvas::setDotGradient(const std::string& gradient)
{	
	gradient_.fromString(gradient);
	initializeGL();
	updateGL();

}
// recalculates the gradient of dataset number i 
void Spectrum3DOpenGLCanvas::recalculateDotGradient_(UnsignedInt i)
{
	canvas_3d_.getDataSet(i).updateRanges(1);
	gradient_.activatePrecalculationMode(canvas_3d_.getDataSet(i).getMinInt(),
																						canvas_3d_.getDataSet(i).getMaxInt(), 
																						UnsignedInt(canvas_3d_.getPref("Preferences:3D:Dot:InterpolationSteps")));
}
// recalculates the log.gradient of dataset number i 
void Spectrum3DOpenGLCanvas::recalculateDotGradientLog_(UnsignedInt i)
{
	canvas_3d_.getDataSet(i).updateRanges(1);
	gradient_.activatePrecalculationMode(log10(canvas_3d_.getDataSet(i).getMinInt()),
																						log10(canvas_3d_.getDataSet(i).getMaxInt()), 
																						UnsignedInt(canvas_3d_.getPref("Preferences:3D:Dot:InterpolationSteps")));
																						
}


///view front toolbar
void Spectrum3DOpenGLCanvas::setBackView()
{
	if(intensity_scale_)
	{
		intensity_scale_ = false;
		initializeGL();
		updateGL();
	}
	
	view_mode_ = VIEW_SELECT;
	
	if(showbackview_)
		{
			xRot_ = xRot_old_;
			yRot_ = yRot_old_;
			zRot_ = zRot_old_;
			initializeGL();
			updateGL();
		}
	
	showbackview_=false;
}

void Spectrum3DOpenGLCanvas::setTopView()
{
	if(intensity_scale_)
		{
			intensity_scale_ = false;
		}
	view_mode_ = VIEW_TOP;
	xRot_old_ = xRot_;
	yRot_old_ = yRot_;
	zRot_old_ = zRot_;	
	xRot_ = 90*16;
	yRot_ = 0;
	zRot_ = 0;	
	showbackview_ = true;
	initializeGL();
 	updateGL();

}

void Spectrum3DOpenGLCanvas::setResetZoomView()
{if(intensity_scale_)
		{
			intensity_scale_ = false;
		}
	canvas_3d_.recalculateRanges_(1,0,2);
	view_mode_ = VIEW_SELECT;
	xRot_ = 0;
	yRot_ = 0;
	zRot_ = 0;
	zoom_ = 1.25;
 	initializeGL();
 	updateGL();
	
}
void Spectrum3DOpenGLCanvas::setIntensityScale(bool on)
{
	intensity_scale_ = on;
	updateIntensityScale();
	view_mode_ = VIEW_SELECT;
	xRot_ = 0;
	yRot_ = 0;
	zRot_ = 0;
	zoom_ = 1.25;
 	initializeGL();
 	updateGL();
}
void Spectrum3DOpenGLCanvas::setZoomView()
{
	if(intensity_scale_)
		{
			intensity_scale_ = false;
		}
	view_mode_ = VIEW_ZOOM;
	zoom_ = 1.0;
	//saving the old angels
	xRot_old_ = xRot_;
	yRot_old_ = yRot_;
	zRot_old_ = zRot_;		
	// setting the topview angels
	xRot_ = 90*16;
	yRot_ = 0;
	zRot_ = 0;
	initializeGL();
	updateGL();
}

void Spectrum3DOpenGLCanvas::setSelectView()
{
	if(intensity_scale_)
	{
			intensity_scale_ = false;
	}
	zoom_ = 1.25;
	view_mode_ = VIEW_SELECT;
	x_1_ = 0.0;
	y_1_ = 0.0;
	x_2_ = 0.0;
	y_2_ = 0.0;
	initializeGL();
	xRot_ = 0;
	yRot_ = 0;
	zRot_ = 0;
	zoom_ = 1.25;
	resizeGL(int(width_),int( heigth_));
	updateGL();
}

bool Spectrum3DOpenGLCanvas::getShowSelect()
{
	if(view_mode_==VIEW_SELECT)
		{
			return true;
		}
	else
		{
			return false;
	}
}

bool Spectrum3DOpenGLCanvas::getShowZoom()
{
	if(view_mode_==VIEW_ZOOM)
		{
			return true;
		}
	else
		{
			return false;
	}
}
}//end of namespace
