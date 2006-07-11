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
// $Maintainer: Cornelia Friedle $
// --------------------------------------------------------------------------

//Qt
#include<qcolor.h>
//Open_GL
#include<qgl.h>
#include<GL/glu.h>
//STL
#include<iostream.h>
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
	setFocusPolicy(QWidget::StrongFocus);
	corner_=100.0;	
	near_=0.0;	
	far_=600.0;
	zoom_=1.25;	
	x_1_=0.0;
	y_1_=0.0;
	x_2_=0.0;
	y_2_=0.0;	
	xrot_=0;
	yrot_ = 0;
	zrot_=0;
	translation_on_ = false;
	trans_x_ =0.0;
	trans_y_ = 0.0;
	setSizePolicy(QSizePolicy::MinimumExpanding,QSizePolicy::MinimumExpanding);
}
	
Spectrum3DOpenGLCanvas::~Spectrum3DOpenGLCanvas()
{
 	grid_rt_.erase(	grid_rt_.begin(),	grid_rt_.end());
	grid_mz_.erase(grid_mz_.begin(),grid_mz_.end());
	grid_intensity_.erase(grid_intensity_.begin(),grid_intensity_.end());
}

void Spectrum3DOpenGLCanvas::calculateGridLines_()
{
	switch(canvas_3d_.intensity_mode_)
	{
	case SpectrumCanvas::IM_SNAP:
		updateIntensityScale();
		grid_intensity_=  AxisTickCalculator::calcGridLines_(int_scale_.min_[0],int_scale_.max_[0],3); 
	
		break;
		
	case SpectrumCanvas::IM_NONE:
		grid_intensity_=  AxisTickCalculator::calcGridLines_(canvas_3d_.overall_data_range_.min_[2],canvas_3d_.overall_data_range_.max_[2],3); 
		break;
		
	case SpectrumCanvas::IM_LOG:
		grid_intensity_log_=  AxisTickCalculator::calcLogGridLines_(log10(canvas_3d_.overall_data_range_.min_[2]),log10(canvas_3d_.overall_data_range_.max_[2]));	
 		break;
		
	case SpectrumCanvas::IM_PERCENTAGE:
	grid_intensity_=  AxisTickCalculator::calcGridLines_(0.0,100.0,3); 
	
			break;
	}
	grid_rt_=  AxisTickCalculator::calcGridLines_(canvas_3d_.visible_area_.min_[0],canvas_3d_.visible_area_.max_[0],3);
	grid_mz_=  AxisTickCalculator::calcGridLines_(canvas_3d_.visible_area_.min_[1],canvas_3d_.visible_area_.max_[1],3);
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
	gluLookAt(0.0, 0.0, 0.0, 
						0.0, 0.0, 0.0, 
						0.0, 1.0, 0.0);
	if(canvas_3d_.getDataSetCount()!=0 && canvas_3d_.recalculate_ )
	{
		switch(canvas_3d_.action_mode_)
			{
			case SpectrumCanvas::AM_ZOOM:

				calculateGridLines_();
				coord_ = makeCoordinates();
				if(canvas_3d_.show_grid_)
				{
					gridlines_ = makeGridLines();
				}
				xrot_ = 90*16;
				yrot_ = 0;
				zrot_ = 0;		
				stickdata_ = makeDataAsTopView();
				axeslabel_ = makeAxesLabel();	
				if(show_zoom_selection_)
				{
					zoomselection_ = makeZoomSelection();
				}
				break;
			case SpectrumCanvas::AM_TRANSLATE:
				calculateGridLines_();
				if(canvas_3d_.show_grid_)
				{
					gridlines_ = makeGridLines();
				}
				coord_ = makeCoordinates();
				ground_ = makeGround();
				x_1_ = 0.0;
				y_1_ = 0.0;
				x_2_ = 0.0;
				y_2_ = 0.0;
				
				if(canvas_3d_.intensity_mode_== SpectrumCanvas::IM_LOG)
				{
					stickdata_ = makeDataAsStickLog();
				}
				else
				{ 
					if(canvas_3d_.intensity_mode_ == SpectrumCanvas::IM_NONE || 
						 canvas_3d_.intensity_mode_ == SpectrumCanvas::IM_SNAP||
						 canvas_3d_.intensity_mode_ == SpectrumCanvas::IM_PERCENTAGE )
					{
						stickdata_ =  makeDataAsStick();
					}
				}
				axeslabel_ = makeAxesLabel();
				break;
			case SpectrumCanvas::AM_MEASURE:
				break;
 			case SpectrumCanvas::AM_SELECT:
				break;
			}
		}
}
void Spectrum3DOpenGLCanvas::resetAngels()
{
	xrot_=0;
	yrot_=0;
	zrot_=0;

}
void Spectrum3DOpenGLCanvas::resetTranslation()
{
	trans_x_ = 0.0;
	trans_y_ = 0.0;
}
void Spectrum3DOpenGLCanvas::paintGL()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	glTranslated(0.0, 0.0,-3.0*corner_);
	glRotated(xrot_ / 16.0, 1.0, 0.0, 0.0);
	glRotated(yrot_ / 16.0, 0.0, 1.0, 0.0);
	glRotated(zrot_/16.0, 0.0, 0.0, 1.0);
	glTranslated(0.0, 0.0,3.0*corner_);
	if(translation_on_)
	{
		glTranslated(trans_x_, trans_y_,0.0);
	}
	
	QColor color(canvas_3d_.getPrefAsString("Preferences:3D:BackgroundColor").c_str());
 	qglClearColor(color);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	if(canvas_3d_.getDataSetCount()!=0)
	{
		switch (canvas_3d_.action_mode_)
		{
	
		case SpectrumCanvas::AM_ZOOM:
			glCallList(stickdata_);	
			if(show_zoom_selection_)
			{
				glCallList(zoomselection_);
			}
			if(grid_rt_.size() !=0 )
			{	
				paintAxesScale();
			}
			if(canvas_3d_.show_grid_)
			{
				glCallList(gridlines_);
			}
			glDisable(GL_DEPTH_TEST);
			glCallList(coord_);
			glEnable(GL_DEPTH_TEST);
			break;
		
		case SpectrumCanvas::AM_TRANSLATE:	
			glCallList(ground_);
			glCallList(stickdata_);	
			glCallList(axeslabel_);
			if(grid_intensity_.size() !=0 )
			{
				paintAxesScale();
			}
			if(canvas_3d_.show_grid_)
			{
				glCallList(gridlines_);
			}
			glDisable(GL_DEPTH_TEST);
			glCallList(coord_);
			glEnable(GL_DEPTH_TEST);
			break;
			case SpectrumCanvas::AM_MEASURE:
				break;
		case SpectrumCanvas::AM_SELECT:
			break;
		}
	}
	
}

void Spectrum3DOpenGLCanvas::paintAxesScale()
{
	GLfloat black[3] = { 0.0, 0.0, 0.0 };
	glColor3fv(black);	//rt-axe
	if(yrot_> 280*16 ^ yrot_<80*16)
		{
			QString result("rt");
			renderText (0.0, -corner_-20.0, -near_-2*corner_+20.0,result);
			if(zoom_<2 && grid_rt_.size()>=2)
			{
				for(UnsignedInt i = 0;i<grid_rt_[0].size();i++)
				{
					QString result = QString("%1").arg(grid_rt_[0][i]);
					renderText (-corner_-result.length()+scaledRT(grid_rt_[0][i]), 
											-corner_-9.0,
											-near_-2*corner_+10.0,
											result);
				}
				for(UnsignedInt i = 0;i<grid_rt_[1].size();i++)
				{
					QString result = QString("%1").arg(grid_rt_[1][i]);
					renderText (-corner_-result.length()+scaledRT(grid_rt_[1][i]), 
											-corner_-9.0,
											-near_-2*corner_+10.0,
											result);
				}
			}
			if(width_>1000 && heigth_>700&& zoom_<1.5 && grid_rt_.size()>=3)
			{
				for(UnsignedInt i = 0;i<grid_rt_[2].size();i++)
					{
						QString result = QString("%1").arg(grid_rt_[2][i]);
						renderText (-corner_-result.length()+scaledRT(grid_rt_[2][i]), 
												-corner_-9.0,
												-near_-2*corner_+10.0,
												result);
					}
			}
		}
	if(yrot_>10*16 && yrot_<190*16 || canvas_3d_.action_mode_ == SpectrumCanvas::AM_ZOOM)
		{
			QString result("mz");
			renderText (-corner_-20.0, 
									-corner_-20.0,
									-near_-3*corner_,
									result);
			
			if(zoom_<2 && grid_mz_.size()>=2)
			{
				for(UnsignedInt i = 0;i<grid_mz_[0].size();i++)
				{
					QString result = QString("%1").arg(grid_mz_[0][i]);
					renderText (-corner_-result.length()-10.0, 
											-corner_-9.0,
											-near_-2*corner_-scaledMZ(grid_mz_[0][i]),
											result);
				}
				for(UnsignedInt i = 0;i<grid_mz_[1].size();i++)
				{
					QString result = QString("%1").arg(grid_mz_[1][i]);
					renderText (-corner_-result.length()-10.0, 
											-corner_-9.0,
											-near_-2*corner_-scaledMZ(grid_mz_[1][i]),
											result);
				}
			}
			if(width_>1000 && heigth_>700 && zoom_<1.5 && grid_mz_.size()>=3)
			{
				for(UnsignedInt i = 0;i<grid_mz_[2].size();i++)
				{
					QString result = QString("%1").arg(grid_mz_[2][i]);
					renderText (-corner_-result.length()-10.0, 
											-corner_-9.0,
											-near_-2*corner_-scaledMZ(grid_mz_[2][i]),
											result);
				}
			}
		}
	
	if(canvas_3d_.action_mode_ != SpectrumCanvas::AM_ZOOM)
		{
					QString result("intensity");
				renderText (-corner_-20.0, 
										corner_+10.0,
										-near_-2*corner_+20.0,
										result);
	
			if(canvas_3d_.intensity_mode_== SpectrumCanvas::IM_LOG)
				{
					if(zoom_<3 && grid_intensity_.size()>=1)
						{
							for(UnsignedInt i = 0;i<grid_intensity_log_[0].size();i++)
								{
									QString result = QString("%1").arg(grid_intensity_log_[0][i]);
									renderText (-corner_-result.length()-3.0, 
															-corner_+scaledIntensity(grid_intensity_log_[0][i]),
															-near_-2*corner_,
															result);
								}
						}
				}
			else
				{
					if(canvas_3d_.intensity_mode_ == SpectrumCanvas::IM_NONE || 
						 canvas_3d_.intensity_mode_ == SpectrumCanvas::IM_SNAP ||
						 canvas_3d_.intensity_mode_ == SpectrumCanvas::IM_PERCENTAGE)
						{
							if(zoom_<2.0 && width_>=800 && 
								 grid_intensity_.size()>=2)
								{
									for(UnsignedInt i = 0;i<grid_intensity_[0].size();i++)
										{
											QString result = QString("%1").arg(grid_intensity_[0][i]);
											renderText (-corner_-result.length()-width_/200.0-5.0, 
																	-corner_+scaledIntensity(grid_intensity_[0][i]),
																	-near_-2*corner_,
																	result);
											



										}
									for(UnsignedInt i = 0;i<grid_intensity_[1].size();i++)
										{
											QString result = QString("%1").arg(grid_intensity_[1][i]);
											renderText (-corner_-result.length()-width_/200.0-5.0, 
																	-corner_+scaledIntensity(grid_intensity_[1][i]),
																	-near_-2*corner_,
																	result);
										}
								}
							if(width_>1000 && heigth_>700&& zoom_<1.5 && grid_intensity_.size()>=3)
							{
								for(UnsignedInt i = 0;i<grid_intensity_[2].size();i++)
									{
										QString result = QString("%1").arg(grid_intensity_[2][i]);
										renderText (-corner_-result.length()-width_/200.0-5.0, 
																-corner_+scaledIntensity(grid_intensity_[2][i]),
																-near_-2*corner_,
																	result);
									}
								}
						}
				}
		}
}
GLuint Spectrum3DOpenGLCanvas::makeGround()
{
	GLuint list = glGenLists(1);
	glNewList(list,GL_COMPILE);
 	QColor color(canvas_3d_.getPrefAsString("Preferences:3D:BackgroundColor").c_str());
	qglColor(color);
	glBegin(GL_QUADS);
	glVertex3d(-corner_, 
						 -corner_-2.0,
						 -near_-2*corner_);
	glVertex3d( -corner_, 
							-corner_-2.0,
							-far_+2*corner_);
	glVertex3d( corner_, 
							-corner_-2.0,
							-far_+2*corner_);
	glVertex3d(corner_, 
						 -corner_-2.0,
						 -near_-2*corner_);
	glEnd();
	glEndList();
	return list;
}


GLuint Spectrum3DOpenGLCanvas::makeZoomSelection()
{
	GLuint list = glGenLists(1);
	glNewList(list,GL_COMPILE);
	//	glEnable(GL_BLEND);
	//	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
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
	glVertex3d(-corner_, 
						 -corner_,
						 -near_-2*corner_);
	glVertex3d( corner_, 
							-corner_,
							-near_-2*corner_);
	//z-achse
	glVertex3d(-corner_, 
						 -corner_,
						 -near_-2*corner_);
	glVertex3d( -corner_, 
							-corner_,
							-far_+2*corner_);
	//y-achse
	glVertex3d(-corner_, 
						 -corner_,
						 -near_-2*corner_);
	glVertex3d( -corner_, 
							corner_,
							-near_-2*corner_);
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
 	//recalculateDotGradient_();		 
	for(UnsignedInt i =0;i<canvas_3d_.getDataSetCount();i++)
	{
		if(canvas_3d_.layer_visible_[i]==true)
		{	
			for (Spectrum3DCanvas::ExperimentType::ConstIterator spec_it = canvas_3d_.getDataSet(i).RTBegin(canvas_3d_.visible_area_.min_[0]); 
					 spec_it != canvas_3d_.getDataSet(i).RTEnd(canvas_3d_.visible_area_.max_[0]); 
					 ++spec_it)
			{
				if (spec_it->getMSLevel()!=1)
				{
					continue;
				}
				for (BaseSpectrum::ConstIterator it = spec_it->MZBegin(canvas_3d_.visible_area_.min_[1]); it!=spec_it->MZEnd(canvas_3d_.visible_area_.max_[1]); ++it)
				{	
					if(it->getIntensity()>= canvas_3d_.disp_ints_[i].first && it->getIntensity()<= canvas_3d_.disp_ints_[i].second)
					{
						glBegin(GL_POINTS);
						if(int(canvas_3d_.getPref("Preferences:3D:Dot:Mode")))
						{
							qglColor(QColor( gradient_.precalculatedColorAt(it->getIntensity())));
							glVertex3d(-corner_+(GLfloat)scaledRT(spec_it->getRetentionTime()), 
												 -corner_,
												 -near_-2*corner_-(GLfloat)scaledMZ(it->getPosition()[0]));
						}
						else
						{
							qglColor(black);			
							glVertex3d(-corner_+(GLfloat)scaledRT(spec_it->getRetentionTime()), 
												 -corner_,
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
	//recalculateDotGradient_();		 
	for(UnsignedInt i =0;i<canvas_3d_.getDataSetCount();i++)
	{
		if(canvas_3d_.isDataSetVisible(i))
		{	
			for (Spectrum3DCanvas::ExperimentType::ConstIterator spec_it = canvas_3d_.getDataSet(i).RTBegin(canvas_3d_.visible_area_.min_[0]); 
					 spec_it != canvas_3d_.getDataSet(i).RTEnd(canvas_3d_.visible_area_.max_[0]); 
					 ++spec_it)
			{
				if (spec_it->getMSLevel()!=1)
				{
					continue;
				}
				for (BaseSpectrum::ConstIterator it = spec_it->MZBegin(canvas_3d_.visible_area_.min_[1]); it!=spec_it->MZEnd(canvas_3d_.visible_area_.max_[1]); ++it)
				{			
					if(it->getIntensity()>= canvas_3d_.disp_ints_[i].first && it->getIntensity()<= canvas_3d_.disp_ints_[i].second)
					{
						glBegin(GL_LINES);

						if(int(canvas_3d_.getPref("Preferences:3D:Dot:Mode")))
						{
							if(canvas_3d_.intensity_mode_ == SpectrumCanvas::IM_PERCENTAGE)
							{
								double intensity = (it->getIntensity() * 100.0) / canvas_3d_.overall_data_range_.max_[2];
								qglColor(QColor( gradient_.precalculatedColorAt(0)));
								glVertex3d(-corner_+(GLfloat)scaledRT(spec_it->getRetentionTime()), 
														 -corner_,
													 -near_-2*corner_-(GLfloat)scaledMZ(it->getPosition()[0]));
								
								qglColor(QColor( gradient_.precalculatedColorAt(intensity)));
								glVertex3d(-corner_+(GLfloat)scaledRT(spec_it->getRetentionTime()),
													 -corner_+(GLfloat)scaledIntensity(intensity),
														 -near_-2*corner_-(GLfloat)scaledMZ(it->getPosition()[0]));
							}
							else
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
						}
						else
						{
							qglColor(black);	
							if(canvas_3d_.intensity_mode_ == SpectrumCanvas::IM_PERCENTAGE)
							{
								double intensity =  it->getIntensity()  * 100.0 / canvas_3d_.overall_data_range_.max_[2];
								glVertex3d(-corner_+(GLfloat)scaledRT(spec_it->getRetentionTime()), 
													 -corner_,
													 -near_-2*corner_-(GLfloat)scaledMZ(it->getPosition()[0]));
								glVertex3d(-corner_+(GLfloat)scaledRT(spec_it->getRetentionTime()),
													 -corner_+(GLfloat)scaledIntensity(intensity),
													 -near_-2*corner_-(GLfloat)scaledMZ(it->getPosition()[0]));
								
							}	
							else
							{
								glVertex3d(-corner_+(GLfloat)scaledRT(spec_it->getRetentionTime()), 
													 -corner_,
													 -near_-2*corner_-(GLfloat)scaledMZ(it->getPosition()[0]));
								glVertex3d(-corner_+(GLfloat)scaledRT(spec_it->getRetentionTime()),
															 -corner_+(GLfloat)scaledIntensity(it->getIntensity()),
													 -near_-2*corner_-(GLfloat)scaledMZ(it->getPosition()[0]));
							}
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
	
	//recalculateDotGradient_();
	for(UnsignedInt i =0;i<canvas_3d_.getDataSetCount();i++)
		{
			if(canvas_3d_.isDataSetVisible(i))
			{		
				for (Spectrum3DCanvas::ExperimentType::ConstIterator spec_it = canvas_3d_.getDataSet(i).RTBegin(canvas_3d_.visible_area_.min_[0]); 
						 spec_it != canvas_3d_.getDataSet(i).RTEnd(canvas_3d_.visible_area_.max_[0]); 
						 ++spec_it)
				{
					if (spec_it->getMSLevel()!=1)
					{
						continue;
					}
					for (BaseSpectrum::ConstIterator it = spec_it->MZBegin(canvas_3d_.visible_area_.min_[1]); it!=spec_it->MZEnd(canvas_3d_.visible_area_.max_[1]); ++it)
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
													 -corner_+(GLfloat)scaledIntensity(log10(it->getIntensity())),
													 -near_-2*corner_-(GLfloat)scaledMZ(it->getPosition()[0]));
							}
							else
							{
								qglColor(black);			
								glVertex3d(-corner_+(GLfloat)scaledRT(spec_it->getRetentionTime()), 
													 -corner_,
													 -near_-2*corner_-(GLfloat)scaledMZ(it->getPosition()[0]));
								glVertex3d(-corner_+(GLfloat)scaledRT(spec_it->getRetentionTime()),
													 -corner_+(GLfloat)scaledIntensity(log10(it->getIntensity())),
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
GLuint Spectrum3DOpenGLCanvas::makeGridLines()
{
	GLuint list = glGenLists(1);
	glNewList(list, GL_COMPILE);
	glBegin(GL_LINES);
	glColor4ub(0, 0, 0, 80);	
	//rt
	if(grid_rt_.size()>=1)
		{
	for(UnsignedInt i = 0;i<grid_rt_[0].size();i++)
	{
		glVertex3d(-corner_+scaledRT(grid_rt_[0][i]), 
							 -corner_,
							  -near_-2*corner_);
		glVertex3d( -corner_+scaledRT(grid_rt_[0][i]), 
								-corner_,
								-far_+2*corner_);
	}
		}
	if(grid_rt_.size()>=2)
		{
			for(UnsignedInt i = 0;i<grid_rt_[1].size();i++)
				{	
					glVertex3d(-corner_+scaledRT(grid_rt_[1][i]), 
										 -corner_,
										 -near_-2*corner_);
					glVertex3d( -corner_+scaledRT(grid_rt_[1][i]), 
								-corner_,
											-far_+2*corner_);
					
				}
		}
	if(grid_rt_.size()>=3)
		{
			for(UnsignedInt i = 0;i<grid_rt_[2].size();i++)
				{	
					glVertex3d(-corner_+scaledRT(grid_rt_[2][i]), 
										 -corner_,
										 -near_-2*corner_);
					glVertex3d( -corner_+scaledRT(grid_rt_[2][i]), 
											-corner_,	
											-far_+2*corner_);
				}
		}
	//mz
	if(grid_mz_.size()>=1)
		{
			for(UnsignedInt i = 0;i<grid_mz_[0].size();i++)
				{
					glVertex3d(-corner_, 
										 -corner_,
										 -near_-2*corner_-scaledMZ(grid_mz_[0][i]));
					glVertex3d( corner_,
											-corner_,
											-near_-2*corner_-scaledMZ(grid_mz_[0][i]));
				}
		}
	if(grid_mz_.size()>=2)
	{
		for(UnsignedInt i = 0;i<grid_mz_[1].size();i++)
		{
			glVertex3d( -corner_, 
									-corner_,
									-near_-2*corner_-scaledMZ(grid_mz_[1][i]));
			glVertex3d( corner_,
									-corner_,
								-near_-2*corner_-scaledMZ(grid_mz_[1][i]));
		}
	}
	if(grid_mz_.size()>=3)
	{
		for(UnsignedInt i = 0;i<grid_mz_[2].size();i++)
		{
			glVertex3d( -corner_, 
									-corner_,
									-near_-2*corner_-scaledMZ(grid_mz_[2][i]));
			glVertex3d( corner_,
									-corner_,
									-near_-2*corner_-scaledMZ(grid_mz_[2][i]));
		}
	}
	glEnd();
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
	if(grid_rt_.size()>=1)
	{
		for(UnsignedInt i = 0;i<grid_rt_[0].size();i++)
		{
			glVertex3d(-corner_+scaledRT(grid_rt_[0][i]), 
								 -corner_,
								 -near_-2*corner_);
			glVertex3d( -corner_+scaledRT(grid_rt_[0][i]), 
									-corner_+4.0,
									-near_-2*corner_);
			}
	}
	if(grid_rt_.size()>=2)
	{
		for(UnsignedInt i = 0;i<grid_rt_[1].size();i++)
		{
		glVertex3d(-corner_+scaledRT(grid_rt_[1][i]), 
							 -corner_,
							 -near_-2*corner_);
		glVertex3d( -corner_+scaledRT(grid_rt_[1][i]), 
								-corner_+3.0,
								-near_-2*corner_);
		}
	}
	if(grid_rt_.size()>=3)
	{
		for(UnsignedInt i = 0;i<grid_rt_[2].size();i++)
		{
			glVertex3d(-corner_+scaledRT(grid_rt_[2][i]), 
								 -corner_,
								 -near_-2*corner_);
			glVertex3d( -corner_+scaledRT(grid_rt_[2][i]), 
									-corner_+2.0,
									-near_-2*corner_);
		}
	}
	
	//MZ	
	if(grid_mz_.size()>=1)
	{
		for(UnsignedInt i = 0;i<grid_mz_[0].size();i++)
			{
				glVertex3d(-corner_, 
									 -corner_,
									 -near_-2*corner_-scaledMZ(grid_mz_[0][i]));
				glVertex3d( -corner_,
										-corner_+4.0,
										-near_-2*corner_-scaledMZ(grid_mz_[0][i]));
			}
	}
	if(grid_mz_.size()>=2)
	{
			for(UnsignedInt i = 0;i<grid_mz_[1].size();i++)
		{
			glVertex3d(-corner_, 
								 -corner_,
							 -near_-2*corner_-scaledMZ(grid_mz_[1][i]));
			glVertex3d( -corner_,
									-corner_+3.0,
									-near_-2*corner_-scaledMZ(grid_mz_[1][i]));
		}
	}	
	if(grid_mz_.size()>=3)
	{
		for(UnsignedInt i = 0;i<grid_mz_[2].size();i++)
		{
			glVertex3d(-corner_, 
								 -corner_,
								 -near_-2*corner_-scaledMZ(grid_mz_[2][i]));
			glVertex3d( -corner_,
									-corner_+2.0,
									-near_-2*corner_-scaledMZ(grid_mz_[2][i]));
		}
	}



	if(canvas_3d_.intensity_mode_== SpectrumCanvas::IM_LOG)
		{
			if(grid_intensity_.size()>=1)
			{
				for(UnsignedInt i = 0;i<grid_intensity_log_[0].size();i++)
				{	
					glVertex3d(-corner_, 
										 -corner_+scaledIntensity(grid_intensity_log_[0][i]),
										 -near_-2*corner_);
					glVertex3d( -corner_+3.0, 
											-corner_+scaledIntensity(grid_intensity_log_[0][i]),
											-near_-2*corner_-3.0);
				}
				for(UnsignedInt i = 0;i<grid_intensity_log_[1].size();i++)
				{
					glVertex3d(-corner_, 
										 -corner_+scaledIntensity(grid_intensity_log_[1][i]),
										 -near_-2*corner_);
					glVertex3d( -corner_+2.0, 
											-corner_+scaledIntensity(grid_intensity_log_[1][i]),
											-near_-2*corner_-2.0);
				}
			}
		}
	else
	{
		if(canvas_3d_.intensity_mode_== SpectrumCanvas::IM_NONE||
			 canvas_3d_.intensity_mode_== SpectrumCanvas::IM_SNAP||
			 canvas_3d_.intensity_mode_== SpectrumCanvas::IM_PERCENTAGE)
		{
			if(grid_intensity_.size()>=1)
				{
					for(UnsignedInt i = 0;i<grid_intensity_[0].size();i++)
						{	
							glVertex3d(-corner_, 
												 -corner_+scaledIntensity(grid_intensity_[0][i]),
														 -near_-2*corner_);
							glVertex3d( -corner_+4.0, 
															-corner_+scaledIntensity(grid_intensity_[0][i]),
															-near_-2*corner_-4.0);
								}
						}
						if(grid_intensity_.size()>=2)
						{
							for(UnsignedInt i = 0;i<grid_intensity_[1].size();i++)
								{
									glVertex3d(-corner_, 
														 -corner_+scaledIntensity(grid_intensity_[1][i]),
														 -near_-2*corner_);
									glVertex3d( -corner_+3.0, 
															-corner_+scaledIntensity(grid_intensity_[1][i]),
															-near_-2*corner_-3.0);
								}
						}
						if(grid_intensity_.size()>=3)
							{
								for(UnsignedInt i = 0;i<grid_intensity_[2].size();i++)
									{ 
										glVertex3d(-corner_, 
															 -corner_+scaledIntensity(grid_intensity_[2][i]),
															 -near_-2*corner_);
										glVertex3d( -corner_+2.0, 
													-corner_+scaledIntensity(grid_intensity_[2][i]),
																-near_-2*corner_-2.0);
									}
							}
				} 
		}
	glEnd();
	glEndList();
	return list; 
}

double Spectrum3DOpenGLCanvas::scaledRT(double rt)
{
	double scaledrt = rt - canvas_3d_.visible_area_.min_[0];
	scaledrt = scaledrt * 2.0 * corner_ /(canvas_3d_.visible_area_.max_[0]-canvas_3d_.visible_area_.min_[0]);
	return scaledrt;
}

double Spectrum3DOpenGLCanvas::scaledInversRT(double rt)
{
	double i_rt =(rt* canvas_3d_.visible_area_.max_[0] -canvas_3d_.visible_area_.min_[0]*rt);
	i_rt = i_rt/200.0;
	i_rt = i_rt + canvas_3d_.visible_area_.min_[0]; 
	return i_rt;
 }

double Spectrum3DOpenGLCanvas::scaledMZ(double mz)
{
	double scaledmz = mz - canvas_3d_.visible_area_.min_[1];
	scaledmz = scaledmz * 2.0 * corner_/(canvas_3d_.visible_area_.max_[1]-canvas_3d_.visible_area_.min_[1])/*dis_mz_*/;
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
	switch(canvas_3d_.intensity_mode_)
	{
	case  SpectrumCanvas::IM_SNAP:
		scaledintensity = intensity -int_scale_.min_[0];
		scaledintensity = ( scaledintensity * 2.0 * corner_)/(int_scale_.max_[0]-int_scale_.min_[0]);
		break;
	case  SpectrumCanvas::IM_NONE:
		scaledintensity = intensity -canvas_3d_.overall_data_range_.min_[2];
		scaledintensity = ( scaledintensity * 2.0 * corner_)/(canvas_3d_.overall_data_range_.max_[2]-canvas_3d_.overall_data_range_.min_[2]);
		break;
	case  SpectrumCanvas::IM_PERCENTAGE:
		scaledintensity = (intensity *100.0)/100.0;
		scaledintensity = ( scaledintensity * 2.0 * corner_)/100.0;
 		break;	
	case SpectrumCanvas::IM_LOG:
		scaledintensity = intensity -log10(canvas_3d_.overall_data_range_.min_[2]);	
		scaledintensity =(scaledintensity * 2.0 * corner_)/(log10(canvas_3d_.overall_data_range_.max_[2])-log10(canvas_3d_.overall_data_range_.min_[2]));
		if(scaledintensity<0)
		{
			scaledintensity =0.0;
		}
		break;
	}
	return scaledintensity;
}


void Spectrum3DOpenGLCanvas::setRotationX(int angle)
{
	normalizeAngle(&angle);
	if (angle != xrot_) 
	{
		xrot_ = angle;
		updateGL();
	}
}

void Spectrum3DOpenGLCanvas::setRotationY(int angle)
{
	normalizeAngle(&angle);
	if (angle != yrot_) 
	{
		yrot_ = angle;
		updateGL();
	}
}

void Spectrum3DOpenGLCanvas::setRotationZ(int angle)
{
	normalizeAngle(&angle);
	if (angle != zrot_)
	{
		zrot_ = angle;
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
	if(canvas_3d_.action_mode_ == SpectrumCanvas::AM_TRANSLATE)
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
	if(canvas_3d_.action_mode_ == SpectrumCanvas::AM_TRANSLATE && !translation_on_)
	{
		int d_x = e->x() - lastMousePos_.x();
		int d_y = e->y() - lastMousePos_.y();
		setRotationX(xrot_ + 8 * d_y);
		setRotationY(yrot_ + 8 * d_x);
		lastMousePos_ = e->pos();
	}
	else
	{
		if(canvas_3d_.action_mode_ == SpectrumCanvas::AM_ZOOM)
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
	if(translation_on_)
	{
		lastMousePos_ = e->pos();
		trans_x_= lastMousePos_.x()-firstMousePos_.x();
		trans_y_ = (heigth_-lastMousePos_.y())-(heigth_ -firstMousePos_.y());
		canvas_3d_.recalculate_ = false;
		canvas_3d_.invalidate_();
	}
}
void Spectrum3DOpenGLCanvas::mousePressEvent ( QMouseEvent * e)
{
	firstMousePos_ = e->pos();
	lastMousePos_ = e->pos();
}

void Spectrum3DOpenGLCanvas::mouseReleaseEvent ( QMouseEvent * e)
{
	canvas_3d_.recalculate_ = true;
	translation_on_ = false;

	if(e->button()==Qt::RightButton)
		{
			emit rightButton(e->globalPos());
		}
	if(canvas_3d_.action_mode_ == SpectrumCanvas::AM_ZOOM)
		{			
			dataToZoomArray(x_1_, y_1_, x_2_, y_2_);
			show_zoom_selection_ = false;
			initializeGL();
			updateGL();
			 
		}
}
void Spectrum3DOpenGLCanvas::keyPressEvent(QKeyEvent * e) 
{
	if(e->key()==Qt::Key_Control)
	{
		translation_on_ = true;
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
		canvas_3d_.visible_area_.min_[0]= scale_x1;
		canvas_3d_.visible_area_.max_[0]= scale_x2;
	} 
 else
	 {
		 canvas_3d_.visible_area_.min_[0]= scale_x2;
		 canvas_3d_.visible_area_.max_[0]= scale_x1;
	 }
 if(scale_y1<=scale_y2)
	{
		canvas_3d_.visible_area_.min_[1]= scale_y1;
		canvas_3d_.visible_area_.max_[1]= scale_y2;
	} 
 else
	 {
		 canvas_3d_.visible_area_.min_[1]= scale_y2;
		 canvas_3d_.visible_area_.max_[1]= scale_y1;
	 } 
}
void Spectrum3DOpenGLCanvas::updateIntensityScale()
{
	int_scale_.min_[0]= canvas_3d_.overall_data_range_.max_[2];
	int_scale_.max_[0]= canvas_3d_.overall_data_range_.min_[2];
	
		for(UnsignedInt i =0;i<canvas_3d_.getDataSetCount();i++)
		{
				for (Spectrum3DCanvas::ExperimentType::ConstIterator spec_it = canvas_3d_.getDataSet(i).RTBegin(canvas_3d_.visible_area_.min_[0]); 
						 spec_it != canvas_3d_.getDataSet(i).RTEnd(canvas_3d_.visible_area_.max_[0]); 
						 ++spec_it)
				{
					for (BaseSpectrum::ConstIterator it = spec_it->MZBegin(canvas_3d_.visible_area_.min_[1]); it!=spec_it->MZEnd(canvas_3d_.visible_area_.max_[1]); ++it)
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
void Spectrum3DOpenGLCanvas::recalculateDotGradient_()
{
	switch(canvas_3d_.intensity_mode_)
	{
	case SpectrumCanvas::IM_SNAP:
		gradient_.activatePrecalculationMode(int_scale_.min_[0],
																				 int_scale_.max_[0], 
																				 UnsignedInt(canvas_3d_.getPref("Preferences:3D:Dot:InterpolationSteps")));
		break;
	case SpectrumCanvas::IM_NONE:
		gradient_.activatePrecalculationMode(canvas_3d_.overall_data_range_.min_[2],
																				 canvas_3d_.overall_data_range_.max_[2], 
																				 UnsignedInt(canvas_3d_.getPref("Preferences:3D:Dot:InterpolationSteps")));
		break;
	case SpectrumCanvas::IM_LOG:
		gradient_.activatePrecalculationMode(log10(canvas_3d_.overall_data_range_.min_[2]),
																						log10(canvas_3d_.overall_data_range_.max_[2]), 
																						UnsignedInt(canvas_3d_.getPref("Preferences:3D:Dot:InterpolationSteps")));
		
		break;
	case SpectrumCanvas::IM_PERCENTAGE:
		gradient_.activatePrecalculationMode(0.0,
																				 100.0,
																				 UnsignedInt(canvas_3d_.getPref("Preferences:3D:Dot:InterpolationSteps")));
		break;
	}
}

}//end of namespace
