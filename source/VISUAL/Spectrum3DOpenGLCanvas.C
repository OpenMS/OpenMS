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


#include <OpenMS/VISUAL/Spectrum3DOpenGLCanvas.h>
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/VISUAL/AxisTickCalculator.h>

#include <QtGui/QWheelEvent>
#include <QtGui/QMouseEvent>
#include <QtGui/QKeyEvent>

using std::cout;
using std::endl;

namespace OpenMS
{
	  
	Spectrum3DOpenGLCanvas::Spectrum3DOpenGLCanvas(QWidget *parent,Spectrum3DCanvas& canvas_3d)
	  : QGLWidget(parent),
	    canvas_3d_(canvas_3d),
	    zoom_mode_(false)
	{
		canvas_3d.rubber_band_.setParent(this);
	  setFocusPolicy(Qt::StrongFocus);
	  corner_=100.0;  
	  near_=0.0;  
	  far_=600.0;
	  zoom_= 1.5; 
	  xrot_=220;
	  yrot_ = 220;
	  zrot_=0;
	  trans_x_ =0.0;
	  trans_y_ = 0.0;
	  setSizePolicy(QSizePolicy::MinimumExpanding,QSizePolicy::MinimumExpanding);
	}
	  
	Spectrum3DOpenGLCanvas::~Spectrum3DOpenGLCanvas()
	{
	}
	
	void Spectrum3DOpenGLCanvas::calculateGridLines_()
	{
	  double dist;
	  switch(canvas_3d_.intensity_mode_)
	  {
		  case SpectrumCanvas::IM_SNAP:
		    updateIntensityScale();
		    AxisTickCalculator::calcGridLines(int_scale_.min_[0],int_scale_.max_[0],3,grid_intensity_,7,5,dist); 
		    
		    break;
		  case SpectrumCanvas::IM_NONE:
		    AxisTickCalculator::calcGridLines(canvas_3d_.overall_data_range_.min_[2],canvas_3d_.overall_data_range_.max_[2],3,grid_intensity_,7,5,dist); 
		    break;
		  case SpectrumCanvas::IM_LOG:
		    double log_min;
		    if(log10(canvas_3d_.overall_data_range_.min_[2])<0)
		    {
		      log_min = 0.0;
		    }
		    else
		    {
		      log_min = log10(canvas_3d_.overall_data_range_.min_[2]);
		    }
		    AxisTickCalculator::calcLogGridLines(log_min,log10(canvas_3d_.overall_data_range_.max_[2]), grid_intensity_log_); 
		    break;
		  case SpectrumCanvas::IM_PERCENTAGE:
		    AxisTickCalculator::calcGridLines(0.0,100.0,3,grid_intensity_,7,5,dist); 
		    break;
	  }
	  
	  //cout << endl;
	  //for (UInt level=0; level!=grid_intensity_.size(); ++level)
	  //{
	  //	cout << "Level: " << level << endl;
	  //	for (UInt tick=0; tick!=grid_intensity_[level].size(); ++tick)
	  //	{
	  //		cout << "  " << grid_intensity_[level][tick] << endl;
	  //	}
	  //}
	  
	  AxisTickCalculator::calcGridLines(canvas_3d_.visible_area_.min_[0],canvas_3d_.visible_area_.max_[0],3,grid_rt_,7,5,dist);
	  AxisTickCalculator::calcGridLines(canvas_3d_.visible_area_.min_[1],canvas_3d_.visible_area_.max_[1],3,grid_mz_,7,5,dist);
	}
	
	
	void Spectrum3DOpenGLCanvas::resizeGL(int w,int h)
	{
    width_ = (float)w;
    heigth_ = (float)h;
    glViewport(0,0,(GLsizei) w, (GLsizei) h);   
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-corner_*zoom_, corner_*zoom_,  -corner_*zoom_, corner_*zoom_ , near_, far_);
    glMatrixMode(GL_MODELVIEW);
	}
	
	void Spectrum3DOpenGLCanvas::initializeGL()
	{
	  QColor color(canvas_3d_.param_.getValue("background_color").toQString());
	  qglClearColor(color);
	  glEnable(GL_DEPTH_TEST);
	  glEnable(GL_BLEND);
	  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	  calculateGridLines_();
	  switch(canvas_3d_.action_mode_)
	  {
	  	case SpectrumCanvas::AM_ZOOM:
	      if(canvas_3d_.getLayerCount()!=0 )
	      {
					if (!zoom_mode_) //translate
					{
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
			      stickdata_ =  makeDataAsStick();
			      axeslabel_ = makeAxesLabel();
			      if(canvas_3d_.legend_shown_)
			      {
			        axeslegend_ = makeLegend();
			      }
					}
					else //zoom
					{
		        if(!canvas_3d_.rubber_band_.isVisible())
		        {
		          coord_ = makeCoordinates();
		          if(canvas_3d_.show_grid_)
		          {
		            gridlines_ = makeGridLines();
		          }
		          xrot_ = 90*16;
		          yrot_ = 0;
		          zrot_ = 0;  
		          zoom_ = 1.25;
		          stickdata_ = makeDataAsTopView();
		          axeslabel_ = makeAxesLabel();   
		          if(canvas_3d_.legend_shown_)
		          {
		            axeslegend_ = makeLegend();
		          }
		        }						
					}
	      }
	      break;
		  default:
		    break;
	  }
	  
	}
	
	void Spectrum3DOpenGLCanvas::setAngels(int xrot, int yrot, int zrot)
	{
		xrot_=xrot;
		yrot_ =yrot;
		zrot_=zrot;
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
		glTranslated(trans_x_, trans_y_,3.0*corner_);
		
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		if(canvas_3d_.getLayerCount()!=0)
		{
			switch (canvas_3d_.action_mode_)
			{
				case SpectrumCanvas::AM_ZOOM:
					if (!zoom_mode_) //translate
					{
						glCallList(ground_);
						glCallList(stickdata_);	
						glCallList(axeslabel_);
						if(canvas_3d_.legend_shown_)
						{
							glCallList(axeslegend_);	
						}
						if(canvas_3d_.show_grid_)
						{
							glCallList(gridlines_);
						}
						//glDisable(GL_DEPTH_TEST);
						glCallList(coord_);
						//glEnable(GL_DEPTH_TEST);
					}
					else //zoom
					{
						glCallList(stickdata_);	
						if(canvas_3d_.legend_shown_)
							{
								glCallList(axeslegend_);	
							}
						glCallList(axeslabel_);
						if(canvas_3d_.show_grid_)
						{	
							glEnable (GL_LINE_STIPPLE);	
							glCallList(gridlines_);
							glDisable (GL_LINE_STIPPLE);	
						}
						//glDisable(GL_DEPTH_TEST);
						glCallList(coord_);
						//glEnable(GL_DEPTH_TEST);
					}
					break;
				default:
					break;			
			}
		}
	}
	
	GLuint Spectrum3DOpenGLCanvas::makeLegend()
	{
		QFont font("Typewriter");
		font.setPixelSize(12);
 		renderText (0.0, 0.0, 0.0, "", font);
		
		GLuint list = glGenLists(1);
		glNewList(list,GL_COMPILE);
		QColor color(canvas_3d_.param_.getValue("axis_color").toQString());
 		qglColor(color);
		QString result("RT");
 		renderText (0.0, 
								-corner_-20.0, 
								-near_-2*corner_+20.0,
								result,
								font);
		font.setPixelSize(10);
		if(grid_rt_.size()>0)
		{
			for(UInt i = 0;i<grid_rt_[0].size();i++)
				{
					result = QString::number(grid_rt_[0][i]);
					renderText (-corner_-result.length()+scaledRT(grid_rt_[0][i]), 
											-corner_-5.0,
											-near_-2*corner_+15.0,
											result,
											font);
				}
		}
		if(zoom_<3.0 && grid_rt_.size()>=2)
			{
				for(UInt i = 0;i<grid_rt_[1].size();i++)
					{
						result = QString::number(grid_rt_[1][i]);
						renderText (-corner_-result.length()+scaledRT(grid_rt_[1][i]), 
												-corner_-5.0,
												-near_-2*corner_+15.0,
												result,
												font);
					}
			}
		if(zoom_<2.0 && grid_rt_.size()>=3)
			{
				for(UInt i = 0;i<grid_rt_[2].size();i++)
					{
						result = QString::number(grid_rt_[2][i]);
						renderText (-corner_-result.length()+scaledRT(grid_rt_[2][i]), 
												-corner_-5.0,
												-near_-2*corner_+15.0,
												result,
												font);
					}
			}
		font.setPixelSize(12);	
		result= "m/z";
		renderText (-corner_-20.0, 
								-corner_-20.0,
								-near_-3*corner_,
								result,
								font);
		font.setPixelSize(10);
		if(grid_mz_.size()>0)
		{
			for(UInt i = 0;i<grid_mz_[0].size();i++)
			{
				result = QString::number(grid_mz_[0][i]);
				renderText (-corner_-15.0, 
										-corner_-5.0,
										-near_-2*corner_-scaledMZ(grid_mz_[0][i]),
										result,
												font);
			}
		}
		if(zoom_<3.0 && grid_mz_.size()>=2)
			{
				
				for(UInt i = 0;i<grid_mz_[1].size();i++)
					{
						result = QString::number(grid_mz_[1][i]);
						renderText (-corner_-15.0, 
												-corner_-5.0,
												-near_-2*corner_-scaledMZ(grid_mz_[1][i]),
												result,
												font);
					}
			}
		if(zoom_<2.0 && grid_mz_.size()>=3)
			{
				for(UInt i = 0;i<grid_mz_[2].size();i++)
					{
						result = QString::number(grid_mz_[2][i]);
						renderText (-corner_-15.0, 
												-corner_-5.0,
												-near_-2*corner_-scaledMZ(grid_mz_[2][i]),
												result,
												font);
					}
			}
		//draw intensity legend if not in zoom mode (bird's eye view)
		if(!(canvas_3d_.action_mode_ == SpectrumCanvas::AM_ZOOM && zoom_mode_) )
			{
				switch (canvas_3d_.intensity_mode_)
					{	
					case SpectrumCanvas::IM_LOG:	
						font.setPixelSize(12);
						result= QString("intensity ");
						renderText (-corner_-20.0, 
												corner_+10.0,
												-near_-2*corner_+20.0,
												result,
												font);
						font.setPixelSize(10);
						if(zoom_<3 && grid_intensity_.size()>=1)
							{
								for(UInt i = 0;i<grid_intensity_log_[0].size();i++)
									{
										result = QString::number(grid_intensity_log_[0][i]);
										renderText (-corner_-result.length()-3.0, 
																-corner_+scaledIntensity(grid_intensity_log_[0][i],canvas_3d_.current_layer_),
																-near_-2*corner_,
																result,
																font);
									}
							}
						break;
					case SpectrumCanvas::IM_PERCENTAGE:
						font.setPixelSize(12);
						result =  QString("intensity %");
						renderText (-corner_-20.0, 
												corner_+10.0,
												-near_-2*corner_+20.0,
												result,
												font);
						font.setPixelSize(10);
						for(UInt i = 0;i<grid_intensity_[0].size();i++)
							{ 
								result = QString::number(grid_intensity_[0][i]);
								renderText (-corner_-result.length()-width_/200.0-5.0, 
														-corner_+(2.0*grid_intensity_[0][i]),
														-near_-2*corner_,
														result,
														font);
							}
						break;
					case SpectrumCanvas::IM_NONE:
					case SpectrumCanvas::IM_SNAP:
						int expo = 0;
						if(grid_intensity_.size()>=1)
						{
							expo = (int)ceil(log10(grid_intensity_[0][0]));
						}
						if(grid_intensity_.size()>=2)
						{
							if(expo>=ceil(log10(grid_intensity_[1][0])))
							{
								expo = (int)ceil(log10(grid_intensity_[1][0]));
							}
						}
						if(grid_intensity_.size()>=3)
						{
							if(expo>=ceil(log10(grid_intensity_[2][0])))
							{
								expo =(int) ceil(log10(grid_intensity_[2][0]));
							}	
						}
						//cout << "Expo: " << expo << endl;
						
						font.setPixelSize(12);
						result =   QString("intensity e+%1").arg((double)expo,0,'f',1);
						renderText (-corner_-20.0, 
												corner_+10.0,
												-near_-2*corner_+20.0,
												result,
												font);
						font.setPixelSize(10);
						if(zoom_<3.0 && grid_intensity_.size()>=2)
							{
								for(UInt i = 0;i<grid_intensity_[0].size();i++)
									{ 
										double intensity = (double)grid_intensity_[0][i]/pow(10,expo);
										result = QString("%1").arg(intensity,0,'f',1);
										renderText (-corner_-result.length()-width_/200.0-5.0, 
																-corner_+scaledIntensity(grid_intensity_[0][i],canvas_3d_.current_layer_),
																-near_-2*corner_,
																result,
																font);
									}
								for(UInt i = 0;i<grid_intensity_[1].size();i++)
									{
										double intensity = (double)grid_intensity_[1][i]/pow(10,expo);
										result = QString("%1").arg(intensity,0,'f',1);
										renderText (-corner_-result.length()-width_/200.0-5.0, 
																-corner_+scaledIntensity(grid_intensity_[1][i],canvas_3d_.current_layer_),
																-near_-2*corner_,
																result,
																font);
									}
							}
						if(width_>800 && heigth_>600&& zoom_<2.0 && grid_intensity_.size()>=3)
							{
								for(UInt i = 0;i<grid_intensity_[2].size();i++)
									{
										double intensity = (double)grid_intensity_[2][i]/pow(10,expo);
										result = QString("%1").arg(intensity,0,'f',1);
										renderText (-corner_-result.length()-width_/200.0-5.0, 
																-corner_+scaledIntensity(grid_intensity_[2][i],canvas_3d_.current_layer_),
																-near_-2*corner_,
																result,
																font);
									}
							}
					
					break;
				}
			}
		glEndList();
		return list;
	}
	
	GLuint Spectrum3DOpenGLCanvas::makeGround()
	{
		GLuint list = glGenLists(1);
		glNewList(list,GL_COMPILE);
		glBegin(GL_QUADS);
	 	QColor color(canvas_3d_.param_.getValue("background_color").toQString());
		qglColor(color);
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
	
	GLuint Spectrum3DOpenGLCanvas::makeCoordinates()
	{
		GLuint list = glGenLists(1);
		glNewList(list,GL_COMPILE);
		glLineWidth(3.0);
		glShadeModel(GL_FLAT);
		glBegin(GL_LINES);
		QColor color(canvas_3d_.param_.getValue("axis_color").toQString());
		
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
	
		for(UInt i =0;i<canvas_3d_.getLayerCount();i++)
		{	
			//cout << "Layer: " << i << endl;
			if(canvas_3d_.getLayer(i).visible)
			{	
				if((Int)canvas_3d_.getLayer(i).param.getValue("dot:shade_mode"))
				{
					glShadeModel(GL_SMOOTH); 
				}
				else
				{
					glShadeModel(GL_FLAT); 
				}	

				//cout << "  - VISIBLE" << endl;
				double min_int = canvas_3d_.getLayer(i).min_int;
				double max_int = canvas_3d_.getLayer(i).max_int;
				
				//cout << "  - RT: " << canvas_3d_.visible_area_.min_[0] << " " << canvas_3d_.visible_area_.max_[0] << endl;
				//cout << "  - MZ: " << canvas_3d_.visible_area_.min_[1] << " " << canvas_3d_.visible_area_.max_[1] << endl;
		
				for (Spectrum3DCanvas::ExperimentType::ConstAreaIterator it = canvas_3d_.getPeakData(i).areaBeginConst(canvas_3d_.visible_area_.min_[0],canvas_3d_.visible_area_.max_[0],canvas_3d_.visible_area_.min_[1],canvas_3d_.visible_area_.max_[1]); 
						 it != canvas_3d_.getPeakData(i).areaEndConst(); 
						 ++it)
				{
					if (it->getIntensity()>=min_int && it->getIntensity()<=max_int)
					{
						
						glBegin(GL_POINTS);
						double intensity = 0;
						switch (canvas_3d_.intensity_mode_)
						{
							case SpectrumCanvas::IM_NONE:
								qglColor(canvas_3d_.getLayer(i).gradient.precalculatedColorAt(it->getIntensity()));
								break;
							case SpectrumCanvas::IM_LOG:
								qglColor(canvas_3d_.getLayer(i).gradient.precalculatedColorAt(log10(it->getIntensity())));
								break;
							case SpectrumCanvas::IM_PERCENTAGE:	
								intensity = it->getIntensity() * 100.0 /canvas_3d_.getMaxIntensity(i);
								qglColor(canvas_3d_.getLayer(i).gradient.precalculatedColorAt(intensity ));
								break;
							case SpectrumCanvas::IM_SNAP:
								qglColor(canvas_3d_.getLayer(i).gradient.precalculatedColorAt(it->getIntensity()));
								break;
						}
						glVertex3d(-corner_+(GLfloat)scaledRT(it.getRT()), 
											 -corner_,
											 -near_-2*corner_-(GLfloat)scaledMZ(it->getMZ()));
						glEnd();		
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
	
		for(UInt i =0;i<canvas_3d_.getLayerCount();i++)
		{	
			//cout << "Layer: " << i << endl;
			if(canvas_3d_.getLayer(i).visible)
			{	

			if((Int)canvas_3d_.getLayer(i).param.getValue("dot:shade_mode"))
			{
				glShadeModel(GL_SMOOTH); 
			}
			else
			{
				glShadeModel(GL_FLAT); 
			}

			glLineWidth(canvas_3d_.getLayer(i).param.getValue("dot:line_width"));

				//cout << "  - VISIBLE" << endl;
				double min_int = canvas_3d_.getLayer(i).min_int;
				double max_int = canvas_3d_.getLayer(i).max_int;
				
				//cout << "  - RT: " << canvas_3d_.visible_area_.min_[0] << " " << canvas_3d_.visible_area_.max_[0] << endl;
				//cout << "  - MZ: " << canvas_3d_.visible_area_.min_[1] << " " << canvas_3d_.visible_area_.max_[1] << endl;
	
				for (Spectrum3DCanvas::ExperimentType::ConstAreaIterator it = canvas_3d_.getPeakData(i).areaBeginConst(canvas_3d_.visible_area_.min_[0],canvas_3d_.visible_area_.max_[0],canvas_3d_.visible_area_.min_[1],canvas_3d_.visible_area_.max_[1]); 
						 it != canvas_3d_.getPeakData(i).areaEndConst(); 
						 ++it)
				{
					if (it->getIntensity()>=min_int && it->getIntensity()<=max_int)
					{
						glBegin(GL_LINES);
						double intensity = 0;
						switch (canvas_3d_.intensity_mode_)
						{
						
						case SpectrumCanvas::IM_PERCENTAGE:	
							
							intensity = it->getIntensity() * 100.0 /canvas_3d_.getMaxIntensity(i);
							qglColor( canvas_3d_.getLayer(i).gradient.precalculatedColorAt(0));
							glVertex3d(-corner_+(GLfloat)scaledRT(it.getRT()), 
												 -corner_,
												 -near_-2*corner_-(GLfloat)scaledMZ(it->getMZ()));
							qglColor( canvas_3d_.getLayer(i).gradient.precalculatedColorAt(intensity ));
							glVertex3d(-corner_+(GLfloat)scaledRT(it.getRT()),
												 -corner_+(GLfloat)scaledIntensity(it->getIntensity(),i),
												 -near_-2*corner_-(GLfloat)scaledMZ(it->getMZ()));
								break;
						
						case SpectrumCanvas::IM_NONE:
						
							qglColor( canvas_3d_.getLayer(i).gradient.precalculatedColorAt(canvas_3d_.overall_data_range_.min_[2]));
							glVertex3d(-corner_+(GLfloat)scaledRT(it.getRT()), 
												 -corner_,
												 -near_-2*corner_-(GLfloat)scaledMZ(it->getMZ()));
							qglColor( canvas_3d_.getLayer(i).gradient.precalculatedColorAt(it->getIntensity()));
							glVertex3d(-corner_+(GLfloat)scaledRT(it.getRT()),
												 -corner_+(GLfloat)scaledIntensity(it->getIntensity(),i),
												 -near_-2*corner_-(GLfloat)scaledMZ(it->getMZ()));
							break;
						
						case SpectrumCanvas::IM_SNAP:
							
							qglColor(canvas_3d_.getLayer(i).gradient.precalculatedColorAt(int_scale_.min_[0]));
							glVertex3d(-corner_+(GLfloat)scaledRT(it.getRT()), 
												 -corner_,
												 -near_-2*corner_-(GLfloat)scaledMZ(it->getMZ()));
							qglColor(canvas_3d_.getLayer(i).gradient.precalculatedColorAt(it->getIntensity()));
							glVertex3d(-corner_+(GLfloat)scaledRT(it.getRT()),
												 -corner_+(GLfloat)scaledIntensity(it->getIntensity(),i),
												 -near_-2*corner_-(GLfloat)scaledMZ(it->getMZ()));
							
							break;
						
						case SpectrumCanvas::IM_LOG:
							qglColor(canvas_3d_.getLayer(i).gradient.precalculatedColorAt(log10(canvas_3d_.overall_data_range_.min_[2])));
							glVertex3d(-corner_+(GLfloat)scaledRT(it.getRT()), 
												 -corner_,
												 -near_-2*corner_-(GLfloat)scaledMZ(it->getMZ()));
							qglColor(canvas_3d_.getLayer(i).gradient.precalculatedColorAt(log10(it->getIntensity())));
							glVertex3d(-corner_+(GLfloat)scaledRT(it.getRT()),
												 -corner_+(GLfloat)scaledIntensity(log10(it->getIntensity()),i),
												 -near_-2*corner_-(GLfloat)scaledMZ(it->getMZ()));
							break;
							
						}
						glEnd();
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
		glEnable (GL_LINE_STIPPLE);
		glLineStipple (1, 0x0101);	
	 	glBegin(GL_LINES);
		glColor4ub(0, 0, 0, 80);	
		//rt
		if(grid_rt_.size()>=1)
			{
		for(UInt i = 0;i<grid_rt_[0].size();i++)
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
				for(UInt i = 0;i<grid_rt_[1].size();i++)
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
				for(UInt i = 0;i<grid_rt_[2].size();i++)
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
				for(UInt i = 0;i<grid_mz_[0].size();i++)
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
			for(UInt i = 0;i<grid_mz_[1].size();i++)
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
			for(UInt i = 0;i<grid_mz_[2].size();i++)
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
		glDisable (GL_LINE_STIPPLE);
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
		QColor color(canvas_3d_.param_.getValue("axis_color").toQString());
		qglColor(color);//x_achse
	
		//RT
		if(grid_rt_.size()>=1)
		{
			for(UInt i = 0;i<grid_rt_[0].size();i++)
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
			for(UInt i = 0;i<grid_rt_[1].size();i++)
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
			for(UInt i = 0;i<grid_rt_[2].size();i++)
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
			for(UInt i = 0;i<grid_mz_[0].size();i++)
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
				for(UInt i = 0;i<grid_mz_[1].size();i++)
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
			for(UInt i = 0;i<grid_mz_[2].size();i++)
			{
				glVertex3d(-corner_, 
									 -corner_,
									 -near_-2*corner_-scaledMZ(grid_mz_[2][i]));
				glVertex3d( -corner_,
										-corner_+2.0,
										-near_-2*corner_-scaledMZ(grid_mz_[2][i]));
			}
		}
	
		switch(canvas_3d_.intensity_mode_)
			{
			case SpectrumCanvas::IM_LOG:
				if(grid_intensity_.size()>=1)
					{
						for(UInt i = 0;i<grid_intensity_log_[0].size();i++)
							{	
								glVertex3d(-corner_, 
													 -corner_+scaledIntensity(grid_intensity_log_[0][i],canvas_3d_.current_layer_),
													 -near_-2*corner_);
								glVertex3d( -corner_+3.0, 
														-corner_+scaledIntensity(grid_intensity_log_[0][i],canvas_3d_.current_layer_),
														-near_-2*corner_-3.0);
							}
						for(UInt i = 0;i<grid_intensity_log_[1].size();i++)
							{
								glVertex3d(-corner_, 
													 -corner_+scaledIntensity(grid_intensity_log_[1][i],canvas_3d_.current_layer_),
													 -near_-2*corner_);
								glVertex3d( -corner_+2.0, 
														-corner_+scaledIntensity(grid_intensity_log_[1][i],canvas_3d_.current_layer_),
														-near_-2*corner_-2.0);
							}
					}
				break;
		
			case SpectrumCanvas::IM_PERCENTAGE:
				if(grid_intensity_.size()>=1)
					{
						for(UInt i = 0;i<grid_intensity_[0].size();i++)
							{	
								glVertex3d(-corner_, 
													 -corner_+(2.0 *grid_intensity_[0][i]),
													 -near_-2*corner_);
								glVertex3d( -corner_+4.0, 
														-corner_+(2.0 * grid_intensity_[0][i]),
														-near_-2*corner_-4.0);
							}
					}
				break;
			case SpectrumCanvas::IM_NONE:
				
			case SpectrumCanvas::IM_SNAP:
				if(grid_intensity_.size()>=1)
					{
						for(UInt i = 0;i<grid_intensity_[0].size();i++)
							{	
								glVertex3d(-corner_, 
													 -corner_+scaledIntensity(grid_intensity_[0][i],canvas_3d_.current_layer_),
													 -near_-2*corner_);
								glVertex3d( -corner_+4.0, 
														-corner_+scaledIntensity(grid_intensity_[0][i],canvas_3d_.current_layer_),
														-near_-2*corner_-4.0);
							}
					}
				if(grid_intensity_.size()>=2)
					{
						for(UInt i = 0;i<grid_intensity_[1].size();i++)
									{
										glVertex3d(-corner_, 
															 -corner_+scaledIntensity(grid_intensity_[1][i],canvas_3d_.current_layer_),
															 -near_-2*corner_);
										glVertex3d( -corner_+3.0, 
																-corner_+scaledIntensity(grid_intensity_[1][i],canvas_3d_.current_layer_),
																-near_-2*corner_-3.0);
									}
					}
				if(grid_intensity_.size()>=3)
					{
						for(UInt i = 0;i<grid_intensity_[2].size();i++)
							{ 
								glVertex3d(-corner_, 
													 -corner_+scaledIntensity(grid_intensity_[2][i],canvas_3d_.current_layer_),
													 -near_-2*corner_);
								glVertex3d( -corner_+2.0, 
														-corner_+scaledIntensity(grid_intensity_[2][i],canvas_3d_.current_layer_),
																	-near_-2*corner_-2.0);
							}
					}
				
		
				break;
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
		//	cout<<"rt"<<rt<<"  "<<"scaledinver"<<i_rt<<endl;
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
		double i_mz =(mz *canvas_3d_.visible_area_.max_[1] - mz *canvas_3d_.visible_area_.min_[1]);
		i_mz = i_mz/200;
		i_mz = i_mz + canvas_3d_.visible_area_.min_[1]; 
		return i_mz;
	}
	
	double Spectrum3DOpenGLCanvas::scaledIntensity(double intensity,int layer_index)
	{
		double scaledintensity= 0;
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
		 		scaledintensity =  intensity * 100.0 /canvas_3d_.getMaxIntensity(layer_index);
				scaledintensity = scaledintensity * 2.0 * corner_/100.0;	
				break;	
			case SpectrumCanvas::IM_LOG:
				double log_min = 0.0;
				if(log10(canvas_3d_.overall_data_range_.min_[2])<0)
				{
					log_min = 0.0;
				}
				else
				{
					log_min = log10(canvas_3d_.overall_data_range_.min_[2]);
				}
				scaledintensity = intensity -log_min;	
				scaledintensity =(scaledintensity * 2.0 * corner_)/(log10(canvas_3d_.overall_data_range_.max_[2])-log_min);
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
	
	void Spectrum3DOpenGLCanvas::setZoomFactor(double zoom, bool repaint)
	{
	 	zoom_ = zoom;
		if(repaint)
		{
			resizeGL((int)width_,(int) heigth_);
			glDraw (); 
		}
	}
	
	///////////////wheel- and MouseEvents//////////////////
	void Spectrum3DOpenGLCanvas::wheelEvent (QWheelEvent* e)
	{
		if(!zoom_mode_) //translate
		{
			double zoom = zoom_ + double(e->delta()/480.0);
			if(zoom>0.0)
			{	
				setZoomFactor( zoom,true);
			}
			else
			{
				zoom = 0.25;
				setZoomFactor( zoom,true);
			}
		}
	}

	void Spectrum3DOpenGLCanvas::keyPressEvent(QKeyEvent* e)
	{
		if (canvas_3d_.action_mode_== SpectrumCanvas::AM_ZOOM && e->key()==Qt::Key_Control)
		{
			zoom_mode_ = true;
			setAngels(220,220,0);
			setZoomFactor(1.5,false);
			canvas_3d_.update_buffer_ = true;
			canvas_3d_.update_(__PRETTY_FUNCTION__);
		  e->accept();
		  return;
		}
		e->ignore();
	}
	
	void Spectrum3DOpenGLCanvas::keyReleaseEvent(QKeyEvent* e)
	{
		if (canvas_3d_.action_mode_== SpectrumCanvas::AM_ZOOM && e->key()==Qt::Key_Control)
		{
			zoom_mode_ = false;
			setAngels(1440,0,0);
			resetTranslation();
			setZoomFactor(1.25,false);
			canvas_3d_.update_buffer_ = true;
			canvas_3d_.update_(__PRETTY_FUNCTION__);
		  e->accept();
		  return;
		}
		e->ignore();
	}

	void Spectrum3DOpenGLCanvas::focusOutEvent(QFocusEvent* e)
	{
		if (canvas_3d_.action_mode_== SpectrumCanvas::AM_ZOOM && zoom_mode_==true)
		{
			zoom_mode_ = false;
			setAngels(1440,0,0);
			resetTranslation();
			setZoomFactor(1.25,false);
			canvas_3d_.update_buffer_ = true;
			canvas_3d_.update_(__PRETTY_FUNCTION__);
		  e->accept();
		  return;
		}
		e->ignore();
	}

	void Spectrum3DOpenGLCanvas::mousePressEvent (QMouseEvent* e)
	{
		mouse_move_begin_ = e->pos();
		mouse_move_end_ = e->pos();	

	  switch(canvas_3d_.action_mode_)
	  {
	  	case SpectrumCanvas::AM_ZOOM:
				if(zoom_mode_ && e->button()==Qt::LeftButton)
				{			
					canvas_3d_.rubber_band_.setGeometry(e->pos().x(),e->pos().y(),0,0);
					canvas_3d_.rubber_band_.show();
					canvas_3d_.update_(__PRETTY_FUNCTION__);
				}
	      break;
		  default:
		    break;
	  }
	}

	void Spectrum3DOpenGLCanvas::mouseMoveEvent(QMouseEvent* e)
	{
	  switch(canvas_3d_.action_mode_)
	  {
	  	case SpectrumCanvas::AM_ZOOM:
				if(zoom_mode_ && e->buttons() & Qt::LeftButton)
				{
					canvas_3d_.rubber_band_.setGeometry(QRect(mouse_move_begin_, e->pos()));
					canvas_3d_.update_(__PRETTY_FUNCTION__);
				}
				else if(!zoom_mode_ && e->buttons() & Qt::LeftButton)
				{
					int d_x = e->x() - mouse_move_end_.x();
					int d_y = e->y() - mouse_move_end_.y();
					setRotationX(xrot_ + 8 * d_y);
					setRotationY(yrot_ + 8 * d_x);
					mouse_move_end_ = e->pos();
					canvas_3d_.update_(__PRETTY_FUNCTION__);
				}
				else if(!zoom_mode_ && e->buttons() & Qt::RightButton)
				{
					mouse_move_end_ = e->pos();
					trans_x_= mouse_move_end_.x()-mouse_move_begin_.x();
					trans_y_ = (heigth_-mouse_move_end_.y())-(heigth_ -mouse_move_begin_.y());
					canvas_3d_.update_(__PRETTY_FUNCTION__);
				}
	      break;
		  default:
		    break;
	  }
	}
	
	void Spectrum3DOpenGLCanvas::mouseReleaseEvent (QMouseEvent* e)
	{
		if(canvas_3d_.action_mode_ == SpectrumCanvas::AM_ZOOM && zoom_mode_ && e->button()==Qt::LeftButton)
		{				
			QRect rect = canvas_3d_.rubber_band_.geometry();
	 		x_1_ = ((rect.topLeft().x()- width_/2) * corner_ *1.25* 2) / width_;
 			y_1_ = -300 + (((rect.topLeft().y()-heigth_/2) * corner_*1.25* 2) / heigth_);
 			x_2_ = ((rect.bottomRight().x()- width_/2) * corner_ *1.25* 2) / width_;
 			y_2_ = -300 + (((rect.bottomRight().y()-heigth_/2) * corner_*1.25* 2) / heigth_);
			dataToZoomArray(x_1_, y_1_, x_2_, y_2_);
			canvas_3d_.rubber_band_.hide();
			canvas_3d_.update_buffer_ = true;
			canvas_3d_.update_(__PRETTY_FUNCTION__);
		}
	}
	
	void Spectrum3DOpenGLCanvas::dataToZoomArray(double x_1, double y_1, double x_2, double y_2)
	{
		double scale_x1 = scaledInversRT(x_1+100.0);
		double scale_x2 = scaledInversRT(x_2+100.0);
		double scale_y1 = scaledInversMZ(-200-y_1);
		double scale_y2 = scaledInversMZ(-200-y_2);
		DRange<2> new_area_;
		if(scale_x1<=scale_x2)
		{
			new_area_.min_[0]= scale_x1;
			new_area_.max_[0]= scale_x2;
		} 
		else
		{
			new_area_.min_[0]= scale_x2;
			new_area_.max_[0]= scale_x1;
		}
		if(scale_y1<=scale_y2)
		{
			new_area_.min_[1]= scale_y1;
			new_area_.max_[1]= scale_y2;
		} 
		else
		{
		 new_area_.min_[1]= scale_y2;
		 new_area_.max_[1]= scale_y1;
		} 
		canvas_3d_.changeVisibleArea_(new_area_, true);
	}
	
	void Spectrum3DOpenGLCanvas::updateIntensityScale()
	{
		int_scale_.min_[0]= canvas_3d_.overall_data_range_.max_[2];
		int_scale_.max_[0]= canvas_3d_.overall_data_range_.min_[2];
		
		for(UInt i =0;i<canvas_3d_.getLayerCount();i++)
		{
			for (SpectrumCanvas::ExperimentType::ConstIterator spec_it = canvas_3d_.getPeakData(i).RTBegin(canvas_3d_.visible_area_.min_[0]); 
					 spec_it != canvas_3d_.getPeakData(i).RTEnd(canvas_3d_.visible_area_.max_[0]); 
					 ++spec_it)
			{
				for (SpectrumCanvas::ExperimentType::SpectrumType::ConstIterator it = spec_it->MZBegin(canvas_3d_.visible_area_.min_[1]); it!=spec_it->MZEnd(canvas_3d_.visible_area_.max_[1]); ++it)
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
	
	void Spectrum3DOpenGLCanvas::recalculateDotGradient_(UInt layer)
	{
		canvas_3d_.getLayer_(layer).gradient.fromString(canvas_3d_.getLayer_(layer).param.getValue("dot:gradient"));
		switch(canvas_3d_.intensity_mode_)
		{
			case SpectrumCanvas::IM_SNAP:
				canvas_3d_.getLayer_(layer).gradient.activatePrecalculationMode(int_scale_.min_[0], int_scale_.max_[0], UInt(canvas_3d_.param_.getValue("dot:interpolation_steps")));
				break;
			case SpectrumCanvas::IM_NONE:
				canvas_3d_.getLayer_(layer).gradient.activatePrecalculationMode(canvas_3d_.overall_data_range_.min_[2], canvas_3d_.overall_data_range_.max_[2], UInt(canvas_3d_.param_.getValue("dot:interpolation_steps")));
				break;
			case SpectrumCanvas::IM_LOG:
				double log_min;
				if(log10(canvas_3d_.overall_data_range_.min_[2])<0)
				{
					log_min = 0.0;
				}
				else
				{
					log_min = log10(canvas_3d_.overall_data_range_.min_[2]);
				}
				canvas_3d_.getLayer_(layer).gradient.activatePrecalculationMode(log_min, log10(canvas_3d_.overall_data_range_.max_[2]), UInt(canvas_3d_.param_.getValue("dot:interpolation_steps")));
				
				break;
			case SpectrumCanvas::IM_PERCENTAGE:
				canvas_3d_.getLayer_(layer).gradient.activatePrecalculationMode(0.0, 100.0, UInt(canvas_3d_.param_.getValue("dot:interpolation_steps")));
				break;
		}
	}

}//end of namespace
