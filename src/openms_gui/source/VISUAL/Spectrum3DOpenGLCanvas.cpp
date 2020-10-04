// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/Spectrum3DOpenGLCanvas.h>

#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/VISUAL/AxisTickCalculator.h>

#include <QMouseEvent>
#include <QKeyEvent>

using std::cout;
using std::endl;
using std::max;

namespace OpenMS
{

  Spectrum3DOpenGLCanvas::Spectrum3DOpenGLCanvas(QWidget * parent, Spectrum3DCanvas & canvas_3d) :
    QOpenGLWidget(parent),
    canvas_3d_(canvas_3d)
  {
    canvas_3d.rubber_band_.setParent(this);

    x_label_ = (String(Peak2D::shortDimensionName(Peak2D::MZ)) + " [" + String(Peak2D::shortDimensionUnit(Peak2D::MZ)) + "]").toQString();
    y_label_ = (String(Peak2D::shortDimensionName(Peak2D::RT)) + " [" + String(Peak2D::shortDimensionUnit(Peak2D::RT)) + "]").toQString();

    //Set focus policy and mouse tracking in order to get keyboard events
    setMouseTracking(true);
    setFocusPolicy(Qt::StrongFocus);

    corner_ = 100.0;
    near_ = 0.0;
    far_ = 600.0;
    zoom_ = 1.5;
    xrot_ = 220;
    yrot_ = 220;
    zrot_ = 0;
    trans_x_ = 0.0;
    trans_y_ = 0.0;
  }

  Spectrum3DOpenGLCanvas::~Spectrum3DOpenGLCanvas()
  {
  }

  void Spectrum3DOpenGLCanvas::calculateGridLines_()
  {
    switch (canvas_3d_.intensity_mode_)
    {
    case SpectrumCanvas::IM_SNAP:
      updateIntensityScale();
      AxisTickCalculator::calcGridLines(0.0, int_scale_.max_[0], grid_intensity_);
      break;

    case SpectrumCanvas::IM_NONE:
      AxisTickCalculator::calcGridLines(0.0, canvas_3d_.overall_data_range_.max_[2], grid_intensity_);
      break;

    case SpectrumCanvas::IM_PERCENTAGE:
      AxisTickCalculator::calcGridLines(0.0, 100.0, grid_intensity_);
      break;

    case SpectrumCanvas::IM_LOG:
      AxisTickCalculator::calcLogGridLines(0.0, log10(1 + max(0.0, canvas_3d_.overall_data_range_.max_[2])), grid_intensity_);
      break;
    }

    AxisTickCalculator::calcGridLines(canvas_3d_.visible_area_.min_[1], canvas_3d_.visible_area_.max_[1], grid_rt_);
    AxisTickCalculator::calcGridLines(canvas_3d_.visible_area_.min_[0], canvas_3d_.visible_area_.max_[0], grid_mz_);
  }

  void Spectrum3DOpenGLCanvas::transformPoint_(GLdouble out[4], const GLdouble m[16], const GLdouble in[4])
  {
  #define M(row,col)  m[col*4+row]
    out[0] = M(0, 0) * in[0] + M(0, 1) * in[1] + M(0, 2) * in[2] + M(0, 3) * in[3];
    out[1] = M(1, 0) * in[0] + M(1, 1) * in[1] + M(1, 2) * in[2] + M(1, 3) * in[3];
    out[2] = M(2, 0) * in[0] + M(2, 1) * in[1] + M(2, 2) * in[2] + M(2, 3) * in[3];
    out[3] = M(3, 0) * in[0] + M(3, 1) * in[1] + M(3, 2) * in[2] + M(3, 3) * in[3];
  #undef M
  }

  GLint Spectrum3DOpenGLCanvas::project_(GLdouble objx, GLdouble objy, GLdouble objz, GLdouble * winx, GLdouble * winy)
  {
    int height= this->height();
    int width = this->width();

    GLdouble in[4], out[4];

    in[0] = objx;
    in[1] = objy;
    in[2] = objz;
    in[3] = 1.0;


    GLdouble model[16];
    GLdouble proj[16];

    glGetDoublev(GL_MODELVIEW_MATRIX, model);
    glGetDoublev(GL_PROJECTION_MATRIX, proj);

    transformPoint_(out, model, in);
    transformPoint_(in, proj, out);

    // transform homogeneous coordinates into normalized device coordinates
    if (in[3] == 0.0) { return GL_FALSE; }
    in[0] /= in[3];
    in[1] /= in[3];
    in[2] /= in[3];

    // viewport transformation (0,0 is in corner of screen not in middle of screen)
    *winx = 0 + (1 + in[0]) * width / 2;
    *winy = 0 + (1 + in[1]) * height / 2;

    return GL_TRUE;
  }  

  void Spectrum3DOpenGLCanvas::renderText_(double x, double y, double z, const QString & text) 
  {
    // project to screen coordinates
    GLdouble textPosX = 0, textPosY = 0;
    project_(x, y, z, &textPosX, &textPosY);
    const int height = this->height();
    textPosY = height - textPosY; // y is inverted

    // cout << x << " " << y << " " << z << " : " << textPosX << " " << textPosY << " " << textPosZ << endl;

    // render text using the current font and color
    painter_->drawText(textPosX, textPosY, text);
  }

  void Spectrum3DOpenGLCanvas::qglColor_(QColor color) 
  {
    glColor4f(color.redF(), color.greenF(), color.blueF(), color.alphaF());
  }

  void Spectrum3DOpenGLCanvas::qglClearColor_(QColor clearColor) 
  {
    glClearColor(clearColor.redF(), clearColor.greenF(), clearColor.blueF(), clearColor.alphaF());
  }

  void Spectrum3DOpenGLCanvas::resizeGL(int w, int h)
  {
    width_ = (float)w;
    height_ = (float)h;
    glViewport(0, 0, (GLsizei) w, (GLsizei) h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-corner_ * zoom_, corner_ * zoom_, -corner_ * zoom_, corner_ * zoom_, near_, far_);
    glMatrixMode(GL_MODELVIEW);
  }

  void Spectrum3DOpenGLCanvas::initializeGL()
  {
    initializeOpenGLFunctions();

    // The following line triggered a bug where the whole screen would turn
    // black during scrolling (specifically it seems that multiple calls to
    // this function causes the issue):
    // glEnable(GL_DEPTH_TEST);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    QColor color(canvas_3d_.param_.getValue("background_color").toQString());
    qglClearColor_(color);
    calculateGridLines_();

    //abort if no layers are displayed
    if (canvas_3d_.getLayerCount() == 0) { return; }

    if (canvas_3d_.action_mode_ == SpectrumCanvas::AM_ZOOM)
    {
      if (!canvas_3d_.rubber_band_.isVisible())
      {
        axes_ = makeAxes_();
        if (canvas_3d_.show_grid_)
        {
          gridlines_ = makeGridLines_();
        }
        xrot_ = 90 * 16;
        yrot_ = 0;
        zrot_ = 0;
        zoom_ = 1.25;

        if (stickdata_ != 0)
        {
          glDeleteLists(stickdata_, 1);
#ifdef DEBUG_TOPPVIEW
          cout << "BEGIN " << OPENMS_PRETTY_FUNCTION << endl;
          std::cout << "Deleting sticklist near makeDataAsTopView" << std::endl;
#endif
        }

        stickdata_ = makeDataAsTopView_();
        axes_ticks_ = makeAxesTicks_();
        //drawAxesLegend_();
      }
    }
    else if (canvas_3d_.action_mode_ == SpectrumCanvas::AM_TRANSLATE)
    {
      if (canvas_3d_.show_grid_) { gridlines_ = makeGridLines_(); }
      axes_ = makeAxes_();
      ground_ = makeGround_();
      x_1_ = 0.0;
      y_1_ = 0.0;
      x_2_ = 0.0;
      y_2_ = 0.0;

      if (stickdata_ != 0)
      {
        glDeleteLists(stickdata_, 1);
#ifdef DEBUG_TOPPVIEW
        cout << "BEGIN " << OPENMS_PRETTY_FUNCTION << endl;
        std::cout << "Deleting sticklist near makeDataAsTopView" << std::endl;
#endif
      }

      stickdata_ =  makeDataAsStick_();
      axes_ticks_ = makeAxesTicks_();
      //drawAxesLegend_();
    }
  }

  void Spectrum3DOpenGLCanvas::resetTranslation()
  {
    trans_x_ = 0.0;
    trans_y_ = 0.0;
  }

  void Spectrum3DOpenGLCanvas::storeRotationAndZoom()
  {
    xrot_tmp_ = xrot_;
    yrot_tmp_ = yrot_;
    zrot_tmp_ = zrot_;
    zoom_tmp_ = zoom_;
  }

  void Spectrum3DOpenGLCanvas::restoreRotationAndZoom()
  {
    xrot_ = xrot_tmp_;
    yrot_ = yrot_tmp_;
    zrot_ = zrot_tmp_;
    zoom_ = zoom_tmp_;
  }

  void Spectrum3DOpenGLCanvas::paintGL()
  {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    glTranslated(0.0, 0.0, -3.0 * corner_);
    glRotated(xrot_ / 16.0, 1.0, 0.0, 0.0);
    glRotated(yrot_ / 16.0, 0.0, 1.0, 0.0);
    glRotated(zrot_ / 16.0, 0.0, 0.0, 1.0);
    glTranslated(trans_x_, trans_y_, 3.0 * corner_);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    if (canvas_3d_.getLayerCount() != 0)
    {
      glCallList(ground_);

      if (canvas_3d_.show_grid_) { glCallList(gridlines_); }

      glCallList(axes_);
      glCallList(axes_ticks_);

      if (canvas_3d_.action_mode_ == SpectrumCanvas::AM_ZOOM 
       || canvas_3d_.action_mode_ == SpectrumCanvas::AM_TRANSLATE)
      {
        glCallList(stickdata_);
      }

      // draw axes legend
      if (this->paintEngine()) // check if the paint device is properly initialized to suppress Qt warning
      {
        painter_ = new QPainter(this);
        if (painter_->isActive())
        {
          drawAxesLegend_(); 
          painter_->end();
        }
        delete(painter_);
      }
    }
  }

  void Spectrum3DOpenGLCanvas::drawAxesLegend_()
  {
    QFont font("Typewriter");
    font.setPixelSize(10);

    QString text;
    qglColor_(Qt::black);

    // Draw x and y axis legend
    if (canvas_3d_.legend_shown_)
    {
      font.setPixelSize(12);
      renderText_(0.0, -corner_ - 20.0, -near_ - 2 * corner_ + 20.0, x_label_);
      renderText_(-corner_ - 20.0, -corner_ - 20.0, -near_ - 3 * corner_, y_label_);
      font.setPixelSize(10);
    }

    // RT tick labels
    {
      if (grid_rt_.size() > 0)
      {
        for (Size i = 0; i < grid_rt_[0].size(); i++)
        {
          text = QString::number(grid_rt_[0][i]);
          renderText_(-corner_ - 15.0, -corner_ - 5.0, -near_ - 2 * corner_ - scaledRT_(grid_rt_[0][i]), text);
        }
      }
      if (zoom_ < 3.0 && grid_rt_.size() >= 2)
      {
        for (Size i = 0; i < grid_rt_[1].size(); i++)
        {
          text = QString::number(grid_rt_[1][i]);
          renderText_(-corner_ - 15.0, -corner_ - 5.0, -near_ - 2 * corner_ - scaledRT_(grid_rt_[1][i]), text);
        }
      }
      if (zoom_ < 2.0 && grid_rt_.size() >= 3)
      {
        for (Size i = 0; i < grid_rt_[2].size(); i++)
        {
          text = QString::number(grid_rt_[2][i]);
          renderText_(-corner_ - 15.0, -corner_ - 5.0, -near_ - 2 * corner_ - scaledRT_(grid_rt_[2][i]), text);
        }
      }
    }

    // m/z tick labels
    {
      if (grid_mz_.size() > 0)
      {
        for (Size i = 0; i < grid_mz_[0].size(); i++)
        {
          text = QString::number(grid_mz_[0][i]);
          renderText_(-corner_ - text.length() + scaledMZ_(grid_mz_[0][i]), -corner_ - 5.0, -near_ - 2 * corner_ + 15.0, text);
        }
      }
      if (zoom_ < 3.0 && grid_mz_.size() >= 2)
      {
        for (Size i = 0; i < grid_mz_[1].size(); i++)
        {
          text = QString::number(grid_mz_[1][i]);
          renderText_(-corner_ - text.length() + scaledMZ_(grid_mz_[1][i]), -corner_ - 5.0, -near_ - 2 * corner_ + 15.0, text);
        }
      }
      if (zoom_ < 2.0 && grid_mz_.size() >= 3)
      {
        for (Size i = 0; i < grid_mz_[2].size(); i++)
        {
          text = QString::number(grid_mz_[2][i]);
          renderText_(-corner_ - text.length() + scaledMZ_(grid_mz_[2][i]), -corner_ - 5.0, -near_ - 2 * corner_ + 15.0, text);
        }
      }
    }

    // draw intensity legend if not in zoom mode
    if (canvas_3d_.action_mode_ != SpectrumCanvas::AM_ZOOM)
    {
      switch (canvas_3d_.intensity_mode_)
      {
      case SpectrumCanvas::IM_LOG:

        if (canvas_3d_.legend_shown_)
        {
          font.setPixelSize(12);
          text = QString("intensity log");
          renderText_(-corner_ - 20.0, corner_ + 10.0, -near_ - 2 * corner_ + 20.0, text);
          font.setPixelSize(10);
        }

        if (zoom_ < 3.0 && grid_intensity_.size() >= 2)
        {
          for (Size i = 0; i < grid_intensity_[0].size(); i++)
          {
            double intensity = (double)grid_intensity_[0][i];
            text = QString("%1").arg(intensity, 0, 'f', 0);
            renderText_(-corner_ - text.length() - width_ / 200.0 - 5.0, -corner_ + scaledIntensity_(pow(10.0, grid_intensity_[0][i]) - 1, canvas_3d_.layers_.getCurrentLayerIndex()), -near_ - 2 * corner_, text);
          }
        }
        break;

      case SpectrumCanvas::IM_PERCENTAGE:
        if (canvas_3d_.legend_shown_)
        {
          font.setPixelSize(12);
          renderText_(-corner_ - 20.0, corner_ + 10.0, -near_ - 2 * corner_ + 20.0, "intensity %");
          font.setPixelSize(10);
        }

        for (Size i = 0; i < grid_intensity_[0].size(); i++)
        {
          text = QString::number(grid_intensity_[0][i]);
          renderText_(-corner_ - text.length() - width_ / 200.0 - 5.0, -corner_ + (2.0 * grid_intensity_[0][i]), -near_ - 2 * corner_, text);
        }
        break;

      case SpectrumCanvas::IM_NONE:
      case SpectrumCanvas::IM_SNAP:
        int expo = 0;
        if (grid_intensity_.size() >= 1)
        {
          expo = (int)ceil(log10(grid_intensity_[0][0]));
        }
        if (grid_intensity_.size() >= 2)
        {
          if (expo >= ceil(log10(grid_intensity_[1][0])))
          {
            expo = (int)ceil(log10(grid_intensity_[1][0]));
          }
        }
        if (grid_intensity_.size() >= 3)
        {
          if (expo >= ceil(log10(grid_intensity_[2][0])))
          {
            expo = (int) ceil(log10(grid_intensity_[2][0]));
          }
        }

        if (canvas_3d_.legend_shown_)
        {
          font.setPixelSize(12);
          text = QString("intensity e+%1").arg((double)expo, 0, 'f', 1);
          renderText_(-corner_ - 20.0, corner_ + 10.0, -near_ - 2 * corner_ + 20.0, text);
          font.setPixelSize(10);
        }


        if (zoom_ < 3.0 && grid_intensity_.size() >= 2)
        {
          for (Size i = 0; i < grid_intensity_[0].size(); i++)
          {
            double intensity = (double)grid_intensity_[0][i] / pow(10.0, expo);
            text = QString("%1").arg(intensity, 0, 'f', 1);
            renderText_(-corner_ - text.length() - width_ / 200.0 - 5.0, -corner_ + scaledIntensity_(grid_intensity_[0][i], canvas_3d_.layers_.getCurrentLayerIndex()), -near_ - 2 * corner_, text);
          }
          for (Size i = 0; i < grid_intensity_[1].size(); i++)
          {
            double intensity = (double)grid_intensity_[1][i] / pow(10.0, expo);
            text = QString("%1").arg(intensity, 0, 'f', 1);
            renderText_(-corner_ - text.length() - width_ / 200.0 - 5.0, -corner_ + scaledIntensity_(grid_intensity_[1][i], canvas_3d_.layers_.getCurrentLayerIndex()), -near_ - 2 * corner_, text);
          }
        }
        if (width_ > 800 && height_ > 600 && zoom_ < 2.0 && grid_intensity_.size() >= 3)
        {
          for (Size i = 0; i < grid_intensity_[2].size(); i++)
          {
            double intensity = (double)grid_intensity_[2][i] / pow(10.0, expo);
            text = QString("%1").arg(intensity, 0, 'f', 1);
            renderText_(-corner_ - text.length() - width_ / 200.0 - 5.0, -corner_ + scaledIntensity_(grid_intensity_[2][i], canvas_3d_.layers_.getCurrentLayerIndex()), -near_ - 2 * corner_, text);
          }
        }
        break;

      }
    }
  }

  GLuint Spectrum3DOpenGLCanvas::makeGround_()
  {
    GLuint list = glGenLists(1);
    glNewList(list, GL_COMPILE);
    glBegin(GL_QUADS);
    QColor color(canvas_3d_.param_.getValue("background_color").toQString());
    qglColor_(color);
    glVertex3d(-corner_, -corner_ - 2.0, -near_ - 2 * corner_);
    glVertex3d(-corner_, -corner_ - 2.0, -far_ + 2 * corner_);
    glVertex3d(corner_, -corner_ - 2.0, -far_ + 2 * corner_);
    glVertex3d(corner_, -corner_ - 2.0, -near_ - 2 * corner_);
    glEnd();
    glEndList();
    return list;
  }

  GLuint Spectrum3DOpenGLCanvas::makeAxes_()
  {
    GLuint list = glGenLists(1);
    glNewList(list, GL_COMPILE);
    glLineWidth(3.0);
    glShadeModel(GL_FLAT);
    glBegin(GL_LINES);
    qglColor_(Qt::black);
    // x
    glVertex3d(-corner_, -corner_, -near_ - 2 * corner_);
    glVertex3d(corner_, -corner_, -near_ - 2 * corner_);
    // z
    glVertex3d(-corner_, -corner_, -near_ - 2 * corner_);
    glVertex3d(-corner_, -corner_, -far_ + 2 * corner_);
    // y
    glVertex3d(-corner_, -corner_, -near_ - 2 * corner_);
    glVertex3d(-corner_, corner_, -near_ - 2 * corner_);
    glEnd();
    glEndList();
    return list;
  }

  GLuint Spectrum3DOpenGLCanvas::makeDataAsTopView_()
  {
    GLuint list = glGenLists(1);
    glNewList(list, GL_COMPILE);
    glPointSize(3.0);

    for (Size i = 0; i < canvas_3d_.getLayerCount(); ++i)
    {
      const LayerData & layer = canvas_3d_.getLayer(i);
      if (layer.visible)
      {
        if ((Int)layer.param.getValue("dot:shade_mode"))
        {
          glShadeModel(GL_SMOOTH);
        }
        else
        {
          glShadeModel(GL_FLAT);
        }

        auto begin_it = layer.getPeakData()->areaBeginConst(canvas_3d_.visible_area_.min_[1], canvas_3d_.visible_area_.max_[1], canvas_3d_.visible_area_.min_[0], canvas_3d_.visible_area_.max_[0]);
        auto end_it = layer.getPeakData()->areaEndConst();

        // count peaks in area
        int count = std::distance(begin_it, end_it);

        int max_displayed_peaks = 10000;
        int step = 1;
        if (count > max_displayed_peaks)
        {
          step = 1 + count / max_displayed_peaks;
        }

        for (auto it = begin_it; it != end_it; ++it)
        {
          if (step > 1)
          {
            for (int i = 0; i < step - 1; ++i)
            {
              ++it;
            }
          }
          if (it == end_it)
          {
            glEndList();
            return list;
          }

          PeakIndex pi = it.getPeakIndex();
          if (layer.filters.passes((*layer.getPeakData())[pi.spectrum], pi.peak))
          {
            glBegin(GL_POINTS);
            double intensity = 0;
            switch (canvas_3d_.intensity_mode_)
            {
            case SpectrumCanvas::IM_NONE:
              qglColor_(layer.gradient.precalculatedColorAt(it->getIntensity()));
              break;

            case SpectrumCanvas::IM_PERCENTAGE:
              intensity = it->getIntensity() * 100.0 / canvas_3d_.getMaxIntensity(i);
              qglColor_(layer.gradient.precalculatedColorAt(intensity));
              break;

            case SpectrumCanvas::IM_SNAP:
              qglColor_(layer.gradient.precalculatedColorAt(it->getIntensity()));
              break;

            case SpectrumCanvas::IM_LOG:
              qglColor_(layer.gradient.precalculatedColorAt(log10(1 + max(0.0, (double)(it->getIntensity())))));
              break;

            }
            glVertex3f(-corner_ + (GLfloat)scaledMZ_(it->getMZ()),
                       -corner_,
                       -near_ - 2 * corner_ - (GLfloat)scaledRT_(it.getRT()));
            glEnd();
          }
        }
      }
    }
    glEndList();
    return list;
  }

  GLuint Spectrum3DOpenGLCanvas::makeDataAsStick_()
  {
    GLuint list = glGenLists(1);
    glNewList(list, GL_COMPILE);

    for (Size i = 0; i < canvas_3d_.getLayerCount(); i++)
    {
      LayerData& layer = canvas_3d_.getLayer(i);
      if (layer.visible)
      {
        recalculateDotGradient_(layer);

        if ((Int)layer.param.getValue("dot:shade_mode"))
        {
          glShadeModel(GL_SMOOTH);
        }
        else
        {
          glShadeModel(GL_FLAT);
        }

        glLineWidth(layer.param.getValue("dot:line_width"));

        auto begin_it = layer.getPeakData()->areaBeginConst(canvas_3d_.visible_area_.min_[1], canvas_3d_.visible_area_.max_[1], canvas_3d_.visible_area_.min_[0], canvas_3d_.visible_area_.max_[0]);
        auto end_it = layer.getPeakData()->areaEndConst();
        // count peaks in area
        int count = std::distance(begin_it, end_it);

        int max_displayed_peaks = 100000;
        int step = 1;
        if (count > max_displayed_peaks)
        {
          step = 1 + count / max_displayed_peaks;
        }

        for (auto it = begin_it; it != end_it; ++it)
        {
          if (step > 1)
          {
            for (int i = 0; i < step - 1; ++i)
            {
              ++it;
            }
          }
          if (it == end_it)
          {
            glEndList();
            return list;
          }

          PeakIndex pi = it.getPeakIndex();
          if (layer.filters.passes((*layer.getPeakData())[pi.spectrum], pi.peak))
          {
            glBegin(GL_LINES);
            double intensity = 0;
            switch (canvas_3d_.intensity_mode_)
            {

            case SpectrumCanvas::IM_PERCENTAGE:

              intensity = it->getIntensity() * 100.0 / canvas_3d_.getMaxIntensity(i);
              qglColor_(layer.gradient.precalculatedColorAt(0.0));
              glVertex3d(-corner_ + (GLfloat)scaledMZ_(it->getMZ()),
                         -corner_,
                         -near_ - 2 * corner_ - (GLfloat)scaledRT_(it.getRT()));
              qglColor_(layer.gradient.precalculatedColorAt(intensity));
              glVertex3d(-corner_ + (GLfloat)scaledMZ_(it->getMZ()),
                         -corner_ + (GLfloat)scaledIntensity_(it->getIntensity(), i),
                         -near_ - 2 * corner_ - (GLfloat)scaledRT_(it.getRT()));
              break;

            case SpectrumCanvas::IM_NONE:

              qglColor_(layer.gradient.precalculatedColorAt(0.0));
              glVertex3d(-corner_ + (GLfloat)scaledMZ_(it->getMZ()),
                         -corner_,
                         -near_ - 2 * corner_ - (GLfloat)scaledRT_(it.getRT()));
              qglColor_(layer.gradient.precalculatedColorAt(it->getIntensity()));
              glVertex3d(-corner_ + (GLfloat)scaledMZ_(it->getMZ()),
                         -corner_ + (GLfloat)scaledIntensity_(it->getIntensity(), i),
                         -near_ - 2 * corner_ - (GLfloat)scaledRT_(it.getRT()));
              break;

            case SpectrumCanvas::IM_SNAP:

              qglColor_(layer.gradient.precalculatedColorAt(0.0));
              glVertex3d(-corner_ + (GLfloat)scaledMZ_(it->getMZ()),
                         -corner_,
                         -near_ - 2 * corner_ - (GLfloat)scaledRT_(it.getRT()));
              qglColor_(layer.gradient.precalculatedColorAt(it->getIntensity()));
              glVertex3d(-corner_ + (GLfloat)scaledMZ_(it->getMZ()),
                         -corner_ + (GLfloat)scaledIntensity_(it->getIntensity(), i),
                         -near_ - 2 * corner_ - (GLfloat)scaledRT_(it.getRT()));

              break;

            case SpectrumCanvas::IM_LOG:

              qglColor_(layer.gradient.precalculatedColorAt(0.0));
              glVertex3d(-corner_ + (GLfloat)scaledMZ_(it->getMZ()),
                         -corner_,
                         -near_ - 2 * corner_ - (GLfloat)scaledRT_(it.getRT()));
              qglColor_(layer.gradient.precalculatedColorAt(log10(1 + max(0.0, (double)(it->getIntensity())))));
              glVertex3d(-corner_ + (GLfloat)scaledMZ_(it->getMZ()),
                         -corner_ + (GLfloat)scaledIntensity_(it->getIntensity(), i),
                         -near_ - 2 * corner_ - (GLfloat)scaledRT_(it.getRT()));

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

  GLuint Spectrum3DOpenGLCanvas::makeGridLines_()
  {
    GLuint list = glGenLists(1);
    glNewList(list, GL_COMPILE);
    glEnable(GL_LINE_STIPPLE);
    glLineStipple(1, 0x0101);
    glBegin(GL_LINES);
    glColor4ub(0, 0, 0, 80);
    // mz
    if (grid_mz_.size() >= 1)
    {
      for (Size i = 0; i < grid_mz_[0].size(); i++)
      {
        glVertex3d(-corner_ + scaledMZ_(grid_mz_[0][i]), -corner_, -near_ - 2 * corner_);
        glVertex3d(-corner_ + scaledMZ_(grid_mz_[0][i]), -corner_, -far_ + 2 * corner_);
      }
    }
    if (grid_mz_.size() >= 2)
    {
      for (Size i = 0; i < grid_mz_[1].size(); i++)
      {
        glVertex3d(-corner_ + scaledMZ_(grid_mz_[1][i]), -corner_, -near_ - 2 * corner_);
        glVertex3d(-corner_ + scaledMZ_(grid_mz_[1][i]), -corner_, -far_ + 2 * corner_);
      }
    }
    if (grid_mz_.size() >= 3)
    {
      for (Size i = 0; i < grid_mz_[2].size(); i++)
      {
        glVertex3d(-corner_ + scaledMZ_(grid_mz_[2][i]), -corner_, -near_ - 2 * corner_);
        glVertex3d(-corner_ + scaledMZ_(grid_mz_[2][i]), -corner_, -far_ + 2 * corner_);
      }
    }
    // rt
    if (grid_rt_.size() >= 1)
    {
      for (Size i = 0; i < grid_rt_[0].size(); i++)
      {
        glVertex3d(-corner_, -corner_, -near_ - 2 * corner_ - scaledRT_(grid_rt_[0][i]));
        glVertex3d(corner_, -corner_, -near_ - 2 * corner_ - scaledRT_(grid_rt_[0][i]));
      }
    }
    if (grid_rt_.size() >= 2)
    {
      for (Size i = 0; i < grid_rt_[1].size(); i++)
      {
        glVertex3d(-corner_, -corner_, -near_ - 2 * corner_ - scaledRT_(grid_rt_[1][i]));
        glVertex3d(corner_, -corner_, -near_ - 2 * corner_ - scaledRT_(grid_rt_[1][i]));
      }
    }
    if (grid_rt_.size() >= 3)
    {
      for (Size i = 0; i < grid_rt_[2].size(); i++)
      {
        glVertex3d(-corner_, -corner_, -near_ - 2 * corner_ - scaledRT_(grid_rt_[2][i]));
        glVertex3d(corner_, -corner_, -near_ - 2 * corner_ - scaledRT_(grid_rt_[2][i]));
      }
    }
    glEnd();
    glDisable(GL_LINE_STIPPLE);
    glEndList();
    return list;
  }

  GLuint Spectrum3DOpenGLCanvas::makeAxesTicks_()
  {
    GLuint list = glGenLists(1);
    glNewList(list, GL_COMPILE);
    glShadeModel(GL_FLAT);
    glLineWidth(2.0);
    glBegin(GL_LINES);
    qglColor_(Qt::black);

    // mz
    if (grid_mz_.size() >= 1)
    {
      for (Size i = 0; i < grid_mz_[0].size(); i++)
      {
        glVertex3d(-corner_ + scaledMZ_(grid_mz_[0][i]), -corner_, -near_ - 2 * corner_);
        glVertex3d(-corner_ + scaledMZ_(grid_mz_[0][i]), -corner_ + 4.0, -near_ - 2 * corner_);
      }
    }
    if (grid_mz_.size() >= 2)
    {
      for (Size i = 0; i < grid_mz_[1].size(); i++)
      {
        glVertex3d(-corner_ + scaledMZ_(grid_mz_[1][i]), -corner_, -near_ - 2 * corner_);
        glVertex3d(-corner_ + scaledMZ_(grid_mz_[1][i]), -corner_ + 3.0, -near_ - 2 * corner_);
      }
    }
    if (grid_mz_.size() >= 3)
    {
      for (Size i = 0; i < grid_mz_[2].size(); i++)
      {
        glVertex3d(-corner_ + scaledMZ_(grid_mz_[2][i]), -corner_, -near_ - 2 * corner_);
        glVertex3d(-corner_ + scaledMZ_(grid_mz_[2][i]), -corner_ + 2.0, -near_ - 2 * corner_);
      }
    }

    // rt
    if (grid_rt_.size() >= 1)
    {
      for (Size i = 0; i < grid_rt_[0].size(); i++)
      {
        glVertex3d(-corner_, -corner_, -near_ - 2 * corner_ - scaledRT_(grid_rt_[0][i]));
        glVertex3d(-corner_, -corner_ + 4.0, -near_ - 2 * corner_ - scaledRT_(grid_rt_[0][i]));
      }
    }
    if (grid_rt_.size() >= 2)
    {
      for (Size i = 0; i < grid_rt_[1].size(); i++)
      {
        glVertex3d(-corner_, -corner_, -near_ - 2 * corner_ - scaledRT_(grid_rt_[1][i]));
        glVertex3d(-corner_, -corner_ + 3.0, -near_ - 2 * corner_ - scaledRT_(grid_rt_[1][i]));
      }
    }
    if (grid_rt_.size() >= 3)
    {
      for (Size i = 0; i < grid_rt_[2].size(); i++)
      {
        glVertex3d(-corner_, -corner_, -near_ - 2 * corner_ - scaledRT_(grid_rt_[2][i]));
        glVertex3d(-corner_, -corner_ + 2.0, -near_ - 2 * corner_ - scaledRT_(grid_rt_[2][i]));
      }
    }

    //Intensity
    switch (canvas_3d_.intensity_mode_)
    {
    case SpectrumCanvas::IM_PERCENTAGE:
      if (grid_intensity_.size() >= 1)
      {
        for (Size i = 0; i < grid_intensity_[0].size(); i++)
        {
          glVertex3d(-corner_, -corner_ + (2.0 * grid_intensity_[0][i]), -near_ - 2 * corner_);
          glVertex3d(-corner_ + 4.0, -corner_ + (2.0 * grid_intensity_[0][i]), -near_ - 2 * corner_ - 4.0);
        }
      }
      break;

    case SpectrumCanvas::IM_NONE:
    case SpectrumCanvas::IM_SNAP:
      if (grid_intensity_.size() >= 1)
      {
        for (Size i = 0; i < grid_intensity_[0].size(); i++)
        {
          glVertex3d(-corner_, -corner_ + scaledIntensity_(grid_intensity_[0][i], canvas_3d_.layers_.getCurrentLayerIndex()), -near_ - 2 * corner_);
          glVertex3d(-corner_ + 4.0, -corner_ + scaledIntensity_(grid_intensity_[0][i], canvas_3d_.layers_.getCurrentLayerIndex()), -near_ - 2 * corner_ - 4.0);
        }
      }
      if (grid_intensity_.size() >= 2)
      {
        for (Size i = 0; i < grid_intensity_[1].size(); i++)
        {
          glVertex3d(-corner_, -corner_ + scaledIntensity_(grid_intensity_[1][i], canvas_3d_.layers_.getCurrentLayerIndex()), -near_ - 2 * corner_);
          glVertex3d(-corner_ + 3.0, -corner_ + scaledIntensity_(grid_intensity_[1][i], canvas_3d_.layers_.getCurrentLayerIndex()), -near_ - 2 * corner_ - 3.0);
        }
      }
      if (grid_intensity_.size() >= 3)
      {
        for (Size i = 0; i < grid_intensity_[2].size(); i++)
        {
          glVertex3d(-corner_, -corner_ + scaledIntensity_(grid_intensity_[2][i], canvas_3d_.layers_.getCurrentLayerIndex()), -near_ - 2 * corner_);
          glVertex3d(-corner_ + 2.0, -corner_ + scaledIntensity_(grid_intensity_[2][i], canvas_3d_.layers_.getCurrentLayerIndex()), -near_ - 2 * corner_ - 2.0);
        }
      }
      break;

    case SpectrumCanvas::IM_LOG:
      if (grid_intensity_.size())
      {
        for (Size i = 0; i < grid_intensity_[0].size(); i++)
        {
          glVertex3d(-corner_, -corner_ + scaledIntensity_(pow(10.0, grid_intensity_[0][i]) - 1, canvas_3d_.layers_.getCurrentLayerIndex()), -near_ - 2 * corner_);
          glVertex3d(-corner_ + 4.0, -corner_ + scaledIntensity_(pow(10.0, grid_intensity_[0][i]) - 1, canvas_3d_.layers_.getCurrentLayerIndex()), -near_ - 2 * corner_ - 4.0);
        }
      }
      break;

    }
    glEnd();
    glEndList();
    return list;
  }

  double Spectrum3DOpenGLCanvas::scaledRT_(double rt)
  {
    double scaledrt = rt - canvas_3d_.visible_area_.min_[1];
    scaledrt = scaledrt * 2.0 * corner_ / (canvas_3d_.visible_area_.max_[1] - canvas_3d_.visible_area_.min_[1]);
    return scaledrt;
  }

  double Spectrum3DOpenGLCanvas::scaledInversRT_(double rt)
  {
    double i_rt = (rt * canvas_3d_.visible_area_.max_[1] - canvas_3d_.visible_area_.min_[1] * rt);
    i_rt = i_rt / 200.0;
    i_rt = i_rt + canvas_3d_.visible_area_.min_[1];
    // cout<<"rt"<<rt<<"  "<<"scaledinver"<<i_rt<<endl;
    return i_rt;
  }

  double Spectrum3DOpenGLCanvas::scaledMZ_(double mz)
  {
    double scaledmz = mz - canvas_3d_.visible_area_.min_[0];
    scaledmz = scaledmz * 2.0 * corner_ / (canvas_3d_.visible_area_.max_[0] - canvas_3d_.visible_area_.min_[0]) /*dis_mz_*/;
    return scaledmz;
  }

  double Spectrum3DOpenGLCanvas::scaledInversMZ_(double mz)
  {
    double i_mz = (mz * canvas_3d_.visible_area_.max_[0] - mz * canvas_3d_.visible_area_.min_[0]);
    i_mz = i_mz / 200;
    i_mz = i_mz + canvas_3d_.visible_area_.min_[0];
    return i_mz;
  }

  double Spectrum3DOpenGLCanvas::scaledIntensity_(float intensity, Size layer_index)
  {
    double scaledintensity = intensity * 2.0 * corner_;
    switch (canvas_3d_.intensity_mode_)
    {
    case SpectrumCanvas::IM_SNAP:
      scaledintensity /= int_scale_.max_[0];
      break;

    case SpectrumCanvas::IM_NONE:
      scaledintensity /= canvas_3d_.overall_data_range_.max_[2];
      break;

    case SpectrumCanvas::IM_PERCENTAGE:
      scaledintensity /= canvas_3d_.getMaxIntensity(layer_index);
      break;

    case SpectrumCanvas::IM_LOG:
      scaledintensity = log10(1 + max(0.0, (double)intensity)) * 2.0 * corner_ / log10(1 + max(0.0, canvas_3d_.overall_data_range_.max_[2]));
      break;
    }
    return scaledintensity;
  }

  void Spectrum3DOpenGLCanvas::normalizeAngle(int* angle)
  {
    while (*angle < 0) { *angle += 360 * 16; }
    while (*angle > 360 * 16) { *angle -= 360 * 16; }
  }

  ///////////////wheel- and MouseEvents//////////////////

  void Spectrum3DOpenGLCanvas::actionModeChange()
  {
    //change from translate to zoom
    if (canvas_3d_.action_mode_ == SpectrumCanvas::AM_ZOOM)
    {
      storeRotationAndZoom();
      xrot_ = 220;
      yrot_ = 220;
      zrot_ = 0;
      canvas_3d_.update_buffer_ = true;
      canvas_3d_.update_(OPENMS_PRETTY_FUNCTION);
    }
    //change from zoom to translate
    else if (canvas_3d_.action_mode_ == SpectrumCanvas::AM_TRANSLATE)
    {
      // if still in selection mode, quit selection mode first:
      if (canvas_3d_.rubber_band_.isVisible())
      {
        computeSelection_();
      }
      restoreRotationAndZoom();
      canvas_3d_.update_buffer_ = true;
      canvas_3d_.update_(OPENMS_PRETTY_FUNCTION);
    }
    update();
  }

  void Spectrum3DOpenGLCanvas::focusOutEvent(QFocusEvent * e)
  {
    canvas_3d_.focusOutEvent(e);
    update();
  }

  void Spectrum3DOpenGLCanvas::mousePressEvent(QMouseEvent * e)
  {
    mouse_move_begin_ = e->pos();
    mouse_move_end_ = e->pos();

    if (canvas_3d_.action_mode_ == SpectrumCanvas::AM_ZOOM && e->button() == Qt::LeftButton)
    {
      canvas_3d_.rubber_band_.setGeometry(QRect(e->pos(), QSize()));
      canvas_3d_.rubber_band_.show();
      canvas_3d_.update_buffer_ = true;
      canvas_3d_.update_(OPENMS_PRETTY_FUNCTION);
    }
    update();
  }

  void Spectrum3DOpenGLCanvas::mouseMoveEvent(QMouseEvent * e)
  {
    if (e->buttons() & Qt::LeftButton)
    {
      if (canvas_3d_.action_mode_ == SpectrumCanvas::AM_ZOOM)
      {
        canvas_3d_.rubber_band_.setGeometry(QRect(mouse_move_begin_, e->pos()).normalized());
        canvas_3d_.update_(OPENMS_PRETTY_FUNCTION);
      }
      else if (canvas_3d_.action_mode_ == SpectrumCanvas::AM_TRANSLATE)
      {
        Int x_angle = xrot_ + 8 * (e->y() - mouse_move_end_.y());
        normalizeAngle(&x_angle);
        xrot_ = x_angle;

        Int y_angle = yrot_ + 8 * (e->x() - mouse_move_end_.x());
        normalizeAngle(&y_angle);
        yrot_ = y_angle;

        //drawAxesLegend_();

        mouse_move_end_ = e->pos();
        canvas_3d_.update_(OPENMS_PRETTY_FUNCTION);
      }
    }
    update();
  }

  void Spectrum3DOpenGLCanvas::mouseReleaseEvent(QMouseEvent * e)
  {
    if (canvas_3d_.action_mode_ == SpectrumCanvas::AM_ZOOM && e->button() == Qt::LeftButton)
    {
      computeSelection_();
    }
    update();
  }

  void Spectrum3DOpenGLCanvas::computeSelection_()
  {
    QRect rect = canvas_3d_.rubber_band_.geometry();
    x_1_ = ((rect.topLeft().x() - width_ / 2) * corner_ * 1.25 * 2) / width_;
    y_1_ = -300 + (((rect.topLeft().y() - height_ / 2) * corner_ * 1.25 * 2) / height_);
    x_2_ = ((rect.bottomRight().x() - width_ / 2) * corner_ * 1.25 * 2) / width_;
    y_2_ = -300 + (((rect.bottomRight().y() - height_ / 2) * corner_ * 1.25 * 2) / height_);
    dataToZoomArray_(x_1_, y_1_, x_2_, y_2_);
    canvas_3d_.rubber_band_.hide();
    canvas_3d_.update_buffer_ = true;
    canvas_3d_.update_(OPENMS_PRETTY_FUNCTION);
  }

  void Spectrum3DOpenGLCanvas::dataToZoomArray_(double x_1, double y_1, double x_2, double y_2)
  {
    double scale_x1 = scaledInversMZ_(x_1 + 100.0);
    double scale_x2 = scaledInversMZ_(x_2 + 100.0);
    double scale_y1 = scaledInversRT_(-200 - y_1);
    double scale_y2 = scaledInversRT_(-200 - y_2);
    DRange<2> new_area_;
    if (scale_x1 <= scale_x2)
    {
      new_area_.min_[0] = scale_x1;
      new_area_.max_[0] = scale_x2;
    }
    else
    {
      new_area_.min_[0] = scale_x2;
      new_area_.max_[0] = scale_x1;
    }
    if (scale_y1 <= scale_y2)
    {
      new_area_.min_[1] = scale_y1;
      new_area_.max_[1] = scale_y2;
    }
    else
    {
      new_area_.min_[1] = scale_y2;
      new_area_.max_[1] = scale_y1;
    }
    canvas_3d_.changeVisibleArea_(new_area_, true, true);
  }

  void Spectrum3DOpenGLCanvas::updateIntensityScale()
  {
    int_scale_.min_[0] = canvas_3d_.overall_data_range_.max_[2];
    int_scale_.max_[0] = canvas_3d_.overall_data_range_.min_[2];

    for (Size i = 0; i < canvas_3d_.getLayerCount(); i++)
    {
      auto rt_begin_it = canvas_3d_.getLayer(i).getPeakData()->RTBegin(canvas_3d_.visible_area_.min_[1]);
      auto rt_end_it = canvas_3d_.getLayer(i).getPeakData()->RTEnd(canvas_3d_.visible_area_.max_[1]);

      for (auto spec_it = rt_begin_it; spec_it != rt_end_it; ++spec_it)
      {
        for (auto it = spec_it->MZBegin(canvas_3d_.visible_area_.min_[0]); it != spec_it->MZEnd(canvas_3d_.visible_area_.max_[0]); ++it)
        {
          if (int_scale_.min_[0] >= it->getIntensity()) { int_scale_.min_[0] = it->getIntensity(); }
          if (int_scale_.max_[0] <= it->getIntensity()) { int_scale_.max_[0] = it->getIntensity(); }
        }
      }
    }
  }

  void Spectrum3DOpenGLCanvas::recalculateDotGradient_(LayerData& layer)
  {
    layer.gradient.fromString(layer.param.getValue("dot:gradient"));
    switch (canvas_3d_.intensity_mode_)
    {
    case SpectrumCanvas::IM_SNAP:
      layer.gradient.activatePrecalculationMode(0.0, int_scale_.max_[0], UInt(canvas_3d_.param_.getValue("dot:interpolation_steps")));
      break;

    case SpectrumCanvas::IM_NONE:
      layer.gradient.activatePrecalculationMode(0.0, canvas_3d_.overall_data_range_.max_[2], UInt(canvas_3d_.param_.getValue("dot:interpolation_steps")));
      break;

    case SpectrumCanvas::IM_PERCENTAGE:
      layer.gradient.activatePrecalculationMode(0.0, 100.0, UInt(canvas_3d_.param_.getValue("dot:interpolation_steps")));
      break;

    case SpectrumCanvas::IM_LOG:
      layer.gradient.activatePrecalculationMode(0.0, log10(1 + max(0.0, canvas_3d_.overall_data_range_.max_[2])), UInt(canvas_3d_.param_.getValue("dot:interpolation_steps")));
      break;
    }
  }

} //end of namespace

