// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <iostream>

//OpenMS
#include <OpenMS/VISUAL/Spectrum3DWidget.h>
#include <OpenMS/VISUAL/Spectrum3DOpenGLCanvas.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum2DGoToDialog.h>

//QT
#include <QtGui/QGridLayout>

using namespace std;

namespace OpenMS
{
  using namespace Internal;
  using namespace Math;

  Spectrum3DWidget::Spectrum3DWidget(const Param & preferences, QWidget * parent) :
    SpectrumWidget(preferences, parent)
  {
    setCanvas_(new Spectrum3DCanvas(preferences, this));

    x_axis_->hide();
    y_axis_->hide();

    // delegate signals from canvas
    connect(canvas(), SIGNAL(showCurrentPeaksAs2D()), this, SIGNAL(showCurrentPeaksAs2D()));
  }

  Spectrum3DWidget::~Spectrum3DWidget()
  {

  }

  void Spectrum3DWidget::recalculateAxes_()
  {
  }

  Histogram<> Spectrum3DWidget::createIntensityDistribution_() const
  {
    //initialize histogram
    DoubleReal min = canvas_->getCurrentMinIntensity();
    DoubleReal max = canvas_->getCurrentMaxIntensity();
    if (min == max)
    {
      min -= 0.01;
      max += 0.01;
    }
    Histogram<> tmp(min, max, (max - min) / 500.0);

    for (ExperimentType::ConstIterator spec_it = canvas_->getCurrentLayer().getPeakData()->begin(); spec_it != canvas_->getCurrentLayer().getPeakData()->end(); ++spec_it)
    {
      if (spec_it->getMSLevel() != 1)
        continue;
      for (ExperimentType::SpectrumType::ConstIterator peak_it = spec_it->begin(); peak_it != spec_it->end(); ++peak_it)
      {
        tmp.inc(peak_it->getIntensity());
      }
    }

    return tmp;
  }

  Histogram<> Spectrum3DWidget::createMetaDistribution_(const String & name) const
  {
    Histogram<> tmp;

    //determine min and max of the data
    Real m_min = (numeric_limits<Real>::max)(), m_max = -(numeric_limits<Real>::max)();
    for (ExperimentType::const_iterator s_it = canvas_->getCurrentLayer().getPeakData()->begin(); s_it != canvas_->getCurrentLayer().getPeakData()->end(); ++s_it)
    {
      if (s_it->getMSLevel() != 1)
        continue;
      //float arrays
      for (ExperimentType::SpectrumType::FloatDataArrays::const_iterator it = s_it->getFloatDataArrays().begin(); it != s_it->getFloatDataArrays().end(); ++it)
      {
        if (it->getName() == name)
        {
          for (Size i = 0; i < it->size(); ++i)
          {
            if ((*it)[i] < m_min)
              m_min = (*it)[i];
            if ((*it)[i] > m_max)
              m_max = (*it)[i];
          }
          break;
        }
      }
      //integer arrays
      for (ExperimentType::SpectrumType::IntegerDataArrays::const_iterator it = s_it->getIntegerDataArrays().begin(); it != s_it->getIntegerDataArrays().end(); ++it)
      {
        if (it->getName() == name)
        {
          for (Size i = 0; i < it->size(); ++i)
          {
            if ((*it)[i] < m_min)
              m_min = (*it)[i];
            if ((*it)[i] > m_max)
              m_max = (*it)[i];
          }
          break;
        }
      }
    }
    if (m_min >= m_max)
      return tmp;

    //create histogram
    tmp.reset(m_min, m_max, (m_max - m_min) / 500.0);
    for (ExperimentType::const_iterator s_it = canvas_->getCurrentLayer().getPeakData()->begin(); s_it != canvas_->getCurrentLayer().getPeakData()->end(); ++s_it)
    {
      if (s_it->getMSLevel() != 1)
        continue;
      //float arrays
      for (ExperimentType::SpectrumType::FloatDataArrays::const_iterator it = s_it->getFloatDataArrays().begin(); it != s_it->getFloatDataArrays().end(); ++it)
      {
        if (it->getName() == name)
        {
          for (Size i = 0; i < it->size(); ++i)
          {
            tmp.inc((*it)[i]);
          }
          break;
        }
      }
      //integer arrays
      for (ExperimentType::SpectrumType::IntegerDataArrays::const_iterator it = s_it->getIntegerDataArrays().begin(); it != s_it->getIntegerDataArrays().end(); ++it)
      {
        if (it->getName() == name)
        {
          for (Size i = 0; i < it->size(); ++i)
          {
            tmp.inc((*it)[i]);
          }
          break;
        }
      }
    }

    return tmp;
  }

  void Spectrum3DWidget::showLegend(bool show)
  {
    canvas()->showLegend(show);
  }

  bool Spectrum3DWidget::isLegendShown() const
  {
    return static_cast<const Spectrum3DCanvas *>(canvas_)->isLegendShown();
  }

  void Spectrum3DWidget::showGoToDialog()
  {
    Spectrum2DGoToDialog goto_dialog(this);
    const DRange<3> & area = canvas()->getDataRange();
    goto_dialog.setRange(area.minY(), area.maxY(), area.minX(), area.maxX());
    goto_dialog.enableFeatureNumber(false);
    if (goto_dialog.exec())
    {
      goto_dialog.fixRange(); // in case user did something invalid
      canvas()->setVisibleArea(SpectrumCanvas::AreaType(goto_dialog.getMinMZ(), goto_dialog.getMinRT(), goto_dialog.getMaxMZ(), goto_dialog.getMaxRT()));
    }
  }

} //namespace
