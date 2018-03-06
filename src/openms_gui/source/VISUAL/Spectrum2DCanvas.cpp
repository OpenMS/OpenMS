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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum2DPrefDialog.h>
#include <OpenMS/VISUAL/ColorSelector.h>
#include <OpenMS/VISUAL/MultiGradientSelector.h>
#include <OpenMS/VISUAL/DIALOGS/FeatureEditDialog.h>
#include <OpenMS/SYSTEM/FileWatcher.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
//STL
#include <algorithm>

//QT
#include <QtGui/QMouseEvent>
#include <QtGui/QPainter>
#include <QtGui/QMenu>
#include <QtGui/QBitmap>
#include <QtGui/QPolygon>
#include <QtCore/QTime>
#include <QtGui/QComboBox>
#include <QtGui/QFileDialog>
#include <QtGui/QMessageBox>

//boost
#include <boost/math/special_functions/fpclassify.hpp>

#define PEN_SIZE_MAX_LIMIT 100    // maximum size of a rectangle representing a point for raw peak data
#define PEN_SIZE_MIN_LIMIT 1      // minimum. This should not be changed without adapting the way dots are plotted 
                                  // (might lead to inconsistencies when switching between drawing modes of
                                  //  paintMaximumIntensities() vs. paintAllIntensities() )

// the two following constants describe the valid range of the 'canvas_coverage_min_' member (adaptable by the user)
#define CANVAS_COVERAGE_MIN_LIMITHIGH 0.5
#define CANVAS_COVERAGE_MIN_LIMITLOW 0.1


using namespace std;

namespace OpenMS
{
  using namespace Internal;

  Spectrum2DCanvas::Spectrum2DCanvas(const Param & preferences, QWidget * parent) :
    SpectrumCanvas(preferences, parent),
    projection_mz_(),
    projection_rt_(),
    selected_peak_(),
    measurement_start_(),
    pen_size_min_(1),
    pen_size_max_(20),
    canvas_coverage_min_(0.2)
  {
    //Parameter handling
    defaults_.setValue("background_color", "#ffffff", "Background color.");
    defaults_.setValue("interpolation_steps", 1000, "Number of interpolation steps for peak gradient pre-calculation.");
    defaults_.setMinInt("interpolation_steps", 1);
    defaults_.setMaxInt("interpolation_steps", 1000);
    defaults_.setValue("dot:gradient", "Linear|0,#eeeeee;1,#ffea00;6,#ff0000;14,#aa00ff;23,#5500ff;100,#000000", "Multi-color gradient for peaks.");
    defaults_.setValue("dot:feature_icon", "circle", "Icon used for features and consensus features.");
    defaults_.setValidStrings("dot:feature_icon", ListUtils::create<String>("diamond,square,circle,triangle"));
    defaults_.setValue("dot:feature_icon_size", 4, "Icon size used for features and consensus features.");
    defaults_.setMinInt("dot:feature_icon_size", 1);
    defaults_.setMaxInt("dot:feature_icon_size", 999);
    defaults_.setValue("mapping_of_mz_to", "y_axis", "Determines which axis is the m/z axis.");
    defaults_.setValidStrings("mapping_of_mz_to", ListUtils::create<String>("x_axis,y_axis"));
    defaultsToParam_();
    setName("Spectrum2DCanvas");
    setParameters(preferences);

    linear_gradient_.fromString(param_.getValue("dot:gradient"));

    projection_mz_.resize(1);
    projection_rt_.resize(1);

    //set preferences and update widgets accordingly
    if (String(param_.getValue("mapping_of_mz_to")) != "x_axis")
    {
      mzToXAxis(false);
    }
    //connect preferences change to the right slot
    connect(this, SIGNAL(preferencesChange()), this, SLOT(currentLayerParametersChanged_()));
  }

  Spectrum2DCanvas::~Spectrum2DCanvas()
  {
  }

  void Spectrum2DCanvas::highlightPeak_(QPainter & painter, const PeakIndex & peak)
  {
    if (!peak.isValid())
      return;

    //determine coordinates;
    QPoint pos;
    if (getCurrentLayer().type == LayerData::DT_FEATURE)
    {
      dataToWidget_(peak.getFeature(*getCurrentLayer().getFeatureMap()).getMZ(), peak.getFeature(*getCurrentLayer().getFeatureMap()).getRT(), pos);
    }
    else if (getCurrentLayer().type == LayerData::DT_PEAK)
    {
      dataToWidget_(peak.getPeak(*getCurrentLayer().getPeakData()).getMZ(), peak.getSpectrum(*getCurrentLayer().getPeakData()).getRT(), pos);
    }
    else if (getCurrentLayer().type == LayerData::DT_CONSENSUS)
    {
      dataToWidget_(peak.getFeature(*getCurrentLayer().getConsensusMap()).getMZ(), peak.getFeature(*getCurrentLayer().getConsensusMap()).getRT(), pos);
    }
    else if (getCurrentLayer().type == LayerData::DT_CHROMATOGRAM)
    {
      const LayerData & layer = getCurrentLayer();
      const ExperimentSharedPtrType exp = layer.getPeakData();

      // create iterator on chromatogram spectrum passed by PeakIndex
      vector<MSChromatogram >::const_iterator chrom_it = exp->getChromatograms().begin();
      chrom_it += peak.spectrum;
      dataToWidget_(chrom_it->getPrecursor().getMZ(), chrom_it->front().getRT(), pos);
    }
    else if (getCurrentLayer().type == LayerData::DT_IDENT)
    {
      //TODO IDENT
    }

    // paint highlighted peak
    painter.save();
    painter.setPen(QPen(Qt::red, 2));

    if (getCurrentLayer().type == LayerData::DT_CHROMATOGRAM) // highlight: chromatogram
    {
      const LayerData & layer = getCurrentLayer();
      const ExperimentSharedPtrType exp = layer.getPeakData();

      vector<MSChromatogram >::const_iterator iter = exp->getChromatograms().begin();
      iter += peak.spectrum;

      painter.drawRect(pos.x() - 5, pos.y() - 5, (int((iter->back().getRT() - iter->front().getRT()) / visible_area_.height() * width())) + 10, 10);
    }
    else // highlight: peak, feature, consensus feature
    {
      painter.drawEllipse(pos.x() - 5, pos.y() - 5, 10, 10);
    }

    //restore painter
    painter.restore();
  }

  PeakIndex Spectrum2DCanvas::findNearestPeak_(const QPoint & pos)
  {
    ///no layers => return invalid peak index
    if (layers_.empty())
      return PeakIndex();

    //Constructing the area corrects swapped mapping of RT and m/z
    AreaType area(widgetToData_(pos - QPoint(5, 5)), widgetToData_(pos + QPoint(5, 5)));

    float max_int = -1 * numeric_limits<float>::max();
    PeakIndex max_pi;

    if (getCurrentLayer().type == LayerData::DT_PEAK)
    {
      for (ExperimentType::ConstAreaIterator i = getCurrentLayer().getPeakData()->areaBeginConst(area.minPosition()[1], area.maxPosition()[1],
                                                                                                 area.minPosition()[0], area.maxPosition()[0]);
           i != getCurrentLayer().getPeakData()->areaEndConst();
           ++i)
      {
        PeakIndex pi = i.getPeakIndex();
        if (i->getIntensity() > max_int && getCurrentLayer().filters.passes((*getCurrentLayer().getPeakData())[pi.spectrum], pi.peak))
        {
          //cout << "new max: " << i.getRT() << " " << i->getMZ() << endl;
          max_int = i->getIntensity();
          max_pi = pi;
        }
      }
    }
    else if (getCurrentLayer().type == LayerData::DT_FEATURE)
    {
      for (FeatureMapType::ConstIterator i = getCurrentLayer().getFeatureMap()->begin();
           i != getCurrentLayer().getFeatureMap()->end();
           ++i)
      {
        if (i->getRT() >= area.minPosition()[1] &&
            i->getRT() <= area.maxPosition()[1] &&
            i->getMZ() >= area.minPosition()[0] &&
            i->getMZ() <= area.maxPosition()[0] &&
            getCurrentLayer().filters.passes(*i))
        {
          if (i->getIntensity() > max_int)
          {
            max_int = i->getIntensity();
            max_pi = PeakIndex(i - getCurrentLayer().getFeatureMap()->begin());
          }
        }
      }
    }
    else if (getCurrentLayer().type == LayerData::DT_CONSENSUS)
    {
      for (ConsensusMapType::ConstIterator i = getCurrentLayer().getConsensusMap()->begin();
           i != getCurrentLayer().getConsensusMap()->end();
           ++i)
      {
        // consensus feature in visible area?
        if (i->getRT() >= area.minPosition()[1] &&
            i->getRT() <= area.maxPosition()[1] &&
            i->getMZ() >= area.minPosition()[0] &&
            i->getMZ() <= area.maxPosition()[0] &&
            getCurrentLayer().filters.passes(*i))
        {
          if (i->getIntensity() > max_int)
          {
            max_int = i->getIntensity();
            max_pi = PeakIndex(i - getCurrentLayer().getConsensusMap()->begin());
          }
        }
      }
    }
    else if (getCurrentLayer().type == LayerData::DT_CHROMATOGRAM)
    {
      const LayerData & layer = getCurrentLayer();

      PeakMap exp;
      exp = *layer.getPeakData();
      float mz_origin = 0.0;

      for (vector<MSChromatogram >::const_iterator iter = exp.getChromatograms().begin(); iter != exp.getChromatograms().end(); ++iter)
      {
        if (iter->empty()) {continue;} // ensure that empty chromatograms are not examined (iter->front = segfault)
        MSChromatogram::ConstIterator cit = iter->begin();

        // for (MSChromatogram::ConstIterator cit = iter->begin(); cit != iter->end(); ++cit)
        // {
        //   cout << "Chrom Values RT/INT: " << cit->getRT() << "/" << " " << cit->getIntensity() << endl;
        //  }

        if (mz_origin != iter->getPrecursor().getMZ())
        {
          mz_origin = iter->getPrecursor().getMZ();
        }

        for (int i = int(iter->front().getRT()); i <= int(iter->back().getRT()); i++)
        {
          if (i >= area.minPosition()[1] &&
              i <= area.maxPosition()[1] &&
              mz_origin >= area.minPosition()[0] &&
              mz_origin <= area.maxPosition()[0])
          {
            return PeakIndex(iter - exp.getChromatograms().begin(), cit - iter->begin());
          }
        }
      }
    }
    else if (getCurrentLayer().type == LayerData::DT_IDENT)
    {
      //TODO IDENT
    }

    return max_pi;
  }

  void Spectrum2DCanvas::paintDots_(Size layer_index, QPainter & painter)
  {
    const LayerData & layer = getLayer(layer_index);

    //update factors (snap and percentage)
    double snap_factor = snap_factors_[layer_index];
    percentage_factor_ = 1.0;
    if (intensity_mode_ == IM_PERCENTAGE)
    {
      if (layer.type == LayerData::DT_PEAK && layer.getPeakData()->getMaxInt() > 0.0)
      {
        percentage_factor_ = overall_data_range_.maxPosition()[2] / layer.getPeakData()->getMaxInt();
      }
      else if (layer.type == LayerData::DT_FEATURE && layer.getFeatureMap()->getMaxInt() > 0.0)
      {
        percentage_factor_ = overall_data_range_.maxPosition()[2] / layer.getFeatureMap()->getMaxInt();
      }
      else if (layer.type == LayerData::DT_CONSENSUS && layer.getConsensusMap()->getMaxInt() > 0.0)
      {
        percentage_factor_ = overall_data_range_.maxPosition()[2] / layer.getConsensusMap()->getMaxInt();
      }
      else if (layer.type == LayerData::DT_CHROMATOGRAM && layer.getConsensusMap()->getMaxInt() > 0.0)
      {
        //TODO CHROM not sure if needed here
      }
    }

    //temporary variables
    Int image_width = buffer_.width();
    Int image_height = buffer_.height();

    if (layer.type == LayerData::DT_PEAK)   //peaks
    {
      // renaming some values for readability
      const ExperimentType & peak_map = *layer.getPeakData();
      const double rt_min = visible_area_.minPosition()[1];
      const double rt_max = visible_area_.maxPosition()[1];
      const double mz_min = visible_area_.minPosition()[0];
      const double mz_max = visible_area_.maxPosition()[0];

      // skip empty peak maps
      if (peak_map.empty())
      {
        return;
      }

      //determine number of pixels for each dimension
      Size rt_pixel_count = image_height;
      Size mz_pixel_count = image_width;
      if (!isMzToXAxis())
      {
        swap(rt_pixel_count, mz_pixel_count);
      }

      //-----------------------------------------------------------------------------------------------
      // Determine number of shown scans (MS1)
      std::vector<Size> rt_indices; // list of visible RT scans in MS1 with at least 2 points
      for (ExperimentType::ConstIterator it = peak_map.RTBegin(rt_min); it != peak_map.RTEnd(rt_max); ++it)
      {
        if (it->getMSLevel() == 1 && it->size() > 1)
        {
          rt_indices.push_back(std::distance(peak_map.begin(), it));
        }
      }
      Size n_ms1_scans = rt_indices.size();

      if (n_ms1_scans > 0)
      {
        // sample #of points at 3 scan locations (25%, 50%, 75%)
        // and take the median value
        Size n_peaks_in_scan(0);
        {
          double quantiles[] = {0.25, 0.50, 0.75};
          std::vector<Size> n_s;
          for (Size i=0; i<sizeof(quantiles)/sizeof(double); ++i)
          {
            const ExperimentType::SpectrumType& spec = peak_map[rt_indices[n_ms1_scans*quantiles[i]]];
            n_s.push_back(std::distance(spec.MZBegin(mz_min), spec.MZEnd(mz_max)) + 1); // +1 to since distance is 0 if only one m/z is shown
          } 
          std::sort(n_s.begin(), n_s.end());
          n_peaks_in_scan = n_s[1]; // median
        }
      
        double ratio_data2pixel_rt = n_ms1_scans / (double)rt_pixel_count;
        double ratio_data2pixel_mz = n_peaks_in_scan / (double)mz_pixel_count;

  #ifdef DEBUG_TOPPVIEW
        std::cerr << rt_min << ":" << rt_max <<"   " << mz_min << ":" << mz_max << "\n";
        std::cerr << n_ms1_scans << "rt " << n_peaks_in_scan << "ms\n";
        std::cerr << rt_pixel_count << " x " << mz_pixel_count << " px\n";
        std::cerr << ratio_data2pixel_rt << " " << ratio_data2pixel_mz << " ratio\n";
  #endif
        // minimum fraction of image expected to be filled with data
        // if not reached, we upscale point size
        bool has_low_pixel_coverage = ratio_data2pixel_rt < canvas_coverage_min_ || ratio_data2pixel_mz < canvas_coverage_min_;

        // Are several peaks expected to be drawn on the same pixel in either RT or m/z?
        // --> thin out and show only maxima
        // Also, we cannot upscale in this mode (since we operate on the buffer directly, i.e. '1 data point == 1 pixel'
        if (!has_low_pixel_coverage && (n_peaks_in_scan > mz_pixel_count || n_ms1_scans > rt_pixel_count))
        {
          paintMaximumIntensities_(layer_index, rt_pixel_count, mz_pixel_count, painter);
        }
        else
        { // this is slower to paint, but allows scaling points
          // when data is zoomed in to single peaks these are visualized as circles

          // compute ideal pen width (from data);
          // since points are rectangular, we take the value of the most "crowded" dimension
          // i.e. so that adjacent points do not overlap      
          double pen_width = std::min(1/ratio_data2pixel_rt, 1/ratio_data2pixel_mz);
          // ... and make sure its within our boundaries
          pen_width = std::max(pen_width, pen_size_min_);
          pen_width = std::min(pen_width, pen_size_max_);
  #ifdef DEBUG_TOPPVIEW
          std::cerr << "pen-width " << pen_width << "\n";
  #endif
          // However: if one dimension is sparse (e.g. only a few, but very long scans), we want to
          //          avoid showing lots of white background by increasing point size
          // This might lead to 'overplotting', but the paint method below can deal with it since
          // it will paint high intensities last.
          adaptPenScaling_(ratio_data2pixel_mz, pen_width);
          adaptPenScaling_(ratio_data2pixel_rt, pen_width);
  #ifdef DEBUG_TOPPVIEW
          std::cerr << "new pen: " << pen_width << "\n";
  #endif

          // few data points expected: more expensive drawing of all data points (circles or points depending on zoom level)
          paintAllIntensities_(layer_index, pen_width, painter);
        }

      } // end of no-scans check

      //-----------------------------------------------------------------
      //draw precursor peaks
      if (getLayerFlag(layer_index, LayerData::P_PRECURSORS))
      {
        paintPrecursorPeaks_(layer_index, painter);
      }
    }
    else if (layer.type == LayerData::DT_FEATURE)   //features
    {
      paintFeatureData_(layer_index, painter);
    }
    else if (layer.type == LayerData::DT_CONSENSUS)  // consensus features
    {
      String icon = layer.param.getValue("dot:feature_icon");
      Size icon_size = layer.param.getValue("dot:feature_icon_size");

      for (ConsensusMapType::ConstIterator i = layer.getConsensusMap()->begin();
           i != layer.getConsensusMap()->end();
           ++i)
      {
        if (i->getRT() >= visible_area_.minPosition()[1] &&
            i->getRT() <= visible_area_.maxPosition()[1] &&
            i->getMZ() >= visible_area_.minPosition()[0] &&
            i->getMZ() <= visible_area_.maxPosition()[0] &&
            layer.filters.passes(*i))
        {
          // determine color
          QColor color;
          if (i->metaValueExists(5))
          {
            color = QColor(i->getMetaValue(5).toQString());
          }
          else
          {
            // use intensity as color
            color = heightColor_(i->getIntensity(), layer.gradient, snap_factor);
          }

          // paint
          QPoint pos;
          dataToWidget_(i->getMZ(), i->getRT(), pos);
          if (pos.x() > 0 && pos.y() > 0 && pos.x() < image_width - 1 && pos.y() < image_height - 1)
          {
            paintIcon_(pos, color.rgb(), icon, icon_size, painter);
          }
        }
      }
    }
    else if (layer.type == LayerData::DT_CHROMATOGRAM) // chromatograms
    {
      const PeakMap exp = *layer.getPeakData();
      //TODO CHROM implement layer filters
      //TODO CHROM implement faster painting

      // paint chromatogram rt start and end as line
      float mz_origin = 0;
      float min_rt = 0;
      float max_rt = 0;

      for (vector<MSChromatogram >::const_iterator iter = exp.getChromatograms().begin(); iter != exp.getChromatograms().end(); ++iter)
      {
        if (mz_origin != iter->getPrecursor().getMZ())
        {
          mz_origin = iter->getPrecursor().getMZ();
          if (!iter->empty())
          {
            min_rt = iter->front().getRT();
            max_rt = iter->back().getRT();
          }
        }

        QPoint posi;
        QPoint posi2;

        dataToWidget_(iter->getPrecursor().getMZ(), min_rt, posi);
        dataToWidget_(iter->getPrecursor().getMZ(), max_rt, posi2);

        painter.drawLine(posi.x(), posi.y(), posi2.x(), posi2.y());
      }
    }
    else if (layer.type == LayerData::DT_IDENT)   // peptide identifications
    {
      paintIdentifications_(layer_index, painter);
    }
  }

  double Spectrum2DCanvas::adaptPenScaling_(double ratio_data2pixel, double& pen_width) const
  {
    // is the coverage ok using current pen width?
    bool has_low_pixel_coverage_withpen = ratio_data2pixel*pen_width < canvas_coverage_min_;
    int merge_factor(1);
    if (has_low_pixel_coverage_withpen)
    { // scale up the sparse dimension until we reach the desired coverage (this will lead to overlap in the crowded dimension)
      double scale_factor = canvas_coverage_min_ / ratio_data2pixel;
      // however, within bounds (no need to check the pen_size_min_, because we can only exceed here, not underestimate)
      scale_factor = std::min(pen_size_max_, scale_factor);
      // The difference between the original pen_width vs. this scale 
      // gives the number of peaks to merge in the crowded dimension
      merge_factor = scale_factor / pen_width;
      // set pen width to the new scale
      pen_width = scale_factor;
    }
    return merge_factor;
  }


  void Spectrum2DCanvas::paintPrecursorPeaks_(Size layer_index, QPainter & painter)
  {
    const LayerData & layer = getLayer(layer_index);
    const ExperimentType & peak_map = *layer.getPeakData();

    QPoint pos_ms1;
    QPoint pos_ms2;
    QPen p;
    p.setColor(Qt::black);
    painter.setPen(p);

    ExperimentType::ConstIterator it_prec = peak_map.end();
    ExperimentType::ConstIterator it_end = peak_map.RTEnd(visible_area_.maxPosition()[1]);
    for (ExperimentType::ConstIterator it = peak_map.RTBegin(visible_area_.minPosition()[1]);
         it != it_end;
         ++it)
    {
      // remember last precursor spectrum (do not call it->getPrecursorSpectrum(), since it can be very slow if no MS1 data is present)
      if (it->getMSLevel() == 1)
      {
        it_prec = it;
      }
      else if (it->getMSLevel() == 2 && !it->getPrecursors().empty())
      { // this is an MS/MS scan
        dataToWidget_(it->getPrecursors()[0].getMZ(), it->getRT(), pos_ms2);   // position of precursor in MS2
        const int x2 = pos_ms2.x();
        const int y2 = pos_ms2.y();

        if (it_prec != peak_map.end())
        {
          dataToWidget_(it->getPrecursors()[0].getMZ(), it_prec->getRT(), pos_ms1);  // position of precursor in MS1
          const int x = pos_ms1.x();
          const int y = pos_ms1.y();
          // diamond shape in MS1
          painter.drawLine(x, y + 3, x + 3, y);
          painter.drawLine(x + 3, y, x, y - 3);
          painter.drawLine(x, y - 3, x - 3, y);
          painter.drawLine(x - 3, y, x, y + 3);

          // rt position of corresponding MS2
          painter.drawLine(x2 - 3, y2, x2 + 3, y2);
          painter.drawLine(x, y, x2, y2);
        }
        else // no preceding MS1
        {
          // rt position of corresponding MS2 (cross)
          painter.drawLine(x2 - 3, y2, x2 + 3, y2);
          painter.drawLine(x2, y2 - 3, x2, y2 + 3);
        }
      }
    }
  }

  void Spectrum2DCanvas::paintAllIntensities_(Size layer_index, double pen_width, QPainter & painter)
  {
    const LayerData & layer = getLayer(layer_index);
    Int image_width = buffer_.width();
    Int image_height = buffer_.height();
    QVector<QPolygon> coloredPoints( (int)layer.gradient.precalculatedSize() );

    const ExperimentType & map = *layer.getPeakData();
    const double rt_min = visible_area_.minPosition()[1];
    const double rt_max = visible_area_.maxPosition()[1];
    const double mz_min = visible_area_.minPosition()[0];
    const double mz_max = visible_area_.maxPosition()[0];

    double snap_factor = snap_factors_[layer_index];

    for (ExperimentType::ConstAreaIterator i = map.areaBeginConst(rt_min, rt_max, mz_min, mz_max);
         i != map.areaEndConst();
         ++i)
    {
      PeakIndex pi = i.getPeakIndex();
      if (layer.filters.passes(map[pi.spectrum], pi.peak))
      {
        QPoint pos;
        dataToWidget_(i->getMZ(), i.getRT(), pos);
        if (pos.x() > 0 && pos.y() > 0 && pos.x() < image_width - 1 && pos.y() < image_height - 1)
        {
          // store point in the array of its color
          Int colorIndex = precalculatedColorIndex_(i->getIntensity(), layer.gradient, snap_factor);
          coloredPoints[ colorIndex ].push_back( pos );
        }
      }
    }
    // draw point arrays from minimum to maximum intensity,
    // avoiding low-intensity points obscuring the high-intensity ones
    painter.save();
    QPen newPointsPen;
    newPointsPen.setWidthF( pen_width );
    for ( Int colorIx = 0; colorIx < coloredPoints.size(); colorIx++ ) {
        const QPolygon& pointsArr = coloredPoints[colorIx];
        if ( pointsArr.size() ) {
            newPointsPen.setColor( layer.gradient.precalculatedColorByIndex( colorIx ) );
            painter.setPen( newPointsPen );
            painter.drawPoints( pointsArr );
        }
    }
    painter.restore();
  }

  void Spectrum2DCanvas::paintMaximumIntensities_(Size layer_index, Size rt_pixel_count, Size mz_pixel_count, QPainter & painter)
  {
    //set painter to black (we operate directly on the pixels for all colored data)
    painter.setPen(Qt::black);
    //temporary variables
    Int image_width = buffer_.width();
    Int image_height = buffer_.height();

    const LayerData & layer = getLayer(layer_index);
    const ExperimentType & map = *layer.getPeakData();
    const double rt_min = visible_area_.minPosition()[1];
    const double rt_max = visible_area_.maxPosition()[1];
    const double mz_min = visible_area_.minPosition()[0];
    const double mz_max = visible_area_.maxPosition()[0];

    double snap_factor = snap_factors_[layer_index];

    //calculate pixel size in data coordinates
    double rt_step_size = (rt_max - rt_min) / rt_pixel_count;
    double mz_step_size = (mz_max - mz_min) / mz_pixel_count;

    // start at first visible RT scan
    Size scan_index = std::distance(map.begin(), map.RTBegin(rt_min));
    //iterate over all pixels (RT dimension)
    for (Size rt = 0; rt < rt_pixel_count; ++rt)
    {
      // interval in data coordinates for the current pixel
      double rt_start = rt_min + rt_step_size * rt;
      double rt_end = rt_start + rt_step_size;
      //cout << "rt: " << rt << " (" << rt_start << " - " << rt_end << ")" << endl;

      // reached the end of data
      if (rt_end >= (--map.end())->getRT()) break;

      //determine the relevant spectra and reserve an array for the peak indices
      vector<Size> scan_indices, peak_indices;
      for (Size i = scan_index; i < map.size(); ++i)
      {
        if (map[i].getRT() >= rt_end)
        {
          scan_index = i;   //store last scan index for next RT pixel
          break;
        }
        if (map[i].getMSLevel() == 1 && map[i].size() > 0)
        {
          scan_indices.push_back(i);
          peak_indices.push_back(map[i].MZBegin(mz_min) - map[i].begin());
        }
      }
      //cout << "  scans: " << scan_indices.size() << endl;

      if (scan_indices.empty())
        continue;

      //iterate over all pixels (m/z dimension)
      for (Size mz = 0; mz < mz_pixel_count; ++mz)
      {
        double mz_start = mz_min + mz_step_size * mz;
        double mz_end = mz_start + mz_step_size;

        //iterate over all relevant peaks in all relevant scans
        float max = -1.0;
        for (Size i = 0; i < scan_indices.size(); ++i)
        {
          Size s = scan_indices[i];
          Size p = peak_indices[i];
          for (; p < map[s].size(); ++p)
          {
            if (map[s][p].getMZ() >= mz_end)
              break;
            if (map[s][p].getIntensity() > max && layer.filters.passes(map[s], p))
            {
              max = map[s][p].getIntensity();
            }
          }
          peak_indices[i] = p;   //store last peak index for next m/z pixel
        }

        //draw to buffer
        if (max >= 0.0)
        {
          QPoint pos;
          dataToWidget_(mz_start + 0.5 * mz_step_size, rt_start + 0.5 * rt_step_size, pos);
          if (pos.y() < image_height && pos.x() < image_width)
          {
            buffer_.setPixel(pos.x(), pos.y(), heightColor_(max, layer.gradient, snap_factor).rgb());
          }
        }
      }
    }
  }

  void Spectrum2DCanvas::paintFeatureData_(Size layer_index, QPainter& painter)
  {
    const LayerData& layer = getLayer(layer_index);
    double snap_factor = snap_factors_[layer_index];
    Int image_width = buffer_.width();
    Int image_height = buffer_.height();

    int line_spacing = QFontMetrics(painter.font()).lineSpacing();
    String icon = layer.param.getValue("dot:feature_icon");
    Size icon_size = layer.param.getValue("dot:feature_icon_size");
    bool show_label = (layer.label != LayerData::L_NONE);
    UInt num = 0;
    for (FeatureMapType::ConstIterator i = layer.getFeatureMap()->begin();
         i != layer.getFeatureMap()->end(); ++i)
    {
      if (i->getRT() >= visible_area_.minPosition()[1] &&
          i->getRT() <= visible_area_.maxPosition()[1] &&
          i->getMZ() >= visible_area_.minPosition()[0] &&
          i->getMZ() <= visible_area_.maxPosition()[0] &&
          layer.filters.passes(*i))
      {
        // determine color
        QColor color;
        if (i->metaValueExists(5))
        {
          color = QColor(i->getMetaValue(5).toQString());
        }
        else
        {
          color = heightColor_(i->getIntensity(), layer.gradient, snap_factor);
        }
        // paint
        QPoint pos;
        dataToWidget_(i->getMZ(), i->getRT(), pos);
        if (pos.x() > 0 && pos.y() > 0 && pos.x() < image_width - 1 && pos.y() < image_height - 1)
        {
          paintIcon_(pos, color.rgb(), icon, icon_size, painter);
        }
        // labels
        if (show_label)
        {
          if (layer.label == LayerData::L_INDEX)
          {
            painter.setPen(Qt::darkBlue);
            painter.drawText(pos.x() + 10, pos.y() + 10, QString::number(num));
          }
          else if ((layer.label == LayerData::L_ID || layer.label == LayerData::L_ID_ALL) && !i->getPeptideIdentifications().empty() && !i->getPeptideIdentifications()[0].getHits().empty())
          {
            painter.setPen(Qt::darkGreen);
            Size maxHits = (layer.label == LayerData::L_ID_ALL) ? i->getPeptideIdentifications()[0].getHits().size() : 1;
            for (Size j = 0; j < maxHits; ++j)
            {
              painter.drawText(pos.x() + 10, pos.y() + 10 + int(j) * line_spacing, i->getPeptideIdentifications()[0].getHits()[j].getSequence().toString().toQString());
            }
          }
          else if (layer.label == LayerData::L_META_LABEL)
          {
            painter.setPen(Qt::darkBlue);
            painter.drawText(pos.x() + 10, pos.y() + 10, i->getMetaValue(3).toQString());
          }
        }
      }
      ++num;
    }
  }

  void Spectrum2DCanvas::paintIcon_(const QPoint & pos, const QRgb & color, const String & icon, Size s, QPainter & p) const
  {
    p.save();
    p.setPen(color);
    p.setBrush(QBrush(QColor(color), Qt::SolidPattern));

    int s_half = (int)s / 2;

    if (icon == "diamond")
    {
      QPolygon pol;
      pol.putPoints(0, 4, pos.x() + s_half, pos.y(),
                    pos.x(), pos.y() + s_half,
                    pos.x() - (int)s_half, pos.y(),
                    pos.x(), pos.y() - (int)s_half);
      p.drawConvexPolygon(pol);
    }
    else if (icon == "square")
    {
      QPolygon pol;
      pol.putPoints(0, 4, pos.x() + s_half, pos.y() + s_half,
                    pos.x() - s_half, pos.y() + s_half,
                    pos.x() - s_half, pos.y() - s_half,
                    pos.x() + s_half, pos.y() - s_half);
      p.drawConvexPolygon(pol);
    }
    else if (icon == "circle")
    {
      p.drawEllipse(QRectF(pos.x() - s_half, pos.y() - s_half, s, s));
    }
    else if (icon == "triangle")
    {
      QPolygon pol;
      pol.putPoints(0, 3, pos.x(), pos.y() + s_half,
                    pos.x() + s_half, pos.y() - (int)s_half,
                    pos.x() - (int)s_half, pos.y() - (int)s_half);
      p.drawConvexPolygon(pol);
    }
    p.restore();
  }

  void Spectrum2DCanvas::paintTraceConvexHulls_(Size layer_index, QPainter & painter)
  {
    painter.setPen(Qt::black);

    const LayerData & layer = getLayer(layer_index);
    for (FeatureMapType::ConstIterator i = layer.getFeatureMap()->begin(); i != layer.getFeatureMap()->end(); ++i)
    {
      if (i->getRT() >= visible_area_.minPosition()[1] &&
          i->getRT() <= visible_area_.maxPosition()[1] &&
          i->getMZ() >= visible_area_.minPosition()[0] &&
          i->getMZ() <= visible_area_.maxPosition()[0] &&
          layer.filters.passes(*i)
          )
      {
        bool hasIdentifications = i->getPeptideIdentifications().size()>0
                               && i->getPeptideIdentifications()[0].getHits().size()>0;
        paintConvexHulls_(i->getConvexHulls(), hasIdentifications, painter);
      }
    }
  }

  void Spectrum2DCanvas::paintFeatureConvexHulls_(Size layer_index, QPainter & painter)
  {
    const LayerData & layer = getLayer(layer_index);
    for (FeatureMapType::ConstIterator i = layer.getFeatureMap()->begin(); i != layer.getFeatureMap()->end(); ++i)
    {
      if (i->getRT() >= visible_area_.minPosition()[1] &&
          i->getRT() <= visible_area_.maxPosition()[1] &&
          i->getMZ() >= visible_area_.minPosition()[0] &&
          i->getMZ() <= visible_area_.maxPosition()[0] &&
          layer.filters.passes(*i))
      {
        //paint hull points
        ConvexHull2D hull = i->getConvexHull();
        ConvexHull2D::PointArrayType ch_points = hull.getHullPoints();
        QPolygon points;
        points.resize((int)ch_points.size());

        UInt index = 0;
        QPoint pos;
        //iterate over hull points

        for (ConvexHull2D::PointArrayType::const_iterator it = ch_points.begin(); it != ch_points.end(); ++it, ++index)
        {
          dataToWidget_(it->getY(), it->getX(), pos);
          points.setPoint(index, pos);
        }
        //cout << "Hull: " << hull << " Points: " << points.size()<<endl;
        bool hasIdentifications = i->getPeptideIdentifications().size()>0
                               && i->getPeptideIdentifications()[0].getHits().size()>0;
        painter.setPen( hasIdentifications ? Qt::darkGreen : Qt::darkBlue );
        painter.drawPolygon(points);
      }
    }
  }

  void Spectrum2DCanvas::paintIdentifications_(Size layer_index, QPainter & painter)
  {
    const LayerData & layer = getLayer(layer_index);
    vector<PeptideIdentification>::const_iterator pep_begin, pep_end;
    if (layer.type == LayerData::DT_FEATURE)
    {
      pep_begin = layer.getFeatureMap()->getUnassignedPeptideIdentifications().begin();
      pep_end = layer.getFeatureMap()->getUnassignedPeptideIdentifications().end();
    }
    else if (layer.type == LayerData::DT_IDENT)
    {
      pep_begin = layer.peptides.begin();
      pep_end = layer.peptides.end();
    }
    else
      return;

    painter.setPen(Qt::darkRed);

    for (; pep_begin != pep_end; ++pep_begin)
    {
      if (!pep_begin->getHits().empty())
      {
        if (!pep_begin->hasRT() ||
            !pep_begin->hasMZ())
        {
          // TODO: show error message here
          continue;
        }
        double rt = pep_begin->getRT();
        if (rt < visible_area_.minPosition()[1] || rt > visible_area_.maxPosition()[1])
          continue;
        double mz = getIdentificationMZ_(layer_index, *pep_begin);
        if (mz < visible_area_.minPosition()[0] || mz > visible_area_.maxPosition()[0])
          continue;

        //draw dot
        QPoint pos;
        dataToWidget_(mz, rt, pos);
        painter.drawLine(pos.x(), pos.y() - 1.0, pos.x(), pos.y() + 1.0);
        painter.drawLine(pos.x() - 1.0, pos.y(), pos.x() + 1.0, pos.y());

        //draw sequence
        String sequence = pep_begin->getHits()[0].getSequence().toString();
        if (pep_begin->getHits().size() > 1)
          sequence += "...";
        painter.drawText(pos.x() + 10.0, pos.y() + 10.0, sequence.toQString());
      }
    }
  }

  void Spectrum2DCanvas::paintConvexHulls_(const vector<ConvexHull2D> & hulls, bool hasIdentifications, QPainter & painter)
  {
    QPolygon points;

    //iterate over all convex hulls
    for (Size hull = 0; hull < hulls.size(); ++hull)
    {
      ConvexHull2D::PointArrayType ch_points = hulls[hull].getHullPoints();
      points.resize((int)ch_points.size());
      UInt index = 0;
      QPoint pos;
      //iterate over hull points
      for (ConvexHull2D::PointArrayType::const_iterator it = ch_points.begin(); it != ch_points.end(); ++it, ++index)
      {
        dataToWidget_(it->getY(), it->getX(), pos);
        points.setPoint(index, pos);
      }
      painter.setPen(QPen(Qt::white, 5, Qt::DotLine, Qt::RoundCap, Qt::RoundJoin));
      painter.drawPolygon(points);
      painter.setPen(QPen( hasIdentifications ? Qt::green : Qt::blue, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
      painter.drawPolygon(points);
    }
  }

  void Spectrum2DCanvas::paintConsensusElements_(Size layer_index, QPainter & p)
  {
    const LayerData & layer = getLayer(layer_index);

    for (ConsensusMapType::ConstIterator i = layer.getConsensusMap()->begin(); i != layer.getConsensusMap()->end(); ++i)
    {
      paintConsensusElement_(layer_index, *i, p, true);
    }
  }

  void Spectrum2DCanvas::paintConsensusElement_(Size layer_index, const ConsensusFeature & cf, QPainter & p, bool use_buffer)
  {
    Int image_width = buffer_.width();
    Int image_height = buffer_.height();

    const LayerData & layer = getLayer(layer_index);

    if (isConsensusFeatureVisible_(cf, layer_index) && layer.filters.passes(cf))
    {
      //calculate position of consensus feature (centroid)
      QPoint consensus_pos;
      dataToWidget_(cf.getMZ(), cf.getRT(), consensus_pos);
      //iterate over elements
      for (ConsensusFeature::HandleSetType::const_iterator element = cf.begin(); element != cf.end(); ++element)
      {
        //calculate position of consensus element
        QPoint pos;
        dataToWidget_(element->getMZ(), element->getRT(), pos);
        //paint line
        p.drawLine(consensus_pos, pos);
        //paint point
        if (pos.x() > 0 && pos.y() > 0 && pos.x() < image_width - 1 && pos.y() < image_height - 1)
        {
          // use buffer only when not highlighting
          if (use_buffer)
          {
            buffer_.setPixel(pos.x(), pos.y(), Qt::black);
            buffer_.setPixel(pos.x() - 1, pos.y(), Qt::black);
            buffer_.setPixel(pos.x() + 1, pos.y(), Qt::black);
            buffer_.setPixel(pos.x(), pos.y() - 1, Qt::black);
            buffer_.setPixel(pos.x(), pos.y() + 1, Qt::black);
          }
          else
          {
            p.drawPoint(pos.x(), pos.y());
            p.drawPoint(pos.x() - 1, pos.y());
            p.drawPoint(pos.x() + 1, pos.y());
            p.drawPoint(pos.x(), pos.y() - 1);
            p.drawPoint(pos.x(), pos.y() + 1);
          }
        }
      }
    }

  }

  bool Spectrum2DCanvas::isConsensusFeatureVisible_(const ConsensusFeature & ce, Size layer_index)
  {
    // check the centroid first
    if (ce.getRT() >= visible_area_.minPosition()[1] &&
        ce.getRT() <= visible_area_.maxPosition()[1] &&
        ce.getMZ() >= visible_area_.minPosition()[0] &&
        ce.getMZ() <= visible_area_.maxPosition()[0])
    {
      return true;
    }

    // if element-flag is set, check if any of the consensus elements is visible
    if (getLayerFlag(layer_index, LayerData::C_ELEMENTS))
    {
      ConsensusFeature::HandleSetType::const_iterator element = ce.getFeatures().begin();
      for (; element != ce.getFeatures().end(); ++element)
      {
        if (element->getRT() >= visible_area_.minPosition()[1] &&
            element->getRT() <= visible_area_.maxPosition()[1] &&
            element->getMZ() >= visible_area_.minPosition()[0] &&
            element->getMZ() <= visible_area_.maxPosition()[0])
        {
          return true;
        }
      }
    }
    return false;
  }

  void Spectrum2DCanvas::intensityModeChange_()
  {
    String gradient_str;
    if (intensity_mode_ == IM_LOG)
    {
      gradient_str = MultiGradient::getDefaultGradientLogarithmicIntensityMode().toString();
    }
    else // linear
    {
      gradient_str = linear_gradient_.toString();
    }
    for (Size i = 0; i < layers_.size(); ++i)
    {
      layers_[i].param.setValue("dot:gradient", gradient_str);
      recalculateDotGradient_(i);
    }
    SpectrumCanvas::intensityModeChange_();
  }

  void Spectrum2DCanvas::recalculateDotGradient_(Size layer)
  {
    getLayer_(layer).gradient.fromString(getLayer_(layer).param.getValue("dot:gradient"));
    if (intensity_mode_ == IM_LOG)
    {
      getLayer_(layer).gradient.activatePrecalculationMode(0.0, std::log1p(overall_data_range_.maxPosition()[2]), param_.getValue("interpolation_steps"));
    }
    else
    {
      getLayer_(layer).gradient.activatePrecalculationMode(0.0, overall_data_range_.maxPosition()[2], param_.getValue("interpolation_steps"));
    }
  }

  void Spectrum2DCanvas::recalculateCurrentLayerDotGradient()
  {
    recalculateDotGradient_(current_layer_);
  }

  void Spectrum2DCanvas::updateProjections()
  {
    //find the last (visible) peak layers
    Size layer_count = 0;
    Size last_layer = 0;
    Size visible_layer_count = 0;
    Size visible_last_layer = 0;
    for (Size i = 0; i < getLayerCount(); ++i)
    {
      if (getLayer(i).type == LayerData::DT_PEAK)
      {
        layer_count++;
        last_layer = i;

        if (getLayer(i).visible)
        {
          visible_layer_count++;
          visible_last_layer = i;
        }
      }
      if (getLayer(i).type == LayerData::DT_CHROMATOGRAM)
      {
        //TODO CHROM
      }
    }

    //try to find the right layer to project
    const LayerData * layer = nullptr;
    //first choice: current layer
    if (layer_count != 0 && getCurrentLayer().type == LayerData::DT_PEAK)
    {
      layer = &(getCurrentLayer());
    }
    //second choice: the only peak layer
    else if (layer_count == 1)
    {
      layer = &(getLayer(last_layer));
    }
    //third choice: the only visible peak layer
    else if (visible_layer_count == 1)
    {
      layer = &(getLayer(visible_last_layer));
    }
    //no layer with peaks: disable projections
    else
    {
      emit toggleProjections();
      return;
    }

    //create projection data
    map<float, float> rt;
    map<int, float> mzint;
    map<int, int> mzcount;
    map<int, float> mzsum;

    UInt peak_count = 0;
    double intensity_max = 0.0;
    double intensity_sum = 0.0;

    // divide visible range into 100 bins (much faster than using a constant, e.g. 0.05, leading to many peaks for large maps without more information)
    float range = visible_area_.maxPosition()[0] - visible_area_.minPosition()[0];
    float mult = 100.0f / (range <= 0 ? 1 : range);

    for (ExperimentType::ConstAreaIterator i = layer->getPeakData()->areaBeginConst(visible_area_.minPosition()[1], visible_area_.maxPosition()[1], visible_area_.minPosition()[0], visible_area_.maxPosition()[0]);
         i != layer->getPeakData()->areaEndConst();
         ++i)
    {
      PeakIndex pi = i.getPeakIndex();
      if (layer->filters.passes((*layer->getPeakData())[pi.spectrum], pi.peak))
      {
        //sum
        ++peak_count;
        intensity_sum += i->getIntensity();
        mzint[int(i->getMZ() * mult)] += i->getIntensity();
        mzcount[int(i->getMZ() * mult)]++;
        mzsum[int(i->getMZ() * mult)] += i->getMZ();

        rt[i.getRT()] += i->getIntensity();
        //max
        intensity_max = max(intensity_max, (double)(i->getIntensity()));
      }
    }

    // write to spectra
    projection_mz_[0].resize(mzint.size() + 2);
    projection_mz_[0][0].setMZ(visible_area_.minPosition()[0]);
    projection_mz_[0][0].setIntensity(0.0);
    projection_mz_[0][1].setMZ(visible_area_.maxPosition()[0]);
    projection_mz_[0][1].setIntensity(0.0);
    projection_rt_[0].resize(rt.size() + 2);
    projection_rt_[0][0].setMZ(visible_area_.minPosition()[1]);
    projection_rt_[0][0].setIntensity(0.0);
    projection_rt_[0][1].setMZ(visible_area_.maxPosition()[1]);
    projection_rt_[0][1].setIntensity(0.0);

    Size i = 2;
    map<int, float>::iterator intit = mzint.begin();
    map<int, int>::iterator cit = mzcount.begin();

    for (map<int, float>::iterator it = mzsum.begin(); it != mzsum.end(); ++it)
    {
      projection_mz_[0][i].setMZ(it->second / cit->second);
      projection_mz_[0][i].setIntensity(intit->second);
      ++intit;
      ++cit;
      ++i;
    }

    i = 2;
    for (map<float, float>::iterator it = rt.begin(); it != rt.end(); ++it)
    {
      projection_rt_[0][i].setMZ(it->first);
      projection_rt_[0][i].setIntensity(it->second);
      ++i;
    }

    ExperimentSharedPtrType projection_mz_sptr(new ExperimentType(projection_mz_));
    ExperimentSharedPtrType projection_rt_sptr(new ExperimentType(projection_rt_));

    if (isMzToXAxis())
    {
      emit showProjectionHorizontal(projection_mz_sptr);
      emit showProjectionVertical(projection_rt_sptr);
    }
    else
    {
      emit showProjectionHorizontal(projection_rt_sptr);
      emit showProjectionVertical(projection_mz_sptr);
    }
    showProjectionInfo(peak_count, intensity_sum, intensity_max);
  }

  bool Spectrum2DCanvas::finishAdding_()
  {
    // unselect all peaks
    selected_peak_.clear();
    measurement_start_.clear();

    current_layer_ = getLayerCount() - 1;

    if (layers_.back().type == LayerData::DT_PEAK)   //peak data
    {
      update_buffer_ = true;
      //Abort if no data points are contained
      if ((currentPeakData_()->size() == 0 || currentPeakData_()->getSize() == 0) && currentPeakData_()->getDataRange().isEmpty())
      {
        layers_.resize(getLayerCount() - 1);
        if (current_layer_ != 0)
        {
          current_layer_ = current_layer_ - 1;
        }
        QMessageBox::critical(this, "Error", "Cannot add a dataset that contains no survey scans. Aborting!");
        return false;
      }
      if ((currentPeakData_()->getSize() == 0) && (!currentPeakData_()->getDataRange().isEmpty()))
      {
        setLayerFlag(LayerData::P_PRECURSORS, true); // show precursors if no MS1 data is contained
      }
    }
    else if (layers_.back().type == LayerData::DT_FEATURE)  //feature data
    {
      getCurrentLayer_().getFeatureMap()->updateRanges();
      setLayerFlag(LayerData::F_HULL, true);

      //Abort if no data points are contained
      if (getCurrentLayer_().getFeatureMap()->size() == 0)
      {
        layers_.resize(getLayerCount() - 1);
        if (current_layer_ != 0)
        {
          current_layer_ = current_layer_ - 1;
        }
        QMessageBox::critical(this, "Error", "Cannot add an empty dataset. Aborting!");
        return false;
      }
    }
    else if (layers_.back().type == LayerData::DT_CONSENSUS)  //consensus feature data
    {
      getCurrentLayer_().getConsensusMap()->updateRanges();

      // abort if no data points are contained
      if (getCurrentLayer_().getConsensusMap()->size() == 0)
      {
        layers_.resize(getLayerCount() - 1);
        if (current_layer_ != 0)
          current_layer_ = current_layer_ - 1;
        QMessageBox::critical(this, "Error", "Cannot add an empty dataset. Aborting!");
        return false;
      }
    }
    else if (layers_.back().type == LayerData::DT_CHROMATOGRAM)  //chromatogram data
    {

      //TODO CHROM
      currentPeakData_()->sortChromatograms(true);
      currentPeakData_()->updateRanges(1);

      update_buffer_ = true;

      // abort if no data points are contained
      if (currentPeakData_()->getChromatograms().empty())
      {
        layers_.resize(getLayerCount() - 1);
        if (current_layer_ != 0)
          current_layer_ = current_layer_ - 1;
        QMessageBox::critical(this, "Error", "Cannot add a dataset that contains no chromatograms. Aborting!");
        return false;
      }
    }
    else if (layers_.back().type == LayerData::DT_IDENT)   // identification data
    {
      // abort if no data points are contained
      if (getCurrentLayer_().peptides.empty())
      {
        layers_.resize(getLayerCount() - 1);
        if (current_layer_ != 0)
          current_layer_ = current_layer_ - 1;
        QMessageBox::critical(this, "Error", "Cannot add an empty dataset. Aborting!");
        return false;
      }
    }

    // warn if negative intensities are contained
    if (getMinIntensity(current_layer_) < 0.0)
    {
      QMessageBox::warning(this, "Warning", "This dataset contains negative intensities. Use it at your own risk!");
    }

    // overall values update
    recalculateRanges_(0, 1, 2);
    if (layers_.size() == 1)
    {
      resetZoom(false); //no repaint as this is done in intensityModeChange_() anyway
    }

    if (getLayerCount() == 2)
    {
      setIntensityMode(IM_PERCENTAGE);
    }
    intensityModeChange_();

    emit layerActivated(this);

    return true;
  }

  void Spectrum2DCanvas::removeLayer(Size layer_index)
  {
    if (layer_index >= getLayerCount())
    {
      return;
    }

    // remove the data
    layers_.erase(layers_.begin() + layer_index);

    // update visible area and boundaries
    DRange<3> old_data_range = overall_data_range_;
    recalculateRanges_(0, 1, 2);

    // only reset zoom if data range has been changed
    if (old_data_range != overall_data_range_)
    {
      resetZoom(false); // no repaint as this is done in intensityModeChange_() anyway
    }

    // update current layer if it became invalid
    if (current_layer_ != 0 && current_layer_ >= getLayerCount())
    {
      current_layer_ = getLayerCount() - 1;
    }

    if (layers_.empty())
    {
      overall_data_range_ = DRange<3>::empty;
      update_buffer_ = true;
      update_(OPENMS_PRETTY_FUNCTION);
      return;
    }

    // unselect all peaks
    selected_peak_.clear();
    measurement_start_.clear();

    intensityModeChange_();

    emit layerActivated(this);
  }

  // change the current layer
  void Spectrum2DCanvas::activateLayer(Size layer_index)
  {
    if (layer_index >= getLayerCount() || layer_index == current_layer_)
    {
      return;
    }

    // unselect all peaks
    selected_peak_.clear();
    measurement_start_.clear();

    current_layer_ = layer_index;
    emit layerActivated(this);

    update_(OPENMS_PRETTY_FUNCTION);
  }

  void Spectrum2DCanvas::recalculateSnapFactor_()
  {
    snap_factors_ = vector<double>(getLayerCount(), 1.0);

    if (intensity_mode_ == IM_SNAP)
    {
      for (Size i = 0; i < getLayerCount(); i++)
      {
        if (getLayer(i).visible)
        {
          double local_max  = -numeric_limits<double>::max();
          if (getLayer(i).type == LayerData::DT_PEAK)
          {
            for (ExperimentType::ConstAreaIterator it = getLayer(i).getPeakData()->areaBeginConst(visible_area_.minPosition()[1], visible_area_.maxPosition()[1], visible_area_.minPosition()[0], visible_area_.maxPosition()[0]);
                 it != getLayer(i).getPeakData()->areaEndConst();
                 ++it)
            {
              PeakIndex pi = it.getPeakIndex();
              if (it->getIntensity() > local_max && getLayer(i).filters.passes((*getLayer(i).getPeakData())[pi.spectrum], pi.peak))
              {
                local_max = it->getIntensity();
              }
            }
          }
          else if (getLayer(i).type == LayerData::DT_FEATURE)         // features
          {
            for (FeatureMapType::ConstIterator it = getLayer(i).getFeatureMap()->begin();
                 it != getLayer(i).getFeatureMap()->end();
                 ++it)
            {
              if (it->getRT() >= visible_area_.minPosition()[1] &&
                  it->getRT() <= visible_area_.maxPosition()[1] &&
                  it->getMZ() >= visible_area_.minPosition()[0] &&
                  it->getMZ() <= visible_area_.maxPosition()[0] &&
                  getLayer(i).filters.passes(*it) &&
                  it->getIntensity() > local_max)
              {
                local_max = it->getIntensity();
              }
            }
          }
          else if (getLayer(i).type == LayerData::DT_CONSENSUS)         // consensus
          {
            for (ConsensusMapType::ConstIterator it = getLayer(i).getConsensusMap()->begin();
                 it != getLayer(i).getConsensusMap()->end();
                 ++it)
            {
              if (it->getRT() >= visible_area_.minPosition()[1] &&
                  it->getRT() <= visible_area_.maxPosition()[1] &&
                  it->getMZ() >= visible_area_.minPosition()[0] &&
                  it->getMZ() <= visible_area_.maxPosition()[0] &&
                  getLayer(i).filters.passes(*it) &&
                  it->getIntensity() > local_max)
              {
                local_max = it->getIntensity();
              }
            }
          }
          else if (getLayer(i).type == LayerData::DT_CHROMATOGRAM)         // chromatogr.
          {
            //TODO CHROM
          }
          else if (getLayer(i).type == LayerData::DT_IDENT)         // identifications
          {
            //TODO IDENT
          }

          if (local_max > 0.0)
          {
            snap_factors_[i] = overall_data_range_.maxPosition()[2] / local_max;
          }
        }
      }
    }
  }

  void Spectrum2DCanvas::updateScrollbars_()
  {
    if (isMzToXAxis())
    {
      emit updateHScrollbar(overall_data_range_.minPosition()[0], visible_area_.minPosition()[0], visible_area_.maxPosition()[0], overall_data_range_.maxPosition()[0]);
      emit updateVScrollbar(overall_data_range_.minPosition()[1], visible_area_.minPosition()[1], visible_area_.maxPosition()[1], overall_data_range_.maxPosition()[1]);
    }
    else
    {
      emit updateVScrollbar(overall_data_range_.minPosition()[0], visible_area_.minPosition()[0], visible_area_.maxPosition()[0], overall_data_range_.maxPosition()[0]);
      emit updateHScrollbar(overall_data_range_.minPosition()[1], visible_area_.minPosition()[1], visible_area_.maxPosition()[1], overall_data_range_.maxPosition()[1]);
    }
  }

  void Spectrum2DCanvas::horizontalScrollBarChange(int value)
  {
    AreaType new_area = visible_area_;
    if (!isMzToXAxis())
    {
      new_area.setMinY(value);
      new_area.setMaxY(value + (visible_area_.height()));
    }
    else
    {
      new_area.setMinX(value);
      new_area.setMaxX(value + (visible_area_.width()));
    }
    //cout << OPENMS_PRETTY_FUNCTION << endl;
    changeVisibleArea_(new_area);
    emit layerZoomChanged(this);
  }

  void Spectrum2DCanvas::verticalScrollBarChange(int value)
  {
    // invert 'value' (since the VERTICAL(!) scrollbar's range is negative -- see SpectrumWidget::updateVScrollbar())
    // this is independent on isMzToXAxis()!
    value *= -1;
    
    AreaType new_area = visible_area_;
    if (!isMzToXAxis())
    {
      new_area.setMinX(value);
      new_area.setMaxX(value + visible_area_.width());
    }
    else
    {
      new_area.setMinY(value);
      new_area.setMaxY(value + (visible_area_.height()));
    }
    //cout << OPENMS_PRETTY_FUNCTION << endl;
    changeVisibleArea_(new_area);
    emit layerZoomChanged(this);
  }

  void Spectrum2DCanvas::paintEvent(QPaintEvent * e)
  {
    //Only fill background if no layer is present
    if (getLayerCount() == 0)
    {
      QPainter painter;
      painter.begin(this);
      painter.fillRect(0, 0, this->width(), this->height(), QColor(param_.getValue("background_color").toQString()));
      painter.end();
      e->accept();
      return;
    }

#ifdef DEBUG_TOPPVIEW
    cout << "BEGIN " << OPENMS_PRETTY_FUNCTION << endl;
    cout << "  Visible area -- m/z: " << visible_area_.minX() << " - " << visible_area_.maxX() << " rt: " << visible_area_.minY() << " - " << visible_area_.maxY() << endl;
    cout << "  Overall area -- m/z: " << overall_data_range_.minPosition()[0] << " - " << overall_data_range_.maxPosition()[0] << " rt: " << overall_data_range_.minPosition()[1] << " - " << overall_data_range_.maxPosition()[1] << endl;
#endif

    //timing
    QTime overall_timer;
    if (show_timing_)
    {
      overall_timer.start();
      if (update_buffer_)
      {
        cout << "Updating buffer:" << endl;
      }
      else
      {
        cout << "Copying buffer:" << endl;
      }
    }

    QPainter painter;
    if (update_buffer_)
    {
      update_buffer_ = false;

      // recalculate snap factor
      recalculateSnapFactor_();

      buffer_.fill(QColor(param_.getValue("background_color").toQString()).rgb());
      painter.begin(&buffer_);
      QTime layer_timer;

      for (Size i = 0; i < getLayerCount(); i++)
      {
        // timing
        if (show_timing_)
        {
          layer_timer.start();
        }

        if (getLayer(i).visible)
        {
          if (getLayer(i).type == LayerData::DT_PEAK)
          {
           // nothing .. currently
          }
          else if (getLayer(i).type == LayerData::DT_FEATURE)
          {
            //cout << "dot feature layer: " << i << endl;
            if (getLayerFlag(i, LayerData::F_HULLS))
            {
              paintTraceConvexHulls_(i, painter);
            }
            if (getLayerFlag(i, LayerData::F_HULL))
            {
              paintFeatureConvexHulls_(i, painter);
            }
            if (getLayerFlag(i, LayerData::F_UNASSIGNED))
            {
              paintIdentifications_(i, painter);
            }
          }
          else if (getLayer(i).type == LayerData::DT_CONSENSUS)
          {
            if (getLayerFlag(i, LayerData::C_ELEMENTS))
            {
              paintConsensusElements_(i, painter);
            }
          }
          else if (getLayer(i).type == LayerData::DT_CHROMATOGRAM)
          {
          }
          else if (getLayer(i).type == LayerData::DT_IDENT)
          {
          }
          // all layers
          paintDots_(i, painter);
        }
        //timing
        if (show_timing_)
        {
          cout << "  -layer " << i << " time: " << layer_timer.elapsed() << " ms" << endl;
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
      painter.drawImage(rects[i].topLeft(), buffer_, rects[i]);
    }

    //draw measurement peak
    if (action_mode_ == AM_MEASURE && measurement_start_.isValid())
    {
      painter.setPen(Qt::black);

      QPoint line_begin, line_end;
      // start of line
      if (selected_peak_.isValid())
      {
        if (getCurrentLayer().type == LayerData::DT_FEATURE)
        {
          dataToWidget_(selected_peak_.getFeature(*getCurrentLayer().getFeatureMap()).getMZ(), selected_peak_.getFeature(*getCurrentLayer().getFeatureMap()).getRT(), line_begin);
        }
        else if (getCurrentLayer().type == LayerData::DT_PEAK)
        {
          dataToWidget_(selected_peak_.getPeak(*getCurrentLayer().getPeakData()).getMZ(), selected_peak_.getSpectrum(*getCurrentLayer().getPeakData()).getRT(), line_begin);
        }
        else if (getCurrentLayer().type == LayerData::DT_CONSENSUS)
        {
          dataToWidget_(selected_peak_.getFeature(*getCurrentLayer().getConsensusMap()).getMZ(), selected_peak_.getFeature(*getCurrentLayer().getConsensusMap()).getRT(), line_begin);
        }
        else if (getCurrentLayer().type == LayerData::DT_CHROMATOGRAM)
        {
          //TODO CHROM
        }
      }
      else
      {
        line_begin = last_mouse_pos_;
      }

      // end of line
      if (getCurrentLayer().type == LayerData::DT_FEATURE)
      {
        dataToWidget_(measurement_start_.getFeature(*getCurrentLayer().getFeatureMap()).getMZ(), measurement_start_.getFeature(*getCurrentLayer().getFeatureMap()).getRT(), line_end);
      }
      else if (getCurrentLayer().type == LayerData::DT_PEAK)
      {
        dataToWidget_(measurement_start_.getPeak(*getCurrentLayer().getPeakData()).getMZ(), measurement_start_.getSpectrum(*getCurrentLayer().getPeakData()).getRT(), line_end);
      }
      else if (getCurrentLayer().type == LayerData::DT_CONSENSUS)
      {
        dataToWidget_(measurement_start_.getFeature(*getCurrentLayer().getConsensusMap()).getMZ(), measurement_start_.getFeature(*getCurrentLayer().getConsensusMap()).getRT(), line_end);
      }
      else if (getCurrentLayer().type == LayerData::DT_CHROMATOGRAM)
      {
        //TODO CHROM
      }
      painter.drawLine(line_begin, line_end);

      highlightPeak_(painter, measurement_start_);
    }

    // draw convex hulls or consensus feature elements
    if (selected_peak_.isValid())
    {
      if (getCurrentLayer().type == LayerData::DT_FEATURE)
      {
        painter.setPen(QPen(Qt::red, 2));
        const Feature& f = selected_peak_.getFeature(*getCurrentLayer().getFeatureMap());
        paintConvexHulls_(f.getConvexHulls(),
                          f.getPeptideIdentifications().size() && f.getPeptideIdentifications()[0].getHits().size(),
                          painter);
      }
      else if (getCurrentLayer().type == LayerData::DT_CONSENSUS && getLayerFlag(current_layer_, LayerData::C_ELEMENTS))
      {
        painter.setPen(QPen(Qt::red, 2));
        paintConsensusElement_(current_layer_, selected_peak_.getFeature(*getCurrentLayer().getConsensusMap()), painter, false);
      }
    }

    if (action_mode_ == AM_MEASURE || action_mode_ == AM_TRANSLATE)
    {
      highlightPeak_(painter, selected_peak_);
    }

    // draw delta for measuring
    if (action_mode_ == AM_MEASURE && measurement_start_.isValid())
    {
      drawDeltas_(painter, measurement_start_, selected_peak_);
    }
    else
    {
      drawCoordinates_(painter, selected_peak_);
    }

    painter.end();
#ifdef DEBUG_TOPPVIEW
    cout << "END   " << OPENMS_PRETTY_FUNCTION << endl;
#endif
    if (show_timing_)
    {
      cout << "  -overall time: " << overall_timer.elapsed() << " ms" << endl << endl;
    }
  }

  void Spectrum2DCanvas::drawCoordinates_(QPainter & painter, const PeakIndex & peak)
  {
    if (!peak.isValid())
      return;

    //determine coordinates;
    double mz = 0.0;
    double rt = 0.0;
    float it = 0.0;
    Int charge = 0;
    double quality = 0.0;
    Size size = 0;
    const Feature* f = nullptr;
    const ConsensusFeature* cf = nullptr;
    ConsensusFeature::HandleSetType sub_features;

    switch (getCurrentLayer().type)
    {
    case LayerData::DT_FEATURE:
    {
      f = &peak.getFeature(*getCurrentLayer().getFeatureMap());
      mz = f->getMZ();
      rt = f->getRT();
      it = f->getIntensity();
      charge  = f->getCharge();
      quality = f->getOverallQuality();
    }
    break;

    case LayerData::DT_PEAK:
    {
      const Peak1D & p = peak.getPeak(*getCurrentLayer().getPeakData());
      const MSSpectrum & s = peak.getSpectrum(*getCurrentLayer().getPeakData());
      mz = p.getMZ();
      rt = s.getRT();
      it = p.getIntensity();
    }
    break;

    case LayerData::DT_CONSENSUS:
    {
      cf = &peak.getFeature(*getCurrentLayer().getConsensusMap());

      mz = cf->getMZ();
      rt = cf->getRT();
      it = cf->getIntensity();
      charge  = cf->getCharge();
      quality = cf->getQuality();
      sub_features = cf->getFeatures();
      size =  sub_features.size();
    }
    break;

    case LayerData::DT_CHROMATOGRAM:
    {
      const LayerData & layer = getCurrentLayer();

      PeakMap exp;
      exp = *layer.getPeakData();

      vector<MSChromatogram >::const_iterator iter = exp.getChromatograms().begin();
      iter += peak.spectrum;

      mz = iter->getPrecursor().getMZ();
      rt = iter->front().getRT();
    }
    break;

    case LayerData::DT_IDENT:
      // TODO implement
      break;

    default:
      break;
    }

    // draw text
    QStringList lines;
    lines.push_back("RT:  " + QLocale::c().toString(rt, 'f', 2)); // adds group separators (consistency with intensity)
    lines.push_back("m/z: " + QLocale::c().toString(mz, 'f', 5)); // adds group separators (consistency with intensity)
    lines.push_back("Int: " + QLocale::c().toString(it, 'f', 2)); // adds group separators (every 1e3), to better visualize large numbers (e.g. 23.009.646.54,3)

    if (getCurrentLayer().type == LayerData::DT_FEATURE || getCurrentLayer().type == LayerData::DT_CONSENSUS)
    {
      lines.push_back("Charge: " + QString::number(charge));
      lines.push_back("Quality: " + QString::number(quality, 'f', 4));
      // peptide identifications
      const PeptideIdentification* pis = nullptr;
      if ( f && f->getPeptideIdentifications().size() > 0 ) {
        pis = &f->getPeptideIdentifications()[0];
      }
      else if ( cf && cf->getPeptideIdentifications().size() > 0 ) {
        pis = &cf->getPeptideIdentifications()[0];
      }
      if ( pis && pis->getHits().size() ) {
          Size nHits = pis->getHits().size();
          for (Size j = 0; j < nHits; ++j)
          {
            lines.push_back( "Peptide" + ( nHits > 1 ? "[" + QString::number(j+1) + "]" : "" ) + ": "
                             + pis->getHits()[j].getSequence().toString().toQString() );
          }
      }
    }

    if (getCurrentLayer().type == LayerData::DT_CONSENSUS)
    {
      lines.push_back("Size: " + QString::number(size));
      for (ConsensusFeature::HandleSetType::const_iterator it = sub_features.begin(); it != sub_features.end(); ++it)
      {
        lines.push_back("Feature m/z: " + QLocale::c().toString(it->getMZ(), 'f', 5) + // adds group separators (consistency with intensity)
                        "  rt: " + QLocale::c().toString(it->getRT(), 'f', 2) +        // adds group separators (consistency with intensity)
                        "   q: " + QString::number(it->getCharge(), 'f', 2) +
                        "  intensity: " + QLocale::c().toString(it->getIntensity(), 'f', 2)); // adds group separators (every 1e3), to better visualize large numbers (e.g. 23.009.646.54,3)
      }
    }

    drawText_(painter, lines);
  }

  void Spectrum2DCanvas::drawDeltas_(QPainter & painter, const PeakIndex & start, const PeakIndex & end)
  {
    if (!start.isValid())
      return;

    //determine coordinates;
    double mz = 0.0;
    double rt = 0.0;
    float it = 0.0;
    float ppm = 0.0;

    if (getCurrentLayer().type == LayerData::DT_FEATURE)
    {
      if (end.isValid())
      {
        mz = end.getFeature(*getCurrentLayer().getFeatureMap()).getMZ() - start.getFeature(*getCurrentLayer().getFeatureMap()).getMZ();
        rt = end.getFeature(*getCurrentLayer().getFeatureMap()).getRT() - start.getFeature(*getCurrentLayer().getFeatureMap()).getRT();
        it = end.getFeature(*getCurrentLayer().getFeatureMap()).getIntensity() / start.getFeature(*getCurrentLayer().getFeatureMap()).getIntensity();
      }
      else
      {
        PointType point = widgetToData_(last_mouse_pos_);
        mz = point[0] - start.getFeature(*getCurrentLayer().getFeatureMap()).getMZ();
        rt = point[1] - start.getFeature(*getCurrentLayer().getFeatureMap()).getRT();
        it = std::numeric_limits<double>::quiet_NaN();
      }
      ppm = (mz / start.getFeature(*getCurrentLayer().getFeatureMap()).getMZ()) * 1e6;
    }
    else if (getCurrentLayer().type == LayerData::DT_PEAK)
    {
      if (end.isValid())
      {
        mz = end.getPeak(*getCurrentLayer().getPeakData()).getMZ() - start.getPeak(*getCurrentLayer().getPeakData()).getMZ();
        rt = end.getSpectrum(*getCurrentLayer().getPeakData()).getRT() - start.getSpectrum(*getCurrentLayer().getPeakData()).getRT();
        it = end.getPeak(*getCurrentLayer().getPeakData()).getIntensity() / start.getPeak(*getCurrentLayer().getPeakData()).getIntensity();
      }
      else
      {
        PointType point = widgetToData_(last_mouse_pos_);
        mz = point[0] - start.getPeak(*getCurrentLayer().getPeakData()).getMZ();
        rt = point[1] - start.getSpectrum(*getCurrentLayer().getPeakData()).getRT();
        it = std::numeric_limits<double>::quiet_NaN();
      }
      ppm = (mz / start.getPeak(*getCurrentLayer().getPeakData()).getMZ()) * 1e6;
    }
    else if (getCurrentLayer().type == LayerData::DT_CONSENSUS)
    {
      if (end.isValid())
      {
        mz = end.getFeature(*getCurrentLayer().getConsensusMap()).getMZ() - start.getFeature(*getCurrentLayer().getConsensusMap()).getMZ();
        rt = end.getFeature(*getCurrentLayer().getConsensusMap()).getRT() - start.getFeature(*getCurrentLayer().getConsensusMap()).getRT();
        it = end.getFeature(*getCurrentLayer().getConsensusMap()).getIntensity() / start.getFeature(*getCurrentLayer().getConsensusMap()).getIntensity();
      }
      else
      {
        PointType point = widgetToData_(last_mouse_pos_);
        mz = point[0] - start.getFeature(*getCurrentLayer().getConsensusMap()).getMZ();
        rt = point[1] - start.getFeature(*getCurrentLayer().getConsensusMap()).getRT();
        it = std::numeric_limits<double>::quiet_NaN();
      }
      ppm = (mz / start.getFeature(*getCurrentLayer().getConsensusMap()).getMZ()) * 1e6;
    }
    else if (getCurrentLayer().type == LayerData::DT_CHROMATOGRAM)
    {
      //TODO CHROM
    }
    else if (getCurrentLayer().type == LayerData::DT_IDENT)
    {
      // TODO IDENT
    }

    //draw text
    QStringList lines;
    lines.push_back("RT delta:  " + QString::number(rt, 'f', 2));
    lines.push_back("m/z delta: " + QString::number(mz, 'f', 6) + " (" + QString::number(ppm, 'f', 1) +" ppm)");
    if (boost::math::isinf(it) || boost::math::isnan(it))
    {
      lines.push_back("Int ratio: n/a");
    }
    else
    {
      lines.push_back("Int ratio: " + QString::number(it, 'f', 2));
    }
    drawText_(painter, lines);
  }

  void Spectrum2DCanvas::mousePressEvent(QMouseEvent * e)
  {
    last_mouse_pos_ = e->pos();

    if (e->button() == Qt::LeftButton)
    {
      if (action_mode_ == AM_MEASURE)
      {
        if (selected_peak_.isValid())
        {
          measurement_start_ = selected_peak_;
        }
        else
        {
          measurement_start_.clear();
        }
      }
      else if (action_mode_ == AM_ZOOM)
      {
        //translate (if not moving features)
        if ( !(getCurrentLayer().type == LayerData::DT_FEATURE) || !selected_peak_.isValid())
        {
          rubber_band_.setGeometry(QRect(e->pos(), QSize()));
          rubber_band_.show();
        }
      }
    }
  }

  void Spectrum2DCanvas::mouseMoveEvent(QMouseEvent * e)
  {
    grabKeyboard();     // (re-)grab keyboard after it has been released by unhandled key
    QPoint pos = e->pos();
    PointType data_pos = widgetToData_(pos);
    emit sendCursorStatus(data_pos[0], data_pos[1]);

    PeakIndex near_peak = findNearestPeak_(pos);

    //highlight current peak and display peak coordinates
    if (action_mode_ == AM_MEASURE || (action_mode_ == AM_TRANSLATE && !(e->buttons() & Qt::LeftButton)))
    {
      //highlight peak
      selected_peak_ = near_peak;
      update_(OPENMS_PRETTY_FUNCTION);

      //show meta data in status bar (if available)
      if (selected_peak_.isValid())
      {
        String status;
        if (getCurrentLayer().type == LayerData::DT_FEATURE || getCurrentLayer().type == LayerData::DT_CONSENSUS)
        {
          //add meta info
          const BaseFeature* f;
          if (getCurrentLayer().type == LayerData::DT_FEATURE)
          {
            f = &selected_peak_.getFeature(*getCurrentLayer().getFeatureMap());
          }
          else
          {
            f = &selected_peak_.getFeature(*getCurrentLayer().getConsensusMap());
          }
          std::vector<String> keys;
          f->getKeys(keys);
          for (Size m = 0; m < keys.size(); ++m)
          {
            status += " " + keys[m] + ": ";
            DataValue dv = f->getMetaValue(keys[m]);
            if (dv.valueType() == DataValue::DOUBLE_VALUE)
            { // use less precision for large numbers, for better readability
              int precision(2);
              if ((double)dv < 10) precision = 5;
              // ... and add 1k separators, e.g. '540,321.99'
              status += QLocale::c().toString((double)dv, 'f', precision);
            }
            else
            {
              status += (String)dv;
            }
          }
        }
        else if (getCurrentLayer().type == LayerData::DT_PEAK)
        {
          //meta info
          const ExperimentType::SpectrumType & s = selected_peak_.getSpectrum(*getCurrentLayer().getPeakData());
          for (Size m = 0; m < s.getFloatDataArrays().size(); ++m)
          {
            if (selected_peak_.peak < s.getFloatDataArrays()[m].size())
            {
              status += s.getFloatDataArrays()[m].getName() + ": " + s.getFloatDataArrays()[m][selected_peak_.peak] + " ";
            }
          }
          for (Size m = 0; m < s.getIntegerDataArrays().size(); ++m)
          {
            if (selected_peak_.peak < s.getIntegerDataArrays()[m].size())
            {
              status += s.getIntegerDataArrays()[m].getName() + ": " + s.getIntegerDataArrays()[m][selected_peak_.peak] + " ";
            }
          }
          for (Size m = 0; m < s.getStringDataArrays().size(); ++m)
          {
            if (selected_peak_.peak < s.getStringDataArrays()[m].size())
            {
              status += s.getStringDataArrays()[m].getName() + ": " + s.getStringDataArrays()[m][selected_peak_.peak] + " ";
            }
          }
        }
        else if (getCurrentLayer().type == LayerData::DT_CHROMATOGRAM) // chromatogram
        {
          //TODO CHROM
        }
        if (status != "")
        {
          emit sendStatusMessage(status, 0);
        }
      }
    }
    else if (action_mode_ == AM_ZOOM)
    {
      //Zoom mode => no peak should be selected
      selected_peak_.clear();
      update_(OPENMS_PRETTY_FUNCTION);
    }

    if (action_mode_ == AM_MEASURE)
    {
      last_mouse_pos_ = pos;
    }
    else if (action_mode_ == AM_ZOOM)
    {
      //if mouse button is held down, enlarge the selection
      if (e->buttons() & Qt::LeftButton)
      {
        rubber_band_.setGeometry(QRect(last_mouse_pos_, pos).normalized());
        rubber_band_.show();         //if the mouse button is pressed before the zoom key is pressed

        update_(OPENMS_PRETTY_FUNCTION);
      }
    }
    else if (action_mode_ == AM_TRANSLATE)
    {
      if (e->buttons() & Qt::LeftButton)
      {
        if (getCurrentLayer().modifiable && getCurrentLayer().type == LayerData::DT_FEATURE && selected_peak_.isValid())       //move feature
        {
          PointType new_data = widgetToData_(pos);
          double mz = new_data[0];
          double rt = new_data[1];

          //restrict the movement to the data range
          mz = max(mz, overall_data_range_.minPosition()[0]);
          mz = min(mz, overall_data_range_.maxPosition()[0]);
          rt = max(rt, overall_data_range_.minPosition()[1]);
          rt = min(rt, overall_data_range_.maxPosition()[1]);

          (*getCurrentLayer_().getFeatureMap())[selected_peak_.peak].setRT(rt);
          (*getCurrentLayer_().getFeatureMap())[selected_peak_.peak].setMZ(mz);

          update_buffer_ = true;
          update_(OPENMS_PRETTY_FUNCTION);
          modificationStatus_(activeLayerIndex(), true);
        }
        else         //translate
        {
          //calculate data coordinates of shift
          PointType old_data = widgetToData_(last_mouse_pos_);
          PointType new_data = widgetToData_(pos);
          //calculate x shift
          double shift = old_data.getX() - new_data.getX();
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
          shift = old_data.getY() - new_data.getY();
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
          //cout << "New area: x " << newLoX <<"-"<< newHiX << " - y "<<newLoY <<"-"<< newHiY << endl;
          //cout << OPENMS_PRETTY_FUNCTION << endl;
          changeVisibleArea_(AreaType(newLoX, newLoY, newHiX, newHiY));
          emit layerZoomChanged(this);
          last_mouse_pos_ = pos;
        }
      }
    }
  }

  void Spectrum2DCanvas::mouseReleaseEvent(QMouseEvent * e)
  {
    if (e->button() == Qt::LeftButton)
    {
      if (action_mode_ == AM_MEASURE)
      {
        if (!selected_peak_.isValid())
        {
          measurement_start_.clear();
        }
        measurement_start_.clear();
        update_(OPENMS_PRETTY_FUNCTION);
      }
      else if (action_mode_ == AM_ZOOM)
      {
        rubber_band_.hide();
        QRect rect = rubber_band_.geometry();
        if (rect.width() != 0 && rect.height() != 0)
        {
          AreaType area(widgetToData_(rect.topLeft()), widgetToData_(rect.bottomRight()));
          //cout << OPENMS_PRETTY_FUNCTION << endl;
          changeVisibleArea_(area, true, true);
          emit layerZoomChanged(this);
        }
      }
    }
  }

  void Spectrum2DCanvas::contextMenuEvent(QContextMenuEvent * e)
  {
    //Abort if there are no layers
    if (layers_.empty())
    {
      return;
    }

    double mz = widgetToData_(e->pos())[0];
    double rt = widgetToData_(e->pos())[1];

    const LayerData & layer = getCurrentLayer();

    QMenu * context_menu = new QMenu(this);

    QAction * a = nullptr;
    QAction * result = nullptr;

    //Display name and warn if current layer invisible
    String layer_name = String("Layer: ") + layer.name;
    if (!layer.visible)
    {
      layer_name += " (invisible)";
    }
    context_menu->addAction(layer_name.toQString())->setEnabled(false);
    context_menu->addSeparator();

    context_menu->addAction("Layer meta data");

    QMenu * settings_menu = new QMenu("Settings");
    settings_menu->addAction("Show/hide grid lines");
    settings_menu->addAction("Show/hide axis legends");
    context_menu->addSeparator();

    context_menu->addAction("Switch to 3D view");

    //-------------------PEAKS----------------------------------
    if (layer.type == LayerData::DT_PEAK)
    {
      //add settings
      settings_menu->addSeparator();
      settings_menu->addAction("Show/hide projections");
      settings_menu->addAction("Show/hide MS/MS precursors");

      //add surrounding survey scans
      //find nearest survey scan
      SignedSize size = getCurrentLayer().getPeakData()->size();
      Int current = getCurrentLayer().getPeakData()->RTBegin(rt) - getCurrentLayer().getPeakData()->begin();
      if (current == size)  // if only one element is present RTBegin points to one after the last element (see RTBegin implementation)
      {
        current = 0;
      }

      SignedSize i = 0;
      while (current + i < size || current - i >= 0)
      {
        if (current + i < size && (*getCurrentLayer().getPeakData())[current + i].getMSLevel() == 1)
        {
          current = current + i;
          break;
        }
        if (current - i >= 0 && (*getCurrentLayer().getPeakData())[current - i].getMSLevel() == 1)
        {
          current = current - i;
          break;
        }
        ++i;
      }
      //search for four scans in both directions
      vector<Int> indices;
      indices.push_back(current);
      i = 1;
      while (current - i >= 0 && indices.size() < 5)
      {
        if ((*getCurrentLayer().getPeakData())[current - i].getMSLevel() == 1)
        {
          indices.push_back(current - i);
        }
        ++i;
      }
      i = 1;
      while (current + i < size && indices.size() < 9)
      {
        if ((*getCurrentLayer().getPeakData())[current + i].getMSLevel() == 1)
        {
          indices.push_back(current + i);
        }
        ++i;
      }
      sort(indices.rbegin(), indices.rend());
      QMenu * ms1_scans = context_menu->addMenu("Survey scan in 1D");
      QMenu * ms1_meta = context_menu->addMenu("Survey scan meta data");
      context_menu->addSeparator();
      for (i = 0; i < (Int)indices.size(); ++i)
      {
        if (indices[i] == current)
        {
          ms1_scans->addSeparator();
        }
        a = ms1_scans->addAction(QString("RT: ") + QString::number((*getCurrentLayer().getPeakData())[indices[i]].getRT()));
        a->setData(indices[i]);
        if (indices[i] == current)
        {
          ms1_scans->addSeparator();
        }

        if (indices[i] == current)
        {
          ms1_meta->addSeparator();
        }
        a = ms1_meta->addAction(QString("RT: ") + QString::number((*getCurrentLayer().getPeakData())[indices[i]].getRT()));
        a->setData(indices[i]);
        if (indices[i] == current)
        {
          ms1_meta->addSeparator();
        }
      }

      // add surrounding fragment scans
      // - We first attempt to look at the position where the user clicked
      // - Next we look within the +/- 5 scans around that position
      // - Next we look within the whole visible area
      QMenu * msn_scans = new QMenu("fragment scan in 1D");
      QMenu * msn_meta = new QMenu("fragment scan meta data");
      DPosition<2> p1 = widgetToData_(e->pos() + QPoint(10, 10));
      DPosition<2> p2 = widgetToData_(e->pos() - QPoint(10, 10));
      double rt_min = min(p1[1], p2[1]);
      double rt_max = max(p1[1], p2[1]);
      double mz_min = min(p1[0], p2[0]);
      double mz_max = max(p1[0], p2[0]);
      bool item_added = collectFragmentScansInArea(rt_min, rt_max, mz_min, mz_max, a, msn_scans, msn_meta);
      if (!item_added)
      {
        // Now simply go for the 5 closest points in RT and check whether there
        // are any scans.
        // NOTE: that if we go for the visible area, we run the
        // risk of iterating through *all* the scans.

        const AreaType & area = getVisibleArea();
        double mz_min_vis = area.minPosition()[0]; 
        double mz_max_vis = area.maxPosition()[0];

        double rt5s_min = getCurrentLayer().getPeakData()->getSpectra()[ indices[indices.size()-1] ].getRT();
        double rt5s_max = getCurrentLayer().getPeakData()->getSpectra()[ indices[0] ].getRT();

        item_added = collectFragmentScansInArea(rt5s_min, rt5s_max, mz_min_vis, mz_max_vis, a, msn_scans, msn_meta);

        if (!item_added)
        {
          // ok, now lets search the whole visible area (may be large!)
          double rt_min_vis = area.minPosition()[1]; 
          double rt_max_vis = area.maxPosition()[1];
          item_added = collectFragmentScansInArea(rt_min_vis, rt_max_vis, mz_min_vis, mz_max_vis, a, msn_scans, msn_meta);
        }
      }
      if (item_added)
      {
        context_menu->addMenu(msn_scans);
        context_menu->addMenu(msn_meta);
        context_menu->addSeparator();
      }

      finishContextMenu_(context_menu, settings_menu);

      //evaluate menu
      if ((result = context_menu->exec(mapToGlobal(e->pos()))))
      {
        if (result->parent() == ms1_scans  || result->parent() == msn_scans)
        {
          emit showSpectrumAs1D(result->data().toInt());
        }
        else if (result->parent() == ms1_meta || result->parent() == msn_meta)
        {
          showMetaData(true, result->data().toInt());
        }
      }
    }
    //-------------------FEATURES----------------------------------
    else if (layer.type == LayerData::DT_FEATURE)
    {
      //add settings
      settings_menu->addSeparator();
      settings_menu->addAction("Show/hide convex hull");
      settings_menu->addAction("Show/hide trace convex hulls");
      settings_menu->addAction("Show/hide numbers/labels");
      settings_menu->addAction("Show/hide unassigned peptide hits");

      // search for nearby features
      DPosition<2> p1 = widgetToData_(e->pos() + QPoint(10, 10));
      DPosition<2> p2 = widgetToData_(e->pos() - QPoint(10, 10));
      double rt_min = min(p1[1], p2[1]);
      double rt_max = max(p1[1], p2[1]);
      double mz_min = min(p1[0], p2[0]);
      double mz_max = max(p1[0], p2[0]);

      QMenu * meta = new QMenu("Feature meta data");
      bool present = false;
      FeatureMapType & features = *getCurrentLayer_().getFeatureMap();
      // feature meta data menu
      for (FeatureMapType::Iterator it = features.begin(); it != features.end(); ++it)
      {
        if (it->getMZ() <= mz_max && it->getMZ() >= mz_min && it->getRT() <= rt_max && it->getRT() >= rt_min)
        {
          present = true;
          a = meta->addAction(QString("RT: ") + QString::number(it->getRT()) + "  m/z:" + QString::number(it->getMZ()) + "  charge:" + QString::number(it->getCharge()));
          a->setData((int)(it - features.begin()));
        }
      }
      if (present)
      {
        context_menu->addMenu(meta);
        context_menu->addSeparator();
      }

      //add modifiable flag
      settings_menu->addSeparator();
      settings_menu->addAction("Toggle edit/view mode");

      finishContextMenu_(context_menu, settings_menu);

      //evaluate menu
      if ((result = context_menu->exec(mapToGlobal(e->pos()))))
      {
        if (result->text().left(3) == "RT:")
        {
          showMetaData(true, result->data().toInt());
        }
      }
    }
    //-------------------CONSENSUS FEATURES----------------------------------
    else if (layer.type == LayerData::DT_CONSENSUS)
    {
      //add settings
      settings_menu->addSeparator();
      settings_menu->addAction("Show/hide elements");

      //search for nearby features
      DPosition<2> p1 = widgetToData_(e->pos() + QPoint(10, 10));
      DPosition<2> p2 = widgetToData_(e->pos() - QPoint(10, 10));
      double rt_min = min(p1[1], p2[1]);
      double rt_max = max(p1[1], p2[1]);
      double mz_min = min(p1[0], p2[0]);
      double mz_max = max(p1[0], p2[0]);

      QMenu * consens_meta = new QMenu("Consensus meta data");
      bool present = false;
      ConsensusMapType & features = *getCurrentLayer_().getConsensusMap();
      //consensus feature meta data menu
      for (ConsensusMapType::Iterator it = features.begin(); it != features.end(); ++it)
      {
        if (it->getMZ() <= mz_max && it->getMZ() >= mz_min && it->getRT() <= rt_max && it->getRT() >= rt_min)
        {
          present = true;

          a = consens_meta->addAction(QString("RT: ") + QString::number(it->getRT()) + "  m/z:" + QString::number(it->getMZ()) + "  charge:" + QString::number(it->getCharge()));
          a->setData((int)(it - features.begin()));
        }
      }
      if (present)
      {
        context_menu->addMenu(consens_meta);
        context_menu->addSeparator();
      }

      finishContextMenu_(context_menu, settings_menu);

      if ((result = context_menu->exec(mapToGlobal(e->pos()))))
      {
        if (result->text().left(3) == "RT:")
        {
          showMetaData(true, result->data().toInt());
        }
      }
    }
    //------------------CHROMATOGRAMS----------------------------------
    else if (layer.type == LayerData::DT_CHROMATOGRAM)
    {
      settings_menu->addSeparator();
      settings_menu->addAction("Show/hide projections");
      settings_menu->addAction("Show/hide MS/MS precursors");

      PeakMap exp;
      exp = *layer.getPeakData();

      int CHROMATOGRAM_SHOW_MZ_RANGE = 10;

      // collect all precursor that fall into the mz rt window
      typedef std::set<Precursor, Precursor::MZLess> PCSetType;
      PCSetType precursor_in_rt_mz_window;
      for (vector<MSChromatogram >::const_iterator iter = exp.getChromatograms().begin(); iter != exp.getChromatograms().end(); ++iter)
      {
        if (mz + CHROMATOGRAM_SHOW_MZ_RANGE >= iter->getPrecursor().getMZ() &&
            mz - CHROMATOGRAM_SHOW_MZ_RANGE <= iter->getPrecursor().getMZ() &&
            rt >= iter->front().getRT() &&
            rt <= iter->back().getRT())
        {
          precursor_in_rt_mz_window.insert(iter->getPrecursor());
        }
      }

      // determine product chromatograms for each precursor
      map<Precursor, vector<Size>, Precursor::MZLess> map_precursor_to_chrom_idx;
      for (PCSetType::const_iterator pit = precursor_in_rt_mz_window.begin(); pit != precursor_in_rt_mz_window.end(); ++pit)
      {
        for (vector<MSChromatogram >::const_iterator iter = exp.getChromatograms().begin(); iter != exp.getChromatograms().end(); ++iter)
        {
          if (iter->getPrecursor() == *pit)
          {
            map_precursor_to_chrom_idx[*pit].push_back(iter - exp.getChromatograms().begin());
          }
        }
      }

      /*
        // dump precursor and associated products
        for (map<Precursor, vector<Size>, Precursor::MZLess >::iterator mit = map_precursor_to_chrom_idx.begin(); mit != map_precursor_to_chrom_idx.end(); ++mit)
        {
          cout << "Precursor: " << mit->first << endl;
          for (vector<Size>::iterator vit = mit->second.begin(); vit != mit->second.end(); ++vit)
          {
            cout << "  Product mz: " << exp.getChromatograms()[*vit].getMZ() << endl;
          }
        }
        */

      QMenu * msn_chromatogram  = nullptr;
      QMenu * msn_chromatogram_meta = nullptr;

      if (!map_precursor_to_chrom_idx.empty())
      {
        msn_chromatogram = context_menu->addMenu("Chromatogram");
        msn_chromatogram_meta = context_menu->addMenu("Chromatogram meta data");
        context_menu->addSeparator();

        for (map<Precursor, vector<Size>, Precursor::MZLess>::iterator mit = map_precursor_to_chrom_idx.begin(); mit != map_precursor_to_chrom_idx.end(); ++mit)
        {
          // Show the peptide sequence if available, otherwise show the m/z and charge only
          QString precursor_string = QString("Precursor m/z: (")  + String(mit->first.getCharge()).toQString() + ") " + QString::number(mit->first.getMZ());
          if (mit->first.metaValueExists("peptide_sequence"))
          {
            precursor_string = QString::number(mit->first.getMZ()) + " : " + String(mit->first.getMetaValue("peptide_sequence")).toQString() + " (" + QString::number(mit->first.getCharge()) + "+)";
          }
          QMenu * msn_precursor = msn_chromatogram->addMenu(precursor_string);  // neuer Eintrag für jeden Precursor

          // Show all: iterate over all chromatograms corresponding to the current precursor and add action containing all chromatograms
          a = msn_precursor->addAction(QString("Show all"));
          QList<QVariant> chroms_idx;
          for (vector<Size>::iterator vit = mit->second.begin(); vit != mit->second.end(); ++vit)
          {
            chroms_idx.push_back((unsigned int)*vit);
          }
          a->setData(chroms_idx);

          // Show single chromatogram: iterate over all chromatograms corresponding to the current precursor and add action for the single chromatogram
          for (vector<Size>::iterator vit = mit->second.begin(); vit != mit->second.end(); ++vit)
          {
            a = msn_precursor->addAction(QString("Chromatogram m/z: ") + QString::number(exp.getChromatograms()[*vit].getMZ())); // Precursor => Chromatogram MZ Werte eintragen
            a->setData((int)(*vit));
          }
        }
      }

      finishContextMenu_(context_menu, settings_menu);

      // show context menu and evaluate result
      if ((result = context_menu->exec(mapToGlobal(e->pos()))))
      {
        if (result->parent()->parent() == msn_chromatogram) // clicked on chromatogram entry (level 2)
        {
          if (result->text() == "Show all")
          {
            std::vector<int> chrom_indices;
            const QList<QVariant> & res = result->data().toList();
            for (Int i = 0; i != res.size(); ++i)
            {
              chrom_indices.push_back(res[i].toInt());
              cout << "chrom_indices: " << res[i].toInt() << std::endl;
            }
            emit showSpectrumAs1D(chrom_indices);
          }
          else   // Show single chromatogram
          {
            //cout << "Chromatogram result " << result->data().toInt() << endl;
            emit showSpectrumAs1D(result->data().toInt());
          }
        }
        else if (result->parent() == msn_chromatogram_meta)
        {
          showMetaData(true, result->data().toInt());
        }
      }
    }

    //common actions of peaks and features
    if (result)
    {
      if (result->text() == "Preferences")
      {
        showCurrentLayerPreferences();
      }
      else if (result->text() == "Show/hide grid lines")
      {
        showGridLines(!gridLinesShown());
      }
      else if (result->text() == "Show/hide axis legends")
      {
        emit changeLegendVisibility();
      }
      else if (result->text() == "Layer" || result->text() == "Visible layer data")
      {
        saveCurrentLayer(result->text() == "Visible layer data");
      }
      else if (result->text() == "As image")
      {
        spectrum_widget_->saveAsImage();
      }
      else if (result->text() == "Show/hide projections")
      {
        emit toggleProjections();
      }
      else if (result->text() == "Show/hide MS/MS precursors")
      {
        setLayerFlag(LayerData::P_PRECURSORS, !getLayerFlag(LayerData::P_PRECURSORS));
      }
      else if (result->text() == "Show/hide convex hull")
      {
        setLayerFlag(LayerData::F_HULL, !getLayerFlag(LayerData::F_HULL));
      }
      else if (result->text() == "Show/hide trace convex hulls")
      {
        setLayerFlag(LayerData::F_HULLS, !getLayerFlag(LayerData::F_HULLS));
      }
      else if (result->text() == "Show/hide unassigned peptide hits")
      {
        setLayerFlag(LayerData::F_UNASSIGNED, !getLayerFlag(LayerData::F_UNASSIGNED));
      }
      else if (result->text() == "Show/hide numbers/labels")
      {
        if (layer.label == LayerData::L_NONE)
        {
          getCurrentLayer_().label = LayerData::L_META_LABEL;
        }
        else
        {
          getCurrentLayer_().label = LayerData::L_NONE;
        }
      }
      else if (result->text() == "Toggle edit/view mode")
      {
        getCurrentLayer_().modifiable = !getCurrentLayer_().modifiable;
      }
      else if (result->text() == "Show/hide elements")
      {
        setLayerFlag(LayerData::C_ELEMENTS, !getLayerFlag(LayerData::C_ELEMENTS));
      }
      else if (result->text() == "Layer meta data")
      {
        showMetaData(true);
      }
      else if (result->text() == "Switch to 3D view")
      {
        emit showCurrentPeaksAs3D();
      }
    }

    e->accept();
  }

  void Spectrum2DCanvas::finishContextMenu_(QMenu * context_menu, QMenu * settings_menu)
  {
    //finish settings menu
    settings_menu->addSeparator();
    settings_menu->addAction("Preferences");

    //create save menu
    QMenu * save_menu = new QMenu("Save");
    save_menu->addAction("Layer");
    save_menu->addAction("Visible layer data");
    save_menu->addAction("As image");

    //add settings menu
    context_menu->addMenu(save_menu);
    context_menu->addMenu(settings_menu);

    //add external context menu
    if (context_add_)
    {
      context_menu->addSeparator();
      context_menu->addMenu(context_add_);
    }
  }

  void Spectrum2DCanvas::showCurrentLayerPreferences()
  {
    Internal::Spectrum2DPrefDialog dlg(this);
    LayerData & layer = getCurrentLayer_();

    ColorSelector * bg_color = dlg.findChild<ColorSelector *>("bg_color");
    QComboBox * mapping = dlg.findChild<QComboBox *>("mapping");
    MultiGradientSelector * gradient = dlg.findChild<MultiGradientSelector *>("gradient");
    QComboBox * feature_icon = dlg.findChild<QComboBox *>("feature_icon");
    QSpinBox * feature_icon_size = dlg.findChild<QSpinBox *>("feature_icon_size");

    bg_color->setColor(QColor(param_.getValue("background_color").toQString()));
    if (isMzToXAxis())
    {
      mapping->setCurrentIndex(0);
    }
    else
    {
      mapping->setCurrentIndex(1);
    }
    gradient->gradient().fromString(layer.param.getValue("dot:gradient"));
    feature_icon->setCurrentIndex(feature_icon->findText(layer.param.getValue("dot:feature_icon").toQString()));
    feature_icon_size->setValue((int)layer.param.getValue("dot:feature_icon_size"));

    if (dlg.exec())
    {
      param_.setValue("background_color", bg_color->getColor().name());
      layer.param.setValue("dot:feature_icon", feature_icon->currentText());
      layer.param.setValue("dot:feature_icon_size", feature_icon_size->value());
      if ((mapping->currentIndex() == 0 && !isMzToXAxis()) || (mapping->currentIndex() == 1 && isMzToXAxis()))
      {
        mzToXAxis(!isMzToXAxis());
      }
      layer.param.setValue("dot:gradient", gradient->gradient().toString());

      emit preferencesChange();
    }
  }

  void Spectrum2DCanvas::currentLayerParametersChanged_()
  {
    recalculateDotGradient_(activeLayerIndex());

    update_buffer_ = true;
    update_(OPENMS_PRETTY_FUNCTION);
  }

  void Spectrum2DCanvas::saveCurrentLayer(bool visible)
  {
    const LayerData & layer = getCurrentLayer();

    //determine proposed filename
    String proposed_name = param_.getValue("default_path");
    if (visible == false && layer.filename != "")
    {
      proposed_name = layer.filename;
    }

    if (layer.type == LayerData::DT_PEAK)   //peak data
    {
      QString selected_filter = "";
      QString file_name = QFileDialog::getSaveFileName(this, "Save file", proposed_name.toQString(), "mzML files (*.mzML);;mzData files (*.mzData);;mzXML files (*.mzXML);;All files (*)", &selected_filter);
      if (!file_name.isEmpty())
      {
        // check whether a file type suffix has been given
        // first check mzData and mzXML then mzML
        // if the setting is at "All files"
        // mzML will be used
        String upper_filename = file_name;
        upper_filename.toUpper();
        if (selected_filter == "mzData files (*.mzData)")
        {
          if (!upper_filename.hasSuffix(".MZDATA"))
          {
            file_name += ".mzData";
          }
        }
        else if (selected_filter == "mzXML files (*.mzXML)")
        {
          if (!upper_filename.hasSuffix(".MZXML"))
          {
            file_name += ".mzXML";
          }
        }
        else
        {
          if (!upper_filename.hasSuffix(".MZML"))
          {
            file_name += ".mzML";
          }
        }

        if (visible)     //only visible data
        {
          ExperimentType out;
          getVisiblePeakData(out);
          addDataProcessing_(out, DataProcessing::FILTERING);
          FileHandler().storeExperiment(file_name, out, ProgressLogger::GUI);
        }
        else         //all data
        {
          FileHandler().storeExperiment(file_name, *layer.getPeakData(), ProgressLogger::GUI);
        }
        modificationStatus_(activeLayerIndex(), false);
      }
    }
    else if (layer.type == LayerData::DT_FEATURE) //features
    {
      QString file_name = QFileDialog::getSaveFileName(this, "Save file", proposed_name.toQString(), "featureXML files (*.featureXML);;All files (*)");
      if (!file_name.isEmpty())
      {
        // add suffix ".featureXML" if not given
        String upper_filename = file_name;
        upper_filename.toUpper();
        if (!upper_filename.hasSuffix(".FEATUREXML"))
        {
          file_name += ".featureXML";
        }
        if (visible)     //only visible data
        {
          FeatureMapType out;
          getVisibleFeatureData(out);
          FeatureXMLFile().store(file_name, out);
        }
        else         //all data
        {
          FeatureXMLFile().store(file_name, *layer.getFeatureMap());
        }
        modificationStatus_(activeLayerIndex(), false);
      }
    }
    else if (layer.type == LayerData::DT_CONSENSUS) //consensus feature data
    {
      QString file_name = QFileDialog::getSaveFileName(this, "Save file", proposed_name.toQString(), "consensusXML files (*.consensusXML);;All files (*)");
      if (!file_name.isEmpty())
      {
        // add suffix ".consensusXML" if not given
        String upper_filename = file_name;
        upper_filename.toUpper();
        if (!upper_filename.hasSuffix(".CONSENSUSXML"))
        {
          file_name += ".consensusXML";
        }

        if (visible)     //only visible data
        {
          ConsensusMapType out;
          getVisibleConsensusData(out);
          ConsensusXMLFile().store(file_name, out);
        }
        else         //all data
        {
          ConsensusXMLFile().store(file_name, *layer.getConsensusMap());
        }
        modificationStatus_(activeLayerIndex(), false);
      }
    }
    else if (layer.type == LayerData::DT_CHROMATOGRAM) //chromatograms
    {
      //TODO CHROM
    }
  }

  void Spectrum2DCanvas::updateLayer(Size i)
  {
    //update nearest peak
    selected_peak_.clear();
    recalculateRanges_(0, 1, 2);
    resetZoom(false);     //no repaint as this is done in intensityModeChange_() anyway
    intensityModeChange_();
    modificationStatus_(i, false);
  }

  void Spectrum2DCanvas::translateVisibleArea_( double mzShiftRel, double rtShiftRel )
  {
    double rtShift = rtShiftRel * visible_area_.height();
    double mzShift = mzShiftRel * visible_area_.width();
    AreaType newArea( visible_area_ );
    // shift the visible area avoiding moving out data range bounds
    if ( mzShift > 0 ) {
        newArea.setMaxX( qMin( overall_data_range_.maxX(), visible_area_.maxX() + mzShift ) );
        newArea.setMinX( qMax( overall_data_range_.minX(), newArea.maxX() - visible_area_.width() ) );
    } else {
        newArea.setMinX( qMax( overall_data_range_.minX(), visible_area_.minX() + mzShift ) );
        newArea.setMaxX( qMin( overall_data_range_.maxX(), newArea.minX() + visible_area_.width() ) );
    }
    if ( rtShift > 0 ) {
        newArea.setMaxY( qMin( overall_data_range_.maxY(), visible_area_.maxY() + rtShift ) );
        newArea.setMinY( qMax( overall_data_range_.minY(), newArea.maxY() - visible_area_.height() ) );
    } else {
        newArea.setMinY( qMax( overall_data_range_.minY(), visible_area_.minY() + rtShift ) );
        newArea.setMaxY( qMin( overall_data_range_.maxY(), newArea.minY() + visible_area_.height() ) );
    }
    //change visible area
    changeVisibleArea_(newArea);
    emit layerZoomChanged(this);
  }

  void Spectrum2DCanvas::translateLeft_(Qt::KeyboardModifiers /*m*/)
  {
    if ( isMzToXAxis() ) translateVisibleArea_( -0.05, 0.0 );
    else translateVisibleArea_( 0.0, -0.05 );
  }

  void Spectrum2DCanvas::translateRight_(Qt::KeyboardModifiers /*m*/)
  {
    if ( isMzToXAxis() ) translateVisibleArea_( 0.05, 0.0 );
    else translateVisibleArea_( 0.0, 0.05 );
  }

  void Spectrum2DCanvas::translateForward_()
  {
    if ( isMzToXAxis() ) translateVisibleArea_( 0.0, 0.05 );
    else translateVisibleArea_( 0.05, 0.0 );
  }

  void Spectrum2DCanvas::translateBackward_()
  {
    if ( isMzToXAxis() ) translateVisibleArea_( 0.0, -0.05 );
    else translateVisibleArea_( -0.05, 0.0 );
  }

  void Spectrum2DCanvas::keyPressEvent(QKeyEvent * e)
  {
    // CTRL+ALT (exactly, not e.g. CTRL+ALT+KEYPAD<X>|SHIFT...)
    // note that Qt::KeypadModifier is also a modifier which gets activated when any keypad key is pressed
    if (e->modifiers() == (Qt::ControlModifier | Qt::AltModifier))
    {
      String status_changed;
      // +Home (MacOSX small keyboard: Fn+ArrowLeft) => increase point size
      if ((e->key() == Qt::Key_Home) && (pen_size_max_ < PEN_SIZE_MAX_LIMIT))
      { 
        ++pen_size_max_;
        status_changed = "Max. dot size increased to '" + String(pen_size_max_) + "'";
      }
      // +End (MacOSX small keyboard: Fn+ArrowRight) => decrease point size
      else if ((e->key() == Qt::Key_End) && (pen_size_max_ > PEN_SIZE_MIN_LIMIT))
      {
        --pen_size_max_;
        status_changed = "Max. dot size decreased to '" + String(pen_size_max_) + "'";
      }
      // +PageUp => increase min. coverage threshold
      else if (e->key() == Qt::Key_PageUp && canvas_coverage_min_ < CANVAS_COVERAGE_MIN_LIMITHIGH)
      {
        canvas_coverage_min_ += 0.05; // 5% steps
        status_changed = "Min. coverage threshold increased to '" + String(canvas_coverage_min_) + "'";
      }
      // +PageDown => decrease min. coverage threshold
      else if (e->key() == Qt::Key_PageDown && canvas_coverage_min_ > CANVAS_COVERAGE_MIN_LIMITLOW)
      {
        canvas_coverage_min_ -= 0.05; // 5% steps
        status_changed = "Min. coverage threshold decreased to '" + String(canvas_coverage_min_) + "'";
      }
      if (!status_changed.empty())
      {
        emit sendStatusMessage(status_changed, 0);
        update_buffer_ = true; // full repaint
        update_(OPENMS_PRETTY_FUNCTION); // schedule repaint
        return;
      }
    }

    // Delete features
    LayerData& layer = getCurrentLayer_();
    if (e->key() == Qt::Key_Delete && getCurrentLayer().modifiable && layer.type == LayerData::DT_FEATURE && selected_peak_.isValid())
    {
      layer.getFeatureMap()->erase(layer.getFeatureMap()->begin() + selected_peak_.peak);
      selected_peak_.clear();
      update_buffer_ = true;
      update_(OPENMS_PRETTY_FUNCTION);
      modificationStatus_(activeLayerIndex(), true);
      return;
    }

    // call parent class
    SpectrumCanvas::keyPressEvent(e);
  }

  void Spectrum2DCanvas::keyReleaseEvent(QKeyEvent * e)
  {
    //zoom if in zoom mode and a valid rectangle is selected
    if (action_mode_ == AM_ZOOM && rubber_band_.isVisible())
    {
      rubber_band_.hide();
      QRect rect = rubber_band_.geometry();
      if (rect.width() != 0 && rect.height() != 0)
      {
        AreaType area(widgetToData_(rect.topLeft()), widgetToData_(rect.bottomRight()));
        changeVisibleArea_(area, true, true);
        emit layerZoomChanged(this);
      }
    }
    else if (action_mode_ == AM_MEASURE)
    {
      measurement_start_.clear();
      update_(OPENMS_PRETTY_FUNCTION);
    }

    // do the normal stuff
    SpectrumCanvas::keyReleaseEvent(e);
  }

  void Spectrum2DCanvas::mouseDoubleClickEvent(QMouseEvent * e)
  {
    LayerData & current_layer = getCurrentLayer_();

    if (current_layer.modifiable && current_layer.type == LayerData::DT_FEATURE)
    {
      Feature tmp;
      if (selected_peak_.isValid())       //edit existing feature
      {
        FeatureEditDialog dialog(this);
        dialog.setFeature((*current_layer.getFeatureMap())[selected_peak_.peak]);
        if (dialog.exec())
        {
          tmp = dialog.getFeature();
          (*current_layer.getFeatureMap())[selected_peak_.peak] = tmp;
        }
      }
      else       //create new feature
      {
        tmp.setRT(widgetToData_(e->pos())[1]);
        tmp.setMZ(widgetToData_(e->pos())[0]);
        FeatureEditDialog dialog(this);
        dialog.setFeature(tmp);
        if (dialog.exec())
        {
          tmp = dialog.getFeature();
          current_layer.getFeatureMap()->push_back(tmp);
        }
      }

      // update gradient if the min/max intensity changes
      if (tmp.getIntensity() < current_layer.getFeatureMap()->getMinInt() || tmp.getIntensity() > current_layer.getFeatureMap()->getMaxInt())
      {
        current_layer.getFeatureMap()->updateRanges();
        recalculateRanges_(0, 1, 2);
        intensityModeChange_();
      }
      else // just repaint to show the changes
      {
        update_buffer_ = true;
        update_(OPENMS_PRETTY_FUNCTION);
      }

      modificationStatus_(activeLayerIndex(), true);
    }
  }

  void Spectrum2DCanvas::mergeIntoLayer(Size i, FeatureMapSharedPtrType map)
  {
    OPENMS_PRECONDITION(i < layers_.size(), "Spectrum2DCanvas::mergeIntoLayer(i, map) index overflow");
    OPENMS_PRECONDITION(layers_[i].type == LayerData::DT_FEATURE, "Spectrum2DCanvas::mergeIntoLayer(i, map) non-feature layer selected");
    //reserve enough space
    layers_[i].getFeatureMap()->reserve(layers_[i].getFeatureMap()->size() + map->size());
    //add features
    for (Size j = 0; j < map->size(); ++j)
    {
      layers_[i].getFeatureMap()->push_back((*map)[j]);
    }
    //update the layer and overall ranges (if necessary)
    RangeManager<2>::PositionType min_pos_old = layers_[i].getFeatureMap()->getMin();
    RangeManager<2>::PositionType max_pos_old = layers_[i].getFeatureMap()->getMax();
    double min_int_old = layers_[i].getFeatureMap()->getMinInt();
    double max_int_old = layers_[i].getFeatureMap()->getMaxInt();
    layers_[i].getFeatureMap()->updateRanges();
    if (min_pos_old > layers_[i].getFeatureMap()->getMin() || max_pos_old < layers_[i].getFeatureMap()->getMax())
    {
      recalculateRanges_(0, 1, 2);
      resetZoom(true);
    }
    if (min_int_old > layers_[i].getFeatureMap()->getMinInt() || max_int_old < layers_[i].getFeatureMap()->getMaxInt())
    {
      intensityModeChange_();
    }
  }

  void Spectrum2DCanvas::mergeIntoLayer(Size i, ConsensusMapSharedPtrType map)
  {
    OPENMS_PRECONDITION(i < layers_.size(), "Spectrum2DCanvas::mergeIntoLayer(i, map) index overflow");
    OPENMS_PRECONDITION(layers_[i].type == LayerData::DT_CONSENSUS, "Spectrum2DCanvas::mergeIntoLayer(i, map) non-consensus-feature layer selected");
    //reserve enough space
    layers_[i].getConsensusMap()->reserve(layers_[i].getFeatureMap()->size() + map->size());
    //add features
    for (Size j = 0; j < map->size(); ++j)
    {
      layers_[i].getConsensusMap()->push_back((*map)[j]);
    }
    //update the layer and overall ranges (if necessary)
    RangeManager<2>::PositionType min_pos_old = layers_[i].getConsensusMap()->getMin();
    RangeManager<2>::PositionType max_pos_old = layers_[i].getConsensusMap()->getMax();
    double min_int_old = layers_[i].getConsensusMap()->getMinInt();
    double max_int_old = layers_[i].getConsensusMap()->getMaxInt();
    layers_[i].getConsensusMap()->updateRanges();
    if (min_pos_old > layers_[i].getConsensusMap()->getMin() || max_pos_old < layers_[i].getConsensusMap()->getMax())
    {
      recalculateRanges_(0, 1, 2);
      resetZoom(true);
    }
    if (min_int_old > layers_[i].getConsensusMap()->getMinInt() || max_int_old < layers_[i].getConsensusMap()->getMaxInt())
    {
      intensityModeChange_();
    }
  }

  void Spectrum2DCanvas::mergeIntoLayer(Size i, vector<PeptideIdentification> & peptides)
  {
    OPENMS_PRECONDITION(i < layers_.size(), "Spectrum2DCanvas::mergeIntoLayer(i, peptides) index overflow");
    OPENMS_PRECONDITION(layers_[i].type == LayerData::DT_IDENT, "Spectrum2DCanvas::mergeIntoLayer(i, peptides) non-identification layer selected");
    // reserve enough space
    layers_[i].peptides.reserve(layers_[i].peptides.size() + peptides.size());
    // insert peptides
    layers_[i].peptides.insert(layers_[i].peptides.end(), peptides.begin(),
                               peptides.end());
    // update the layer and overall ranges
    recalculateRanges_(0, 1, 2);
    resetZoom(true);
  }

    bool Spectrum2DCanvas::collectFragmentScansInArea(double rt_min, double rt_max, double mz_min, double mz_max, QAction* a, QMenu * msn_scans, QMenu * msn_meta)
    {

      bool item_added = false;
      for (ExperimentType::ConstIterator it = getCurrentLayer().getPeakData()->RTBegin(rt_min); it != getCurrentLayer().getPeakData()->RTEnd(rt_max); ++it)
      {
        double mz = 0.0;
        if (!it->getPrecursors().empty())
        {
          mz = it->getPrecursors()[0].getMZ();
        }

        if (it->getMSLevel() > 1 && mz >= mz_min && mz <= mz_max)
        {
          a = msn_scans->addAction(QString("RT: ") + QString::number(it->getRT()) + " mz: " + QString::number(mz));
          a->setData((int)(it - getCurrentLayer().getPeakData()->begin()));
          a = msn_meta->addAction(QString("RT: ") + QString::number(it->getRT()) + " mz: " + QString::number(mz));
          a->setData((int)(it - getCurrentLayer().getPeakData()->begin()));
          item_added = true;
        }
      }
      return item_added;
    }

} //namespace OpenMS
