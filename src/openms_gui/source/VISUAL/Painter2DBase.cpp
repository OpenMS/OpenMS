// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------


#include <OpenMS/VISUAL/Painter2DBase.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <OpenMS/VISUAL/LayerDataChrom.h>
#include <OpenMS/VISUAL/LayerDataConsensus.h>
#include <OpenMS/VISUAL/LayerDataIdent.h>
#include <OpenMS/VISUAL/LayerDataIonMobility.h>
#include <OpenMS/VISUAL/LayerDataFeature.h>
#include <OpenMS/VISUAL/LayerDataPeak.h>

#include <OpenMS/VISUAL/Plot2DCanvas.h>

#include <QColor>
#include <QPainter>
#include <QPoint>

using namespace std;

namespace OpenMS
{


  void Painter2DBase::highlightElement(QPainter* /*painter*/, Plot2DCanvas* /*canvas*/, const PeakIndex /*element*/)
  {
    // does nothing by default
  }

  void Painter2DBase::paintConvexHull_(QPainter& painter, Plot2DCanvas* canvas, const ConvexHull2D& hull, bool has_identifications)
  {
    QPolygon points;
    ConvexHull2D::PointArrayType ch_points = hull.getHullPoints();
    points.resize((int)ch_points.size());
    UInt index = 0;
    QPoint pos;
    // iterate over hull points
    for (ConvexHull2D::PointArrayType::const_iterator it = ch_points.begin(); it != ch_points.end(); ++it, ++index)
    {
      Peak2D ms_peak({it->getX(), it->getY()}, 0); // assume that CH of a feature is RT and m/z
      points.setPoint(index, canvas->dataToWidget_(canvas->unit_mapper_.map(ms_peak)));
    }
    painter.setPen(QPen(Qt::white, 5, Qt::DotLine, Qt::RoundCap, Qt::RoundJoin));
    painter.drawPolygon(points);
    painter.setPen(QPen(has_identifications ? Qt::green : Qt::blue, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
    painter.drawPolygon(points);
  }

  void Painter2DBase::paintConvexHulls_(QPainter& painter, Plot2DCanvas* canvas, const std::vector<ConvexHull2D>& hulls, bool has_identifications)
  {
    for (const auto& hull : hulls)
    {
      paintConvexHull_(painter, canvas, hull, has_identifications);
    }
  }

  void Painter2DBase::paintPeptideIDs_(QPainter* painter, Plot2DCanvas* canvas, const IPeptideIds::PepIds& ids, int layer_index)
  {
    painter->setPen(Qt::darkRed);
    bool show_labels = canvas->getLayerFlag(layer_index, LayerDataBase::I_LABELS);

    for (const auto& id : ids)
    {
      if (!id.getHits().empty() || show_labels)
      {
        if (!id.hasRT() || !id.hasMZ())
        {
          // TODO: show error message here
          continue;
        }
        double rt = id.getRT();
        if (!canvas->visible_area_.getAreaUnit().containsRT(rt))
        {
          continue;
        }
        double mz = canvas->getIdentificationMZ_(layer_index, id);
        if (!canvas->visible_area_.getAreaUnit().containsMZ(mz))
        {
          continue;
        }
        // draw dot
        QPoint pos = canvas->dataToWidget_(canvas->unit_mapper_.map(id));
        painter->drawLine(pos.x(), pos.y() - 1.0, pos.x(), pos.y() + 1.0);
        painter->drawLine(pos.x() - 1.0, pos.y(), pos.x() + 1.0, pos.y());

        // draw sequence
        String sequence;
        if (show_labels)
        {
          sequence = id.getMetaValue("label");
        }
        else
        {
          sequence = id.getHits()[0].getSequence().toString();
        }
        if (sequence.empty() && !id.getHits().empty())
        {
          sequence = id.getHits()[0].getMetaValue("label");
        }
        if (id.getHits().size() > 1)
          sequence += "...";
        painter->drawText(pos.x() + 10, pos.y() + 10, sequence.toQString());
      }
    }
  }

  // ###############################
  // ###### 2D Peak
  // ###############################

  Painter2DPeak::Painter2DPeak(const LayerDataPeak* parent) : layer_(parent)
  {
  }

  void Painter2DPeak::paint(QPainter* painter, Plot2DCanvas* canvas, int layer_index)
  {
    // renaming some values for readability
    const auto& peak_map = *layer_->getPeakData();

    // skip empty peak maps
    if (peak_map.empty())
    {
      return;
    }

    const auto [rt_min, rt_max] = canvas->visible_area_.getAreaUnit().RangeRT::getNonEmptyRange();
    const auto [mz_min, mz_max] = canvas->visible_area_.getAreaUnit().RangeMZ::getNonEmptyRange();
    const auto [im_min, im_max] = canvas->visible_area_.getAreaUnit().RangeMobility::getNonEmptyRange();

    // do we currently show an IM frame (with IM + m/z) as units?
    const bool is_IM_frame = peak_map.isIMFrame();

    auto is_visible_scan = [&](MSExperiment::ConstIterator it_scan) {
      if (it_scan->size() <= 1) return false;
      // an IM scan? (where we do not care about MS level)
      if (is_IM_frame) { return true; }
      // for 'standard' RT, m/z data
      return it_scan->getMSLevel() == 1;
    };

    //-----------------------------------------------------------------------------------------------
    // Determine number of shown scans (MS1)
    std::vector<Size> scan_indices; // list of visible RT/IM scans with at least 2 points
    const auto rt_end = peak_map.RTEnd(rt_max);
    for (auto it = peak_map.RTBegin(rt_min); it != rt_end; ++it)
    {
      if (is_visible_scan(it) && Math::contains(it->getDriftTime(), im_min, im_max))
      {
        scan_indices.push_back(std::distance(peak_map.begin(), it));
      }
    }
    Size n_ms1_scans = scan_indices.size();

    if (n_ms1_scans > 0)
    {
      // sample #of points at 3 scan locations (25%, 50%, 75%)
      // and take the median value
      Size n_peaks_in_scan(0);
      {
        static constexpr std::array<double, 3> quantiles = {0.25, 0.50, 0.75};
        std::array<Size, quantiles.size()> n_s;
        for (Size i = 0; i < quantiles.size(); ++i)
        {
          const auto& spec = peak_map[scan_indices[n_ms1_scans * quantiles[i]]];
          if (!spec.isSorted())
            throw Exception::BaseException();
          n_s[i] = std::distance(spec.MZBegin(mz_min), spec.MZEnd(mz_max));
        }
        std::sort(n_s.begin(), n_s.end());
        n_peaks_in_scan = n_s[1]; // median
      }

      Size rt_pixel_count, mz_pixel_count;
      { // obtain number of pixels in RT/IM and MZ dimension
        auto tmp = canvas->getPixelRange().getAreaUnit();
        rt_pixel_count = tmp.RangeRT::isEmpty() ? tmp.getMaxMobility() : tmp.getMaxRT();
        mz_pixel_count = tmp.getMaxMZ();
      }
      double ratio_data2pixel_rt = n_ms1_scans / (double)rt_pixel_count;
      double ratio_data2pixel_mz = n_peaks_in_scan / (double)mz_pixel_count;

      // minimum fraction of image expected to be filled with data
      // if not reached, we upscale point size
      bool has_low_pixel_coverage = ratio_data2pixel_rt < canvas->canvas_coverage_min_ || ratio_data2pixel_mz < canvas->canvas_coverage_min_;

      // Are several peaks expected to be drawn on the same pixel in either RT or m/z?
      // --> thin out and show only maxima
      // Also, we cannot upscale in this mode (since we operate on the buffer directly, i.e. '1 data point == 1 pixel'
      if (!has_low_pixel_coverage && (n_peaks_in_scan > mz_pixel_count || n_ms1_scans > rt_pixel_count))
      {
        paintMaximumIntensities_(*painter, canvas, layer_index, rt_pixel_count, mz_pixel_count);
      }
      else
      { // this is slower to paint, but allows scaling points
        // when data is zoomed in to single peaks these are visualized as circles

        // compute ideal pen width (from data);
        // since points are rectangular, we take the value of the most "crowded" dimension
        // i.e. so that adjacent points do not overlap
        double pen_width = std::min(1 / ratio_data2pixel_rt, 1 / ratio_data2pixel_mz);
        // ... and make sure its within our boundaries
        pen_width = std::max(pen_width, canvas->pen_size_min_);
        pen_width = std::min(pen_width, canvas->pen_size_max_);
#ifdef DEBUG_TOPPVIEW
        std::cerr << "pen-width " << pen_width << "\n";
#endif
        // However: if one dimension is sparse (e.g. only a few, but very long scans), we want to
        //          avoid showing lots of white background by increasing point size
        // This might lead to 'overplotting', but the paint method below can deal with it since
        // it will paint high intensities last.
        canvas->adaptPenScaling_(ratio_data2pixel_mz, pen_width);
        canvas->adaptPenScaling_(ratio_data2pixel_rt, pen_width);
#ifdef DEBUG_TOPPVIEW
        std::cerr << "new pen: " << pen_width << "\n";
#endif

        // few data points expected: more expensive drawing of all data points (circles or points depending on zoom level)
        paintAllIntensities_(*painter, canvas, layer_index, pen_width);
      }

    } // end of no-scans check

    //-----------------------------------------------------------------
    // draw precursor peaks
    if (canvas->getLayerFlag(layer_index, LayerDataBase::P_PRECURSORS))
    {
      paintPrecursorPeaks_(*painter, canvas);
    }
  }

  void Painter2DPeak::paintAllIntensities_(QPainter& painter, Plot2DCanvas* canvas, Size layer_index, double pen_width)
  {
    QVector<QPolygon> coloredPoints((int)layer_->gradient.precalculatedSize());

    const double snap_factor = canvas->snap_factors_[layer_index];
    const auto& map = *layer_->getPeakData();
    const auto& area = canvas->visible_area_.getAreaUnit();
    const auto end_area = map.areaEndConst();
    // for IM data, use whatever is there. For RT/mz data, use MSlevel 1
    const UInt MS_LEVEL = (! map.empty() && map.isIMFrame()) ? map[0].getMSLevel() : 1;
    for (auto i = map.areaBeginConst(area, MS_LEVEL); i != end_area; ++i)
    {
      PeakIndex pi = i.getPeakIndex();
      if (layer_->filters.passes(map[pi.spectrum], pi.peak))
      {
        auto from = canvas->unit_mapper_.map(i);
        QPoint pos = canvas->dataToWidget_(from);
        // store point in the array of its color
        Int colorIndex = canvas->precalculatedColorIndex_(i->getIntensity(), layer_->gradient, snap_factor);
        coloredPoints[colorIndex].push_back(pos);
      }
    }

    // draw point arrays from minimum to maximum intensity,
    // avoiding low-intensity points obscuring the high-intensity ones
    painter.save();
    QPen newPointsPen;
    newPointsPen.setWidthF(pen_width);
    for (Int colorIx = 0; colorIx < coloredPoints.size(); colorIx++)
    {
      const QPolygon& pointsArr = coloredPoints[colorIx];
      if (pointsArr.size())
      {
        newPointsPen.setColor(layer_->gradient.precalculatedColorByIndex(colorIx));
        painter.setPen(newPointsPen);
        painter.drawPoints(pointsArr);
      }
    }
    painter.restore();
  }

  /// abstract away the RT or IM dimension (so we can iterate over both using the same code in paintMaximumIntensities_())
  struct DimInfo
  {
    DimInfo(const MSExperiment& map, const RangeBase& dim, const DimMapper<2>& mapper) : exp_(map), dim_(dim), mapper_(mapper)
    {
    }
    const MSExperiment& exp_;
    const RangeBase& dim_;
    const DimMapper<2>& mapper_;
    /// Maximum RT or IM value for the whole layer
    virtual double getMaximum() const = 0;
    /// get the RT or IM from a spectrum
    virtual double getValue(const MSSpectrum& spec) const = 0;
    /// Iterator to first scan with this RT/IM value
    virtual MSExperiment::ConstIterator getFirstScan(double value) const = 0;
    /// Map a pair of RT/mz (or IM/mz) to the XY-plane
    virtual DimMapper<2>::Point mapToPoint(double value, double mz) const = 0;
  };
  struct DimInfoRT : DimInfo
  {
    using DimInfo::DimInfo; // inherit C'tor
    double getMaximum() const override
    {
      return exp_.getMaxRT();
    }
    double getValue(const MSSpectrum& spec) const override
    {
      return spec.getRT();
    }
    MSExperiment::ConstIterator getFirstScan(double value) const override
    {
      return exp_.RTBegin(value);
    }
    DimMapper<2>::Point mapToPoint(double value, double mz) const override
    {
      return mapper_.map(Peak2D({value, mz}, 0));
    }
  };
  struct DimInfoIM : DimInfo {
    using DimInfo::DimInfo; // inherit C'tor
    double getMaximum() const override
    {
      return exp_.getMaxMobility();
    }
    double getValue(const MSSpectrum& spec) const override
    {
      return spec.getDriftTime();
    }
    MSExperiment::ConstIterator getFirstScan(double value) const override
    {
      return exp_.IMBegin(value);
    }
    DimMapper<2>::Point mapToPoint(double value, double mz) const override
    {
      // there is no datastructure 
      return mapper_.map(MobilityPeak2D({value, mz}, 0));
    }
  };



  void Painter2DPeak::paintMaximumIntensities_(QPainter& painter, Plot2DCanvas* canvas, Size layer_index, Size rt_pixel_count, Size mz_pixel_count)
  {
    // set painter to black (we operate directly on the pixels for all colored data)
    painter.setPen(Qt::black);
    const double snap_factor = canvas->snap_factors_[layer_index];
    const auto& map = *layer_->getPeakData();
    const auto& area = canvas->visible_area_.getAreaUnit();

    // for IM data, use whatever is there. For RT/mz data, use MSlevel 1
    const UInt MS_LEVEL = (! map.empty() && map.isIMFrame()) ? map[0].getMSLevel() : 1;

    auto RT_or_IM_paint = [&](const DimInfo& mapper) {
      // note: the variables are named, assuming we have an RT+mz canvas.
      //       However, by using 'DimInfo' this could well be an IM+mz canvas (i.e. an IM Frame)
      
      const double rt_min = mapper.dim_.getMin();
      const double rt_max = mapper.dim_.getMax();
      const double mz_min = area.getMinMZ();
      const double mz_max = area.getMaxMZ();

      // calculate pixel size in data coordinates
      double rt_step_size = (rt_max - rt_min) / rt_pixel_count;
      double mz_step_size = (mz_max - mz_min) / mz_pixel_count;

      // start at first visible RT scan
      Size scan_index = std::distance(map.begin(), mapper.getFirstScan(rt_min));
      vector<Size> scan_indices, peak_indices;
      // iterate over all pixels (RT dimension)
      for (Size rt = 0; rt < rt_pixel_count; ++rt)
      {
        // interval in data coordinates for the current pixel
        double rt_start = rt_min + rt_step_size * rt;
        double rt_end = rt_start + rt_step_size;
        // cout << "rt: " << rt << " (" << rt_start << " - " << rt_end << ")" << endl;

        // reached the end of data
        if (rt_end >= mapper.getMaximum())
        {
          break;
        }
        // determine the relevant spectra and reserve an array for the peak indices
        scan_indices.clear();
        peak_indices.clear();
        for (Size i = scan_index; i < map.size(); ++i)
        {
          const auto& spec = map[i];
          if (mapper.getValue(spec) >= rt_end)
          {
            scan_index = i; // store last scan index for next RT pixel
            break;
          }
          if (spec.getMSLevel() == MS_LEVEL && ! spec.empty())
          {
            scan_indices.push_back(i);
            peak_indices.push_back(spec.MZBegin(mz_min) - spec.begin());
          }
        }

        if (scan_indices.empty())
        {
          continue;
        }

        // iterate over all pixels (m/z dimension)
        for (Size mz = 0; mz < mz_pixel_count; ++mz)
        {
          double mz_start = mz_min + mz_step_size * mz;
          double mz_end = mz_start + mz_step_size;

          // iterate over all relevant peaks in all relevant scans
          float max = -1.0;
          for (Size i = 0; i < scan_indices.size(); ++i)
          {
            Size s = scan_indices[i];
            Size p = peak_indices[i];
            const auto& spec = map[s];
            for (; p < spec.size(); ++p)
            {
              if (spec[p].getMZ() >= mz_end)
              {
                break;
              }
              if (spec[p].getIntensity() > max && layer_->filters.passes(spec, p))
              {
                max = spec[p].getIntensity();
              }
            }
            peak_indices[i] = p; // store last peak index for next m/z pixel
          }

          // draw to buffer
          if (max >= 0.0)
          {
            QPoint pos = canvas->dataToWidget_(mapper.mapToPoint(rt_start + 0.5 * rt_step_size, mz_start + 0.5 * mz_step_size));
            canvas->buffer_.setPixel(pos.x(), pos.y(), canvas->heightColor_(max, layer_->gradient, snap_factor).rgb());
          }
        }
      }
    }; // end lambda

    if (map.isIMFrame())
    {
      RT_or_IM_paint(DimInfoIM(map, area.getRangeForDim(MSDim::IM), canvas->unit_mapper_));
    }
    else
    {
      RT_or_IM_paint(DimInfoRT(map, area.getRangeForDim(MSDim::RT), canvas->unit_mapper_));
    }
  }

  void Painter2DPeak::paintPrecursorPeaks_(QPainter& painter, Plot2DCanvas* canvas)
  {
    const auto& peak_map = *layer_->getPeakData();

    QPen p;
    p.setColor(Qt::black);
    painter.setPen(p);

    auto it_prec = peak_map.end(); // MS1 spectrum of parent ion
    auto it_end = peak_map.RTEnd(canvas->visible_area_.getAreaUnit().getMaxRT());
    for (auto it = peak_map.RTBegin(canvas->visible_area_.getAreaUnit().getMinRT()); it != it_end; ++it)
    {
      // remember last precursor spectrum (do not call it->getPrecursorSpectrum(), since it can be very slow if no MS1 data is present)
      if (it->getMSLevel() == 1)
      {
        it_prec = it;
      }
      else if (it->getMSLevel() == 2 && !it->getPrecursors().empty())
      { // this is an MS/MS scan
        
        // position of precursor in MS2 (only works for 2D views with RT, m/z), not for ion mobility (IM, m/z) views.
        try
        {
          const auto data_xy_ms2 = canvas->unit_mapper_.map(Peak2D({it->getRT(), it->getPrecursors()[0].getMZ()}, {}));
          const QPoint pos_px_ms2 = canvas->dataToWidget_(data_xy_ms2); 
          const int x2 = pos_px_ms2.x();
          const int y2 = pos_px_ms2.y();

          if (it_prec != peak_map.end())
          {
            // position of precursor in MS1
            const auto data_xy_ms1 = canvas->unit_mapper_.map(Peak2D({it_prec->getRT(), it->getPrecursors()[0].getMZ()}, {}));
            const QPoint pos_px_ms1 = canvas->dataToWidget_(data_xy_ms1);
            const int x = pos_px_ms1.x();
            const int y = pos_px_ms1.y();
            // diamond shape in MS1
            drawDiamond({x, y}, &painter, 6);

            // rt position of corresponding MS2
            painter.drawLine(x, y, x2, y2);
          }
          else // no preceding MS1
          {
            // rt position of corresponding MS2 (cross)
            drawCross({x2, y2}, &painter, 6);
          }
        } // end try
        catch (...)
        {
          // paint nothing, since the coordinate system is wrong
        }
      }
    }
  }

  Painter2DChrom::Painter2DChrom(const LayerDataChrom* parent) : layer_(parent)
  {
  }

  void Painter2DChrom::paint(QPainter* painter, Plot2DCanvas* canvas, int /*layer_index*/)
  {
    const PeakMap& exp = *layer_->getChromatogramData();
    // TODO CHROM implement layer filters

    // paint chromatogram rt start and end as line
    float mz_origin = 0;

    QPoint posi;
    QPoint posi2;
    for (const auto& chrom : exp.getChromatograms())
    {
      if (chrom.empty())
      {
        continue;
      }
      mz_origin = chrom.getPrecursor().getMZ();
      posi = canvas->dataToWidget_(canvas->unit_mapper_.map(Peak2D {{chrom.front().getRT(), mz_origin}, 0}));
      posi2 = canvas->dataToWidget_(canvas->unit_mapper_.map(Peak2D {{chrom.back().getRT(), mz_origin}, 0}));
      painter->drawLine(posi.x(), posi.y(), posi2.x(), posi2.y());
    }
  }

  Painter2DIonMobility::Painter2DIonMobility(const LayerDataIonMobility* parent) : layer_(parent)
  {
  }

  void Painter2DIonMobility::paint(QPainter* /*painter*/, Plot2DCanvas* /*canvas*/, int /*layer_index*/)
  {
  }

  Painter2DFeature::Painter2DFeature(const LayerDataFeature* parent) : layer_(parent)
  {
  }

  void Painter2DFeature::paint(QPainter* painter, Plot2DCanvas* canvas, int layer_index)
  {
    if (canvas->getLayerFlag(layer_index, LayerDataBase::F_HULLS))
    {
      paintTraceConvexHulls_(painter, canvas);
    }
    if (canvas->getLayerFlag(layer_index, LayerDataBase::F_HULL))
    {
      paintFeatureConvexHulls_(painter, canvas);
    }
    if (canvas->getLayerFlag(layer_index, LayerDataBase::F_UNASSIGNED))
    {
      paintPeptideIDs_(painter, canvas, layer_->getPeptideIds(), layer_index);
    }


    const double snap_factor = canvas->snap_factors_[layer_index];

    int line_spacing = QFontMetrics(painter->font()).lineSpacing();
    const auto icon = toShapeIcon(layer_->param.getValue("dot:feature_icon").toString());
    Size icon_size = layer_->param.getValue("dot:feature_icon_size");
    bool show_label = (layer_->label != LayerDataBase::L_NONE);
    UInt num = 0;
    for (const auto& f : *layer_->getFeatureMap())
    {
      if (canvas->visible_area_.getAreaUnit().containsRT(f.getRT()) && canvas->visible_area_.getAreaUnit().containsMZ(f.getMZ()) && layer_->filters.passes(f))
      {
        // determine color
        QColor color;
        if (f.metaValueExists(5))
        {
          color = QColor(f.getMetaValue(5).toQString());
        }
        else
        {
          color = canvas->heightColor_(f.getIntensity(), layer_->gradient, snap_factor);
        }
        // paint
        QPoint pos = canvas->dataToWidget_(canvas->unit_mapper_.map(f));
        drawIcon(pos, color.rgb(), icon, icon_size, *painter);
        // labels
        if (show_label)
        {
          if (layer_->label == LayerDataBase::L_INDEX)
          {
            painter->setPen(Qt::darkBlue);
            painter->drawText(pos.x() + 10, pos.y() + 10, QString::number(num));
          }
          else if ((layer_->label == LayerDataBase::L_ID || layer_->label == LayerDataBase::L_ID_ALL) && !f.getPeptideIdentifications().empty() && !f.getPeptideIdentifications()[0].getHits().empty())
          {
            painter->setPen(Qt::darkGreen);
            Size maxHits = (layer_->label == LayerDataBase::L_ID_ALL) ? f.getPeptideIdentifications()[0].getHits().size() : 1;
            for (Size j = 0; j < maxHits; ++j)
            {
              painter->drawText(pos.x() + 10, pos.y() + 10 + int(j) * line_spacing, f.getPeptideIdentifications()[0].getHits()[j].getSequence().toString().toQString());
            }
          }
          else if (layer_->label == LayerDataBase::L_META_LABEL)
          {
            painter->setPen(Qt::darkBlue);
            painter->drawText(pos.x() + 10, pos.y() + 10, f.getMetaValue(3).toQString());
          }
        }
      }
      ++num;
    }
  }

  void Painter2DFeature::highlightElement(QPainter* painter, Plot2DCanvas* canvas, const PeakIndex /*element*/)
  {
    painter->setPen(QPen(Qt::red, 2));
    const Feature& f = canvas->selected_peak_.getFeature(*layer_->getFeatureMap());
    paintConvexHulls_(*painter, canvas, f.getConvexHulls(), f.getPeptideIdentifications().size() && f.getPeptideIdentifications()[0].getHits().size());
  }

  void Painter2DFeature::paintTraceConvexHulls_(QPainter* painter, Plot2DCanvas* canvas)
  {
    painter->setPen(Qt::black);
    const auto& area = canvas->visible_area_.getAreaUnit();
    for (const auto& f : *layer_->getFeatureMap())
    {
      if (area.containsRT(f.getRT()) && area.containsMZ(f.getMZ()) && layer_->filters.passes(f))
      {
        bool hasIdentifications = !f.getPeptideIdentifications().empty() && !f.getPeptideIdentifications()[0].getHits().empty();
        paintConvexHulls_(*painter, canvas, f.getConvexHulls(), hasIdentifications);
      }
    }
  }

  void Painter2DFeature::paintFeatureConvexHulls_(QPainter* painter, Plot2DCanvas* canvas)
  {
    const auto& area = canvas->visible_area_.getAreaUnit();
    for (const auto& f : *layer_->getFeatureMap())
    {
      if (area.containsRT(f.getRT()) && area.containsMZ(f.getMZ()) && layer_->filters.passes(f))
      {
        paintConvexHull_(*painter, canvas, f.getConvexHull(), !f.getPeptideIdentifications().empty() && !f.getPeptideIdentifications()[0].getHits().empty());
      }
    }
  }

  Painter2DConsensus::Painter2DConsensus(const LayerDataConsensus* parent) : layer_(parent)
  {
  }

  void Painter2DConsensus::paint(QPainter* painter, Plot2DCanvas* canvas, int layer_index)
  {
    if (canvas->getLayerFlag(layer_index, LayerDataBase::C_ELEMENTS))
    {
      paintConsensusElements_(painter, canvas, layer_index);
    }

    const double snap_factor = canvas->snap_factors_[layer_index];
    const auto icon = toShapeIcon(layer_->param.getValue("dot:feature_icon").toString());
    Size icon_size = layer_->param.getValue("dot:feature_icon_size");

    const auto area = canvas->visible_area_.getAreaUnit();
    for (const auto& cf : *layer_->getConsensusMap())
    {
      if (area.containsRT(cf.getRT()) && area.containsMZ(cf.getMZ()) && layer_->filters.passes(cf))
      {
        // determine color
        QColor color;
        if (cf.metaValueExists(5))
        {
          color = cf.getMetaValue(5).toQString();
        }
        else
        {
          // use intensity as color
          color = canvas->heightColor_(cf.getIntensity(), layer_->gradient, snap_factor);
        }

        // paint
        auto pos_unit = canvas->unit_mapper_.map(cf);
        drawIcon(canvas->dataToWidget_(pos_unit), color.rgb(), icon, icon_size, *painter);
      }
    }
  }

  void Painter2DConsensus::highlightElement(QPainter* painter, Plot2DCanvas* canvas, const PeakIndex element)
  {
    painter->setPen(QPen(Qt::red, 2));
    paintConsensusElement_(painter, canvas, canvas->getCurrentLayerIndex(), element.getFeature(*layer_->getConsensusMap()));
  }

  void Painter2DConsensus::paintConsensusElements_(QPainter* painter, Plot2DCanvas* canvas, Size layer_index)
  {
    for (const auto& cf : *layer_->getConsensusMap())
    {
      paintConsensusElement_(painter, canvas, layer_index, cf);
    }
  }

  void Painter2DConsensus::paintConsensusElement_(QPainter* painter, Plot2DCanvas* canvas, Size layer_index, const ConsensusFeature& cf)
  {
    // Is CF or any of its handles visible?
    if (!isConsensusFeatureVisible_(canvas, cf, layer_index) || !layer_->filters.passes(cf))
    {
      return;
    }

    // calculate position of consensus feature (centroid)
    QPoint consensus_pos = canvas->dataToWidget_(canvas->unit_mapper_.map(cf));
    // iterate over elements
    for (const FeatureHandle& element : cf)
    {
      // calculate position of consensus element
      QPoint pos = canvas->dataToWidget_(canvas->unit_mapper_.map(element));
      // paint line
      painter->drawLine(consensus_pos, pos);
      painter->drawPoint(pos.x(), pos.y());
      painter->drawPoint(pos.x() - 1, pos.y());
      painter->drawPoint(pos.x() + 1, pos.y());
      painter->drawPoint(pos.x(), pos.y() - 1);
      painter->drawPoint(pos.x(), pos.y() + 1);
    }
  }

  bool Painter2DConsensus::isConsensusFeatureVisible_(const Plot2DCanvas* canvas, const ConsensusFeature& cf, Size layer_index)
  {
    const auto& area = canvas->visible_area_.getAreaUnit();
    // check the centroid first
    if (area.containsRT(cf.getRT()) && area.containsMZ(cf.getMZ()))
    {
      return true;
    }

    // if element-flag is set, check if any of the consensus elements is visible
    if (canvas->getLayerFlag(layer_index, LayerDataBase::C_ELEMENTS))
    {
      for (const auto& ce : cf.getFeatures())
      {
        if (area.containsRT(ce.getRT()) && area.containsMZ(ce.getMZ()))
        {
          return true;
        }
      }
    }
    return false;
  }

  Painter2DIdent::Painter2DIdent(const LayerDataIdent* parent) : layer_(parent)
  {
  }

  void Painter2DIdent::paint(QPainter* painter, Plot2DCanvas* canvas, int layer_index)
  {
    paintPeptideIDs_(painter, canvas, layer_->getPeptideIds(), layer_index);
  }
} // namespace OpenMS
