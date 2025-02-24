// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

//OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/MATH/STATISTICS/Histogram.h>
#include <OpenMS/VISUAL/EnhancedTabBarWidgetInterface.h>
#include <OpenMS/VISUAL/PlotCanvas.h>

class QCloseEvent;
class QGridLayout;
class QMimeData;
class QScrollBar;

namespace OpenMS
{

  class AxisWidget;

  /**
      @brief Base class for spectrum widgets

      This class is the base class for the different MDI window types in the
      TOPPView application. For each type of spectrum view (such as 1D view, 2D
      view, 3D view etc.), there must exist a corresponding class derived from
      this class.

      In TOPPView, each PlotWidget holds an enclosed PlotCanvas with
      which it is paired (e.g. a Plot1DWidget holds a Plot1DCanvas)
      which can retrieved with the canvas() function. While the PlotCanvas
      does the actual drawing, the PlotWidget holds information about the
      axes (axis widgets), scrolling (scrollbars) etc. The PlotWidget uses
      a grid layout (QGridLayout) with a default 3x3 grid where the upper right
      corner of the grid is the canvas and the spaces left and below the canvas
      are for the axes widget and scrollbars.

      To integrate a new spectrum view (i.e. classes derived from
      PlotWidget and PlotCanvas) into the TOPPView application,
      a class must be derived from this class which holds an
      instance of the PlotCanvas class as a child widget.

      @todo Add support to store the displayed data as SVG image (HiWi)
  */
  class OPENMS_GUI_DLLAPI PlotWidget :
    public QWidget,
    public EnhancedTabBarWidgetInterface
  {
    Q_OBJECT

public:
    static const char RT_AXIS_TITLE[];
    static const char MZ_AXIS_TITLE[];
    static const char INTENSITY_AXIS_TITLE[];
    static const char IM_MS_AXIS_TITLE[];
    static const char IM_ONEKZERO_AXIS_TITLE[];

    /** @name Type definitions */
    //@{

    /// Main data type (experiment)
    typedef LayerDataBase::ExperimentType ExperimentType;
    /// Main data type (features)
    typedef LayerDataBase::FeatureMapType FeatureMapType;
    /// Spectrum type
    typedef ExperimentType::SpectrumType SpectrumType;
    //@}

    /// Default constructor
    PlotWidget(const Param & preferences, QWidget * parent = nullptr);
    /// Destructor
    ~PlotWidget() override;

    /**
        @brief Returns a pointer to canvas object

        This method is overwritten for 1D, 2D, 3D to make the class specific members accessible.
        
        The canvas object is set with the setCanvas_() method.
        This is usually done in the constructor.
    */
    virtual PlotCanvas* canvas() const = 0;

    /// Returns a pointer to the x-axis axis widget.
    virtual inline AxisWidget * xAxis()
    {
      return x_axis_;
    }

    /// Returns a pointer to the y-axis axis widget.
    virtual inline AxisWidget * yAxis()
    {
      return y_axis_;
    }

    /// Get the mouse action mode
    Int getActionMode() const;

    /// Returns if the axis labels are shown
    virtual bool isLegendShown() const;

    /// Shows/hides axis labels
    virtual void showLegend(bool show);

    /// Sets the intensity mode of the PlotCanvas
    void setIntensityMode(PlotCanvas::IntensityModes mode);

    /// Hides x-axis and y-axis
    virtual void hideAxes();

    /// Saves the widget's content as image file
    virtual void saveAsImage();

signals:
    /// Emits a status message that should be displayed for @p time ms. If @p time is 0 the message should be displayed until the next message is emitted.
    void sendStatusMessage(std::string, OpenMS::UInt);
    /// Emitted when the cursor position changes (for displaying e.g. in status bar)
    void sendCursorStatus(const String& x_value, const String& y_value);
    /// Message about the destruction of this widget
    void aboutToBeDestroyed(int window_id);
    /// Shows the main preferences dialog
    void openPreferences();
    /// Signal that is emitted, when a drag-and-drop action ends on this widget
    void dropReceived(const QMimeData* data, QWidget* source, int id);

public slots:
    /// Shows statistics about the data (count, min, max, avg of intensity, charge, quality and meta data)
    void showStatistics();
    /// Shows the intensity distribution of the current layer
    void showIntensityDistribution(const Math::Histogram<>& dist);
    /// Shows the meta data distribution of value @p name of the current layer
    void showMetaDistribution(const String& name, const Math::Histogram<>& dist);
    /// Updates the axes by setting the right labels and calling recalculateAxes_();
    void updateAxes();
    /**
        @brief Updates the horizontal scrollbar

        @param min The overall minimum of the range
        @param disp_min The displayed minimum
        @param disp_max The displayed maximum
        @param max The overall maximum of the range
    */
    void updateHScrollbar(float min, float disp_min, float disp_max, float max);
    /**
        @brief Updates the vertical scrollbar

        @param min The overall minimum of the range
        @param disp_min The displayed minimum
        @param disp_max The displayed maximum
        @param max The overall maximum of the range
    */
    void updateVScrollbar(float min, float disp_min, float disp_max, float max);
    /// Shows a goto dialog
    virtual void showGoToDialog() = 0;
    /// Toggles the axis legend visibility
    void changeLegendVisibility();

    /**
     * \brief Set a new mapper for the canvas and axis. Internally, all dependent components are updated (e.g. projections in 2D View)
     * \param mapper The new mapper for translating between units and axis
     */
    virtual void setMapper(const DimMapper<2>& mapper) = 0;

protected:
    /// @name Reimplemented Qt events
    //@{
    void closeEvent(QCloseEvent * e) override;
    //@}

    /**
        @brief Adds the canvas, axes and scrollbars to the layout

        @p row and @p col define the position of the canvas.
        Axes and scrollbars are added to the left and bottom of the canvas.
    */
    void setCanvas_(PlotCanvas * canvas, UInt row = 0, UInt col = 2);
    /// Switch between different intensity modes
    virtual void intensityModeChange_();
    /// recalculates the Axis ticks
    virtual void recalculateAxes_() = 0;

    ///@name reimplemented Qt events
    //@{
    void dragEnterEvent(QDragEnterEvent * event) override;
    void dragMoveEvent(QDragMoveEvent * event) override;
    void dropEvent(QDropEvent * event) override;
    /// make our subclassed QWidget listen to things like stylesheet changes
    void paintEvent(QPaintEvent * /*event*/) override;
    //@}

    /// Pointer to the canvas widget
    PlotCanvas* canvas_;
    /// Main layout
    QGridLayout* grid_;
    /// Vertical axis
    AxisWidget* y_axis_;
    /// Horizontal axis
    AxisWidget* x_axis_;
    /// Horizontal scrollbar
    QScrollBar* x_scrollbar_;
    /// Vertical scrollbar
    QScrollBar* y_scrollbar_;
  };
}

