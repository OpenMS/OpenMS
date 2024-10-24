// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/FileHandler.h>

// declare Qt classes OUTSIDE of namespace OpenMS!
class QPainter;
class QPoint;
class QPointF;
class QRectF;
class QWidget;

#include <QColor>
#include <QFont>
#include <QtCore/qcontainerfwd.h> // for QStringList

#include <array>

namespace OpenMS
{

  class FileTypeList;

  /**
    Namespace which holds static GUI-related helper functions.
  */
  namespace GUIHelpers
  {        
    /// Open a folder in file explorer
    /// Will show a message box on failure
    OPENMS_GUI_DLLAPI void openFolder(const QString& folder);

    /// Open a dialog to select a filename to save data to.

    OPENMS_GUI_DLLAPI QString getSaveFilename(QWidget* parent,
                                              const QString& caption,
                                              const QString& dir,
                                              const FileTypeList& supported_file_types, 
                                              bool add_all_filter,
                                              const FileTypes::Type fallback_extension);


    /// Open TOPPView (e.g. from within TOPPAS) as a detached process (i.e. will continue running when this process ends)
    /// @return true if process started successfully
    OPENMS_GUI_DLLAPI bool startTOPPView(QStringList args);

    /// Open a certain URL (in a browser)
    /// Will show a message box on failure
    OPENMS_GUI_DLLAPI void openURL(const QString& target);

    /**
       @brief draw a multi-line text at coordinates XY using a specific font and color

       Internally used getTextDimension() to figure out the size of the text-block/background which needs to be painted.

       @param painter Where to draw
       @param text Each item is a new line
       @param where Coordinates where to start drawing (upper left corner of text)
       @param col_fg Optional text color; if invalid (=default) will use the current painter's color
       @param col_bg Optional background color of bounding rectangle; if invalid (=default) no background will be painted
       @param font Font to use; will use Courier by default
    */
    OPENMS_GUI_DLLAPI void drawText(QPainter& painter, const QStringList& text, const QPoint& where, const QColor& col_fg = QColor("invalid"), const QColor& col_bg = QColor("invalid"),
                                   const QFont& font = QFont("Courier"));


    /**
      @brief Obtains the bounding rectangle of a text (useful to determine overlaps etc)
    
    */
    OPENMS_GUI_DLLAPI QRectF getTextDimension(const QStringList& text, const QFont& font, int& line_spacing);


    /// Returns the point in the @p list that is nearest to @p origin
    OPENMS_GUI_DLLAPI QPointF nearestPoint(const QPointF& origin, const QList<QPointF>& list);

    /**
     * \brief Find the point on a rectangle where a ray/line from a point @p p to its center would intersect at
     * \param rect Rectangle which intersects with the line from @p p to its center
     * \param p A point outside the rectangle
     * \return The intersection point or the center() of @p rect if @p p is inside the rectangle
     */
    OPENMS_GUI_DLLAPI QPointF intersectionPoint(const QRectF& rect, const QPointF& p);


    /**
      @brief A heuristic: Given a set of levels (rows), try to add items at to topmost row which does not overlap an already placed item in this row (according to its x-coordinate)

      If a collision occurs, try the row below.
      If no row is collision-free, pick the one with the smallest overlap.

      Only positions beyond the largest x-coordinate observed so far are considered (i.e. gaps are not filled).

      X-coordinates should always be positive (a warning is issued otherwise).
    */
    class OPENMS_GUI_DLLAPI OverlapDetector
    {
    public:
      /// C'tor: number of @p levels must be >=1
      /// @throw Exception::InvalidSize if levels <= 0
      explicit OverlapDetector(int levels);

      /// try to put an item which spans from @p x_start to @p x_end in the topmost row possible
      /// @return the smallest row index (starting at 0) which has none (or the least) overlap
      size_t placeItem(double x_start, double x_end);

    private:
      std::vector<double> rows_; ///< store the largest x_end for each row
    };

    /**
      @brief RAII class to disable the GUI and set a busy cursor and go back to the original state when this class is destroyed
    */
    class OPENMS_GUI_DLLAPI GUILock
    {
    public:
      /// C'tor receives the widget to lock
      /// @param gui QWidget to lock(including all children); can be nullptr (nothing will be locked)
      GUILock(QWidget* gui);

      /// no copy/assignment allowed
      GUILock(const GUILock& rhs) = delete;
      GUILock(GUILock&& rhs) = delete;
      GUILock& operator=(const GUILock& rhs) = delete;

      /// D'tor: unlocks the GUI (does nothing if already unlocked)
      ~GUILock();
      
      /// manually lock the GUI (does nothing if already locked)
      void lock();
      /// manually unlock the GUI (does nothing if already unlocked)
      void unlock();

    private:
      QWidget* locked_widget_{ nullptr };
      bool currently_locked_{ false };
      bool was_enabled_{ true };
    };

    /// color palette for certain purposes
    /// Currently, only a set of distinct colors is supported
    class ColorBrewer
    {
    public:
      struct Distinct
      {
        enum NAMES
        {
          Red,
          Blue,
          Green,
          Brown,
          Purple,
          LightGrey,
          LightGreen,
          LightBlue,
          Cyan,
          Orange,
          Yellow,
          Tan,
          Pink,
          DarkGrey,
          SIZE_OF_NAMES
        };

        const std::array<QColor, NAMES::SIZE_OF_NAMES> values = { { Qt::red,
                                                                    Qt::blue,
                                                                    Qt::green,
                                                                    QColor(129, 74, 25) /*brown*/,
                                                                    QColor(129, 38, 192) /*purple*/,
                                                                    Qt::lightGray,
                                                                    QColor(129,197,122) /*lightGreen*/,
                                                                    QColor(157,175,255) /*lightBlue*/,
                                                                    Qt::cyan,
                                                                    QColor(255,146,51) /*orange*/,
                                                                    Qt::yellow,
                                                                    QColor(233,222,187) /*tan*/,
                                                                    QColor(255,205,243) /*pink*/,
                                                                    Qt::darkGray } };
      };

      /// get a certain color. If @p index is larger than the maximum color, modulo operator will applied (cycling through colors)
      template<class COLOR_CLASS>
      static QColor getColor(uint32_t index)
      {
        // cycle if necessary
        if (index >= COLOR_CLASS::NAMES::SIZE_OF_NAMES) index = index % COLOR_CLASS::NAMES::SIZE_OF_NAMES;
        return COLOR_CLASS().values[index];
      }
    }; // ColorBrewer

    OPENMS_GUI_DLLAPI StringList convert(const QStringList& in);
    OPENMS_GUI_DLLAPI QStringList convert(const StringList& in);

  }; // GUIHelpers
}
