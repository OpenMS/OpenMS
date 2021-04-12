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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

// declare Qt classes OUTSIDE of namespace OpenMS!
class QString; 
class QStringList;
class QPainter;
class QPoint;

#include <QColor>
#include <QCursor>
#include <QFont>

#include <array>

namespace OpenMS
{
  /**
    Namespace which holds static GUI-related helper functions.
  */
  namespace GUIHelpers
  {
    
    /// Open a folder in file explorer
    /// Will show a message box on failure
    void openFolder(const QString& folder);

    /// Open TOPPView (e.g. from within TOPPAS)
    void startTOPPView(const QStringList& args);

    /// Open a certain URL (in a browser)
    /// Will show a message box on failure
    void openURL(const QString& target);

    /**
       @brief draw a multi-line text at coordinates XY using a specific font and color

       Internally used getTextDimension() to figure out the size of the text-block/background which needs to be painted.

       @param painter Where to draw
       @param text Each item is a new line
       @param where Coordinates where to start drawing (upper left corner of text)
       @param col_fg Optional text color; if invalid (=default) will use the current painter's color
       @param col_bg Optional background color of bounding rectangle; if invalid (=default) no background will be painted
       @param Optional font; will use Courier by default
    */
    void drawText(QPainter& painter, const QStringList& text, const QPoint& where, const QColor col_fg = QColor("invalid"), const QColor col_bg = QColor("invalid"), const QFont& f = QFont("Courier"));


    /**
      @brief Obtains the bounding rectangle of a text (useful to determine overlaps etc)
    
    */
    QRectF getTextDimension(const QStringList& text, const QFont& font, int& line_spacing);


    /**
      @brief Given a set of levels (rows), try to add items at to topmost row which does not overlap an already placed item in this row (according to its x-coordinate)

      If a collision occurs, try the row below.
      If no row is collision-free, pick the one with the smallest overlap.

      X coordinates should always be positive.
    */
    class OverlapDetector
    {
    public:
      /// C'tor: number of @p levels must be >=1
      explicit OverlapDetector(int levels);

      /// try to put an item which spans from @p x_start to @p x_end in the topmost row possible
      /// @return the smallest row index (starting at 0) which has none (or the least) overlap
      size_t placeItem(double x_start, double x_end);

    private:
      std::vector<double> rows_;
    };

    /**
      @brief RAII class to disable the GUI and set a busy cursor and go back to the orignal state when this class is destroyed
    */
    class GUILock
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

    StringList convert(const QStringList& in);
    QStringList convert(const StringList& in);

  }; // GUIHelpers
}
