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
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_ANNOTATION_ANNOTATIONS1DCONTAINER_H
#define OPENMS_VISUAL_ANNOTATION_ANNOTATIONS1DCONTAINER_H

#include <OpenMS/DATASTRUCTURES/DPosition.h>

#include <list>

#include <QtGui/QPen>

class QPoint;
class QObject;
class QRectF;
class QPainter;

namespace OpenMS
{
  class Annotation1DItem;

  /// Container for annotations to content of Spectrum1DCanvas
  class Annotations1DContainer :
    public std::list<Annotation1DItem *>
  {
public:
    /// Default constructor
    Annotations1DContainer();

    /// Copy constructor
    Annotations1DContainer(const Annotations1DContainer & rhs);

    /// Assignment operator
    Annotations1DContainer & operator=(const Annotations1DContainer & rhs);

    /// Destructor
    virtual ~Annotations1DContainer();

    /// Iterator for the 1D annotations
    typedef std::list<Annotation1DItem *>::iterator Iterator;

    /// Const iterator for the 1D annotations
    typedef std::list<Annotation1DItem *>::const_iterator ConstIterator;

    /// Type of the Points
    typedef DPosition<2> PointType;

    /// Coordinate type
    typedef double CoordinateType;

    /** @brief Returns a pointer to the item at @p pos, or 0, if not existent

            If more than one item's bounding box encloses @p pos , the one in the
            foreground is returned.
    */
    Annotation1DItem * getItemAt(const QPoint & pos) const;

    /// Selects the item at @p pos on the canvas, if it exists.
    void selectItemAt(const QPoint & pos);

    /// Deselects the item at @p pos on the canvas, if it exists.
    void deselectItemAt(const QPoint & pos);

    /// Selects all items
    void selectAll();

    /// Deselects all items
    void deselectAll();

    /// Removes the selected items
    void removeSelectedItems();

    /// Returns the selected items
    std::vector<Annotation1DItem*> getSelectedItems();

    /// Sets the pen_
    void setPen(const QPen & pen);

    /// Returns the pen_
    const QPen & getPen() const;

    /// Sets the selected_pen_
    void setSelectedPen(const QPen & pen);

    /// Returns the selected_pen_
    const QPen & getSelectedPen() const;

    /// The pen used to draw items
    QPen pen_;

    /// The pen used to draw selected items
    QPen selected_pen_;
  };

} // namespace

#endif
