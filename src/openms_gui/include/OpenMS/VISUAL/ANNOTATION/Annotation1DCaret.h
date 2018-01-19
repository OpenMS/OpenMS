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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_ANNOTATION_ANNOTATION1DCARET_H
#define OPENMS_VISUAL_ANNOTATION_ANNOTATION1DCARET_H

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DItem.h>

#include <QtGui/QColor>
#include <QStaticText> 
#include <QTextDocument>

#include <vector>

namespace OpenMS
{
  /** @brief An annotation item which paints a set of carets on the canvas.

      Most useful to visualize (theoretical) isotope distributions (one caret per isotope position).
      Additionally, a text annotation can be provided.

      @see Annotation1DItem
  */
  class Annotation1DCaret :
    public Annotation1DItem
  {
public:

    typedef Annotation1DItem::PointType PointType;
    typedef std::vector<PointType> PositionsType;

    /// Constructor
    Annotation1DCaret(const PositionsType& poly_positions, const QString& text, const QColor& colour, const QColor& connection_line_color);

    /// Copy constructor
    Annotation1DCaret(const Annotation1DCaret& rhs);

    /// Destructor
    ~Annotation1DCaret() override;

    // Docu in base class
    void ensureWithinDataRange(Spectrum1DCanvas* const canvas) override;

    // Docu in base class
    void draw(Spectrum1DCanvas* const canvas, QPainter& painter, bool flipped = false) override;

    // Docu in base class
    void move(const PointType& delta) override;

    /// Returns the positions of the lines (in MZ / intensity coordinates)
    const PositionsType& getCaretPositions() const;

    /// Sets the position of the label (in MZ / intensity coordinates)
    void setPosition(const PointType& position);

    /// Returns the position of the annotated peak (in MZ / intensity coordinates)
    const PointType& getPosition() const;

    /// Set the colour of the carets (colour of text must be set using html)
    void setColor(const QColor& color);

    /// Returns the colour of the carets
    const QColor& getColor() const;

    /// The text to display (optional).
    /// Rendered using QStaticText, so HTML formatting is allowed.
    void setRichText(const QString& text);

protected:

    /// The positions of points (in MZ/intensity coordinates)
    /// Ensure positions are sorted by m/z when assigning
    PositionsType caret_positions_;

    /// The position of the label (in MZ/intensity coordinates)
    PointType position_;

    /// The colour of the label
    QColor color_;

    /// The colour of the (optional) dashed line connecting peak and label
    QColor connection_line_color_;

    /// Holds the (rich) text
    QStaticText st_;

  };
} // namespace OpenMS

#endif
