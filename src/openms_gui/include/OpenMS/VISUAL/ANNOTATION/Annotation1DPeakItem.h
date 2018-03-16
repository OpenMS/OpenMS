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

#ifndef OPENMS_VISUAL_ANNOTATION_ANNOTATION1DPEAKITEM_H
#define OPENMS_VISUAL_ANNOTATION_ANNOTATION1DPEAKITEM_H

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DItem.h>

#include <QtGui/QColor>

namespace OpenMS
{
  /** @brief A peak annotation item
            @see Annotation1DItem
    */
  class Annotation1DPeakItem :
    public Annotation1DItem
  {
public:
    /// Constructor
    Annotation1DPeakItem(const PointType& peak_position, const QString& text, const QColor& color);

    /// Copy constructor
    Annotation1DPeakItem(const Annotation1DPeakItem& rhs);

    /// Destructor
    ~Annotation1DPeakItem() override;

    /// Docu in base class
    void ensureWithinDataRange(Spectrum1DCanvas* const canvas) override;

    /// Docu in base class
    void draw(Spectrum1DCanvas* const canvas, QPainter& painter, bool flipped = false) override;

    /// Docu in base class
    void move(const PointType& /*delta*/) override;

    /// Returns the position of the label (peak) (in MZ/intensity coordinates)
    const PointType& getPeakPosition() const;

    /// Sets the position of the label (in MZ/intensity coordinates)
    void setPosition(const PointType& position);

    /// Returns the position of the annotated peak (in MZ/intensity coordinates)
    const PointType& getPosition() const;

    /// Set the color of the label
    void setColor(const QColor& color);

    /// Returns the color of the label
    const QColor& getColor() const;
protected:
    /// The position of the anchor (peak) (in MZ / intensity coordinates)
    PointType peak_position_;

    /// The position of the label (in MZ / intensity coordinates)
    PointType position_;

    /// The color of the label
    QColor color_;
  };
} // namespace OpenMS

#endif
