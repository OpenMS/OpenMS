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
// $Maintainer: Hendrik Weisser $
// $Authors: Johannes Junker, Chris Bielow, Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_TOPPASSPLITTERVERTEX_H
#define OPENMS_VISUAL_TOPPASSPLITTERVERTEX_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/VISUAL/TOPPASVertex.h>

namespace OpenMS
{
  /**
      @brief A special vertex that allows to split a list of inputs.

      Tools that produce lists of output files (several files in each processing round, e.g. map alignment tools) cannot directly provide input for tools that only take a single input file in TOPPAS. This "Splitter" node provides the necessary glue, by splitting a list of output files into several rounds of single input files.

      See the @ref OpenMS::TOPPASMergerVertex "Collector" node for the opposite operation.

      @ingroup TOPPAS_elements
  */
  class OPENMS_GUI_DLLAPI TOPPASSplitterVertex :
    public TOPPASVertex
  {
    Q_OBJECT

public:

    /// Default constructor
    TOPPASSplitterVertex();
    /// Copy constructor
    TOPPASSplitterVertex(const TOPPASSplitterVertex& rhs);
    /// Destructor
    ~TOPPASSplitterVertex() override;
    /// Assignment operator
    TOPPASSplitterVertex& operator=(const TOPPASSplitterVertex& rhs);
    /// returns "SplitterVertex"
    String getName() const override;
    /// check if upstream nodes are finished and call downstream nodes
    void run() override;
    // documented in base class
    void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget) override;
    // documented in base class
    QRectF boundingRect() const override;
    // documented in base class
    QPainterPath shape() const override;
    // documented in base class
    void markUnreachable() override;

protected:

    ///@name reimplemented Qt events
    //@{
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* e) override;
    //@}

  };
}

#endif // OPENMS_VISUAL_TOPPASSPLITTERVERTEX_H
