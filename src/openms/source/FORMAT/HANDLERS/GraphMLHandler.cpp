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
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/GraphMLHandler.h>
#include <ostream>

using namespace std;

namespace OpenMS
{
  namespace Internal
  {
    GraphMLHandler::GraphMLHandler(std::vector<Size>& nodes, std::map<std::pair<Size, Size>, Size>& edges, const String &filename)
    :
    XMLHandler(filename, "2.0"),
    nodes_(nodes),
    edges_(edges),
    filepath_(filename)
    {
    }

    void GraphMLHandler::writeTo(std::ostream & os)
    {
      // write header
      writeYFileHeader(os);

      // write keys
      os << "  <key for=\"node\" id=\"Ng\" yfiles.type=\"nodegraphics\"/>\n" // for text on nodes
         << "  <key for=\"edge\" id=\"Eg\" yfiles.type=\"edgegraphics\"/>\n"; // for edge style (undirected)

      // writing graph
      os << "  <graph edgedefault=\"undirected\" id=\"G\">\n";

      // writing nodes
      for (auto& node_index : nodes_)
      {
        os << "    <node id=\"n" << to_string(node_index) << "\">\n"
              "      <data key=\"Ng\">\n"
              "        <y:ShapeNode>\n"
              "          <y:NodeLabel>" << to_string(node_index) << "</y:NodeLabel>\n"
              "        </y:ShapeNode>\n"
              "      </data>\n"
              "    </node>\n";
      }

      // writing edges
      Size counter = 0;
      for (auto& edge_pair : edges_)
      {
        os << "    <edge id=\"e" << to_string(counter) << "\" source=\"n"<< edge_pair.first.first <<"\" target=\"n" << edge_pair.first.second <<"\">\n"
              "      <data key=\"Eg\">\n"
              "        <y:PolyLineEdge>\n"
              "          <y:Arrows source=\"none\" target=\"none\"/>\n"
//              "          <y:EdgeLabel>" << edge_pair.second << "</y:EdgeLabel>"
              "        </y:PolyLineEdge>\n"
              "      </data>\n"
              "    </edge>\n";
        ++counter;
      }

      os << "  </graph>";

      // ending clause for header
      os << "</graphml>\n" ;
    }

    void GraphMLHandler::writeYFileHeader(std::ostream & os) const
    {
      // based on yEd
      os << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n"
            "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\" xmlns:java=\"http://www.yworks.com/xml/yfiles-common/1.0/java\" "
            "xmlns:sys=\"http://www.yworks.com/xml/yfiles-common/markup/primitives/2.0\" xmlns:x=\"http://www.yworks.com/xml/yfiles-common/markup/2.0\" "
            "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns:y=\"http://www.yworks.com/xml/graphml\" xmlns:yed=\"http://www.yworks.com/xml/yed/3\" "
            "xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns http://www.yworks.com/xml/schema/graphml/1.1/ygraphml.xsd\">" << "\n";
    }
  }
}