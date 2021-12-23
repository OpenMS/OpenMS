// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_INDEX_ESA_DRAWING_H
#define SEQAN_HEADER_INDEX_ESA_DRAWING_H

namespace SEQAN_NAMESPACE_MAIN
{

template <typename TFile, typename TText, typename TESASpec>
void write(TFile & file, 
	   Index<TText, IndexEsa<TESASpec> > & stree,
	   DotDrawing) 
{
//IOREV _nodoc_
	SEQAN_CHECKPOINT
	typedef Index<TText, IndexEsa<TESASpec> > TIndex;
	
	streamPut(file, "digraph G {\n");
	streamPut(file, '\n');
	streamPut(file, "/* Graph Attributes */\n");
	streamPut(file, "graph [rankdir = LR];\n");
	streamPut(file, '\n');
	streamPut(file, "/* Node Attributes */\n");
	streamPut(file, "node [shape = ellipse, fillcolor = lightgrey, style = filled, fontname = \"Times-Italic\"];\n");
	streamPut(file, '\n');
	streamPut(file, "/* Edge Attributes */\n");
	streamPut(file, "edge [fontname = \"Times-Italic\", arrowsize = 0.75, fontsize = 16];\n");
	streamPut(file, '\n');

	streamPut(file, "/* Edges */\n");
	typedef typename Iterator<TIndex, TopDown<ParentLinks<Preorder> > >::Type TIterator;
	typedef typename Iterator<TIndex, TopDown<> >::Type TIteratorSimple;
	TIterator it(stree);

	for(;!atEnd(it);++it) 
	{
		// dump node
       		streamPut(file, "\"[");
 		streamPut(file, (int)value(it).range.i1);
		streamPut(file, ':');
		streamPut(file, (int)value(it).range.i2);
       		streamPut(file, ")\"");
       		if (!isRightTerminal(it))
			streamPut(file, " [style = dashed]");
       		streamPut(file, ";\n");

		// dump edge from parent (if not root)
		if (!isRoot(it)) {
			TIteratorSimple src(container(it), nodeUp(it));

			streamPut(file, "\"[");
			streamPut(file, (int)value(src).range.i1);
			streamPut(file, ':');
			streamPut(file, (int)value(src).range.i2);
			streamPut(file, ")\"");

			streamPut(file, " -> ");

			streamPut(file, "\"[");
			streamPut(file, (int)value(it).range.i1);
			streamPut(file, ':');
			streamPut(file, (int)value(it).range.i2);
			streamPut(file, ")\"");

			streamPut(file, " [label = \"");
			streamPut(file, parentEdgeLabel(it));
			streamPut(file, "\"];\n");
		}
	}
	streamPut(file, '\n');

	streamPut(file, "}\n");
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
