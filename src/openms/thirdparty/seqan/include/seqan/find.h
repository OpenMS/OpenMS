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

#ifndef SEQAN_HEADER_FIND_H
#define SEQAN_HEADER_FIND_H

//____________________________________________________________________________
// prerequisites

#include <cmath>

#include <deque>

#include <seqan/sequence.h>
#include <seqan/modifier.h>
#include <seqan/score.h>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <seqan/map.h>

//____________________________________________________________________________

#include <seqan/find/find_base.h>
#include <seqan/find/find_pattern_base.h>

//____________________________________________________________________________
// exact pattern matching
#include <seqan/find/find_simple.h>
#include <seqan/find/find_horspool.h>
#include <seqan/find/find_shiftand.h>
#include <seqan/find/find_shiftor.h>
#include <seqan/find/find_bndm.h>
#include <seqan/find/find_bom.h>

//____________________________________________________________________________
//complex pattern matching
#include <seqan/find/find_wild_shiftand.h>

//____________________________________________________________________________
//multiple pattern search
#include <seqan/find/find_ahocorasick.h>
#include <seqan/find/find_multiple_shiftand.h>
#include <seqan/find/find_set_horspool.h>

//#include <seqan/find/find_multi.h> //wegwerfen
#include <seqan/find/find_wumanber.h> 
#include <seqan/find/find_multiple_bfam.h>

//____________________________________________________________________________
// approximate pattern matching
#include <seqan/find/find_begin.h>

#include <seqan/find/find_score.h>
#include <seqan/find/find_myers_ukkonen.h>
#include <seqan/find/find_abndm.h>
#include <seqan/find/find_pex.h>

#include <seqan/find/find_hamming_simple.h>

//____________________________________________________________________________

#endif //#ifndef SEQAN_HEADER_...
