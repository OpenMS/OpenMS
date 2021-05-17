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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Define SeqAn version.
// ==========================================================================

#ifndef SEQAN_VERSION_H_
#define SEQAN_VERSION_H_

/**
.Macro.SEQAN_VERSION_MAJOR
..cat:Versioning
..summary:Major SeqAn revision number.
..signature:SEQAN_VERSION_MAJOR
..example:For SeqAn version "1.3", this value is $1$, for "2.5.4", it is $2$.
..include:seqan/version.h

.Macro.SEQAN_VERSION_MINOR
..cat:Versioning
..summary:Minor SeqAn revision number.
..signature:SEQAN_VERSION_MINOR
..example:For SeqAn version "1.3", this value is $3$, for "1.5.4", it is $5$.
..include:seqan/version.h

.Macro.SEQAN_VERSION_PATCH
..cat:Versioning
..summary:SeqAn patch revision number.
..signature:SEQAN_VERSION_PATCH
..example:For SeqAn version "1.3", this value is $0$, for "1.3.4", it is $4$.
..include:seqan/version.h

.Macro.SEQAN_VERSION_PRE_RELEASE
..cat:Versioning
..summary:Flag ($0$/$1$) to indicate whether this is a pre-release (i.e. SVN version).
..signature:SEQAN_VERSION_PRE_RELEASE
..include:seqan/version.h
*/ 

#define SEQAN_VERSION_MAJOR 1

#define SEQAN_VERSION_MINOR 4

#define SEQAN_VERSION_PATCH 1

#define SEQAN_VERSION_PRE_RELEASE 0

#endif  // SEQAN_VERSION_H_
