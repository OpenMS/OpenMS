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
// Author: Konrad Ludwig Moritz Rudolph <konrad.rudolph@fu-berlin.de>
// ==========================================================================

//SEQAN_xNO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_HEADER_PIZZACHILI_API_H
#define SEQAN_HEADER_PIZZACHILI_API_H

#include <seqan/basic.h>

namespace SEQAN_NAMESPACE_MAIN {

namespace impl {
    typedef unsigned char uchar_t;
    typedef unsigned long ulong_t;
    typedef void* index_t;
    typedef int error_t;
} // namespace impl

struct InvalidPizzaChiliSpec;

template <typename TSpec>
struct PizzaChiliCodeProvider {
    typedef InvalidPizzaChiliSpec Type;
};

/**
.Tag.Pizza & Chili Index Tags
..summary:Tag specifying the Pizza & Chili library to use.
..remarks:More information for all the index libraries can be found in the
@http://pizzachili.dcc.uchile.cl|original documentation@ (or the
@http://pizzachili.di.unipi.it|Italian mirror@).
..cat:Index
..tag.PizzaChiliAF:The alphabet-friendly FM index.
..tag.PizzaChiliCcsa:The compressed compact suffix array index.
..tag.PizzaChiliFM: The FM (full-text in minute space) index.
..tag.PizzaChiili_RSA:The repair suffix array index.
...remarks:The index cannot be saved and loaded.
..tag.PizzaChiliSA: The simple suffix array index.
...remarks:The index cannot be saved and loaded.
..tag.PizzaChiliSada: the compressed suffix array index.
...remarks:The index cannot be saved and loaded.
..see:Spec.Pizza & Chili Index
..see:Spec.Pizza & Chili String
..include:seqan/index.h
*/

// We need to declare these explicitly instead through macro expansion in order
// for them to be included in the forward generated declarations.

struct PizzaChiliAF_;
typedef Tag<PizzaChiliAF_> const PizzaChiliAF;

struct PizzaChiliCcsa_;
typedef Tag<PizzaChiliCcsa_> const PizzaChiliCcsa;

struct PizzaChiliFM_;
typedef Tag<PizzaChiliFM_> const PizzaChiliFM;

struct PizzaChiliLZ_;
typedef Tag<PizzaChiliLZ_> const PizzaChiliLZ;

struct PizzaChiliRsa_;
typedef Tag<PizzaChiliRsa_> const PizzaChiliRsa;

struct PizzaChiliRlfm_;
typedef Tag<PizzaChiliRlfm_> const PizzaChiliRlfm;

struct PizzaChiliSA_;
typedef Tag<PizzaChiliSA_> const PizzaChiliSA;

struct PizzaChiliSada_;
typedef Tag<PizzaChiliSada_> const PizzaChiliSada;

struct PizzaChiliSsa_;
typedef Tag<PizzaChiliSsa_> const PizzaChiliSsa;

struct PizzaChiliTest_;
typedef Tag<PizzaChiliTest_> const PizzaChiliTest;

#define SEQAN_MAKE_PIZZACHILI_PROVIDER(name) \
    class PizzaChiliApi##name { \
    public: \
        static char* error_index(impl::error_t e); \
        static int build_index( \
            impl::uchar_t* text, \
            impl::ulong_t length, \
            char* build_options, \
            impl::index_t* index \
            ); \
        static int save_index(impl::index_t index, char* filename); \
        static int load_index(char* filename, impl::index_t* index); \
        static int free_index(impl::index_t index); \
        static int index_size(impl::index_t index, impl::ulong_t* size); \
        static int count( \
            impl::index_t index, \
            impl::uchar_t* pattern, \
            impl::ulong_t length, \
            impl::ulong_t* numocc \
        ); \
        static int locate( \
            impl::index_t index, \
            impl::uchar_t* pattern, \
            impl::ulong_t length, \
            impl::ulong_t** occ, \
            impl::ulong_t* numocc \
        ); \
        static int get_length(impl::index_t index, impl::ulong_t* length); \
        static int extract( \
            impl::index_t index, \
            impl::ulong_t from, \
            impl::ulong_t to, \
            impl::uchar_t** snippet, \
            impl::ulong_t* snippet_length \
        ); \
        static int display( \
            impl::index_t index, \
            impl::uchar_t* pattern, \
            impl::ulong_t length, \
            impl::ulong_t numc,  \
            impl::ulong_t* numocc, \
            impl::uchar_t** snippet_text, \
            impl::ulong_t** snippet_length \
        ); \
        static int init_ds_ssort(int adist, int bs_ratio); \
    }; \
    \
    /*struct PizzaChili##name_; \
    typedef Tag<PizzaChili##name_> const PizzaChili##name;*/ \
    \
    template <> \
    struct PizzaChiliCodeProvider<PizzaChili##name> { \
        typedef PizzaChiliApi##name Type; \
    };

SEQAN_MAKE_PIZZACHILI_PROVIDER(AF)
SEQAN_MAKE_PIZZACHILI_PROVIDER(Ccsa)
SEQAN_MAKE_PIZZACHILI_PROVIDER(FM)
SEQAN_MAKE_PIZZACHILI_PROVIDER(LZ)
SEQAN_MAKE_PIZZACHILI_PROVIDER(Rsa)
SEQAN_MAKE_PIZZACHILI_PROVIDER(Rlfm)
SEQAN_MAKE_PIZZACHILI_PROVIDER(SA)
SEQAN_MAKE_PIZZACHILI_PROVIDER(Sada)
SEQAN_MAKE_PIZZACHILI_PROVIDER(Ssa)
SEQAN_MAKE_PIZZACHILI_PROVIDER(Test)

#undef SEQAN_MAKE_PIZZACHILI_PROVIDER

} // namespace SEQAN_NAMESPACE_MAIN

#endif // SEQAN_HEADER_PIZZACHILI_API_H
