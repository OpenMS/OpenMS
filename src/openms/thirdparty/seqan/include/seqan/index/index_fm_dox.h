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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================

/*!
 * @defgroup FMIndexRankDictionarySpec FMIndex RankDictionary Specialisations
 * @brief Tags that can be chosen to specify a certain @link RankDictionary @endlink.
 *
 * @tag FMIndexRankDictionarySpec#WT
 * @brief Tag that specifies the @link FMIndex @endlink to use a wavelet tree as the occurrence table.
 *
 * @tag FMIndexRankDictionarySpec#SBM
 * @brief Tag that specifies the @link FMIndex @endlink to use a StringSet of rank support bis strings as the occurrence table.
 *
 */

/*!
 * @defgroup FMIndexCompressionSpec FMIndex Compression Specialisations
 * @brief Tags that can be chosen to specify if the index stored the text or not.
 *
 * @tag FMIndexCompressionSpec#CompressText
 * @brief Tag to select a FM index variant that can be used such that it is not
 *        necessary to store the text after index construction. This index is
 *        very space efficient.
 *
 * @tag FMIndexCompressionSpec#void
 * @brief The text is kept and not cleared. This FM Index version is faster but more memory is required.
 */

/*!
 * @defgroup FMIndexFibres FM Index Fibres
 * 
 * @brief Tag to select a specific fibre of a @link FMIndex @endlink.
 * 
 * @section Remarks
 * 
 * These tags can be used to get @link Fibre Fibres @endlink of a FM index.
 * 
 * @see Fibre
 * @see getFibre
 * 
 * @tag FMIndexFibres#FibreText
 * 
 * @brief The original text of the index.
 * 
 * @tag FMIndexFibres#FibreSaLfTable
 * 
 * @brief The lf table as well as the compressed suffix array.
 * 
 * @section Remarks
 * 
 * This tag can only be used with the functions @link indexRequire @endlink or
 * @link indexSupplied @endlink.
 * 
 * @tagFMIndexFibres#FibreSA
 * 
 * @brief The compressed suffix array of the text.
 * 
 * @tag FMIndexFibres#FibrePrefixSumTable
 * 
 * @brief The prefix sum table of the index.
 * 
 * @tag FMIndexFibres#FibreLfTable
 * 
 * @brief The lf table.
 */

/*!
 * @class FMIndex
 * 
 * @extends Index
 * 
 * @headerfile seqan/index.h
 * 
 * @brief An index based on the Burrows-Wheeler transform.
 * 
 * @signature Index<TText, FMIndex<TOccSpec, TSpec> >
 * 
 * @tparam TOccSpec Occurrence table specialisation.The tags are really
 *                  shortcuts for the different @link SentinelRankDictionary
 *                  @endlinks Types: @link FMIndexRankDictionarySpec#WT @endlink, @link FMIndexRankDictionarySpec#SBM @endlink, Default: @link FMIndexRankDictionarySpec#WT @endlink
 *
 * @tparam TSpec FM index specialisation. Types: @link FMIndexCompressionSpec#CompressText @endlink, @link FMIndexCompressionSpec#void @endlink, Default: @link FMIndexCompressionSpec#void @endlink
 *
 * @tparam TText The text type. Types: @link String @endlink, @link StringSet @endlink
 */

/*!
 * @fn FMIndex#begin
 * 
 * @brief Returns an iterator pointing to the root not of the virtual prefix
 *        trie of the reversed text of the index.
 * 
 * @signature Iterator begin(index, tag)
 * 
 * @param index The index to be traversed.
 * @param tag The specialisation of the iterator to be returned by the function.
 *            Types: VSTree Iterator
 * 
 * @return TReturn Returns an iterator pointing to the root not of the virtual
 *                 prefix trie of the reversed text of the the index. Types:
 *                 nolink:<tt>The result of Iterator<Index<TText, TIndexSpec>,
 *                 TSpec >::Type</tt>
 */

/*!
 * @fn FMIndex#toSuffixPosition
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This function computes the position of a specified position in the
 *        suffix array (additionally containing entries for the sentinels. The
 *        returned position correspond to the suffix array of the original text
 *        without sentinels.
 * 
 * @signature toSuffixPosition(fmIndex, pos, offset)
 * 
 * @param pos The position in the suffix array of the fm index (with sentinels).
 *            Types: @link UnsignedIntegerConcept @endlink
 * @param fmIndex The FM index.
 * @param offset The number of sequences in the original text. Types:
 *               @link UnsignedIntegerConcept @endlink
 */

/*!
 * @defgroup CompressedSAFibres  CompressedSA Fibres
 * 
 * @brief Tag to select a specific fibre of a @link CompressedSA @endlink.
 * 
 * @see Fibre
 * @see getFibre
 * 
 * @tag CompressedSAFibres#FibreSparseString
 * 
 * @brief The sparse string.
 */

/*!
 * @class CompressedSA
 * 
 * @headerfile seqan/index.h
 * 
 * @brief A suffix array storing only a few suffix array entries and computing
 *        the remaining on demand.
 * 
 * @signature CompressedSA<TSparseString, TLfTable, TSpec>
 * 
 * @tparam TSpec Possibility to specialise a compressed suffix array. Default:
 *               void.
 * @tparam TLfTable The lfTable containg an occurrence table and a prefix sum
 *                  table. Types: @link LfTable @endlink
 * @tparam TSparseString The string containing specific suffix array entries.
 *                       Types: @link SparseString @endlink
 * 
 * @section Remarks
 * 
 * The compressed suffix array can only be used together with the @link LfTable @endlink.
 */

/*!
 * @fn CompressedSA#clear
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Resets an object.
 * 
 * @signature clear(compressedSA)
 * 
 * @param compressesSA The compressed suffix array to be cleared.
 */

/*!
 * @fn CompressedSA#empty
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Checks whether or not a compressed suffix array contains any elements.
 * 
 * @signature bool empty (compressedSA)
 * 
 * @param compressesSA The compressed suffix array to be cleared.
 *
 * @return bool Returns true if the compressed suffix array is empty and false otherwise.
 */

/*!
 * @fn CompressedSA#createCompressedSa
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions creates a compressed suffix array with a specified
 *        compression factor.
 * 
 * @signature void createCompressedSa(compressedSA, completeSA,
 *            compressionFactor [,offset])
 * 
 * @param compressedSA The compressed suffix array.
 *
 * @param compressionFactor The compression factor. A compression factor of x
 *                          means that the compressed suffix array specifically
 *                          stores a value for every x values in the complete
 *                          suffix array. Types: @link UnsignedIntegerConcept @endlink
 * @param completeSA A complete suffix array containing all values. Types: @link String @endlink
 */

/*!
 * @fn CompressedSA#getFibre
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns a specific fibre of a compressed suffix array.
 * 
 * @signature getFibre(compressedSA, fibreTag)
 * 
 * @param fibreTag A tag that identifies the @link Fibre @endlink. Types:
 *                 @link CompressedSAFibres CompressedSA Fibres @endlink
 * @param compressedSA The container holding the fibre.
 * 
 * @return TReturn A reference to the @link Fibre @endlink object.
 */

/*!
 * @fn CompressedSA#length
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the number of elements in the compressed suffix array.
 * 
 * @signature TSize length(compressedSA)
 * 
 * @param compressesSA The compressed suffix array.
 *
 * @return TSize The number of elements in the compressed suffix array. Types: The result of @link Size @endlink of the
 * compressed suffix array.
 */
/*!
 * @fn CompressedSA#resize
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Resets the number of elements in the compressed suffix array.
 * 
 * @signature TSize resize(compressedSA, newLenght)
 * 
 * @param compressesSA The compressed suffix array.
 * @param newLength The number of elements which should be stored in the compressed suffix array. Types: @link
 * UnsignedIntegerConcept @endlink.
 *
 * @return TSize The number of elements in the compressed suffix array. Types: The result of @link Size @endlink of the
 * compressed suffix array.
 *
 * @section Note If the new length is smaller than the actual one then the last <bb>x<bb> items of the compressed suffix array
 * are deleted with x = oldLength - newLength.
 */

/*!
 * @fn CompressedSA#open
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions opens a compressed suffix array from disk.
 * 
 * @signature open(compressedSA, fileName [, mode])
 * 
 * @param mode The combination of flags defining how the file should be
 *             opened.To open a file read-only, write-only or to read and write
 *             use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *             <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *             <tt>OPEN_CREATE</tt>.To append a file if existing add
 *             <tt>OPEN_APPEND</tt>.To circumvent problems, files are always
 *             opened in binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *             OPEN_APPEND</tt>
 * @param compressedSA The compressed suffix array to be opened.
 * @param fileName C-style character string containing the file name.
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on success.
 */


/*!
 * @fn CompressedSA#save
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions saves a compressed suffix array to disk.
 * 
 * @signature save(compressedSA, fileName [, mode])
 * 
 * @param mode The combination of flags defining how the file should be
 *             opened.To open a file read-only, write-only or to read and write
 *             use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *             <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *             <tt>OPEN_CREATE</tt>.To append a file if existing add
 *             <tt>OPEN_APPEND</tt>.To circumvent problems, files are always
 *             opened in binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *             OPEN_APPEND</tt>
 * @param compressedSA The compressed suffix array to be opened.
 * @param fileName C-style character string containing the file name.
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on success.
 */

/*!
 * @fn CompressedSA#setLfTable
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Set the LfTable of the compressed suffix array.
 * 
 * @signature setLfTable(compressedSa, lfTable)
 * 
 * @param lfTable The LfTable Types: LfTable
 * @param compressedSA The compressed suffix array.
 */

/*!
 * @fn CompressedSA#value
 * 
 * @brief Returns the value stored at a specified position in the compressed
 *        suffix-array.
 * 
 * @signature Value value(compressedSA, pos)
 * 
 * @param pos Position at which to access the suffix array. Types:
 *            @link UnsignedIntegerConcept @endlink
 * @param compressedSA The compressed suffix array to access.
 * 
 * @section Remarks
 * 
 * Note that the compressed suffix array is read only. Therefore a const
 * reference is return by this function.
 */

/*!
 * @defgroup LFTableFibres LF Table Fibres
 * 
 * @brief Tag to select a specific fibre of a @link LfTable @endlink.
 * 
 * @section Remarks
 * 
 * These tags can be used to get @link Fibre Fibres @endlink of a @link LfTable @endlink.
 * 
 * @see Fibre
 * @see getFibre
 * 
 * @tag LFTableFibres#FMTablePrefixSumTable
 * 
 * @brief The prefix sum table of the lf table.
 * 
 * @tag LFTableFibres#FibreOccTable
 * 
 * @brief The occurrence table of the lf table.
 */

/*!
 *!
 * @class LfTable
 * 
 * @headerfile seqan/Index.h
 * 
 * @brief LfTable is an object storing all necessary information for the LF-
 *        mapping.
 * 
 * @signature LfTable<TOccTable, TPrefixSumTable>
 * 
 * @tparam TPrefixSumTable The specialisation tag. Default: @link String @endlink
 * @tparam TOccTable The occurrence table data structure. Types: @link RankDictionary @endlink, @link SentinelRankDictionary @endlink.
 */

/*!
 * @fn LfTable#clear
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Resets the LF table.
 * 
 * @signature clear(lfTable)
 * 
 * @param lfTable The LF table to be cleared.
 */

/*!
  * @fn LfTable#empty
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Clears the LF table.
 * 
 * @signature empty(lfTable)
 * 
 * @param lfTable The LF table to be checked.
 * 
 * @return TReturn <tt>true</tt> if the LF table is empty, <tt>false</tt>
 *                 otherwise. Types: nolink:<tt>bool</tt>
 */

/*!
 * @fn createLfTable
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Creates the LF table
 * 
 * @signature createLfTable(lfTable, text)
 * 
 * @param lfTable The LF table to be constructed.
 * @param text The underlying text Types: @link String @endlink.
 * 
 * @return TReturn <tt>true</tt> on successes, <tt>false</tt> otherwise.
 */

/*!
 * @fn LfTable#getFibre
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns a specific fibre of a LF table.
 * 
 * @signature getFibre(lfTable, fibreTag)
 * 
 * @param fibreTag A tag that identifies the @link Fibre @endlink. Types: @link LFTableFibres @endlink
 * @param container The container holding the fibre.
 * 
 * @return TReturn A reference to the @link Fibre @endlink object.
 */

/*!
 * @fn LfTable#lfMapping
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the position of an character at a specified position of L in
 *        F. L corresponds to the last column of the sorted cyclic rotations of
 *        the original text, while F correspond to the first column.
 * 
 * @signature lfMapping(lfTable, pos)
 * 
 * @param lfTable The @link LfTable @endlink holding the occurrence and prefix
 *                sum table.
 *
 * @param pos The position in L. Types: @link UnsignedIntegerConcept @endlink
 * 
 * @return TReturn Returns the position of the character L[c] in F. The returned
 *                 position is of the same type as pos. Types: The type of the position.
 */

/*!
 * @fn LfTable#open
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions loads a LF table from disk.
 * 
 * @signature open(lfTable, fileName [, openMode])
 * 
 * @param openMode The combination of flags defining how the file should be
 *                 opened.To open a file read-only, write-only or to read and
 *                 write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                 <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *                 <tt>OPEN_CREATE</tt>.To append a file if existing add
 *                 <tt>OPEN_APPEND</tt>.To circumvent problems, files are always
 *                 opened in binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                 OPEN_APPEND</tt>
 * @param lfTable The lfTable. Types: LfTable
 * @param fileName C-style character string containing the file name.
 * 
 * @return TReturn A nolink:<tt>bool</tt> which is <tt>true</tt> on success.
 */

/*!
 * @fn LfTable#save
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions saves a LF table to disk.
 * 
 * @signature save(lfTable, fileName [, openMode])
 * 
 * @param openMode The combination of flags defining how the file should be
 *                 opened.To open a file read-only, write-only or to read and
 *                 write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                 <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *                 <tt>OPEN_CREATE</tt>.To append a file if existing add
 *                 <tt>OPEN_APPEND</tt>.To circumvent problems, files are always
 *                 opened in binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                 OPEN_APPEND</tt>
 * @param lfTable The dictionary. Types: LfTable
 * @param fileName C-style character string containing the file name.
 * 
 * @return TReturn A nolink:<tt>bool</tt> which is <tt>true</tt> on success.
 */

/*!
 * @defgroup PrefixSumTableFibres PrefixSumTable Fibres
 * 
 * @brief Tag to select a specific fibre of a @link PrefixSumTableFibres @endlink.
 * 
 * @see Fibre
 * @see getFibre
 * 
 * @tag PrefixSumTableFibres#FibreEntries
 * 
 * @brief The entries in the prefix sum table.
 */

/*!
 * @class PrefixSumTable
 * 
 * @headerfile seqan/Index.h
 * 
 * @brief The prefix-sum table is a data structure which stores for each
 *        character the number of smaller lexicographic smaller characters in a
 *        given text.
 * 
 * @signature PrefixSumTable<TChar, TSpec>
 * 
 * @tparam TSpec A specialisation tag. Default: void
 * @tparam TChar The character type
 */

/*!
 * @fn PrefixSumTable#clear
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Resets the prefix sum table.
 * 
 * @signature clear(prefixSumTable)
 * @param prefixSumTable The prefix sum table to be cleared.
 */

/*!
 * @fn PrefixSumTable#createPrefixSumTable
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Creates the prefix sum table
 * 
 * @signature createPrefixSumTable(prefixSumTable, text)
 * 
 * @param text The underlying text. Types: @link String @endlink, @link StringSet @endlink
 * @param prefixSumTable The prefix sum table to be constructed.
 */

/*!
 * @fn PrefixSumTable#getAlphabetSize
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the number of different characters in the prefix sum table.
 * 
 * @signature getAlphabetSize(prefixSumTable)
 * 
 * @param prefixSumTable A prefix sum table.
 */

/*!
 * @fn PrefixSumTable#getCharacterPosition
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the position of a given character within the prefix sum table.
 * 
 * @signature getCharacterPosition(prefixSumTable, character)
 * 
 * @param character A character.
 * @param prefixSumTable A prefix sum table.
 */

/*!
 * @fn PrefixSumTable#getCharacter
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the character of a given position within the prefix sum table.
 * 
 * @signature getCharacter(prefixSumTable, pos)
 * 
 * @param pos A position. Types: @link UnsignedIntegerConcept @endlink
 * @param prefixSumTable A prefix sum table.
 */

/*!
 * @fn PrefixSumTable#getPrefixSum
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the prefix sum of a given position.
 * 
 * @signature getPrefixSum(prefixSumTable, pos)
 * 
 * @param pos A position. Types: @link UnsignedIntegerConcept @endlink
 * @param prefixSumTable A prefix sum table.
 */

/*!
 * @fn PrefixSumTable#getValue
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the prefix sum of a given position.
 * 
 * @signature getValue(prefixSumTable, pos)
 * 
 * @param pos A position. Types: @link UnsignedIntegerConcept @endlink
 * @param prefixSumTable A prefix sum table.
 */

/*!
 * @fn PrefixSumTable#getFibre
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns a specific fibre of a prefix-sum table.
 * 
 * @signature getFibre(prefixSumTable, fibreTag)
 * 
 * @param fibreTag A tag that identifies the @link Fibre @endlink. Types:
 *                 @link PrefixSumTableFibres PrefixSumTanble Fibres @endlink.
 * @param prefixSumTable The prefix sum table.
 * 
 * @return TReturn A reference to the @link Fibre @endlink object.
 */

/*!
 * @fn PrefixSumTable#length
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the number of different characters in the prefix-sum table.
 * 
 * @signature length(lfTable)
 * 
 * @param lfTable The prefix-sum table.
 * 
 * @return TReturn Returns the number of different characters in the prefix-sum
 *                 table.If the type of the characters of the prefix-sum table
 *                 consists of more than 8 bit only the characters actually
 *                 occurring in the original text are accounted for when calling
 *                 length. Types: @link Size @endlink of the prefix-sum table.
 */

/*!
 * @fn PrefixSumTable#prefixSum
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns a reference to the entry of the prefix sum table of a given
 *        position.
 * 
 * @signature prefixSum(prefixSumTable, pos)
 * 
 * @param pos A position. Types: @link UnsignedIntegerConcept @endlink.
 * @param prefixSumTable A prefix sum table.
 */

/*!
 * @fn PrefixSumTable#resize
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Resize the prefix-sum table to be able to store more or less
 *        characters.
 * 
 * @signature resize(prefixSumTable, size [,value, resizeTag])
 * 
 * @param resizeTag Specifies the strategy that is applied if the capacity of
 *                  <tt>object</tt> is less than <tt>newLength</tt>. (optional)
 *                  Types: @link OverflowStrategy @endlink Default: Specified by @link
 *                  DefaultOverflowExplicit @endlink.
 * @param prefixSumTable A prefix sum table. Types: PrefixSumTable
 * @param value The value to be used to initialize the new storage.
 * @param size The new size. Types: @link UnsignedIntegerConcept @endlink
 */

/*!
 * @fn PrefixSumTable#setPrefixSum
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns a reference to the entry of the prefix-sum table of a given
 *        position.
 * 
 * @signature setPrefixSum(prefixSumTable, value, pos)
 * 
 * @param pos A position. Types: @link UnsignedIntegerConcept @endlink
 * @param value A specified value to be inserted.
 * @param prefixSumTable A prefix sum table.
 */

/*!
  * @fn PrefixSumTable#open
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions loads a prefix-sum table from disk.
 * 
 * @signature open(prefixSumTable, fileName [, openMode])
 * 
 * @param openMode The combination of flags defining how the file should be
 *                 opened.To open a file read-only, write-only or to read and
 *                 write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                 <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *                 <tt>OPEN_CREATE</tt>.To append a file if existing add
 *                 <tt>OPEN_APPEND</tt>.To circumvent problems, files are always
 *                 opened in binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                 OPEN_APPEND</tt>
 * @param prefixSumTable The prefix-sum table.
 * @param fileName C-style character string containing the file name.
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on success.
 */

/*!
 * @fn PrefixSumTable#save
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions saves a prefix-sum table to disk.
 * 
 * @signature save(prefixSumTable, fileName [, openMode])
 * 
 * @param openMode The combination of flags defining how the file should be
 *                 opened.To open a file read-only, write-only or to read and
 *                 write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                 <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *                 <tt>OPEN_CREATE</tt>.To append a file if existing add
 *                 <tt>OPEN_APPEND</tt>.To circumvent problems, files are always
 *                 opened in binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                 OPEN_APPEND</tt>
 * @param prefixSumTable The prefix-sum table.
 * @param fileName C-style character string containing the file name.
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on success.
 */

/*!
 * @class RankDictionary
 * 
 * @headerfile seqan/index.h
 * 
 * @brief A rank dictionary is a data structure to store the rank of an element
 *        of a sequence at every position of the sequence.
 * 
 * @signature RankDictionary<TSpec>
 * 
 * @tparam TSpec The rank dictionary specialisation. Types: WaveletTree,
 *               SequenceBitMask Default: @link WaveletTree @endlink
 */

/*!
 * @class WaveletTree
 *
 * @extends RankDictionary
 * 
 * @headerfile seqan/index.h
 * 
 * @brief A wavelet tree is a tree like binary encoding of a text.
 * 
 * @signature template <typename TValue>
 *            RankDictionary<WaveletTree<TValue> >
 * 
 * @tparam TValue The value type of the wavelet tree.
 * 
 * @section Remarks
 * 
 * The nodes of a wavelet tree consist of a bit string as well as a character c.
 * In each level of the tree, characters smaller than c are represented as a 0
 * while character greater or equal to c are represented with a 1. The
 * characters represented by a 0 form the string to be represented by the left
 * subtree while characters represented by a 1 form the string of the right
 * subtree. Therefore, only the bit string of the root node represents all
 * characters while all other nodes represent subsets.
 */

/*!
 * @defgroup WaveletTreeFibres WaveletTree Fibres
 * 
 * @brief Tag to select a specific fibre (e.g. table, object, ...) of a @link
 *        WaveletTree @endlink.
 * 
 * @section Remarks
 * 
 * These tags can be used to get @link Fibre Fibres @endlink of a @link
 * WaveletTree @endlink.
 * 
 * @see Fibre
 * @see getFibre
 * 
 * @tag WaveletTreeFibres#FibreTreeStructure
 * 
 * @brief The wavelet tree structure of the wavelet tree.
 * 
 * @tag WaveletTreeFibres#FibreBitStrings
 * 
 * @brief The string set containing a bit string for each node.
 */

/*!
 * @class SequenceBitMask
 *
 * @extends RankDictionary
 * 
 * @headerfile seqan/index.h
 * 
 * @brief The string set bit string dictionary is a string set of rank support
 *        bit strings for constant time acces of the rank of a specified
 *        character at a specified position.
 * 
 * @signature template <typename TValue>
 *            RankDictionary<SequenceBitMask<TValue> >
 * 
 * @tparam TValue The value type of the .
 * 
 * @section Remarks
 * 
 * This data structure is optimized for very small alphabets, such as @link Dna @endlink
 * or  @link Dna5 @endlink.  Consider using a @link WaveletTree @endlink if
 * your alphabet size is larger.
 */

/*!
 * @class SequenceBitMaskFibres SequenceBitMask Fibres
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Tag to select a specific fibre of a SequenceBitMask.
 *
 * @tag SequenceBitMaskFibres#FibreBitStrings The string set containing a bit string for each character.
 * 
 * @see Fibre
 * @see getFibre
 */

/*!
 * @mfn RankDictionary#Fibre
 *
 * @brief Returns the type of a specified fibre of a @link RankDictionary @endlink.
 *
 * signature Fibre<RankDictionary, FibreSpec>::Type
 *
 * @tparam FibreSpec The Fibre of interest. Types: @link WaveletTreeFibres @endlink, @link SequenceBitMaskFibres @endlink.
 *
 */

/*!
 * @fn RankDictionary#clear
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Resets the rank dictionary.
 * 
 * @signature clear(dictionary)
 * 
 * @param dictionary The rank dictionary to be cleared.
 */

/*!
 * @fn RankDictionary#empty
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns whether or not the rank dictionary is empty.
 * 
 * @signature empty(dictionary)
 * 
 * @param dictionary The rank dictionary to be checked.
 * 
 * @return TReturn <tt>true</tt> if the dictionary is empty, <tt>false</tt>
 *                 otherwise.
 */

/*!
 * @fn RankDictionary#getValue
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the character of a specified position.
 * 
 * @signature getValue(dictionary, pos)
 * 
 * @param pos The position. Types: @link UnsignedIntegerConcept @endlink.
 * @param dictionary The dictionary.
 */

/*!
  * @fn RankDictionary#getFibre
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns a specific fibre of a dictionary.
 * 
 * @signature getFibre(dictionary, fibreTag)
 * 
 * @param fibreTag A tag that identifies the @link Fibre @endlink. Types:
 *                 @link SequenceBitMaskFibres SequenceBitMask Fibres @endlink, @link WaveletTreeFibres WaveletTree Fibres @endlink.
 * @param dictionary The dictionary holding the fibre.
 * 
 * @return TReturn A reference to the @link Fibre @endlink object.
 */

/*!
 * @fn RankDictionary#countOccurrences
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the rank (number of occurrences) of a specified character up
 *        to a specified position.
 * 
 * @signature countOccurrences(dictionary, character, pos)
 * 
 * @param character The character of interest.
 * @param pos The position (which is also included in the rank computation).
 * @param dictionary The dictionary. 
 *
 * @return TReturn The rank (number of occurrences) of a specified character up
 *                 to a specified position. Types: <tt>unsigned</tt>
 */

/*!
 * @fn RankDictionary#createRankDictionary
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions creates the dictionary.
 * 
 * @signature createRankDictionary(dictionary, text)
 * 
 * @param text A text to be transfered into a wavelet tree. Types: @link String #endlink
 * @param dictionary The dictionary.
 */

/*!
 * @fn RankDictionary#open
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions loads a rank dictionary from disk.
 * 
 * @signature open(dictionary, fileName [, openMode])
 * 
 * @param openMode The combination of flags defining how the file should be
 *                 opened.To open a file read-only, write-only or to read and
 *                 write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                 <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *                 <tt>OPEN_CREATE</tt>.To append a file if existing add
 *                 <tt>OPEN_APPEND</tt>.To circumvent problems, files are always
 *                 opened in binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                 OPEN_APPEND</tt>
 * @param dictionary The dictionary.
 * @param fileName C-style character string containing the file name.
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on success.
 */

/*!
 * @fn RankDictionary#save
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions saves a dictionary to disk.
 * 
 * @signature save(dictionary, fileName [, openMode])
 * 
 * @param openMode The combination of flags defining how the file should be
 *                 opened.To open a file read-only, write-only or to read and
 *                 write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                 <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *                 <tt>OPEN_CREATE</tt>.To append a file if existing add
 *                 <tt>OPEN_APPEND</tt>.To circumvent problems, files are always
 *                 opened in binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                 OPEN_APPEND</tt>
 * @param dictionary The dictionary.
 * @param fileName C-style character string containing the file name.
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on success.
 */

/*!
 * @defgroup RankSupportBitStringFibres RankSupportBitString Fibres
 * 
 * @brief Tag to select a specific fibre (e.g. table, object, ...) of a @link
 *        RankSupportBitString @endlink.
 * 
 * @section Remarks
 * 
 * These tags can be used to get @link Fibre.Fibres @endlink of a rank support
 * bit string.
 * 
 * @see Fibre
 * @see getFibre
 * 
 * @tag RankSupportBitStringFibres#FibreSuperBlocks
 * 
 * @brief The super block string.
 * 
 * @tag RankSupportBitStringFibres#FibreBits
 * 
 * @brief The bit string.
 * 
 * @tag RankSupportBitStringFibres#FibreBlocks
 * 
 * @brief The block string.
 */

/*!
 * @class RankSupportBitString
 * 
 * @headerfile seqan/index.h
 * 
 * @brief A bit string supporting rank queries in constant time.
 * 
 * @signature RankSupportBitString<TSpec>
 * 
 * @tparam TSpec Specialisation tag. Default: void
 * 
 * @section Remarks
 * 
 * The constant rank query time is achieved by evaluating precomputed
 * subsolutions. In order to do so, the bit string is divided into blocks of
 * length l. A super block string stores for each block of l blocks the number
 * of bits set from the beginning. In addition a block string stores the number
 * of bits set in each block from the start of the last super block block.
 * Therefore it is possible to compute the result of a rank query in constant
 * time by adding information from the bit, block and super block string.
 */

/*!
 * @fn RankSupportBitString#appendValue
 * 
 * @headerfile seqan/sequence.h
 * 
 * @brief Appends a bit to a @link RankSupportBitString @endlink.
 * 
 * @signature appendValue(bitString, bit)
 * 
 * @param target A container. Types: RankSupportBitString
 * @param bit Value that is appended to <tt>target</tt>.If the value is
 *              different from 0 it is interpreted as 1. Types: @link UnsignedIntegerConcept @endlink, bool
 */

/*!
 * @fn RankSupportBitString#clear
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Resets an rank-support-bit string.
 * 
 * @signature clear(bitString)
 * 
 * @param bitString The bit string to be cleared.
 */

/*!
 * @fn RankSupportBitString#getRank
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the rank (the number of bits set from the start of the bit
 *        string) of a specified position.
 * 
 * @signature getRank(bitString, pos)
 * 
 * @param bitString The bit string. Types: @link RankSupportBitString @endlink
 * @param pos Position of a bit.
 * 
 * @return TReturn @link Value @endlink of @link Fibre @endlink of the rank-support-bit string.
 */

/*!
 * @fn RankSupportBitString#empty
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Checks whether or not a rank-support-bit string contains any elements.
 * 
 * @signature bool empty(bitString)
 * 
 * @param bitString The rank-support-bit string to be checked.
 *
 * @return bool Returns true if the rank-support-bit string is empty and false otherwise.
 */

/*!
 * @fn RankSupportBitString#isBitSet
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns whether the bit with the given index is set to 1.
 * 
 * @signature isSetBit(bitString, pos)
 * 
 * @param bitString The bit string.
 * @param pos Position of the bit. Types: @link UnsignedIntegerConcept @endlink
 * 
 * @return TReturn Returns whether a specified bit is set or not.
 */

/*!
 * @fn RankSupportBitString#getFibre
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns a specific fibre of a rank-support-bit string.
 * 
 * @signature getFibre(bitString, fibreTag)
 * 
 * @param fibreTag A tag that identifies the @link Fibre @endlink. Types:
 *                 @link RankSupportBitStringFibres RankSupportBitString Fibres @endlink.
 * @param bitString The rank-support-bit string holding the fibre.
 * 
 * @return TReturn A reference to the @link Fibre @endlink object.
 */

/*!
 * @fn RankSupportBitString#length
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the number of bits in the rank-support-bit string.
 * 
 * @signature length(bitString)
 * 
 * @param bitString The rank-support-bit string.
 * 
 * @return TReturn Types: @link Value @endlink of @link Fibre @endlink of the rank-support-bit string.
 */

/*!
 * @fn RankSupportBitString#resize
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Resets the number of bits in the rank-support-bit string.
 * 
 * @signature TSize resize(bitString, newLenght)
 * 
 * @param bitString The rank-support-bit string.
 * @param newLength The number of elements which should be stored in the compressed suffix array. Types: @link
 * UnsignedIntegerConcept @endlink.
 *
 * @return TSize The number of elements in the rank-support-bit string. Types: The result of @link Size @endlink of the
 * rank-support-bit string.
 *
 * @section Note If the new length is smaller than the actual one then the last <bb>x<bb> items of the rank-support-bit string
 * are deleted with x = oldLength - newLength.
 */

/*!
 * @fn RankSupportBitString#setBitTo
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Set the bit with the given position to the given value.
 * 
 * @signature setBitTo(bitString, pos, value)
 * 
 * @param pos Position of the bit. Types: @link UnsignedIntegerConcept @endlink
 * @param bitString The bit string. 
 * @param bit The value of the bit. Note that values different from 0 are
 *            interpreted as 1.
 * 
 * @return TReturn <tt>void</tt>
 * 
 * @section Examples
 * 
 * @see RankSupportBitString#isBitSet
 */

/*!
 * @fn RankSupportBitString#open
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions opens a @link RankSupportBitString @endlink from disk.
 * 
 * @signature open(bitString, fileName [, openMode])
 * 
 * @param openMode The combination of flags defining how the file should be
 *                 opened.To open a file read-only, write-only or to read and
 *                 write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                 <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *                 <tt>OPEN_CREATE</tt>.To append a file if existing add
 *                 <tt>OPEN_APPEND</tt>.To circumvent problems, files are always
 *                 opened in binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                 OPEN_APPEND</tt>
 * @param bitString The bit string to be opened.
 * @param fileName C-style character string containing the file name.
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on success.
 */

/*!
 * @fn RankSupportBitString#save
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions saves a @link RankSupportBitString @endlink to disk.
 * 
 * @signature save(bitString, fileName [, openMode])
 * 
 * @param openMode The combination of flags defining how the file should be
 *                 opened.To open a file read-only, write-only or to read and
 *                 write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                 <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *                 <tt>OPEN_CREATE</tt>.To append a file if existing add
 *                 <tt>OPEN_APPEND</tt>.To circumvent problems, files are always
 *                 opened in binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                 OPEN_APPEND</tt>
 * @param bitString The bit string to be saved.
 * @param fileName C-style character string containing the file name.
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on success.
 */

/*!
 * @class RightArrayBinaryTree
 * 
 * @headerfile seqan/index.h
 * 
 * @brief A special format to encode the structure of a wavelet tree. The
 *        structure is very space efficient because only one position is stored
 *        which encodes where the left and right subtree of a given node exist.
 * 
 * @signature RightArrayBinaryTree<TValue, TSpec>
 * 
 * @tparam TValue The value type - type of the stored characters.
 *
 * @tparam TSpec The wavelet tree structure specialisation. Default: void.
 */

/*!
 * @defgroup RightArrayBinaryTreeFibres RightArrayBinaryTree Fibres
 * 
 * @brief Tag to select a specific fibre (e.g. table, object, ...) of a @link
 *        RightArrayBinaryTree @endlink.
 * 
 * @section Remarks
 * 
 * These tags can be used to get @link Fibre Fibres @endlink of a @link RightArrayBinaryTree @endlink.
 * 
 * @see Fibre
 * @see getFibre
 * 
 * @tag RightArrayBinaryTreeFibres#FibreTreeStructureEncoding
 * 
 * @brief The string encoding the wavelet tree structure.
 */

/*!
 * @fn RightArrayBinaryTree#clear
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Resets a right-array-binary tree.
 * 
 * @signature clear(rightArrayBinaryTree)
 * 
 * @param rightArrayBinaryTree The right-array-binary tree to be cleared.
 */

/*!
 * @fn RightArrayBinaryTree#createRightArrayBinaryTree
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Computes the right-array-binary tree of a text.
 * 
 * @signature createRightArrayBinaryTree(rightArrayBinaryTree, text)
 * 
 * @param text A text. Types: @link String @endlink
 * @param rightArrayBinaryTree A wavelet tree structure.
 */

/*!
 * @fn RightArrayBinaryTree#empty
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Checks whether or not a right-array-binary tree contains any elements.
 * 
 * @signature bool empty(rightArrayBinaryTree)
 * 
 * @param rightArrayBinaryTree The right-array-binary tree to be checked.
 *
 * @return bool Returns true if the rank-support-bit string is empty and false otherwise.
 */

/*!
 * @fn RightArrayBinaryTree#getFibre
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns a specific fibre of a right-array-binary tree.
 * 
 * @signature getFibre(rightArrayBinaryTree, fibreTag)
 * 
 * @param fibreTag A tag that identifies the @link Fibre @endlink. Types:
 *                 @link RightArrayBinaryTreeFibres RightArrayBinaryTree Fibres @endlink.
 * @param container The container holding the fibre.
 * 
 * @return TReturn A reference to the @link Fibre @endlink object.
 */

/*!
 * @fn RightArrayBinaryTree#open
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions loads a @link RightArrayBinaryTree @endlink from disk.
 * 
 * @signature open(rightArrayBinaryTree, fileName [, openMode])
 * 
 * @param openMode The combination of flags defining how the file should be
 *                 opened.To open a file read-only, write-only or to read and
 *                 write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                 <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *                 <tt>OPEN_CREATE</tt>.To append a file if existing add
 *                 <tt>OPEN_APPEND</tt>.To circumvent problems, files are always
 *                 opened in binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                 OPEN_APPEND</tt>
 * @param rightArrayBinaryTree The rightArrayBinaryTree. 
 *
 * @param fileName C-style character string containing the file name.
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on success.
 */

/*!
 * @fn RightArrayBinaryTree#save
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions saves a @link RightArrayBinaryTree @endlink to disk.
 * 
 * @signature save(rightArrayBinaryTree, fileName [, openMode])
 * 
 * @param openMode The combination of flags defining how the file should be
 *                 opened.To open a file read-only, write-only or to read and
 *                 write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                 <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *                 <tt>OPEN_CREATE</tt>.To append a file if existing add
 *                 <tt>OPEN_APPEND</tt>.To circumvent problems, files are always
 *                 opened in binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                 OPEN_APPEND</tt>
 * @param rightArrayBinaryTree The rightArrayBinaryTree. 
 *
 * @param fileName C-style character string containing the file name.
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on success.
 */

/*!
 * @class RightArrayBinaryTreeIterator RightArrayBinaryTree Iterator
 * 
 * @extends Iter
 * 
 * @headerfile seqan/index.h
 * 
 * @brief An iterator for @link RightArrayBinaryTree @endlink.
 * 
 * @signature Iter<RightArrayBinaryTree, TSpec >
 * 
 * @tparam TSpec Specialisation Tag. Types: TopDownIterator
 * @tparam RightArrayBinaryTree The @link RightArrayBinaryTree @endlink. Types:
 *                              @link RightArrayBinaryTree @endlink
 */

/*!
 * @fn RightArrayBinaryTree#begin
 * 
 * @headerfile seqan/index.h
 * 
 * @brief The begin (root) of a @link RightArrayBinaryTree @endlink.
 * 
 * @signature Iterator begin(rightArrayBinaryTree, iterSpec)
 * 
 * @param rightArrayBinaryTree The right-array-binary tree.
 *
 * @param iterSpec A specialisation tag. Types: TopDown<>, TopDown<ParentLinks<> >
 * 
 * @return TReturn An iterator to the first item in <tt>object</tt>.
 *                 Metafunctions: Metafunction.Iterator
 */

/*!
 * @fn RightArrayBinaryTreeIterator#container
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Container of an iterator.
 * 
 * @signature Container container(iterator)
 * 
 * @param iterator An iterator.
 * 
 * @return TReturn The container that <tt>iterator</tt> traverses.
 */

/*!
 * @fn RightArrayBinaryTree#end
 * 
 * @headerfile seqan/index.h
 * 
 * @brief The end (rigthmost laef) of a @link RightArrayBinaryTree @endlink.
 * 
 * @signature Iterator end(rightArrayBinaryTree, iterSpec)
 * 
 * @param rightArrayBinaryTree The right-array-binary tree.
 *
 * @param iterSpec A specialisation tag. Types: TopDown<>, TopDown<ParentLinks<> >
 * 
 * @return TReturn An iterator to the first item in <tt>object</tt>.
 *                 Metafunctions: Metafunction.Iterator
 */

/*!
 * @fn RightArrayBinaryTreeIterator#getCharacter
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This function returns the pivot character of the node the iterator
 *        currently points to.
 * 
 * @signature getCharacter(iterator)
 * 
 * @param iterator The iterator.
 */

/*!
 * @fn RightArrayBinaryTreeIterator#getLeftChildPos
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the position in @link RightArrayBinaryTree @endlink of the
 *        left child node.
 * 
 * @signature getLeftChildPos(iterator)
 * 
 * @param iterator The iterator.
 */

/*!
 * @fn RightArrayBinaryTreeIterator#getSubTreeSize
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the number of vertices in the subtree starting at the position
 *        an iterator points to.
 * 
 * @signature getSubTreeSize(iterator)
 * 
 * @param iterator The iterator.
 */

/*!
 * @fn RightArrayBinaryTreeIterator#getPosition
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the position of the iterator in the host.
 * 
 * @signature getPosition(iterator)
 * 
 * @param iterator The iterator.
 */

/*!
 * @fn RightArrayBinaryTreeIterator#getRightChildPos
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the position in @link RightArrayBinaryTree @endlink of the
 *        right child node.
 * 
 * @signature getLeftChildPos(iterator)
 * 
 * @param iterator The iterator.
 */

/*!
 * @fn RightArrayBinaryTreeIterator#goDown
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Iterates down the leftmost edge in a @link RightArrayBinaryTree @endlink.
 * 
 * @signature bool goDown(iterator)
 *
 * @param iterator The iterator
 * 
 * @return TReturn <tt>true</tt> if an edge to go down exists,
 *                 otherwise <tt>false</tt>. Types: bool
 */

/*!
 * @fn RightArrayBinaryTreeIterator#goLeftChild
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Sets the iterator to the left child of the current node if it exists
 *        and returns true, otherwise the iterator does not change position and
 *        the function returns false.
 * 
 * @signature bool goLeftChild(iterator)
 * 
 * @param iterator The iterator
 * 
 * @return TReturn <tt>true</tt> if the edge to go down exists,
 *                 otherwise <tt>false</tt>.
 */

/*!
 * @fn RightArrayBinaryTreeIterator#goRight
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Iterates to the next sibling in a @link RightArrayBinaryTree @endlink.
 * 
 * @signature goRight(iterator)
 * 
 * @param iterator The iterator
 * 
 * @return TReturn <tt>true</tt> if the iterator could be moved, otherwise
 *                 <tt>false</tt>. Types: nolink:bool
 */

/*!
 * @fn RightArrayBinaryTreeIterator#goRightChild
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Sets the iterator to the right child of the current node if it exists
 *        and returns true, otherwise the iterator does not change position and
 *        the function returns false.
 * 
 * @signature bool goRightChild(iterator)
 * 
 * @param iterator The iterator.
 * 
 * @return TReturn <tt>true</tt> if the edge to go down exists,
 *                 otherwise <tt>false</tt>.
 */

/*!
 * @fn RightArrayBinaryTreeIterator#goUp
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Iterates up one edge to the parent in a @link RightArrayBinaryTree @endlink.
 * 
 * @signature goUp(iterator)
 * 
 * @param iterator The iterator.
 * 
 * @return TReturn <tt>true</tt> if the iterator could be moved, otherwise
 *                 <tt>false</tt>. Types: nolink:bool
 */

/*!
 * @fn RightArrayBinaryTreeIterator#isLeaf
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Tests whether a given node is a leaf or not.
 * 
 * @signature isLeaf(iterator)
 * 
 * @param iterator The iterator.
 * 
 * @return TReturn True if the node is a leaf.
 */

/*!
 * @fn RightArrayBinaryTreeIterator#setCharacter
 * 
 * @headerfile seqan/index.h
 * 
 * @brief The function sets the character of the
 *        node the iterator points to to character.
 * 
 * @signature void setCharacter(iterator, character)
 * 
 * @param character The character to be assigned to a node.
 * @param iterator The iterator.
 */

/*!
 * @fn RightArrayBinaryTreeIterator#isRoot
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Test whether a iterator points to the root node.
 * 
 * @signature bool isRoot(iterator)
 * 
 * @param iterator The iterator.
 * 
 * @return TReturn <tt>true</tt> if <tt>iterator</tt> points to the root of the
 *                 tree, otherwise <tt>false</tt>. Types: nolink:bool
 */

/*!
 * @defgroup SentinelRankDictionaryFibres SentinelRankDictionary Fibres
 * 
 * @brief Tag to select a specific fibre of a @link
 *        SentinelRankDictionary @endlink.
 *  
 * @see Fibre
 * @see getFibre
 * 
 * @tag SentinelRankDictionaryFibres#FibreSentinelPosition
 * 
 * @brief The bit string encoding the position of the sentinel sign.
 * 
 * @tag SentinelRankDictionaryFibres#FibreRankDictionary
 * 
 * @brief The rank dictionary.
 */

/*!
 * @class SentinelRankDictionary
 * 
 * @headerfile seqan/index.h
 * 
 * @brief A rank dictionary, additional storing sentinel character which are not
 *        accounted for in a rank query.
 * 
 * @signature SentinelRankDictionary<TRankDictionary, TSpec>
 * 
 * @tparam TSpec Specialisation
 * @tparam TRankDictionary The rank dictionary of a text.
 */

/*!
 * @fn SentinelRankDictionary#clear
 * @headerfile seqan/index.h
 * 
 * @brief Clears the dictionary.
 * 
 * @signature clear(dictionary)
 * 
 * @param dictionary The rank dictionary to be cleared.
 */

/*!
 * @fn SentinelRankDictionary#sentinelPosition
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns whether a specified position is a sentinel position.
 * 
 * @signature sentinelPosition(dictionary, pos)
 * 
 * @param pos The position. Types: @link UnsignedIntegerConcept @endlink
 * @param dictionary The dictionary.
 */

/*!
 * @fn SentinelRankDictionary#empty
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns whether or not the dictionary is empty.
 * 
 * @signature empty(dictionary)
 * 
 * @param dictionary The rank dictionary to be checked. 
 */

/*!
 * @fn SentinelRankDictionary#getValue
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the character of a specified position.
 * 
 * @signature getCharacter(dictionary, pos)
 * 
 * @param pos The position
 * @param dictionary The rank dictionary.
 */

/*!
 * @fn SentinelRankDictionary#getFibre
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns a specific fibre of a dictionary.
 * 
 * @signature getFibre(dictionary, fibreTag)
 * 
 * @param fibreTag A tag that identifies the @link Fibre @endlink. Types:
 *                 @link SentinelRankDictionaryFibres SentinelRankDictionary Fibres
 * @param dictionary The dictionary holding the fibre.
 * 
 * @return TReturn A reference to the @link Fibre @endlink object.
 */

/*!
 * @fn SentinelRankDictionary#countOccurrences
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the number of occurrences of a specified character from the
 *        start to a specified position.
 * 
 * @signature countOccurrences(dictionary, character, pos)
 * 
 * @param character The character.
 * @param pos The position (which is included in the counting).
 * @param dictionary The dictionary.
 */

/*!
 * @fn SentinelRankDictionary#getSentinelSubstitute
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the character used to substitute the sentinel sign.
 * 
 * @signature getSentinelSubstitute(dictionary)
 * 
 * @param dictionary The dictionary.
 */

/*!
 * @fn SentinelRankDictionary#setSentinelSubstitute
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Sets the character used to substitute the sentinel sign.
 * 
 * @signature setSentinelSubstitute(dictionary, character)
 * 
 * @param character The sentinel substitute.
 * @param dictionary The dictionary.
 */

/*!
 * @fn SentinelRankDictionary#setSentinelPosition
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Sets the sentinel position..
 * 
 * @signature setSentinelPosition(dictionary, pos)
 * 
 * @param pos The sentinel position. Types: @link UnsignedIntegerConcept @endlink
 * @param dictionary The dictionary.
 */

/*!
 * @fn SentinelRankDictionary#createSentinelRankDictionary
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions creates the dictionary structure.
 * 
 * @signature void createSentinelRankDictionary(dictionary, text)
 * 
 * @param text A text to be transfered into a dictionary. Types: @link String @endlink
 * @param dictionary The dictionary. 
 */

/*!
 * @fn SentinelRankDictionary#save
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions saves a dictionary to disk.
 * 
 * @signature save(dictionary, fileName [, openMode])
 * 
 * @param openMode The combination of flags defining how the file should be
 *                 opened.To open a file read-only, write-only or to read and
 *                 write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                 <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *                 <tt>OPEN_CREATE</tt>.To append a file if existing add
 *                 <tt>OPEN_APPEND</tt>.To circumvent problems, files are always
 *                 opened in binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                 OPEN_APPEND</tt>
 * @param dictionary The dictionary. Types: SentinelRankDictionary
 * @param fileName C-style character string containing the file name.
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on success.
 */

/*!
 * @fn SentinelRankDictionary#open
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions loads a dictionary from disk.
 * 
 * @signature open(dictionary, fileName [, openMode])
 * 
 * @param openMode The combination of flags defining how the file should be
 *                 opened.To open a file read-only, write-only or to read and
 *                 write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                 <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *                 <tt>OPEN_CREATE</tt>.To append a file if existing add
 *                 <tt>OPEN_APPEND</tt>.To circumvent problems, files are always
 *                 opened in binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                 OPEN_APPEND</tt>
 * @param dictionary The dictionary. Types: SentinelRankDictionary
 * @param fileName C-style character string containing the file name.
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on success.
 */

/*!
 * @class SparseString
 * 
 * @headerfile seqan/index.h
 * 
 * @brief A string storing only a fraction of the values of the original
 *        string.
 * 
 * @signature SparseString<TValueString, TSpec>
 * 
 * @tparam TSpec The specialisation tag. Default: void.
 * @tparam TValueString The type of the string containing the values. Types: @link String @endlink
 */

/*!
 * @defgroup SparseStringFibres FM Index Fibres
 * 
 * @brief Tag to select a specific fibre of a @link FMIndex @endlink.
 * 
 * @section Remarks
 * 
 * These tags can be used to get @link Fibre Fibres @endlink of a FM index.
 * 
 * @see Fibre
 * @see getFibre
 * 
 * @tag SparseStringFibres#FibreValueString
 * 
 * @brief The String containing the stored values.
 * 
 * @tag SparseStringFibres#FibreIndicatorString
 * 
 * @brief The string storing for each position if a value different from a default value is stored.
 */

/*!
 * @fn SparseString#clear
 *
 * @headerfile seqan/index.h
 * 
 * @brief Clears the @link SparseString @endlink.
 * 
 * @signature clear(sparseString)
 * 
 * @param sparseString The  @link SparseString @endlink to be cleared.
 */

/*!
 * @fn SparseString#empty
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns whether or not the @link SparseString @endlink is empty.
 * 
 * @signature empty(sparseString)
 * 
 * @param sparseString The @link SparseString @endlink to be checked. 
 */

/*!
 * @fn SparseString#getValue
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the value of a @link SparseString @endlink.
 * 
 * @signature getValue(sparseString, pos)
 * 
 * @param sparseString The @link SparseString @endlink. 
 * @param pos The position at which a value should be assign to the sparse string. 
 *        Types: @link UnsignedIntergerConcept @endlink
 *
 * @return TValue The type @link GetValue @endlink of @link SparseString @endlink is returned.
 */

/*!
 * @fn SparseString#value
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the value of a @link SparseString @endlink.
 * 
 * @signature value(sparseString, pos)
 * 
 * @param sparseString The @link SparseString @endlink. 
 * @param pos The position at which a value should be assign to the sparse string. 
 *        Types: @link UnsignedIntergerConcept @endlink
 *
 * @return TValue The type @link Reference @endlink of @link SparseString @endlink is returned.
 */

/*!
 * @fn SparseString#getFibre
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns a specific fibre of a @link SparseString @endlink.
 * 
 * @signature getFibre(sparseString, fibreTag)
 * 
 * @param fibreTag A tag that identifies the @link Fibre @endlink. Types:
 *                 @link SparseStringFibres SparseString Fibres @endlink
 * @param sparseString The sparseString holding the fibre.
 * 
 * @return TReturn A reference to the @link Fibre @endlink object.
 */

/*!
 * @fn SparseString#length
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the number of elements in the @link SparseString @endlink.
 * 
 * @signature TSize length(sparseString)
 * 
 * @param sparseString The sparse string suffix array.
 *
 * @return TSize The number of elements in the sparse string array. Types: The result of @link Size @endlink of the
 * sparse string.
 */

/*!
 * @fn SparseString#resize
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Resets the number of elements in the compressed suffix array.
 * 
 * @signature TSize resize(sparseString, newLenght)
 * 
 * @param sparseString The sparse string.
 * @param newLength The number of elements which should be stored in the  sparse string. Types: @link
 * UnsignedIntegerConcept @endlink.
 *
 * @return TSize The number of elements in the  sparse string. Types: The result of @link Size @endlink of the
 * sparse string.
 *
 * @section Note If the new length is smaller than the actual one then the last <bb>x<bb> items of the compressed suffix array
 * are deleted with x = oldLength - newLength.
 */

/*!
 * @fn SparseString#open
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions open a sparse string from disk.
 * 
 * @signature open(string, fileName [, openMode])
 * 
 * @param openMode The combination of flags defining how the file should be
 *                 opened.To open a file read-only, write-only or to read and
 *                 write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                 <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *                 <tt>OPEN_CREATE</tt>.To append a file if existing add
 *                 <tt>OPEN_APPEND</tt>.To circumvent problems, files are always
 *                 opened in binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                 OPEN_APPEND</tt>
 * @param string The string to be saved. Types: SparseString
 * @param fileName C-style character string containing the file name.
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on success.
 */
/*!
 * @fn SparseString#save
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions saves a sparse string to disk.
 * 
 * @signature open(string, fileName [, openMode])
 * 
 * @param openMode The combination of flags defining how the file should be
 *                 opened.To open a file read-only, write-only or to read and
 *                 write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                 <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *                 <tt>OPEN_CREATE</tt>.To append a file if existing add
 *                 <tt>OPEN_APPEND</tt>.To circumvent problems, files are always
 *                 opened in binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                 OPEN_APPEND</tt>
 * @param string The string to be saved. Types: SparseString
 * @param fileName C-style character string containing the file name.
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on success.
 */
