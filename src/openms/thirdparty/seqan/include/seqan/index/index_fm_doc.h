/**
.Tag.WT
..summary:Tag that specifies the @Spec.FMIndex@ to use a wavelet tree as the occurrence table.
..cat:Index
*/
/**
.Tag.FM Index Fibres
..summary:Tag to select a specific fibre of a @Spec.FMIndex@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a FM index.
..cat:Index

..tag.FibrePrefixSumTable:The prefix sum table of the index.
..tag.FibreSA:The compressed suffix array of the text.
..tag.FibreText:The original text of the index.
..tag.FibreLfTable:The lf table.
..tag.FibreSaLfTable:The lf table as well as the compressed suffix array.
...remarks:This tag can only be used with the functions @Function.indexRequire@ or @Function.indexSupplied@.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index_fm.h
*/
/**
.Tag.CompressText
..cat:Index
..summary:Tag to select a FM index variant that can be used such that it is 
not necessary to store the text after index construction. This index is very
space efficient.
*/
/**
.Spec.FMIndex:
..summary:An index based on the Burrows-Wheeler transform.
..cat:Index
..general:Class.Index
..signature:Index<TText, FMIndex<TOccSpec, TSpec> >
..param.TText:The text type.
...type:Class.String
...type:Class.StringSet
..param.TOccSpec:Occurrence table specialisation. 
...type:Tag.WT
...type:Tag.SBM
...remarks:The tags are really shortcuts for the different @Class.SentinelRankDictionary@s
...default:Tag.WT
..param.TSpec:FM index specialisation.
...type:Tag.CompressText
...default:void
..include:seqan/index.h
*/
/**
.Function.FMIndex#getFibre:
..summary:Returns a specific fibre of a fm index.
..signature:getFibre(index, fibreTag)
..class:Spec.FMIndex
..cat:Index
..param.index:The index holding the fibre.
...type:Spec.FMIndex
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.FM Index Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
*/
/**
.Function.FMIndex#toSuffixPosition
..class:Spec.FMIndex
..summary:This function computes the position of a specified position in the suffix array (additionally containing 
entries for the sentinels. The returned position correspond to the suffix array of the original text without sentinels.
..signature:toSuffixPosition(fmIndex, pos, offset)
..param.fmIndex:The FM index.
...type:Spec.FMIndex
..param.pos:The position in the suffix array of the fm index (with sentinels).
...type:Concept.UnsignedIntegerConcept
..param.offset:The number of sequences in the original text.
...type:Concept.UnsignedIntegerConcept
..include:seqan/index.h
*/
/**
.Function.FMIndex#indexCreate
..summary:Creates a specific @Metafunction.Fibre@.
..signature:indexCreate(index, fibreTag)
..param.index:The index to be created.
...type:Spec.FMIndex
..param.fibreTag:The fibre of the index to be computed.
...type:Tag.FM Index Fibres.tag.FibreSaLfTable
..remarks:If you call this function on the compressed text version of the FM index
you will get an error message: "Logic error. It is not possible to create this index without a text."
*/
/**
.Function.FMIndex#indexSupplied:
..summary:Returns whether a specific @Metafunction.Fibre@ is present.
..param.fibreTag:
...type:Tag.FM Index Fibres
*/
/**
.Tag.CompressedSA Fibres
..summary:Tag to select a specific fibre of a @Class.CompressedSA@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a sparse string.
..cat:Index

..tag.FibreSparseString:The sparse string.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index.h
*/
/**
.Class.CompressedSA:
..cat:Index
..summary:A suffix array storing only a few suffix array entries and computing the remaining on demand.
..signature:CompressedSA<TSparseString, TLfTable, TSpec>
..param.TSparseString:The string containing specific suffix array entries.
...type:Class.SparseString
..param.TLfTable:The lfTable containg an occurrence table and a prefix sum table.
...type:Class.LfTable
..param.TSpec:Possibility to specialize a compressed suffix array.
...default:void.
..remarks:The compressed suffix array can only be used with the FM index.
..include:seqan/index.h
*/
/**
.Function.clear.param.object.type:Class.CompressedSA
.Function.clear.class:Class.CompressedSA
*/
/**
.Function.empty.param.object.type:Class.CompressedSA
.Function.empty.class:Class.CompressedSA
*/
/**
.Function.createCompressedSa
..summary:This functions creates a compressed suffix array with a specified compression factor.
..signature:void createCompressedSa(compressedSA, completeSA, compressionFactor [,offset])
..class:Class.CompressedSA
..param.compressedSA:The compressed suffix array
...type:Class.CompressedSA
..param.completeSA:A complete suffix array containing all values
..param.compressionFactor:The compression factor.
...type:Concept.UnsignedIntegerConcept
...remarks:A compression factor of x means that the compressed suffix array specifically stores a value for every x values in the complete suffix array.
..param:offset:Number of elements at the beginning which should contain the default value.
..include:seqan/index.h
*/
/**
.Function.CompressedSA#getFibre:
..summary:Returns a specific fibre of a container.
..signature:getFibre(compressedSA, fibreTag)
..class:Class.CompressedSA
..cat:Index
..param.compressedSA:The compressed suffix array holding the fibre.
...type:Class.CompressedSA
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.CompressedSA Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
..example.code:
Index< String<char> > index_esa("tobeornottobe");

String<char> & text = getFibre(indexEsa, EsaText());
*/
/**
.Function.length.param.object.type:Class.CompressedSA
.Function.length.class:Class.CompressedSA
*/
/**
.Function.resize.param.object.type:Class.CompressedSA
.Function.resize.class:Class.CompressedSA
*/
/**
.Function.CompressedSA#open
..class:Class.CompressedSA
..summary:This functions opens a compressed suffix array from disk.
..signature:open(compressedSA, fileName [, openMode])
..param.compressedSA:The compresses suffix array to be opened.
...type:Class.CompressedSA
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/template <typename TSparseString, typename TLfTable, typename TSpec>
/**
.Function.CompressedSA#save
..class:Class.CompressedSA
..summary:This functions saves a compressed suffix array to disk.
..signature:save(compressedSA, fileName [, openMode])
..param.compressedSA:The compressed suffix array to be saved.
...type:Class.CompressedSA
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/
/**
.Function.setLfTable
..class:Class.CompressedSA
..summary:Set the LfTable of the compressed suffix array.
..signature:setLfTable(CompressedSA<TSparseString, TLfTable, TSpec> compressedSa, TLfTable & lfTable)
..param.CompressedSA<TSparseString, TLfTable, TSpec>:The compressed suffix array.
...type:Class.CompressedSA
..param.lfTable
...type:Class.LfTable
..include:seqan/index.h
*/
/**
.Function.CompressedSA#value
..class:Class.CompressedSA
..summary:Returns the value stored at a specified position in the compressed suffix-array.
..signature:value(compressedSA, pos)
..param.compressedSA:The compressed suffix array to access.
...type:Class.CompressedSA
..param.pos:Position at which to access the suffix array.
...type:Concept.UnsignedIntegerConcept
..remarks:Note that the compressed suffix array is read only. Therefore a const reference is return by
this function.
*/
///.Function.begin.param.object.type:Class.CompressedSA
///.Function.begin.class:Class.CompressedSA
///.Function.end.param.object.type:Class.CompressedSA
///.Function.end.class:Class.CompressedSA
/**
.Tag.LF Table Fibres
..summary:Tag to select a specific fibre of a @Spec.FMIndex@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a FM index.
..cat:Index

..tag.FibreOccTable:The occurrence table of the lf table.
..tag.FMTablePrefixSumTable:The prefix sum table of the lf table.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index_fm.h
*/
/**
.Class.LfTable:
..cat:Index
..summary:LfTable is an object storing all necessary information for the LF-mapping.
..signature:LfTable<TOccTable, TPrefixSumTable>
..param.TOccTable:The occurrence table data structure.
...type:Class.SentinelRankDictionary
..param.TPrefixSumTable:The specialisation tag.
...default:String
..include:seqan/Index.h
*/
/**
.Function.LfTable#clear
..class:Class.LfTable
..summary:Clears the LF table.
..signature:clear(lfTable)
..param.lfTable:The LF table to be cleared.
...type:Class.LfTable
..include:seqan/index.h
*/
/**
.Function.LfTable#empty
..class:Class.LfTable
..summary:Clears the LF table.
..signature:empty(lfTable)
..param.lfTable:The LF table to be cleared.
...type:Class.LfTable
..returns:$true$ if the LF table is empty, $false$ otherwise.
...type:nolink:$bool$
..include:seqan/index.h
*/
/**
.Function.createLfTable
..class:Class.LfTable
..summary:Creates the LF table
..signature:createLfTable(lfTable, text)
..param.lfTable:The LF table to be constructed.
...type:Class.LfTable.
..param.text:The underlying text
...type:Class.String
..returns:$true$ on successes, $false$ otherwise.
..include:seqan/index.h
*/
/**
.Function.LfTable#getFibre:
..summary:Returns a specific fibre of a container.
..signature:getFibre(container, fibreTag)
..class:Class.LfTable
..cat:Index
..param.container:The container holding the fibre.
...type:Class.LfTable
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.LF Table Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
*/
/**
.Function.lfMapping:
..class:Class.LfTable
..summary:Returns the position of an character at a specified position of L in F. L corresponds to the last column of 
the sorted cyclic rotations of the original text, while F correspond to the first column.
..cat:Index
..signature:lfMapping(lfTable, pos)
..param.lfTable:The @Class.LfTable@ holding the occurrence and prefix sum table.
...type:Class.LfTable
..param.pos:The position in L
..returns:Returns the position of the character L[c] in F. The returned position is of the same type as pos.
..include:seqan/index.h
*/
/**
.Function.LfTable#open
..class:Class.LfTable
..summary:This functions loads a LF table from disk.
..signature:open(lfTable, fileName [, openMode])
..param.lfTable:The lfTable.
...type:Class.LfTable
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A nolink:$bool$ which is $true$ on success.
..include:seqan/index.h
*/
/**
.Function.LfTable#save
..class:Class.LfTable
..summary:This functions saves a LF table to disk.
..signature:save(lfTable, fileName [, openMode])
..param.lfTable:The dictionary.
...type:Class.LfTable
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A nolink:$bool$ which is $true$ on success.
..include:seqan/index.h
*/
/**
.Tag.PrefixSumTable Fibres
..summary:Tag to select a specific fibre of a @Class.CompressedSA@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a sparse string.
..cat:Index

..tag.FibreEntries:The entries in the prefix sum table.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index.h
*/
/**
.Class.PrefixSumTable:
..cat:Index
..summary:The prefix-sum table is a data structure which stores for each character the number of smaller lexicographic 
smaller characters in a given text.
..signature:PrefixSumTable<TChar, TSpec>
..param.TChar:The character type
..param.TSpec:A specialisation tag.
...default:void
..include:seqan/Index.h
*/
/**
.Function.PrefixSumTable#clear
..class:Class.PrefixSumTable
..summary:Clears the prefix sum table.
..signature:clear(prefixSumTable)
..param.prefixSumTable:The prefix sum table to be cleared.
...type:Class.LfTable
..include:seqan/index.h
*/
/**
.Function.createPrefixSumTable
..class:Class.PrefixSumTable
..summary:Creates the prefix sum table
..signature:createPrefixSumTable(prefixSumTable, text)
..param.prefixSumTable:The prefix sum table to be constructed.
...type:Class.PrefixSumTable
..param.text:The underlying text.
...type:Class.String
...type:Class.StringSet
..include:seqan/index.h
*/
/**
.Function.PrefixSumTable#getAlphabetSize
..class:Class.PrefixSumTable
..summary:Returns the number of different characters in the prefix sum table.
..signature:getAlphabetSize(prefixSumTable)
..param.prefixSumTable:A prefix sum table.
...type:Class.PrefixSumTable.
..include:seqan/index.h
*/
/**
.Function.PrefixSumTable#getCharacterPosition
..class:Class.PrefixSumTable
..summary:Returns the position of a given character within the prefix sum table.
..signature:getCharacterPosition(prefixSumTable, character)
..param.prefixSumTable:A prefix sum table.
...type:Class.PrefixSumTable
..param.character:A character.
..include:seqan/index.h
*/
/**
.Function.PrefixSumTable#getCharacter
..class:Class.PrefixSumTable
..summary:Returns the character of a given position within the prefix sum table.
..signature:getCharacter(prefixSumTable, pos)
..param.prefixSumTable:A prefix sum table.
...type:Class.PrefixSumTable.
..param.pos:A position.
...type:Concept.UnsignedIntegerConcept
..include:seqan/index.h
*/
/**
.Function#getPrefixSum
..class:Class.PrefixSumTable
..summary:Returns the prefix sum of a given position. 
..signature:getPrefixSum(prefixSumTable, pos)
..param.prefixSumTable:A prefix sum table.
...type:Class.PrefixSumTable
..param.pos:A position.
...type:Concept.UnsignedIntegerConcept
..include:seqan/index.h
*/
/**
.Function.PrefixSumTable#getValue
..class:Class.PrefixSumTable
..summary:Returns the prefix sum of a given position. 
..signature:getValue(prefixSumTable, pos)
..param.prefixSumTable:A prefix sum table.
...type:Class.PrefixSumTable
..param.pos:A position.
...type:Concept.UnsignedIntegerConcept
..include:seqan/index.h
*/
/**
.Function.PrefixSumTable#getFibre:
..class:Class.PrefixSumTable
..summary:Returns a specific fibre of a prefix-sum table.
..signature:getFibre(prefixSumTable, fibreTag)
..class:Class.PrefixSumTable
..cat:Index
..param.prefixSumTable:The container holding the fibre.
...type:Class.PrefixSumTable
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.PrefixSumTable Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
*/
/**
.Function.PrefixSumTable#length
..class:Class.PrefixSumTable
..summary:Returns the number of different characters in the prefix-sum table.
..signature:length(lfTable)
..param.lfTable:The prefix-sum table.
...type:Class.LfTable
..returns:Returns the number of different characters in the prefix-sum table.
...type:Metafunction.Size
...remarks:If the type of the characters of the prefix-sum table consists of more than 8 bit only the characters
actually occurring in the original text are accounted for when calling length.
..include:seqan/index.h
*/
/**
.Function.PrefixSumTable#prefixSum
..class:Class.PrefixSumTable
..summary:Returns a reference to the entry of the prefix sum table of a given position. 
..signature:prefixSum(prefixSumTable, pos)
..param.prefixSumTable:A prefix sum table.
...type:Class.PrefixSumTable
..param.pos:A position.
...type:Concept.UnsignedIntegerConcept
..include:seqan/index.h
*/
/**
.Function.PrefixSumTable#resize
..class:Class.PrefixSumTable
..summary:Resize the prefix sum table to be able to store more or less characters. 
..signature:resize(prefixSumTable, size [,value, resizeTag])
..param.prefixSumTable:A prefix sum table.
...type:Class.PrefixSumTable
..param.size:The new size.
...type:Concept.UnsignedIntegerConcept
..param.resizeTag: Specifies the strategy that is applied if the capacity of $object$ is less than $newLength$. (optional)
...type:Tag.Overflow Strategy
...default:Specified by @Metafunction.DefaultOverflowExplicit@.

..param.value:The value to be used to initialize the new storage.
..include:seqan/index.h
*/
/**
.Function.setPrefixSum
..class:Class.PrefixSumTable
..summary:Returns a reference to the entry of the prefix-sum table of a given position. 
..signature:setPrefixSum(prefixSumTable, value, pos)
..param.prefixSumTable:A prefix sum table.
...type:Class.PrefixSumTable
..param.value:A specified value to be inserted.
..param.pos:A position.
...type:Concept.UnsignedIntegerConcept
..include:seqan/index.h
*/
/**
.Function.PrefixSumTable#open
..class:Class.PrefixSumTable
..summary:This functions loads a prefix-sum table from disk.
..signature:open(prefixSumTable, fileName [, openMode])
..param.prefixSumTable:The prefix-sum table.
...type:Class.PrefixSumTable
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/
/**
.Function.PrefixSumTable#save
..class:Class.PrefixSumTable
..summary:This functions saves a prefix-sum table to disk.
..signature:save(prefixSumTable, fileName [, openMode])
..param.prefixSumTable:The prefix-sum table.
...type:Class.PrefixSumTable
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/
/**
.Tag.SBM
..summary:Tag that specifies the @Spec.FMIndex@ to use a StringSet of rank support bis strings as the occurrence table.
..cat:Index
*/
/**
.Spec.SequenceBitMask Fibres
..cat:Index
..summary:Tag to select a specific fibre of a SequenceBitMask.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a SequenceBitMask.

..DISABLED.tag.FibreBitStrings:The string set containing a bit string for each node.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index.h
*/
/**
.Spec.SequenceBitMask:
..cat:Index
..general:Class.RankDictionary
..summary:The string set bit string dictionary is a string set of rank support bit strings for constant time acces
of the rank of a specified character at a specified position.
..signature:SequenceBitMask<TValue>
..param.TValue:The value type of the .
..include:seqan/index.h
..remarks:This data structure is optimized for very small alphabets, such as @Spec.Dna@ or @Spec.Dna5@. Consider using a @Spec.WaveletTree@ if your alphabet size is larger.
*/
///.Function.RankDictionary#getFibre.param.fibreTag.type:Spec.SequenceBitMask Fibres
/**
.Tag.WaveletTree Fibres
..summary:Tag to select a specific fibre (e.g. table, object, ...) of a @Spec.WaveletTree@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a @Spec.WaveletTree@.
..cat:Spec.WaveletTree

..tag.FibreBitStrings:The string set containing a bit string for each node.

..tag.FibreTreeStructure:The wavelet tree structure of the wavelet tree.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index.h
*/
/**
.Class.RankDictionary:
..cat:Index
..summary:A rank dictionary is a data structure to store the rank of an element of a sequence at every position of the 
sequence.
..signature:RankDictionary<TSpec>
..param.TSpec:The rank dictionary specialisation.
...type:Spec.WaveletTree
...type:Spec.SequenceBitMask
...default:@Spec.WaveletTree@
..include:seqan/index.h
*/
/**
.Spec.WaveletTree:
..general:Class.RankDictionary
..cat:Index
..summary:A wavelet tree is a tree like binary encoding of a text.
..signature:WaveletTree<TValue>
..param.TValue:The value type of the wavelet tree.
..include:seqan/index.h
..remarks:The nodes of a wavelet tree consist of a bit string as well as a character c. In each level of the tree, 
characters smaller than c are represented as a 0 while character greater or equal to c are represented with a 1.
The characters represented by a 0 form the string to be represented by the left subtree while characters represented
by a 1 form the string of the right subtree. Therefore, only the bit string of the root node represents all characters while all other nodes represent subsets.
*/
/**
.Function.RankDictionary#clear
..class:Class.RankDictionary
..summary:Clears the rank dictionary.
..signature:clear(dictionary)
..param.dictionary:The rank dictionary to be cleared.
...type:Class.RankDictionary
..include:seqan/index.h
*/
/**
.Function.RankDictionary#empty
..class:Class.RankDictionary
..summary:Returns whether or not the rank dictionary is empty.
..signature:empty(dictionary)
..param.dictionary:The rank dictionary to be checked.
...type:Class.RankDictionary
..returns:$true$ if the dictionary is empty, $false$ otherwise.
..include:seqan/index.h
*/
/**
.Function.RankDictionary#getValue
..summary:Returns the character of a specified position.
..signature:getValue(dictionary, pos)
..class:Class.RankDictionary
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..param.pos:The position
..include:seqan/index.h
*/
/**
.Function.RankDictionary#getFibre:
..summary:Returns a specific fibre of a dictionary.
..signature:getFibre(dictionary, fibreTag)
..class:Class.RankDictionary
..cat:Index
..param.dictionary:The dictionary holding the fibre.
...type:Class.RankDictionary
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.WaveletTree Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
*/
/**
.Function.RankDictionary#countOccurrences:
..summary:Returns the rank (number of occurrences) of a specified character up to a specified position. 
..signature:countOccurrences(dictionary, character, pos)
..class:Class.RankDictionary
..cat:Index
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..param.character:The character of interest.
..param.pos:The position (which is also included in the rank computation).
..returns:The rank (number of occurrences) of a specified character up to a specified position. 
...type:nolink:$unsigned$
..include:seqan/index.h
*/
/**
.Function.createRankDictionary
..class:Class.RankDictionary
..summary:This functions creates the dictionary.
..signature:createRankDictionary(dictionary, text)
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..param.text:A text to be transfered into a wavelet tree.
...type:Class.String
..include:seqan/index.h
*/
/**
.Function.RankDictionary#open
..class:Class.RankDictionary
..summary:This functions loads a dictionary from disk.
..signature:open(dictionary, fileName [, openMode])
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/
/**
.Function.RankDictionary#save
..class:Class.RankDictionary
..summary:This functions saves a dictionary to disk.
..signature:open(dictionary, fileName [, openMode])
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/
///.Metafunction.Fibre.param.TContainer.type:Class.RankDictionary
///.Metafunction.Fibre.param.TSpec.type:Tag.WaveletTree Fibres
/**
.Tag.RankSupportBitString Fibres
..summary:Tag to select a specific fibre (e.g. table, object, ...) of a @Class.RankSupportBitString@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a rank support bit string.
..cat:Index

..tag.FibreBits:The bit string.
..tag.FibreBlocks:The block string.
..tag.FibreSuperBlocks:The super block string.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index.h
*/
/**
.Class.RankSupportBitString:
..summary:A bit string supporting rank queries in constant time.
..cat:Index
..signature:RankSupportBitString<TSpec>
..param.TSpec:Specialisation tag.
...default:void
..remarks:The constant rank query time is achieved by evaluating precomputed subsolutions. In order to do so, the bit string is divided into blocks of length l. A super block string stores for each block of l blocks the number of bits set from the beginning. In addition a block string stores the number of bits set in each block from the start of the last super block block. Therefore it is possible to compute the result of a rank query in constant time by adding information from the bit, block and super block string.
..include:seqan/index.h
*/
/**
.Function.RankSupportBitString#appendValue:
..signature:appendValue(ranKSupportBitString, value)
..cat:Index
..summary:Appends a value to a container.
..param.ranKSupportBitString:A @Class.RankSupportBitString@.
...type:Class.RankSupportBitString
..param.value:Value that is appended to $target$.
...type:Concept.UnsignedIntegerConcept
...type:nolink:bool
...remarks:If the value is different from 0 it is interpreted as 1.
..include:seqan/sequence.h
*/
/**
.Function.clear.param.object.type:Class.RankSupportBitString
.Function.clear.class:Class.RankSupportBitString
*/
/**
.Function.getRank
..class:Class.RankSupportBitString
..summary:Returns the rank (the number of bits set from the start of the bit string) of a specified position.
..signature:getRank(bitString, pos)
..param.bitString:The bit string.
...type:Class.RankSupportBitString
..param.pos:Position of a bit.
..returns:Value type of the super block fibre (default unsigned long).
..include:seqan/index.h
..example.code:
*/
/**
.Function.empty.param.object.type:Class.RankSupportBitString
.Function.empty:Class.RankSupportBitString
*/
/**
.Function.isSetBit
..class:Class.RankSupportBitString
..summary:Returns whether the bit with the given index is set to 1.
..signature:isSetBit(bitString, pos)
..param.bitString:The bit string.
...type:Class.RankSupportBitString
..param.pos:Position of the bit.
..returns:Returns whether a specified bit is set or not.
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";

RankSupportBitString<> bitString;
resize(bitString, length(genome));
...
mark all 'a's
...

for (unsigned i = 0; i < length(bitString); ++i)
    if(isSetBit(bitString, i))
        std::cout << "a found at: " << i << std::endl;
*/
/**
.Function.RankSupportBitString#getFibre:
..class:Class.RankSupportBitString
..summary:Returns a specific fibre of a container.
..signature:getFibre(container, fibreTag)
..class:Class.RankSupportBitString
..cat:Index
..param.container:The container holding the fibre.
...type:Class.RankSupportBitString
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.RankSupportBitString Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
..example.code:
Index< String<char> > index_esa("tobeornottobe");

String<char> & text = getFibre(indexEsa, EsaText());
*/
/**
.Function.length.param.object.type:Class.RankSupportBitString
.Function.length.class:Class.RankSupportBitString
*/
/**
.Function.resize.param.object.type:Class.RankSupportBitString
.Function.resize.class:Class.RankSupportBitString
*/
/**
.Function.setBitTo
..class:Class.RankSupportBitString
..signature:setBitTo(bitString, pos, bit)
..param.bitString:The bit string.
...type:Class.RankSupportBitString
..param.pos:Position of the bit.
..param.bit:The value of the bit.
...remarks:Note that values different from 0 are interpreted as 1.
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";

RankSupportBitString<> bitString;
resize(bitString, length(genome));

for (unsigned i = 0; i < length(genome); ++i)
    if(genome[i] < Dna5('c'))
        setBitTo(bitString, 1);

updateRanks_(bitString);
*/
/**
.Function.RankSupportBitString#open
..class:Class.RankSupportBitString
..summary:This functions saves a @Class.RankSupportBitString@ to disk.
..signature:open(bitString, fileName [, openMode])
..param.bitString:The bit string to be saved.
...type:Class.RankSupportBitString
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/
/**
.Function.RankSupportBitString#save
..class:Class.RankSupportBitString
..summary:This functions saves a @Class.RankSupportBitString@ to disk.
..signature:save(bitString, fileName [, openMode])
..param.bitString:The bit string to be saved.
...type:Class.RankSupportBitString
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/
/**
.Tag.RightArrayBinaryTree Fibres
..summary:Tag to select a specific fibre (e.g. table, object, ...) of a @Class.RightArrayBinaryTree@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a RightArrayBinaryTree.
..cat:RightArrayBinaryTree
..tag.FibreTreeStructureEncoding:The string encoding the wavelet tree structure.
..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index.h
*/
/**
.Class.RightArrayBinaryTree:
..cat:WaveletTree
..summary:A special format to encode the structure of a wavelet tree. The structure is very space efficient because only one position is stored which encodes where the left and right subtree of a given node exist.
..signature:RightArrayBinaryTree<TValue, TSpec>
..param.TSpec:The value type, that is the type of the stored characters.
..param.TSpec:The wavelet tree structure specialisation.
...default:void.
..include:seqan/index.h
*/
/**
.Function.clear.param.object.type:Class.RightArrayBinaryTree
.Function.clear.class:Class.RightArrayBinaryTree
*/
/**
.Function.createRightArrayBinaryTree
..class:Class.RightArrayBinaryTree
..summary:Computes the wavelet tree structure of a text.
..signature:createRightArrayBinaryTree(waveletTreeStructure, text)
..param.waveletTreeStructure:A wavelet tree structure.
...type:Class.RightArrayBinaryTree
..param.text:A text.
...type:Class.String
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";

RightArrayBinaryTree<Dna5> waveletTreeStructure;
computeRightArrayBinaryTree(genome);
*/
/**
.Function.empty.param.object.type:Class.RightArrayBinaryTree
.Function.empty.class:Class.RightArrayBinaryTree
*/
/**
.Function.RightArrayBinaryTree#getFibre:
..class:Class.RightArrayBinaryTree
..summary:Returns a specific fibre of a container.
..signature:getFibre(container, fibreTag)
..class:Class.RightArrayBinaryTree
..cat:Index
..param.container:The container holding the fibre.
...type:Class.RightArrayBinaryTree
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.RightArrayBinaryTree Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
*/
/**
.Function.RightArrayBinaryTree#open
..class:Class.RightArrayBinaryTree
..summary:This functions loads a @Class.RightArrayBinaryTree@ from disk.
..signature:open(rightArrayBinaryTree, fileName [, openMode])
..param.rightArrayBinaryTree:The rightArrayBinaryTree.
...type:Class.RightArrayBinaryTree
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/
/**
.Function.RightArrayBinaryTree#save
..class:Class.RightArrayBinaryTree
..summary:This functions saves a @Class.RightArrayBinaryTree@ to disk.
..signature:save(rightArrayBinaryTree, fileName [, openMode])
..param.rightArrayBinaryTree:The rightArrayBinaryTree.
...type:Class.RightArrayBinaryTree
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/
///.Metafunction.Fibre.param.TContainer.type:Class.RightArrayBinaryTree
///.Metafunction.Fibre.param.TSpec.type:Tag.RightArrayBinaryTree Fibres
/**
.Spec.RightArrayBinaryTree Iterator:
..summary:An iterator for @Class.RightArrayBinaryTree@.
..cat:Iter
..general:Class.Iter
..signature:Iter<RightArrayBinaryTree, TSpec >
..param.RightArrayBinaryTree:The @Class.RightArrayBinaryTree@.
...type:Class.RightArrayBinaryTree
..param.TSpec:Specialisation Tag.
...type:Spec.TopDown Iterator
..include:seqan/index.h
*/
/**
.Function.getCharacter
..class:Spec.RightArrayBinaryTree Iterator
..summary:This function returns the pivot character of the node the iterator currently points to.
..signature:getCharacter(it)
..param.it:The iterator.
...type:Spec.RightArrayBinaryTree Iterator
..include:seqan/index.h
*/
/**
.Function.getLeftChildPos
..class:Spec.RightArrayBinaryTree Iterator
..summary:Returns the position in @Class.RightArrayBinaryTree@ of the left child vertex.
..signature:getLeftChildPos(it)
..param.it:The iterator.
...type:Spec.RightArrayBinaryTree Iterator.
..include:seqan/index.h
*/
/**
.Function.getSubTreeSize
..class:Spec.RightArrayBinaryTree Iterator
..summary:Returns the number of vertices in the subtree starting at the position an iterator points to.
..signature:getSubTreeSize(it)
..param.it:The iterator.
...type:Spec.RightArrayBinaryTree Iterator.
..include:seqan/index.h
*/
/**
.Function.RightArrayBinaryTree Iterator#getPosition
..class:Spec.RightArrayBinaryTree Iterator
..summary:Returns the position of the iterator in the host.
..signature:getPosition(it)
..param.it:The iterator.
...type:Spec.RightArrayBinaryTree Iterator.
..include:seqan/index.h
*/
/**
.Function.getRightChildPos
..class:Spec.RightArrayBinaryTree Iterator
..summary:Returns the position in @Class.RightArrayBinaryTree@ of the right child vertex.
..signature:getLeftChildPos(it)
..param.it:The iterator.
...type:Spec.RightArrayBinaryTree Iterator.
..include:seqan/index.h
*/
/**
.Function.goLeftChild
..class:Spec.RightArrayBinaryTree Iterator
..summary:Sets the iterator to the left child of the current node if it exists and returns true, otherwise the iterator does not change position and the function returns false.
..signature:bool goLeftChild(iterator)
..param.iterator:An iterator of a wavelet tree structure.
...type:Spec.RightArrayBinaryTree Iterator
..remarks:$goLeftChild(iterator)$ goes down the left edge if it exist.
..returns:$true$ if the edge or path to go down exists, otherwise $false$.
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";
RightArrayBinaryTree<Dna5> waveletTreeStructure(genome);

Iterator<RightArrayBinaryTree<Dna5>, TopDown<> >::Type it;
it = begin(waveletTreeStructure); // go to root node

goLeftChild(it); // go to left child of root node
*/
/**
.Function.RightArrayBinaryTree Iterator#goRight
..class:Spec.RightArrayBinaryTree Iterator
..summary:Iterates to the next sibling in a tree.
..signature:bool goRight(iterator)
..param.iterator:
...type:Spec.RightArrayBinaryTree Iterator
*/
/**
.Function.goRightChild
..class:Spec.RightArrayBinaryTree Iterator
..summary:Sets the iterator to the right child of the current node if it exists and returns true, otherwise the iterator does not change position and the function returns false.
..signature:bool goRightChild(iterator)
..param.iterator:An iterator of a wavelet tree structure.
...type:Spec.RightArrayBinaryTree Iterator
..remarks:$goRightChild(iterator)$ goes down the right edge if it exist.
..returns:$true$ if the edge or path to go down exists, otherwise $false$.
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";
RightArrayBinaryTree<Dna5> waveletTreeStructure(genome);

Iterator<RightArrayBinaryTree<Dna5>, TopDown<> >::Type it;
it = begin(waveletTreeStructure); // go to root node

goRightChild(it); // go to right child of root node
*/
/**
.Function.goUp.param.iterator.type:Spec.TopDownHistory Iterator
.Function.goUp.class:Spec.RightArrayBinaryTree Iterator
*/
/**
.Function.isLeaf
..class:Spec.RightArrayBinaryTree Iterator
..param.iterator.type:Spec.RightArrayBinaryTree Iterator
*/
/**
.Function.RightArrayBinaryTree Iterator#setCharacter
..class:Spec.RightArrayBinaryTree Iterator
..signature:bool setCharacter(iterator, character)
..param.iterator:An iterator of a wavelet tree structure.
...type:Spec.RightArrayBinaryTree Iterator
..param.character:The character to be assigned to a node.
..summary:$setCharacter(iterator, character)$ sets the character of the node the iterator points to to character.
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";
RightArrayBinaryTree<Dna5> waveletTreeStructure(genome);

Iterator<RightArrayBinaryTree<Dna5>, TopDown<> >::Type it;
it = begin(waveletTreeStructure); // go to root node

goRightChild(it); // go to right child of root node
setCharacter(it,'T'); // sets the character of the root's
                      // right child to 'T'
*/
/**
.Function#isRoot
..class:Spec.RightArrayBinaryTree Iterator
..param.iterator.type:Spec.RightArrayBinaryTree Iterator
*/
///.Function.begin.param.object.type:Class.RightArrayBinaryTree
///.Function.begin.class:Spec.RightArrayBinaryTree Iterator
///.Function.container.param.iterator.type:Class.RightArrayBinaryTree
///.Function.container.class:Spec.RightArrayBinaryTree Iterator
///.Function.end.param.object.type:Class.RightArrayBinaryTree
///.Function.end.class:Spec.RightArrayBinaryTree Iterator
///.Function.goDown.param.iterator.type:Spec.RightArrayBinaryTree Iterator
///.Function.goDown.class:Spec.RightArrayBinaryTree Iterator
/**
.Tag.SentinelRankDictionary Fibres
..summary:Tag to select a specific fibre (e.g. table, object, ...) of a @Class.SentinelRankDictionary@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a SentinelRankDictionary.
..cat:SentinelRankDictionary

..tag.FibreRankDictionary:The rank dictionary.

..tag.FibreSentinelPosition:The bit string encoding the position of the sentinel sign.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index.h
*/
/**
.Class.SentinelRankDictionary:
..cat:Index
..summary:A rank dictionary, additional storing sentinel character which are not accounted for in a rank querry.
..signature:SentinelRankDictionary<TRankDictionary, TSpec>
..param.TRankDictionary:The rank dictionary of a text.
...type:Class.RankDictionary
..param.TSpec:Specialisation
..include:seqan/index.h
*/
/**
.Function.SentinelRankDictionary#clear
..class:Class.SentinelRankDictionary
..summary:Clears the dictionary.
..signature:clear(dictionary)
..param.dictionary:The rank dictionary to be cleared.
...type:Class.SentinelRankDictionary
..include:seqan/index.h
*/
/**
.Function.sentinelPosition
..class:Class.SentinelRankDictionary
..summary:Returns whether a specified position is a sentinel position.
..signature:sentinelPosition(dictionary, pos)
..param.dictionary:The dictionary.
...type:Class.SentinelRankDictionary
..param.pos:The position.
...type:Concept.UnsignedIntegerConcept
..include:seqan/index.h
*/
/**
.Function.SentinelRankDictionary#empty
..class:Class.SentinelRankDictionary
..summary:Returns whether or not the dictionary is empty.
..signature:empty(dictionary)
..param.dictionary:The rank dictionary to be checked.
...type:Class.SentinelRankDictionary
..include:seqan/index.h
*/
/**
.Function.SentinelRankDictionary#getValue
..class:Class.SentinelRankDictionary
..summary:Returns the character of a specified position.
..signature:getCharacter(dictionary, pos)
..param.dictionary:The rank dictionary.
...type:Class.RankDictionary
..param.pos:The position
..include:seqan/index.h
..example.code:
*/
/**
.Function.SentinelRankDictionary#getFibre:
..class:Class.SentinelRankDictionary
..summary:Returns a specific fibre of a dictionary.
..signature:getFibre(dictionary, fibreTag)
..class:Class.RankDictionary
..cat:Index
..param.dictionary:The dictionary holding the fibre.
...type:Class.RankDictionary
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.SentinelRankDictionary Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
*/
/**
.Function.SentinelRankDictionary#countOccurrences
..class:Class.SentinelRankDictionary
..summary:Returns the number of occurrences of a specified character from the start
to a specified position.
..signature:countOccurrences(dictionary, character, pos)
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..param.character:The character.
..param.pos:The position (which is included in the counting).
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";
RankDictionary<String<Dna5> > dictionary(genome);

std::cout << countOccurrences(dictionary, 'a', 3) << std::endl; // 1
std::cout << countOccurrences(dictionary, 'a', 4) << std::endl; // 2
*/
/**
.Function.SentinelRankDictionary#getSentinelSubstitute
..class:Class.SentinelRankDictionary
..summary:Returns the character used to substitute the sentinel sign.
..signature:getSentinelSubstitute(dictionary)
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..include:seqan/index.h
*/
/**
.Function.SentinelRankDictionary#setSentinelSubstitute
..class:Class.SentinelRankDictionary
..summary:Sets the character used to substitute the sentinel sign.
..signature:setSentinelSubstitute(dictionary, character)
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..param.character:The sentinel substitute.
..include:seqan/index.h
*/
/**
.Function.SentinelRankDictionary#setSentinelPosition
..class:Class.SentinelRankDictionary
..summary:Sets the sentinel position..
..signature:setSentinelPosition(dictionary, pos)
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..param.pos:The sentinel position.
..include:seqan/index.h
*/
/**
.Function.SentinelRankDictionary#createSentinelRankDictionary
..class:Class.SentinelRankDictionary
..summary:This functions creates the dictionary structure.
..signature:void createSentinelRankDictionary(dictionary, text)
..param.dictionary:The dictionary.
...type:Class.RankDictionary.
..param.text:A text to be transfered into a dictionary.
...type:Class.String
..include:seqan/index.h
*/
/**
.Function.SentinelRankDictionary#open
..class:Class.SentinelRankDictionary
..summary:This functions loads a dictionary from disk.
..signature:open(dictionary, fileName [, openMode])
..param.dictionary:The dictionary.
...type:Class.SentinelRankDictionary
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/
/**
.Function.SentinelRankDictionary#save
..class:Class.SentinelRankDictionary
..summary:This functions saves a dictionary to disk.
..signature:save(dictionary, fileName [, openMode])
..param.dictionary:The dictionary.
...type:Class.SentinelRankDictionary
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/
///.Metafunction.Fibre.param.TSpec.type:Tag.SentinelRankDictionary Fibres
/**
.Tag.SparseString Fibres
..summary:Tag to select a specific fibre of a @Class.SparseString@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a sparse string.
..cat:Index

..tag.FibreValueString:The string storing values different from a default value.
..tag.FibreIndicatorString:The bit string indicating if the value of a given position is different from the default
value.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index.h
*/
/**
.Class.SparseString:
..cat:Index
..summary:A string storing only a fraction of the values of the original string..
..signature:SparseString<TValueString, TSpec>
..param.TValueString:The string containing the values.
..param.TSpec:The specialisation tag.
...default:void.
..include:seqan/String.h
*/
/**
.Function.SparseString#getFibre:
..summary:Returns a specific fibre of a container.
..signature:getFibre(sparseString, fibreTag)
..class:Class.SparseString
..cat:Index
..param.sparseString:The compressed suffix array holding the fibre.
...type:Class.SparseString
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.SparseString Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
*/
/**
.Function.SparseString#open
..class:Class.SparseString
..summary:This functions saves a sparse string to disk.
..signature:open(string, fileName [, openMode])
..param.string:The string to be saved.
...type:Class.SparseString
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/
/**
.Function.SparseString#save
..class:Class.SparseString
..summary:This functions saves a sparse string to disk.
..signature:save(string, fileName [, openMode])
..param.string:The string to be saved.
...type:Class.SparseString
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/
///.Function.clear.param.object.type:Class.SparseString
///.Function.clear.class:Class.SparseString
///.Function.empty.param.object.type:Class.SparseString
///.Function.empty.class:Class.SparseString
///.Function.assignValue.param.container.type:Class.SparseString
///.Function.assignValue.class:Class.SparseString
///.Function.getValue.param.container.type:Class.SparseString
///.Function.getValue.class:Class.SparseString
///.Function.value.param.container.type:Class.SparseString
///.Function.value.class:Class.SparseString
///.Function.length.param.object.type:Class.SparseString
///.Function.length.class:Class.SparseString
///.Function.resize.param.object.type:Class.SparseString
///.Function.resize.class:Class.SparseString
///.Spec.VSTree Iterator.param.TContainer.type:Spec.FMIndex
///.Function.begin.param.object.type:Spec.FMIndex
