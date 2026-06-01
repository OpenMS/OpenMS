Overview
======

MS-GF+ (aka MSGF+ or MSGFPlus) performs peptide identification by scoring
MS/MS spectra against peptides derived from a protein sequence database.
It supports the HUPO PSI standard input file (mzML) and saves results in
the mzIdentML format, though results can easily be transformed to TSV.
ProteomeXchange supports Complete data submissions using MS-GF+ search results.

MS-GF+ is optimized for a variety of spectral types, i.e., combinations
of fragmentation method, instrument, enzyme, and experimental protocols.
It supports a variety of input file formats, including mzML, mzXML,
Mascot Generic File (mgf), MS2 files, Micromass Peak List files (pkl),
and Concatenated DTA files (_dta.txt).

Requirements
======

Java Runtime 8 or higher (use 64-bit Java)\
At least 2GB of memory (recommended to use 4GB); larger FASTA files require more memory

There are some issues running with newer version of Java (Java 11 and newer), due to the deprecation and removal of some libraries. 

Downloads / Updates
======

* https://github.com/MSGFPlus/msgfplus/releases

*Version number notes*

As of [January 20, 2016 (commit 375d462)](https://github.com/MSGFPlus/msgfplus/commit/375d462e30cbe460b699091a7d6ba52bc192aba1) the version numbering scheme changed.
Previously the version number was the SVN commit number; git does not have simple commit numbers, so MSGFPlus was changed to a date-based version numbering scheme.

An example: v10282 became v2016.01.20

Installation
======

Unzip MSGFPlus.zip\
Place MSGFPlus.jar in any folder

Usage Information
======

Type `java -jar MSGFPlus.jar` for command line arguments.

To convert an mzid output file into a tsv file, run `java -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv`

Alternatively, use the Mzid-To-Tsv-Converter, which is a faster converter that supports larger result files.
It is a C# application that works on Windows or on Linux using mono.
Download the Mzid-To-Tsv-Converter <a href="https://github.com/PNNL-Comp-Mass-Spec/Mzid-To-Tsv-Converter/releases">from GitHub</a>.

For detailed documentation, see the "docs" subfolder, or visit:
* [GitHub project help pages](https://msgfplus.github.io/msgfplus/)
* [GitHub repo HTML help pages - same as above, but may have issues](https://htmlpreview.github.io/?https://github.com/MSGFPlus/msgfplus/blob/master/docs/index.html)
* (previously at https://bix-lab.ucsd.edu/pages/viewpage.action?pageId=13533355)

Contact Information
======

PNNL Proteomics [proteomics@pnnl.gov]\
Sangtae Kim [sangtae.kim (at) gmail.com]

Publications
======

"MS-GF+ makes progress towards a universal database search tool for proteomics,"\
Sangtae Kim and Pavel A Pevzner,
Nat Commun. 2014 Oct 31; 5:5277. doi: 10.1038/ncomms6277.\
https://pubmed.ncbi.nlm.nih.gov/25358478/

"Spectral Probabilities and Generating Functions of Tandem Mass Spectra: A Strike against Decoy Databases",\
Sangtae Kim, Nitin Gupta, and Pavel A Pevzner,
J Proteome Res. 2008 Aug; 7(8):3354-63. doi: 10.1021/pr8001244.\
https://pubmed.ncbi.nlm.nih.gov/18597511/

Source
======

https://github.com/MSGFPlus/msgfplus
