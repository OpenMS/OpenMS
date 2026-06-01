2023/03/31

Comet version "2023.01 rev. 2".
https://uwpr.github.io/Comet/releases/release_202301.html

release 2023.01 rev. 2 (2023.01.2), release date 2023/03/30

mzML/mzXML files without the optional scan index would not be searched because their spectra could not be read. Support for this functionality was implemented in v2019.01.4 thru v2022.01.1 but was lost as of v2022.01.2 when MSToolkit code was updated from that library’s repository. This functionality is re-implemented in this release. Thanks to J. Wang for reporting the issue.
A parameter “export_additional_pepxml_scores” has been implemented. This is an optional/hidden parameter in that it is not written by default in the comet.params files; it needs to be added manually. When this parameter is present and it’s value set to “1”, additional search scores (lnrSp, deltLCn, lnExpect, and IonFrac) are reported in the pep.xml output. Feature request by J. Scheid, OpenMS group.
release 2023.01 rev. 1 (2023.01.1), release date 2023/03/15

For pep.xml output when “num_output_lines = 1”, the deltaCn scores will always be reported as “1.0”. This occurs when only the top hit is reported and has been this behavior since the original release of Comet. This update will correctly report the deltaCn value for this case. Thanks to J. Scheid for reporting this issue.
Comet and MSToolkit had input file name/path limits of 512 and 256 characters, respectively. Any input file strings longer than 256 would cause an error in reading the input spectra data. The file name buffer has been expanded to 4096 to mitigate this issue. Thanks to B. Connolly for reporting this issue.
The Linux/Ubuntu and Mac runners for compiling the release binaries via GitHub Actions have been changed from “ubuntu-latest” and “macos-latest” to “ubuntu-20.04” and “macos-11”. This is to give these binaries wider, backwards compatibility with older OS’s. Thanks to T. Sachsenberg for reporting this issue.
release 2023.01 rev. 0 (2023.01.0), release date 2023/01/31

Address issue where a peptide with a low xcorr identified in a very poor/sparse spectrum would be assigned a good E-value. This is due to the majority of matched xcorr scores being “0” so any poor scoring match looks like a good outlier. Thanks to D. Shteynberg for reporting the issue.
Add “scale_fragmentNL” parameter entry which scales (multiplies) the neutral loss mass value by the number of modified residues in the fragment. Feature requested by A. Keller.
Add contributions of fragment neutral loss peaks in preliminary (Sp) score; previously they only applied to the cross-correlation score.
Fix bug where the fragment neutral loss peak was not analyzed if the primary fragment peak was not matched.
Fix bug where Comet files to analyze a variable modification if more than 19 residues are specified for that mod. Thanks to D. Tabb for reporting the issue.
Fix minor typo in command line help. Thanks to M. Riffle for the pull request.
