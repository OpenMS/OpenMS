Examples
=============

This folder contains example scripts using pyOpenMS:

* `convertToMGF.py`
    * A simple script to convert any spectra file readable by OpenMS (for
      example `mzXML` or `mzML`)  to the `MGF` format. Python is only used to
      read in the spectra while writing of MGF is done natively in Python. MGF
      is used an example of any file format that is not currenly available
      through the OpenMS C++ library, showing how such a sample implementation
      could look like.
- `peakpicker_scipyFFT.py`
    * A peakpicker using fast fourier transform (FFT) provided in `scipy` to
      smooth the provided spectrum and then applies the OpenMS
      `PeakPickerHiRes` on the result. The input spectrum and the picked peaks
      are then plotted using `matplotlib`, completing the integration of Python
      and OpenMS on multiple levels.
- `pymol_example.py`
    * An example of mapping identified peptide sequences from an MS/MS
      experiment onto a 3D crystal strucutre of a protein using Python with
      `pymol`.
