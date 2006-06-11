#!/bin/sh

## Simple Tandem MS Identification Pipeline

# Convert raw data to mzData format.
FileConverter -in $1 -out id.mzData

# Extract MS/MS spectra only.
FileFilter -in id.mzData -out raw.mzData -level 2

# Noise filtering, baseline subtraction
NoiseFilter -in raw.mzData -out rawf.mzData -ini ID.ini

# Construct stick spectra
PeakPicker -in rawf.mzData -out peaks.mzData -ini ID.ini

# Identify spectra.
MascotAdapter -in peaks.mzData -out id.xml -ini ID.ini

# Filter out reliable identifications
IDFilter -in id.xml -out result.xml -ini ID.ini
