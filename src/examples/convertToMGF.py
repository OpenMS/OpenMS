"""
--------------------------------------------------------------------------
                  OpenMS -- Open-Source Mass Spectrometry
--------------------------------------------------------------------------
Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
ETH Zurich, and Freie Universitaet Berlin 2002-2013.

This software is released under a three-clause BSD license:
 * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
 * Neither the name of any author or any participating institution
   may be used to endorse or promote products derived from this software
   without specific prior written permission.
For a full list of authors, refer to the file AUTHORS.
--------------------------------------------------------------------------
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
###########################################################################
## Example script to convert any file to MGF
###########################################################################

# Read input
import pyopenms, sys
if len(sys.argv) <= 2:
    print ("Usage: convertToMGF.py inputfile outputfile")
    sys.exit()

msdata = pyopenms.MSExperiment()
pyopenms.FileHandler().loadExperiment(sys.argv[1], msdata)
outfile = open(sys.argv[2], "w")

# Create header
outfile.write("COM=Testfile\n")
outfile.write("ITOL=1\n")
outfile.write("ITOLU=Da\n")
outfile.write("CLE=Trypsin\n")
outfile.write("CHARGE=1,2,3\n")

# Iterate through all spectra, skip all MS1 spectra and then write mgf format
nr_ms2_spectra = 0
for spectrum in msdata:
    if spectrum.getMSLevel() == 1:
        continue
    nr_ms2_spectra += 1
    outfile.write("\nBEGIN IONS\n")
    outfile.write("TITLE=%s\n" % spectrum.getNativeID())
    outfile.write("RTINSECONDS=%s\n" % spectrum.getRT())
    try:
        outfile.write("PEPMASS=%s\n" % spectrum.getPrecursors()[0].getMZ())
        ch = spectrum.getPrecursors()[0].getCharge()
        if ch > 0:
            outfile.write("CHARGE=%s\n" % ch)
    except IndexError:
        outfile.write("PEPMASS=unknown\n")
    for peak in spectrum:
        outfile.write("%s %s\n" % (peak.getMZ(), peak.getIntensity() ))
    outfile.write("END IONS\n")

if nr_ms2_spectra == 0:
    print("Did not find any MS2 spectra in your input, thus the output file is empty!")

