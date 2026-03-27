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
# Example of displaying a Streptococcal protein and showing peptide coverage
import pymol
import pyopenms

print("Starting pymol")
pymol.finish_launching()
pepxml_file = "data/4D8B.pepXML"


def get_peptides_protein_seq():
    # Read in PepXML
    protein_ids = []
    peptide_ids = []
    pyopenms.PepXMLFile().load(pepxml_file, protein_ids, peptide_ids)
    peptides = [pid.getHits()[0].getSequence().toString()
                for pid in peptide_ids]

    # Sequence could be from FASTA file (or just provided here)
    sequence = "".join(
        """
    MNKKKLGIRLLSLLALGGFVLANPVFADQNFARNEKEAKDSAITFIQKSAAIKAGARSAE
    DIKLDKVNLGGELSGSNMYVYNISTGGFVIVSGDKRSPEILGYSTSGSFDANGKENIASF
    MESYVEQIKENKKLDTTYAGTAEIKQPVVKSLLNSKGIHYNQGNPYNLLTPVIEKVKPGE
    QSFVGQHAATGCVATATAQIMKYHNYPNKGLKDYTYTLSSNNPYFNHPKNLFAAISTRQY
    NWNNILPTYSGRESNVQKMAISELMADVGISVDMDYGPSSGSAGSSRVQRALKENFGYNQ
    SVHQINRSDFSKQDWEAQIDKELSQNQPVYYQGVGKVGGHAFVIDGADGRNFYHVNWGWG
    GVSDGFFRLDALNPSALGTGGGAGGFNGYQSAVVGIKP
    """.split()
    )
    return peptides, sequence


print("Plotting Protein 4D8B")
pymol.cmd.window("show")
pymol.cmd.fetch("4D8B")
pymol.cmd.hide("everything")
pymol.cmd.show("cartoon")

peptides, protein_sequence = get_peptides_protein_seq()
# select which residues to color (set of peptides, on protein sequence)
res_str = ""
for peptide in peptides:
    found_at = protein_sequence.find(peptide)
    res_str += " res %s-%s" % (found_at, found_at + len(peptide))

print("Superimposing MS-identified peptides .. please wait")
pymol.cmd.color("bluewhite")
pymol.cmd.bg_color("white")
pymol.cmd.color("wheat", res_str)

# Do raytracing
pymol.cmd.set("orthoscopic")
pymol.cmd.set("ray_trace_mode", 1)
pymol.cmd.ray(1920, 1440)

print("Printing output image")
pymol.cmd.png("4D8B_coverage_white_5_ray_ortho_mode1.png", dpi=600)
