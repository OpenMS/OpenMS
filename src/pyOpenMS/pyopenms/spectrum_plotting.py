"""Convert MSSpectrum to spectrum_utils MsmsSpectrum.

Specified in:
Wout Bittremieux. “spectrum_utils: A Python package for mass spectrometry data processing and visualization.”
"""

# Code adopted from:
# Wout Bittremieux. “spectrum_utils: A Python package for mass spectrometry data processing and visualization.”
# Analytical Chemistry 92 (1) 659-661 (2020) doi:10.1021/acs.analchem.9b04884.


#                               Apache License
#                         Version 2.0, January 2004
#                      http://www.apache.org/licenses/
#
# TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION
#
# 1. Definitions.
#
#    "License" shall mean the terms and conditions for use, reproduction,
#    and distribution as defined by Sections 1 through 9 of this document.
#
#    "Licensor" shall mean the copyright owner or entity authorized by
#    the copyright owner that is granting the License.
#
#    "Legal Entity" shall mean the union of the acting entity and all
#    other entities that control, are controlled by, or are under common
#    control with that entity. For the purposes of this definition,
#    "control" means (i) the power, direct or indirect, to cause the
#    direction or management of such entity, whether by contract or
#    otherwise, or (ii) ownership of fifty percent (50%) or more of the
#    outstanding shares, or (iii) beneficial ownership of such entity.
#
#    "You" (or "Your") shall mean an individual or Legal Entity
#    exercising permissions granted by this License.
#
#    "Source" form shall mean the preferred form for making modifications,
#    including but not limited to software source code, documentation
#    source, and configuration files.
#
#    "Object" form shall mean any form resulting from mechanical
#    transformation or translation of a Source form, including but
#    not limited to compiled object code, generated documentation,
#    and conversions to other media types.
#
#    "Work" shall mean the work of authorship, whether in Source or
#    Object form, made available under the License, as indicated by a
#    copyright notice that is included in or attached to the work
#    (an example is provided in the Appendix below).
#
#    "Derivative Works" shall mean any work, whether in Source or Object
#    form, that is based on (or derived from) the Work and for which the
#    editorial revisions, annotations, elaborations, or other modifications
#    represent, as a whole, an original work of authorship. For the purposes
#    of this License, Derivative Works shall not include works that remain
#    separable from, or merely link (or bind by name) to the interfaces of,
#    the Work and Derivative Works thereof.
#
#    "Contribution" shall mean any work of authorship, including
#    the original version of the Work and any modifications or additions
#    to that Work or Derivative Works thereof, that is intentionally
#    submitted to Licensor for inclusion in the Work by the copyright owner
#    or by an individual or Legal Entity authorized to submit on behalf of
#    the copyright owner. For the purposes of this definition, "submitted"
#    means any form of electronic, verbal, or written communication sent
#    to the Licensor or its representatives, including but not limited to
#    communication on electronic mailing lists, source code control systems,
#    and issue tracking systems that are managed by, or on behalf of, the
#    Licensor for the purpose of discussing and improving the Work, but
#    excluding communication that is conspicuously marked or otherwise
#    designated in writing by the copyright owner as "Not a Contribution."
#
#    "Contributor" shall mean Licensor and any individual or Legal Entity
#    on behalf of whom a Contribution has been received by Licensor and
#    subsequently incorporated within the Work.
#
# 2. Grant of Copyright License. Subject to the terms and conditions of
#    this License, each Contributor hereby grants to You a perpetual,
#    worldwide, non-exclusive, no-charge, royalty-free, irrevocable
#    copyright license to reproduce, prepare Derivative Works of,
#    publicly display, publicly perform, sublicense, and distribute the
#    Work and such Derivative Works in Source or Object form.
#
# 3. Grant of Patent License. Subject to the terms and conditions of
#    this License, each Contributor hereby grants to You a perpetual,
#    worldwide, non-exclusive, no-charge, royalty-free, irrevocable
#    (except as stated in this section) patent license to make, have made,
#    use, offer to sell, sell, import, and otherwise transfer the Work,
#    where such license applies only to those patent claims licensable
#    by such Contributor that are necessarily infringed by their
#    Contribution(s) alone or by combination of their Contribution(s)
#    with the Work to which such Contribution(s) was submitted. If You
#    institute patent litigation against any entity (including a
#    cross-claim or counterclaim in a lawsuit) alleging that the Work
#    or a Contribution incorporated within the Work constitutes direct
#    or contributory patent infringement, then any patent licenses
#    granted to You under this License for that Work shall terminate
#    as of the date such litigation is filed.
#
# 4. Redistribution. You may reproduce and distribute copies of the
#    Work or Derivative Works thereof in any medium, with or without
#    modifications, and in Source or Object form, provided that You
#    meet the following conditions:
#
#    (a) You must give any other recipients of the Work or
#        Derivative Works a copy of this License; and
#
#    (b) You must cause any modified files to carry prominent notices
#        stating that You changed the files; and
#
#    (c) You must retain, in the Source form of any Derivative Works
#        that You distribute, all copyright, patent, trademark, and
#        attribution notices from the Source form of the Work,
#        excluding those notices that do not pertain to any part of
#        the Derivative Works; and
#
#    (d) If the Work includes a "NOTICE" text file as part of its
#        distribution, then any Derivative Works that You distribute must
#        include a readable copy of the attribution notices contained
#        within such NOTICE file, excluding those notices that do not
#        pertain to any part of the Derivative Works, in at least one
#        of the following places: within a NOTICE text file distributed
#        as part of the Derivative Works; within the Source form or
#        documentation, if provided along with the Derivative Works; or,
#        within a display generated by the Derivative Works, if and
#        wherever such third-party notices normally appear. The contents
#        of the NOTICE file are for informational purposes only and
#        do not modify the License. You may add Your own attribution
#        notices within Derivative Works that You distribute, alongside
#        or as an addendum to the NOTICE text from the Work, provided
#        that such additional attribution notices cannot be construed
#        as modifying the License.
#
#    You may add Your own copyright statement to Your modifications and
#    may provide additional or different license terms and conditions
#    for use, reproduction, or distribution of Your modifications, or
#    for any such Derivative Works as a whole, provided Your use,
#    reproduction, and distribution of the Work otherwise complies with
#    the conditions stated in this License.
#
# 5. Submission of Contributions. Unless You explicitly state otherwise,
#    any Contribution intentionally submitted for inclusion in the Work
#    by You to the Licensor shall be under the terms and conditions of
#    this License, without any additional terms or conditions.
#    Notwithstanding the above, nothing herein shall supersede or modify
#    the terms of any separate license agreement you may have executed
#    with Licensor regarding such Contributions.
#
# 6. Trademarks. This License does not grant permission to use the trade
#    names, trademarks, service marks, or product names of the Licensor,
#    except as required for reasonable and customary use in describing the
#    origin of the Work and reproducing the content of the NOTICE file.
#
# 7. Disclaimer of Warranty. Unless required by applicable law or
#    agreed to in writing, Licensor provides the Work (and each
#    Contributor provides its Contributions) on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
#    implied, including, without limitation, any warranties or conditions
#    of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A
#    PARTICULAR PURPOSE. You are solely responsible for determining the
#    appropriateness of using or redistributing the Work and assume any
#    risks associated with Your exercise of permissions under this License.
#
# 8. Limitation of Liability. In no event and under no legal theory,
#    whether in tort (including negligence), contract, or otherwise,
#    unless required by applicable law (such as deliberate and grossly
#    negligent acts) or agreed to in writing, shall any Contributor be
#    liable to You for damages, including any direct, indirect, special,
#    incidental, or consequential damages of any character arising as a
#    result of this License or out of the use or inability to use the
#    Work (including but not limited to damages for loss of goodwill,
#    work stoppage, computer failure or malfunction, or any and all
#    other commercial damages or losses), even if such Contributor
#    has been advised of the possibility of such damages.
#
# 9. Accepting Warranty or Additional Liability. While redistributing
#    the Work or Derivative Works thereof, You may choose to offer,
#    and charge a fee for, acceptance of support, warranty, indemnity,
#    or other liability obligations and/or rights consistent with this
#    License. However, in accepting such obligations, You may act only
#    on Your own behalf and on Your sole responsibility, not on behalf
#    of any other Contributor, and only if You agree to indemnify,
#    defend, and hold each Contributor harmless for any liability
#    incurred by, or claims asserted against, such Contributor by reason
#    of your accepting any such warranty or additional liability.
#
# END OF TERMS AND CONDITIONS

from typing import Optional, Callable
from pyopenms import MSSpectrum, StringDataArray


def parse_annotation(annotation: str):
    """ Parse an OpenMS fragment annotation string.

    :param annotation: The annotation string to be parsed.

    :returns: A dictionary containing the parsed annotation.
    """
    import re
    ion_type, neutral_loss, adduct, charge = None, None, None, 0
    regex = '([abcxyz]+\d+)-?(\w*)([\+\-]?)'
    match = re.search(regex, annotation)
    adduct_regex = '(\[\S+\])\s?(\S+)?'
    adduct_match = re.search(adduct_regex, annotation)

    if match is not None:
        ion_type = match.group(1)
        neutral_loss = match.group(2) if match.group(2) != '' else None
        charge = match.group(3).count('+') - match.group(3).count('-')
    if annotation.find(':') != -1:  # compound
        ion_type = '_' + annotation
    if adduct_match is not None:
        ion_type = 'f' if ion_type is None else ion_type
        adduct = adduct_match.group(1)
        neutral_loss = adduct_match.group(2) if adduct_match.group(2) != '' else None
    if ion_type is None:
        ion_type = '?' + annotation
    if ion_type != '?' and charge == 0:
        charge = 1

    analyte_number = None
    isotope = 0
    mz_delta = None
    return {
        'analyte_number': analyte_number,
        'ion_type': ion_type,
        'neutral_loss': neutral_loss,
        'isotope': isotope,
        'charge': charge,
        'adduct': adduct,
        'mz_delta': mz_delta
    }


def spectrum_to_spectrum_utils(
        spectrum: "MSSpectrum", annot_func: Optional[Callable] = lambda x: parse_annotation(x)
):
    """ Convert an OpenMS MSSpectrum to a spectrum_utils MsmsSpectrum.

    :param spectrum: The OpenMS MSSpectrum to be plotted.
    Reads annotations from the first StringDataArray if it has the same length as the number of peaks.

    :param annot_func: Fragment annotation function parsing the annotation string. Is required to return a dictionary
    with the following keys: 'analyte_number', 'ion_type', 'neutral_loss', 'isotope', 'charge', 'adduct', 'mz_delta'.

    :returns: The spectrum_utils MsmsSpectrum.

    """
    from spectrum_utils.spectrum import MsmsSpectrum
    from spectrum_utils.fragment_annotation import FragmentAnnotation, PeakInterpretation
    if len(spectrum.getPrecursors()) > 0:
        precursor_mz = spectrum.getPrecursors()[0].getMZ()
        precursor_charge = spectrum.getPrecursors()[0].getCharge()
    else:
        precursor_mz = -1
        precursor_charge = -1

    utils_spectrum = MsmsSpectrum(
        identifier=spectrum.getNativeID(),
        precursor_mz=precursor_mz,
        precursor_charge=precursor_charge,
        mz=list(spectrum.get_peaks()[0]),
        intensity=list(spectrum.get_peaks()[1])
    )
    peak_interpretations = []
    for string_data in spectrum.getStringDataArrays()[0]:
        annotation = str(string_data)[1:].replace('"', '').replace("'", '')
        pi = PeakInterpretation()
        if annotation != '':
            annotation_data = annot_func(annotation)
            f = FragmentAnnotation(**annotation_data)
            pi.fragment_annotations.append(f)
        peak_interpretations.append(pi)
    utils_spectrum._annotation = peak_interpretations
    return utils_spectrum


def peptide_hit_to_spectrum_utils(
        peptide_hit: "PeptideHit", annot_func: Optional[Callable] = lambda x: parse_annotation(x)
):
    """ Convert an OpenMS peptideHit to a spectrum_utils MsmsSpectrum.

    :param peptide_hit: The OpenMS peptideHit to be plotted.
    Reads annotations from the first StringDataArray if it has the same length as the number of peaks.

    :param annot_func: Fragment annotation function parsing the annotation string. Is required to return a dictionary
    with the following keys: 'analyte_number', 'ion_type', 'neutral_loss', 'isotope', 'charge', 'adduct', 'mz_delta'.

    :returns: The spectrum_utils MsmsSpectrum.

    """
    intensities, mass_to_charges, annotations = [], [], []

    for peak in peptide_hit.getPeakAnnotations():
        intensities.append(peak.intensity)
        mass_to_charges.append(peak.mz)
        annotations.append(str(peak.annotation))

    spectrum = MSSpectrum()
    spectrum.set_peaks([mass_to_charges, intensities])

    # Sort the peaks according to ascending mass-to-charge ratio
    spectrum.sortByPosition()

    sda = StringDataArray()
    sda.setName("Peak annotation")
    for element in annotations:
        sda.push_back(element)

    spectrum.setStringDataArrays([sda])
    return spectrum_to_spectrum_utils(spectrum, annot_func)


def show_annotations(annotation: "FragmentAnnotation", ion_types: str='abcxyz', ions: Optional[bool]=True,
                     charges: Optional[bool]=False, losses: Optional[bool]=False, adducts: Optional[bool]=False):
    """ toggle which annotations to show.

    :param annotation: The annotation string to be parsed.

    :param ion_types: The ion types to show. Default: 'abcxyz'.

    :param ions: Whether to show the ion type. Default: True.

    :param charges: Whether to show the charge. Default: False.

    :param losses: Whether to show the neutral loss. Default: False.

    :param adducts: Whether to show the adduct. Default: False.

    :returns: The annotation string.

    """
    output = ''
    if ions:
        output += annotation.ion_type if annotation.ion_type[0] in ion_types else annotation.ion_type[1:]
    if charges:
        output += f'{annotation.charge:+d}' if annotation.charge is not None else ''
    if losses:
        output += annotation.neutral_loss if annotation.neutral_loss is not None else ''
    if adducts:
        output += annotation.adduct if annotation.adduct is not None else ''
    return output

