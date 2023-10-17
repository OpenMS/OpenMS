"""Plot MSSpectrum & MSChromatogram data

The MSSpectrum plotting function are adapted from:
Wout Bittremieux. “spectrum_utils: A Python package for mass spectrometry data processing and visualization.”
Plot a single spectrum with plot_spectrum or two with mirror_plot_spectrum, using matplotlib.
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

import math
from typing import Dict, Optional, Tuple, Union, List, Set
import itertools

colors = {'a': '#388E3C', 'b': '#1976D2', 'c': '#00796B',
          'x': '#7B1FA2', 'y': '#D32F2F', 'z': '#F57C00',
          'p': '#512DA8', '?': '#212121', 'f': '#212121', None: '#212121'}
zorders = {'a': 3, 'b': 4, 'c': 3, 'x': 3, 'y': 4, 'z': 3,
           'p': 3, '?': 2, 'f': 5, None: 1}



# TODO switch to forward declarations via from __future__ import annotations
# when py 3.7 is minimum. Then you can use types instead of strings
def plot_chromatogram(c: "MSChromatogram"):
    """Plot chromatogram peaks.

    :param c: The chromatogram to be plotted.
    :type c: MSChromatogram
    """

    import matplotlib.pyplot as plt
    x, y = c.get_peaks()
    plt.plot(x, y)
    plt.xlabel("Retention time")
    plt.ylabel("Intensity")
    plt.show()


def _annotate_ion(mz: float, intensity: float, annotation: Optional[str],
                  color_ions: bool, annotate_ions: bool, matched: Optional[bool],
                  annotation_kws: Dict[str, object], ax) -> Tuple[str, int]:
    """Annotate a specific fragment peak.

    :param mz: The peak's m/z value (position of the annotation on the x axis).
    :type mz: float

    :param intensity: The peak's intensity (position of the annotation on the y axis).
    :type intensity: float

    :param annotation: The annotation that will be plotted.
    :type annotation: str, optional

    :param color_ions: Flag whether to color the peak annotation or not.
    :type color_ions: bool

    :param annotate_ions: Flag whether to annotation the peak or not.
    :type annotate_ions: bool

    :param annotation_kws:  Keyword arguments for `ax.text` to customize peak annotations.
    :type annotation_kws: Dict[str, object]

    :param ax: Axes instance on which to plot the annotation.
    :type ax: plt.Axes


    :return: A tuple of the annotation's color as a hex string and the annotation's zorder.
    :rtype: Tuple[str, int]
    """

    # No annotation -> Just return peak styling information.
    if annotation is None:
        return colors.get(None), zorders.get(None)
    # Else: Add the textual annotation.
    ion_type = annotation[0]
    if ion_type == '[': # precursor ion
        ion_type = 'p'
    if ion_type not in colors and color_ions:
        raise ValueError('Ion type not supported')

    color = (colors[ion_type] if color_ions else
             colors[None])
    zorder = (1 if ion_type not in zorders else zorders[ion_type])

    if matched is not None and not matched:
        color = '#aaaaaa'

    if annotate_ions:
        annotation_pos = intensity
        if annotation_pos > 0:
            annotation_pos += 0.02
        kws = annotation_kws.copy()
        del kws['zorder']
        ax.text(mz, annotation_pos, str(annotation), color=color,
                zorder=zorder, **kws)

    return color, zorder


def plot_spectrum(spectrum: "MSSpectrum", color_ions: bool = True,
                  annotate_ions: bool = True, matched_peaks: Optional[Set] = None, annot_kws: Optional[Dict] = None,
                  mirror_intensity: bool = False, grid: Union[bool, str] = True, ax=None):
    """Plot an MS/MS spectrum.

    :param spectrum: The spectrum to be plotted.
        Reads annotations from the first StringDataArray if it has the same length as the number of peaks.
    :type spectrum: MSSpectrum

    :param color_ions: Flag indicating whether to color annotated fragment ions. The default is True.
    :type color_ions: bool, optional

    :param annotate_ions: Flag indicating whether to annotate fragment ions. The default is True.
    :type annotate_ions: bool, optional

    :param matched_peaks: Indices of matched peaks in a spectrum alignment.
    :type matched_peaks: Optional[Set], optional

    :param annot_kws: Keyword arguments for `ax.text` to customize peak annotations.
    :type annot_kws: Optional[Dict], optional

    :param mirror_intensity: Flag indicating whether to flip the intensity axis or not.
    :type mirror_intensity : bool, optional

    :param grid: Draw grid lines or not. Either a boolean to enable/disable both major
        and minor grid lines or 'major'/'minor' to enable major or minor grid lines respectively.
    :type grid: Union[bool, str], optional

    :param ax: Axes instance on which to plot the spectrum. If None the current Axes instance is used.
    :type ax : Optional[plt.Axes], optional


    :return The matplotlib Axes instance on which the spectrum is plotted.
    :rtype: plt.Axes
    """

    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker

    if ax is None:
        ax = plt.gca()

    mz, intensity = spectrum.get_peaks()
    min_mz = max(0, math.floor(mz[0] / 100 - 1) * 100)
    max_mz = math.ceil(mz[-1] / 100 + 1) * 100
    ax.set_xlim(min_mz, max_mz)
    ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1.))
    y_max = 1.15 if annotate_ions else 1.05
    ax.set_ylim(*(0, y_max) if not mirror_intensity else (-y_max, 0))

    max_intensity = intensity.max()
    if max_intensity == 0: max_intensity = 1
    if len(spectrum.getStringDataArrays()) > 0 and len(list(spectrum.getStringDataArrays()[0])) == len(mz):
        annotations = [ion.decode() for ion in spectrum.getStringDataArrays()[0]]
    else:
        annotations = itertools.repeat(None)
    annotation_kws = {
        'horizontalalignment': 'left' if not mirror_intensity else 'right',
        'verticalalignment': 'center', 'rotation': 90,
        'rotation_mode': 'anchor', 'zorder': 5}
    if annot_kws is not None:
        annotation_kws.update(annot_kws)
    for i_peak, (peak_mz, peak_intensity, peak_annotation) in enumerate(zip(mz, intensity, annotations)):
        peak_intensity = peak_intensity / max_intensity
        if mirror_intensity:
            peak_intensity *= -1

        matched = matched_peaks is not None and i_peak in matched_peaks

        color, zorder = _annotate_ion(
            peak_mz, peak_intensity, peak_annotation, color_ions, annotate_ions,
            matched, annotation_kws, ax)

        ax.plot([peak_mz, peak_mz], [0, peak_intensity], color=color, zorder=zorder)

    ax.xaxis.set_minor_locator(mticker.AutoLocator())
    ax.yaxis.set_minor_locator(mticker.AutoLocator())
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    if grid in (True, 'both', 'major'):
        ax.grid(visible=True, which='major', color='#9E9E9E', linewidth=0.2)
    if grid in (True, 'both', 'minor'):
        ax.grid(visible=True, which='minor', color='#9E9E9E', linewidth=0.2)
    ax.set_axisbelow(True)

    ax.tick_params(axis='both', which='both', labelsize='small')
    y_ticks = ax.get_yticks()
    ax.set_yticks(y_ticks[y_ticks <= 1.])

    ax.set_xlabel('m/z', style='italic')
    ax.set_ylabel('Intensity')

    return ax


def mirror_plot_spectrum(spec_top: "MSSpectrum", spec_bottom: "MSSpectrum", alignment: Optional[List] = None,
                         spectrum_top_kws: Optional[Dict] = None, spectrum_bottom_kws: Optional[Dict] = None,
                         ax=None):
    """Mirror plot two MS/MS spectra.

    :param spec_top: The spectrum to be plotted on the top.
        Reads annotations from the first StringDataArray if it has the same length as the number of peaks.
    :type spec_top: MSSpectrum

    :param spec_bottom: The spectrum to be plotted on the bottom.
        Reads annotations from the first StringDataArray if it has the same length as the number of peaks.
    :type spec_bottom: MSSpectrum

    :param alignment: List of aligned peak pairs.
    :type alignment: Optional[List], optional

    :param spectrum_top_kws: Keyword arguments for `Plotting.plot_spectrum` of top spectrum.
    :type spectrum_top_kws: Optional[Dict], optional

    :param spectrum_bottom_kws: Keyword arguments for `Plotting.plot_spectrum` of bottom spectrum.
    :type spectrum_bottom_kws: Optional[Dict], optional

    :param ax: Axes instance on which to plot the spectrum. If None the current Axes instance is used.
    :type ax: Optional[plt.Axes], optional


    :return: The matplotlib Axes instance on which the spectra are plotted.
    :rtype: plt.Axes
    """

    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    if ax is None:
        ax = plt.gca()

    if spectrum_top_kws is None:
        spectrum_top_kws = {}

    if spectrum_bottom_kws is None:
        spectrum_bottom_kws = {}

    if alignment is not None:
        matched_peaks_bottom, matched_peaks_top = set(zip(*alignment))
    else:
        matched_peaks_bottom, matched_peaks_top = None, None

    # Top spectrum.
    plot_spectrum(spec_top, mirror_intensity=False, ax=ax, matched_peaks=matched_peaks_top, **spectrum_top_kws)
    y_max = ax.get_ylim()[1]
    # Mirrored bottom spectrum.
    plot_spectrum(spec_bottom, mirror_intensity=True, ax=ax, matched_peaks=matched_peaks_bottom, **spectrum_bottom_kws)
    y_min = ax.get_ylim()[0]
    ax.set_ylim(y_min, y_max)

    ax.axhline(0, color='#9E9E9E', zorder=10)

    # Update axes so that both spectra fit.
    spec_top_mz, sp_top_intensity = spec_top.get_peaks()
    spec_bottom_mz, sp_bottom_intensity = spec_bottom.get_peaks()
    min_mz = max([0, math.floor(spec_top_mz[0] / 100 - 1) * 100,
                  math.floor(spec_bottom_mz[0] / 100 - 1) * 100])
    max_mz = max([math.ceil(spec_top_mz[-1] / 100 + 1) * 100,
                  math.ceil(spec_bottom_mz[-1] / 100 + 1) * 100])
    ax.set_xlim(min_mz, max_mz)
    ax.yaxis.set_major_locator(mticker.AutoLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(
        lambda x, pos: f'{abs(x):.0%}'))

    return ax
