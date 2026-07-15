---
orphan: true
---
Subtracting a Baseline from a Spectrum
=====================================

First, load the spectrum to be analyzed in TOPPView. To use the described tools, open the tutorial data via the
File-menu (**File** > **Open example file**, then select `peakpicker_tutorial_1.mzML`). The BaselineFilter can be called via
the Tools-menu (**Tools** > **Apply TOPP tool**), then select **BaselineFilter** as TOPPtool (red rectangle). You can choose,
between different types of filters (green rectangle), the one mainly used is TopHat. The other important parameter is
the length of the structuring element (blue rectangle). The default value is `3` Thomson. Press **Ok** to start the baseline
subtraction.

![TOPPView Tools Baseline](/_images/tutorials/topp/TOPPView_tools_baseline.png)

The following image shows:
- A part of the spectrum after baseline filtering as a green line.
- The original raw data as a blue line.

![TOPPView Tools Baseline Filtered](/_images/tutorials/topp/TOPPView_tools_baseline_filtered.png)
