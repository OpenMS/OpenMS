---
orphan: true
---
Smoothing Raw Data
==================

To smooth raw data, call one of the available NoiseFilters via the Tools-menu, (select **Tools** > **Apply TOPP tool**), then select **NoiseFilterSGolay** or **NoiseFilterGaussian** as TOPP tool (green rectangle). The parameters for the filter type can be adapted (blue rectangle). For the `Savitzky-Golay` filter, set the **frame_length** and the **polynomial_order** fitted.
For the Gaussian filter, the gaussian width and the ppm tolerance for a flexible gaussian width depending on the `m/z`
value can be adapted. Press **Ok** to run the selected `NoiseFilter`.

![](/_images/tutorials/topp/TOPPView_tools_noisefilter.png)

The following image shows a part of the spectrum after smoothing as red line with the un-smoothed data in green.

![](/_images/tutorials/topp/TOPPView_tools_noisefilter_filtered.png)
