Calibration
===========

OpenMS offer two calibration methods: an internal and an external calibration. Both can handle peak data as well as
profile data. To calibrate profile data, a peak picking step is necessary, the important parameters can be set via the
ini-file. If you have already picked data, don't forget the `-peak_data` flag.

The external calibration (**TOFCalibration**) is used to convert flight times into m/z values with the help of external
calibrant spectra containing e.g. a polymer like polylysine. For the calibrant spectra, the calibration constants the
machine uses need to be known as well as the expected masses. Then a quadratic function is fitted to convert the flight
times into m/z-values.

The internal calibration (**InternalCalibration**) uses reference masses in the spectra to correct the m/z-values
using a linear function.

In a typical setting one would first pick the TOF-data, then perform the TOFCalibration and then the InternalCalibration:

```bash
PeakPickerWavelet -in raw_tof.mzML -out picked_tof.mzML -ini pp.ini
TOFCalibration -in picked_tof.mzML -out picked.mzML -ext_calibrants ext_cal.mzML
               -ref_masses ext_cal_masses
               -tof_const tof_conv_consts -peak_data
InternalCalibration -in picked.mzML -out picked_calibrated.mzML
                    -ref_masses internal_calibrant_masses -peak_data
```
