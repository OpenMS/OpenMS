#ifndef __PYTHON_MS_DATA_CONSUMER_HPP__
#define __PYTHON_MS_DATA_CONSUMER_HPP__

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>

// see ../pxds/PythonMSDataConsumer.pxd for Cython def
class PythonMSDataConsumer :
  virtual public OpenMS::Interfaces::IMSDataConsumer
{

    typedef OpenMS::PeakMap::SpectrumType SpectrumType;
    typedef OpenMS::PeakMap::ChromatogramType ChromatogramType;

    // typedefs for function ptr (helper fxn to convert C++ type to Python object)
    // see ../addons/MzMLFile.pyx
    typedef PyObject* (*SpectrumToPythonWrapper) (const SpectrumType &);
    typedef PyObject* (*ChromatogramToPythonWrapper) (const ChromatogramType &);
    typedef PyObject* (*ExperimentalSettingsToPythonWrapper) (const OpenMS::ExperimentalSettings &);

    private:

        PyObject *py_consumer_;

        SpectrumToPythonWrapper wrap_spectrum_;
        ChromatogramToPythonWrapper wrap_chromatogram_;
        ExperimentalSettingsToPythonWrapper wrap_experimental_settings_;

    public:

        /// Constructor
        PythonMSDataConsumer(PyObject *py_consumer,
                             SpectrumToPythonWrapper wrap_spectrum,
                             ChromatogramToPythonWrapper wrap_chromatogram,
                             ExperimentalSettingsToPythonWrapper wrap_experimental_settings) :
          py_consumer_(py_consumer),
          wrap_spectrum_(wrap_spectrum),
          wrap_chromatogram_(wrap_chromatogram),
          wrap_experimental_settings_(wrap_experimental_settings)
        {
           Py_INCREF(py_consumer_);
        };

        /// Destructor
        ~PythonMSDataConsumer()
        {
           Py_DECREF(py_consumer_);
        };


        /// Consume spectrum (call Python method "consumeSpectrum" of the py_consumer_ object from C++)
        virtual void consumeSpectrum(SpectrumType & spec)
        {
            PyObject * py_spec = wrap_spectrum_(spec);
            PyObject * method_name = PyUnicode_FromString("consumeSpectrum");
            PyObject * r = PyObject_CallMethodObjArgs(py_consumer_, method_name, py_spec, NULL);
            Py_DECREF(py_spec);
            Py_DECREF(method_name);
            // NULL indicates python exception:
            if (r == NULL)
                throw "exception"; // not sense needed, as cython evaluates python strack trace
            Py_DECREF(r);
        };

        /// Consume chromatogram (call Python method "consumeChromatogram" of the py_consumer_ object from C++)
        virtual void consumeChromatogram(ChromatogramType & chrom)
        {
            PyObject * py_chrom = wrap_chromatogram_(chrom);
            PyObject * method_name = PyUnicode_FromString("consumeChromatogram");
            PyObject * r = PyObject_CallMethodObjArgs(py_consumer_, method_name, py_chrom, NULL);
            Py_DECREF(py_chrom);
            Py_DECREF(method_name);
            // NULL indicates python exception:
            if (r == NULL)
                throw "exception"; // not sense needed, as cython evaluates python strack trace
            Py_DECREF(r);
        };

        virtual void setExpectedSize(OpenMS::Size expectedSpectra,
                                     OpenMS::Size expectedChromatograms)
        {
            PyObject * expected_spectra = PyInt_FromSize_t(expectedSpectra);
            PyObject * expected_chromatograms = PyInt_FromSize_t(expectedChromatograms);
            PyObject * method_name = PyUnicode_FromString("setExpectedSize");
            PyObject * r = PyObject_CallMethodObjArgs(py_consumer_, method_name, expected_spectra,
                                                      expected_chromatograms, NULL);
            Py_DECREF(expected_spectra);
            Py_DECREF(expected_chromatograms);
            Py_DECREF(method_name);
            // NULL indicates python exception:
            if (r == NULL)
                throw "exception"; // not sense needed, as cython evaluates python strack trace
            Py_DECREF(r);
        };

        virtual void setExperimentalSettings(const OpenMS::ExperimentalSettings & exp_settings)
        {
            PyObject * py_exp_settings = wrap_experimental_settings_(exp_settings);
            PyObject * method_name = PyUnicode_FromString("setExperimentalSettings");
            PyObject * r = PyObject_CallMethodObjArgs(py_consumer_,
                                                      method_name, py_exp_settings, NULL);

            Py_DECREF(py_exp_settings);
            Py_DECREF(method_name);
            // NULL indicates python exception:
            if (r == NULL)
                throw "exception"; // not sense needed, as cython evaluates python strack trace
            Py_DECREF(r);
        };
};

#endif
