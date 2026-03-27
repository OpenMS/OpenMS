Support
=======

Feature Requests
****************

pyOpenMS is an evolving project. We are happy to learn about missing features you would like to
see in pyOpenMS.

Feel free to open an issue on `GitHub issue <https://github.com/OpenMS/OpenMS/issues>`_
 to request a novel feature.

Troubleshooting
***************

If you encounter a problem running pyOpenMS feel free to contact
us by opening a `GitHub issue <https://github.com/OpenMS/OpenMS/issues>`_
or `via our Homepage <https://openms.de/communication/>`_ today.

Please provide information on what pyOpenMS version you have installed.
You can check if importing pyOpenMS works and print your pyOpenMS version with:

.. code-block:: python
    :linenos:

    import pyopenms as oms

    print("Version: " + oms.VersionInfo.getVersion())
    print("OpenMP: " + str(oms.OpenMSBuildInfo.isOpenMPEnabled()))
    print("Build type: " + oms.OpenMSBuildInfo.getBuildType())
    print("Architecture: " + oms.OpenMSOSInfo.getBinaryArchitecture())

