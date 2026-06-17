Frequently Asked Questions
==========================

.. contents:: The following is a list of common questions asked about pyOpenMS.

How can I wrap a new method with pyOpenMS?
**************************************************

Add a ``.def(...)`` call for the method in the appropriate ``src/pyOpenMS/bindings/bind_<domain>.cpp``
file using nanobind syntax. See the `wrapping guide <../community/wrapping_workflows_new_classes.html#how-to-wrap-new-classes>`_ for details.


How can I wrap a new class with pyOpenMS?
*************************************************

Add a ``nb::class_<OpenMS::ClassName>(m, "ClassName", "docstring")`` block in the appropriate
``src/pyOpenMS/bindings/bind_<domain>.cpp`` file and use the `procedure outlined <../community/wrapping_workflows_new_classes.html#how-to-wrap-new-classes>`_. 


Can I use multiple output parameters?
*************************************

Python does not support passing primitive types (``int``, ``double``, etc.) by reference, therefore ``void calculate(double &)`` will not work.
