Frequently Asked Questions
==========================

.. contents:: The following is a list of common questions asked about pyOpenMS.

How can I wrap a new method with pyOpenMS?
**************************************************

Add an entry to ``src/pyOpenMS/pxds/CLASS_NAME.pxd`` with the signature of your new method(s).


How can I wrap a new class with pyOpenMS?
*************************************************

Create a new file ``src/pyOpenMS/pxds/CLASS_NAME.pxd`` and use the `procedure outlined <../community/wrapping_workflows_new_classes.html#how-to-wrap-new-classes>`_. 


Can I use multiple output parameters?
*************************************

Python does not support passing primitive types (``int``, ``double``, etc.) by reference, therefore ``void calculate(double &)`` will not work.
