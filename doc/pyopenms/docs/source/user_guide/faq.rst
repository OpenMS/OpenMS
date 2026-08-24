Frequently Asked Questions
==========================

.. contents:: The following is a list of common questions asked about pyOpenMS.

I modified a spectrum from my experiment, but the change is gone. Why?
**********************************************************************

pyOpenMS containers use value semantics: indexing, iteration, and getters
return independent copies, so ``exp[0].setRT(5.0)`` edits only the copy.
To change data inside a container, read it, edit it, and put it back:

.. code-block:: python

    spec = exp[0]        # read it (a copy)
    spec.setRT(5.0)      # edit it
    exp[0] = spec        # put it back

For zero-copy access — reading large data without duplicating it, or
editing in place — use the *view* accessors, whose names end in ``_view``
(``exp.spectrum_view(i)``, ``exp.iter_spectrum_views()``); edits through a
view land directly, but the view is only valid until the container is
resized or sorted. The rule of thumb: anything called ``get_*`` returns a
copy you own; ``_view``/``_views``/``_struct`` methods alias. See the
`MS data <ms_data.html>`_ chapter for details.


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
