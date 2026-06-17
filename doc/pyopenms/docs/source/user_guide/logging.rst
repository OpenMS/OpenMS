Logging 
=======

Since pyopenms is mainly a C(++) extension and it is not easy to redirect the complete C++ logging mechanism to python,
it does not (yet) support full control over its logging with the built-in `python-based logging module <https://docs.python.org/3/library/logging.html>`_.

Instead, you currently need to control logging on the python side AND on the C extension side. Since our
python-only part is rather small and does not use much logging, we will skip this part for now.

Logging in the wrapped C extension part can be controlled by the following:

.. code-block:: python

   import pyopenms as pyoms

   pyoms.LogConfigHandler().setLogLevel("ERROR")

Unfortunately, there are still some limitations: `OpenMS issue 6827 <https://github.com/OpenMS/OpenMS/issues/6827>`_
Controlling logging inside the C extension is only supported globally per process. That means, once you change the log level,
it will stay changed for all subsequent calls, until you change it back or the program terminates.
Behaviour in multithreaded regions is untested.
    