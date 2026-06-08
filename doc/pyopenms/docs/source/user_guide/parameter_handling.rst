Parameter Handling 
==================

Parameter handling in OpenMS and pyOpenMS is usually implemented through inheritance
from ``DefaultParamHandler`` and allow access to parameters through the :py:class:`~.Param` object. This
means, the classes implement the methods ``getDefaults``, ``getParameters`` and ``setParameters``,
to access to the default parameters, the current parameters and to set new parameters, respectively.
The class :py:class:`~.TheoreticalSpectrumGenerator` is just one example of many which makes use of parameter handling via 
``DefaultParamHandler``.

The :py:class:`~.Param` object is the central data structure here. It can be manipulated through the :py:meth:`~.Param.setValue`
and :py:meth:`~.Param.getValue` methods. The :py:meth:`~.Param.exists` method can be used to check for existence of a key and should
always be used if a param value might be missing, since accessing a missing value via :py:meth:`~.Param.getValue` 
will result in a RuntimeError exception.

Using the :py:meth:`~.Param.getDescription` method, it is possible to get a descriptive help for each parameter value in an
interactive session without consulting the documentation.

.. code-block:: python
    :linenos:

    import pyopenms as oms

    p = oms.Param()
    p.setValue("param1", 4.0, "This is value 1")
    p.setValue("param2", 5.0, "This is value 2")
    p.setValue(
        "param3",
        [b"H:+:0.6", b"Na:+:0.2", b"K:+:0.2"],
        "This is value 3 (StringList)",
    )
    print(p[b"param1"])
    p[b"param1"] += 3  # add three to the parameter value
    print(p[b"param1"])
    print(p[b"param3"])


The parameters can also be accessed as 

.. code-block:: pycon

    >>> p.asDict()
    {'param2': 4.0, 'param1': 7.0}
    >>> p.values()
    [4.0, 7.0]
    >>> p.keys()
    ['param1', 'param2']
    >>> p.items()
    [('param1', 7.0), ('param2', 4.0)]
    >>> p.exists("param1")
    True


The param object can be copied and merged into another param object:
 
.. code-block:: python
    :linenos:

    # print the key and value pairs stored in a Param object
    def printParamKeyAndValues(p):
        if p.size():
            for i in p.keys():
                print("Key:", i, "Value:", p[i])
        else:
            print("no data available")


    new_p = p  # new deep copy of p

    # we will add 4 more keys to the new_p
    new_p.setValue("param2", 9.0, "This is value 9")
    new_p.setValue("example1", 6.0, "This is value 6")
    new_p.setValue("example2", 8.0, "This is value 8")
    new_p.setValue("example3", 10.0, "This is value 10")

    # names "example1", "example2" , "example3" keys will added to p, but "param2" will update the value
    p.merge(new_p)
    print(" print the key and values pairs stored in a Param object p ")
    printParamKeyAndValues(p)

In a param object, the keys can be removed by key name or prefix:

.. code-block:: python
    :linenos:

    # We now call the remove method with the key of the entry we want to delete ("example3")
    new_p.remove("example3")
    print("Key and value pairs after removing the entry with key: example3")
    printParamKeyAndValues(new_p)

    # We now want to delete all keys with prefix "exam"
    new_p.removeAll("exam")
    print(
        "Key and value pairs after removing all entries with keys starting with: exam"
    )
    printParamKeyAndValues(new_p)

    # we can compare Param objects for identical content
    if p == new_p:  # check p is equal to new_p
        new_p.clear()  # Example: delete all keys from new_p

    print("Keys and values after deleting all entries.")
    printParamKeyAndValues(new_p)  # All keys of new_p deleted

For the algorithms that inherit from ``DefaultParamHandler``, you can list all parameters along with their 
description by using, for instance, the following simple function.

.. code-block:: python
    :linenos:

    # print all parameters with description
    def printParams(p):
        if p.size():
            for i in p.keys():
                print(
                    "Param:", i, "Value:", p[i], "Description:", p.getDescription(i)
                )
        else:
            print("no data available")

    # print all parameters in GaussFilter class
    gf = oms.GaussFilter()
    printParams(gf.getParameters())

.. code-block:: output

    Param: b'gaussian_width' Value: 0.2 Description: Use a gaussian filter width which has approximately the same width as your mass peaks (FWHM in m/z).
    Param: b'ppm_tolerance' Value: 10.0 Description: Gaussian width, depending on the m/z position.
    The higher the value, the wider the peak and therefore the wider the gaussian.
    Param: b'use_ppm_tolerance' Value: false Description: If true, instead of the gaussian_width value, the ppm_tolerance is used. The gaussian is calculated in each step anew, so this is much slower.
    Param: b'write_log_messages' Value: false Description: true: Warn if no signal was found by the Gauss filter algorithm.

To print a simple key-value list, you can use ``asDict()``, as shown above:

.. code-block:: python
    :linenos:
    
    gf = oms.GaussFilter()
    gf.getParameters().asDict()
    
    
Types of Parameter Values
************************************************

A :py:class:`~.Param` object can hold many parameters of mixed value type. Above, we have seen floating point values, e.g.

.. code-block:: python
    :linenos:
    
    new_p.setValue("param2", 9.0, "This is value 9")
    
Other possible values include ``int``, ``float``, ``bytes``, ``str``, ``List[int]``, ``List[float]``, ``List[bytes]`` (aka StringList).
E.g.

.. code-block:: python
    :linenos:
    
    p = oms.Param()
    p.setValue("p_float", 4.0, "This is a float")
    p.setValue("p_int", 5, "This is an integer")
    p.setValue("p_string", "myvalue", "This is a string")
    p.setValue("p_stringlist", [b"H:+:0.6", b"Na:+:0.2", b"K:+:0.2"], "This is a StringList")
    p.setValue("p_floatlist", [1.0, 2.0, 3.0], "This is a list of floats")
    p.setValue("p_intlist", [1, 2, 3], "This is a list of integers")
    
    
Restrictions(=Validity) of Parameter Values
******************************************************* 
    
For certain types of values, pyOpenMS supports restrictions,
e.g. for single strings only a restricted set of values may be allowed.
Also, for floats/ints only a restricted interval of numbers may be valid.

Usually, these restrictions are set by the OpenMS algorithm/class which hands out the parameters.
Then, if you provide invalid values via ``setParameters``, the algorithm will throw an exception.

It is usually interesting to inspect the restrictions to know what methods a class supports, e.g. see below for an example
using a GaussFilter and the Normalizer.

In theory, you can create your own restrictions. Usually this is done when defining the algorithm in C++ and is out of scope here.

E.g.

.. code-block:: python
    :linenos:

    gf = oms.GaussFilter()
    gfp = gf.getParameters()
    gfp.getValidStrings("use_ppm_tolerance")  ## yields [b'true', b'false']
    
    gfp.setValue(b"use_ppm_tolerance", "maybe") ## is invalid but setValue does not complain
    ##  ... until you actually set the parameters:
    try:
      gf.setParameters(gfp)   ## --> throws a RuntimeError
    except RuntimeError as e:
      print(f"RuntimeError: {str(e)}")
      ## prints `GaussFilter: Invalid string parameter value 'maybe' for parameter 'use_ppm_tolerance' given! Valid values are: 'true,false'.`
         
    
    nor = oms.Normalizer()
    norp = nor.getParameters()
    norp.getValidStrings("method")  ## yields [b'to_one', b'to_TIC']
    norp.setValue("method", "to_TIC") ## pick the 'to_TIC' method
    nor.setParameters(norp)
    # ... now run the Normalizer ...
 

Unfortunately, it is not possible to retrieve the valid ranges for floats and ints, if they have been set via the pyOpenMS API (yet).
However, one can look at either the documentation of the class in pyOpenMS docs. There will be a link to the C++ version which contains the
restrictions (if any) of all parameters of a class.
Alternatively, you can simply write the parameters to an INI file (also called :py:class:`~.ParamXMLFile`), which is a special XML file format which OpenMS uses to store parameters.

E.g.

.. code-block:: python
    :linenos:
    
    pphr = oms.PeakPickerHiRes()

    px = oms.ParamXMLFile()
    px.store("tmp.ini", pphr.getParameters())  ## store PeakPickerHiRes params (or any Param object you like)

    ## either look at the file in Python, or open it in an Editor of your choice
    print(open('tmp.ini').read())    

The INI file looks something like this (shortened):

.. code-block:: xml
    :linenos:
    
    <?xml version="1.0" encoding="ISO-8859-1"?>
    <PARAMETERS version="1.7.0" xsi:noNamespaceSchemaLocation="https://raw.g.../Param_1_7_0.xsd" xmlns:xsi="...">
      <ITEM name="signal_to_noise" value="0.0" type="double" description="Minimal signal... SNT estimation!)" required="false" advanced="false" restrictions="0.0:" />
      <ITEM name="spacing_difference_gap" value="4.0" type="double" description="The extension ... chromatograms." required="false" advanced="true" restrictions="0.0:" />
    ...

Any parameter which has restrictions on its value (strings, ints and floats) will have a ``restrictions`` attribute.
In the above example, the restriction on the ``signal_to_noise`` parameter are ``restrictions="0.0:"``, i.e. only the lower bound is restricted to 0.0. The upper bound can be any value larger than 0.

    
      