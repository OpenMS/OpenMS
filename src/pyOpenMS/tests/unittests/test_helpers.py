#!/usr/bin/env python
# -*- coding: utf-8  -*-

## ----------------------------------------------------------------------------
## $Maintainer: $
## $Authors: Hannes Roest, Timo Sachsenberg, axelwalter,
##           Samuel Wein, Uwe Schmitt, Joshua Charkow,
##           Nikos Patikas, Chris Bielow, Julianus Pfeuffer,
##           Oliver Alka, Stephan Aiche $
## ----------------------------------------------------------------------------
"""
Shared test helpers for pyOpenMS unit tests.

This module contains helper functions used across multiple test files
after splitting test000.py. These helpers provide common testing utilities
for validating interfaces and generating test data.
"""

from __future__ import print_function

import pyopenms
import sys
import time
import numpy as np
from functools import wraps
from typing import Tuple

# Python 2/3 compatibility
try:
    long
except NameError:
    long = int


def _testStrOutput(input_str):
    """
    Test helper to validate string output type across Python versions.
    
    Args:
        input_str: String to validate
        
    Asserts:
        - Python 2: input_str is unicode
        - Python 3: input_str is str
    """
    if sys.version_info[0] < 3:
        assert isinstance(input_str, unicode)
    else:
        assert isinstance(input_str, str)


def report(f):
    """
    Decorator to print test function name before execution.
    
    Args:
        f: Test function to decorate
        
    Returns:
        Wrapped function that prints "run <function_name>" before execution
    """
    @wraps(f)
    def wrapper(*a, **kw):
        print("run ", f.__name__)
        f(*a, **kw)
    return wrapper


def _testMetaInfoInterface(what):
    """
    Test helper to validate MetaInfoInterface implementation.
    
    Tests the following MetaInfoInterface methods:
    - getKeys()
    - getMetaValue()
    - setMetaValue()
    - metaValueExists()
    - removeMetaValue()
    - clearMetaInfo()
    
    Args:
        what: Object implementing MetaInfoInterface
        
    Note:
        Preserves essential meta values (like "rank" for PeptideHit) before
        testing and restores them after clearMetaInfo().
    """
    # Store essential meta values (like rank for PeptideHit) before testing
    essential_keys = [b"rank"]
    preserved_values = {}
    for key in essential_keys:
        if what.metaValueExists(key):
            preserved_values[key] = what.getMetaValue(key)

    what.setMetaValue("key", 42)
    what.setMetaValue("key2", 42)

    keys = []
    what.getKeys(keys)
    assert len(keys) and all(isinstance(k, bytes) for k in keys)
    # Check that our set keys exist and have correct values (keys are unordered)
    assert b"key" in keys
    assert b"key2" in keys
    assert what.getMetaValue(b"key") == 42
    assert what.getMetaValue(b"key2") == 42

    assert what.metaValueExists("key")
    what.removeMetaValue("key")

    keys = []
    what.getKeys(keys)
    assert b"key2" in keys
    assert what.getMetaValue(b"key2") == 42

    what.clearMetaInfo()
    
    # Restore essential meta values after clearMetaInfo
    for key, value in preserved_values.items():
        what.setMetaValue(key, value)
    
    keys = []
    what.getKeys(keys)
    # Check that only essential keys remain (if any were preserved)
    expected_len = len(preserved_values)
    assert len(keys) == expected_len


def _testUniqueIdInterface(what):
    """
    Test helper to validate UniqueIdInterface implementation.
    
    Tests the following UniqueIdInterface methods:
    - hasInvalidUniqueId()
    - hasValidUniqueId()
    - ensureUniqueId()
    - getUniqueId()
    - clearUniqueId()
    - setUniqueId()
    
    Args:
        what: Object implementing UniqueIdInterface
    """
    assert what.hasInvalidUniqueId()
    assert not what.hasValidUniqueId()
    assert what.ensureUniqueId()
    assert isinstance(what.getUniqueId(), (int, long))
    assert what.getUniqueId() > 0
    assert not what.hasInvalidUniqueId()
    assert what.hasValidUniqueId()

    what.clearUniqueId()
    assert what.getUniqueId() == 0
    assert what.hasInvalidUniqueId()
    assert not what.hasValidUniqueId()

    assert what.ensureUniqueId()
    assert isinstance(what.getUniqueId(), (int, long))
    assert what.getUniqueId() > 0
    assert not what.hasInvalidUniqueId()
    assert what.hasValidUniqueId()

    what.setUniqueId(1234)
    assert what.getUniqueId() == 1234


def _testProgressLogger(ff):
    """
    Test helper to validate ProgressLogger functionality.
    
    Tests the following ProgressLogger methods:
    - setLogType()
    - getLogType()
    - startProgress()
    - setProgress()
    - endProgress()
    
    Args:
        ff: Object implementing ProgressLogger interface
    """
    ff.setLogType(pyopenms.LogType.NONE)
    assert ff.getLogType() == pyopenms.LogType.NONE
    ff.startProgress(0, 3, "label")
    ff.setProgress(0)
    ff.setProgress(1)
    ff.setProgress(2)
    ff.setProgress(3)
    ff.endProgress()


def _testParam(p):
    """
    Test helper to validate Param class functionality.
    
    Tests the following Param methods:
    - asDict()
    - size()
    - keys()
    - values()
    - items()
    - getValue() / __getitem__
    - setValue() / __setitem__
    - getDescription()
    - getTags()
    - setSectionDescription()
    - getSectionDescription()
    - get()
    - exists()
    - addTag()
    - addTags()
    - clearTags()
    - insert()
    - copy()
    - update()
    
    Args:
        p: Param object to test
    """
    assert p == p

    dd = p.asDict()
    assert len(dd) == p.size()
    assert isinstance(dd, dict)

    for k in p.keys():
        value = p[k]
        p[k] = value
        p.update(p)
        p.update(p.asDict())
        assert p[k] == value
        desc = p.getDescription(k)
        tags = p.getTags(k)
        p.setValue(k, value, desc, tags)
        p.setValue(k, value, desc)
        assert p.exists(k)
        # only set the section description if there are actually two or more sections
        if len(k.split(b":")) < 2:
            continue
        f = k.split(b":")[0]
        p.setSectionDescription(f, k)
        # FIXME: keys inside maps are not yet properly decoded
        assert p.getSectionDescription(f) == k.decode()

        assert p.get(k) is not None

    assert len(p.values()) == len([p[k] for k in p.keys()])
    assert sorted(p.items()) == sorted((k, p[k]) for k in p.keys())

    assert not p.exists("asdflkj01231321321v")
    p.addTag(k, "a")
    p.addTags(k, [b"", b"c"])
    assert sorted(p.getTags(k)) == [b"", b"a", b"c"]
    p.clearTags(k)
    assert p.getTags(k) == []

    pn = pyopenms.Param()
    pn.insert("master:", p)
    assert pn.exists(b"master:" + k)

    p1 = pn.copy("master:", True)
    assert p1 == p

    p1.update(p)
    p1.update(p, 0)
    p1.update(p, 1)


def generate_random_ranges(exp: pyopenms.MSExperiment, 
                          n_ranges: int,
                          rt_width: float,
                          mz_width: float) -> np.ndarray:
    """
    Generate random ranges within the experiment bounds.
    
    Used for testing extraction and aggregation performance of MSExperiment.
    
    Args:
        exp: MSExperiment object
        n_ranges: Number of ranges to generate
        rt_width: Width of RT window
        mz_width: Width of m/z window
        
    Returns:
        numpy array of shape (n_ranges, 4) with [mz_min, mz_max, rt_min, rt_max]
    """
    # Set the seed for reproducibility
    np.random.seed(4711)
    
    min_mz = exp.getMinMZ()
    max_mz = exp.getMaxMZ()
    min_rt = exp.getMinRT()
    max_rt = exp.getMaxRT()    
    print(f"spectra min_mz: {min_mz}, max_mz: {max_mz}, min_rt: {min_rt}, max_rt: {max_rt}")

    # Generate random centers
    mz_centers = np.random.uniform(min_mz + mz_width/2, max_mz - mz_width/2, n_ranges)
    rt_centers = np.random.uniform(min_rt + rt_width/2, max_rt - rt_width/2, n_ranges)
    
    # Create ranges array
    ranges = np.zeros((n_ranges, 4))
    ranges[:, 0] = mz_centers - mz_width/2  # mz_min
    ranges[:, 1] = mz_centers + mz_width/2  # mz_max
    ranges[:, 2] = rt_centers - rt_width/2  # rt_min
    ranges[:, 3] = rt_centers + rt_width/2  # rt_max
    
    return ranges


def extraction_performance_test(exp: pyopenms.MSExperiment, ms_level: int) -> Tuple[float, float]:
    """
    Run performance test for both aggregateFromMatrix and extractXICsFromMatrix.
    
    Tests extraction performance by generating 1,000,000 random ranges and
    measuring the time taken for aggregate and XIC extraction operations.
    
    Args:
        exp: MSExperiment object
        ms_level: MS level to extract from
        
    Returns:
        Tuple of (aggregate_time, xic_time) in seconds
    """
    # Generate 1,000,000 random ranges
    print("\nGenerating random ranges...")
    ranges = generate_random_ranges(exp, 
                                   n_ranges=1_000_000,
                                   rt_width=60.0,
                                   mz_width=0.01)
    
    # Convert to MatrixDouble
    ranges_matrix = pyopenms.MatrixDouble.fromNdArray(ranges)
    
    # Time aggregateFromMatrix
    print("Running aggregateFromMatrix...")
    start_time = time.time()
    _ = exp.aggregateFromMatrix(ranges_matrix, ms_level, b"sum")
    aggregate_time = time.time() - start_time
    
    # Time extractXICsFromMatrix
    print("Running extractXICsFromMatrix...")
    start_time = time.time()
    _ = exp.extractXICsFromMatrix(ranges_matrix, ms_level, b"sum")
    xic_time = time.time() - start_time
    
    return aggregate_time, xic_time
