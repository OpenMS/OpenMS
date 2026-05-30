#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test that pyOpenMS doesn't affect Python's locale settings.

This test ensures that importing pyOpenMS doesn't change the Python
locale, which would affect other Python libraries.
"""

import unittest
import locale


class TestLocale(unittest.TestCase):
    """Test that pyOpenMS import doesn't affect Python locale."""

    def test_locale_preserved_after_import(self):
        """Test that locale is preserved after importing pyOpenMS."""
        # Get the locale before importing pyOpenMS
        locale_before = locale.getlocale()
        
        # Import pyOpenMS
        import pyopenms
        
        # Get the locale after importing pyOpenMS
        locale_after = locale.getlocale()
        
        # Check that the locale hasn't changed
        # Note: Some systems might have (None, None) as the default locale,
        # so we just check that it's the same as before
        self.assertEqual(
            locale_before, 
            locale_after,
            f"Locale changed from {locale_before} to {locale_after} after importing pyOpenMS"
        )
        
        # Also ensure it's not (None, None) if we started with a valid locale
        if locale_before != (None, None):
            self.assertNotEqual(
                locale_after,
                (None, None),
                "Locale was reset to (None, None) after importing pyOpenMS"
            )


if __name__ == '__main__':
    unittest.main()
