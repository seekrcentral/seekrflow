"""
Unit and regression test for the seekrflow package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import seekrflow


def test_seekrflow_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "seekrflow" in sys.modules
