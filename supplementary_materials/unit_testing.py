#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 12:10:56 2026

@author: bani4ki
"""
#Unit testing of minimum two components

import unittest
from my_program import calculate_gc


class TestCalculateGC(unittest.TestCase):

    def test_normal_sequence(self):
        self.assertEqual(calculate_gc("GCGCAAAT"), 50.0)

    def test_no_gc(self):
        self.assertEqual(calculate_gc("AAAAAA"), 0.0)

    def test_all_gc(self):
        self.assertEqual(calculate_gc("GCGCGC"), 100.0)


if __name__ == "__main__":
    unittest.main()