# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 14:08:34 2017

@author: Nan
"""

from intervaltree import Interval, IntervalTree
t = IntervalTree()
t[4:7] = "001"
t[5:6] = "002"
print list(t.search(5, 6))[0].data