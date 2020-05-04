#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 16:20:10 2020

@author: ceballos
"""

import numpy as np

def ecdf(data):
    """Compute the Empirical Cumulative Distribution Function
         or from statsmodels.distributions.empirical_distribution import ECDF
    """
    x = np.sort(data)
    n = x.size
    y = np.arange(1, n + 1) / n
    return x, y
