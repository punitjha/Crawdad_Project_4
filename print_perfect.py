#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 10:47:01 2017

@author: mintuser
"""

import matplotlib.pyplot as plt
import numpy as np

def print_perfect (over_mat):
    for row in range(over_mat.shape[0]):
        for col in range(over_mat.shape[1]):
            print('{:11.8f}'.format(over_mat[row,col]),end=' ')
        print()