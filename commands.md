#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 14:38:24 2022

@author: François Clapuyt
"""

# Test

in_dtm="/home/elic/fclapuyt/SedInConnect_2.3/inputs/Otemma/Otemma_filledWL2019.tif"
in_cs=2
out_w="/home/elic/fclapuyt/SedInConnect_2.3/outputs/Otemma/otemma_weight.tif"
out_ic="/home/elic/fclapuyt/SedInConnect_2.3/outputs/Otemma/otemma_ic.tif"

CavalliConnectivityout(in_dtm, in_cs, out_w, out_ic)