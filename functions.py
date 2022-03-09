#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 11:38:46 2022

@author: franz
"""

def __init__(self)
def initUI(self)
    
# to enable drag and drop   
def eventFilter(self, object, event)
def AddTarget(self, state)
def AddSinks(self, state)
def WInput(self, state)
def buttonClicked(self)
def selectFile1(self)#input DTM
def saveRough(self)#save RI
def saveWeight(self)#isave W
def selectFile2(self)#input weight raster

        
# sinks detection and function
def sinks(self, dtm_f, sink_shp)

# Roughness index calculation
def rw_cavalli(self, dtm_r, f_s, ri_out, w_out, sk_flag = 0)
        
# CONNECTIVITY TO THE OUTLET
def CavalliConnectivityout (self, dtm_f, c_s, w_f, out_ic, sink_flag = 0)

# CONNECTIVITY TO TARGETS
def CavalliConnectivitytg (self, dtm_f, c_s, tg_f, w_f, out_ic_tg, sink_flag_tg = 0)

    
def OK_rw_conn_outlet(self)

