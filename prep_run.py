# -*- coding: utf-8 -*-
import os as os

vec=range(148,152)+range(170,174)+range(178,186)+[204,237]
for ii in vec:
    print 'processing run '+str(ii)
    os.system('python ./run'+str(ii)+'/preprocess/rect_ana.py')
    os.system('python ./run'+str(ii)+'/preprocess/initial.py')
    os.system('python ./run'+str(ii)+'/preprocess/sponge.py')

