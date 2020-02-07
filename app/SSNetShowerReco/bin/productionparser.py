from __future__ import print_function
from __future__ import division
import json
import numpy as np
from math import sqrt


def FidVol(x,y,z):

        if x < 15:        return False
        if x > 256-15:    return False
        if y < -116.5+15: return False
        if y > 116.5-25:  return False
        if z < 15:        return False
        if z > 1000:      return False

        return True

filelist = open("filelist.txt").read().splitlines()

for i,f in enumerate(filelist):
    print(f)
    with open(f,"r") as read_file:
        data = json.load(read_file)
        for p in data['entries']:
            if( len(p['vertex_pos'])>=0):
                print('Run: ' ,p['run'], ' Subrun: ' ,p['subrun'], 'Event: ' ,p['event'])
                print('Vertex position: ', p['vertex_pos'])
                print('Shower Energies: ', p['shower_energies'])
                print('Sum Qs: ', p['shower_sumQs'])
                print('Length ', p['shower_shlengths'])
                print('True Energy: ', p['true_electron_energy'])
                print('True Pos: ', p['true_vertex_sce'])
                print('')
