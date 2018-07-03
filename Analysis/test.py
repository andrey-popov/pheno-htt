#!/usr/bin/env python

"""Tests machinery from parton cross sections to limits."""

from math import isclose

import hmssm
from spectrum import RecoMtt
import statscan


if __name__ == '__main__':
    
    mA = 500.
    tanbeta = 5.
    
    significance_ref = 1.5227446521132153
    cls_ref = 0.13491426645087826
    
    
    parton_xsec = hmssm.XSecHMSSM(mA, tanbeta, 'hMSSM_13TeV.root')
    reco_mtt = RecoMtt(parton_xsec, resolution=0.1)
    calc = statscan.StatCalc(reco_mtt, 'ttbar_res10.root', 150e3)
    
    with statscan.suppress_stdout():
        significance = calc.significance()
        cls = calc.cls()
    
    
    print('Significance\n  Current:   {}\n  Reference: {}'.format(significance, significance_ref))
    print('CLs\n  Current:   {}\n  Reference: {}'.format(cls, cls_ref))
    
    if isclose(significance, significance_ref) and isclose(cls, cls_ref):
        print('\033[1;32mTest passed\033[0m')
    else:
        print('\033[1;31mTest failed\033[0m')
