#!/usr/bin/env python

"""Plots mtt in signal for varying degree of smearing."""

import os

import numpy as np

import matplotlib as mpl
mpl.use('agg')
from matplotlib import pyplot as plt

import lhapdf

import hmssm
from spectrum import RecoMtt


if __name__ == '__main__':
    
    mpl.rc('xtick', top=True, direction='in')
    mpl.rc('ytick', right=True, direction='in')
    mpl.rc('axes', labelsize='large')
    mpl.rc('axes.formatter', limits=[-2, 4], use_mathtext=True)
    
    fig_dir = 'fig'
    
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    
    lhapdf.setVerbosity(0)
    
    
    mA, tanbeta = 500., 1.
    parton_xsec = hmssm.XSecHMSSM(mA, tanbeta, 'hMSSM_13TeV.root')
    reco_mtt = RecoMtt(parton_xsec)
    
    
    fig = plt.figure()
    axes = fig.add_subplot(111)
    
    mtt = np.linspace(300., 700., num=500)
    xsec_nosmear = np.empty_like(mtt)
    
    for i in range(len(mtt)):
        xsec_nosmear[i] = reco_mtt.xsec_no_smear(mtt[i])
    
    axes.plot(mtt, xsec_nosmear, c='black', label='No smearing')
    
    for resolution in [1e-2, 0.05, 0.1, 0.2]:
        reco_mtt.resolution = resolution
        xsec = np.empty_like(mtt)
        
        for i in range(len(mtt)):
            xsec[i] = reco_mtt.xsec(mtt[i])
        
        axes.plot(mtt, xsec, label='Smearing {:g}%'.format(resolution * 100.))
    
    axes.axhline(0., c='black', ls='dashed', lw=0.8)
    
    axes.legend()
    axes.margins(x=0.)
    axes.set_xlabel('$m_{t\\bar t}$ [GeV]')
    axes.set_ylabel('Differential cross section [pb/GeV]')
    
    axes.text(
        1., 1., '$m_A = {:g}$ GeV, $\\tan\\beta = {:g}$'.format(mA, tanbeta),
        ha='right', va='bottom', transform=axes.transAxes
    )
    
    fig.savefig(os.path.join(fig_dir, 'resolution.pdf'))
    plt.close(fig)
