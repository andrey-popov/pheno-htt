#!/usr/bin/env python

"""Plots significance and exclusion in (mA, tan(beta)) plane.

Data are read from a CSV file.
"""

import argparse
import os

import numpy as np
from scipy.interpolate import RectBivariateSpline

import matplotlib as mpl
mpl.use('agg')
from matplotlib import pyplot as plt


if __name__ == '__main__':
    
    argParser = argparse.ArgumentParser(epilog=__doc__)
    argParser.add_argument('dataFile', help='CSV file with significance')
    argParser.add_argument(
        '-l', '--label', default=';',
        help='Label consisting of two parts separated by semicolon'
    )
    argParser.add_argument(
        '-o', '--output', default='fig/significance.pdf',
        help='Name for output figure file'
    )
    args = argParser.parse_args()
    
    labels = args.label.split(';')
    
    if len(labels) != 2:
        raise RuntimeError('Failed to parse label "{}".'.format(args.label))
    
    
    mpl.rc('xtick', top=True, direction='in')
    mpl.rc('ytick', right=True, direction='in')
    mpl.rc('axes', labelsize='large')
    mpl.rc('axes.formatter', limits=[-2, 4], use_mathtext=True)
    
    figDir = os.path.dirname(args.output)
    
    if figDir and not os.path.exists(figDir):
        os.makedirs(figDir)
    
    
    # Read CSV file
    data = []
    
    with open(args.dataFile) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parsedRow = [float(token) for token in line.split(',')]
            data.append(tuple(parsedRow))
    
    
    # Make sure data are sorted (first in mA, then in tan(beta)) and
    # convert to a NumPy array
    data.sort()
    data = np.asarray(data)
    
    
    # Fits executed while evaluating the significance of CLs could have
    # failed, leading to infinite and NaN values.  This seems to only
    # happen for very large significances.  Replace the problematic
    # values with large dummy values for significance and zero CLs.
    # Also replace valid but large values of significance with the same
    # placeholders, which is needed to avoid instabilities in the
    # interpolation below.
    maxSignificance = 10.
    
    for row in data:
        if not np.isfinite(row[2]) or row[2] > maxSignificance:
            row[2] = maxSignificance + 1
        if not np.isfinite(row[3]):
            row[3] = 0.
    
    
    # Deduce the number of points in tan(beta) by finding the first row
    # in which mA differs from data[0, 0].  This allows to reconstruct
    # the grid used for mA and tan(beta), which are denoted as x and y.
    numPointsY = np.searchsorted(data[:, 0], data[0, 0], 'right')
    pointsX = data[::numPointsY, 0]
    pointsY = data[:numPointsY, 1]
    
    
    # Upsample the significance by interpolating between evaluated
    # points.  This is needed to construct smooth contours.
    significanceInterp = RectBivariateSpline(
        pointsX, pointsY,
        np.reshape(data[:, 2], (-1, numPointsY))
    )
    
    upsampleFactor = 10
    pointsXUpsampled = np.linspace(pointsX[0], pointsX[-1], num=len(pointsX) * upsampleFactor)
    pointsYUpsampled = np.linspace(pointsY[0], pointsY[-1], num=len(pointsY) * upsampleFactor)
    significanceUpsampled = np.transpose(significanceInterp(pointsXUpsampled, pointsYUpsampled))
    
    
    # Similarly, upsample CLs
    clsInterp = RectBivariateSpline(
        pointsX, pointsY,
        np.reshape(data[:, 3], (-1, numPointsY))
    )
    clsUpsampled = np.transpose(clsInterp(pointsXUpsampled, pointsYUpsampled))
    
    
    # Significance will be plotted with imshow, which will treat each
    # bin as a pixel of an image.  Find the span of the image in
    # physical coordinates along each axis.  Take into account that
    # imshow does not clip boundary "pixels", so the span must be
    # extended by half a pixel size on each side.
    stepX = (pointsXUpsampled[-1] - pointsXUpsampled[0]) / len(pointsXUpsampled)
    stepY = (pointsYUpsampled[-1] - pointsYUpsampled[0]) / len(pointsYUpsampled)
    extent = (
        pointsXUpsampled[0] - stepX / 2, pointsXUpsampled[-1] + stepX / 2,
        pointsYUpsampled[0] - stepY / 2, pointsYUpsampled[-1] + stepY
    )
    
    
    # Plot significance together with characteristic contours
    fig = plt.figure()
    axes = fig.add_subplot(111)
    
    colormap = plt.get_cmap('viridis')
    # colormap.set_over('white')
    
    image = axes.imshow(
        significanceUpsampled, aspect='auto', origin='lower', interpolation='bicubic',
        extent=extent, cmap=colormap, vmin=0., vmax=maxSignificance
    )
    fig.colorbar(image, fraction=0.05, pad=0.02, label='Significance [$\\sigma$]', spacing='proportional')
    
    contourSet = axes.contour(
        pointsXUpsampled, pointsYUpsampled, significanceUpsampled,
        [1., 3., 5.], colors='white', zorder=1.5
    )
    axes.clabel(contourSet, fmt='%g $\\sigma$')
    
    contourCLs = axes.contour(
        pointsXUpsampled, pointsYUpsampled, clsUpsampled,
        [0.05], colors='red', zorder=1.5
    )
    axes.clabel(contourCLs, fmt='95%% CL excl.')
    
    axes.set_xlim(pointsXUpsampled[0], pointsXUpsampled[-1])
    axes.set_ylim(pointsYUpsampled[0], pointsYUpsampled[-1])
    
    axes.set_xlabel('$m_A$ [GeV]')
    axes.set_ylabel('$\\tan\/\\beta$')
    
    axes.text(0., 1.005, labels[0], ha='left', va='bottom', transform=axes.transAxes)
    axes.text(1., 1.005, labels[1], ha='right', va='bottom', transform=axes.transAxes)
    
    fig.savefig(args.output)
