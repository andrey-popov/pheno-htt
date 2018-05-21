#!/usr/bin/env python

"""Computes convolution of PDF for an s-channel process.

Differential cross section in sHat can be computed as
  d[sigma(sHat)] / d[sHat]
    = Int[sigmaHat(sHat) f(x1) f(x2) delta(x1 x2 s - sHat) d[x1] d[x2]]
    = sigmaHat(sHat) / s * Int[f(x1) f(sHat / (s x1)) d[x1] / x1]
Integration limits in the last expression are x1 = sHat / s .. 1, and it
is convenient to express the integrand as a function of ln(x1).

This script computes the integral from the last expression (with the 1/s
factor included) using gluon PDF.  The factorization scale is set at
sqrt(sHat) / 2.  The computation is done for various values of sHat.

Results are saved in a NumPy file.  They are represented with an array
of shape (2, n), whose first row containts values of sHat and the second
row gives values of the integral at these sHat.
"""

import argparse
import math
import warnings

import numpy as np
from scipy.integrate import quad

import lhapdf


def pdf_prod(lnX, pdf, scale2, sHatFrac):
    """Compute f(x) * f(sHat / (s * x))."""
    
    x = math.exp(lnX)
    return pdf.xfxQ2(21, x, scale2) * pdf.xfxQ2(21, sHatFrac / x, scale2) / sHatFrac


if __name__ == '__main__':
    
    argParser = argparse.ArgumentParser(epilog=__doc__)
    argParser.add_argument(
        '--pdf', default='PDF4LHC15_nlo_30_pdfas',
        help='PDF set'
    )
    argParser.add_argument(
        '--sqrt-s', dest='sqrtS', type=float, default=13e3,
        help='Collision energy'
    )
    argParser.add_argument(
        '--min-sqrt-s-hat', dest='minSqrtSHat', type=float, default=100.,
        help='Minimal value for sqrt(s-hat)'
    )
    argParser.add_argument(
        '--max-sqrt-s-hat', dest='maxSqrtSHat', type=float, default=5000.,
        help='Maximal value for sqrt(s-hat)'
    )
    argParser.add_argument(
        '-n', type=int, default=5000,
        help='Number of points in s-hat to compute'
    )
    argParser.add_argument(
        '-o', '--output', default='sHatPDF.npy',
        help='Name for output .npy file'
    )
    args = argParser.parse_args()
    
    s = args.sqrtS ** 2
    
    lhapdf.setVerbosity(0)
    
    
    results = np.empty((2, args.n))
    results[0] = np.geomspace(
        args.minSqrtSHat ** 2, args.maxSqrtSHat ** 2, num=args.n
    )
    
    pdf = lhapdf.mkPDF(args.pdf, 0)
    
    for i, sHat in enumerate(results[0]):
        sHatFrac = sHat / s
        scale2 = sHat / 4
        
        # Report all warnings in the integration below
        with warnings.catch_warnings():
            warnings.simplefilter('always')
            
            integral = quad(
                pdf_prod, math.log(sHatFrac), 0.,
                args=(pdf, scale2, sHatFrac), epsrel=1e-5, epsabs=0., limit=100
            )
            results[1, i] = integral[0] / s
    
    
    # Save results in a numpy file
    np.save(args.output, results)
