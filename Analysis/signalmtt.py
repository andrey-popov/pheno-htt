import math

import numpy as np

from scipy.integrate import simps
from scipy.interpolate import interp1d

import lhapdf


class SignalMtt:
    """A class to evaluate diff. cross section in mtt for H/A -> tt.
    
    Applies selection efficiency and smearing to parton-level cross
    section.
    """
    
    def __init__(
        self, partonCalc, resolution=0.2, sHatPDFFile='sHatPDF.npy',
        pdfLabel='PDF4LHC15_nlo_30_pdfas'
    ):
        """Initialize.
        
        Arguments:
            partonCalc:  An object to compute parton-level cross
                sections.
            resolution:  Relative resolution in mtt.
            sHatPDFFile:  NumPy file with array that provides PDF
                convolution as a function of parton-level s.
            pdfLabel:  PDF set from which strong coupling constant will
                be read.
        """
        
        self.partonCalc = partonCalc
        self.resolution = resolution
        self.pdf = lhapdf.mkPDF(pdfLabel, 0)
        
        sHatPDF = np.load(sHatPDFFile)
        self.sHatPDFInterp = interp1d(sHatPDF[0], sHatPDF[1], copy=False, assume_sorted=True)
        
        # Branching ratio for targeted decays.  Set to l+jets, l = e/mu.
        self.targetBranching = 8 / 27
        
        # Scale factor to vary renormalization scale
        self.muRScaleFactor = 1.
        
        # Cross section without smearing precomputed on an adaptive
        # grid.  Needed to compute smeared cross section efficiently.
        # A NumPy array with shape (2, N) that contains values of mtt
        # and the cross section.
        self.xSecNoSmearGrid = None
        
        # An interpolation function for self.xSecNoSmearGrid
        self.xSecNoSmearInterp = None
    
    
    def build_xsec_grid(self, relTolerance=0.005):
        """Approximate cross section without smearing with a grid.
        
        Approximate the cross section with an adaptive grid so that
        integration of it can be replaced with a summation over the
        grid.  The grid is constructed so that at every two consecutive
        points the absolute change of the cross section is smaller than
        a predetermined tolerance.
        
        Whenever parameters affecting the cross section without smearing
        change, this method must be rerun.
        
        Arguments:
            relTolerance:  The fraction of the overall span of the cross
                section to be used as the tolerance.
        
        Return value:
            None.
        
        Start from a coarse uniform grid and iteratively add points at
        centres of segments until the change of the cross section over
        each segment becomes small.  Constructed grid with precomputed
        cross section is stored as self.xSecNoSmearGrid.
        """
        
        mttRange = (340., 2000.)
        
        # Initial uniform grid
        numPoints = max(100, round((mttRange[1] - mttRange[0]) / self.partonCalc.variationScale))
        mttValues = list(np.linspace(mttRange[0], mttRange[1], num=numPoints))
        xsecValues = [0.] * len(mttValues)
        
        for i in range(len(mttValues)):
            xsecValues[i] = self.xsec_no_smear(mttValues[i])
        
        
        tolerance = relTolerance * (max(xsecValues) - min(xsecValues))
        
        
        # Iteratively adjust the grid adding more points where needed
        i = 0
        
        while i < len(mttValues) - 1:
            
            # Check how well the function is approximated with a linear
            # extrapolation.  To do it, compute the vertical distance
            # between the interpolated and the actual value of the
            # function at the centre of the segment.
            meanMtt = (mttValues[i] + mttValues[i + 1]) / 2
            xsecMeanMtt = self.xsec_no_smear(meanMtt)
            deviation = abs(xsecMeanMtt - (xsecValues[i] + xsecValues[i + 1]) / 2)
            
            # If the overall change of the function over the current
            # segment is larger than the tolerance or the function
            # deviates to much from a linear interpolation, add the
            # middle point to the grid.  The condition imposed on the
            # deviation of the middle point also means that the largest
            # change in the function on all three points is less than
            # the tolerance.
            if abs(xsecValues[i + 1] - xsecValues[i]) > tolerance or deviation > tolerance / 2:
                mttValues.insert(i + 1, meanMtt)
                xsecValues.insert(i + 1, xsecMeanMtt)
            else:
                i += 1
        
        self.xSecNoSmearGrid = np.array([mttValues, xsecValues])
        self.xSecNoSmearInterp = None
    
    
    def scale(self, mtt):
        """Evaluate renormalization scale, given parton-level mtt."""
        
        return mtt / 2 * self.muRScaleFactor
    
    
    def selection_efficiency(self, mtt, subprocess):
        """Compute efficiency of event selection.
        
        Arguments:
            mtt:  Parton-level mtt, in GeV.
            subprocess:  String 'Res' or 'Int' to choose resonant part
                or interference.
        
        Return value:
            Efficiency.
        
        Computed efficiency does not include efficiencies of lepton
        identification or b-tagging, but it includes the branching
        ratio for targeted decays.
        """
        
        # Fitted parameters of a log-cubic approximation
        if subprocess == 'Res':
            coeffs = [-0.1166857, 2.19917117, -13.58990087, 27.78332692]
        elif subprocess == 'Int':
            coeffs = [-0.05877867, 1.03660773, -5.91517033, 11.11336388]
        else:
            raise RuntimeError('Do not recognize subprocess "{}".'.format(subprocess))
        
        return np.polyval(coeffs, math.log(mtt)) * self.targetBranching
    
    
    def xsec(self, mtt, numSigma=3):
        """Compute differential cross section.
        
        Effects of reconstruction are approximated by smearing the
        parton-level mtt (with selection efficiency applied).
        
        A grid of precomputed values of cross section without smearing
        is utilized to improve the performance.  If any parameters
        affected the parton-level cross section or the selection
        efficiency have changed, the grid must be rebuilt using method
        build_xsec_grid.
        
        Arguments:
            mtt:  Value of smeared mtt, in GeV.
            numSigma:  Defines truncation for Gaussian kernel.
        
        Return value:
            Full differential cross section in mtt, in pb / GeV.
        """
        
        # Precompute cross section without smearing if not done yet
        if self.xSecNoSmearGrid is None:
            self.build_xsec_grid()
        
        
        # Find the integration window
        halfWidth = numSigma * self.resolution * mtt
        start = np.searchsorted(self.xSecNoSmearGrid[0], mtt - halfWidth, 'left')
        end = np.searchsorted(self.xSecNoSmearGrid[0], mtt + halfWidth, 'right')
        
        
        # Perform a convolution of the cross section without smearing
        # with a Gaussian kernel.  Will use different algorithms
        # depending on how many precomputed points the integration
        # window contains.
        if (end - start) / numSigma > 20:
            # Compute weights for the Gaussian kernel
            weights = np.empty(end - start)
            
            for i in range(end - start):
                curMtt = self.xSecNoSmearGrid[0, start + i]
                sigma = self.resolution * curMtt
                weights[i] = 1 / sigma * math.exp(-0.5 * ((mtt - curMtt) / sigma) ** 2)
            
            weights /= math.sqrt(2 * math.pi)
            
            # Correct for truncation in the kernel.  This is not exact
            # because the variance of the kernel is not constant.
            weights /= math.erf(numSigma / math.sqrt(2))
            
            return simps(
                self.xSecNoSmearGrid[1, start:end] * weights,
                self.xSecNoSmearGrid[0, start:end], 'first'
            )
        
        else:
            # The grid of precomputed points is too coarse for the
            # assumed resolution.  The cross section without smearing
            # can be approaximated using linear interpolation as it does
            # not change rapidly within the window, but more poitns are
            # needed for the convolution with the Gaussian kernel.
            
            # Construct the interpolation object if not done yet
            if self.xSecNoSmearInterp is None:
                self.xSecNoSmearInterp = interp1d(
                    self.xSecNoSmearGrid[0], self.xSecNoSmearGrid[1],
                    copy=False, assume_sorted=True, bounds_error=False, fill_value=0.
                )
            
            x = np.linspace(mtt - halfWidth, mtt + halfWidth, num=101)
            y = np.empty_like(x)
            
            for i in range(len(x)):
                sigma = self.resolution * x[i]
                weight = 1 / sigma * math.exp(-0.5 * ((mtt - x[i]) / sigma) ** 2)
                y[i] = self.xSecNoSmearInterp(x[i]) * weight
            
            y /= math.sqrt(2 * math.pi) * math.erf(numSigma / math.sqrt(2))
            
            return simps(y, x)
    
    
    def xsec_no_smear(self, mtt):
        """Compute differential cross section without smearing.
        
        Apply convolution with PDF and selection efficiencies, but not
        the smearing.
        
        Arguments:
            mtt:  Parton-level mtt, in GeV.
        
        Return value:
            Differential cross section in mtt, in pb / GeV, up to but
            not including smearing.
        """
        
        sHat = mtt ** 2
        alphaS = self.pdf.alphasQ(self.scale(mtt))
        
        # Evaluate parts that are different for the two subprocesses
        res = self.partonCalc.xsec_res(sHat, alphaS) * self.selection_efficiency(mtt, 'Res')
        int_ = self.partonCalc.xsec_int(sHat, alphaS) * self.selection_efficiency(mtt, 'Int')
        
        # Convolute with PDF.  The last two terms appear from
        # translation of d[sigma] / d[sqrt(sHat)] to d[sigma] / d[sHat].
        xsec = (res + int_) * self.sHatPDFInterp(sHat) * 2 * mtt
        
        return xsec
