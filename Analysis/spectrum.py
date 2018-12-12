"""
Exports an abstract base class to describe parton-level cross section
for gg -> S -> tt and a class to emulate reconstructed mtt distribution
starting from the parton-level cross section.
"""

import abc
import math

import numpy as np

from scipy.integrate import simps
from scipy.interpolate import interp1d

import lhapdf


lhapdf.setVerbosity(0)


class PartonXSec(abc.ABC):
    
    """Base class to compute cross sections for gg -> S -> tt.
    
    Provides an interface for computation of parton-level cross sections
    for the resonant production of gg -> S -> tt and the interference
    with the SM gg -> tt background.  The cross sections are evaluated
    as functions of the square root of Mandelstam s variable, which is
    equivalent to the mtt invariant mass, and the strong coupling
    constant.
    
    Several auxiliary methods, such as computation of the partial width
    for the H -> tt decay, are also provided.
    
    Class attributes:
        mt:  Mass of the top quark (GeV).
        gF:  Fermi constant (1/GeV^2).
    
    A subclass must override methods to compute the two cross sections
    and set the typical scale at which the cross section can change
    substantially.
    """
    
    mt = 173.  # GeV
    gF = 1.166390e-05  # GeV^(-2)
    
    
    def __init__(self):
        self.var_scale_ = None
    
    
    @staticmethod
    def beta(s):
        """Compute velocity of each of the top quarks.
        
        Arguments:
            s:  Mandelstam s variable, in GeV^2.
        
        Return value:
            The velocity beta.
        """
        
        return math.sqrt(1 - 4 * PartonXSec.mt ** 2 / s)
    
    
    @staticmethod
    def loop_ampl_fermion(cp, s, mf=None):
        """Compute fermionic loop amplitude.
        
        Amplitudes associated with the fermion loop attached to a Higgs
        boson are taken from Section 2.2.1 of [1].
        [1] Branco et al., https://arxiv.org/abs/1106.0034
        
        Arguments:
            cp:  CP state of the Higgs boson, 'A' or 'H'.
            s:  Mandelstam s variable, in GeV^2.
            mf:  Mass of the fermion, in GeV.  Alternatively, can be
                None.  In this case mass of the top quark is used.
        
        Return value:
            Computed amplitude.  It is a complex number.
        """
        
        PartonXSec._check_cp(cp)
        
        if mf is None:
            mf = PartonXSec.mt
        
        tau = s / (2 * mf) ** 2
        
        if tau > 1:
            beta = math.sqrt(1 - 1 / tau)
            f = -0.25 * (math.log((1 + beta) / (1 - beta)) - math.pi * 1j) ** 2
        else:
            f = math.asin(math.sqrt(tau)) ** 2
        
        if cp == 'H':
            A = 2 * (tau + (tau - 1) * f) / tau ** 2
        else:
            A = 2 * f / tau
        
        return A
    
    
    @staticmethod
    def to_pb(xsec):
        """Convert cross section from GeV^(-2) to pb."""
        
        return xsec / 2.56819e-9
    
    
    @property
    def var_scale(self):
        """Scale at which cross section changes substantially.
        
        This is a typical scale of a change in sqrt(s), in GeV, over
        which the cross section can change substantially.  When the
        spectrum is evaluated on a uniform grid with the step size set
        to this scale, it is supposed to represent all important
        features.
        
        The scale must be set in a subclass.  Otherwise, this method
        raises an exception.
        """
        
        if self.var_scale_ is None:
            raise NotImplementedError('Property var_scale must be set in a subclass.')
        else:
            return self.var_scale_
    
    
    @var_scale.setter
    def var_scale(self, value):
        self.var_scale_ = value
    
    
    @staticmethod
    def width_tt(cp, mass, s=None, g=1.):
        """Compute partial width for the decay to tt.
        
        Use formulae for energy-dependent widths from [1].
        [1] Dicus et al., http://arxiv.org/abs/hep-ph/9404359
        
        Arguments:
            cp:  CP state of the Higgs boson, 'A' or 'H'.
            mass:  Mass of the Higgs boson, in GeV.
            s:  Mandelstam s variable, in GeV^2.  Alternatively, can be
                None.  In this case compute the width at the scale
                mass^2.
            g:  Reduced coupling to top quarks.
        
        Return value:
            Partial width for the decay H -> tt, in GeV.
        """
        
        PartonXSec._check_cp(cp)
        
        if cp == 'H':
            beta_power = 3
        else:
            beta_power = 1
        
        if s is None:
            s = mass ** 2
        
        if s < 4 * PartonXSec.mt ** 2:
            return 0.
        
        return g ** 2 * 3 * PartonXSec.gF * PartonXSec.mt ** 2 / (4 * math.pi * math.sqrt(2)) \
            * PartonXSec.beta(s) ** beta_power * s / mass
    
    
    def xsec(self, sqrt_s, alpha_s):
        """Compute cross section for gg -> S -> tt.
        
        Arguments:
            sqrt_s:  Square root of Mandelstam s variable, in GeV.
            alpha_s:  Value of the strong coupling constant.
        
        Return value:
            Computed cross section, in pb.
        """
        
        return self.xsec_res(sqrt_s, alpha_s) + self.xsec_int(sqrt_s, alpha_s)
    
    
    @abc.abstractmethod
    def xsec_res(self, sqrt_s, alpha_s):
        """Compute cross section for resonant part in gg -> S -> tt.
        
        Arguments:
            sqrt_s:  Square root of Mandelstam s variable, in GeV.
            alpha_s:  Value of the strong coupling constant.
        
        Return value:
            Computed cross section, in pb.
        
        Method must be implemented in a subclass.
        """
        
        raise NotImplementedError
    
    
    @abc.abstractmethod
    def xsec_int(self, sqrt_s, alpha_s):
        """Compute cross section for interference in gg -> S -> tt.
        
        Arguments:
            sqrt_s:  Square root of Mandelstam s variable, in GeV.
            alpha_s:  Value of the strong coupling constant.
        
        Return value:
            Computed cross section, in pb.
        
        Method must be implemented in a subclass.
        """
        
        raise NotImplementedError
    
    
    @staticmethod
    def _check_cp(cp):
        """Check if given CP state is supported.
        
        Raise an exception if given label of CP state is not 'A' or 'H'.
        
        Arguments:
            cp:  Label of the CP state.
        
        Return value:
            None.
        """
        
        if cp not in {'A', 'H'}:
            raise RuntimeError('Cannot recognize CP state "{}".'.format(cp))


class RecoMtt:
    """Class to emulate reconstructed mtt distribution for H -> tt.
    
    This class starts from a parton-level cross section provided with an
    instance of PartonXSec.  It performs convolution with PDF, applies
    selection efficiency, and smears the mtt distribution to mimic
    experimental resolution.
    """
    
    def __init__(
        self, parton_xsec, resolution=0.2, shat_pdf_file='sHatPDF.npy',
        pdflabel='PDF4LHC15_nlo_30_pdfas'
    ):
        """Initialize.
        
        Arguments:
            parton_xsec:  Object to compute parton-level cross sections.
                A subclass of PartonXSec is expected.
            resolution:  Relative resolution in mtt.
            shat_pdf_file:  NumPy file with array that provides PDF
                convolution as a function of parton-level s.
            pdflabel:  PDF set from which strong coupling constant will
                be read.
        """
        
        self.parton_xsec_ = parton_xsec
        self.resolution = resolution
        self.pdf = lhapdf.mkPDF(pdflabel, 0)
        
        shat_pdf = np.load(shat_pdf_file)
        self.shat_pdf_interp = interp1d(shat_pdf[0], shat_pdf[1], copy=False, assume_sorted=True)
        
        # Branching ratio for targeted decays.  Set to l+jets, l = e/mu.
        self.target_branching = 8 / 27
        
        # Efficiency of additional selection.  Set according to [1].
        # [1] https://github.com/andrey-popov/pheno-htt/issues/1
        self.add_sel_eff = 0.3
        
        # Scale factor to vary renormalization scale
        self.muR_scale_factor = 1.
        
        # For reasons of performance, the cross section before smearing
        # is precomputed on an adaptive grid.  This attribute will be
        # set to a NumPy array of shape (2, N) that contains values of
        # mtt and corresponding cross section.
        self.xsec_nosmear_grid = None
        
        # An interpolation function for self.xsec_nosmear_grid
        self.xsec_nosmear_interp = None
    
    
    def build_xsec_grid(self, rel_tolerance=0.005):
        """Approximate cross section before smearing with a grid.
        
        Approximate the cross section with an adaptive grid so that
        integration of it can be replaced with a summation over the
        grid.  The grid is constructed so that at every two consecutive
        points the absolute change of the cross section is smaller than
        a predetermined tolerance.
        
        Whenever parameters affecting the cross section before smearing
        change, this method must be rerun.
        
        Arguments:
            rel_tolerance:  The fraction of the overall span of the
                cross section to be used as the tolerance.
        
        Return value:
            None.
        
        Start from a coarse uniform grid and iteratively add points at
        centres of segments until the change of the cross section over
        each segment becomes small.  Constructed grid with precomputed
        cross section is stored as self.xsec_nosmear_grid.
        """
        
        mtt_range = (340., 2000.)
        
        # Initial uniform grid
        num_points = max(100, round((mtt_range[1] - mtt_range[0]) / self.parton_xsec.var_scale))
        mtt_values = list(np.linspace(mtt_range[0], mtt_range[1], num=num_points))
        xsec_values = [0.] * len(mtt_values)
        
        for i in range(len(mtt_values)):
            xsec_values[i] = self.xsec_no_smear(mtt_values[i])
        
        
        tolerance = rel_tolerance * (max(xsec_values) - min(xsec_values))
        
        
        # Iteratively adjust the grid adding more points where needed
        i = 0
        
        while i < len(mtt_values) - 1:
            
            # Check how well the function is approximated with a linear
            # extrapolation.  To do it, compute the vertical distance
            # between the interpolated and the actual value of the
            # function at the centre of the segment.
            mean_mtt = (mtt_values[i] + mtt_values[i + 1]) / 2
            xsec_mean_mtt = self.xsec_no_smear(mean_mtt)
            deviation = abs(xsec_mean_mtt - (xsec_values[i] + xsec_values[i + 1]) / 2)
            
            # If the overall change of the function over the current
            # segment is larger than the tolerance or the function
            # deviates to much from a linear interpolation, add the
            # middle point to the grid.  The condition imposed on the
            # deviation of the middle point also means that the largest
            # change in the function on all three points is less than
            # the tolerance.
            if abs(xsec_values[i + 1] - xsec_values[i]) > tolerance or deviation > tolerance / 2:
                mtt_values.insert(i + 1, mean_mtt)
                xsec_values.insert(i + 1, xsec_mean_mtt)
            else:
                i += 1
        
        self.xsec_nosmear_grid = np.array([mtt_values, xsec_values])
        self.xsec_nosmear_interp = None
    
    
    @property
    def muR_scale_factor(self):
        """Return scale factor used for renormalization scale."""
        
        return self.muR_scale_factor_
    
    
    @muR_scale_factor.setter
    def muR_scale_factor(self, scale_factor):
        """Update scale factor for renormalization scale.
        
        This invalidates cache of cross section before smearing.
        """
        
        self.muR_scale_factor_ = scale_factor
        self.xsec_nosmear_grid = None
    
    
    @property
    def parton_xsec(self):
        """Return object used to compute parton-level cross section."""
        
        return self.parton_xsec_
    
    
    @parton_xsec.setter
    def parton_xsec(self, parton_xsec_):
        """Update object used to compute parton-level cross section.
        
        This invalidates cache of cross section before smearing.
        """
        
        self.parton_xsec_ = parton_xsec_
        self.xsec_nosmear_grid = None
    
    
    def scale(self, mtt):
        """Evaluate renormalization scale, given parton-level mtt."""
        
        return mtt / 2 * self.muR_scale_factor
    
    
    def selection_efficiency(self, mtt, subprocess):
        """Compute efficiency of event selection.
        
        Arguments:
            mtt:  Parton-level mtt, in GeV.
            subprocess:  String 'Res' or 'Int' to choose resonant part
                or interference.
        
        Return value:
            Efficiency.
        
        Computed efficiency includes accounts for the branching ratio of
        targeted decays, lepton identification, and b-tagging.
        """
        
        # Fitted parameters of a log-cubic approximation
        if subprocess == 'Res':
            coeffs = [-0.1166857, 2.19917117, -13.58990087, 27.78332692]
        elif subprocess == 'Int':
            coeffs = [-0.05877867, 1.03660773, -5.91517033, 11.11336388]
        else:
            raise RuntimeError('Do not recognize subprocess "{}".'.format(subprocess))
        
        return np.polyval(coeffs, math.log(mtt)) * self.target_branching * self.add_sel_eff
    
    
    def xsec(self, mtt, num_sigma=3):
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
            num_sigma:  Defines truncation for Gaussian kernel.
        
        Return value:
            Full differential cross section in mtt, in pb / GeV.
        """
        
        # Precompute cross section without smearing if not done yet
        if self.xsec_nosmear_grid is None:
            self.build_xsec_grid()
        
        
        # Find the integration window
        half_width = num_sigma * self.resolution * mtt
        start = np.searchsorted(self.xsec_nosmear_grid[0], mtt - half_width, 'left')
        end = np.searchsorted(self.xsec_nosmear_grid[0], mtt + half_width, 'right')
        
        
        # Perform a convolution of the cross section without smearing
        # with a Gaussian kernel.  Will use different algorithms
        # depending on how many precomputed points the integration
        # window contains.
        if (end - start) / num_sigma > 20:
            # Compute weights for the Gaussian kernel
            weights = np.empty(end - start)
            
            for i in range(end - start):
                cur_mtt = self.xsec_nosmear_grid[0, start + i]
                sigma = self.resolution * cur_mtt
                weights[i] = 1 / sigma * math.exp(-0.5 * ((mtt - cur_mtt) / sigma) ** 2)
            
            weights /= math.sqrt(2 * math.pi)
            
            # Correct for truncation in the kernel.  This is not exact
            # because the variance of the kernel is not constant.
            weights /= math.erf(num_sigma / math.sqrt(2))
            
            return simps(
                self.xsec_nosmear_grid[1, start:end] * weights,
                self.xsec_nosmear_grid[0, start:end], 'first'
            )
        
        else:
            # The grid of precomputed points is too coarse for the
            # assumed resolution.  The cross section without smearing
            # can be approaximated using linear interpolation as it does
            # not change rapidly within the window, but more poitns are
            # needed for the convolution with the Gaussian kernel.
            
            # Construct the interpolation object if not done yet
            if self.xsec_nosmear_interp is None:
                self.xsec_nosmear_interp = interp1d(
                    self.xsec_nosmear_grid[0], self.xsec_nosmear_grid[1],
                    copy=False, assume_sorted=True, bounds_error=False, fill_value=0.
                )
            
            x = np.linspace(mtt - half_width, mtt + half_width, num=101)
            y = np.empty_like(x)
            
            for i in range(len(x)):
                sigma = self.resolution * x[i]
                weight = 1 / sigma * math.exp(-0.5 * ((mtt - x[i]) / sigma) ** 2)
                y[i] = self.xsec_nosmear_interp(x[i]) * weight
            
            y /= math.sqrt(2 * math.pi) * math.erf(num_sigma / math.sqrt(2))
            
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
        
        shat = mtt ** 2
        alpha_s = self.pdf.alphasQ(self.scale(mtt))
        
        # Evaluate parts that are different for the two subprocesses
        res = self.parton_xsec.xsec_res(mtt, alpha_s) * self.selection_efficiency(mtt, 'Res')
        int_ = self.parton_xsec.xsec_int(mtt, alpha_s) * self.selection_efficiency(mtt, 'Int')
        
        # Convolute with PDF.  The last two terms appear from
        # translation of d[sigma] / d[sqrt(shat)] to d[sigma] / d[shat].
        xsec = (res + int_) * self.shat_pdf_interp(shat) * 2 * mtt
        
        return xsec
