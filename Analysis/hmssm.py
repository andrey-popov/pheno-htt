import math

import ROOT

from spectrum import PartonXSec


class XSecHMSSM(PartonXSec):
    """A class to compute cross sections for gg -> H -> tt in hMSSM.
    
    Expressions for the cross sections have been adapted from [1], with
    a difference of using fixed widths.
    
    [1] Dicus et al., http://arxiv.org/abs/hep-ph/9404359
    """
    
    def __init__(self, mA, tanbeta, paramfile, components={'ARes', 'AInt', 'HRes', 'HInt'}):
        """Initialize from mA and tan(beta).
        
        Arguments:
            mA, tanbeta:  Core parameters of hMSSM: mass of the CP-odd
                Higgs boson, in GeV, and tan(beta).
            paramfile:  Path to file from which dependent parameters of
                the model will be extracted.  Must follow same format as
                in [1].
            components:  Default set of components of the cross section
                to be included.
        
        [1] https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWGMSSMNeutral#ROOT_histograms_MSSM_benchmark_s
        """
        
        self.mA = mA
        self.tanbeta = tanbeta
        
        
        # Set dependent parameters of the theory
        paramfile = ROOT.TFile(paramfile)
        
        self.mH = paramfile.Get('m_H').Interpolate(mA, tanbeta)
        self.wA = paramfile.Get('width_A').Interpolate(mA, tanbeta)
        self.wH = paramfile.Get('width_H').Interpolate(mA, tanbeta)
        self.gA = 1 / tanbeta
        self.gH = paramfile.Get('rescale_gt_H').Interpolate(mA, tanbeta)
        
        paramfile.Close()
        
        
        # Manually set naive scale factors for higher-order corrections.
        # The NLO k-factors for gg -> A/H are around 2, and the NNLO
        # k-factor for SM tt is around 1.6.  Use the geometric mean of
        # the two for the interference as suggested in [1].
        # [1] Hespel et al., https://arxiv.org/abs/1606.04149
        self.kARes = self.kHRes = 2.
        self.kAInt = self.kHInt = math.sqrt(2. * 1.6)
        
        
        # Store default set of components to evaluate
        self._check_components(components)
        self.components = components
        
        
        # Set scale over which the cross section changes
        self.var_scale = min(self.wA, self.wH)
    
    
    def xsec(self, sqrt_s, alpha_s, components=None):
        """Compute cross section for gg -> H -> tt.
        
        Arguments:
            sqrt_s:  Square root of Mandelstam s variable, in GeV.
            alpha_s:  Value of the strong coupling constant.
        
        Return value:
            Computed cross section, in pb.
        """
        
        if components is None:
            components = self.components
        else:
            self._check_components(components)
        
        xsec = 0.
        
        if 'ARes' in components:
            xsec += self.xsec_odd_res(sqrt_s, alpha_s)
        
        if 'AInt' in components:
            xsec += self.xsec_odd_int(sqrt_s, alpha_s)
        
        if 'HRes' in components:
            xsec += self.xsec_even_res(sqrt_s, alpha_s)
        
        if 'HInt' in components:
            xsec += self.xsec_even_int(sqrt_s, alpha_s)
        
        return xsec
    
    
    def xsec_res(self, sqrt_s, alpha_s):
        """Compute cross section for resonant part in gg -> H -> tt.
        
        Arguments:
            sqrt_s:  Square root of Mandelstam s variable, in GeV.
            alpha_s:  Value of the strong coupling constant.
        
        Return value:
            Computed cross section, in pb.
        
        Implements abstract method of the superclass.
        """
        
        return self.xsec_even_res(sqrt_s, alpha_s) + self.xsec_odd_res(sqrt_s, alpha_s)
    
    
    def xsec_int(self, sqrt_s, alpha_s):
        """Compute cross section for interference in gg -> H -> tt.
        
        Arguments:
            sqrt_s:  Square root of Mandelstam s variable, in GeV.
            alpha_s:  Value of the strong coupling constant.
        
        Return value:
            Computed cross section, in pb.
        
        Implements abstract method of the superclass.
        """
        
        return self.xsec_even_int(sqrt_s, alpha_s) + self.xsec_odd_int(sqrt_s, alpha_s)
    
    
    def xsec_even_res(self, sqrt_s, alpha_s):
        """Compute cross section for resonant gg -> H -> tt."""
        
        s = sqrt_s ** 2
        
        if s < 4 * self.mt ** 2:
            return 0.
        
        beta = math.sqrt(1 - 4 * self.mt ** 2 / s)
        y = math.log((1 + beta) / (1 - beta))
        
        a = 3 * (alpha_s * self.gF * self.mt ** 3) ** 2 * beta ** 3 / (1024 * math.pi ** 3)
        b = 16 + 8 * beta ** 2 * (math.pi ** 2 - y ** 2) + beta ** 4 * (math.pi ** 2 + y ** 2) ** 2
        denom = (s - self.mH ** 2) ** 2 + (self.mH * self.wH) ** 2
        
        return self.to_pb(self.kHRes * self.gH ** 4 * a * b / denom)
    
    
    def xsec_even_int(self, sqrt_s, alpha_s):
        """Compute cross section for interference in gg -> H -> tt."""
        
        s = sqrt_s ** 2
        
        if s < 4 * self.mt ** 2:
            return 0.
        
        beta = math.sqrt(1 - 4 * self.mt ** 2 / s)
        y = math.log((1 + beta) / (1 - beta))
        
        a = -alpha_s ** 2 * self.gF * self.mt ** 4 * beta ** 2 / \
            (32 * math.pi * math.sqrt(2) * s) * y
        b = (s - self.mH ** 2) * (4 + beta ** 2 * (math.pi ** 2 - y ** 2)) + \
            2 * math.pi * beta ** 2 * self.mH * self.wH * y
        denom = (s - self.mH ** 2) ** 2 + (self.mH * self.wH) ** 2
        
        return self.to_pb(self.kHInt * self.gH ** 2 * a * b / denom)
    
    
    def xsec_odd_res(self, sqrt_s, alpha_s):
        """Compute cross section for resonant gg -> A -> tt."""
        
        s = sqrt_s ** 2
        
        if s < 4 * self.mt ** 2:
            return 0.
        
        beta = math.sqrt(1 - 4 * self.mt ** 2 / s)
        y = math.log((1 + beta) / (1 - beta))
        
        a = 3 * (alpha_s * self.gF * self.mt ** 3) ** 2 * beta / (1024 * math.pi ** 3)
        b = (math.pi ** 2 + y ** 2) ** 2
        denom = (s - self.mA ** 2) ** 2 + (self.mA * self.wA) ** 2
        
        return self.to_pb(self.kARes * self.gA ** 4 * a * b / denom)
    
    
    def xsec_odd_int(self, sqrt_s, alpha_s):
        """Compute cross section for interference in gg -> A -> tt."""
        
        s = sqrt_s ** 2
        
        if s < 4 * self.mt ** 2:
            return 0.
        
        beta = math.sqrt(1 - 4 * self.mt ** 2 / s)
        y = math.log((1 + beta) / (1 - beta))
        
        a = -alpha_s ** 2 * self.gF * self.mt ** 4 / (32 * math.pi * math.sqrt(2) * s) * y
        b = (s - self.mA ** 2) * (math.pi ** 2 - y ** 2) + 2 * math.pi * self.mA * self.wA * y
        denom = (s - self.mA ** 2) ** 2 + (self.mA * self.wA) ** 2
        
        return self.to_pb(self.kAInt * self.gA ** 2 * a * b / denom)
    
    
    def _check_components(self, components):
        """Check if all provided component labels are recognized.
        
        Arguments:
            components:  Iterable with labels of components of the cross
                section to be checked.
        
        Return value:
            True if all labels are recognized, raise an exception
            otherwise.
        """
        
        allowedComponents = {'ARes', 'AInt', 'HRes', 'HInt'}
        
        for component in components:
            if component not in allowedComponents:
                raise RuntimeError('Do not recognize component "{}".'.format(component))
        
        return True
