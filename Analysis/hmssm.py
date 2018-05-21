import math

import ROOT


class PartonXSec:
    """A class to compute cross sections for gg -> H/A -> tt in hMSSM.
    
    Provides methods to compute parton-level cross sections for the
    resonant production of CP-even or odd states as well as interference
    with SM gg -> tt background.  Also provides a method to compute the
    total cross section including all these components or some of them.
    The cross sections depend on the s Mandelstam variable and the
    strong coupling constant.  They are computed at the leading order
    and rescaled with naive k-factors.  Expressions for the cross
    sections have been adapted from [1], but fixed widths are used.
    
    There are class attributes: mass of the top quark and Fermi
    constant.
    
    [1] Dicus et al., http://arxiv.org/abs/hep-ph/9404359
    """
    
    mt = 173.  # GeV
    gF = 1.166390e-05  # GeV^(-2)
    
    
    def __init__(self, mA, tanBeta, paramFile, components={'ARes', 'AInt', 'HRes', 'HInt'}):
        """Initialize from mA and tan(beta).
        
        Arguments:
            mA, tanBeta:  Core parameters of hMSSM: mass of the CP-odd
                Higgs boson, in GeV, and tan(beta).
            paramFile:  Path to file from which dependent parameters of
                the model will be extracted.  Must follow same format as
                in [1].
            components:  Default set of components of the cross section
                to be included.
        
        [1] https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWGMSSMNeutral#ROOT_histograms_MSSM_benchmark_s
        """
        
        self.mA = mA
        self.tanBeta = tanBeta
        
        
        # Set dependent parameters of the theory
        paramFile = ROOT.TFile(paramFile)
        
        self.mH = paramFile.Get('m_H').Interpolate(mA, tanBeta)
        self.wA = paramFile.Get('width_A').Interpolate(mA, tanBeta)
        self.wH = paramFile.Get('width_H').Interpolate(mA, tanBeta)
        self.gA = 1 / tanBeta
        self.gH = paramFile.Get('rescale_gt_H').Interpolate(mA, tanBeta)
        
        paramFile.Close()
        
        
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
        
        
        # Typical scale, in GeV, on which cross section can change
        # substantially
        self.variationScale = min(self.wA, self.wH)
    
    
    @staticmethod
    def to_pb(xsec):
        """Convert cross section from GeV^(-2) to pb."""
        
        return xsec / 2.56819e-9
    
    
    def xsec(self, s, alphaS, components=None):
        """Compute combined cross section for gg -> H/A -> tt.
        
        Arguments:
            s:  Value of parton-level s variable, in GeV^2, at which to
                evaluate the cross section.
            alphaS:  Value of the strong coupling constant.
            components:  Set of components of the cross section to be
                included.  If None (default), use the set given at the
                initialization.
        
        Return value:
            Computed cross section, in pb.
        """
        
        if components is None:
            components = self.components
        else:
            self._check_components(components)
        
        xsec = 0.
        
        if 'ARes' in components:
            xsec += self.xsec_odd_res(s, alphaS)
        
        if 'AInt' in components:
            xsec += self.xsec_odd_int(s, alphaS)
        
        if 'HRes' in components:
            xsec += self.xsec_even_res(s, alphaS)
        
        if 'HInt' in components:
            xsec += self.xsec_even_int(s, alphaS)
        
        return xsec
    
    
    def xsec_res(self, s, alphaS):
        """Compute cross section for resonant gg -> H/A -> tt."""
        
        return self.xsec_even_res(s, alphaS) + self.xsec_odd_res(s, alphaS)
    
    
    def xsec_int(self, s, alphaS):
        """Compute cross section for interference in gg -> H/A -> tt."""
        
        return self.xsec_even_int(s, alphaS) + self.xsec_odd_int(s, alphaS)
    
    
    def xsec_even_res(self, s, alphaS):
        """Compute cross section for resonant gg -> H -> tt."""
        
        if s < 4 * self.mt ** 2:
            return 0.
        
        beta = math.sqrt(1 - 4 * self.mt ** 2 / s)
        y = math.log((1 + beta) / (1 - beta))
        
        a = 3 * (alphaS * self.gF * self.mt ** 3) ** 2 * beta ** 3 / (1024 * math.pi ** 3)
        b = 16 + 8 * beta ** 2 * (math.pi ** 2 - y ** 2) + beta ** 4 * (math.pi ** 2 + y ** 2) ** 2
        denom = (s - self.mH ** 2) ** 2 + (self.mH * self.wH) ** 2
        
        return self.to_pb(self.kHRes * self.gH ** 4 * a * b / denom)
    
    
    def xsec_even_int(self, s, alphaS):
        """Compute cross section for interference in gg -> H -> tt."""
        
        if s < 4 * self.mt ** 2:
            return 0.
        
        beta = math.sqrt(1 - 4 * self.mt ** 2 / s)
        y = math.log((1 + beta) / (1 - beta))
        
        a = -alphaS ** 2 * self.gF * self.mt ** 4 * beta ** 2 / \
            (32 * math.pi * math.sqrt(2) * s) * y
        b = (s - self.mH ** 2) * (4 + beta ** 2 * (math.pi ** 2 - y ** 2)) + \
            2 * math.pi * beta ** 2 * self.mH * self.wH * y
        denom = (s - self.mH ** 2) ** 2 + (self.mH * self.wH) ** 2
        
        return self.to_pb(self.kHInt * self.gH ** 2 * a * b / denom)
    
    
    def xsec_odd_res(self, s, alphaS):
        """Compute cross section for resonant gg -> A -> tt."""
        
        if s < 4 * self.mt ** 2:
            return 0.
        
        beta = math.sqrt(1 - 4 * self.mt ** 2 / s)
        y = math.log((1 + beta) / (1 - beta))
        
        a = 3 * (alphaS * self.gF * self.mt ** 3) ** 2 * beta / (1024 * math.pi ** 3)
        b = (math.pi ** 2 + y ** 2) ** 2
        denom = (s - self.mA ** 2) ** 2 + (self.mA * self.wA) ** 2
        
        return self.to_pb(self.kARes * self.gA ** 4 * a * b / denom)
    
    
    def xsec_odd_int(self, s, alphaS):
        """Compute cross section for interference in gg -> A -> tt."""
        
        if s < 4 * self.mt ** 2:
            return 0.
        
        beta = math.sqrt(1 - 4 * self.mt ** 2 / s)
        y = math.log((1 + beta) / (1 - beta))
        
        a = -alphaS ** 2 * self.gF * self.mt ** 4 / (32 * math.pi * math.sqrt(2) * s) * y
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
