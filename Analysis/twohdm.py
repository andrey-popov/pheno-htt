import math

from spectrum import PartonXSec


class XSecTwoHDM(PartonXSec):
    """A class to compute cross sections for gg -> S -> tt in 2HDM.
    
    Expressions for the cross sections have been adapted from [1], with
    a difference of using fixed widths.
    
    [1] Dicus et al., http://arxiv.org/abs/hep-ph/9404359
    """
    
    def __init__(self, mA, wA, gA, mH, wH, gH):
        """Initialize from properties of neutral Higgs bosons.
        
        Arguments:
            mA, mH:  Masses of CP-odd and even states, in GeV.
            wA, wH:  Total widths of CP-odd and even states, in GeV.
            gA, gH:  Couplings of CP-odd and even states to top quarks.
                Expressed as scale factors to be applied to standard
                Yukawa couplings.
        """
        
        self.mA = mA
        self.wA = wA
        self.gA = gA
        
        self.mH = mH
        self.wH = wH
        self.gH = gH
        
        
        # Manually set naive scale factors for higher-order corrections.
        # The NLO k-factors for gg -> A/H are around 2, and the NNLO
        # k-factor for SM tt is around 1.6.  Use the geometric mean of
        # the two for the interference as suggested in [1].
        # [1] Hespel et al., https://arxiv.org/abs/1606.04149
        self.kA_res = self.kH_res = 2.
        k_bkg = 1.6
        self.kA_int = self.kH_int = math.sqrt(self.kA_res * k_bkg)
        
        # Set scale over which the cross section changes
        self.var_scale = min(self.wA, self.wH)
    
    
    def xsec_res(self, sqrt_s, alpha_s):
        """Compute cross section for resonant part in gg -> S -> tt.
        
        Arguments:
            sqrt_s:  Square root of Mandelstam s variable, in GeV.
            alpha_s:  Value of the strong coupling constant.
        
        Return value:
            Computed cross section, in pb.
        
        Implements abstract method of the superclass.
        """
        
        return self.xsec_even_res(sqrt_s, alpha_s) + self.xsec_odd_res(sqrt_s, alpha_s)
    
    
    def xsec_int(self, sqrt_s, alpha_s):
        """Compute cross section for interference in gg -> S -> tt.
        
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
        
        beta = self.beta(s)
        y = math.log((1 + beta) / (1 - beta))
        
        a = 3 * (alpha_s * self.gF * self.mt ** 3) ** 2 * beta ** 3 / (1024 * math.pi ** 3)
        b = 16 + 8 * beta ** 2 * (math.pi ** 2 - y ** 2) + beta ** 4 * (math.pi ** 2 + y ** 2) ** 2
        denom = (s - self.mH ** 2) ** 2 + (self.mH * self.wH) ** 2
        
        return self.to_pb(self.kH_res * self.gH ** 4 * a * b / denom)
    
    
    def xsec_even_int(self, sqrt_s, alpha_s):
        """Compute cross section for interference in gg -> H -> tt."""
        
        s = sqrt_s ** 2
        
        if s < 4 * self.mt ** 2:
            return 0.
        
        beta = self.beta(s)
        y = math.log((1 + beta) / (1 - beta))
        
        a = -alpha_s ** 2 * self.gF * self.mt ** 4 * beta ** 2 / \
            (32 * math.pi * math.sqrt(2) * s) * y
        b = (s - self.mH ** 2) * (4 + beta ** 2 * (math.pi ** 2 - y ** 2)) + \
            2 * math.pi * beta ** 2 * self.mH * self.wH * y
        denom = (s - self.mH ** 2) ** 2 + (self.mH * self.wH) ** 2
        
        return self.to_pb(self.kH_int * self.gH ** 2 * a * b / denom)
    
    
    def xsec_odd_res(self, sqrt_s, alpha_s):
        """Compute cross section for resonant gg -> A -> tt."""
        
        s = sqrt_s ** 2
        
        if s < 4 * self.mt ** 2:
            return 0.
        
        beta = self.beta(s)
        y = math.log((1 + beta) / (1 - beta))
        
        a = 3 * (alpha_s * self.gF * self.mt ** 3) ** 2 * beta / (1024 * math.pi ** 3)
        b = (math.pi ** 2 + y ** 2) ** 2
        denom = (s - self.mA ** 2) ** 2 + (self.mA * self.wA) ** 2
        
        return self.to_pb(self.kA_res * self.gA ** 4 * a * b / denom)
    
    
    def xsec_odd_int(self, sqrt_s, alpha_s):
        """Compute cross section for interference in gg -> A -> tt."""
        
        s = sqrt_s ** 2
        
        if s < 4 * self.mt ** 2:
            return 0.
        
        beta = self.beta(s)
        y = math.log((1 + beta) / (1 - beta))
        
        a = -alpha_s ** 2 * self.gF * self.mt ** 4 / (32 * math.pi * math.sqrt(2) * s) * y
        b = (s - self.mA ** 2) * (math.pi ** 2 - y ** 2) + 2 * math.pi * self.mA * self.wA * y
        denom = (s - self.mA ** 2) ** 2 + (self.mA * self.wA) ** 2
        
        return self.to_pb(self.kA_int * self.gA ** 2 * a * b / denom)
