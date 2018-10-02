import math

import ROOT

from twohdm import XSecTwoHDM


def mh_alpha(mA, tanbeta, mZ=91.15348, mh=125.0):
    """Compute mass of CP-even state and angle alpha in hMSSM.
    
    Arguments:
        mA, tanbeta:  Core parameters of hMSSM: mass of the CP-odd
            Higgs boson, in GeV, and tan(beta).
        mZ:  Mass of the Z boson, GeV.
        mh:  Mass of the lighter CP-even Higgs boson, GeV.
    
    Return value:
        Tuple consisting of computed mass of the heavier CP-even state,
        GeV, and angle alpha.
    
    Use Eq. 5 in [1].
    [1] Djouadi et al., https://arxiv.org/abs/1307.5205
    """
    
    tb2 = tanbeta ** 2
    sb2 = tb2 / (1 + tb2)
    cb2 = 1 / (1 + tb2)
    mH2 = (
        (mA ** 2 + mZ ** 2 - mh ** 2) * (mZ ** 2 * cb2 + mA ** 2 * sb2) -
        mA **2 * mZ ** 2 * ((1 - tb2) / (1 + tb2)) ** 2
    ) / (
        mZ ** 2 * cb2 + mA ** 2 * sb2 - mh ** 2
    )
    tanalpha = -(mZ ** 2 + mA ** 2) * math.sqrt(cb2 * sb2) / \
        (mZ ** 2 * cb2 + mA ** 2 * sb2 - mh ** 2)
    
    return math.sqrt(mH2), math.atan(tanalpha)



class XSecHMSSM(XSecTwoHDM):
    """A class to compute cross sections for gg -> H -> tt in hMSSM.
    
    This is a specialization of the class implementing 2HDM.  For the
    purpose of testing, it is possible to compute cross sections with
    only one CP state and/or only resonant part or interference.
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
        
        paramfile = ROOT.TFile(paramfile)
        
        # Compute dependent parameters of the theory
        mH, alpha = mh_alpha(mA, tanbeta)
        
        wA = paramfile.Get('width_A').Interpolate(mA, tanbeta)
        wH = paramfile.Get('width_H').Interpolate(mA, tanbeta)
        
        gA = 1 / tanbeta
        gH = math.sin(alpha) / (tanbeta / math.sqrt(1 + tanbeta ** 2))
        
        
        # Set parameters of 2HDM
        super().__init__(mA, wA, gA, mH, wH, gH)
        self.tanbeta = tanbeta
        
        
        # Override dummy k-factors.  The k-factor for the SM tt
        # backround is set to 2 [1], and for the interference use the
        # geomentric mean of the k-factors for the background and the
        # resonant part
        # [1] https://github.com/andrey-popov/pheno-htt/issues/2
        # [2] Hespel et al., https://arxiv.org/abs/1606.04149
        k_bkg = 2.
        self.kA_res = paramfile.Get('kA_NNLO_13TeV').Interpolate(mA, tanbeta)
        self.kH_res = paramfile.Get('kH_NNLO_13TeV').Interpolate(mA, tanbeta)
        self.kA_int = math.sqrt(self.kA_res * k_bkg)
        self.kH_int = math.sqrt(self.kH_res * k_bkg)
        
        paramfile.Close()
        
        
        # Store default set of components to evaluate
        self._check_components(components)
        self.components = components
    
    
    def xsec(self, sqrt_s, alpha_s, components=None):
        """Compute cross section for gg -> H -> tt.
        
        Arguments:
            sqrt_s:  Square root of Mandelstam s variable, in GeV.
            alpha_s:  Value of the strong coupling constant.
            components:  Components to use in the computation.  By
                default, use the set provided at initialization.
        
        Return value:
            Computed cross section, in pb.
        
        Reimplements method of a superclass.
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
