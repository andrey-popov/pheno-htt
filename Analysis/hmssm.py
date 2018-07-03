import math

import ROOT

from twohdm import XSecTwoHDM


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
        
        # Compute dependent parameters of the theory
        paramfile = ROOT.TFile(paramfile)
        
        mH = paramfile.Get('m_H').Interpolate(mA, tanbeta)
        wA = paramfile.Get('width_A').Interpolate(mA, tanbeta)
        wH = paramfile.Get('width_H').Interpolate(mA, tanbeta)
        gA = 1 / tanbeta
        gH = paramfile.Get('rescale_gt_H').Interpolate(mA, tanbeta)
        
        paramfile.Close()
        
        
        # Set parameters of 2HDM
        super().__init__(mA, wA, gA, mH, wH, gH)
        self.tanbeta = tanbeta
        
        
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
