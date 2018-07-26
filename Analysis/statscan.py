"""
Exports classes to scan over parameters of a model.
"""

import contextlib
import os
import sys
from uuid import uuid4

import numpy as np
from scipy.integrate import trapz
from scipy.interpolate import RectBivariateSpline

import matplotlib as mpl
mpl.use('agg')
from matplotlib import pyplot as plt

import ROOT
from ROOT.RooStats import HistFactory


mpl.rc('xtick', top=True, direction='in')
mpl.rc('ytick', right=True, direction='in')
mpl.rc('axes', labelsize='large')
mpl.rc('axes.formatter', limits=[-2, 4], use_mathtext=True)

ROOT.Math.MinimizerOptions.SetDefaultMinimizer('Minuit2', 'Minimize')
ROOT.RooStats.UseNLLOffset(True)


@contextlib.contextmanager
def suppress_stdout():
    """Context manager to suppress stdout.
    
    Works for Python and external libraries alike.  Code from [1].
    [1] https://stackoverflow.com/a/17954769/966461
    """
    
    fd = sys.stdout.fileno()
    
    def _redirect_stdout(to):
        sys.stdout.close()
        os.dup2(to.fileno(), fd)
        sys.stdout = os.fdopen(fd, 'w')
    
    with os.fdopen(os.dup(fd), 'w') as old_stdout:
        with open(os.devnull, 'w') as f:
            _redirect_stdout(to=f)
        try:
            yield
        finally:
            _redirect_stdout(to=old_stdout)


class Grid:
    """Class to facilitate scanning over a 2D grid.
    
    For each node of the grid, the significance and the CLs value can be
    attached.  They are stored in 2D arrays with coordinates in the
    order (x, y).  The grid and values of significance and CLs can be
    saved in an .npz file and read back from it.
    """
    
    def __init__(self, x, y):
        """Initialize grid from coordinates along each axis.
        
        Arguments:
            x, y:  1D array_like of coordinates in ascending order.
        """
        
        self.x = np.asarray(x)
        self.y = np.asarray(y)
        
        self.significance = np.zeros((len(x), len(y)))
        self.cls = np.zeros_like(self.significance)
    
    
    @classmethod
    def fromfile(cls, filename):
        """Construct from an .npz file.
        
        Arguments:
            filename:  Path to an .npz file created by method Grid.save.
        
        Return value:
            Instance of Grid constructed from the file.
        """
        
        data = np.load(filename)
        obj = cls(data['x'], data['y'])
        obj.significance = data['significance']
        obj.cls = data['cls']
        
        return obj
    
    
    def __iter__(self):
        """Iterate over the grid.
        
        Yields global index of the current node and its x and y
        coordinates.
        """
        
        for iy in range(len(self.y)):
            for ix in range(len(self.x)):
                global_index = ix + iy * len(self.x)
                yield (global_index, self.x[ix], self.y[iy])
    
    
    def save(self, filename):
        """Save to an .npz file.
        
        Write the grid and attached values of significance and CLs.
        
        Arguments:
            filename:  Name of output file.
        
        Return value:
            None.
        """
        
        np.savez(filename, x=self.x, y=self.y, significance=self.significance, cls=self.cls)
    
    
    def set(self, global_index, significance, cls):
        """Set significance and CLs for a node.
        
        The node is identified by a global index as returned by
        __iter__.
        
        Arguments:
            global_index:  Global index identifying the node.
            significance, cls:  Values of significance and CLs to be
                attached.
        
        Return value:
            None.
        """
        
        ix = global_index % len(self.x)
        iy = global_index // len(self.x)
        
        self.significance[ix, iy] = significance
        self.cls[ix, iy] = cls


class PlotScan:
    """Class to plot significance and CLs obtained in a scan.
    
    Produces a figure with significance of the signal, represented with
    a colour map and contours, and 95% CL exclusion shown with a
    contour.  Provided values of significance and CLs on a grid are
    up-sampled for visual aesthetics.  Resulting figure is returned to
    the user for customization and saving to a file.
    """
    
    def __init__(self, grid, max_significance=10.):
        """Initialize from Grid object.
        
        Arguments:
            grid:  Instance of Grid.
            max_significance:  Maximal significance.  Larger values are
                clipped.
        """
        
        self.x = grid.x
        self.y = grid.y
        self.significance = grid.significance
        self.cls = grid.cls
        
        self.max_significance = max_significance
        
        # Sometimes fits fail when the significance is very large.
        # Replace the results with placeholders in this case.  Also
        # clip valid but large values of significance, which is needed
        # to prevent instability in interpolation.
        for val in np.nditer(self.significance, op_flags=['readwrite']):
            if not np.isfinite(val) or val > self.max_significance:
                val[...] = self.max_significance
        
        for val in np.nditer(self.cls, op_flags=['readwrite']):
            if not np.isfinite(val):
                val[...] = 0.
        
        self.fig = None
        self.axes = None
    
    
    def draw(self):
        """Construct the figure.
        
        Arguments:
            None.
        
        Return value:
            Constructed figure and axes.  They are also added to self as
            attributes.
        """
        
        # Up-sample the arrays to produce smooth colour map and
        # contours
        x_up, y_up, significance_up, cls_up = self._upsample(10)
        
        
        self.fig = plt.figure()
        self.fig.patch.set_alpha(0.)  # Transparent background
        self.axes = self.fig.add_subplot(111)
        
        
        # Plot significance with a colour map.  Coordinates of grid
        # nodes are shifted since pcolormesh expects coordinates of bin
        # edges rather than centres.
        xx_c, yy_c = np.meshgrid(self._centre_bins(x_up), self._centre_bins(y_up))
        colormap = plt.get_cmap('viridis')
        
        image = self.axes.pcolormesh(
            xx_c, yy_c, significance_up.T,
            cmap=colormap, vmin=0., vmax=self.max_significance
        )
        self.fig.colorbar(
            image, fraction=0.05, pad=0.02,
            label='Significance [$\\sigma$]', spacing='proportional'
        )
        
        
        # Draw contours
        xx, yy = np.meshgrid(x_up, y_up)
        
        contours = self.axes.contour(
            xx, yy, significance_up.T,
            [1., 3., 5.], colors='white', zorder=1.5
        )
        self.axes.clabel(contours, fmt='%g $\\sigma$')
        
        contour_cls = self.axes.contour(
            xx, yy, cls_up.T,
            [0.05], colors='red', zorder=1.5
        )
        self.axes.clabel(contour_cls, fmt='95%% CL excl.')
        
        return self.fig, self.axes
    
    
    @staticmethod
    def _centre_bins(a, clip=True):
        """Create new binning whose bins are centered at given points.
        
        Arguments:
            a:  1D array_like that defines positions of bin centres.
            clip:  Determines whether first and last bins extend beyond
                the full range of given array or get clipped.
        
        Return value:
            Array defining the new binning.  Its length is larger than
            that of input array by one.
        """
        
        anew = np.empty(len(a) + 1, dtype=a.dtype)
        anew[1:-1] = (a[:-1] + a[1:]) / 2
        
        if clip:
            anew[0] = a[0]
            anew[-1] = a[-1]
        else:
            anew[0] = a[0] - (a[1] - a[0]) / 2
            anew[-1] = a[-1] + (a[-1] - a[-2]) / 2
        
        return anew
    
    
    def _upsample(self, up_factor=10, degree=3):
        """Up-sample the grid.
        
        Produce a denser grid using spline interpolation.
        
        Arguments:
            up_factor:  Up-sampling factor.
            degree:  Degree for the splines.
        
        Return value:
            Up-sampled x and y coordinates, significance, and CLs
            values.
        """
        
        # Uniformally split each bin in the grid into up_factor segments
        x_up = np.empty((len(self.x) - 1) * up_factor + 1)
        y_up = np.empty((len(self.y) - 1) * up_factor + 1)
        
        for src, up in [(self.x, x_up), (self.y, y_up)]:
            for i in range(len(src) - 1):
                up[i * up_factor : (i + 1) * up_factor] = np.linspace(
                    src[i], src[i + 1], num=up_factor, endpoint=False
                )
            
            up[-1] = src[-1]
        
        
        # Up-sample significance and CLs values using spline
        # interpolation
        significance_interp = RectBivariateSpline(
            self.x, self.y, self.significance,
            kx=min(degree, len(self.x) - 1), ky=min(degree, len(self.y) - 1)
        )
        significance_up = significance_interp(x_up, y_up)
        
        cls_interp = RectBivariateSpline(
            self.x, self.y, self.cls,
            kx=min(degree, len(self.x) - 1), ky=min(degree, len(self.y) - 1)
        )
        cls_up = cls_interp(x_up, y_up)
        
        return (x_up, y_up, significance_up, cls_up)


class StatCalc:
    """Class to perform statistical evaluation for H -> tt.
    
    Distributions of mtt in signal are computed with the help of an
    instance of spectrum.RecoMtt, while distributions for the SM tt
    background are read from a file.  Uncertainty in the renormalization
    scale is included for the signal, and a number of uncertainties are
    considered for the background.  This class allows to compute the
    expected significance if the signal is present and the expected CLs
    value if it does not exist.
    """
    
    def __init__(self, signal_distr, bkgfile, lumi):
        """Initialize from signal and background distributions.
        
        Arguments:
            signal_distr:  Object to compute event density in mtt
                distribution for signal, in pb/GeV.  An instance of
                spectrum.RecoMtt is expected.
            bkgfile:  Name of ROOT file with distributions for the
                background, in pb.
            lumi:  Target integrated luminosity, 1/pb.
        
        Binning of background histograms is used for signal
        distributions as well.
        """
        
        self.signal_distr_calc = signal_distr
        self.lumi = lumi
        
        self._signal_templates = None
        self._workspace = None
        
        
        bkgfile = ROOT.TFile(bkgfile)
        
        # Read binning from the file with backgrounds
        hist = bkgfile.Get('TT')
        self.binning = np.asarray([
            hist.GetBinLowEdge(bin) for bin in range(1, hist.GetNbinsX() + 2)
        ])
        
        
        # Read all background templates converting them to trivial
        # uniform binning, as needed for HistFactory [1].
        # [1] https://root-forum.cern.ch/t/binned-ml-fit-with-histfactory/28867
        self._bkg_templates = {}
        
        for key in bkgfile.GetListOfKeys():
            self._bkg_templates[key.GetName()] = self._strip_binning(key.ReadObj())
        
        bkgfile.Close()


    def cls(self):
        """Compute CLs value on b-only Asimov data set.
        
        Return value:
            CLs value.
        
        Build RooStats models if they are not available.
        """
        
        if not self._workspace:
            self._build_model()
        
        
        # Asymptotic calculator for Asimov b-only data set.  For an
        # exclusion the null hypothesis is that r >= 1 and the
        # alternative hypothesis is r = 0.
        calc = ROOT.RooStats.AsymptoticCalculator(
            self._workspace.data('asimovData_bOnly'), self._bOnlyModel, self._fullModel
        )
        
        # Configure the calculator to use the q_mu tilde statistics from
        # [1], which is intended for upper limits with additional
        # condition r > 0.
        # [1] https://arxiv.org/abs/1007.1727
        calc.SetOneSided(True)
        calc.SetQTilde(True)
        
        
        # Perform the test.  The Asimov data set provided to the
        # constructor of calc is fitted.  Because of the data set used,
        # the "observed" CLs coincides with median expected one.
        testResult = calc.GetHypoTest()
        testResult.SetBackgroundAsAlt(True)
        
        return testResult.CLs()
            

    def significance(self):
        """Compute significance on s+b Asimov data set.
        
        Quantifies the incompatibility between the s+b Asimov data set
        and post-fit SM expectation.
        
        Return value:
            Expected significance.
        
        Build RooStats models if they are not available.
        """
        
        if not self._workspace:
            self._build_model()
        
        
        # Asymptotic calculator for the significance with Asimov s+b
        # data set as the data
        calc = ROOT.RooStats.AsymptoticCalculator(
            self._workspace.data('asimovData'), self._fullModel, self._bOnlyModel
        )
        
        # Configure the calculator to use the q0 statistics from [1],
        # which is intended for discovery
        # [1] https://arxiv.org/abs/1007.1727
        calc.SetOneSidedDiscovery(True)
        
        
        # Perform the test.  The Asimov data set provided to the
        # constructor of calc is fitted.  Because of the data set used,
        # the "observed" significance coincides with median expected
        # one.
        testResult = calc.GetHypoTest()
        
        return testResult.Significance()
    
    
    def update_signal(self, signal_distr):
        """Update object to compute signal distributions.
        
        Arguments:
            signal_distr:  Object to compute event density in mtt
                distribution for signal, in pb/GeV.  An instance of
                spectrum.RecoMtt is expected.
        
        Return value:
            None.
        
        As a side effect, invalidate precomputed signal templates and
        the workspace.
        """
        
        self.signal_distr_calc = signal_distr
        self._signal_templates = None
        self._workspace = None


    def _build_model(self):
        """Build RooFit workspace with statistical model.
        
        Construct signal templates if they are not available.
        """
        
        if not self._signal_templates:
            self._build_signal_templates()
        
        
        measurement = HistFactory.Measurement('Measurement')
        
        measurement.AddPOI('r')
        measurement.AddPreprocessFunction('negR', '-r', 'r[1,0,10]')
        measurement.AddConstantParam('Lumi')
        measurement.SetLumi(self.lumi)
        
        
        # Define all samples and corresponding systematic variations.
        # Clone all histograms because HistFactory takes ownership of
        # them for no reason.
        channel = HistFactory.Channel('Channel')
        
        sgn_pos = HistFactory.Sample('SgnPos')
        sgn_pos.SetHisto(self._signal_templates['SgnPos'].Clone(uuid4().hex))
        sgn_pos.AddNormFactor('r', 1., 0., 10., True)
        sgn_pos.SetNormalizeByTheory(True)
        
        syst = HistFactory.HistoSys('RenormScale')
        syst.SetHistoHigh(self._signal_templates['SgnPos_RenormScaleUp'].Clone(uuid4().hex))
        syst.SetHistoLow(self._signal_templates['SgnPos_RenormScaleDown'].Clone(uuid4().hex))
        sgn_pos.AddHistoSys(syst)
        
        sgn_neg = HistFactory.Sample('SgnNeg')
        sgn_neg.SetHisto(self._signal_templates['SgnNeg'].Clone(uuid4().hex))
        sgn_neg.AddNormFactor('negR', -1., -10., 0., True)
        sgn_neg.SetNormalizeByTheory(True)
        
        syst = HistFactory.HistoSys('RenormScale')
        syst.SetHistoHigh(self._signal_templates['SgnNeg_RenormScaleUp'].Clone(uuid4().hex))
        syst.SetHistoLow(self._signal_templates['SgnNeg_RenormScaleDown'].Clone(uuid4().hex))
        sgn_neg.AddHistoSys(syst)
        
        bkg = HistFactory.Sample('TT')
        bkg.SetHisto(self._bkg_templates['TT'].Clone(uuid4().hex))
        bkg.AddOverallSys('TTRate', 0.9, 1.1)
        
        systNames = ['MttScale', 'RenormScale', 'FactorScale', 'FSR', 'MassT', 'PDFAlphaS']
        systNames += ['PDF{}'.format(i) for i in range(1, 31)]
        
        for systName in systNames:
            syst = HistFactory.HistoSys(systName)
            syst.SetHistoHigh(self._bkg_templates['TT_{}Up'.format(systName)].Clone(uuid4().hex))
            syst.SetHistoLow(self._bkg_templates['TT_{}Down'.format(systName)].Clone(uuid4().hex))
            bkg.AddHistoSys(syst)
        
        
        for sample in [sgn_pos, sgn_neg, bkg]:
            channel.AddSample(sample)
        
        measurement.AddChannel(channel)
        
        
        # An s+b Asimov data set will be created automatically.  Add
        # also a b-only data set.
        bOnlyData = HistFactory.Asimov('asimovData_bOnly')
        bOnlyData.SetParamValue('r', 0.)
        measurement.AddAsimovDataset(bOnlyData)
        
        
        # Construct the workspace as well as s+b and b-only models
        with suppress_stdout():
            # Suppress stdout here since the command below will print
            # out the whole workspace
            self._workspace = \
                HistFactory.HistoToWorkspaceFactoryFast.MakeCombinedModel(measurement)
        
        self._fullModel = self._workspace.obj('ModelConfig')
        
        # Create a snapshot that fixes the POI in the full model
        poi = self._fullModel.GetParametersOfInterest().first()
        poi.setVal(1.)
        self._fullModel.SetSnapshot(ROOT.RooArgSet(poi))
        
        # Create the b-only model by fixing POI to zero in the full
        # model
        self._bOnlyModel = self._fullModel.Clone(self._fullModel.GetName() + '_bOnly')
        poi.setVal(0.)
        self._bOnlyModel.SetSnapshot(ROOT.RooArgSet(poi))


    def _build_signal_templates(self):
        """Construct templates for signal.
        
        Templates are split into contributions with positive and
        negative bin contents, and the sign of the latter ones is
        flipped.  As a result, they should be scaled with negative
        signal strength in the statistical model.  Construct also
        templates representing a variation of the renormalization scale
        by factor 2.  Produced histograms have trivial equidistant
        binning.
        """
        
        self._signal_templates = {}
        
        for muR_scale_factor, postfix in [
            (1., ''), (2., '_RenormScaleUp'), (0.5, '_RenormScaleDown')
        ]:
            self.signal_distr_calc.muR_scale_factor = muR_scale_factor
            
            
            # In each bin given by the binning, integrate the
            # differential cross section
            xsec_integral = np.empty(len(self.binning) - 1)
            num_samples = 10
            
            for ibin in range(len(xsec_integral)):
                
                x = np.linspace(self.binning[ibin], self.binning[ibin + 1], num=num_samples)
                y = np.empty_like(x)
                
                for isample in range(num_samples):
                    y[isample] = self.signal_distr_calc.xsec(x[isample])
                
                xsec_integral[ibin] = trapz(y, x)
            
            
            # Put results into histograms, separately for positive and
            # negative counts
            hist_pos = ROOT.TH1D(
                'SgnPos' + postfix, '',
                len(self.binning) - 1, 0., len(self.binning) - 1
            )
            hist_neg = hist_pos.Clone('SgnNeg' + postfix)
            
            for ibin in range(len(xsec_integral)):
                value = xsec_integral[ibin]
                
                if value >= 0.:
                    hist_pos.SetBinContent(ibin + 1, value)
                else:
                    hist_neg.SetBinContent(ibin + 1, value)
                
                hist_pos.SetBinError(ibin + 1, 0.)
                hist_neg.SetBinError(ibin + 1, 0.)
            
            hist_neg.Scale(-1.)
            
            for hist in [hist_pos, hist_neg]:
                hist.SetDirectory(None)
                self._signal_templates[hist.GetName()] = hist
    
    
    @staticmethod
    def _strip_binning(hist):
        """Copy ROOT histogram stripping binning information.
        
        Create a new histogram that has the same bin contents and
        errors, but uses a trivial equidistant binning.  The new
        histogram is not attached to any directory.
        
        Arguments:
            hist:  Source ROOT histogram.
        
        Return value:
            Newly created ROOT histogram with trivial equidistant
            binning.
        
        Ignore under- and overflow bins.
        """
        
        num_bins = hist.GetNbinsX()
        
        newhist = ROOT.TH1D(uuid4().hex, '', num_bins, 0., num_bins)
        newhist.SetDirectory(None)
        newhist.SetName(hist.GetName())
        
        for bin in range(1, num_bins + 1):
            newhist.SetBinContent(bin, hist.GetBinContent(bin))
            newhist.SetBinError(bin, hist.GetBinError(bin))
        
        return newhist
