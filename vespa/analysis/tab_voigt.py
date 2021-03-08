# Python modules
import os
import io as StringIO
import base64
import datetime

# 3rd party modules
import wx
import numpy as np

# Our modules
import vespa.analysis.constants as constants
import vespa.analysis.block_fit_voigt as block_fit_voigt
import vespa.analysis.tab_base as tab_base
import vespa.analysis.prefs as prefs_module
import vespa.analysis.util_menu as util_menu
import vespa.analysis.util_import as util_import
import vespa.analysis.plot_panel_voigt as plot_panel_voigt  
import vespa.analysis.util_analysis_config as util_analysis_config
import vespa.analysis.dialog_dataset_browser as dialog_dataset_browser
import vespa.analysis.dynamic_list_voigt_metabolite as dynamic_list_voigt_metabolite
import vespa.analysis.figure_layouts as figure_layouts 
import vespa.analysis.auto_gui.voigt as voigt

# TODO - move these into COMMON
import vespa.simulation.dialog_mixed_metabolite_output as dialog_mixout
from vespa.simulation.constants import ThirdPartyExportTypes

import vespa.common.mrs_prior as mrs_prior
import vespa.common.dialog_experiment_browser as dialog_experiment_browser
import vespa.common.wx_gravy.util as wx_util
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.util.misc as util_misc
import vespa.common.util.generic_spectral as util_spectral

from vespa.analysis.constants import FitLineshapeModel

from matplotlib.backends.backend_pdf import PdfPages



# # _PYWAVELETS_EXPLANATION_URL is the URL at which we explain PyWavelets
# # availability (which is slightly complicated).
# _PYWAVELETS_EXPLANATION_URL = "https://vespa-mrs.github.io/vespa.io/development/applications/dev_analysis/technical/PyWavelets.html?highlight=PyWavelets"
#
# # _PYWAVELETS_UNAVAILABLE_MSG is what the user sees when PyWavelets isn't
# # installed and he tries to use the baseline wavelet filter.
# _PYWAVELETS_UNAVAILABLE_MSG = """The PyWavelets module doesn't appear \
# to be installed. Wavelet filtering will be Python based.
#
# For more information, see:
# %s
#
# Do you want to open that Web page now?
# """ % _PYWAVELETS_EXPLANATION_URL


def _create_default_prior():
    # This creates & returns a Prior object populated with some default
    # metabs. It's exists just so that something shows up when users first
    # open this tab.

    metabolites = { "n-acetylaspartate" : { "spins" : 3, 
                                            "ppm"   : 2.0082,   
                                            "area"  : 3.0,
                                            "phase" : 0.0,
                                          },
                    "creatine"          : { "spins" : 3, 
                                            "ppm"   : 3.0241,   
                                            "area"  : 3.0,
                                            "phase" : 0.0,
                                          },
                    "cr2"               : { "spins" : 2, 
                                            "ppm"   : 3.90241,   
                                            "area"  : 2.0,
                                            "phase" : 0.0,
                                          },
                    "choline"           : { "spins"  : 9, 
                                            "ppm"   : 3.21,   
                                            "area"  : 9.0,
                                            "phase" : 0.0,
                                          },
                    "myo_inositol"      : { "spins" : 1, 
                                            "ppm"   : 3.51,   
                                            "area"  : 1.0,
                                            "phase" : 0.0,
                                          },
                    "h2o"               : { "spins"  : 2, 
                                            "ppm"   : 4.69,   
                                            "area"  : 2.0,
                                            "phase" : 0.0,
                                          },
                  }
    
#     metabolites = { "n-acetylaspartate" : { "spins" : 3, 
#                                             "ppm"   : 2.0082,   
#                                             "area"  : 3.0,
#                                             "phase" : 0.0,
#                                           },
#                     "naag"              : { "spins" : 3, 
#                                             "ppm"   : 2.042,   
#                                             "area"  : 3.0,
#                                             "phase" : 0.0,
#                                           },
#                     "creatine"          : { "spins" : 3, 
#                                             "ppm"   : 3.0241,   
#                                             "area"  : 3.0,
#                                             "phase" : 0.0,
#                                           },
#                     "pcr"               : { "spins" : 3, 
#                                             "ppm"   : 3.0295,   
#                                             "area"  : 3.0,
#                                             "phase" : 0.0,
#                                           },
#                     "cr2"               : { "spins" : 3, 
#                                             "ppm"   : 3.90241,   
#                                             "area"  : 3.0,
#                                             "phase" : 0.0,
#                                           },
#                     "pcr2"               : { "spins" : 3, 
#                                             "ppm"   : 3.90295,   
#                                             "area"  : 3.0,
#                                             "phase" : 0.0,
#                                           },
#                     "choline"           : { "spins"  : 9, 
#                                             "ppm"   : 3.21,   
#                                             "area"  : 9.0,
#                                             "phase" : 0.0,
#                                           },
#                     "pcho"              : { "spins"  : 9, 
#                                             "ppm"   : 3.206,   
#                                             "area"  : 9.0,
#                                             "phase" : 0.0,
#                                           },
#                     "gpc"               : { "spins"  : 9, 
#                                             "ppm"   : 3.212,   
#                                             "area"  : 9.0,
#                                             "phase" : 0.0,
#                                           },
#                     "glu"               : { "spins" : 3, 
#                                             "ppm"   : 2.35,   
#                                             "area"  : 3.0,
#                                             "phase" : 0.0,
#                                           },
#                     "gln"               : { "spins" : 3, 
#                                             "ppm"   : 2.43,   
#                                             "area"  : 3.0,
#                                             "phase" : 0.0,
#                                           },
#                     "tau"               : { "spins" : 3, 
#                                             "ppm"   : 3.4,   
#                                             "area"  : 3.0,
#                                             "phase" : 0.0,
#                                           },
#                     "glc"               : { "spins" : 3, 
#                                             "ppm"   : 3.7,   
#                                             "area"  : 3.0,
#                                             "phase" : 0.0,
#                                           },
#                     "h2o"               : { "spins"  : 2, 
#                                             "ppm"   : 4.69,   
#                                             "area"  : 2.0,
#                                             "phase" : 0.0,
#                                           },
#                   }

    deflated_metabolites = [ ]
    for name, metabolite in metabolites.items():
        d = { "name" : name,
              "spins" : metabolite["spins"],
              "dims" : [0, 0, 0],
              "group" : [0],
              "ppms" : [metabolite["ppm"]],
              "areas" : [metabolite["area"]],
              "phases" : [metabolite["phase"]],
            }
        deflated_metabolites.append(d)

             
    d = { "source" : "default",
          "source_id" : "default",
          "comment" : "This is a typical 1H singlet prior basis set.",
          "nucleus" : "1H",
          "seqte"  : 0.07, 
          "prior_metabolites" : deflated_metabolites,
        }

    return mrs_prior.Prior(d)


def _configure_combo(control, choices, selection=''):
        lines = [choices[key] for key in list(choices.keys())]
        control.Clear()        
        control.AppendItems( lines )
        if selection in lines:
            control.SetStringSelection(selection)
        else: 
            control.SetStringSelection(lines[0])


def _paired_event(obj_min, obj_max):
        val_min = obj_min.GetValue()
        val_max = obj_max.GetValue()
        pmin = min(val_min, val_max)
        pmax = max(val_min, val_max)
        obj_min.SetValue(pmin)
        obj_max.SetValue(pmax)
        return pmin, pmax

#------------------------------------------------------------------------------
#
#  Tab VOIGT 
#
#------------------------------------------------------------------------------

class TabVoigt(tab_base.Tab, voigt.PanelVoigtUI):
    
    # self-identify tab to notebook, value does not matter, its presence is sufficient.
    IS_VOIGT = True

    def __init__(self, tab_dataset, top, block):
        
        voigt.PanelVoigtUI.__init__(self, tab_dataset.NotebookDataset)
        tab_base.Tab.__init__(self, tab_dataset, top, prefs_module.PrefsVoigt)

        self.top = top               # application frame

        dataset = self.dataset

        # need to set this to a dummy until after initialize_controls
        self.block = block_fit_voigt.BlockFitVoigt()

        # Plotting is disabled during some of init. That's because the plot
        # isn't ready to plot, but the population of some controls 
        # (e.g. spin controls on the water filter panel) fires their 
        # respective change event which triggers a call to plot(). This
        # appears to happen only under Windows; it might be a Windows-specific
        # bug.
        # In any case, skipping some calls to plot() will speed things up. =)  
        self._plotting_enabled = False 
        self.plot_results = None
        
        # With large basis sets the automatic initial value update can be 
        # tedious, these values (and a check box on the GUI) give me some
        # control over that frustration 
        self.auto_update_initial_values = True
        self.update_button_hit          = False

        self.initialize_controls()

        self.block = block
        
        do_default = False
        if self.block is None:
            do_default = True
        elif not self.block.set.prior.metabolites: 
            # this was the original test for setting up default, but in Py3 it gorked is block==None
            do_default =True

        if do_default:
            # Create some default prior metabs and select them in the list.
            # This ensures that something will display in the plot when users 
            # first open this tab.
            self.block.set.prior = _create_default_prior()
            self.block.set.prior.calculate_full_basis_set(self.block.set.prior_ppm_start, 
                                                          self.block.set.prior_ppm_end, 
                                                          dataset)
            self.block.set.prior_list = block.set.prior.basis_set_names
            self.dynamic_metabolite_list.set_new_values(False)
            self.dynamic_metabolite_list.select_all()
            self.on_dynamic_metabolite_list()
        self.block.check_parameter_dimensions(self.dataset)
            
        self.populate_controls()

        # FIXME - this is temporary until I decide to release LMFIT stuff        
        # self.PanelOptimizeConstraintChecks.Hide()
        # FIXME - end
        
        self._plotting_enabled = True

        #------------------------------------------------------------
        # Setup the canvases
        self.process_and_plot(entry='plot_refresh')
        
        # If the sash position isn't recorded in the INI file, we use the  self.block.set.prior_list
        # arbitrary-ish value of 550.
        if not self._prefs.sash_position:
            self._prefs.sash_position = 550
        
        # Under OS X, wx sets the sash position to 10 (why 10?) *after*
        # this method is done. So setting the sash position here does no
        # good. We use wx.CallAfter() to (a) set the sash position and
        # (b) fake an EVT_SPLITTER_SASH_POS_CHANGED.
        wx.CallAfter(self.SplitterWindow.SetSashPosition, self._prefs.sash_position, True)
        wx.CallAfter(self.on_splitter)


#     @property
#     def block(self):
#         return self.dataset.blocks["fit"]


    #=======================================================
    #
    #           GUI Setup Handlers 
    #
    #=======================================================

    def initialize_controls(self):
        """ 
        This methods goes through the widgets and sets up certain sizes
        and constraints for those widgets. This method does not set the 
        value of any widget except insofar that it is outside a min/max
        range as those are being set up. 
        
        Use populate_controls() to set the values of the widgets from
        a data object.
        """
        
        # calculate a few useful values
        
        dataset = self.dataset
        dim0, dim1, dim2, dim3 = dataset.spectral_dims
        sw      = dataset.sw
        maxppm  = dataset.pts2ppm(0)
        minppm  = dataset.pts2ppm(dim0-1)
        ppmlim  = (minppm, maxppm)
        dim0lim = (0,dim0-1)
        
        
        # The many spin controls on various tabs need configuration of 
        # their size, # of digits displayed, increment and min/max. 
        
        
        wx_util.configure_spin(self.FloatPriorPpmStart, 70, 3, 1.0, ppmlim)
        wx_util.configure_spin(self.FloatPriorPpmEnd,   70, 3, 1.0, ppmlim)
        
        wx_util.configure_spin(self.FloatB0ShiftValue, 70, 3, 0.25, (-sw,sw))

        wx_util.configure_spin(self.FloatInitialBaselineLowessWidth, 70, 3, 1, (0.1,1000.0))
        wx_util.configure_spin(self.FloatBaselineInitIgnoreWidth,    70, 3, 1, (3,200))

        wx_util.configure_spin(self.FloatLinewidthValue, 70, 3, 0.5, 
                       (constants.FitLineWidth.MIN,
                        constants.FitLineWidth.MAX))
        wx_util.configure_spin(self.FloatInitialLinewidthStart, 70, 3, 0.5, ppmlim)
        wx_util.configure_spin(self.FloatInitialLinewidthEnd,   70, 3, 0.5, ppmlim)
        wx_util.configure_spin(self.FloatInitialLinewidthMultiplier, 70, 3, 0.1, (0.001,100.0))
        wx_util.configure_spin(self.FloatInitialPhase0Value,    70, 3, 1.0, (-360,360.0))
        wx_util.configure_spin(self.FloatInitialPhase1Value,    70, 3, 1.0, (-5000,5000.0))
        wx_util.configure_spin(self.SpinInitialKOLinewidthMinimum, 50, None, None, (10,100))
        wx_util.configure_spin(self.SpinInitialKOPoints,        50, None, None, (0,int(dim0/2)))

        wx_util.configure_spin(self.FloatBaselineSmoothingWidth, 70, 3, 5,
                       (constants.FitBaselineLowessWindowSize.MIN,
                        constants.FitBaselineLowessWindowSize.MAX))
        wx_util.configure_spin(self.FloatBaselineUnderestimate,  70, 3, 5,
                       (constants.FitBaselineUnderestimation.MIN,
                        constants.FitBaselineUnderestimation.MAX))
        wx_util.configure_spin(self.FloatBaselineUnderestimateLast,  70, 3, 5,
                       (constants.FitBaselineUnderestimationLast.MIN,
                        constants.FitBaselineUnderestimationLast.MAX))
        wx_util.configure_spin(self.SpinBaselineUnderestimateSteps, 70, None, None, (2,10))
        wx_util.configure_spin(self.SpinBaselineSplineNknots,    70, None, None, (1,100))
        wx_util.configure_spin(self.FloatBaselineSplineSpacing,  70, 3, 1, (1,dim0-1))    
        wx_util.configure_spin(self.SpinBaselineSplineOrder,     70, None, None,
                       (constants.FitBaselineBsplineOrder.MIN,
                        constants.FitBaselineBsplineOrder.MAX))
        wx_util.configure_spin(self.FloatBaselineWaveletMinDyad, 70, 3, 1, (1,1000))    

        wx_util.configure_spin(self.FloatMmolSingleDatasetStartArea,  70, 5, 0.1, (0.00001,10000.0))
        wx_util.configure_spin(self.FloatMmolSingleDatasetStartFudge, 70, 5, 0.1, (0.0000001,1000000.0))
        wx_util.configure_spin(self.FloatMmolSingleDatasetPpmStart,   70, 5, 0.1, ppmlim)
        wx_util.configure_spin(self.FloatMmolSingleDatasetPpmEnd,     70, 5, 0.1, ppmlim)
        wx_util.configure_spin(self.FloatMmolSingleDatasetLimitMax,   70, 5, 0.1, (0.00001,10000.0))
        wx_util.configure_spin(self.FloatMmolSingleDatasetLimitMin,   70, 5, 0.1, (0.00001,10000.0))

        wx_util.configure_spin(self.SpinOptimizeGlobalIterations,   70, None, None,
                       (constants.FitOptimizationGlobalIterations.MIN,
                        constants.FitOptimizationGlobalIterations.MAX))
        wx_util.configure_spin(self.FloatOptimizeStopTolerance,     70, 3, None,
                       (constants.FirtOptimizationStopTolerance.MIN,
                        constants.FirtOptimizationStopTolerance.MAX))
        self.FloatOptimizeStopTolerance.multiplier = 10
        wx_util.configure_spin(self.SpinOptimizeMaxIterations,      50, None, None,
                       (constants.FitOptimizationAlgorithmIterations.MIN,
                        constants.FitOptimizationAlgorithmIterations.MAX))
        wx_util.configure_spin(self.FloatOptimizeBoundsAreaMax,   70, 3, 10,
                       (constants.FitOptimizationAmplitudeBoundPercentMax.MIN,
                        constants.FitOptimizationAmplitudeBoundPercentMax.MAX))
        wx_util.configure_spin(self.FloatOptimizeBoundsAreaMin,   70, 3, 10,
                       (constants.FitOptimizationAmplitudeBoundPercentMin.MIN,
                        constants.FitOptimizationAmplitudeBoundPercentMin.MAX))
        wx_util.configure_spin(self.FloatOptimizeBoundsAreaMaxSmall,   70, 3, 10,
                       (constants.FitOptimizationAmplitudeBoundPercentMax.MIN,
                        constants.FitOptimizationAmplitudeBoundPercentMax.MAX))
        wx_util.configure_spin(self.FloatOptimizeBoundsAreaMinSmall,   70, 3, 10,
                       (constants.FitOptimizationAmplitudeBoundPercentMin.MIN,
                        constants.FitOptimizationAmplitudeBoundPercentMin.MAX))
        # Range is # of degrees in a circle
        wx_util.configure_spin(self.FloatOptimizeBoundsRangePhase0,  70, 3, 10, (1, 360))
        wx_util.configure_spin(self.FloatOptimizeBoundsRangePpm, 70, 3, 1,
                       (constants.FitOptimizationFrequency.MIN,
                        constants.FitOptimizationFrequency.MAX))
        wx_util.configure_spin(self.FloatOptimizeBoundsRangePpmSmall, 70, 3, 1,
                       (constants.FitOptimizationFrequency.MIN,
                        constants.FitOptimizationFrequency.MAX))
        wx_util.configure_spin(self.FloatOptimizeBoundsRangePhase1, 70, 3, 100,
                       (constants.FitOptimizationPhase1.MIN,
                        constants.FitOptimizationPhase1.MAX))
        wx_util.configure_spin(self.FloatOptimizeBoundsMinLinewidth, 70, 3, 0.01,
                       (constants.FitOptimizationTaTb.MIN,
                        constants.FitOptimizationTaTb.MAX))
        wx_util.configure_spin(self.FloatOptimizeBoundsMaxLinewidth, 70, 3, 0.05,
                       (constants.FitOptimizationTaTb.MIN,
                        constants.FitOptimizationTaTb.MAX))
        wx_util.configure_spin(self.FloatOptimizeWeightsScaleFactor, 70, 5, None,
                       (constants.FitOptimizationAreaWeight.MIN,
                        constants.FitOptimizationAreaWeight.MAX))
        self.FloatOptimizeWeightsScaleFactor.multiplier = 10
        wx_util.configure_spin(self.FloatOptimizeWeightsWidthFactor, 70, 3, 1,
                       (constants.FitOptimizationLocalMultiplier.MIN,
                        constants.FitOptimizationLocalMultiplier.MAX))
        wx_util.configure_spin(self.FloatOptimizeWeightsWaterStart,  70, 3, 0.5, ppmlim)
        wx_util.configure_spin(self.FloatOptimizeWeightsWaterEnd,    70, 3, 0.5, ppmlim)
        wx_util.configure_spin(self.FloatOptimizeWeightsLipidStart,  70, 3, 0.5, ppmlim)
        wx_util.configure_spin(self.FloatOptimizeWeightsLipidEnd,    70, 3, 0.5, ppmlim)
        wx_util.configure_spin(self.FloatOptimizeWeightsSmallPeakFactor, 70, 3, 0.5, (0.001,1000.0))
        
        wx_util.configure_spin(self.FloatConfidenceAlpha,            70, 5, 0.1,
                       (constants.FitOptimizationConfidenceAlpha.MIN,
                        constants.FitOptimizationConfidenceAlpha.MAX))
        wx_util.configure_spin(self.FloatCramerRaoPpmStart, 70, 3, 0.5, ppmlim)
        wx_util.configure_spin(self.FloatCramerRaoPpmEnd,   70, 3, 0.5, ppmlim)
        
        # Now clear and repopulate all Combo widget pull down menus
        #
        # Note. found odd bug (at least on windows) when dealing with ComboBox ...
        # if I define three entries in wxGlade, then redefine them here, I get an odd
        # 'box overlay' in my widget that sometimes goes away when I make a pick and 
        # sometimes it just stays and blocks my view of the widget below. The work around
        # was to NOT define the three entries in wxGlage (leave it empty) and just 
        # define them here. Who knew.
        
#        _configure_combo(self.ComboPriorMaskSource,         constants.FitPriorMaskSource.choices)
        _configure_combo(self.ComboInitialB0ShiftMethod,     constants.FitInitialB0ShiftMethod.choices)
        _configure_combo(self.ComboInitialBaselineMethod,    constants.FitInitialBaselineMethod.choices)
        _configure_combo(self.ComboInitialSmallPeakAreas,    constants.FitInitialSmallPeakAreas.choices)
        _configure_combo(self.ComboInitialSmallPeakFreqs,    constants.FitInitialSmallPeakFreqs.choices)
        _configure_combo(self.ComboInitialLinewidthMethod,   constants.FitInitialLinewidthMethod.choices)
        _configure_combo(self.ComboInitialPhaseMethod,       constants.FitInitialPhaseMethod.choices)
        _configure_combo(self.ComboLineshapeModel,           FitLineshapeModel.choices)
        _configure_combo(self.ComboBaselineMethod,           constants.FitBaselineMethod.choices)
        _configure_combo(self.ComboBaselineUnderestimateMethod, constants.FitBaselineUnderestimateMethod.choices)
        _configure_combo(self.ComboMmolModel,                constants.FitMacromoleculeMethod.choices)
        _configure_combo(self.ComboMmolSingleDatasetInitVal, constants.FitMacromoleculeMethodInitVal.choices)
        _configure_combo(self.ComboOptimizeMethod,           constants.FitOptimizeMethod.choices)
        _configure_combo(self.ComboOptimizeWeightsMethod,    constants.FitOptimizeWeightsMethod.choices)
        
        lines = ['1', '2', '4', '8', '16', '32', '64']
        self.ComboBaselineWaveletScale.Clear()        
        self.ComboBaselineWaveletScale.AppendItems( lines )

        #------------------------------------------------------------
        # set up View plot panel 
        #------------------------------------------------------------

        self.view = plot_panel_voigt.PlotPanelVoigt(self.PanelViewVoigt, 
                                                    self,
                                                    self._tab_dataset,
                                                    naxes=4,
                                                    reversex=True,
                                                    zoom='span', 
                                                    reference=True,
                                                    middle=True,
                                                    zoom_button=1,
                                                    middle_button=3,
                                                    refs_button=2,
                                                    do_zoom_select_event=True,
                                                    do_zoom_motion_event=True,
                                                    do_refs_select_event=True,
                                                    do_refs_motion_event=True,
                                                    do_middle_select_event=True,
                                                    do_middle_motion_event=True,
                                                    do_scroll_event=True,
                                                    props_zoom=dict(alpha=0.2, facecolor='yellow'),
                                                    props_cursor=dict(alpha=0.2, facecolor='gray'),
                                                    xscale_bump=0.0,
                                                    yscale_bump=0.05,
                                                    data = [],
                                                    prefs=self._prefs,
                                                    dataset=self.dataset,
                                                    )
        
        self.view.figure.set_size_inches(2.0,4.0)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.view, 1, wx.LEFT | wx.TOP | wx.EXPAND)
        self.PanelViewVoigt.SetSizer(sizer)
        self.view.Fit()    

        self.view.change_naxes(self._prefs.n_plots)

        #------------------------------------------------------------
        # set up Results output HTML widget
        #------------------------------------------------------------

        html_sizer = self.LabelHtmlPlaceholder.GetContainingSizer()
        parent = self.LabelHtmlPlaceholder.GetParent()
        self.LabelHtmlPlaceholder.Destroy()
        self.results_html_ctrl = wx.html.HtmlWindow(parent)
        html_sizer.Add(self.results_html_ctrl, 1, wx.EXPAND|wx.ALIGN_TOP)


        #------------------------------------------------------------
        # set up Metabolite dynamic list
        #------------------------------------------------------------

        # The list grid sizer is marked so we can find it at run-time
        self.MetaboliteGridSizer = self.LabelMetabolites.GetContainingSizer()
        _inner_notebook = self.LabelMetabolites.GetParent()
        self.LabelMetabolites.Destroy()
        
        # Add headings to the first row of the grid sizer.
        self.MetaboliteGridSizer.Clear()
        self.MetaboliteGridSizer.SetRows(1)
        headings = ( "Metabolites", 
                     "Area Scale\nFactor", 
                     "Peak PPM\nLocation",
                     "Peak Search\nWidth [ppm]",
                     "Use DB\nPPM Value",
                     "Fixed Ta\nValue [sec]",
                     "Peak Search\nPhase0 [deg]")
        
        for heading in headings:
            label = wx.StaticText(_inner_notebook, label=heading, style=wx.ALIGN_CENTRE)
            self.MetaboliteGridSizer.Add(label, 0, wx.EXPAND|wx.ALIGN_CENTER_VERTICAL)

        self.dynamic_metabolite_list = \
            dynamic_list_voigt_metabolite.DynamicListVoigtMetabolite( self.PanelMetaboliteLines,
                                                                      self,
                                                                      self.MetaboliteGridSizer,
                                                                      dataset,
                                                                      self.on_dynamic_metabolite_list)


#        #------------------------------------------------------------
#        # set up Macromolecule dynamic list
#        #------------------------------------------------------------
#        
#        # The list grid sizer is marked so we can find it at run-time
#        self.MacromoleculeGridSizer = self.LabelMacromolecules.GetContainingSizer()
#        _inner_notebook = self.LabelMacromolecules.GetParent()
#        self.LabelMacromolecules.Destroy()
#       
#        # Add headings to the first row of the grid sizer.
#        self.MacromoleculeGridSizer.Clear()
#        self.MacromoleculeGridSizer.SetRows(1)
#        headings = (None, "Peaks\n[ppm]", "+/- Range\n[PPM]", "Area", 
#                    "+/- Range",  "Linewidth\n[Hz]", "+/- Range\n[Hz]")
#        
#        for heading in headings:
#            if heading:
#                label = wx.StaticText(_inner_notebook, label=heading, style=wx.ALIGN_CENTRE)
#                self.MacromoleculeGridSizer.Add(label, 0, wx.EXPAND|wx.ALIGN_CENTER_VERTICAL)
#            else:
#                self.MacromoleculeGridSizer.AddSpacer( 1 )        
        

    def populate_controls(self, preset=False):
        """ 
        Populates the widgets with relevant values from the data object. 
        It's meant to be called when a new data object is loaded.
        
        This function trusts that the data object it is given doesn't violate
        any rules. Whatever is in the data object gets slapped into the 
        controls, no questions asked. 
        
        This function is, however, smart enough to enable and disable 
        other widgets depending on settings.
        """
        dataset = self.dataset
        voxel   = self._tab_dataset.voxel
        fit     = self.block

        hpp = dataset.sw / dataset.spectral_dims[0]
        max_lw_sec   = fit.set.optimize_bounds_max_linewidth
        min_lw_hz, _ = util_spectral.voigt_width(max_lw_sec, max_lw_sec, dataset)
        min_lw_sec   = fit.set.optimize_bounds_min_linewidth
        max_lw_hz, _ = util_spectral.voigt_width(min_lw_sec, min_lw_sec, dataset)
        
        #------------------------------------------------------------
        # set up Voigt widgets
        #------------------------------------------------------------
        
        self.TextPriorInformationSource.SetValue(fit.set.prior.source_id)
        self.FloatPriorPpmStart.SetValue(fit.set.prior_ppm_start)
        self.FloatPriorPpmEnd.SetValue(fit.set.prior_ppm_end)
        self.CheckPriorCalculateCombinations.SetValue(fit.set.prior_calculate_combinations)

#        val = constants.FitPriorMaskSource.choices[fit.set.prior_mask_source]
#        self.ComboPriorMaskSource.SetStringSelection( val )

#        self.SpinPriorXrangeStart.SetValue(fit.set.prior_xrange[0])
#        self.SpinPriorXrangeEnd.SetValue(  fit.set.prior_xrange[1])
#        self.SpinPriorYrangeStart.SetValue(fit.set.prior_yrange[0])
#        self.SpinPriorYrangeEnd.SetValue(  fit.set.prior_yrange[1])
#        self.SpinPriorZrangeStart.SetValue(fit.set.prior_zrange[0])
#        self.SpinPriorZrangeEnd.SetValue(  fit.set.prior_zrange[1])
#        self.CheckPriorIgnoreMask.SetValue(fit.set.prior_ignore_mask)
#
#        if dataset.spectral_dims[2] > 1 or dataset.spectral_dims[3] > 1:
#            self.PanelMaskPrior.Show()
#        else:
#            self.PanelMaskPrior.Hide()

        val = constants.FitInitialB0ShiftMethod.choices[fit.set.initial_b0_shift_method]
        self.ComboInitialB0ShiftMethod.SetStringSelection( val )
        self.FloatB0ShiftValue.SetValue(self.dataset.get_frequency_shift(self._tab_dataset.voxel))
        val = constants.FitInitialBaselineMethod.choices[fit.set.initial_baseline_method]
        self.ComboInitialBaselineMethod.SetStringSelection( val )
        self.FloatInitialBaselineLowessWidth.SetValue(fit.set.initial_baseline_lowess_width)
        self.FloatBaselineInitIgnoreWidth.SetValue(fit.set.initial_baseline_lowess_ignore_width)

        self.CheckInitialCrChoSeparation.SetValue(fit.set.initial_cr_cho_separation)
        self.CheckInitialPeakSearchAbs.SetValue(fit.set.initial_peak_search_abs)
        val = constants.FitInitialSmallPeakAreas.choices[fit.set.initial_small_peak_areas]
        self.ComboInitialSmallPeakAreas.SetStringSelection( val )
        val = constants.FitInitialSmallPeakFreqs.choices[fit.set.initial_small_peak_freqs]
        self.ComboInitialSmallPeakFreqs.SetStringSelection( val )
        val = constants.FitInitialLinewidthMethod.choices[fit.set.initial_linewidth_method]
        self.ComboInitialLinewidthMethod.SetStringSelection( val )
        self.FloatLinewidthValue.SetValue(fit.set.initial_linewidth_value)
        self.FloatInitialLinewidthStart.SetValue(fit.set.initial_linewidth_range_start)
        self.FloatInitialLinewidthEnd.SetValue(fit.set.initial_linewidth_range_end)
        self.FloatInitialLinewidthMultiplier.SetValue(fit.set.initial_linewidth_fudge)
        val = constants.FitInitialPhaseMethod.choices[fit.set.initial_phase_method]
        self.ComboInitialPhaseMethod.SetStringSelection( val )
        self.FloatInitialPhase0Value.SetValue(self.dataset.get_phase_0(self._tab_dataset.voxel))
        self.FloatInitialPhase1Value.SetValue(self.dataset.get_phase_1(self._tab_dataset.voxel))
        self.CheckApplyKoFilter.SetValue(fit.set.initial_apply_ko_filter)
        self.SpinInitialKOLinewidthMinimum.SetValue(fit.set.initial_ko_linewidth_minimum)
        self.SpinInitialKOPoints.SetValue(fit.set.initial_ko_points)

        show = bool(fit.set.initial_baseline_method)
        self.PanelBaselineInitLowess.Show(show)

        if fit.set.initial_linewidth_method == constants.FitInitialLinewidthMethod.MANUAL:
            self.PanelInitialLinewidth.Hide()
        else:
            self.PanelInitialLinewidth.Show()
        self.PanelInitialValues.Layout()

        # metabolite dynamic list set below
        self.ComboLineshapeModel.SetStringSelection(fit.set.lineshape_model)

        method = constants.FitBaselineMethod.choices[fit.set.baseline_method]
        self.ComboBaselineMethod.SetStringSelection(method)
        # Fake a selection in the baseline combobox
        self.on_baseline_method()
        
        method2 = constants.FitBaselineUnderestimateMethod.choices[fit.set.baseline_underestimate_method]
        self.ComboBaselineUnderestimateMethod.SetStringSelection(method2)
        
        self.CheckBaselineSmoothingFlag.SetValue(       fit.set.baseline_smoothing_flag)
        self.CheckBaselineSkipLastSmooth.SetValue(      fit.set.baseline_skip_last_smooth)
        self.FloatBaselineSmoothingWidth.SetValue(      fit.set.baseline_smoothing_width)
        self.FloatBaselineUnderestimate.SetValue(       fit.set.baseline_underestimate)
        self.FloatBaselineUnderestimateLast.SetValue(   fit.set.baseline_underestimate_last)
        self.SpinBaselineUnderestimateSteps.SetValue(   fit.set.baseline_underestimate_steps)
        self.SpinBaselineSplineNknots.SetValue(         fit.set.baseline_spline_nknots)
        self.FloatBaselineSplineSpacing.SetValue(       fit.set.baseline_spline_spacing)
        self._set_spacing(fit.set.baseline_spline_spacing * hpp)
        self.SpinBaselineSplineOrder.SetValue(          fit.set.baseline_spline_order)
        self.ComboBaselineWaveletScale.SetStringSelection(  str(fit.set.baseline_wavelet_scale))
        self.FloatBaselineWaveletMinDyad.SetValue(      fit.set.baseline_wavelet_min_dyad)

        if fit.set.baseline_method == constants.FitBaselineMethod.BSPLINE_FIXED_KNOT:
            self.FloatBaselineSplineSpacing.Enable()
            self.SpinBaselineSplineNknots.Disable()
        elif fit.set.baseline_method == constants.FitBaselineMethod.BSPLINE_VARIABLE_KNOT:
            self.FloatBaselineSplineSpacing.Disable()
            self.SpinBaselineSplineNknots.Enable()

        if fit.set.baseline_underestimate_method == constants.FitBaselineUnderestimateMethod.FIRST_ITERATION_ONLY:
            self.LabelBaselineUnderestimateLast.Hide()
            self.FloatBaselineUnderestimateLast.Hide()
            self.LabelBaselineUnderestimateSteps.Hide()
            self.SpinBaselineUnderestimateSteps.Hide()
        elif fit.set.baseline_underestimate_method == constants.FitBaselineUnderestimateMethod.MULTI_STEP_UNDERESTIMATION:
            self.LabelBaselineUnderestimateLast.Show()
            self.FloatBaselineUnderestimateLast.Show()
            self.LabelBaselineUnderestimateSteps.Show()
            self.SpinBaselineUnderestimateSteps.Show()

        label = constants.FitMacromoleculeMethod.choices[fit.set.macromol_model]
        self.ComboMmolModel.SetStringSelection(label)
        fname = ''
        if fit.set.macromol_single_basis_dataset:
            block = fit.set.macromol_single_basis_dataset.blocks["raw"]
            fname = block.data_source
        self.TextMmolSingleDataset.SetValue(fname)
        if fit.set.macromol_model == constants.FitMacromoleculeMethod.SINGLE_BASIS_DATASET:
            self.PanelMmolSingleDataset.Show()
        else:
            self.PanelMmolSingleDataset.Hide()

        label = constants.FitMacromoleculeMethodInitVal.choices[fit.set.macromol_single_basis_dataset_initval]
        self.ComboMmolSingleDatasetInitVal.SetStringSelection(label)
        self.set_mmol_single()
        
        self.FloatMmolSingleDatasetStartArea.SetValue(fit.set.macromol_single_basis_dataset_start_area)
        self.FloatMmolSingleDatasetStartFudge.SetValue(fit.set.macromol_single_basis_dataset_start_fudge)
        self.FloatMmolSingleDatasetPpmStart.SetValue(fit.set.macromol_single_basis_dataset_ppm_start)
        self.FloatMmolSingleDatasetPpmEnd.SetValue(fit.set.macromol_single_basis_dataset_ppm_end)
        self.FloatMmolSingleDatasetLimitMax.SetValue(fit.set.macromol_single_basis_dataset_limit_max)
        self.FloatMmolSingleDatasetLimitMin.SetValue(fit.set.macromol_single_basis_dataset_limit_min)



#        # macromolecule dynamic list set below
#        self.ComboMmolMethod.SetStringSelection(fit.set.macromolecule_method)
#        if fit.set.macromolecule_reference_line in fit.set.prior.names:
#            self.ComboMmolReferenceLine.SetStringSelection(fit.set.macromolecule_reference_line)
#        elif 'n-acetylaspartate' in fit.set.prior.names:
#            self.ComboMmolReferenceLine.SetStringSelection('n-acetylaspartate')
#        else:
#            self.ComboMmolReferenceLine.SetStringSelection(fit.set.prior.names[0])  
#        if fit.set.macromolecule_method == 'None':
#            self.PanelMmolParameters.Hide()
#        else:
#            self.PanelMmolParameters.Show()

        
        val = constants.FitOptimizeMethod.choices[fit.set.optimize_method]
        self.ComboOptimizeMethod.SetStringSelection(    val )
        self.CheckOptimizeScalingFlag.SetValue(         fit.set.optimize_scaling_flag)
        self.SpinOptimizeGlobalIterations.SetValue(     fit.set.optimize_global_iterations)
        self.FloatOptimizeStopTolerance.SetValue(       fit.set.optimize_stop_tolerance)
        self.SpinOptimizeMaxIterations.SetValue(        fit.set.optimize_max_iterations)

        val = constants.FitOptimizeWeightsMethod.choices[fit.set.optimize_weights_method]
        self.ComboOptimizeWeightsMethod.SetStringSelection( val )
        self.FloatOptimizeWeightsScaleFactor.SetValue(  fit.set.optimize_weights_scale_factor)
        self.FloatOptimizeWeightsWidthFactor.SetValue(  fit.set.optimize_weights_width_factor)
        self.CheckOptimizeWeightsWaterFlag.SetValue(      fit.set.optimize_weights_water_flag)
        self.FloatOptimizeWeightsWaterStart.SetValue(     fit.set.optimize_weights_water_start)
        self.FloatOptimizeWeightsWaterEnd.SetValue(       fit.set.optimize_weights_water_end)
        self.CheckOptimizeWeightsLipidFlag.SetValue(      fit.set.optimize_weights_lipid_flag)
        self.FloatOptimizeWeightsLipidStart.SetValue(     fit.set.optimize_weights_lipid_start)
        self.FloatOptimizeWeightsLipidEnd.SetValue(       fit.set.optimize_weights_lipid_end)
        self.FloatOptimizeWeightsSmallPeakFactor.SetValue(fit.set.optimize_weights_small_peak_factor)

        if fit.set.optimize_weights_method == constants.FitOptimizeWeightsMethod.LOCAL_WEIGHTING:
            self.PanelLocalWeights.Show()
        else:
            self.PanelLocalWeights.Hide()
        
        self.FloatOptimizeBoundsAreaMax.SetValue(       fit.set.optimize_bounds_area_max)
        self.FloatOptimizeBoundsAreaMin.SetValue(       fit.set.optimize_bounds_area_min)
        self.FloatOptimizeBoundsAreaMaxSmall.SetValue(  fit.set.optimize_bounds_area_max_small)
        self.FloatOptimizeBoundsAreaMinSmall.SetValue(  fit.set.optimize_bounds_area_min_small)
        self.CheckEnableBoundsAreaSmall.SetValue(       fit.set.optimize_enable_bounds_area_small)
        self.FloatOptimizeBoundsRangePhase0.SetValue(   fit.set.optimize_bounds_range_phase0)
        self.FloatOptimizeBoundsRangePpm.SetValue(      fit.set.optimize_bounds_range_ppm)
        self.FloatOptimizeBoundsRangePpmSmall.SetValue( fit.set.optimize_bounds_range_ppm_small)
        self.CheckEnableBoundsPpmSmall.SetValue(        fit.set.optimize_enable_bounds_ppm_small)
        self.FloatOptimizeBoundsRangePhase1.SetValue(   fit.set.optimize_bounds_range_phase1)
        self.FloatOptimizeBoundsMinLinewidth.SetValue(  fit.set.optimize_bounds_min_linewidth)
        self._set_max_linewidth(max_lw_hz)
        self.FloatOptimizeBoundsMaxLinewidth.SetValue(  fit.set.optimize_bounds_max_linewidth)
        self._set_min_linewidth(min_lw_hz)
        
        self.CheckOptimizeConstrainPpmNaaNaag.SetValue( fit.set.optimize_constrain_ppm_naa_naag)
        self.CheckOptimizeConstrainPpmCrPcr.SetValue(   fit.set.optimize_constrain_ppm_cr_pcr)
        self.CheckOptimizeConstrainPpmGpcPcho.SetValue(    fit.set.optimize_constrain_ppm_gpc_pcho)
        self.CheckOptimizeConstrainPpmCr2Pcr2.SetValue( fit.set.optimize_constrain_ppm_cr2_pcr2)
        self.CheckOptimizeConstrainPpmGluGln.SetValue( fit.set.optimize_constrain_ppm_glu_gln)
        self.CheckOptimizeConstrainPpmTauGlc.SetValue( fit.set.optimize_constrain_ppm_tau_glc)

        self.CheckConfidenceIntervalsFlag.SetValue(fit.set.confidence_intervals_flag)
        self.FloatConfidenceAlpha.SetValue(fit.set.confidence_alpha)
        self.CheckConfidenceAreaFlagFlag.SetValue(fit.set.confidence_area_flag)
        self.CheckCofidencePpmFlag.SetValue(fit.set.confidence_ppm_flag)
        self.CheckConfidenceLinewidthFlag.SetValue(fit.set.confidence_linewidth_flag)
        self.CheckConfidencePhaseFlag.SetValue(fit.set.confidence_phase_flag)
        self.CheckCramerRaoFlag.SetValue(fit.set.cramer_rao_flag)
        self.FloatCramerRaoPpmStart.SetValue(fit.set.cramer_rao_ppm_start)
        self.FloatCramerRaoPpmEnd.SetValue(fit.set.cramer_rao_ppm_end)
        
        # This value is only used for the GUI and is not saved with XML output
        self.CheckAutoInitialValues.SetValue(self.auto_update_initial_values)
        
        # match old settings to new if they exist, but be sure that all other
        # parameters are updated if old values are NOT part of new one
        self.dynamic_metabolite_list.set_new_values(respect_current=(not preset), preset=preset)
        vals = self.dynamic_metabolite_list.get_values()
        names, area_scales, peak_ppms, search_ppms, db_ppms, fix_t2, search_ph0 = vals 
                
        fit.set.prior_list       = names
        fit.set.prior_area_scale = area_scales
        fit.set.prior_peak_ppm   = peak_ppms
        fit.set.prior_search_ppm = search_ppms
        fit.set.prior_db_ppm     = db_ppms
        fit.set.prior_fix_t2     = fix_t2
        fit.set.prior_search_ph0 = search_ph0

#                
#        self.dynamic_macromolecule_list = DynamicMacromoleculeList( self.PanelMmolLines,
#                                                                    self.MacromoleculeGridSizer,
#                                                                    fit.set.macromolecule,
#                                                                    dataset)


    def set_mmol_single(self):
        """ various widgets need to be active or disabled depending on the initval method """
        set = self.block.set
        if set.macromol_single_basis_dataset_initval == constants.FitMacromoleculeMethodInitVal.MANUAL:
            self.FloatMmolSingleDatasetStartArea.Enable(True)
            self.FloatMmolSingleDatasetStartFudge.Enable(False)
            self.FloatMmolSingleDatasetPpmStart.Enable(False)
            self.FloatMmolSingleDatasetPpmEnd.Enable(False)
        elif set.macromol_single_basis_dataset_initval == constants.FitMacromoleculeMethodInitVal.ONEPOINT:
            self.FloatMmolSingleDatasetStartArea.Enable(False)
            self.FloatMmolSingleDatasetStartFudge.Enable(True)
            self.FloatMmolSingleDatasetPpmStart.Enable(True)
            self.FloatMmolSingleDatasetPpmEnd.Enable(False)
        elif set.macromol_single_basis_dataset_initval == constants.FitMacromoleculeMethodInitVal.REGION:
            self.FloatMmolSingleDatasetStartArea.Enable(False)
            self.FloatMmolSingleDatasetStartFudge.Enable(True)
            self.FloatMmolSingleDatasetPpmStart.Enable(True)
            self.FloatMmolSingleDatasetPpmEnd.Enable(True)
        
        
        

    #=======================================================
    #
    #           Global and Menu Event Handlers  
    #
    #=======================================================

    def on_activation(self):
        tab_base.Tab.on_activation(self)

        voxel   = self._tab_dataset.voxel
        phase_0 = self.dataset.get_phase_0(voxel)
        phase_1 = self.dataset.get_phase_1(voxel)

        self.FloatB0ShiftValue.SetValue(self.dataset.get_frequency_shift(voxel))
        self.FloatInitialPhase0Value.SetValue(phase_0)
        self.FloatInitialPhase1Value.SetValue(phase_1)

        self.view.set_phase_0(phase_0, absolute=True)
        self.view.set_phase_1(phase_1, absolute=True)
        
        
    def on_menu_view_option(self, event):
        event_id = event.GetId()
        
        if self._prefs.handle_event(event_id):
            if event_id in (util_menu.ViewIdsVoigt.ZERO_LINE_SHOW,
                            util_menu.ViewIdsVoigt.ZERO_LINE_TOP,
                            util_menu.ViewIdsVoigt.ZERO_LINE_MIDDLE,
                            util_menu.ViewIdsVoigt.ZERO_LINE_BOTTOM,
                            util_menu.ViewIdsVoigt.XAXIS_SHOW,
                           ):
                self.view.update_axes()
                self.view.canvas.draw()

            elif event_id in (util_menu.ViewIdsVoigt.XAXIS_PPM,
                              util_menu.ViewIdsVoigt.XAXIS_HERTZ,
                             ):
                self.view.update(no_draw=True)
                self.view.set_phase_0(0.0)
                self.view.canvas.draw()

            elif event_id in (util_menu.ViewIdsVoigt.AREA_CALC_PLOT_A,
                              util_menu.ViewIdsVoigt.AREA_CALC_PLOT_B,
                              util_menu.ViewIdsVoigt.AREA_CALC_PLOT_C,
                              util_menu.ViewIdsVoigt.AREA_CALC_PLOT_D,
                             ):
                area, rms = self.view.calculate_area()
                if self._prefs.area_calc_plot_a:
                    index = 0
                elif  self._prefs.area_calc_plot_b:
                    index = 1
                elif  self._prefs.area_calc_plot_c:
                    index = 2
                elif  self._prefs.area_calc_plot_d:
                    index = 3
                self.top.statusbar.SetStatusText(self.build_area_text(area[index], rms[index]), 3)

            elif event_id in (util_menu.ViewIdsVoigt.N_PLOTS_1,
                              util_menu.ViewIdsVoigt.N_PLOTS_2,
                              util_menu.ViewIdsVoigt.N_PLOTS_3,
                              util_menu.ViewIdsVoigt.N_PLOTS_4,
                             ):
                self.view.change_naxes(self._prefs.n_plots)



    def on_menu_view_output(self, event):
        # EVENT_TO_TYPE maps menu items to file types
        EVENT_TO_TYPE = { util_menu.ViewIdsVoigt.VIEW_TO_PNG : "png",
                          util_menu.ViewIdsVoigt.VIEW_TO_SVG : "svg",
                          util_menu.ViewIdsVoigt.VIEW_TO_EPS : "eps",
                          util_menu.ViewIdsVoigt.VIEW_TO_PDF : "pdf",
                        }
                        
        type_ = EVENT_TO_TYPE[event.GetId()]
        
        ini_name = "voigt_output_as_image"
        default_path = util_analysis_config.get_path(ini_name)

        filetype_filter = "%s files (*.%s)|*.%s" % (type_.upper(), type_, type_)
        
        filename = common_dialogs.save_as(default_path=default_path,
                                          filetype_filter=filetype_filter)

        if filename:
            figure = self.view.figure
            try:
                figure.savefig( filename,
                                dpi=300, 
                                facecolor='w', 
                                edgecolor='w',
                                orientation='portrait', 
                                papertype='letter', 
                                format=None,
                                transparent=False)
            except IOError:
                msg = """I can't write the file "%s".""" % filename
                common_dialogs.message(msg, style=common_dialogs.E_OK)
            else:
                path, _ = os.path.split(filename)
                util_analysis_config.set_path(ini_name, path)


    def on_menu_view_results(self, event):
        """
        Lots of options for saving results to files.
        
        Two main flavors : Text (csv) and Images (PDF and PNG)
        
        Text output can be for current voxel, or all voxels in the dataset. This
        is a pretty simplistic method. It outputs whatever is in the results 
        array, zeros or whatever.
        
        Image output come in three layouts - LCModel-like, 2 plots just like 
        what are on the GUI, and 4 plots just like what are on the GUI. All three
        layouts can be saved to either PDF or PNG.
        
        """
        
        # EVENT_TO_TYPE maps menu items to file types
        EVENT_TO_TYPE = { util_menu.ViewIdsVoigt.RESULTS_CSV_CURRENT   : ("text","current","csv"),
                          util_menu.ViewIdsVoigt.RESULTS_CSV_ALL       : ("text","all","csv"),
                          util_menu.ViewIdsVoigt.RESULTS_LCM_PDF       : ("image","lcm","pdf"),
                          util_menu.ViewIdsVoigt.RESULTS_LCM_PDF_MULTI : ("image","lcm_multi","pdf"),
                          util_menu.ViewIdsVoigt.RESULTS_LCM_PNG       : ("image","lcm","png"),
                          util_menu.ViewIdsVoigt.RESULTS_BRP512_PNG    : ("image","brp512","png"),
                          util_menu.ViewIdsVoigt.RESULTS_BRP1024_PNG   : ("image","brp1024","png"),
                          util_menu.ViewIdsVoigt.RESULTS_2PLOT_PDF     : ("image","2plot","pdf"),
                          util_menu.ViewIdsVoigt.RESULTS_2PLOT_PNG     : ("image","2plot","png"),
                          util_menu.ViewIdsVoigt.RESULTS_4PLOT_PDF     : ("image","4plot","pdf"),
                          util_menu.ViewIdsVoigt.RESULTS_4PLOT_PNG     : ("image","4plot","png"),
                        }

        flavor, layout, type_ = EVENT_TO_TYPE[event.GetId()]

        if flavor == 'text':

            ini_name = "voigt_output_as_csv"
            default_path_name = util_analysis_config.get_path(ini_name)
            default_path, default_fname = os.path.split(default_path_name)
            
            if type_ == 'csv':
                filetype_filter = "CSV files (*.csv)|*.csv|TXT files (*.txt)|*.txt"
            else:
                filetype_filter = "%s files (*.%s)|*.%s" % (type_.upper(), type_, type_)
            
            filename = common_dialogs.save_as(default_path     = default_path,
                                              filetype_filter  = filetype_filter,
                                              default_filename = default_fname)
            if filename:
                try:
                    if layout == 'current':
                        voxels = [self._tab_dataset.voxel,]
                    else:
                        voxels = self.dataset.all_voxels
                    
                    self.output_results_to_csv(voxels, filename)
                    
                except IOError:
                    msg = """I can't write the file "%s".""" % filename
                    common_dialogs.message(msg, style=common_dialogs.E_OK)
                else:
                    # We saved results, so we write the path to the INI file.
                    util_analysis_config.set_path(ini_name, filename)
            #else:
                # User clicked cancel on the "save as" dialog

            
        elif flavor == 'image':       # image layout result
        
            ini_name = "voigt_results_output_as_image"
            default_path = util_analysis_config.get_path(ini_name)
    
            filetype_filter = "%s files (*.%s)|*.%s" % (type_.upper(), type_, type_)
            
            filename = common_dialogs.save_as(default_path    = default_path,
                                              filetype_filter = filetype_filter)

            voxel = self._tab_dataset.voxel

            if layout == 'lcm':
                fig_call = figure_layouts.lcm_like
                nobase = False
                fixphase = True
                dpi = 300
            elif layout == 'lcm_multi':
                fig_call = figure_layouts.lcm_multipage_pdf
                nobase = False
                fixphase = True
                dpi = 300
            elif layout == '2plot':
                fig_call = figure_layouts.analysis_plot2
                nobase = True
                fixphase = True
                dpi = 300
            elif layout == '4plot':
                fig_call = figure_layouts.debug_plot4
                nobase = True
                fixphase = True
                dpi = 300
            elif layout == 'brp512':  
                fig_call = figure_layouts.analysis_brp_generic
                nobase = True
                fixphase = True
                dpi = 64            # for 8"x8" figure -> 512x512 resolution
            elif layout == 'brp1024':  
                fig_call = figure_layouts.analysis_brp_generic
                nobase = True
                fixphase = True
                dpi = 128           # for 8"x8" figure -> 1024x1024 resolution
            
            if filename:
                
                viffpath = 'Analysis - Fit Tab'
                vespa_version = util_misc.get_vespa_version()
                timestamp = ''      # fig_call will provide
                minplot   = 0.1
                maxplot   = 4.9
                fontname  = 'Courier New'  #'DejaVu Sans'  #'Ariel'  #'Courier New'
    
                figure = fig_call(self.dataset, 
                                  viffpath=viffpath, 
                                  vespa_version=vespa_version,
                                  timestamp=timestamp,
                                  fontname=fontname,
                                  minplot=minplot,
                                  maxplot=maxplot,
                                  nobase=nobase,
                                  extfig=None,
                                  fixphase=fixphase,
                                  verbose=False, 
                                  debug=False,
                                  quantvals=False,
                                  voxel=voxel)
                
                if len(figure) == 1:
                
                    try:
                        figure[0].savefig( filename,
                                        dpi=dpi, 
                                        pad_inches=0.1,
                                        facecolor=figure[0].get_facecolor(), 
                                        edgecolor='none')
                        
                    except IOError:
                        msg = """I can't write the file "%s".""" % filename
                        common_dialogs.message(msg, style=common_dialogs.E_OK)
                    else:
                        path, _ = os.path.split(filename)
                        util_analysis_config.set_path(ini_name, path)

                else:
                    
                    try:
                        if os.path.splitext(filename)[-1].lower() != '.pdf':
                            filename += '.pdf'
                    
                        # Create the PdfPages object to which we will save the pages:
                        # The with statement makes sure that the PdfPages object is closed properly at
                        # the end of the block, even if an Exception occurs.
                        with PdfPages(filename) as pdf:
                            for fig in figure:
                                pdf.savefig(fig)
                        
                            # We can also set the file's metadata via the PdfPages object:
                            today = datetime.date.today()
                            d = pdf.infodict()
                            d['Title'] = 'Vespa Multi-page LCM-like Output'
                            d['Author'] = 'Brian J. Soher'
                            d['Subject'] = 'Vespa results output'
                            d['Keywords'] = 'PdfPages multipage Vespa output lcm-like'
                            d['CreationDate'] = datetime.datetime(today.year, today.month, today.day)
                            d['ModDate'] = datetime.datetime.today()        
                        # print "Testing completed successfully ... returning."                    
                    
                    except IOError:
                        msg = """I can't write the multi-page PDF file "%s".""" % filename
                        common_dialogs.message(msg, style=common_dialogs.E_OK)
                    else:
                        path, _ = os.path.split(filename)
                        util_analysis_config.set_path(ini_name, path)

        else:
            msg = 'this flavor of results output is not recognized'


    def on_menu_view_debug(self, event):
        """
        Lots of options for debugging fitting results to screen or to files.
        
        """
        EVENT_TO_TYPE = { util_menu.ViewIdsVoigt.DEBUG_POPUP_CURRENT : ("text","current","popup"),
                          util_menu.ViewIdsVoigt.DEBUG_CSV_CURRENT   : ("text","current","csv"),
                        }

        flavor, layout, type_ = EVENT_TO_TYPE[event.GetId()]

        if flavor == 'text':
            
            if type_ == 'popup':
                
                raw              = self.dataset.blocks["raw"]
                voxel            = self._tab_dataset.voxel
                data_source      = raw.get_data_source(voxel)
                dataset_filename = self.dataset.dataset_filename

                lines = self.dataset.debug_results_as_popup(voxel,
                                                            source=data_source,
                                                            dsetname=dataset_filename,
                                                          )
                lines = "\n".join(lines)
                
                wx_util.display_text_as_file(lines)
                
            elif type_ == 'csv':

                ini_name = "voigt_debug_as_csv"
                default_path_name = util_analysis_config.get_path(ini_name)
                default_path, default_fname = os.path.split(default_path_name)
                
                filetype_filter = "CSV files (*.csv)|*.csv|TXT files (*.txt)|*.txt"
                
                filename = common_dialogs.save_as(default_path     = default_path,
                                                  filetype_filter  = filetype_filter,
                                                  default_filename = default_fname)
                if filename:
                    try:
                        if layout == 'current':
                            voxels = [self._tab_dataset.voxel,]
                        else:
                            voxels = self.dataset.all_voxels
                        
                        self.debug_results_to_csv(voxels, filename)
                        
                    except IOError:
                        msg = """I can't write the file "%s".""" % filename
                        common_dialogs.message(msg, style=common_dialogs.E_OK)
                    else:
                        # We saved results, so we write the path to the INI file.
                        util_analysis_config.set_path(ini_name, filename)
                #else:
                    # User clicked cancel on the "save as" dialog

        elif flavor == 'image':       # image layout result
            bob = 10
        else:
            msg = 'this flavor of results output is not recognized'
                        

    def on_menu_plot_x(self, event):
        event_id = event.GetId()

        if self._prefs.handle_event(event_id):
            do_break = False
            for i in range(4):
                if event_id in (util_menu.PlotXIds.DATA_TYPE_REAL[i],
                                util_menu.PlotXIds.DATA_TYPE_IMAGINARY[i],
                                util_menu.PlotXIds.DATA_TYPE_MAGNITUDE[i],
                                util_menu.PlotXIds.DATA_TYPE_SUMMED[i],
                               ):
                    if event_id == util_menu.PlotXIds.DATA_TYPE_REAL[i]:
                        self.view.set_data_type_real(index=[i])
                    elif event_id == util_menu.PlotXIds.DATA_TYPE_IMAGINARY[i]:
                        self.view.set_data_type_imaginary(index=[i])
                    elif event_id == util_menu.PlotXIds.DATA_TYPE_MAGNITUDE[i]:
                        self.view.set_data_type_magnitude(index=[i])
                    elif event_id == util_menu.PlotXIds.DATA_TYPE_SUMMED[i]:
                        self.view.set_data_type_summed(index=[i])
                    self.view.update(no_draw=True)
                    self.view.set_phase_0(0.0)
                    self.view.canvas.draw()
                    do_break = True

                if event_id in (util_menu.PlotXIds.RAW_DATA[i],
                                util_menu.PlotXIds.FITTED_DATA[i],
                                util_menu.PlotXIds.BASELINE[i],
                                util_menu.PlotXIds.COMBO_RAW_MINUS_FIT[i],
                                util_menu.PlotXIds.COMBO_RAW_MINUS_BASE[i],
                                util_menu.PlotXIds.COMBO_RAW_MINUS_FIT_MINUS_BASE[i],
                                util_menu.PlotXIds.COMBO_RAW_AND_FIT[i],
                                util_menu.PlotXIds.COMBO_RAW_AND_BASE[i],
                                util_menu.PlotXIds.COMBO_RAW_AND_FIT_PLUS_BASE[i],
                                util_menu.PlotXIds.COMBO_RAW_MINUS_BASE_AND_FIT[i],
                                util_menu.PlotXIds.COMBO_FIT_PLUS_BASE[i],
                                util_menu.PlotXIds.COMBO_FIT_AND_BASE[i],
                                util_menu.PlotXIds.COMBO_RAW_AND_INIT_MODEL[i],
                                util_menu.PlotXIds.COMBO_RAW_AND_WT_ARR[i],
                               ):
                    self.process_and_plot()
                    do_break = True
                    
                if do_break: 
                    break
                    
                   


    #=======================================================
    #
    #           Widget Event Handlers  
    #
    #=======================================================
    
    def on_splitter(self, event=None):
        # This is sometimes called programmatically, in which case event is None
        self._prefs.sash_position = self.SplitterWindow.GetSashPosition()

    
    def on_prior_information_from_database(self, event): 
        # Get the user to select an experiment
        dialog = dialog_experiment_browser.DialogExperimentBrowser(self, self.top.db)
        dialog.ShowModal()
        experiment_id = dialog.selected_experiment_id
        dialog.Destroy()

        if experiment_id:
            
            experiment = self.top.db.fetch_experiment(experiment_id)
            format     = ThirdPartyExportTypes.ANALYSIS_DIRECT
            dialog = dialog_mixout.DialogMixedMetaboliteOutput(self, 
                                                               experiment, 
                                                               self.dataset, 
                                                               format)
            dialog.ShowModal()
            prior = dialog.final_prior
            dialog.Destroy()
            
            if prior is not None:
                preview = self.top.db.fetch_experiment_preview(experiment_id)
                # Now load the new prior into the block & GUI
                self._apply_selected_prior(prior, preview.name)
                
                
    def on_prior_information_from_file(self, event): 
        ini_name = "load_prior"
        default_path = util_analysis_config.get_path(ini_name)
        filetype_filter = "Priors (*.xml,*.xml.gz,*.viff,*.vif)|*.xml;*.xml.gz;*.viff;*.vif"

        filename = common_dialogs.pickfile(default_path=default_path,
                                           filetype_filter=filetype_filter)
        
        if filename:
            msg = ""
            try:
                importer = util_import.PriorImporter(filename)
            except IOError:
                msg = """I can't read the file "%s".""" % filename
            except SyntaxError:
                msg = """The file "%s" isn't a valid prior export file.""" % filename
            else:
                # Time to rock and roll!
                wx.BeginBusyCursor()
                priors = importer.go()
                wx.EndBusyCursor()    

                if not priors:
                    msg = """The file "%s" didn't contain any exported priors.""" % filename
                
            if msg:
                common_dialogs.message(msg, "Analysis - Load Prior Information")
            else:
                # We grab the first prior and ignore the rest
                prior = priors[0]
                prior.source = 'file'
                prior.source_id = filename

                path, _ = os.path.split(filename)
                util_analysis_config.set_path(ini_name, path)
                
                # Now load the new prior into the block & GUI
                self._apply_selected_prior(prior, filename)


    def on_prior_ppm_start(self, event): 
        min, max = _paired_event(self.FloatPriorPpmStart,
                                 self.FloatPriorPpmEnd)
        self.block.set.prior_ppm_start = min
        self.block.set.prior_ppm_end   = max
        self.top.Freeze()
        wx.BeginBusyCursor()
        self.block.set.prior.calculate_full_basis_set(min, max, self.dataset )
        self.dynamic_metabolite_list.set_new_values()
        self.on_dynamic_metabolite_list()
        wx.EndBusyCursor()
        self.top.Thaw()


    def on_prior_ppm_end(self, event): 
        min, max = _paired_event(self.FloatPriorPpmStart,
                                 self.FloatPriorPpmEnd)
        self.block.set.prior_ppm_start = min
        self.block.set.prior_ppm_end   = max
        self.top.Freeze()
        wx.BeginBusyCursor()
        self.block.set.prior.calculate_full_basis_set(min, max, self.dataset )
        self.dynamic_metabolite_list.set_new_values()
        self.on_dynamic_metabolite_list()
        wx.EndBusyCursor()
        self.top.Thaw()

    def on_prior_calculate_combinations(self, event):
        val = event.GetEventObject().GetValue()
        self.block.set.prior_calculate_combinations = val
        

#    def on_prior_mask_editor(self, event): # wxGlade: PanelVoigtUI.<event_handler>
#        print "Event handler `on_prior_mask_editor' not implemented!"
#        event.Skip()
#
#    def on_prior_xrange_start(self, event): 
#        val_min = self.SpinPriorXrangeStart.GetValue()
#        val_max = self.SpinPriorXrangeEnd.GetValue()
#        param_min = min(val_min, val_max)
#        param_max = max(val_min, val_max)
#        self.block.set.prior_xrange = [param_min,param_max]
#        self.SpinPriorXrangeStart.SetValue(param_min)
#        self.SpinPriorXrangeEnd.SetValue(param_max)
#                      
#    def on_prior_xrange_end(self, event): 
#        val_min = self.SpinPriorXrangeStart.GetValue()
#        val_max = self.SpinPriorXrangeEnd.GetValue()
#        param_min = min(val_min, val_max)
#        param_max = max(val_min, val_max)
#        self.block.set.prior_xrange = [param_min,param_max]
#        self.SpinPriorXrangeStart.SetValue(param_min)
#        self.SpinPriorXrangeEnd.SetValue(param_max)
#        
#    def on_prior_yrange_start(self, event):
#        val_min = self.SpinPriorYrangeStart.GetValue()
#        val_max = self.SpinPriorYrangeEnd.GetValue()
#        param_min = min(val_min, val_max)
#        param_max = max(val_min, val_max)
#        self.block.set.prior_yrange = [param_min,param_max]
#        self.SpinPriorYrangeStart.SetValue(param_min)
#        self.SpinPriorYrangeEnd.SetValue(param_max)
#        
#    def on_prior_yrange_end(self, event):
#        val_min = self.SpinPriorYrangeStart.GetValue()
#        val_max = self.SpinPriorYrangeEnd.GetValue()
#        param_min = min(val_min, val_max)
#        param_max = max(val_min, val_max)
#        self.block.set.prior_yrange = [param_min,param_max]
#        self.SpinPriorYrangeStart.SetValue(param_min)
#        self.SpinPriorYrangeEnd.SetValue(param_max)
#
#    def on_prior_zrange_start(self, event):
#        val_min = self.SpinPriorZrangeStart.GetValue()
#        val_max = self.SpinPriorZrangeEnd.GetValue()
#        param_min = min(val_min, val_max)
#        param_max = max(val_min, val_max)
#        self.block.set.prior_zrange = [param_min,param_max]
#        self.SpinPriorZrangeStart.SetValue(param_min)
#        self.SpinPriorZrangeEnd.SetValue(param_max)
#
#    def on_prior_zrange_end(self, event): 
#        val_min = self.SpinPriorZrangeStart.GetValue()
#        val_max = self.SpinPriorZrangeEnd.GetValue()
#        param_min = min(val_min, val_max)
#        param_max = max(val_min, val_max)
#        self.block.set.prior_zrange = [param_min,param_max]
#        self.SpinPriorZrangeStart.SetValue(param_min)
#        self.SpinPriorZrangeEnd.SetValue(param_max)
#
#    def on_prior_mask_source(self, event): 
#        val = event.GetEventObject().GetSelection()
#        val = constants.FitPriorMaskSource.choices.byindex(val)[0]
#        self.block.set.prior_mask_source = val
#
#    def on_prior_ignore_mask(self, event): 
#        val = event.GetEventObject().GetValue()
#        self.block.set.prior_ignore_mask = val

    def on_b0_shift_method(self, event): 
        index = event.GetEventObject().GetSelection()
        value = list(constants.FitInitialB0ShiftMethod.choices.keys())[index]
        self.block.set.initial_b0_shift_method = value
        self.process_and_plot()

    def on_b0_shift_value(self, event): 
        val = event.GetEventObject().GetValue()
        voxel = self._tab_dataset.voxel
        orig = self.dataset.get_frequency_shift(voxel)
        # the method needs to be updated *before* set_freq_shift call
        # so we do not use the pubsub method here
#        key   = constants.FitInitialB0ShiftMethod.MANUAL
#        strng = constants.FitInitialB0ShiftMethod.choices[key]
#        self.ComboInitialB0ShiftMethod.SetStringSelection( strng )
#        self.block.set.initial_b0_shift_method = key
        self._tab_dataset.set_frequency_shift(val-orig, voxel)

    def on_initial_baseline_method(self, event): 
        index = event.GetEventObject().GetSelection()
        value = list(constants.FitInitialBaselineMethod.choices.keys())[index]
        self.block.set.initial_baseline_method = value
        self.PanelBaselineInitLowess.Show(bool(value))
        
        self.PanelInitialValues.Layout()
        self.PanelInitialValues.Refresh()
        self.process_and_plot()

    def on_initial_baseline_lowess_width(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.initial_baseline_lowess_width = val
        self.process_and_plot()

    def on_initial_baseline_ignore_width(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.initial_baseline_lowess_ignore_width = val
        self.process_and_plot()

    def on_initial_cr_cho_separation(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.initial_cr_cho_separation = val
        self.process_and_plot()

    def on_initial_peak_search_abs(self, event):
        val = event.GetEventObject().GetValue()
        self.block.set.initial_peak_search_abs = val
        self.process_and_plot()

    def on_initial_small_peak_areas(self, event):
        index = event.GetEventObject().GetSelection()
        value = list(constants.FitInitialSmallPeakAreas.choices.keys())[index]
        self.block.set.initial_small_peak_areas = value
        self.process_and_plot()

    def on_initial_small_peak_freqs(self, event):
        index = event.GetEventObject().GetSelection()
        value = list(constants.FitInitialSmallPeakFreqs.choices.keys())[index]
        self.block.set.initial_small_peak_freqs = value
        self.process_and_plot()

    def on_initial_linewidth_method(self, event): 
        index = event.GetEventObject().GetSelection()
        value = list(constants.FitInitialLinewidthMethod.choices.keys())[index]
        self.block.set.initial_linewidth_method = value
        show = (value != constants.FitInitialLinewidthMethod.MANUAL)
        self.PanelInitialLinewidth.Show(show)
        
        self.PanelInitialValues.Layout()
        self.PanelInitialValues.Refresh()
        self.process_and_plot()
        
    def on_linewidth_value(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.initial_linewidth_value = val
        key   = constants.FitInitialLinewidthMethod.MANUAL
        strng = constants.FitInitialLinewidthMethod.choices[key]
        self.ComboInitialLinewidthMethod.SetStringSelection( strng )
        self.block.set.initial_linewidth_method = key
        self.PanelInitialLinewidth.Show(False)
        self.PanelInitialValues.Layout()
        self.PanelInitialValues.Refresh()
        
        self.process_and_plot()

    def on_initial_linewidth_start(self, event): 
        min, max = _paired_event(self.FloatInitialLinewidthStart,
                                 self.FloatInitialLinewidthEnd)
        self.block.set.initial_linewidth_range_start = min
        self.block.set.initial_linewidth_range_end   = max
        self.process_and_plot()

    def on_initial_linewidth_end(self, event): 
        min, max = _paired_event(self.FloatInitialLinewidthStart,
                                 self.FloatInitialLinewidthEnd)
        self.block.set.initial_linewidth_range_start = min
        self.block.set.initial_linewidth_range_end   = max
        self.process_and_plot()

    def on_initial_linewidth_multiplier(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.initial_linewidth_fudge = val
        self.process_and_plot()
        
    def on_initial_phase_method(self, event): 
        index = event.GetEventObject().GetSelection()
        value = list(constants.FitInitialPhaseMethod.choices.keys())[index]
        self.block.set.initial_phase_method = value
        self.process_and_plot()

    def on_initial_phase0_value(self, event): 
        val = event.GetEventObject().GetValue()
        voxel = self._tab_dataset.voxel
        orig = self.dataset.get_phase_0(voxel)
        poll_label = [self._tab_dataset.indexAB[0]]
        self.top.notebook_datasets.global_poll_phase(poll_label, val-orig, voxel, do_zero=True)
        self.process_and_plot()
        
    def on_initial_phase1_value(self, event): 
        val = event.GetEventObject().GetValue()
        voxel = self._tab_dataset.voxel
        orig = self.dataset.get_phase_1(voxel)
        poll_label = [self._tab_dataset.indexAB[0]]
        self.top.notebook_datasets.global_poll_phase(poll_label, val-orig, voxel, do_zero=False)
        self.process_and_plot()

    def on_apply_ko_filter(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.initial_apply_ko_filter = val
        self.process_and_plot()
                
    def on_initial_ko_linewidth_minimum(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.initial_ko_linewidth_minimum = val
        self.process_and_plot()

    def on_initial_ko_points(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.initial_ko_points = val
        self.process_and_plot()

    def on_lineshape_model(self, event): 
        index = event.GetEventObject().GetSelection()
        value = list(FitLineshapeModel.choices.keys())[index]
        self.block.set.lineshape_model = value
        self.process_and_plot()
        

#    def on_mmol_method(self, event): 
#        val = event.GetEventObject().GetSelection()
#        val = constants.FitMacromoleculeMethod.choices.byindex(val)[0]
#        self.block.set.macromolecule_method = val
#        if val == 'None':
#            self.PanelMmolParameters.Hide()
#        else:
#            self.PanelMmolParameters.Show()
#
#    def on_mmol_reference_line(self, event): 
#        val = event.GetEventObject().GetStringSelection()
#        self.block.set.lineshape_model = val
#
#    def on_mmol_select_all(self, event):
#        self.dynamic_macromolecule_list.select_all()
#
#    def on_mmol_deselect_all(self, event):
#        self.dynamic_macromolecule_list.deselect_all()
#
#    def on_mmol_add(self, event): 
#        self.dynamic_auto_prior_list.add_row(constants.VOIGT_NEW_MACROMOLECULE_ROW, update=True)
#
#    def on_mmol_delete(self, event): 
#        self.dynamic_macromolecule_list.remove_checked_rows()
#
#    def on_mmol_restore_defaults(self, event): 
#        self.dynamic_macromolecule_list.reset_to_default()

    def on_baseline_method(self, event=None):
        # This handler is sometimes called programmatically, in which case
        # event is None, so don't assume it exists.
        index = self.ComboBaselineMethod.GetSelection()
        method = list(constants.FitBaselineMethod.choices.keys())[index]
        self.block.set.baseline_method = method
        
        if method:
            self.PanelBaselineParameters.Show()
            if method == constants.FitBaselineMethod.WAVELET_FILTER_BASIC:
                # # PyWavelets is an optional dependency, so we might have to
                # # tell the user that it's not installed.
                # if not constants.PYWAVELETS_AVAILABLE:
                #     answer = common_dialogs.message(_PYWAVELETS_UNAVAILABLE_MSG,
                #                                     "Baseline Method",
                #                                     common_dialogs.Q_YES_NO)
                #     if answer == wx.YES:
                #         webbrowser.open(_PYWAVELETS_EXPLANATION_URL, 1)
                    
                self.PanelWaveletParameters.Show()
                self.PanelBsplineParameters.Hide()
            else:
                self.PanelWaveletParameters.Hide()
                self.PanelBsplineParameters.Show()
                if method == constants.FitBaselineMethod.BSPLINE_FIXED_KNOT:
                    self.FloatBaselineSplineSpacing.Enable()
                    self.SpinBaselineSplineNknots.Disable()
                elif method == constants.FitBaselineMethod.BSPLINE_VARIABLE_KNOT:
                    self.FloatBaselineSplineSpacing.Disable()
                    self.SpinBaselineSplineNknots.Enable()
        else:
            self.PanelBaselineParameters.Hide()
        self.PanelBaseline.Layout()
        self.PanelBaseline.Refresh()


    def on_baseline_smoothing_flag(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.baseline_smoothing_flag = val

    def on_baseline_skip_last_smooth(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.baseline_skip_last_smooth = val

    def on_baseline_smoothing_width(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.baseline_smoothing_width = val

    def on_baseline_underestimate_method(self, event): 
        index = self.ComboBaselineUnderestimateMethod.GetSelection()
        method = list(constants.FitBaselineUnderestimateMethod.choices.keys())[index]
        self.block.set.baseline_underestimate_method = method

        if method == constants.FitBaselineUnderestimateMethod.FIRST_ITERATION_ONLY:
            self.LabelBaselineUnderestimateLast.Hide()
            self.FloatBaselineUnderestimateLast.Hide()
            self.LabelBaselineUnderestimateSteps.Hide()
            self.SpinBaselineUnderestimateSteps.Hide()
        else:
            self.LabelBaselineUnderestimateLast.Show()
            self.FloatBaselineUnderestimateLast.Show()
            self.LabelBaselineUnderestimateSteps.Show()
            self.SpinBaselineUnderestimateSteps.Show()
        self.PanelBaseline.Layout()
        self.PanelBaseline.Refresh()

    def on_baseline_underestimate(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.baseline_underestimate = val

    def on_baseline_underestimate_last(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.baseline_underestimate_last = val

    def on_baseline_underestimate_steps(self, event): 
        val = event.GetEventObject().GetValue()
        if val > self.block.set.optimize_global_iterations:
            val = self.block.set.optimize_global_iterations
            event.GetEventObject().SetValue(val)
        self.block.set.baseline_underestimate_steps = val

    def on_baseline_spline_nknots(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.baseline_spline_nknots = val

    def on_baseline_spline_spacing(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.baseline_spline_spacing = val
        self._set_spacing(val * self.dataset.spectral_hpp)

    def on_baseline_spline_order(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.baseline_spline_order = val

    def on_baseline_wavelet_scale(self, event): 
        val = event.GetEventObject().GetStringSelection()
        self.block.set.baseline_wavelet_scale = val

    def on_baseline_wavelet_min_dyad(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.baseline_wavelet_min_dyad = val

    def on_mmol_model(self, event):
        index = self.ComboMmolModel.GetSelection()
        method = list(constants.FitMacromoleculeMethod.choices.keys())[index]
        self.block.set.macromol_model = method

        if method:
            self.PanelMmolSingleDataset.Show()
            if method == constants.FitMacromoleculeMethod.SINGLE_BASIS_DATASET:
                self.PanelMmolSingleDataset.Show()
            else:
                self.PanelMmolSingleDataset.Hide()
        else:
            self.PanelMmolSingleDataset.Hide()
        self.PanelMacromol.Layout()
        self.PanelMacromol.Refresh()

    def on_mmol_single_dataset(self, event):
        # Allows the user to select a dataset with single macromolecule basis function
        dialog = dialog_dataset_browser.DialogDatasetBrowser(self.top.datasets, omit_self=self.dataset)
        dialog.ShowModal()
        mmol_dataset = dialog.dataset
        dialog.Destroy()
        
        if mmol_dataset:
            if mmol_dataset.id == self.dataset.id:
                # can not use self as a macromol basis function
                return
            self.dataset.set_associated_dataset_mmol(mmol_dataset)

        else:
            # allow user to have NO mmol basis set 
            self.block.set.macromol_single_basis_dataset       = None
            self.block.set.macromol_single_basis_dataset_id    = ''
            self.block.set.macromol_single_basis_dataset_fname = ''
        self.TextMmolSingleDataset.SetValue(self.block.set.macromol_single_basis_dataset_fname)

    def on_mmol_single_dataset_initval(self, event):
        index = self.ComboMmolSingleDatasetInitVal.GetSelection()
        method = list(constants.FitMacromoleculeMethodInitVal.choices.keys())[index]
        self.block.set.macromol_single_basis_dataset_initval = method
        self.set_mmol_single()      # enable/disable widgets
        self.PanelMacromol.Layout()
        self.PanelMacromol.Refresh()
                
    def on_mmol_single_dataset_start_area(self, event):
        val = event.GetEventObject().GetValue()
        self.block.set.macromol_single_basis_dataset_start_area = val
        self.process_and_plot()

    def on_mmol_single_dataset_start_fudge(self, event):
        val = event.GetEventObject().GetValue()
        self.block.set.macromol_single_basis_dataset_start_fudge = val
        self.process_and_plot()

    def on_mmol_single_dataset_ppm_start(self, event):
        min, max = _paired_event(self.FloatMmolSingleDatasetPpmStart, 
                                 self.FloatMmolSingleDatasetPpmEnd)
        self.block.set.macromol_single_basis_dataset_ppm_start = min
        self.block.set.macromol_single_basis_dataset_ppm_end   = max
        self.process_and_plot()

    def on_mmol_single_dataset_ppm_end(self, event):
        min, max = _paired_event(self.FloatMmolSingleDatasetPpmStart, 
                                 self.FloatMmolSingleDatasetPpmEnd)
        self.block.set.macromol_single_basis_dataset_ppm_start = min
        self.block.set.macromol_single_basis_dataset_ppm_end   = max
        self.process_and_plot()

    def on_mmol_single_dataset_limit_max(self, event):
        min, max = _paired_event(self.FloatMmolSingleDatasetLimitMin, 
                                 self.FloatMmolSingleDatasetLimitMax)
        self.block.set.macromol_single_basis_dataset_limit_min = min
        self.block.set.macromol_single_basis_dataset_limit_max = max

    def on_mmol_single_dataset_limit_min(self, event):
        min, max = _paired_event(self.FloatMmolSingleDatasetLimitMin, 
                                 self.FloatMmolSingleDatasetLimitMax)
        self.block.set.macromol_single_basis_dataset_limit_min = min
        self.block.set.macromol_single_basis_dataset_limit_max = max



    def on_optimize_method(self, event): 
        index = event.GetEventObject().GetSelection()
        value = list(constants.FitOptimizeMethod.choices.keys())[index]
        self.block.set.optimize_method = value

    def on_optimize_scaling_flag(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.optimize_scaling_flag = val

    def on_optimize_global_iterations(self, event):
        val = event.GetEventObject().GetValue()
        
        if self.block.set.baseline_underestimate_steps > val:
            # underestimation steps has to be LT global iterations and GT 2
            nstep = 2 if val < 2 else val
            self.block.set.baseline_underestimate_steps = nstep
            self.SpinBaselineUnderestimateSteps.SetValue(nstep)
        
        self.block.set.optimize_global_iterations = val

    def on_optimize_stop_tolerance(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.optimize_stop_tolerance = val

    def on_optimize_max_iterations(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.optimize_max_iterations = val

    def on_optimize_bounds_area_max(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.optimize_bounds_area_max = val

    def on_optimize_bounds_area_min(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.optimize_bounds_area_min = val

    def on_optimize_bounds_area_max_small(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.optimize_bounds_area_max_small = val

    def on_optimize_bounds_area_min_small(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.optimize_bounds_area_min_small = val

    def on_optimize_enable_bounds_area_small(self, event): 
        """ alternative area bounds for 'small' metabs  """
        val = event.GetEventObject().GetValue()
        self.block.set.optimize_enable_bounds_area_small = val

    def on_optimize_bounds_range_ppm(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.optimize_bounds_range_ppm = val

    def on_optimize_bounds_range_ppm_small(self, event): 
        """ alternative ppm bounds for 'small' metabs to enforce db values more closely """
        val = event.GetEventObject().GetValue()
        self.block.set.optimize_bounds_range_ppm_small = val

    def on_optimize_enable_bounds_ppm_small(self, event): 
        """ can turn this 'small' option on/off """
        val = event.GetEventObject().GetValue()
        self.block.set.optimize_enable_bounds_ppm_small = val

    def on_optimize_bounds_range_phase0(self, event):
        val = event.GetEventObject().GetValue()
        self.block.set.optimize_bounds_range_phase0 = val

    def on_optimize_bounds_range_phase1(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.optimize_bounds_range_phase1 = val

    def on_optimize_bounds_min_linewidth(self, event): 
        # careful here, there's an inverse relationship
        min, max = _paired_event( self.FloatOptimizeBoundsMinLinewidth,
                                  self.FloatOptimizeBoundsMaxLinewidth)
        self.block.set.optimize_bounds_max_linewidth = min
        self.block.set.optimize_bounds_min_linewidth = max
        max_lw_sec = self.block.set.optimize_bounds_max_linewidth
        min_lw_sec = self.block.set.optimize_bounds_min_linewidth
        min_lw_hz, _  = util_spectral.voigt_width(max_lw_sec, max_lw_sec, self.dataset)
        max_lw_hz, _  = util_spectral.voigt_width(min_lw_sec, min_lw_sec, self.dataset)
        self._set_min_linewidth(max_lw_hz)
        self._set_max_linewidth(min_lw_hz)

    def on_optimize_bounds_max_linewidth(self, event): 
        # careful here, there's an inverse relationship
        min, max = _paired_event(self.FloatOptimizeBoundsMinLinewidth,
                                 self.FloatOptimizeBoundsMaxLinewidth)
        self.block.set.optimize_bounds_max_linewidth = min
        self.block.set.optimize_bounds_min_linewidth = max
        max_lw_sec = self.block.set.optimize_bounds_max_linewidth
        min_lw_sec = self.block.set.optimize_bounds_min_linewidth
        min_lw_hz, _  = util_spectral.voigt_width(max_lw_sec, max_lw_sec, self.dataset)
        max_lw_hz, _  = util_spectral.voigt_width(min_lw_sec, min_lw_sec, self.dataset)
        self._set_min_linewidth(max_lw_hz)
        self._set_max_linewidth(min_lw_hz)

    def on_optimize_weights_method(self, event): 
        index = event.GetEventObject().GetSelection()
        value = list(constants.FitOptimizeWeightsMethod.choices.keys())[index]
        self.block.set.optimize_weights_method = value
        if self.block.set.optimize_weights_method == \
            constants.FitOptimizeWeightsMethod.LOCAL_WEIGHTING:
            self.PanelLocalWeights.Show()
            self.PanelOptimization.Layout()
        else:
            self.PanelLocalWeights.Hide()
            self.PanelOptimization.Layout()
        self.process_and_plot()

    def on_optimize_weights_scale_factor(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.optimize_weights_scale_factor = val
        self.process_and_plot()

    def on_optimize_weights_width_factor(self, event):
        val = event.GetEventObject().GetValue()
        self.block.set.optimize_weights_width_factor = val
        self.process_and_plot()

    def on_optimize_weights_water_flag(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.optimize_weights_water_flag = val
        self.process_and_plot()

    def on_optimize_weights_water_start(self, event):
        min, max = _paired_event(self.FloatOptimizeWeightsWaterStart,
                                 self.FloatOptimizeWeightsWaterEnd)
        self.block.set.optimize_weights_water_start = min
        self.block.set.optimize_weights_water_end = max
        self.process_and_plot()
        
    def on_optimize_weights_water_end(self, event):
        min, max = _paired_event(self.FloatOptimizeWeightsWaterStart,
                                 self.FloatOptimizeWeightsWaterEnd)
        self.block.set.optimize_weights_water_start = min
        self.block.set.optimize_weights_water_end = max
        self.process_and_plot()
        
    def on_optimize_weights_lipid_flag(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.optimize_weights_lipid_flag = val
        self.process_and_plot()
        
    def on_optimize_weights_lipid_start(self, event): 
        min, max = _paired_event(self.FloatOptimizeWeightsLipidStart,
                                 self.FloatOptimizeWeightsLipidEnd)
        self.block.set.optimize_weights_lipid_start = min
        self.block.set.optimize_weights_lipid_end = max
        self.process_and_plot()
        
    def on_optimize_weights_lipid_end(self, event): 
        min, max = _paired_event(self.FloatOptimizeWeightsLipidStart,
                                 self.FloatOptimizeWeightsLipidEnd)
        self.block.set.optimize_weights_lipid_start = min
        self.block.set.optimize_weights_lipid_end = max
        self.process_and_plot()
        
    def on_optimize_weights_small_peak_factor(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.optimize_weights_small_peak_factor = val
        self.process_and_plot()

    def on_optimize_constrain_ppm_naa_naag(self, event):
        val = event.GetEventObject().GetValue()
        unique_abbrs = self.dataset.prior_list_unique    # this is a property in dataset
        if 'naa' not in unique_abbrs or 'naag' not in unique_abbrs:
            val = False
        event.GetEventObject().SetValue(val)
        self.block.set.optimize_constrain_ppm_naa_naag = val

    def on_optimize_constrain_ppm_cr_pcr(self, event):
        val = event.GetEventObject().GetValue()
        unique_abbrs = self.dataset.prior_list_unique    # this is a property in dataset
        if 'cr' not in unique_abbrs or 'pcr' not in unique_abbrs:
            val = False
        event.GetEventObject().SetValue(val)
        self.block.set.optimize_constrain_ppm_cr_pcr = val

    def on_optimize_constrain_ppm_gpc_pcho(self, event):
        val = event.GetEventObject().GetValue()
        unique_abbrs = self.dataset.prior_list_unique    # this is a property in dataset
        if 'gpc' not in unique_abbrs or 'pcho' not in unique_abbrs:
            val = False
        event.GetEventObject().SetValue(val)
        self.block.set.optimize_constrain_ppm_gpc_pcho = val

    def on_optimize_constrain_ppm_cr2_pcr2(self, event):
        val = event.GetEventObject().GetValue()
        unique_abbrs = self.dataset.prior_list_unique    # this is a property in dataset
        if 'cr2' not in unique_abbrs or 'pcr2' not in unique_abbrs:
            val = False
        event.GetEventObject().SetValue(val)
        self.block.set.optimize_constrain_ppm_cr2_pcr2 = val

    def on_optimize_constrain_ppm_glu_gln(self, event):
        val = event.GetEventObject().GetValue()
        unique_abbrs = self.dataset.prior_list_unique    # this is a property in dataset
        if 'glu' not in unique_abbrs or 'gln' not in unique_abbrs:
            val = False
        event.GetEventObject().SetValue(val)
        self.block.set.optimize_constrain_ppm_glu_gln = val

    def on_optimize_constrain_ppm_tau_glc(self, event):
        val = event.GetEventObject().GetValue()
        unique_abbrs = self.dataset.prior_list_unique    # this is a property in dataset
        if 'tau' not in unique_abbrs or 'glc' not in unique_abbrs:
            val = False
        event.GetEventObject().SetValue(val)
        self.block.set.optimize_constrain_ppm_tau_glc = val

    def on_confidence_intervals_flag(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.confidence_intervals_flag = val

    def on_confidence_alpha(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.confidence_alpha = val

    def on_confidence_area_flag(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.confidence_area_flag = val

    def on_confidence_ppm_flag(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.confidence_ppm_flag = val

    def on_confidence_linewidth_flag(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.confidence_linewidth_flag = val

    def on_confidence_phase_flag(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.confidence_phase_flag = val

    def on_cramer_rao_flag(self, event): 
        val = event.GetEventObject().GetValue()
        self.block.set.cramer_rao_flag = val

    def on_cramer_rao_ppm_start(self, event): 
        min_, max_ = _paired_event(self.FloatCramerRaoPpmStart,
                                   self.FloatCramerRaoPpmEnd)
        self.block.set.cramer_rao_ppm_start = min_
        self.block.set.cramer_rao_ppm_end = max_

    def on_cramer_rao_ppm_end(self, event): 
        min_, max_ = _paired_event(self.FloatCramerRaoPpmStart,
                                   self.FloatCramerRaoPpmEnd)
        self.block.set.cramer_rao_ppm_start = min_
        self.block.set.cramer_rao_ppm_end = max_


    def on_output_results_to_html(self, event):
        ini_name = "voigt_output_as_html"
        default_path = util_analysis_config.get_path(ini_name)
        default_filename = "analysis_voigt_results.html"
        filetype_filter = "HTML (*.html)|*.html"
        
        filename = common_dialogs.save_as(default_path=default_path,
                                          filetype_filter=filetype_filter,
                                          default_filename=default_filename)
        if filename:
            # Build the image as a PNG.
            # FIXME PS - matplotlib's PNG support seems shaky under OS X so we 
            # use SVG instead. 
            # SVG works OK (all recent browsers support it, except IE < 9), 
            # but it's not ideal. Safari's rendering and resizing for SVGs is
            # not nearly as optimized as it is for PNGs. Plus it would be 
            # nice to behave the same on all platforms.
            if util_misc.get_platform() == "osx":
                format = "svg"
                mime_type = "image/svg+xml"
            else:
                format = "png"
                mime_type = "image/png"

            # Get matplotlib to write the image to our StringIO object
            fake_file = StringIO.StringIO()
        
            self.view.figure.savefig(fake_file, dpi=300, format=format)
            
            image_data = base64.b64encode(fake_file.getvalue())
        
            fake_file.close()


            voxel = self._tab_dataset.voxel
            raw   = self.dataset.blocks["raw"]

            data_source = raw.get_data_source(voxel)
            
            #html = self.block.results_as_html(voxel, self.plot_results['fitted_lw'],
            html = self.dataset.fit_results_as_html(voxel, self.plot_results['fitted_lw'],
                                                      self.plot_results['minmaxlw'][0],
                                                      self.plot_results['minmaxlw'][1], 
                                                      data_source,
                                                      (mime_type, image_data)
                                                     )

            open(filename, "wb").write(html)
            
            # We saved results, so we write the path to the INI file.
            path, _ = os.path.split(filename)
            util_analysis_config.set_path(ini_name, path)

        #else:
            # User clicked cancel on the "save as" dialog


    def on_output_results_to_csv(self, event):
        """ 
        Send current voxel fitting results to CSV file, save the 
        last successful CSV file as the default for the next call.
        
        """

        voxels = [self._tab_dataset.voxel,]
        
        ini_name = "voigt_output_as_csv"
        default_path_name = util_analysis_config.get_path(ini_name)
        default_path, default_fname = os.path.split(default_path_name)
        filetype_filter = "CSV (*.csv)|*.csv"
        
        filename = common_dialogs.save_as(default_path=default_path,
                                          filetype_filter=filetype_filter,
                                          default_filename=default_fname)
        if filename:

            try:
                self.output_results_to_csv(voxels, filename)

            except IOError:
                msg = """I can't write the file "%s".""" % filename
                common_dialogs.message(msg, style=common_dialogs.E_OK)
            else:
                # We saved results, so we write the path to the INI file.
                util_analysis_config.set_path(ini_name, filename)
            
        #else:
            # User clicked cancel on the "save as" dialog


    def output_results_to_csv(self, voxels, filename):
        """
        There are multiple events that can trigger sending one or all voxel
        fitting results to a text/csv file. Having this standalone method
        available allows me to ensure that each event triggers the same output
        method. 
        
        Inputs:
        
          voxels - list of lists with [x,y,z] voxel values to be saved
                    ex. [[2,0,0],] or [[1,0,0],[2,0,0]]
          
          filename - string, current path/filename to be used for CSV output
        
        """
            
        decor1 = self._prefs.csv_qa_metab_labels
            
        # Create output header and results strings, check element count. 
        # If the file exists, check that the element count is the same in 
        # in the last line as for this results line. If it is, just write
        # out the results string. If different length, output both the 
        # header and results strings.

        items = []
        raw   = self.dataset.blocks["raw"]


        for voxel in voxels:

            data_source = raw.get_data_source(voxel)
            dataset_filename = self.dataset.dataset_filename
            
            # FIXME - bjs, need to update current 'fitted_lw' and 'min/maxlw' 
            #         values if we save more than one voxel. For now, the 
            #         existing values for current voxel are used.
            
            #val, hdr = self.block.results_as_csv(voxel,
            val, hdr = self.dataset.fit_results_as_csv( voxel, 
                                                        self.plot_results['fitted_lw'],
                                                        self.plot_results['minmaxlw'][0],
                                                        self.plot_results['minmaxlw'][1], 
                                                        data_source,
                                                        dataset_filename,
                                                        decor1=decor1)
            nhdr = len(hdr)
            val = ",".join(val)
            hdr = ",".join(hdr)
            val += "\n"
            hdr += "\n"

            items.append([hdr,val])
             
        hdr_flag = True
        if os.path.isfile(filename):
            with open(filename, 'r+') as f:
                data = f.readlines()
                if len(data)>1:
                    last = data[-1]
                    nlast = len(last.split(','))
                    if nlast == nhdr:
                        hdr_flag = False

        with open(filename, 'a') as f:
            for item in items:
                hdr, val = item       
                if hdr_flag:
                    f.write(hdr)
                    hdr_flag = False
                f.write(val)
                
                
    def debug_results_to_csv(self, voxels, filename):
        """
        Used to debug fitting results by sending them to a text/csv file. Having 
        this standalone method available allows me to ensure that more than one  
        event can trigger this same output method. 
        
        Inputs:
        
          voxels - list of lists with [x,y,z] voxel values to be saved
                    ex. [[2,0,0],] or [[1,0,0],[2,0,0]]
          
          filename - string, current path/filename to be used for CSV output
        
        """
            
        decor1 = self._prefs.csv_qa_metab_labels
            
        # Create output header and results strings, check element count. 
        # If the file exists, check that the element count is the same in 
        # in the last line as for this results line. If it is, just write
        # out the results string. If different length, output both the 
        # header and results strings.

        items = []
        raw   = self.dataset.blocks["raw"]


        for voxel in voxels:

            data_source = raw.get_data_source(voxel)
            dataset_filename = self.dataset.dataset_filename
            
            # FIXME - bjs, need to update current 'fitted_lw' and 'min/maxlw' 
            #         values if we save more than one voxel. For now, the 
            #         existing values for current voxel are used.
            
            val, hdr = self.dataset.debug_results_as_csv(   voxel, 
                                                            data_source,
                                                            dataset_filename,
                                                            decor1=decor1)
            nhdr = len(hdr)
            val = ",".join(val)
            hdr = ",".join(hdr)
            val += "\n"
            hdr += "\n"

            items.append([hdr,val])
             
        hdr_flag = True
        if os.path.isfile(filename):
            with open(filename, 'r+') as f:
                data = f.readlines()
                if len(data)>1:
                    last = data[-1]
                    nlast = len(last.split(','))
                    if nlast == nhdr:
                        hdr_flag = False

        with open(filename, 'a') as f:
            for item in items:
                hdr, val = item       
                if hdr_flag:
                    f.write(hdr)
                    hdr_flag = False
                f.write(val)                
    

    def on_update_initial_values(self, event): 
        self.update_button_hit = True
        self.process_and_plot()
        
    def on_auto_initial_values(self, event):
        val = event.GetEventObject().GetValue()
        self.auto_update_initial_values = val
        self.process_and_plot()


    def on_fit_spectrum(self, event): 
        # check that dependent parameters are consistent
        
        # check that if mmol model is on, there is a dataset selected
        index = self.ComboMmolModel.GetSelection()
        method = list(constants.FitMacromoleculeMethod.choices.keys())[index]

        if (method == 'single_basis_dataset') and (self.block.set.macromol_single_basis_dataset is None):
            # inconsistent, reset widget and mmol parameters
            self.PanelMmolSingleDataset.Hide()
            self.PanelMacromol.Layout()
            self.PanelMacromol.Refresh()
            self.TextMmolSingleDataset.SetValue('')
            
            method = list(constants.FitMacromoleculeMethod.choices.keys())[0]   # None method
            self.block.set.macromol_model = method
            self.block.set.macromol_single_basis_dataset = None
            self.block.set.macromol_single_basis_dataset_id = ''
            self.block.set.macromol_single_basis_dataset_fname = ''
        
        self.process_and_plot(entry='full_fit')


    def on_batch_fit_all(self, event): 
        
        # check that if mmol model is on, there is a dataset selected
        index = self.ComboMmolModel.GetSelection()
        method = list(constants.FitMacromoleculeMethod.choices.keys())[index]

        if (method == 'single_basis_dataset') and (self.block.set.macromol_single_basis_dataset is None):
            # inconsistent, reset widget and mmol parameters
            self.PanelMmolSingleDataset.Hide()
            self.PanelMacromol.Layout()
            self.PanelMacromol.Refresh()
            self.TextMmolSingleDataset.SetValue('')
            
            method = list(constants.FitMacromoleculeMethod.choices.keys())[0]   # None method
            self.block.set.macromol_model = method
            self.block.set.macromol_single_basis_dataset       = None
            self.block.set.macromol_single_basis_dataset_id    = ''
            self.block.set.macromol_single_basis_dataset_fname = ''
        
        self._tab_dataset.batch_process_all(statusbar=self.top.statusbar)   # bjs - maybe could do direct call to self._dataset.batch_fit_all() ???
        
        self.process_and_plot(entry='voxel_change')
        self.view.canvas.draw()
                

    def on_dynamic_metabolite_list(self, event=None):
        # This is a fudged event called from the actual event that occurs 
        # inside the dynamic list class but its also invoked programmatically
        # (i.e. by our code, not by wx). In the latter case, event is None.

        # update Metabolite Basis set (subset) and start values
        vals = self.dynamic_metabolite_list.get_values()
        names, item2, item3, item4, item5, item6, item7 = vals
        self.block.set.prior_list           = names
        self.block.set.prior_area_scale     = item2
        self.block.set.prior_peak_ppm       = item3
        self.block.set.prior_search_ppm     = item4
        self.block.set.prior_db_ppm         = item5
        self.block.set.prior_fix_t2         = item6
        self.block.set.prior_search_ph0     = item7

        self.block.check_parameter_dimensions(self.dataset)

        # This is a little hack to work around a shortcoming in the 
        # scrolled panel implementation under OS X & Windows. When the panel
        # receives new content that causes the virtual size to exceed the
        # actual size, scrollbars don't appear until the window receives
        # a size event. This is a cheap way of forcing a size event -- make
        # the panel 1 pixel bigger or smaller. 
        width, height = self.PanelMetaboliteLines.GetSize()
        delta = (1 if height % 2 else -1)
        self.PanelMetaboliteLines.SetSize( (width, height + delta) )
        
        self.process_and_plot()
        


    #=======================================================
    #
    #           Public Methods
    #
    #=======================================================

    def on_voxel_change(self, voxel):
        """
        this just updates widgets that vary based on the voxel number
        selection. We do not update plot here because that is only done
        for the active tab in the inner notebook, and is called from 
        the _dataset.on_spin_update() routine 

        """
        self.FloatB0ShiftValue.SetValue(self.dataset.get_frequency_shift(voxel))
        self.FloatInitialPhase0Value.SetValue(self.dataset.get_phase_0(voxel))
        self.FloatInitialPhase1Value.SetValue(self.dataset.get_phase_1(voxel))
    

    def update_html_results_tab(self):
        """ 
        Update display of results in the html window in Results tab.
        
        At the moment we are only set up to have multiple files loaded into
        the first index, ie. SVS files loaded into the window, thus we use
        the voxel[0] index to select the filename.
        """
        voxel = self._tab_dataset.voxel
        raw   = self.dataset.blocks["raw"]

        #html = self.block.results_as_html(voxel, self.plot_results['fitted_lw'],
        html = self.dataset.fit_results_as_html(voxel, self.plot_results['fitted_lw'],
                                                  self.plot_results['minmaxlw'][0],
                                                  self.plot_results['minmaxlw'][1], 
                                                  raw.get_data_source(voxel))

        self.results_html_ctrl.SetPage(html)
        self.results_html_ctrl.SetBackgroundColour(self.GetBackgroundColour())        


    def process_and_plot(self, entry='initial_only'):
        """
        The process(), plot() and process_and_plot() methods are standard in
        all processing tabs. They are called to update the data in the plot
        results dictionary, the plot_panel in the View side of the tab or both.

        """    
        if self.auto_update_initial_values or self.update_button_hit or entry != 'initial_only':
            
            self.update_button_hit = False
            
            tab_base.Tab.process_and_plot(self, entry)
    
            if self._plotting_enabled:        
                self.process(entry=entry)
                self.plot()


    def process(self, entry='initial_only'):
        """
        Data processing results are stored into the Block inside the Chain,
        but the View results are returned as a dictionary from the Chain.run()
        method. The plot routine takes its inputs from this dictionary.

        Entry values available for this tab are (see chain_voigt module for 
        details of each entry point):
            "initial_only", "full_fit"        
        """
        tab_base.Tab.process(self, entry)

        if self._plotting_enabled: 
            dataset = self.dataset
            voxel   = self._tab_dataset.voxel
            
            self.plot_results = self.block.chain.run([voxel], 
                                                     entry=entry, 
                                                     statusbar=self.top.statusbar)
            
            self.update_html_results_tab()

            # this result is properly done, the block is updated in the chain, and
            # here we update the widgets in the tab
            self.FloatLinewidthValue.SetValue(self.block.set.initial_linewidth_value)
            
            self.FloatMmolSingleDatasetStartArea.SetValue(self.plot_results['mmol_area'])
            
            # new values for b0, ph0 and ph1 are set by the chain before returning, but
            # we need to update GUI here, so we call the 'tab_dataset' method that view
            # uses so much.  It want to know delta change to what is already set in the
            # block, so our hack here just a delta of 0 to refresh Spectral and Fit GUI
            # to what is already set in the block
            
            self._tab_dataset.set_frequency_shift(0.0, voxel, auto_calc=True)
            self._tab_dataset.set_phase_0(        0.0, voxel, auto_calc=True)
            self._tab_dataset.set_phase_1(        0.0, voxel, auto_calc=True)
            

    def plot(self):
        tab_base.Tab.plot(self)

        if self._plotting_enabled: 
            dataset = self.dataset
            tdfd    = dataset.blocks["spectral"]

            if not tdfd.set.fft: 
                return

            voxel = self._tab_dataset.voxel
            ph0   = self.dataset.get_phase_0(voxel)
            ph1   = self.dataset.get_phase_1(voxel)

            # Set up aliases and do helpful calculations
            results = self.plot_results
            dim0    = dataset.spectral_dims[0]
            nmet    = len(self.block.set.prior_list)

            # take max for use in scaling wtarr
            ymax = max(np.abs(results['data']))

            # update the wtarr array in overlay1 attribute init display to False
            self.view.set_overlay1(results['weight_array'].copy()*ymax*0.5)  
            self.view.set_show_overlay1(False)
            
            # loop through the 4 axes and set data in each as needed
            data = []
            for i in range(4):

                # Get copies of results arrays. For the spectra that are 1D 
                # arrays, I have edited their shape so as to be able to 
                # concatenate them as needed in the plot types below
                freq  = results['data'].copy()
                base  = results['fit_baseline'].copy()
                base.shape = 1, base.shape[0]
                yfit  = results['yfit'].copy()
                if len(yfit.shape) == 1:
                    yfit.shape = 1,yfit.shape[0]
                yini  = results['yini'].copy()
                if len(yini.shape) == 1:
                    yini.shape = 1,yini.shape[0]
                yinib = results['init_baseline'].copy()
                yinib.shape = 1, yinib.shape[0]
                yfits = results['yfit'].copy() if len(yfit.shape)==1 else np.sum(yfit, axis=0)
                
                line1     = np.zeros(dim0, complex)
                line2     = np.zeros(dim0, complex)
                color1    = self._prefs.line_color_raw
                color2    = self._prefs.line_color_raw

                if self._prefs.plotx[i].raw_data:           
                    line1  = freq
                    line1 = {'data' : line1, 
                             'line_color_real'      : self._prefs.line_color_real,
                             'line_color_imaginary' : self._prefs.line_color_imaginary,
                             'line_color_magnitude' : self._prefs.line_color_magnitude }
                elif self._prefs.plotx[i].fitted_data:
                    line1  = yfit
                    color1 = self._prefs.line_color_fit
                elif self._prefs.plotx[i].baseline:
                    line1  = base
                    color1 = self._prefs.line_color_base
                elif self._prefs.plotx[i].combo_raw_minus_fit:
                    line1  = freq - yfits
                    line1 = {'data' : line1, 
                             'line_color_real'      : self._prefs.line_color_real,
                             'line_color_imaginary' : self._prefs.line_color_imaginary,
                             'line_color_magnitude' : self._prefs.line_color_magnitude }
                elif self._prefs.plotx[i].combo_raw_minus_base:
                    line1  = freq - base
                    line1 = {'data' : line1, 
                             'line_color_real'      : self._prefs.line_color_real,
                             'line_color_imaginary' : self._prefs.line_color_imaginary,
                             'line_color_magnitude' : self._prefs.line_color_magnitude }
                elif self._prefs.plotx[i].combo_raw_minus_fit_minus_base:
                    line1  = freq-base-yfits
                    line1 = {'data' : line1, 
                             'line_color_real'      : self._prefs.line_color_real,
                             'line_color_imaginary' : self._prefs.line_color_imaginary,
                             'line_color_magnitude' : self._prefs.line_color_magnitude }
                elif self._prefs.plotx[i].combo_raw_and_fit:
                    line1  = freq
                    line1 = {'data' : line1, 
                             'line_color_real'      : self._prefs.line_color_real,
                             'line_color_imaginary' : self._prefs.line_color_imaginary,
                             'line_color_magnitude' : self._prefs.line_color_magnitude }
                    line2  = yfit
                    color2 = self._prefs.line_color_fit
                elif self._prefs.plotx[i].combo_raw_and_base:
                    line1  = freq
                    line1 = {'data' : line1, 
                             'line_color_real'      : self._prefs.line_color_real,
                             'line_color_imaginary' : self._prefs.line_color_imaginary,
                             'line_color_magnitude' : self._prefs.line_color_magnitude }
                    line2  = base
                    color2 = self._prefs.line_color_base
                elif self._prefs.plotx[i].combo_raw_and_fit_plus_base:
                    line1  = freq
                    line1 = {'data' : line1, 
                             'line_color_real'      : self._prefs.line_color_real,
                             'line_color_imaginary' : self._prefs.line_color_imaginary,
                             'line_color_magnitude' : self._prefs.line_color_magnitude }
                    line2  = np.concatenate((yfit,base))
                    color2 = self._prefs.line_color_fit
                elif self._prefs.plotx[i].combo_raw_minus_base_and_fit:
                    line1  = (freq-base)
                    line1 = {'data' : line1, 
                             'line_color_real'      : self._prefs.line_color_real,
                             'line_color_imaginary' : self._prefs.line_color_imaginary,
                             'line_color_magnitude' : self._prefs.line_color_magnitude }
                    line2  = yfit
                    color2 = self._prefs.line_color_fit
                elif self._prefs.plotx[i].combo_fit_plus_base:
                    line1  = yfits+base
                    color1 = self._prefs.line_color_fit
                elif self._prefs.plotx[i].combo_fit_and_base:
                    line1  = yfit
                    line2  = base
                    color1 = self._prefs.line_color_fit
                    color2 = self._prefs.line_color_base
                elif self._prefs.plotx[i].combo_raw_and_init_model:
                    line1  = freq
                    line1 = {'data' : line1, 
                             'line_color_real'      : self._prefs.line_color_real,
                             'line_color_imaginary' : self._prefs.line_color_imaginary,
                             'line_color_magnitude' : self._prefs.line_color_magnitude }
                    line2  = np.concatenate((yini, yinib))
                    color2 = self._prefs.line_color_init
                elif self._prefs.plotx[i].combo_raw_and_wt_arr:
                    line1  = freq
                    line1 = {'data' : line1, 
                             'line_color_real'      : self._prefs.line_color_real,
                             'line_color_imaginary' : self._prefs.line_color_imaginary,
                             'line_color_magnitude' : self._prefs.line_color_magnitude }
                    self.view.set_show_overlay1(True, index=i)

                if isinstance(line1, dict):
                    dat1 = line1
                else:
                    dat1 = {'data' : line1, 
                            'line_color_real'      : color1,
                            'line_color_imaginary' : color1,
                            'line_color_magnitude' : color1 }

                if isinstance(line2, dict):
                    dat2 = line2
                else:
                    dat2 = {'data' : line2, 
                            'line_color_real'      : color2,
                            'line_color_imaginary' : color2,
                            'line_color_magnitude' : color2 }
                
                data.append([dat1,dat2])

            self.view.set_data(data)
            self.view.update(no_draw=True, set_scale=not self._scale_intialized)
            if not self._scale_intialized:
                self._scale_intialized = True

            self.view.set_phase_0(ph0, absolute=True, no_draw=True)
            self.view.set_phase_1(ph1, absolute=True, no_draw=True)
            
            self.view.canvas.draw()

            # Calculate the new area after phasing
            area, rms = self.view.calculate_area()
            if self._prefs.area_calc_plot_a:
                index = 0
            elif  self._prefs.area_calc_plot_b:
                index = 1
            elif  self._prefs.area_calc_plot_c:
                index = 2
            elif  self._prefs.area_calc_plot_d:
                index = 3
            self.top.statusbar.SetStatusText(self.build_area_text(area[index], rms[index]), 3)


    #=======================================================
    #
    #           Internal Helper Functions 
    #
    #=======================================================

    def _apply_selected_prior(self, prior, source):
        """
        Given a newly-selected prior & source info, updates the prior and GUI. 
        Meant to be called after user loads prior info from the database or a file.
        
        """
        self.block.set.prior = prior
        prior.calculate_full_basis_set(self.block.set.prior_ppm_start,
                                       self.block.set.prior_ppm_end,
                                       self.dataset)

        self.TextPriorInformationSource.SetValue(source)

        self.top.Freeze()
        self.dynamic_metabolite_list.set_new_values()
        
        # Force the scrolled panel to adjust to its new conent.
        self.PanelMetaboliteLines.Layout()
        
        self.top.Thaw()
        self.on_dynamic_metabolite_list()


    def _set_max_linewidth(self, value):
        # A convenience/consistency function for populating the 
        # max linewdith control.
        self.StaticMaxLinewidth.SetLabel("%.3f" % value)


    def _set_min_linewidth(self, value):
        # A convenience/consistency function for populating the 
        # min linewdith control.
        self.StaticMinLinewidth.SetLabel("%.3f" % value)


    def _set_spacing(self, value):
        # A convenience/consistency function for populating the 
        # baseline/b-spline spacing label.
        self.StaticHzSpacing.SetLabel("%.3f" % value)





   

################################################################################
###### Dynamic Macromolecule List Class ########################################
#
#class DynamicMacromoleculeList(object):
#
#    def __init__(self, _inner_notebook, GridSizer, generic_basis, dataset):
#        
#        self._tab_dataset        = _inner_notebook
#        self.generic_basis = generic_basis
#        self.dataset       = dataset
#        self.block         = dataset.blocks["fit"]
#        
#        # We follow the wx CamelCaps naming convention for this wx object.
#        self.GridSizer = GridSizer
#                
#        self.list_lines = []
#
#        self.maxppm = self.dataset.pts2ppm(0)
#        self.minppm = self.dataset.pts2ppm(self._dataset.dims[0]-1)
#
#        self.default = self.generic_basis.get_rows()
#
#        for row_vals in default:
#            self.add_row(row_vals)
#
#    
#    def add_row(self, row_vals, update=False):
#        """
#        Adds a row to the end of the list. 
#        """
#
#        # create widgets to go into the line
#        list_line = { }
#
#        checkbox = wx.CheckBox(self._tab_dataset)
#        
#        value_ppm   = FloatSpin(self._tab_dataset, agwStyle=FS_LEFT)
#        value_area  = FloatSpin(self._tab_dataset, agwStyle=FS_LEFT)
#        value_lwhz  = FloatSpin(self._tab_dataset, agwStyle=FS_LEFT)
#        limit_ppm  = FloatSpin(self._tab_dataset, agwStyle=FS_LEFT)
#        limit_area = FloatSpin(self._tab_dataset, agwStyle=FS_LEFT)
#        limit_lwhz = FloatSpin(self._tab_dataset, agwStyle=FS_LEFT)
#        
#        # keep a copy of panel and widgets to access later
#        line = { "checkbox"     : checkbox, 
#                 "value_ppm"    : ppm, 
#                 "value_area"   : area, 
#                 "value_lwhz"   : lwhz,
#                 "ppm_range"    : limit_ppm, 
#                 "area_range"   : limit_area, 
#                 "lwhz_range"   : limit_lwhz,
#               }
#
#        # Add the controls to the grid sizer
#        self.GridSizer.Add(line["checkbox"], 0, wx.ALIGN_CENTER_VERTICAL)
#        for key in ("value_ppm", "limit_ppm", "value_area", "limit_area", "value_lwhz", "limit_lwhz"):
#            self.GridSizer.Add(line[key], 0, wx.EXPAND)
#
#        # Configure the controls I just created
#
#        # All of the floatspins have the same size. 
#        floatspin_size = wx.Size(70, -1)
#
#        # Note. On these Spin and FloatSpin widgets, if the value you want to
#        #    set is outside the wxGlade standard range, you should make the 
#        #    call to reset the range first and then set the value you want.
#
#        wx_util.configure_spin(value_ppm,  70, 2, 0.25,(self.minppm,self.maxppm))
#        wx_util.configure_spin(value_area, 70, 3, 0.25,(0.001,100000.0))
#        wx_util.configure_spin(value_lwhz, 70, 2, 1.0, (0.001,10000.0))
#        wx_util.configure_spin(limit_ppm,  70, 2, 1.0, (self.minppm,self.maxppm))
#        wx_util.configure_spin(limit_area, 70, 3, 0.5, (0.001,100000.0))
#        wx_util.configure_spin(limit_lwhz, 70, 2, 5.0, (0.001,10000.0))
#
#        self.list_lines.append(line)
#        
#        self._tab_dataset.Layout()
#        self._tab_dataset.Fit()
#        
#        if update:
#            self.update_generic_basis()
#        
#        
#    def remove_checked_rows(self):
#        # gather indices of all checked boxes
#        checklist = []
#        
#        for i, line in enumerate(self.list_lines):
#            if line["checkbox"].GetValue():
#                checklist.append(i)
#        
#        # remove in reverse order so we don't invalidate later
#        # indices by removing the ones preceding them in the list
#        checklist.reverse()
#        
#        for i in checklist:
#            # Each line is a dict of controls + mixture info
#            for item in self.list_lines[i].values():
#                if hasattr(item, "Destroy"):
#                    # It's a wx control
#                    item.Destroy()
#                
#            del self.list_lines[i]
#            
#        # Reduce the # of rows in the grid sizer
#        rows = self.GridSizer.GetRows()
#        self.GridSizer.SetRows(rows - len(checklist))
#        self.GridSizer.Layout()
#        self.update_generic_basis()
#        
#    def reset_to_default(self):
#        if not self.default:
#            return
#        for line in self.list_lines:
#            for item in line.values():
#                if hasattr(item, "Destroy"):
#                    # It's a wx control
#                    item.Destroy()
#        self.list_lines = []
#        for row_vals in default:
#            self.add_row(row_vals)
#        self.update_generic_basis()
#        
#    def get_values(self):
#        return [self.get_line_values(line) for line in self.list_lines]
#
#    def get_line_values(self, line):
#        # Returns a dict containing  values of the controls in the line.
#        return { "checkbox"    : line["checkbox"].GetValue(),
#                 "value_ppm"   : line["value_ppm"].GetValue(),
#                 "value_area"  : line["value_area"].GetValue(),
#                 "value_phase" : 0.0,
#                 "value_lwhz"  : line["value_lwhz"].GetValue(),
#                 "limit_ppm"   : line["limit_ppm"].GetValue(),
#                 "limit_area"  : line["limit_area"].GetValue(),
#                 "limit_phase" : 0.0,
#                 "limit_lwhz"  : line["limit_lwhz"].GetValue()
#               }
#
#    def select_all(self):
#        for line in self.list_lines:
#            line["checkbox"].SetValue(True)
#
#    def deselect_all(self):
#        for line in self.list_lines:
#            line["checkbox"].SetValue(False)
#
#    def update_generic_basis(self):
#        vals = self.get_values()
#        self.generic_basis.set_values(vals)

#     
# 
#     
#     def kiss_off_filter(self):
#         """
#         Estimates macromolecular peaks by subtracting the data after a certain number
#         of kiss off points from the original data.
# 
#         """
#         data = self.linktinfo.fitt_data
#         time = self.linktinfo.fitt_time
#         freq = np.fft.fft(time) / len(time)
#         yini = np.sum(self.linktinfo.fitt_yini, axis=1)
#         import pylab
#         
#         # Creates freq with kiss off points
#         nkopts = data.macromolecule_kiss_off_points
#         npts = data.dims[0]
#         temp = time.copy()
#         temp[0:npts-nkopts] = temp[nkopts:npts]
#         temp[npts-nkopts:npts] = 0
#         if nkopts % 2:
#             freq_kopts = np.fft.fft(-temp) / len(temp)
#         else:
#             freq_kopts = np.fft.fft(temp) / len(temp)
# 
#         # Do phasing for each spectrum
#         pivot  = data.ppm2pts(data.tdfd.set.phase_1_pivot)
#         phase0 = self.dataset.get_phase_0()*numpy.pi/180
#         phase1 = self.dataset.get_phase_1()*numpy.pi/180 * (numpy.arange(npts)-pivot)/npts
#         phase  = numpy.exp(1j*(phase0+phase1))
#         freq  *= phase
#         freq_kopts *= phase
# 
#         # Subtract kiss off point spectrum from initial spectrum and do smoothing
#         temp = freq - abs(freq_kopts)
#         base = numpy.zeros(npts)
#         ptran = 5
#         for i in range(ptran,npts-ptran):
#             base[i] = numpy.mean(temp[i-ptran:i+ptran])
#         base = numpy.where(base > 0, base, 0)
#         pylab.plot(base, 'r')
#         
#         ppm = []
#         amp = []
#         lw  = []
#         basis = []
#         max_base = max(base)
#         count = 0
#         
#         # Estimate the amplitude, location, and linewidth of the tallest peak
#         # in base.  Subtract the peak from base until the residual spectrum is
#         # below 20% of the original base.
#         while max(base) > 0.2 * max_base:
#             output = util_voigt.find_mm_peaks(data, base)
#             if output[3] > data.macromolecule_linewidth_minimum:
#                 ppm.append(output[1])
#                 amp.append(output[2])
#                 lw.append(output[3])
#             pylab.plot(output[0], 'k')
#             base -= output[0]
#             base = numpy.where(base > 0, base, 0)
#             count += 1
#             if count > 10:
#                 print "Error: Maximum iterations"
#                 break
# 
#         macromolecules_in_model = len(ppm)
# 
#         # Update the macromolecule number and widget
#         data.macromolecules_in_model = macromolecules_in_model
#         self.spinModelLines.SetValue(macromolecules_in_model)
#         self.update_macromolecule_number()
# 
#         # Set the area, ppm, and linewidth for the calculated macromolecules
#         for i in range(macromolecules_in_model):
#             data.macromolecule_peak_ppm[i] = ppm[i]
#             data.macromolecule_fit_flags[i] = True
#             data.macromolecule_values_ppm[i] = ppm[i]
#             data.macromolecule_values_area[i] = amp[i]
#             data.macromolecule_values_linewidth[i] = lw[i]
#             data.macromolecule_limits_ppm[i] = 0.10
#             data.macromolecule_limits_area[i] = amp[i]/2.0
#             data.macromolecule_limits_linewidth[i] = lw[i]/2.0
#             self.checkboxMacromolecules[i].SetValue(True)
#             self.gridMacromolecules.SetCellValue(i, 0, '%.2f'%(ppm[i]))
#             self.gridMacromolecules.SetCellValue(i, 1, '0.10')
#             self.gridMacromolecules.SetCellValue(i, 2, '%.2f'%(amp[i]))
#             self.gridMacromolecules.SetCellValue(i, 3, '%.2f'%(amp[i]/2.0))
#             self.gridMacromolecules.SetCellValue(i, 4, '%.2f'%(lw[i]))
#             self.gridMacromolecules.SetCellValue(i, 5, '%.2f'%(lw[i]/2.0))
# 
#         # Recalculate the basis
#         for i in range(macromolecules_in_model):
#             data.prior_mmbasis[:,i], data.macromolecule_peak_ppm[i] = \
#                 util_voigt.make_prior_basis([1], 
#                                               data.macromolecule_values_ppm[i], 
#                                               0.0, data)
# 
#         pylab.plot(base, 'b')
#         pylab.plot(freq, 'g')
# #        pylab.show()
         


# def find_mm_peaks(data, base):
#     """
#     =========
#     Arguments
#     =========
#     **data:**  [list][float] fitting data structure that contains all fitting parameters
#     **base:**  [list][float] input spectrum
# 
#     ===========
#     Description
#     ===========
#     Finds the location of macromolecular peaks by attempting to fit a
#     Lorentzian lineshape to the maximum amplitude of the input function.
#     It estimates the location and amplitude of the peak.
# 
#     ======
#     Syntax
#     ======
#     ::
# 
#       spectrum, ppm, amp, linewidth = find_mm_peaks( data, base )
# 
#     """
#     # First, find the ppm of the tallest peak and create initial basis
#     pts = np.where(base == max(base))[0][0]
# 
#     i = pts
#     while base[i] > 0.9 * base[pts]:
#         i -= 1
#     left = i
#     i = pts
#     while base[i] > 0.9 * base[pts]:
#         i += 1
#     right = i
#     pts = int((left+right)/2)
# 
#     # Create a basis function for a peak at the estimated location
#     ppm = data.pts2ppm(pts, )
#     basis = np.zeros(data.dims[0], complex)
#     basis[0:data.time.dims[0]] = calculate_generic_basis_line([1], ppm, 0.0, data, nopkppm=True)
# 
#     linewidth = 10
#     tmm = 1.0 / data.time.sw * np.arange(data.dims[0])
#     tamm  = 1.0 / (linewidth * 0.34 * np.pi)
#     lshapemm = np.exp(-(tmm/tamm)**2)
# 
#     # Create a spectrum from the estimated basis function with amplitude adjustment
#     spectrum = basis * lshapemm
#     spectrum = np.fft.fft(spectrum) / len(spectrum)
#     spectrum = spectrum.real
#     spectrum -= min(spectrum)
#     spectrum *= max(base) / max(spectrum)
# 
#     # Check is minimized to find the best fit
#     check = np.sum(abs(base-spectrum))
# 
#     # Gradually increase the linewidth by 1 Hz up to a maximum of 110 Hz
#     for i in range(100):
#         linewidth += 1
# 
#         # Recalculate the spectrum with new linewidth value
#         tamm  = 1.0 / (linewidth * 0.34 * np.pi)
#         lshapemm = np.exp(-(tmm/tamm)**2)
#         spectrum = basis * lshapemm
#         spectrum = np.fft.fft(spectrum) / len(spectrum)
#         spectrum = spectrum.real
#         spectrum -= min(spectrum)
# 
#         # Adjust amplitude to maximum of input function
#         amp = max(base) / max(spectrum)
#         spectrum *= amp
# 
#         # Recalculate difference between input function and estimated spectrum
#         # If the difference is greater than the previous linewidth value,
#         # end the loop and return the results.
#         tempcheck = np.sum(abs(base-spectrum))
#         if tempcheck > check:
#             break
#         check = tempcheck
# 
#     return spectrum, ppm, amp, linewidth


# FIXME PS - August 2012 - lorgauss_mm() isn't called from anywhere; it
# should be removed.


# def lorgauss_mm(a, fitt_data, pderflg = False, nobase = False, indiv = False, finalw = False):
#     """
#     =========
#     Arguments
#     =========
#     **a:**  [array][float]

#       xxxx

#     **fitt_data:**  [object][Dataset]

#       xxxx

#     **pderflg:**  [keyword][bool]

#       xxxx

#     **nobase:**  [keyword][bool]

#       xxxx

#     **indiv:**  [keyword][bool]

#       xxxx

#     **finalw:**  [keyword][bool]

#       xxxx

#     ===========
#     Description
#     ===========
#     Returns the parameterized metabolite model function. Includes a parameterized
#     macromolecular baseline function as part of the model.

#     A contains either: [[am],[fr],[Ta],[Tb],ph0,ph1]          - LorGauss real/complex
#                        [[am],[fr],[Ta],[Tb] ]                 - LorGauss magnitude
#                        [[am],[fr],[Ta],[Tb],ph0,ph1,[coef]]   - Spline Lineshape

#     metabolite basis function ampls and freqs are taken from the prior info ptr
#     in doof, so the values in [am] and [fr] are relative multipliers
#     and additives respectively.  That is why there is only
#     one value for each compound in each array

#     If the relfreq flag is ON, then [fr] is a single value that
#     is added to each peak freq equivalently.  Ie. the whole
#     spectrum can shift, but relative ppm separations between
#     all metabolites are maintained exactly.  If the flag is OFF,
#     then metabs may shift independently from one another,
#     however, within groups of peaks belonging to the same
#     metabolite, relative ppm separtaions are maintained.

#     Ta and Tb in this function are now arrays. It is possible for each metabolite
#     or macromolecule line to have its own linewidth now. But, in practice the
#     metabolites will all have the same parameters while MM lines will have separate
#     ones.  So the first value in the Ta and Tb arrays will be distributed across
#     all metabolites and the remainder will be used for MM lines.

#     am    - peak amplitude
#     fr    - peak frequency offsets in PPM
#     Ta    - T2 decay constant in sec
#     Tb    - T2 star decay const in sec
#     ph0/1 - zero/first order phase in degrees

#     coef  - are the spline coefs for the lineshape, knot locations are in doof
#     ======
#     Syntax
#     ======
#     ::

#       yfit = lorgauss_mm(a, fitt_data,  pderflg= False,
#                                         nobase = False,
#                                         indiv  = False,
#                                         finalw = False)

#     """

#     nmm = fitt_data.macromolecule_number
#     if fitt_data.macromolecule_number <= 0:
#         return lorgauss(a, fitt_data, pderflg, nobase, indiv, finalw)

#     # Set some constants and flags ---
#     DDTOR   = np.pi / 180.0
#     nmet    = fitt_data.metabolite_number
#     npts    = fitt_data.raw.dims[0]
#     nptszf  = round(npts * dataset.tdfd.set.zero_fill_multiplier)
#     sw      = 1.0*fitt_data.raw.sw
#     td      = 1.0/sw
#     piv     = dataset.ppm2pts(dataset.phase_1_pivot, acq=True)
#     arr1    = np.zeros(npts,float) + 1.0
#     f       = np.zeros((nptszf,nmet+nmm),complex)

#     range2d = np.zeros((npts,nmet))
#     num     = np.arange(npts*nmet)
#     for i in range(nmet):
#         range2d[:,i] = num[i*npts:(i+1)*npts]
#     t       = (range2d % npts) * td

#     range2d = np.zeros((npts,nmm))
#     num     = np.arange(npts*nmm)
#     for i in range(nmm):
#         range2d[:,i] = num[i*npts:(i+1)*npts]
#     tmm     = (range2d % npts) * td  # for T2, always decreases with time

#     # get prior max peak ppm vals for metabs which are flagged ON ---
#     indx    = np.where(fitt_data.metabolite_fit_flags)[0]
#     cnt     = np.size(indx)
#     mpeaks  = dataset.tdfd.prior_peak_ppm[indx]

#     # Decrypt a - depending on which MM model we use ---
#     am   = a[0:nmet]
#     fr   = a[nmet:nmet*2]
#     ta   = a[nmet*2]
#     tb   = a[nmet*2+1]
#     ph0  = a[nmet*2+2]
#     ph1  = a[nmet*2+3]
#     na   = nmet*2+4

#     ph0mm   = ph0
#     mmflag  = np.where(fitt_data.macromolecule_fit_flags)[0]
#     mmpkamp = fitt_data.macromolecule_values_area[mmflag]
#     mmpkppm = fitt_data.macromolecule_peak_ppm[mmflag]
#     mmpkwid = fitt_data.macromolecule_values_linewidth[mmflag]

#     if fitt_data.macromolecule_method_peaks == constants.FittingMacromoleculePeak.GROUPED_PEAKS:
#         # Voigt & grouped MM area/freq
#         ammm = np.array([a[na]]*nmm) * mmpkamp[0]
#         frmm = a[na+1] - fitt_data.ppm2hz(mmpkppm[0])*2.0*np.pi
#         frmm = np.zeros(nmm)+frmm

#         if fitt_data.macromolecule_method_linewidths == constants.FittingMacromoleculeLinewidths.LUMPED:
#             tamm0 = 1.0/(ta * 0.34 * np.pi)    # same as for metab, but convert to Hz
#         elif fitt_data.macromolecule_method_linewidths == constants.FittingMacromoleculeLinewidths.INDEPENDENT:
#             tamm0 = a[na+2]                # indep of metab - but only 1 basis fn

#     elif fitt_data.macromolecule_method_peaks == constants.FittingMacromoleculePeak.INDIVIDUAL_PEAKS:
#         # Voigt & individual MM area/freq
#         ammm = a[na:na+nmm] * mmpkamp
#         frmm = a[na+nmm:na+nmm*2] - fitt_data.ppm2hz(mmpkppm)*2.0*np.pi

#         if fitt_data.macromolecule_method_linewidths == constants.FittingMacromoleculeLinewidths.LUMPED:
#             tamm0 = a[na+nmm*2]            # indep of metab - and only 1 lw param
#         elif fitt_data.macromolecule_method_linewidths == constants.FittingMacromoleculeLinewidths.INDEPENDENT:
#             tamm0 = a[na+nmm*2:na+nmm*3]    # indep of metab - and multiple lw param

#     tamm0 = tamm0 + mmpkwid    # add to baseline widths - bumps tamm to number of elements of nmm
#     tamm  = 1.0 / (tamm0 * 0.34 * np.pi)    # convert the +/- Hz range to Tb decay

#     # Set up Lineshape - first for metab, then for mmolec and combine ---
#     expo     = t/ta + (t/tb)**2
#     lshape   = np.exp(-expo)
#     tamm     = np.outer(arr1,tamm) # outer is matrix multiplication
#     lshapemm = np.exp(-(tmm/tamm)**2)

#     # if FID, then for correct area, first pt must be div by 2.0 ---
#     # apply params to metab basis
#     if nmet:
#         tmp = fitt_data.basis_mets.copy()
#         fre = fr - fitt_data.ppm2hz(mpeaks)*2.0*np.pi    # in Radians here
#         fre = np.exp( 1j * np.outer(arr1,fre) * t) # outer is matrix multiplication
#         amp = np.outer(arr1, am)
#         ph0 = np.outer(arr1, np.exp(1j * (np.zeros(nmet)+ph0)))
#         tmp *= amp * fre * ph0 * lshape

#     # apply params to mmolec basis
#     tmpmm = fitt_data.basis_mmol.copy()
#     fremm = np.exp( 1j * np.outer(arr1,frmm) * tmm)
#     ampmm = np.outer(arr1, ammm)
#     ph0mm = np.outer(arr1, np.exp(1j * (np.zeros(nmm)+ph0mm)))
#     tmpmm *= ampmm * fremm * ph0mm * lshapemm

#     # combine both basis sets and apply lineshape
#     for i in range(nmet):
#         f[0:npts,i] = tmp[:,i]
#     for i in range(nmet,nmet+nmm):
#         f[0:npts,i] = tmpmm[:,i-nmet]

#     f[0,] = f[0,] / 2.0

#     # Calc Phase1 here ---
#     phase1 = np.exp( 1j * (a[nmet*2+3]*DDTOR*(np.arange(nptszf,dtype=float)-piv)/nptszf))

#     # Calculate Partial Derivatives here ---
#     if pderflg:
#         pder = np.zeros((nptszf, np.size(a)), complex)

#         if nmet != 1:
#             pall = np.sum(f,axis=1)   # all lines added
#         else:
#             pall = f

#         pind = f
#         tt         = np.zeros(nptszf,float)
#         tt[0:npts] = np.arange(npts,dtype=float) * td

#         for i in range(nmet):   # Calc the Ampl and Freq pders
#             pder[:,i]      = (np.fft.fft(pind[:,i] / a[i])/nptszf) * phase1
#             pder[:,i+nmet] = (np.fft.fft(tt * 1j * pind[:,i])/nptszf) * phase1

#         pder[:,nmet*2+0]  = (np.fft.fft(     tt     * pall/(a[nmet*2+0]**2))/nptszf) * phase1
#         pder[:,nmet*2+1]  = (np.fft.fft(2.0*(tt**2) * pall/(a[nmet*2+1]**3))/nptszf) * phase1

#         pder[:,nmet*2+2]  = (np.fft.fft(1j*pall)/nptszf) * phase1 * nptszf
#         pder[:,nmet*2+3]  = (np.fft.fft(   pall)/nptszf) * (1j*DDTOR*(np.arange(nptszf,dtype=float)-piv)/nptszf) * phase1

#         return pder

#     # Do the FFT here ---
#     if indiv:   # return individual lines
#         if nmet+nmm != 1:
#             for i in range(nmet+nmm):
#                 f[:,i] = (np.fft.fft(f[:,i])/nptszf) * phase1
#         else:
#             f = (np.fft.fft(f[:,0])/nptszf) * phase1
#     else:  # return summed spectrum
#         if (nmet+nmm) != 1:
#             f = np.sum(f,axis=1)
#             f = (np.fft.fft(f)/nptszf) * phase1
#         else:
#             f = (np.fft.fft(f[:,0])/nptszf) * phase1

#     # Add in baseline unless nobase is true ---
#     if not nobase:
#         if f.ndim > 1:
#             for i in range(len(f)): f[:,i] = f[:,i] + fitt_data.fit_baseline
#         else:
#             f = f + fitt_data.fit_baseline

#     # ComplexData put in one long array to comply to ---
#     #  Curvefit (constrained_levenberg_marquardt) constraints
#     return f


