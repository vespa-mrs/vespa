# Python modules
import os
import io as StringIO
import base64
import datetime

# 3rd party modules
import wx
from wx.lib.mixins.listctrl import CheckListCtrlMixin, ColumnSorterMixin

# Our modules
import vespa.analysis.tab_base as tab_base
import vespa.analysis.constants as constants
import vespa.analysis.util_menu as util_menu
import vespa.analysis.prefs as prefs_module
import vespa.analysis.figure_layouts as figure_layouts
import vespa.analysis.util_analysis_config as util_analysis_config
import vespa.analysis.dialog_dataset_browser as dialog_dataset_browser
import vespa.analysis.auto_gui.watref as watref

import vespa.common.wx_gravy.util as wx_util
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.util.misc as util_misc

from matplotlib.backends.backend_pdf import PdfPages

#------------------------------------------------------------------------------

_HLSVD_RESULTS_DISPLAY_SIZE = 6

class CheckListCtrl(wx.ListCtrl, CheckListCtrlMixin, ColumnSorterMixin):
    def __init__(self, _inner_notebook, tab):
        style = wx.LC_REPORT | wx.LC_HRULES | wx.LC_VRULES
        wx.ListCtrl.__init__(self, _inner_notebook, -1, style=style)
        CheckListCtrlMixin.__init__(self)
        ColumnSorterMixin.__init__(self, _HLSVD_RESULTS_DISPLAY_SIZE)
        self.itemDataMap = {}
        self._tab_dataset = _inner_notebook
        self.tab = tab
        self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnItemActivated)

    def GetListCtrl(self):
        return self

    def OnItemActivated(self, evt):
        self.ToggleItem(evt.m_itemIndex)

    # this is called by the base class when an item is checked/unchecked
    def OnCheckItem(self, index, flag):
        self.tab.on_check_item(self, index, flag)


#------------------------------------------------------------------------------

def _configure_combo(control, choices, selection=''):
        lines = list(choices.values())
        control.SetItems(lines)
        if selection in lines:
            control.SetStringSelection(selection)
        else:
            control.SetStringSelection(lines[0])


#------------------------------------------------------------------------------
#
#  Tab WATREF
#
#------------------------------------------------------------------------------

class TabWatref(tab_base.Tab, watref.PanelWatrefUI):
    
    # self-identify tab to notebook, value does not matter, its presence is sufficient.
    IS_WATREF = True

    def __init__(self, tab_dataset, top, block):

        watref.PanelWatrefUI.__init__(self, tab_dataset.NotebookDataset)
        tab_base.Tab.__init__(self, tab_dataset, top, prefs_module.PrefsWatref)

        self.top      = top               # application frame
        self.block    = block             # processing object

        self.cursor_span_picks_lines = False

        # Plotting is disabled during some of init. That's because the plot
        # isn't ready to plot, but the population of some controls
        # (e.g. spin controls on the water filter panel) fires their
        # respective change event which triggers a call to plot(). This
        # is a Windows-specific bug.
        # See http://trac.wxwidgets.org/ticket/14583
        # In any case, skipping some calls to plot() will speed things up. =)
        self._plotting_enabled = False
        self.plot_results = None

        self.initialize_controls()
        self.populate_controls()

        self._plotting_enabled = True

        #------------------------------------------------------------
        # Setup the canvas
        self.process_and_plot()

        #------------------------------------------------------------
        # PubSub subscriptions

        #------------------------------------------------------------
        # Set window events

        # If the sash position isn't recorded in the INI file, we use the
        # arbitrary-ish value of 400.
        if not self._prefs.sash_position_main:
            self._prefs.sash_position_main = 600

        # Under OS X, wx sets the sash position to 10 (why 10?) *after*
        # this method is done. So setting the sash position here does no
        # good. We use wx.CallAfter() to (a) set the sash position and
        # (b) fake an EVT_SPLITTER_SASH_POS_CHANGED.
        wx.CallAfter(self.SplitterWindow.SetSashPosition, self._prefs.sash_position_main, True)
        wx.CallAfter(self.on_splitter)




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
        # Here I set up the spin controls. Some are standard wxSpinCtrls and
        # some are floatspins (from wx.lib.agw.floatspin or wx_common).

        self.SpinWaterAverages.SetValue(1)
        wx_util.configure_spin(self.SpinWaterAverages,      60, None, None, (1,10000))
        self.SpinMetaboliteAverages.SetValue(1)
        wx_util.configure_spin(self.SpinMetaboliteAverages, 60, None, None, (1,10000))
        self.FloatSequenceTe.SetValue(30.0)
        wx_util.configure_spin(self.FloatSequenceTe,       60, 3, 1.0, (0.001,1000.0))

        self.FloatTissueContentGm.SetValue(60.0)
        wx_util.configure_spin(self.FloatTissueContentGm,  70, 2, 1.0, (0.0,100.0))
        self.FloatTissueContentWm.SetValue(40.0)
        wx_util.configure_spin(self.FloatTissueContentWm,  70, 2, 1.0, (0.0,100.0))
        self.FloatTissueContentCsf.SetValue(0.0) 
        wx_util.configure_spin(self.FloatTissueContentCsf, 70, 2, 1.0, (0.0,100.0))
        self.FloatWaterContentGm.SetValue(0.78)
        wx_util.configure_spin(self.FloatWaterContentGm,   70, 2, 0.1, (0.001,1.5))
        self.FloatWaterContentWm.SetValue(0.65)
        wx_util.configure_spin(self.FloatWaterContentWm,   70, 2, 0.1, (0.001,1.5))
        self.FloatWaterContentCsf.SetValue(0.97)
        wx_util.configure_spin(self.FloatWaterContentCsf,  70, 2, 0.1, (0.001,1.5))
        self.FloatWaterT2Gm.SetValue(110.0)
        wx_util.configure_spin(self.FloatWaterT2Gm,        70, 3, 1.0, (0.001,1000.0))
        self.FloatWaterT2Wm.SetValue(80.0)
        wx_util.configure_spin(self.FloatWaterT2Wm,        70, 3, 1.0, (0.001,1000.0))
        self.FloatWaterT2Csf.SetValue(350.0)
        wx_util.configure_spin(self.FloatWaterT2Csf,       70, 3, 1.0, (0.001,1000.0))
        self.FloatMetaboliteT2.SetValue(160.0)
        wx_util.configure_spin(self.FloatMetaboliteT2,     70, 3, 1.0, (0.001,1000.0))

        #------------------------------------------------------------
        # set up Results output HTML widget
        #------------------------------------------------------------

        html_sizer = self.LabelHtmlPlaceholderWatref.GetContainingSizer()
        parent = self.LabelHtmlPlaceholderWatref.GetParent()
        self.LabelHtmlPlaceholderWatref.Destroy()
        self.results_html_ctrl = wx.html.HtmlWindow(parent)
        html_sizer.Add(self.results_html_ctrl, 1, wx.EXPAND|wx.ALIGN_TOP)


    def populate_controls(self, preset=False):
        """
        Populates the widgets with relevant values from the data object.
        It's meant to be called when a new data object is loaded.

        This function trusts that the data object it is given doesn't violate
        any rules. Whatever is in the data object gets slapped into the
        controls, no questions asked.

        """
        dataset = self.dataset
        voxel   = self._tab_dataset.voxel

        self.TextReferenceDataset.SetValue(self.block.set.watref_filename)

        self.SpinWaterAverages.SetValue(self.block.set.water_averages)
        self.SpinMetaboliteAverages.SetValue(self.block.set.metabolite_averages)
        self.FloatSequenceTe.SetValue(self.block.set.sequence_te)
        self.FloatTissueContentGm.SetValue(self.block.set.tissue_content_gm)
        self.FloatTissueContentWm.SetValue(self.block.set.tissue_content_wm)
        self.FloatTissueContentCsf.SetValue(self.block.set.tissue_content_csf)
        self.FloatWaterContentGm.SetValue(self.block.set.water_content_gm)
        self.FloatWaterContentWm.SetValue(self.block.set.water_content_wm)
        self.FloatWaterContentCsf.SetValue(self.block.set.water_content_csf)
        self.FloatWaterT2Gm.SetValue(self.block.set.water_t2_gm)
        self.FloatWaterT2Wm.SetValue(self.block.set.water_t2_wm)
        self.FloatWaterT2Csf.SetValue(self.block.set.water_t2_csf)
        self.FloatMetaboliteT2.SetValue(self.block.set.metabolite_t2)
        
        self.PanelWaterCorrection.Enable(self.block.set.apply_water_correction)
        self.PanelMetaboliteCorrection.Enable(self.block.set.apply_metabolite_correction)



    #=======================================================
    #
    #           Global and Menu Event Handlers
    #
    #=======================================================

    def on_activation(self):
        tab_base.Tab.on_activation(self)
        voxel = self._tab_dataset.voxel
        self.process_and_plot()


    def on_menu_view_option(self, event):
          pass

    def on_menu_view_output(self, event):
        event_id = event.GetId()
        #        if event_id == util_menu.ViewIdsWatref.OUTPUT_RESULTS_TEXT:
        #            pass
        if event_id == util_menu.ViewIdsWatref.OUTPUT_CURRENT_VOXEL_CSV:
            self._output_csv(event, all_voxels=False)
        elif event_id == util_menu.ViewIdsWatref.OUTPUT_ALL_VOXELS_CSV:
            self._output_csv(event, all_voxels=True)

    def on_menu_view_results(self, event):
        """
        Options for saving results to files.

        Two main flavors : Text (csv) and Images (PDF and PNG)

        Text output can be for current voxel, or all voxels in the dataset. This
        is a pretty simplistic method. It outputs whatever is in the results
        array, zeros or whatever.

        Image output come in one layout - LCModel-like

        """

        # map menu items to file types
        EVENT_TO_TYPE = { util_menu.ViewIdsWatref.RESULTS_CSV_CURRENT   : ("text","current","csv"),
                          util_menu.ViewIdsWatref.RESULTS_CSV_ALL       : ("text","all","csv"),
                          util_menu.ViewIdsWatref.RESULTS_LCM_PDF       : ("image","lcm","pdf"),
                          util_menu.ViewIdsWatref.RESULTS_LCM_PDF_MULTI : ("image","lcm_multi","pdf"),
                          util_menu.ViewIdsWatref.RESULTS_LCM_PNG       : ("image","lcm","png"),
                        }

        flavor, layout, type_ = EVENT_TO_TYPE[event.GetId()]

        if flavor == 'text':

            ini_name = "watref_output_as_csv"
            default_path_name = util_analysis_config.get_path(ini_name)
            default_path, default_fname = os.path.split(default_path_name)

            if type_ == 'csv':
                filetype_filter = "CSV files (*.csv)|*.csv|TXT files (*.txt)|*.txt"
            else:
                filetype_filter = "%s files (*.%s)|*.%s" % (type_.upper(), type_, type_)

            filename = common_dialogs.save_as(default_path=default_path,
                                              filetype_filter=filetype_filter,
                                              default_filename=default_fname)
            if filename:
                try:
                    if layout == 'current':
                        voxels = [self._tab_dataset.voxel, ]
                    else:
                        voxels = self.dataset.all_voxels

                    self.output_results_to_csv(voxels, filename)

                except IOError:
                    msg = """I can't write the file "%s".""" % filename
                    common_dialogs.message(msg, style=common_dialogs.E_OK)
                else:
                    # We saved results, so we write the path to the INI file.
                    util_analysis_config.set_path(ini_name, filename)
            # else:
            # User clicked cancel on the "save as" dialog


        elif flavor == 'image':  # image layout result

            ini_name = "watref_results_output_as_image"
            voxel    = self._tab_dataset.voxel
            dpath    = util_analysis_config.get_path(ini_name)
            filtr    = "%s files (*.%s)|*.%s" % (type_.upper(), type_, type_)
            filename = common_dialogs.save_as(default_path=dpath, filetype_filter=filtr)

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

            if filename:

                viffpath = 'Analysis - Quant Tab'
                vespa_version = util_misc.get_vespa_version()
                timestamp = ''  # fig_call will provide
                minplot = 0.1
                maxplot = 4.9
                fontname = 'Courier New'  # 'DejaVu Sans'  #'Ariel'  #'Courier New'

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
                                  quantvals=True,
                                  voxel=voxel)

                if len(figure) == 1:
                    try:
                        figure[0].savefig(filename,
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



    #=======================================================
    #
    #           Widget Event Handlers
    #
    #=======================================================

    def on_tab_changed(self, event):
        voxel = self._tab_dataset.voxel


    def on_splitter(self, event=None):
        # This is sometimes called programmatically, in which case event is None
        self._prefs.sash_position_main = self.SplitterWindow.GetSashPosition()


    def on_reference_browse(self, event):
        # Allows the user to select an ECC dataset.
 
        dialog = dialog_dataset_browser.DialogDatasetBrowser(self.top.datasets)
        dialog.ShowModal()
        ref_dataset = dialog.dataset
        dialog.Destroy()

        if ref_dataset:
            
            self.dataset.set_associated_dataset_quant(ref_dataset)
            
#             block    = ref_dataset.blocks["raw"]
#             self.block.set.watref_dataset    = ref_dataset
#             self.block.set.watref_dataset_id = ref_dataset.id
#             self.block.set.watref_filename   = block.data_source
            
            self.TextReferenceDataset.SetValue(ref_dataset.blocks["raw"].data_source)
            self.process_and_plot()

    def on_sequence_te(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.sequence_te = value
        self.process_and_plot()

    def on_water_averages(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.water_averages = value
        self.process_and_plot()

    def on_metabolite_averages(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.metabolite_averages = value
        self.process_and_plot()

    def on_apply_water_correction(self, event):
        val = event.GetEventObject().GetValue()
        self.block.set.apply_water_correction = val
        if val:
            self.PanelWaterCorrection.Enable()
        else:
            self.PanelWaterCorrection.Disable()
        self.process_and_plot()

    def on_tissue_content_gm(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.tissue_content_gm = value
        self.process_and_plot()

    def on_tissue_content_wm(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.tissue_content_wm = value
        self.process_and_plot()

    def on_tissue_content_csf(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.tissue_content_csf = value
        self.process_and_plot()

    def on_water_content_gm(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.water_content_gm = value
        self.process_and_plot()

    def on_water_content_wm(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.water_content_wm = value
        self.process_and_plot()

    def on_water_content_csf(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.water_content_csf = value
        self.process_and_plot()

    def on_water_t2_gm(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.water_t2_gm = value
        self.process_and_plot()

    def on_water_t2_wm(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.water_t2_wm = value
        self.process_and_plot()

    def on_water_t2_csf(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.water_t2_csf = value
        self.process_and_plot()

    def on_apply_metabolite_correction(self, event):
        val = event.GetEventObject().GetValue()
        self.block.set.apply_metabolite_correction = val
        if val:
            self.PanelMetaboliteCorrection.Enable()
        else:
            self.PanelMetaboliteCorrection.Disable()
        self.process_and_plot()

    def on_metabolite_t2(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.metabolite_t2 = value
        self.process_and_plot()


    def on_output_html(self, event):
        ini_name = "watref_output_as_html"
        default_path = util_analysis_config.get_path(ini_name)
        default_filename = "analysis_watref_results.html"
        filetype_filter = "HTML (*.html)|*.html"
        
        filename = common_dialogs.save_as(default_path=default_path,
                                          filetype_filter=filetype_filter,
                                          default_filename=default_filename)
        if filename:

            voxel   = self._tab_dataset.voxel
            raw     = self.dataset.blocks["raw"]
            fit     = self.dataset.blocks["fit"]
            tab_fit = self._tab_dataset._tabs["fit"]
            
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
        
            tab_fit.view.figure.savefig(fake_file, dpi=300, format=format)
            
            image_data = base64.b64encode(fake_file.getvalue())
        
            fake_file.close()

            fitted_lw = fit.chain.fitted_lw
            lw_min, lw_max = fit.chain.minmaxlw

            data_source = raw.get_data_source(voxel)
            
            #html = self.block.results_as_html(voxel, fit, fitted_lw,
            html = self.dataset.quant_results_as_html(voxel, fitted_lw, 
                                                      lw_min, lw_max, 
                                                      data_source,
                                                      (mime_type, image_data)
                                                     )

            open(filename, "wb").write(html)
            
            # We saved results, so we write the path to the INI file.
            path, _ = os.path.split(filename)
            util_analysis_config.set_path(ini_name, path)

        #else:
            # User clicked cancel on the "save as" dialog


#     def on_output_csv_orig(self, event):
#         ini_name = "watref_output_as_csv"
#         default_path_name = util_analysis_config.get_path(ini_name)
#         default_path, default_fname = os.path.split(default_path_name)
#         filetype_filter = "CSV (*.csv)|*.csv"
#         
#         filename = common_dialogs.save_as(default_path=default_path,
#                                           filetype_filter=filetype_filter,
#                                           default_filename=default_fname)
# 
#         if filename:
#             
#             # Create output header and results strings, check element count. 
#             # If the file exists, check that the element count is the same in 
#             # in the last line as for this results line. If it is, just write
#             # out the results string. If different length, output both the 
#             # header and results strings.
# 
#             decor1 = self._prefs.csv_qa_metab_labels
# 
#             voxel = self._tab_dataset.voxel
#             raw   = self.dataset.blocks["raw"]
#             data_source = raw.get_data_source(voxel)
#             dataset_filename = self.dataset.dataset_filename
# 
#             fitted_lw = fit.chain.fitted_lw
#             lw_min, lw_max = fit.chain.minmaxlw
#             
#             val, hdr = self.dataset.quant_results_as_csv(voxel, fitted_lw, 
#                                                          lw_min, lw_max,
#                                                          data_source,
#                                                          dataset_filename,
#                                                          decor1=decor1)
#             nhdr = len(hdr)
#             val = ",".join(val)
#             hdr = ",".join(hdr)
#             val += "\n"
#             hdr += "\n"
#              
#             hdr_flag = True
#             if os.path.isfile(filename):
#                 with open(filename, 'r+') as f:
#                     data = f.readlines()
#                     if len(data)>1:
#                         last = data[-1]
#                         nlast = len(last.split(','))
#                         if nlast == nhdr:
#                             hdr_flag = False
#                         
#             with open(filename, 'a') as f:
#                 if hdr_flag:
#                     f.write(hdr)
#                 f.write(val)
# 
#             # We saved results, so we write the path to the INI file.
#             util_analysis_config.set_path(ini_name, filename)
# 
#         #else:
#             # User clicked cancel on the "save as" dialog
            

    def on_output_csv(self, event):
        """ 
        Send all voxel fitting results to CSV file, save the 
        last successful CSV file as the default for the next call.
        
        """
        self._output_csv(event, all_voxels=False)
        
#         voxels = [self._tab_dataset.voxel,]
#         
#         ini_name = "watref_output_as_csv"
#         default_path_name = util_analysis_config.get_path(ini_name)
#         default_path, default_fname = os.path.split(default_path_name)
#         filetype_filter = "CSV (*.csv)|*.csv"
#         
#         filename = common_dialogs.save_as(default_path=default_path,
#                                           filetype_filter=filetype_filter,
#                                           default_filename=default_fname)
# 
#         if filename:
# 
#             try:
#                 self.output_results_to_csv(voxels, filename)
# 
#             except IOError:
#                 msg = """I can't write the file "%s".""" % filename
#                 common_dialogs.message(msg, style=common_dialogs.E_OK)
#             else:
#                 # We saved results, so we write the path to the INI file.
#                 util_analysis_config.set_path(ini_name, filename)
#             
#         #else:
#             # User clicked cancel on the "save as" dialog            


    def _output_csv(self, event, all_voxels=False):
        """ 
        Send all voxel fitting results to CSV file, save the 
        last successful CSV file as the default for the next call.
        
        """
        if not all_voxels:
            voxels = [self._tab_dataset.voxel,]
        else:
            voxels = self.dataset.all_voxels
        
        ini_name = "watref_output_as_csv"
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
        fit   = self.dataset.blocks["fit"]

        for voxel in voxels:

            data_source = raw.get_data_source(voxel)
            dataset_filename = self.dataset.dataset_filename
            
            fitted_lw = fit.chain.fitted_lw
            lw_min, lw_max = fit.chain.minmaxlw

            # FIXME - bjs, need to update current 'fitted_lw' and 'min/maxlw' 
            #         values if we save more than one voxel. For now, the 
            #         existing values for current voxel are used.

            #val, hdr = self.block.results_as_csv(voxel, fit,
            val, hdr = self.dataset.quant_results_as_csv(voxel,  
                                                         lw = fit.chain.fitted_lw, 
                                                         lwmin = fit.chain.minmaxlw[0], 
                                                         lwmax = fit.chain.minmaxlw[1], 
                                                         source=data_source,
                                                         dsetname=dataset_filename,
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



    def on_batch_process_all_voxels(self, event): 
        
        # check that if mmol model is on, there is a dataset selected

        block_fit = self.dataset.blocks['fit']
        
        method = block_fit.set.macromol_model 
        if (method == 'single_basis_dataset') and (block_fit.set.macromol_single_basis_dataset is None):
            # inconsistent, reset widget and mmol parameters
            method = list(constants.FitMacromoleculeMethod.choices.keys())[0]   # None method
            block_fit.set.macromol_model = method
            block_fit.set.macromol_single_basis_dataset = None
            block_fit.set.macromol_single_basis_dataset_id = ''
        
        self._tab_dataset.batch_process_all(statusbar=self.top.statusbar)   # bjs - maybe could do direct call to self._dataset.batch_fit_all() ???
        
        self.process_and_plot(entry='voxel_change')
                  

    #=======================================================
    #
    #           Public Methods
    #
    #=======================================================

    def update_html_results_tab(self):
        """ 
        Update display of results in the html window.
        
        At the moment we are only set up to have multiple files loaded into
        the first index, ie. SVS files loaded into the window, thus we use
        the voxel[0] index to select the filename.
        """
        voxel   = self._tab_dataset.voxel
        raw     = self.dataset.blocks["raw"]
        fit     = self.dataset.blocks["fit"]
        tab_fit = self._tab_dataset._tabs["fit"]

        #html = self.block.results_as_html(voxel, fit,
        html = self.dataset.quant_results_as_html(voxel, 
                                                  tab_fit.plot_results['fitted_lw'],
                                                  tab_fit.plot_results['minmaxlw'][0],
                                                  tab_fit.plot_results['minmaxlw'][1], 
                                                  raw.get_data_source(voxel))

        self.results_html_ctrl.SetPage(html)
        self.results_html_ctrl.SetBackgroundColour(self.GetBackgroundColour())        


    def process_and_plot(self, entry='all'):
        """
        The process(), plot() and process_and_plot() methods are standard in
        all processing tabs. They are called to update the data in the plot
        results dictionary, the plot_panel in the View side of the tab or both.

        """
        tab_base.Tab.process_and_plot(self, entry)

        if self._plotting_enabled:
            self.process(entry=entry)
            self.plot()


    def process(self, entry='all'):
        """
        Data processing results are stored into the Block inside the Chain,
        but the View results are returned as a dictionary from the Chain.run()
        method. The plot routine takes its inputs from this dictionary.

        """
        tab_base.Tab.process(self, entry)

        if self._plotting_enabled:

            dataset = self.dataset
            voxel   = self._tab_dataset.voxel

            self.plot_results = self.block.chain.run([voxel], entry=entry)

            label1  = 'Relaxation Corrected Water Concentration [M] = '
            label1 += self.plot_results['wat_relax_corr'] 
            label2  = 'Metabolite Relaxation Correction = '
            label2 += self.plot_results['met_relax_corr']
            label3  = 'Acquisition Averages Correction = '
            label3 += self.plot_results['averages_corr']

            self.TextRelaxationCorrectedWaterConcentration.SetLabel(label1)
            self.TextMetaboliteRelaxationCorrection.SetLabel(label2)
            self.TextAveragesCorrection.SetLabel(label3)


    def plot(self):
        ''' write something nice here '''
        self.update_html_results_tab()


    def on_voxel_change(self, voxel):
        # this just updates widgets that vary based on the voxel number
        # selection. We do not update plot here because that is only done
        # for the active tab in the inner notebook.
        pass


    #=======================================================
    #
    #           Internal Helper Functions
    #
    #=======================================================









