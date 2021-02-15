# Python modules


# Our modules 
import vespa.common.constants as constants

"""When Vespa needs to create an app's INI file, this is the content it
uses. This could be part of constants.py, but the long triple-quoted
strings make it 
really hard to read so it's easier if this is isolated in its own module.

This module contains just one object which is a dict called 
DEFAULT_INI_FILE_CONTENT. The dict has a key for each INI file ("vespa", 
"simulation", etc.) and associated with that key is the default content
for that INI file. 
"""

# The dict contents are mostly strings but at least one (simulation) 
# uses string interpolation. 



DEFAULT_INI_FILE_CONTENT = {

###############################      Vespa

    "vespa" : """
# The Vespa config file.

[database]
filename=vespa.sqlite

[general]
last_export_path=
# If cpu_limit is set to a positive integer, Vespa will use no more than that
# many CPUs simultaneously. When cpu_limit is set to 0 or left blank, Vespa
# will use as many CPUs as it thinks are available.
cpu_limit=

""",



###############################      Pulse 

    "pulse" : """
# The Pulse config file.

# Colors are described in matplotlib's terms. Matplotlib understands standard
# color names like "red", "black", "blue", etc.
# One can also specify the color using an HTML hex string like #eeefff.
# If you use a string like that, you must put it in quotes, otherwise the
# hash mark will get interpreted as a comment marker. For example, to set
# a background to pale yellow:
#    bgcolor = "#f3f3bb"

[general]
last_file_open_path = 
warn_when_discarding_results = yes

[main]
left = 40
top = 40
width = 1000
height = 740

[main_prefs]
bgcolor = "#ffffff"
data_type_real = True
data_type_real_imaginary = False
line_color_imaginary = red
line_color_magnitude = purple
line_color_phase_degrees = orange
line_color_real = black
line_width = 1.0
sash_position = 500
xaxis_show = True
zero_line_color = goldenrod
zero_line_show = False
zero_line_style = solid

[dialog_editor_transform]
left = 20
top = 20
width = 1200
height = 850

""",



###############################      Simulation 
    "simulation" : """
# The Simulation config file.

# Colors are described in matplotlib's terms. Matplotlib understands standard
# color names like "red", "black", "blue", etc.
# One can also specify the color using an HTML hex string like #eeefff.
# If you use a string like that, you must put it in quotes, otherwise the
# hash mark will get interpreted as a comment marker. For example, to set
# a background to pale yellow:
#    bgcolor = "#f3f3bb"

[general]
last_file_open_path = 
default_isotope = %s
default_b0 = %s

[main]
left = 40
top = 40
width = 1200
height = 800

[main_prefs]
sash_position = 400
bgcolor = "#ffffff"
contour_axes_show = False
contour_grayscale = False
contour_levels = 1
contour_plot_show = False
data_type_imaginary = False
data_type_magnitude = False
data_type_real = True
integral_plot_show = False
integral_xaxis_show = False
integral_yaxis_show = False
line_color_imaginary = red
line_color_magnitude = purple
line_color_real = black
line_shape_gaussian = True
line_shape_lorentzian = False
line_shape_stick = False
line_width = 1.0
xaxis_hertz = True
xaxis_ppm = True
xaxis_show = False
zero_line_bottom = False
zero_line_color = goldenrod
zero_line_middle = True
zero_line_show = False
zero_line_style = solid
zero_line_top = False

[dialog_pulse_sequence_editor]
left = 20
top = 20
width = 1250
height = 850

""" % (constants.DEFAULT_ISOTOPE, constants.DEFAULT_B0),

###############################      Analysis

    "analysis" : """
# The Analysis config file.

# Colors are described in matplotlib's terms. Matplotlib understands standard
# color names like "red", "black", "blue", etc.
# One can also specify the color using an HTML hex string like #eeefff.
# If you use a string like that, you must put it in quotes, otherwise the
# hash mark will get interpreted as a comment marker. For example, to set
# a background to pale yellow:
#    bgcolor = "#f3f3bb"

[main]
left = 40
top = 40
width = 1200
height = 800
maximized = False

[basic_prefs]
bgcolor = "#ffffff"
line_color_imaginary = red
line_color_individual = green
line_color_magnitude = purple
line_color_real = black
line_color_summed = black
line_width = 1.0
sash_position = 550
xaxis_hertz = False
xaxis_ppm = True
xaxis_show = True
zero_line_bottom = True
zero_line_color = darkgoldenrod
zero_line_middle = False
zero_line_show = False
zero_line_style = solid
zero_line_top = False
csv_qa_metab_labels = False

[prep_fidsum_prefs]
area_calc_plot_a = False
area_calc_plot_b = True
bgcolor = "#ffffff"
data_type_imaginary = False
data_type_magnitude = False
data_type_real = True
line_color_imaginary = red
line_color_magnitude = purple
line_color_real = black
line_width = 1.0
sash_position = 522
xaxis_hertz = False
xaxis_ppm = True
xaxis_show = False
zero_line_bottom = True
zero_line_color = darkgoldenrod
zero_line_middle = False
zero_line_show = False
zero_line_style = solid
zero_line_top = False
zero_line_plot_bottom = False
zero_line_plot_color = yellow
zero_line_plot_middle = False
zero_line_plot_show = True
zero_line_plot_style = solid
zero_line_plot_top = True


[spectral_prefs]
area_calc_plot_a = True
area_calc_plot_b = False
area_calc_plot_c = False
bgcolor = "#ffffff"
data_type_imaginary = False
data_type_magnitude = False
data_type_real = True
line_color_imaginary = red
line_color_magnitude = purple
line_color_real = black
line_color_svd = green
line_width = 1.0
plot_c_function_a_minus_b = True
plot_c_function_a_plus_b = False
plot_c_function_b_minus_a = False
plot_c_function_none = False
sash_position_main = 400
sash_position_svd = 487
xaxis_hertz = False
xaxis_ppm = True
xaxis_show = False
zero_line_bottom = True
zero_line_color = darkgoldenrod
zero_line_middle = False
zero_line_show = False
zero_line_style = solid
zero_line_top = False

[voigt_prefs]
area_calc_plot_a = True
area_calc_plot_b = False
area_calc_plot_c = False
area_calc_plot_d = False
bgcolor = "#ffffff"
csv_qa_metab_labels = False
line_color_base = purple
line_color_fit = green
line_color_imaginary = red
line_color_init = green
line_color_magnitude = purple
line_color_real = black
line_color_raw = black
line_color_weight = darkgoldenrod
line_width = 1.0
n_plots_1 = False
n_plots_2 = False
n_plots_3 = False
n_plots_4 = True
sash_position = 550
xaxis_hertz = False
xaxis_ppm = True
xaxis_show = False
zero_line_bottom = True
zero_line_color = darkgoldenrod
zero_line_middle = False
zero_line_show = False
zero_line_style = solid
zero_line_top = False

[voigt_plot_a]
baseline = False
bgcolor = black
combo_fit_and_base = False
combo_fit_plus_base = False
combo_raw_and_base = False
combo_raw_and_fit = False
combo_raw_and_fit_plus_base = False
combo_raw_and_init_model = True
combo_raw_and_wt_arr = False
combo_raw_minus_base = False
combo_raw_minus_base_and_fit = False
combo_raw_minus_fit = False
combo_raw_minus_fit_minus_base = False
data_type_imaginary = False
data_type_magnitude = False
data_type_real = True
fitted_data = False
raw_data = False

[voigt_plot_b]
baseline = False
bgcolor = black
combo_fit_and_base = False
combo_fit_plus_base = False
combo_raw_and_base = True
combo_raw_and_fit = False
combo_raw_and_fit_plus_base = False
combo_raw_and_init_model = False
combo_raw_and_wt_arr = False
combo_raw_minus_base = False
combo_raw_minus_base_and_fit = False
combo_raw_minus_fit = False
combo_raw_minus_fit_minus_base = False
data_type_imaginary = False
data_type_magnitude = False
data_type_real = True
fitted_data = False
raw_data = False

[voigt_plot_c]
baseline = False
bgcolor = black
combo_fit_and_base = False
combo_fit_plus_base = False
combo_raw_and_base = False
combo_raw_and_fit = False
combo_raw_and_fit_plus_base = True
combo_raw_and_init_model = False
combo_raw_and_wt_arr = False
combo_raw_minus_base = False
combo_raw_minus_base_and_fit = False
combo_raw_minus_fit = False
combo_raw_minus_fit_minus_base = False
data_type_imaginary = False
data_type_magnitude = False
data_type_real = True
fitted_data = False
raw_data = False

[voigt_plot_d]
baseline = False
bgcolor = black
combo_fit_and_base = False
combo_fit_plus_base = False
combo_raw_and_base = False
combo_raw_and_fit = False
combo_raw_and_fit_plus_base = False
combo_raw_and_init_model = False
combo_raw_and_wt_arr = False
combo_raw_minus_base = False
combo_raw_minus_base_and_fit = False
combo_raw_minus_fit = False
combo_raw_minus_fit_minus_base = True
data_type_imaginary = False
data_type_magnitude = False
data_type_real = True
fitted_data = False
raw_data = False

[giso_prefs]
area_calc_plot_a = True
area_calc_plot_b = False
area_calc_plot_c = False
area_calc_plot_d = False
bgcolor = "#ffffff"
line_color_base = purple
line_color_fit = green
line_color_imaginary = red
line_color_init = green
line_color_magnitude = purple
line_color_real = black
line_color_raw = black
line_color_weight = darkgoldenrod
line_width = 1.0
n_plots_1 = False
n_plots_2 = False
n_plots_3 = False
n_plots_4 = True
sash_position = 550
xaxis_hertz = False
xaxis_ppm = True
xaxis_show = False
zero_line_bottom = True
zero_line_color = darkgoldenrod
zero_line_middle = False
zero_line_show = False
zero_line_style = solid
zero_line_top = False

[giso_plot_a]
baseline = False
bgcolor = black
combo_fit_and_base = False
combo_fit_plus_base = False
combo_raw_and_base = False
combo_raw_and_fit = False
combo_raw_and_fit_plus_base = False
combo_raw_and_init_model = True
combo_raw_and_wt_arr = False
combo_raw_minus_base = False
combo_raw_minus_base_and_fit = False
combo_raw_minus_fit = False
combo_raw_minus_fit_minus_base = False
data_type_imaginary = False
data_type_magnitude = False
data_type_real = True
fitted_data = False
raw_data = False

[giso_plot_b]
baseline = False
bgcolor = black
combo_fit_and_base = False
combo_fit_plus_base = False
combo_raw_and_base = True
combo_raw_and_fit = False
combo_raw_and_fit_plus_base = False
combo_raw_and_init_model = False
combo_raw_and_wt_arr = False
combo_raw_minus_base = False
combo_raw_minus_base_and_fit = False
combo_raw_minus_fit = False
combo_raw_minus_fit_minus_base = False
data_type_imaginary = False
data_type_magnitude = False
data_type_real = True
fitted_data = False
raw_data = False

[giso_plot_c]
baseline = False
bgcolor = black
combo_fit_and_base = False
combo_fit_plus_base = False
combo_raw_and_base = False
combo_raw_and_fit = False
combo_raw_and_fit_plus_base = True
combo_raw_and_init_model = False
combo_raw_and_wt_arr = False
combo_raw_minus_base = False
combo_raw_minus_base_and_fit = False
combo_raw_minus_fit = False
combo_raw_minus_fit_minus_base = False
data_type_imaginary = False
data_type_magnitude = False
data_type_real = True
fitted_data = False
raw_data = False

[giso_plot_d]
baseline = False
bgcolor = black
combo_fit_and_base = False
combo_fit_plus_base = False
combo_raw_and_base = False
combo_raw_and_fit = False
combo_raw_and_fit_plus_base = False
combo_raw_and_init_model = False
combo_raw_and_wt_arr = False
combo_raw_minus_base = False
combo_raw_minus_base_and_fit = False
combo_raw_minus_fit = False
combo_raw_minus_fit_minus_base = True
data_type_imaginary = False
data_type_magnitude = False
data_type_real = True
fitted_data = False
raw_data = False



[watref_prefs]
csv_qa_metab_labels = False
sash_position_main = 400

""",


###############################      DataSim

    "datasim" : """
# The DataSim config file.

# Colors are described in matplotlib's terms. Matplotlib understands standard
# color names like "red", "black", "blue", etc.
# One can also specify the color using an HTML hex string like #eeefff.
# If you use a string like that, you must put it in quotes, otherwise the
# hash mark will get interpreted as a comment marker. For example, to set
# a background to pale yellow:
#    bgcolor = "#f3f3bb"


[main]
left = 40
top = 40
width = 1200
height = 800

[main_prefs]
bgcolor = "#ffffff"
data_type_imaginary = False
data_type_magnitude = False
data_type_real = True
line_color_baseline = blue
line_color_imaginary = red
line_color_magnitude = purple
line_color_metabolite = blue
line_color_real = black
line_width = 1.0
plot_view_all = True
plot_view_final = False
sash_position = 556
xaxis_hertz = False
xaxis_ppm = True
xaxis_show = True
zero_line_bottom = True
zero_line_color = goldenrod
zero_line_middle = False
zero_line_show = False
zero_line_style = solid
zero_line_top = False

""" 
,

###############################      Analysis Menu Additions

    "analysis_import_menu_additions" : r"""
# The INI file for additions to Analysis' File/Import menu.

# If you want Analysis to be able to import a file format that it doesn't 
# currently understand, you can write a reader for it and reference it here
# The reference here will make it appear on Analysis' File/Import menu along
# with the other 3rd party formats.
#
# See vespa/analysis/src/fileio/template.py for more information on how to write 
# a class that can read you custom format. You can also see template.py online:
# http://scion.duhs.duke.edu/vespa/project/browser/trunk/analysis/src/fileio/template.py
#
# There should be one menu item for each file format, and one section in this
# file for each menu item you want to add.
# 
# Section names should be valid Python module names (all lower case and no
# spaces will work fine) and they shouldn't conflict with the name of any 
# built-in Python modules. Also note that "dicom" is reserved; you may not use
# that name.
# 
# Each section should contain four entries --
#    - path -
#    The fully qualified path to the file that implements your reader (which
#    is based on template.py, as mentioned above).
#    - class_name - 
#    The name of the implementation class in your reader file.
#    - menu_item_text - 
#    The text that will appear on the menu for this item. Usually this is just 
#    a single word or two. If you want to get fancy, you can define an 
#    accelerator for it that will allow you to activate your menu item with 
#    a single keystroke. Append "\tCtrl+Y" to your menu text and you'll be 
#    able to activate it with Ctrl+Y, or Cmd+Y on the Mac. You can use any
#    letter you like (not just 'Y'), just make sure yours doesn't conflict with
#    another menu item!
#    - ini_file_name - 
#    We use this text in our INI file to record the last path from which you
#    opened a file of this type. Just about any string will do, but for 
#    simplicity we strongly recommend all lower case, and no spaces or 
#    punctuation besides underscores or dashes.
#
# Here's an example that reads the fictitious Acme format.
# 
#    [acme]
#    path=/home/wile_e_coyote/acme.py
#    class_name=RawReaderAcme
#    menu_item_text=Acme
#    ini_file_name=import_acme
#
# Here's an example that reads the fictitious XYZ format. Note that its menu
# item uses an accelerator.
#
#    [xyz]
#    path=C:\Users\philip\xyz_reader.py
#    class_name=RawReaderXyz
#    menu_item_text=XYZ\tCtrl+Y
#    ini_file_name=import_xyz


""" 
,}

