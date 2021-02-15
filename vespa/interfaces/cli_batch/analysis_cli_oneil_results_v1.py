# Python modules
import os
import sys

# 3rd party modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, NullFormatter

# Our modules
import vespa.analysis.util_import as util_import
import vespa.common.util.ppm as util_ppm
import vespa.common.util.misc as util_misc

SUPPORTED = ['wbnaa', 'siemens dicom']

DESC =  \
"""Command line interface to process MRS data in Vespa-Analysis. 
 Data filename, preset file name, data type string and CSV output 
 file name values are all required for this command to function 
 properly.
  
 Note. You may have to enclose data/preset/output strings in double 
 quotation marks for them to process properly if they have  
 spaces or other special characters embedded in them.
"""

#mpl.rc('text', usetex=True) 
#mpl.rcParams['text.usetex']=True
#mpl.rcParams['text.latex.unicode']=True


def analysis_cli_oneil_results(dataset, csvfile, viffpath='', 
                                                 vespa_version='',
                                                 timestamp='',
                                                 verbose=False, debug=False):
    
    if viffpath:
        viffname = os.path.basename(viffpath)
    else:
        viffname = 'none'
        
    raw = dataset.blocks["raw"]
    data_source = raw.get_data_source(dataset.all_voxels)
    data_source = os.path.basename(data_source)

    dim0, dim1, dim2, dim3 = dataset.spectral_dims
    sw    = dataset.sw

    voxel = dataset.all_voxels
    key = 'fit'
    if key in list(dataset.blocks.keys()):
        block = dataset.blocks[key]
        results = block.chain.run(voxel, entry='output_refresh')
    else:
        msg = """This dataset has no 'fit' block, returning.""" 
        print(msg, file=sys.stderr)
        print(msg, file=sys.stdout)
        sys.exit(-1)

    freq  = results['data'].copy()
    base  = results['fit_baseline'].copy()
    base.shape = 1, base.shape[0]
    yfit  = results['yfit'].copy()
    yfits = results['yfit'].copy() if len(yfit.shape)==1 else np.sum(yfit, axis=0) 
    if len(yfit.shape) == 1:
        yfit.shape = 1,yfit.shape[0]
    
    # Print Control Setting
    outbase  = 'D:\\Users\\bsoher\\myplot'
    fontname = 'Courier New'    # 'Consolas' 'Calibri' 'Courier New' 'Times New Roman'
    savetype = 'pdf'            # 'svg' 'eps' 'pdf' 'png' 'raw' 'rgba' 'ps' 'pgf' etc.
    minplot, maxplot = 0.1, 4.9   # in ppm

    xvals = [dataset.pts2ppm(val) for val in range(dim0)]
    minppm, maxppm = xvals[-1], xvals[0]
    
    minplot = minplot if minplot >= minppm else minppm
    maxplot = maxplot if maxplot <= maxppm else maxppm
    
    imin = int(dataset.ppm2pts(maxplot))
    imax = int(dataset.ppm2pts(minplot))
    tmin = np.round((minplot+0.5))                  # integer ppm just above minppm
    tmax = np.round((maxplot-0.5))                  # integer ppm just under maxppm
    
    # Create the figure
    fig = plt.figure(figsize=(11,8.5))
    plt.subplots_adjust(hspace=0.001)
    
    fsupported = plt.gcf().canvas.get_supported_filetypes()
    if savetype not in list(fsupported.keys()):
        msg = r"Output file format '%s' not supported by current Matplotlib backend, Returning." % savetype
        print(msg, file=sys.stderr)
        print(msg, file=sys.stdout)
        sys.exit(-1)

    outname = outbase+'.'+savetype
    nullfmt = NullFormatter()           # no labels

    left, bottom    = 0.05, 0.05        # set up for 8.5x11 landscape printout
    w1, w2          = 0.55, 0.35
    h1, h2, h3      = 0.07, 0.61, 0.07
    hpad, vpad      = 0.02, 0.001    
    
    rect1 = [left,         bottom+h1+h2, w1, h3]
    rect2 = [left,         bottom+h1,    w1, h2]
    rect3 = [left,         bottom,       w1, h1]    # xmin, ymin, dx, and dy 
    rect4 = [left+w1+hpad, bottom, w2, h1+h2+h3]
    
    dat1 = freq[imin:imax] - (yfits+base[0,:])[imin:imax]
    dat1 = dat1.real
    min1, max1 = min(dat1),max(dat1)
    delt1 = (max1 - min1)*0.8
    min1, max1 = min1 - delt1, max1 + delt1
#    tmin1, tmax1 = int(min1+0.5), int(max1-0.5)
#    step1 = int((tmax1 - tmin1)/4)
    
    ax1 = fig.add_axes(rect1)               # xmin, ymin, dx, and dy 
    ax1.xaxis.set_major_formatter(nullfmt)  # no x labels, have to go before plot()
    ax1.plot(xvals[imin:imax],dat1)
    ax1.set_ylabel('some numbers', fontsize=8.0) 
    plt.yticks([0.0,max1])
    plt.ylim(min1, max1)
    plt.xticks(np.arange(tmin, tmax, 1.0))
    plt.xlim(maxplot, minplot)

    ax2 = fig.add_axes(rect2)               # xmin, ymin, dx, and dy
    ax2.xaxis.set_major_formatter(nullfmt)  # no x labels, have to go before plot()
    ax2.plot(xvals[imin:imax],freq[imin:imax])
    ax2.plot(xvals[imin:imax],(yfits+base[0,:])[imin:imax])
    ax2.plot(xvals[imin:imax],(base[0,:])[imin:imax])
    ax2.set_ylabel('some more numbers', fontsize=8.0) 
    plt.yticks(np.arange(-2.0, 13.0, 4.0))
    plt.ylim(-3, 14)
    plt.xticks(np.arange(tmin, tmax, 1.0))
    plt.xlim(maxplot, minplot)
    
    ax3 = fig.add_axes(rect3)               # xmin, ymin, dx, and dy
    ax3.plot(xvals[imin:imax],(base[0,:])[imin:imax])
    ax3.set_ylabel('original numbers',     fontsize=8.0) 
    ax3.set_xlabel('Chemical Shift [ppm]', fontsize=8.0)
    plt.yticks(np.arange(-2.0, 13.0, 4.0))
    plt.ylim(-3, 14)
    plt.xticks(np.arange(tmin, tmax, 1.0))
    plt.xlim(maxplot, minplot)

    for ax in [ax1,ax2,ax3]:
        ax.xaxis.set_major_locator(MultipleLocator(1))        # this is for even distrib across plot
        ax.xaxis.set_minor_locator(MultipleLocator(0.2))
        ax.yaxis.set_major_locator(MultipleLocator(4.0))
#        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))  # yticks() above is finer grained
#        ax.yaxis.set_minor_locator(MultipleLocator(2.0))
        ax.grid(which='major', axis='x', linewidth=0.50, linestyle='-', color='0.75')
        ax.grid(which='minor', axis='x', linewidth=0.25, linestyle=':', color='0.75')
        ax.grid(which='major', axis='y', linewidth=0.25, linestyle=':', color='0.75')
        ax.grid(which='minor', axis='y', linewidth=0.25, linestyle=':', color='0.75')        

    ax4 = fig.add_axes(rect4, axisbg='g')       # xmin, ymin, dx, and dy
    ax4.xaxis.set_major_formatter(nullfmt)
    ax4.yaxis.set_major_formatter(nullfmt)
    ax4.axis('off')
    nrow = 10
    clust_data = np.round(np.random.random((nrow,3)),3)
    collabel=("col 1", "col 2", "col 3")
    rowlabel=[str(i+1) for i in range(nrow)]
    the_table = ax4.table(cellText=clust_data, 
                          cellLoc='left',
                          colLoc='left',
                          colWidths=[0.3,0.3,0.3],
                          colLabels=collabel,
                          rowLoc='center',
                          #rowLabels=rowlabel,
                          fontsize=8.0, 
                          loc='upper center')  
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(7.0)  
    
    table_props = the_table.properties()
    cheight = table_props['children'][0].get_height()   # all start with same default
    keys = list(table_props['celld'].keys())
    for key in keys:                            # use cell dict here to test for col/row labels
        cell = table_props['celld'][key]
        cell.set_height(1.1*cheight)
        cell.get_text().set_fontname(fontname)
        if key[0] == 0:
            if key[1] == 0:
                cell.visible_edges = 'BLT'
            elif key[1] == 2:
                cell.visible_edges = 'BRT'
            else:
                cell.visible_edges = 'BT'
            cell.set_linewidth(1.0)
            cell.set_linestyle('-')
            cell.get_text().set_fontweight('bold')
        else:
            if key[1] == 0:
                cell.visible_edges = 'BL'
            elif key[1] == 2:
                cell.visible_edges = 'BR'
            else:
                cell.visible_edges = 'B'
            cell.set_linewidth(0.25)
            cell.set_linestyle('-')
    
    # Retrieve an element of a plot and set properties
    for tick in ax3.xaxis.get_ticklabels():
        tick.set_fontsize(8.0)
        tick.set_fontname(fontname)
        tick.set_color('gray')
        tick.set_weight('bold')    

    for ax in [ax1,ax2,ax3]:
        ax.yaxis.label.set_fontname(fontname)
        for tick in ax.yaxis.get_ticklabels():
            tick.set_fontsize(8.0)
            tick.set_fontname(fontname)
            tick.set_color('gray')
            tick.set_weight('normal')    

    msg = "Vespa-Analysis Version: %s                            Processing Timestamp: %s" % (vespa_version, timestamp)
    plt.figtext(0.03, 0.95, msg, 
                            wrap=True, 
                            horizontalalignment='left', 
                            fontsize=8, 
                            fontname=fontname)
    msg = "VIFF File   : %s \nData Source : %s" % (viffname, data_source)
    plt.figtext(0.03, 0.90, msg, 
                            wrap=True, 
                            horizontalalignment='left', 
                            fontsize=10, 
                            fontname=fontname)

    
    fig.canvas.draw()
    
    fig.savefig(outname, pad_inches=0.5)#, bbox_inches=(6,8))  #'tight')
     
    bob = 10
    bob += 1





#     # Save results to CSV file --------------------------------------
#  
#     if verbose: print """Saving results to CSV file "%s". """ % csvfile
#      
#     fit = dataset.blocks["fit"]
#     data_source = dataset.blocks["raw"].get_data_source(voxel)
#      
#     val, hdr = fit.results_as_csv(voxel[0], fit.chain.fitted_lw,
#                                             fit.chain.minmaxlw[0],
#                                             fit.chain.minmaxlw[1], 
#                                             data_source, outxml)
#     nhdr = len(hdr)
#     val = ",".join(val)
#     hdr = ",".join(hdr)
#     val += "\n"
#     hdr += "\n"
#       
#     hdr_flag = True
#     if os.path.isfile(csvfile):
#         with open(csvfile, 'r+') as f:
#             data = f.readlines()
#             if len(data)>1:
#                 last = data[-1]
#                 nlast = len(last.split(','))
#                 if nlast == nhdr:
#                     hdr_flag = False
#                  
#     with open(csvfile, 'a') as f:
#         if hdr_flag:
#             f.write(hdr)
#         f.write(val)


def _open_viff(datafile):
    
    datasets = []

    filename = datafile
    timestamp = ''

    msg = ""
    try:
        importer = util_import.DatasetCliImporter(filename)
    except IOError:
        msg = """I can't read the file "%s".""" % filename
    except SyntaxError:
        msg = """The file "%s" isn't valid Vespa Interchange File Format.""" % filename


    if msg:        
        print(msg, file=sys.stderr)
        print(msg, file=sys.stdout)
        sys.exit(-1)
    else:
        # Time to rock and roll!
        dsets, timestamp = importer.go()

        for item in dsets:
            datasets.append(item)


    if datasets:

        for dataset in datasets:
            if dataset.id == datasets[-1].id:
                dataset.dataset_filename = filename
                # dataset.filename is an attribute set only at run-time
                # to maintain the name of the VIFF file that was read in
                # rather than deriving a filename from the raw data
                # filenames with *.xml appended. But we need to set this
                # filename only for the primary dataset, not the associated
                # datasets. Associated datasets will default back to their
                # raw filenames if we go to save them for any reason
            else:
                dataset.dataset_filename = ''

    return datasets, timestamp

    


def main():

    verbose = True

    datatype = 'siemens dicom'

    # Processing of SVS_EDIT_OFF files
    STARTDIR = 'D:\\Users\\bsoher\\projects\\2017_oneil_ucla_lipid_contam_svs\\MRS_RAW_DATA\\_all_svs_edit_off' # \\fitted_pass1'
    csvfile         = STARTDIR+'\\bjs_csv_output_file_off.txt'

#     # Processing of SVS_EDIT_DIFF files
#     STARTDIR = 'D:\\Users\\bsoher\\projects\\2017_oneil_ucla_lipid_contam_svs\\MRS_RAW_DATA\\_all_svs_edit_diff\\fitted_pass1'
#     csvfile         = STARTDIR+'\\bjs_csv_output_file_diff.txt'


#     # this gets all files *.IMA in all subdirectories of STARTDIR
#     imafiles = []
#     for dirpath, dirnames, filenames in os.walk(STARTDIR):
#         for filename in [f for f in filenames if f.endswith(".xml")]:
#             imafiles.append(os.path.join(dirpath, filename))
#             print os.path.join(dirpath, filename)

    vespa_version = util_misc.get_vespa_version()

    i = 0

    imafiles = ['D:\\Users\\bsoher\\temp\\dmx\\test_5076.0004.0002.xml',]

    for datafile in imafiles:

        # Test input arguments for consistency --------------------------
        
        msg = ''
        if not os.path.isfile(datafile):
            msg = """Main DATAFILE does not exist "%s".""" % datafile 
            
        if msg:        
            print(msg, file=sys.stderr)
            print(msg, file=sys.stdout)
            sys.exit(-1)
            
        if not os.path.isfile(csvfile):
            if verbose:
                pass
                print("""Output CSV file will be created - %s""" % csvfile)
        
        
        # Load Main Dataset --------------------------
        
#        if verbose: print """%s - Load Data into a Dataset object - %s""" % (str(i), datafile)    
        dataset, timestamp = _open_viff(datafile)
        dataset = dataset[-1]
        print(str(i)+' : '+' - '+datafile)
          
        analysis_cli_oneil_results(dataset, csvfile, viffpath=datafile, 
                                                     vespa_version=vespa_version,
                                                     timestamp=timestamp,
                                                     verbose=verbose, 
                                                     debug=False)
        
        i += 1
        if i >= 1: break    # debug statement to exit after one file processed
        
    bob = 10
    bob += 1
    
    
    


        
if __name__ == '__main__':
    main()        











#     # no labels
#     nullfmt = NullFormatter()         # no labels
#
#     ax4.xaxis.set_major_formatter(nullfmt)  # have to go before plot()
#     ax4.yaxis.set_major_formatter(nullfmt)  # have to go before plot()
#     ax4.plot(freq[1000:1800])
#     plt.xticks(np.arange(100, 800, 200))
#     plt.xlim(0, 800)


#     # Gets rid of xaxis label for ax1 and ax2
#     xticklabels = ax1.get_xticklabels() + ax2.get_xticklabels()
#     plt.setp(xticklabels, visible=False)
    

#     # Annotate example with mathtext and how to set color with 3 term RGB vector
#     mpl_grey_rvb = (51./255., 51./255., 51./255.)
#     mpl_grey_rvb = (1./255., 1./255., 1./255.)
#     tmp2 =r"$\mathrm{Roman}\ , \ \mathit{Italic}\ , \ \mathtt{Typewriter} \, \ \mathrm{or}\ \mathcal{CALLIGRAPHY}$"
#     plt.annotate(tmp2,  xy=(110.0, 7.0),
#                         xycoords='data', 
#                         color=mpl_grey_rvb,
#                         fontsize=7)        