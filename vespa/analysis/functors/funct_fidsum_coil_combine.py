# Python imports


# 3rd party imports
import numpy as np

# Vespa imports
from vespa.analysis.algos.suspect_channel_combination import whiten, svd_weighting, combine_channels


COILCOMBINE_MENU_ITEMS = ['Siemens',
                          'CMRR',
                          'CMRR-Sequential',
                          'CMRR-Hybrid',
                          'SVD (suspect)',
                          'External Dataset',
                          'External Dataset with Offset']



def siemens(raw):
    """
    Combine code from IceSpectroEdit program, a la Mark Brown

    Given a list of scans and the extracted parameter dictionary we process
    the data in as similar a way to Siemens ICE program IceSpectroEdit as
    we can. My way of removing oversampling seems to not match 100% but the
    other steps are in the order and perform very similarly to the ICE program.

    Note. Since this got moved out of the Import step, we are assuming that
    the following processing steps have already been applied:
      1. oversampling removal
      2. global scaling to reasonable numerical range
      3. (optional) complex conjugate for proper display

    Input data is array with (1, ncoils, nfids, npts) dimensions
    Output data array has    (1,      1, nfids, npts) dimension

    Output ndarray collapses ncoils dimension by combining each group of ncoils
    FIDs into a weighted/phased summed FID.

    Other Output:

    all_weight   ndarray[nfids,ncoils], float, weights calculated/used to combine
                   each group of ncoil FIDs into a single FID

    all_phases   ndarray[nfids,ncoils], complex, phases calculated/used to combine
                   each group of ncoil FIDs into a single FID

    """
    nrep, ncoil, nfid, npts = raw.shape
    dat_comb = np.zeros([nfid, npts], dtype=np.complex)
    weight   = np.zeros([nfid, ncoil], dtype=np.float)
    phases   = np.zeros([nfid, ncoil], dtype=np.complex)

    flag_norm_to_sum = False  # default for now

    # Parse each group of FIDs for N channels for each average as separate
    # from all other scans. Perform the following steps:
    # - accumulate weighting factors and phases, all channels
    # - collate all channels as numpy arrays
    # - calc final weights for each channel
    # - apply weight and phase corrections to all channels

    for i in range(nfid):
        chans = []
        for j in range(ncoil):
            chan = raw[0, j, i, :].copy()
            weight[i,j] = np.abs(chan[0])
            phases[i,j] = np.conjugate(chan[0]) / weight[i,j]  # normalized complex conj to cancel phase
            chans.append(chan)

        # normalize weighting function based on spectro data

        tmp = np.sum([val * val for val in weight[i]])  # sum squared values
        if tmp == 0.0: tmp = 1.0
        if flag_norm_to_sum:
            lamda = np.sum(weight) / tmp        # sum of sensitivities
        else:
            lamda = 1.0 / np.sqrt(tmp)          # sqrt of sum of squared sensitivities

        weight[i] = [wt * lamda for wt in weight[i]]

        # apply weighting and phase corrections and sum
        for j, chan in enumerate(chans):
            dat_comb[i, :] += chan * weight[i,j] * phases[i,j]

    print_combine_stats(weight, phases, method='Siemens')

    return normalize_shape(dat_comb), weight, phases


def cmrr_standard(raw, delta=0.0):
    """
    Derived from Matlab code from Gulin Oz and Dinesh Deelchand, UMinn.
    Coil weights and phases are calculated once from the first scan,
    then these values are applied globally to all subsequent scans.

    Input data is array with (1, ncoil, nfid, npts) dimensions
    Output data array has    (1,     1, nfid, npts) dimension

    delta should be in radians (bjs - ??)

    Output array collapses ncoil dimension, combining ncoil FIDs into
    single eighted/phased summed FID.

    Other Output:
        all_weight   ndarray[nfid,ncoil], float, weights calculated/used to combine
                       each group of ncoil FIDs into a single FID
        all_phases   ndarray[nfid,ncoil], complex, phases calculated/used to combine
                       each group of ncoil FIDs into a single FID

    """
    nrep, ncoil, nfid, npts = raw.shape

    xaxis  = list(range(npts))
    flag_norm_to_sum = False  # default for now

    weight = np.zeros([ncoil,], dtype=np.float)
    phases = np.zeros([ncoil,], dtype=np.complex)


    # --------------------------------------------------------------------------
    # Calc weights and phases from first scan only
    # - weight is amplitude of zero order polynomial coefficient using 9th order
    #   polynomial fit of FID in time domain (based on  Uzay's script)

    for j in range(ncoil):
        chan = raw[0, j, 0, :].copy()
        weight[j] = np.abs(chan[0])
        phases[j] = np.conjugate(chan[0]) / weight[j]  # normalized complex conj to cancel phase
        coeffs = np.polyfit(xaxis, np.abs(chan), 9)
        weight[j] = coeffs[-1]                          # last entry is amplitude - zero order coeff

    # --------------------------------------------------------------------------
    # normalize weighting function based on spectro data

    tmp = np.sum([val * val for val in weight])  # sum squared values
    if tmp == 0.0: tmp = 1.0
    if flag_norm_to_sum:
        lamda = np.sum(weight) / tmp        # sum of sensitivities
    else:
        lamda = 1.0 / np.sqrt(tmp)          # sqrt of sum of squared sensitivities

    weight = np.array([val * lamda for val in weight])
    phases = phases + delta

    # --------------------------------------------------------------------------
    # Apply weights and phases from first scan to all scans

    all_weight = np.tile(weight, nfid)
    all_phases = np.tile(phases, nfid)
    all_weight.shape = (nfid, ncoil)
    all_phases.shape = (nfid, ncoil)

    dat_comb  = raw.copy()
    dat_comb *= weight.reshape(1, ncoil, 1, 1) * phases.reshape(1, ncoil, 1, 1)
    dat_comb  = dat_comb.sum(axis=(0,1))

    print_combine_stats(all_weight, all_phases, method='CMRR')

    return normalize_shape(dat_comb), all_weight, all_phases


def cmrr_sequential(raw):
    """ 
    Combine method hybrid of Siemens and Gulin Oz CMRR
    
    Derived from Matlab code from Dinesh. The coil weights and phases are 
    calculated from each scan in the scan list as in CMRR code, but are then
    applied only to its own scan as in the Siemens code.

    Input data is array with (1, ncoils, nfids, npts) dimensions
    Output data array has    (1,      1, nfids, npts) dimension

    Output ndarray collapses ncoils dimension by combining each group of ncoils
    FIDs into a weighted/phased summed FID. 
    
    Other Output:
    
    all_weight   ndarray[nfids,ncoils], float, weights calculated/used to combine 
                   each group of ncoil FIDs into a single FID 

    all_phases   ndarray[nfids,ncoils], complex, phases calculated/used to combine 
                   each group of ncoil FIDs into a single FID 
   
    """
    nrep, ncoil, nfid, npts = raw.shape
    dat_comb   = np.ndarray([nfid,npts], dtype=np.complex)
    all_weight = np.ndarray([nfid,ncoil], dtype=np.float)
    all_phases = np.ndarray([nfid,ncoil], dtype=np.complex)

    xaxis = list(range(npts))
    flag_norm_to_sum = False                    # default for now

    for i in range(nfid):

        # determine weighting and phz for each coil
        #   zero-order phase correction
        #   correct for phase based on 1st point in 1st wref fid

        # for each average, calc phase and weights to correct for coil geometry
        chans  = []
        weight = []
        phases = []
        
        for j in range(ncoil):
            chan = raw[0,j,i,:].copy()
            
            magn = np.abs(chan[0])
            phas = np.conjugate(chan[0])/magn        # normalized complex conj to cancel phase        
            chan = phas * chan                       # Note. applying phase here NOT below as in Siemens
        
            # amplitude of zero order phased fid in time domain
            #   using 9th order polynomial fit (based on  Uzay's script)
            coeffs = np.polyfit(xaxis, np.absolute(chan), 9)
            
            weight.append(coeffs[-1])       # last entry is amplitude - zero order coeff
            phases.append(phas)
            chans.append(chan)
        
        # normalize weighting function based on spectro data 
        tmp = np.sum([val*val for val in weight])   # sum squared values
        if tmp == 0.0: tmp = 1.0
        if flag_norm_to_sum:
            lamda = np.sum(weight) / tmp    # sum of sensitivities
        else:
            lamda = 1.0 / np.sqrt(tmp)      # sqrt of sum of squared sensitivities

        weight = [val*lamda for val in weight]

        all_weight[i,:] = weight
        all_phases[i,:] = phases
        
        # apply weighting ... phase corrections done above
        for j,chan in enumerate(chans):
            chans[j] = chan * weight[j]
        
        # sum corrected FIDs from each coil into one combined FID
        dat_comb[i,:] = np.sum(chans, axis=0) 

    print_combine_stats(all_weight, all_phases, method='CMRR_Sequential')

    return normalize_shape(dat_comb), all_weight, all_phases


def cmrr_hybrid(raw):
    raw_combined1, weights1, phases1 = siemens(raw)
    raw_combined2, weights2, phases2 = cmrr_standard(raw)

    phases1 = np.angle(np.sum(phases1, axis=0), deg=True)
    phases2 = np.angle(np.sum(phases2, axis=0), deg=True)

    a = phases2 - phases1
    vals = (a + 180) % 360 - 180
    delta = np.mean(vals) * np.pi / 180.0

    r = cmrr_standard(raw, delta=delta)
    return r


def svd_suspect(raw):
    """
    Based on Suspect (suspect.processing.channel_combination.svd_weighting)
    - change weights normalization to -> sqrt of sum of squared sensitivities
      to match up with other routines
    - adapted the input/ouput requirements to work for Analysis

    """
    nrep, ncoil, nfid, npts = raw.shape
    data = np.squeeze(np.mean(raw.copy(), axis=2))  # [ncoil,npts] shape

    u, s, v = np.linalg.svd(data, full_matrices=False)

    # we truncate SVD to rank 1, v[0] is our FID -> use v[0,0] to phase signal
    weights = u[:,0].conjugate()
    phases  = np.angle(v[0, 0])

    norm = np.sqrt(np.sum(np.array([np.abs(item)*np.abs(item) for item in weights])))
    norm_weights_phases = weights * np.exp(-1j * phases) / norm

    # --------------------------------------------------------------------------
    # Apply weights and phases from first scan to all scans

    dat_comb  = raw.copy()
    dat_comb *= norm_weights_phases.reshape(ncoil, 1, 1)
    dat_comb  = dat_comb.sum(axis=(0,1))

    all_weight = np.tile(np.abs(norm_weights_phases), nfid)
    all_phases = np.tile(np.angle(norm_weights_phases), nfid)
    all_weight.shape = (nfid, ncoil)
    all_phases.shape = (nfid, ncoil)

    return normalize_shape(dat_comb), all_weight, all_phases


def external_dataset(chain, delta=0.0):
    """ 
    Coil combine method that uses values calculated in another dataset.
    
    Input data is array with (1, ncoils, nfids, npts) dimensions
    Output data array has    (1,      1, nfids, npts) dimension

    Output ndarray collapses ncoils dimension by combining each group of ncoils
    FIDs into a weighted/phased summed FID. 
    
    Other Output:
    
    all_weight   ndarray[nfids,ncoils], float, weights calculated/used to combine 
                   each group of ncoil FIDs into a single FID 

    all_phases   ndarray[nfids,ncoils], complex, phases calculated/used to combine 
                   each group of ncoil FIDs into a single FID 
    
    In some cases we want to use the coil_combine results from a single high
    SNR water FID to combine multiple low SNR metabolite FIDs. In this case 
    we create weight and phase arrays that are nfids copies of the 1,ncoils 
    array of weight and phase values from the external data set.
    
    """
    nrep, ncoil, nfid, npts = chain.raw.shape
    dat_comb   = np.ndarray([nfid,npts],  dtype=np.complex)
    all_weight = chain.weights.copy()
    all_phases = chain.phases.copy() * np.exp(1j*delta)

    # expand/copy weights/phases if nfids dimension differs from raw data here
    if all_weight.shape[0] != nfid:
        all_weight = np.tile(all_weight[0], nfid)
        all_weight.shape = (nfid, ncoil)
    if all_phases.shape[0] != nfid:
        all_phases = np.tile(all_phases[0], nfid)
        all_phases.shape = (nfid, ncoil)

    for i in range(nfid):
        # apply pre-calc phase and weights to correct for coil geometry
        chans  = []
        for j in range(ncoil):
            data    = chain.raw[0,j,i,:].copy()
            weights = all_weight[i,j]
            phases  = all_phases[i,j]
            chans.append( data * weights * phases ) 
        
        # sum corrected FIDs from each coil into one combined FID
        dat_comb[i,:] = np.sum(chans, axis=0) 

    print_combine_stats(all_weight, all_phases, method='External Dataset')
        
    return normalize_shape(dat_comb), all_weight, all_phases 


def external_dataset_with_offset(chain):
    raw_combined1, weights1, phases10 = external_dataset(chain)
    raw_combined2, weights2, phases20 = siemens(chain)

    phases1 = np.angle(np.sum(phases10, axis=0), deg=True)
    phases2 = np.angle(np.sum(phases20, axis=0), deg=True)

    # find minimum angle between the two methods, note that depending on
    # lead/lag and whether the angles span the dislocation at 0/359 degrees.
    a = phases1 - phases2
    vals = (a + 180) % 360 - 180

    # get rid of outliers if we have enough coils to do so
    # - this keeps a few wrong numbers from skewing the mean offset
    ncoil = len(vals)
    if ncoil >= 32:
        nout = 4
    elif ncoil >= 8:
        nout = 2
    else:
        nout = 0

    if nout:
        cmean = np.mean(vals)
        cdiff = np.abs(cmean - vals)
        for n in range(nout):
            indx = np.argmax(cdiff)
            cdiff = np.delete(cdiff, indx)
            vals = np.delete(vals, indx)

    # then we take the mean value and use it as the overall offset
    delta = np.mean(vals)
    delta = delta * np.pi / 180.0

    r = external_dataset(chain, delta=-delta)

    return r


def coil_combine_none(chain):
    """ here we just copy data from first coil channel """

    dat_comb = chain.raw[0,0,:,:].copy()
    return normalize_shape(dat_comb)


def normalize_shape(data):
    """returns 'data' with 4 dimensional shape array (x, ncoils, nfids, npts) """
    while len(data.shape) < 4:
        data.shape = [1, ] + list(data.shape)
    return data




#------------------------------------------------------------------------------
# Helper Functions

def print_combine_stats(all_weight, all_phases, method=''):
    #pass
    return
    # take mean values for all fids on one coil 
    mean_weight = np.mean(all_weight, axis=0)
    mean_phases = np.angle(np.mean(all_phases, axis=0), deg=True) 
    stdv_weight = np.std(all_weight, axis=0)
    stdv_phases = np.std(np.angle(all_phases, deg=True), axis=0)

#    print '\n'
#    print method+" Weights mean      = \n", mean_weight
#    print method+" Weights ratios = \n", ["{0:0.2f}".format(i/max(abs(mean_weight))) for i in mean_weight]
#    print method+" Weights % stdv    = ", 100.0*stdv_weight/mean_weight
#    print method+" Phases [deg] mean = \n", mean_phases 
#    print method+" Phases [deg] mean = \n", ["{0:0.2f}".format(i) for i in mean_phases]
#    print method+" Phases % stdv     = ", 100.0*stdv_phases/mean_phases    



'''
function varargout = weightsmod(varargin)
% function varargout = weightsmod(varargin)
% determine and apply weighting Factor for each coil
% method = 1, determine weighting and phz for each coil
% method = 2, apply weights and phz to data
%
% Dinesh Deelchand, 18 March 2013

if (nargin < 2)
    method = 1;
    fidw = varargin{1};
elseif (nargin == 3)
    method = 2;
    fiduse = varargin{1};
    wk = varargin{2};
    refphz =  varargin{3};
else
    errorbox('Number of input parameters is not valid in weightsmod','warn')
    return
end

if (method == 1)
    %% determine weighting and phz for each coil
    % zero-order phase correction
    % correct for phase based on 1st point in 1st wref fid
    [np,nbCoils] = size(fidw);
    fidwphz = complex(zeros(size(fidw)));
    refphz = zeros(nbCoils,1);
    for iy = 1:nbCoils
        refphz(iy) = angle(fidw(1,iy));
        fidwphz(:,iy) = phase_adjust(fidw(:,iy),refphz(iy));
    end
    
    % weighting factor
    fidecc = squeeze(fidwphz(:,:));
    
    % amplitude of fid in time domain
    % using 9th order polynomial fit (based on  Uzay's script)
    warning off MATLAB:polyfit:RepeatedPointsOrRescale
    pts =(1:np)';
    amplfid = zeros(nbCoils,1);
    for ical = 1:nbCoils
        wktemp = polyfit(pts,abs(double(fidecc(:,ical))),9);
        amplfid(ical) = wktemp(end); % last entry is amplitude
        
        % show fit
        %f=polyval(wktemp,pts);
        %figure, plot(abs(double(fidecc(:,ical))));
        %hold on, plot(f,'r');
    end
    
    % weighting function
    wk = amplfid/sqrt(sum(amplfid.^2));
    
    % return variables
    varargout{1} = wk;
    varargout{2} = refphz;
    
elseif (method == 2)
    %% apply weights and phz to data
    if (ndims(fiduse) == 2) %2D size
        
        fiduse_phz = complex(zeros(size(fiduse)));
        nbCoils = size(fiduse,3);
        
        % zero-order phase correction
        for iy = 1:nbCoils
            fiduse_phz(:,iy) = phase_adjust(fiduse(:,iy),refphz(iy));
        end
        
        % apply weightings
        fidweighted = complex(zeros(size(fiduse_phz)));
        for ical=1:length(wk)
            fidweighted(:,ical) = wk(ical).*fiduse_phz(:,ical);
        end
        
        %sum all channels
        sumfidweighted = sum(fidweighted,2);
        
    elseif (ndims(fiduse) == 3)   %3D size
        
        fiduse_phz = complex(zeros(size(fiduse)));
        nbCoils = size(fiduse,3);
        nt = size(fiduse,2);
        
        % zero-order phase correction
        for iy = 1:nbCoils
            for ix = 1:nt
                fiduse_phz(:,ix,iy) = phase_adjust(fiduse(:,ix,iy),refphz(iy));
            end
        end
        
        % apply weightings
        fidweighted = complex(zeros(size(fiduse_phz)));
        for ical=1:length(wk)
            fidweighted(:,:,ical) = wk(ical).*fiduse_phz(:,:,ical);
        end
        
        %sum all channels
        sumfidweighted = sum(fidweighted,3);
    end
    
    % return variables
    varargout{1} = sumfidweighted ;
end
return


function outfid = phase_adjust(trc, phase)
% function to apply phase angles to complex data, e.g. to use
% variable receiver phases in VNMR

rtrc = real(trc);
itrc = imag(trc);

cosp = cos(phase); % real_cor
sinp = sin(phase); % imag_cor

rout = rtrc.*cosp + itrc.*sinp;
iout = itrc.*cosp - rtrc.*sinp;

outfid = complex(rout, iout);
return

'''


