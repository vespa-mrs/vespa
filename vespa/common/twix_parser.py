"""
This is a simple Python parser for Siemens twix files created using version
VB software. As of VERSION 0.2.0, I am only maintaining this for Python 3.6 and
later. It contains two classes and one function. The classes are TwixMdh and
TwixScan. The former represents the measurement data header of a single scan
from a twix file. The latter is the data header plus the scan data.

The function read() accepts the name of a twix file and returns a list of
TwixScan instances and EVPs.

Brief sample usage is at the end of this file.

This file itself is executable; invoke like so:
   python twix.py  your_filename.dat


The class structure for twix objects was derived from many sources including
the ICE User's Guide, the Siemens presentation "ICE in IDEA for VD11A - News
on RAID, MDH and Tools", and Assaf Tal, Uzay Emir, and Jonathan Polimeni whose
code informed ours

  Classes and inline functions and other info about the MDH header was taken
  from the Siemens IDEA environment for software VB19 from:

  n4/pkg/MrServers/MrMeasSrv/SeqIF/MDH/mdh.h
  n4/pkg/MrServers/MrMeasSrv/SeqIF/MDH/MdhProxy.h

  Twix files layout are described (briefly!) on p138 of IceUsersGuide.pdf

Twix <= VB19 - Single RAID measurements per twix file
-------------------------------------------------------------------------------------------------

|  *  *  *  *  One Measurement  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  |
.                                                                                           .
.                                                                                           .
|  *  *  *  *   Header  *  *  *  *  |           Scan 1            | Scan 2 |  ...  | Scan N |
<---------    Header Size   --------> *                             *
.                                   .  *                              *
.                                   .   *                               *
| Header Size |   Header raw data   |   |    sMDH     | Measurement data |
<- uint32_t ->  (e.g. seqdata, evp)     <- 128 byte ->


Enjoy,
Brian Soher and Philip Semanchuk on behalf of the Vespa project
https://github.com/vespa-mrs/vespa

==============================
 NEWS
==============================

2012-06-27  Released version 0.1.0, developed and tested only under VB17
2017-01-23  Added more documentation, since it was already in the twix
             multi-raid parser file
2018-07-11  Bumped version to 0.1.1 since we've added the TwixRaid class
             to keep the API in line with the newer VE11 ++ class we had
             to write for MultiRaid twix files.
2020-08-18  Released version 0.2.0, this version has been updated for Python 3
             and has only been tested for Python 3.6 and later.
2020-09-12  Released version 0.3.0, Sort functionality moved from twix_sort.py
             into twix_parser.py and twix_parser_multi_raid.py modules. This
             adds Numpy as a dependency.
"""

# Python modules
import os
import sys
import time
import struct

# Third party modules
import numpy as np


VERSION = "0.3.0"



_MDH_ACQEND = 0x01              # mask for eval_info_mask first bit, indicates if last acquisition
_NDMA_FLAG_BYTES = 2            # There are 2 bytes of DMA flags
_NEVAL_INFO_FLAG_BYTES = 8      # There are 8 bytes of eval info flags
_NICE_PARAMETERS = 4            # There are 4 ICE params (2 bytes each = 8 bytes)
_NFREE_PARAMETERS = 4           # 4 free params for pulseq use (2 bytes each = 8 bytes)


SCAN_INDICES = ['ide', 'idd', 'idc', 'idb', 'ida',
                'segment', 'phase',
                'echo', 'partition',
                'slice', 'acquisition',
                'line', 'repetition',
                'set', 'channel_id' ]

MULTI_INDICES = ['ide', 'idd', 'idc', 'idb', 'ida',
                 'segment', 'phase',
                 'echo', 'partition',
                 'slice', 'acquisition',
                 'line', 'repetition', 'set' ]      # no channel_id index, inherent now


class TwixMdh(object):
    """ Object for a Twix MDH (measurement data header), contains no data. """

    def __init__(self):
        """ Inline comments are names of variables in Siemens' mdh.h """

        # Some code treats DMA length as a single 32-bit int and some treats
        # it as a 16-bit in plus two bytes of flags. I'm guessing the former
        # is the easygoing interpretation and that the latter is more correct.
        self.dma_length = 0                     # ulDMALength
        self.dma_flags = [0] * _NDMA_FLAG_BYTES
        self.measurement_userid = 0             # lMeasUID
        self.scan_count = 0                     # ulScanCounter
        self.timestamp = 0                      # ulTimeStamp - 2.5ms ticks since 00:00
        self.pmu_timestamp = 0                  # ulPMUTimeStamp - 2.5ms ticks since last trigger

        self.eval_info_mask = [0] * _NEVAL_INFO_FLAG_BYTES  # aulEvalInfoMask

        self.samples_in_scan = 0                # Number of complex data points.

        # Loop counters. These are 14 unsigned short values which index 
        # the measurement (in imaging; not used in spectroscopy).

        self.used_channels = 0                  # ushUsedChannels
        self.line = 0                           # ushLine
        self.acquisition = 0                    # ushAcquisition
        self.slice = 0                          # ushSlice
        self.partition = 0                      # ushPartition
        self.echo = 0                           # ushEcho
        self.phase = 0                          # ushPhase
        self.repetition = 0                     # ushRepetition
        self.set = 0                            # ushSet
        self.segment = 0                        # ushSeg

        self.ida = 0                            # ushIda - these 5 are ICE dimension indices
        self.idb = 0                            # ushIdb
        self.idc = 0                            # ushIdc
        self.idd = 0                            # ushIdd
        self.ide = 0                            # ushIde

        self.pre = 0                            # ushPre - Cut-off data
        self.post = 0                           # ushPost - Cut-off data

        self.kspace_center_column = 0           # ushKSpaceCentreColumn
        self.coil_select = 0                    # ushCoilSelect
        self.readout_off_center = 0.0           # fReadOutOffcentre
        self.time_since_last_rf = 0             # ulTimeSinceLastRF - timestamp since last RF pulse
        self.kspace_center_line = 0             # ushKSpaceCentreLineNo - number of K-space center line
        self.kspace_center_partition = 0        # ushKSpaceCentrePartitionNo

        # 4 ushorts (2 bytes each = 8 total) are reserved for ICE program params
        # 4 ushorts (2 bytes each = 8 total) are free for pulse sequence use
        self.ice_parameters = [0] * _NICE_PARAMETERS
        self.free_parameters = [0] * _NFREE_PARAMETERS

        # Slice data
        self.saggital = 0.0                     # flSag
        self.coronal = 0.0                      # flCor
        self.transverse = 0.0                   # flTra
        self.quaternion = [0.0] * 4             # aflQuaternion

        # docs say, "Channel ID must be the last parameter". Unsure what's meant.
        self.channel_id = 0                     # ulChannelId

        # negative, absolute PTAB position in [0.1 mm]
        self.ptab_pos_neg = 0                   # ushPTABPosNeg


    @property
    def is_last_acquisition(self):
        """ True if scan is last acquisition (NB. it contains no FID data), else False """
        return bool(self.eval_info_mask[0] & _MDH_ACQEND)


    @property
    def clock_time(self):
        """
        Returns the scan's timestamp as a string representing the clock
        time when the scan was taken. The string is in ISO format and 
        includes microseconds -- HH:MM:SS,uSEC

        This function assumes that the scanner's timestamp is in UTC and 
        converts it to local time. That assumption might be incorrect. Caveat 
        emptor.

        self.timestamp is the # of 2.5ms ticks since midnight, but when is
        midnight? Is that in local time or UTC, since many Linux machines
        use UTC for their system clock?

        """
        # Convert to seconds since midnight.
        timestamp = self.timestamp * 2.5 / 10e4
        microseconds = int(round((timestamp - int(timestamp)) * 10e5))
        timestamp = time.localtime(timestamp)

        # "%H:%M:%S" is ISO "extended" time format (HH:MM:SS)
        # http://en.wikipedia.org/wiki/ISO_8601#Times
        timestamp = time.strftime("%H:%M:%S", timestamp)

        # According to Wikipedia, "Decimal fractions may also be added...
        # A decimal mark, either a comma or a dot...with a preference for a 
        # comma according to ISO 8601:2004) is used as a separator..."
        # http://en.wikipedia.org/wiki/ISO_8601#Times
        timestamp += (",%d" % microseconds) 

        return timestamp


    def __str__(self):

        lines = [ ]
        lines.append("DMA length:      %d (0x%x)" % (self.dma_length, 
                                                self.dma_length))

        lines.append("-- DMA Flags --")
        for i, byte_ in enumerate(self.dma_flags):
            bits = _bit_string(byte_, 8)
            lines.append("   byte[%d]:      0x%02x == %s" % (i, byte_, bits))

        lines.append("Meas userid:     %d" % self.measurement_userid)
        lines.append("Scan count:      %d" % self.scan_count)
        lines.append("Timestamp:       %d" % self.timestamp)
        lines.append("PMU timestamp:   %d" % self.pmu_timestamp)
        lines.append("Clock time:      %s" % self.clock_time)
        lines.append("-- Eval info mask --")
        for i, mask in enumerate(self.eval_info_mask):
            bits = _bit_string(mask, 8)
            lines.append("   mask[%d]:      0x%02x == %s" % (i, mask, bits))

        lines.append("Samples:         %d" % self.samples_in_scan)
        lines.append("Channels used:   %d" % self.used_channels)
        lines.append("Line:            %d" % self.line)
        lines.append("Acquisition:     %d" % self.acquisition)
        lines.append("Slice:           %d" % self.slice)
        lines.append("Partition:       %d" % self.partition)
        lines.append("Echo:            %d" % self.echo)
        lines.append("Phase:           %d" % self.phase)
        lines.append("Repetition:      %d" % self.repetition)
        lines.append("Set:             %d" % self.set)
        lines.append("Segment:         %d" % self.segment)

        lines.append("ICE dim A:       %d" % self.ida)
        lines.append("ICE dim B:       %d" % self.idb)
        lines.append("ICE dim C:       %d" % self.idc)
        lines.append("ICE dim D:       %d" % self.idd)
        lines.append("ICE dim E:       %d" % self.ide)

        lines.append("Pre:             %d" % self.pre)
        lines.append("Post:            %d" % self.post)

        lines.append("K space center column: %d" % self.coil_select)
        lines.append("Coil select:           %d" % self.kspace_center_column)
        lines.append("Readout off center:    %f" % self.readout_off_center)
        lines.append("Time since last RF:    %d" % self.time_since_last_rf)
        lines.append("K space center line:   %d" % self.kspace_center_line)
        lines.append("K space center part:   %d" % self.kspace_center_partition)

        lines.append("-- ICE parameters --")
        for i, parameter in enumerate(self.ice_parameters):
            lines.append("   param[%d]:  %d (0x%x)" % (i, parameter, parameter))

        lines.append("-- Free parameters --")
        for i, parameter in enumerate(self.free_parameters):
            lines.append("   param[%d]:  %d (0x%x)" % (i, parameter, parameter))

        lines.append("Saggital:      %f" % self.saggital)
        lines.append("Coronal:       %f" % self.coronal)
        lines.append("Transverse:    %f" % self.transverse)
        lines.append("Quaternion:    %f %f %f %f" % tuple(self.quaternion))

        lines.append("Channel ID:    %d " % self.channel_id)
        lines.append("PTAB Position: %d " % self.ptab_pos_neg)

        return '\n'.join(lines)


    def populate_from_file(self, infile):
        """
        Given an open file or file-like object (like a StringIO instance) 
        that's positioned at the first byte of an MDH, populates this TwixMdh 
        object from the file. The file pointer is advanced to the end of the 
        header.

        """
        self.dma_length = _read_ushort(infile)
        self.dma_flags = [_read_byte(infile) for i in range(_NDMA_FLAG_BYTES)]

        self.measurement_userid = _read_int(infile)
        self.scan_count         = _read_uint(infile)
        self.timestamp          = _read_uint(infile)
        self.pmu_timestamp      = _read_uint(infile)

        self.eval_info_mask = [_read_byte(infile) for i in range(_NEVAL_INFO_FLAG_BYTES)]

        self.samples_in_scan = _read_ushort(infile) 
        self.used_channels  = _read_ushort(infile)
        self.line           = _read_ushort(infile)
        self.acquisition    = _read_ushort(infile)
        self.slice          = _read_ushort(infile)
        self.partition      = _read_ushort(infile)
        self.echo           = _read_ushort(infile)
        self.phase          = _read_ushort(infile)
        self.repetition     = _read_ushort(infile)
        self.set            = _read_ushort(infile)
        self.segment        = _read_ushort(infile)
        # ICE dimension indices
        self.ida = _read_ushort(infile)
        self.idb = _read_ushort(infile)
        self.idc = _read_ushort(infile)
        self.idd = _read_ushort(infile)
        self.ide = _read_ushort(infile)

        # Cut-off data
        self.pre  = _read_ushort(infile)
        self.post = _read_ushort(infile)

        self.kspace_center_column   = _read_ushort(infile)
        self.coil_select            = _read_ushort(infile)
        self.readout_off_center     = _read_float(infile)
        self.time_since_last_rf     = _read_uint(infile)
        self.kspace_center_line     = _read_ushort(infile)
        self.kspace_center_partition = _read_ushort(infile)

        self.ice_parameters  = [_read_ushort(infile) for i in range(_NICE_PARAMETERS)]
        self.free_parameters = [_read_ushort(infile) for i in range(_NFREE_PARAMETERS)]

        # Slice data
        self.saggital   = _read_float(infile)
        self.coronal    = _read_float(infile)
        self.transverse = _read_float(infile)
        self.quaternion = [_read_float(infile) for i in range(4)]

        self.channel_id = _read_ushort(infile)
        self.ptab_pos_neg = _read_ushort(infile)


class TwixScan(TwixMdh):
    """
    This is the same as a TwixMdh object, but with the addition of the 
    actual scan data (stored as an ordinary Python list).

    """
    def __init__(self):
        TwixMdh.__init__(self)
        self.data = []


    def __str__(self):
        s = TwixMdh.__str__(self)
        return ("Data points:     %d\n" % len(self.data)) + s


    def populate_from_file(self, infile, data_file=None):
        """
        Given an open file or file-like object (like a StringIO instance) 
        that's positioned at the first byte of an MDH, populates this TwixScan 
        object's data and metadata from the file. The file pointer is advanced 
        to the end of the header.

        By default, the data is read from the same file as the MDH. When the 
        optional data_file parameter is supplied, this method reads the data
        from that file instead. This is useful when the MDH and data are in 
        separate files.

        """
        TwixMdh.populate_from_file(self, infile)

        if not data_file:
            data_file = infile

        # Data is in complex #s that are stored as (real, imag) pairs.
        data = _read_float(data_file, 2 * self.samples_in_scan)

        self.data = _collapse_complexes(data)



class TwixRaid(object):
    """
    Twix files through VB software only copied a single measurement. From VD11
    on, Twix Multi-RAID expanded this constraint.  So, this module does not need
    the complexity of Multi-RAID. 
    
    Normally, we would just recommend the built-in read() method of this module
    but we provide the TwixRaid class to have a similar functionality to
    TwixMultiRaid to ease the plight of our users and provide different access 
    to the underlying data from the scans and headers.
    
    """
    def __init__(self):
    
        self.header_size     = 0
        self.evps            = None
        self.scans           = None
        self.free_parameters = []
        self.ice_parameters  = []

        # the following are Sort functionality attributes

        self.indices_list    = []       # fill with self.create_ice_indices()
        self.indices_unique  = []       # fill with self.check_unique_indices()
        self.dims            = []       # fill with self.get_dims()


    @property
    def current(self):
        return self         # no-op here to match VD/VE call in TwixMultiRaid


    ########  General methods  ################################################

    def __str__(self):
        s = self.scans.__unicode__(self)
        return  s


    def get_free_parameters(self):
        free_params = []
        for scan in self.scans:
            free_params.append(scan.free_parameters)
        return free_params


    def get_ice_parameters(self):
        ice_params = []
        for scan in self.scans:
            ice_params.append(scan.ice_parameters)
        return ice_params


    def populate_from_file(self, filename):
        """
        Given the name of a Twix file, reads the global header and all the
        scans in the file. 
    
        Returns a two-tuple of (scans, EVPs) where scans is a list of TwixScan 
        objects (one for each scan in the file) and EVPs is a list of 2-tuples.
        The EVP two-tuples are (name, data) and contain the name and data 
        associated with the ASCII EVP chunks in the global header. They're returned
        in the order in which they appeared. 
    
        To summarize, the return value is a 2-tuple containing two lists, all the 
        scans and all the evp headers:
        ( 
            [scan1, scan2, ... scanN], 
            [ (evp_name1, evp_data1), (evp_name2, evp_data2) (etc.) ]
        )
        """
        file_size = os.path.getsize(filename)
    
        infile = open(filename, 'rb')
    
        # Global header size is stashed in the first 4 bytes. This size includes
        # the size itself. In other words, the size can also be interpreted as 
        # an absolute offset from 0 that points to the the start of the first scan.
        # In fact, the doc says, "Skip data of this length for reading the real 
        # measurement data sorted as in former versions."
        header_size = _read_uint(infile)
    
        # Next is # of EVPs (whatever they are)
        nevps = _read_uint(infile)
    
        evps = [ ]
        while nevps:
            # Each EVP contains the name (as a NULL-terminated C string) followed by
            # the number of bytes (characters) of content followed by the content
            # itself. There's no accomodation made for multibyte characters or
            # anything other than ASCII, AFAICT. Or maybe since all the examples in
            # the manual are Windows-based, we should assume a character set 
            # of windows-1252?
            name = _read_cstr(infile)
            size = _read_uint(infile)
            data = _read_byte(infile, size)
            data = ''.join(map(chr, data))
            evps.append( (name, data) )
            nevps -= 1
    
        # There's a few bytes of padding between the EVPs and the start of the 
        # scans. Here we jump over that.
        infile.seek(header_size, os.SEEK_SET)       # start from beginning of file
    
        # Read the scans one by one until we hit the last (which should be flagged)
        # or run off the end of the file.
        scans = [ ]
        more_scans = True
        index = 0
        while more_scans:
            index += 1
            scan = TwixScan()
            scan.populate_from_file(infile)
    
            if scan.is_last_acquisition:
                more_scans = False
            else:
                scans.append(scan)
    
                # PS - I'm not sure if this is necessary. In the samples I have,
                # the scan.is_last_acquisition flag is set appropriately so I 
                # never run off the end of the file.
                if infile.tell() >= file_size:
                    more_scans = False
                    
        self.evps  = evps
        self.scans = scans
        self.free_parameters = self.get_free_parameters()
        self.ice_parameters  = self.get_ice_parameters()

        self.indices_list    = self.create_ice_indices()
        self.indices_unique  = self.check_unique_indices()
        self.dims            = self.get_dims()

            
    ######## Sort functionality methods #######################################
    #
    # On successful return from populate_from_file() we automatically create a
    # list of all ICE dimensions for each scan and run a check on these dims to
    # ensure that they are unique.
    #
    # If dims are unique, we can (if asked) return a numpy array will all 14
    # ICE dimensions and let the user make full use of the ICE loop counter
    # organization of the data. The array would look like:
    #
    # [ide, idd, idc, idb, ... , partition, slice, line, chan, col]
    #
    # If dims are unique, we can also return a numpy array in this order:
    #   [ repetition, channel_id, set (averages), spectral points ]
    # but we require that all other dims be equal to 1.
    #
    # If we don't care if dims unique or not, we can sort data into two types
    # of numpy arrays. The first is in [scan order, spectral points] and the
    # second is in [ channel_id, set (averages), spectral points] order. The
    # second version requires that the number of scans divides evenly into the
    # number of channels.

    def get_ice_index(self, iscan):
        """
        Return tuple of ICE loop values in SCAN_INDICES order given an index
        in the self.indices_list attribute. This method is typically used by
        the get_data_numpy() call, so we add an Ellipsis object to the end
        of the index list and convert the list into a tuple so it can be used
        to index a location in the numpy output array for the scan data.

        """
        if iscan > len(self.indices_list):
            msg = "Scan index outside range of indices_list"
            raise ValueError(msg)

        indx = list(self.indices_list[iscan])
        indx.append(Ellipsis)                   # simplify indexing in 14 dim arrays
        indx = tuple(indx)
        return indx

    def create_ice_indices(self):
        """ return list of the dims of each Scan """
        dims = []

        for scan in self.scans:
            vals = []
            for attr in SCAN_INDICES:
                vals.append(getattr(scan, attr))
            dims.append(vals)

        return dims

    def get_dims(self):
        """ return list of the max value in each dim, plus samples per scan """

        # transpose the list of lists and check max val in each dimension
        dims = [max(item) + 1 for item in zip(*self.indices_list)]

        dims.append(self.scans[0].samples_in_scan)

        return dims

    def check_unique_indices(self):
        """ uses set to determine if all ice dims are unique """
        tmp = [tuple(item) for item in self.indices_list]
        tmp = set(tmp)
        if len(tmp) == len(self.indices_list):
            return True
        else:
            return False

    def get_data_numpy_ice_indexed(self):
        """
        Use the ICE loop indices to sort each scan into a numpy array with
        all 16 dimensions.  This option can not be used if we don't have a
        unique set of ICE indices. We use check_unique_indices() to determine
        if this is the case.  If everything is OK, we return a numpy array,
        if it is not, we return None.

        """
        if not self.indices_unique:
            msg = "ICE indices are not unique, can not return in a Loop Counter dimensioned numpy array."
            raise ValueError(msg)

        nparr = np.zeros(self.dims, np.complex64)

        for i, scan in enumerate(self.scans):
            indx = self.get_ice_index(i)
            print("List, index = " + str(self.indices_list[i]) + "  " + str(indx))
            nparr[indx] = np.array(scan.data)

        return nparr

    def get_data_numpy_scan_order(self):
        """ return numpy array of all FID data in scan order - no indexing """

        nscan = len(self.scans)
        npts = self.scans[0].samples_in_scan

        nparr = np.zeros([nscan, npts], np.complex64)

        nparr = np.zeros([len(self.scans), self.scans[0].samples_in_scan], np.complex64)

        for i, scan in enumerate(self.scans):
            nparr[i, :] = np.array(scan.data)

        return nparr

    def get_data_numpy_channel_scan_order(self):
        """
        Sort all scans into numpy array in [channel_id, fid, spectral points] order.

        This does not require unique dimensions. It does require that total
        number of scans divided by number of channels is an integer. It is
        assumed that all channels for a given FID follow one after the other
        in the scan list.

        """
        npts = self.scans[0].samples_in_scan
        nscan = len(self.scans)
        ncha = self.scans[0].used_channels

        if (nscan % ncha) > 0:
            msg = "Number of scans not divisible by number of channels, inconsistent array size, returning.."
            raise ValueError(msg)

        nfid = int(nscan / ncha)

        nparr = np.zeros([ncha, nfid, npts], np.complex64)

        for j in range(nfid):
            for i in range(ncha):
                scan = self.scans[i + j * ncha]
                icha = scan.channel_id
                nparr[icha, j, :] = np.array(scan.data)

        return nparr

    def get_data_numpy_scan_channel_order(self):
        """
        Sort all scans into numpy array in [fid, channel_id, spectral points] order.

        This does not require unique dimensions. It does require that total
        number of scans divided by number of channels is an integer. It is
        assumed that all channels for a given FID follow one after the other
        in the scan list.

        """
        npts = self.scans[0].samples_in_scan
        nscan = len(self.scans)
        ncha = self.scans[0].used_channels

        if (nscan % ncha) > 0:
            msg = "Number of scans not divisible by number of channels, inconsistent array size, returning.."
            raise ValueError(msg)

        nfid = int(nscan / ncha)

        nparr = np.zeros([ncha, nfid, npts], np.complex64)

        for j in range(nfid):
            for i in range(ncha):
                scan = self.scans[i + j * ncha]
                icha = scan.channel_id
                nparr[icha, j, :] = np.array(scan.data)

        return nparr


    def get_data_numpy_channel_scan_order_with_prep(self):
        """
        Sort all scans into numpy array in [fid, channel_id, spectral points] order.

        This does not require unique dimensions. It does require that total
        number of scans divided by number of channels is an integer. It is
        assumed that all channels for a given FID follow one after the other
        in the scan list.

        """
        npts = self.scans[0].samples_in_scan
        nscan = len(self.scans)
        ncha = self.scans[0].used_channels

        nrep = self.dims[SCAN_INDICES.index('repetition')]  # number of repetitions, may be 1
        nset = self.dims[SCAN_INDICES.index('set')]  # number of FIDs for spectroscopy
        ncha2 = self.dims[SCAN_INDICES.index('channel_id')]  # number of channels, may be 1


        if (nscan % ncha) > 0:
            msg = "Number of scans not divisible by number of channels, inconsistent array size, returning.."
            raise ValueError(msg)

        nfid = int(nscan / ncha)

        header = _parse_protocol_data(self.evps[3][1])
        nprep_hdr = int(header['sSpecPara.lPreparingScans']) if 'sSpecPara.lPreparingScans' in header.keys() else 0
        nprep = nfid - nrep * nset
        if nprep_hdr != nprep: print('Warning: get_data_numpy_scan_channel_order_with_prep() - nprep_hdr != nprep')

        # I have seen in CMRR sLASER data that the scan.set value does not increment
        # while the prep scans are taken. We get around this empirically. We know
        # how many Averages are being taken, multiply by NCha and subtract from the
        # total number of scans and divide by NCha and we know how many prep averages
        # there will be.  I calc the iset index below using this sort of math. Just
        # in case someday it start to increment, I make allowances for that case too.

        prep_arr = None
        if nprep != 0:
            prep_arr = np.zeros([ncha, nprep, npts], np.complex64)
            for iscan, scan in enumerate(self.scans[0:int(nprep * ncha)]):
                ifid = iscan // ncha
                icha = scan.channel_id
                prep_arr[icha, ifid, :] = np.array(scan.data)

        nparr = np.zeros([ncha, nfid-nprep, npts], np.complex64)

        for iscan, scan in enumerate(self.scans[int(nprep * ncha):]):
            ifid = iscan // ncha
            icha = scan.channel_id
            nparr[icha, ifid, :] = np.array(scan.data)

        return nparr, prep_arr


    def get_data_numpy_rep_coil_avg_npts_order(self, return_prep=False):
        return self.get_data_numpy_rep_channel_set_order(return_prep=return_prep)

    def get_data_numpy_rep_channel_set_order(self, return_prep=False):
        """
        Use the ICE loop indices to sort each scan into a numpy array with
        all 16 dimensions.  This option can not be used if we don't have a
        unique set of ICE indices. We use check_unique_indices() to determine
        if this is the case.  If everything is OK, we return a numpy array,
        if it is not, we return None.

        """
        items = ['ide', 'idd', 'idc', 'idb', 'ida', 'segment', 'phase', 'echo', 'partition', 'slice', 'acquisition',
                 'line']
        for item in items:
            indx = SCAN_INDICES.index(item)
            if self.dims[indx] > 1:
                msg = "ICE dimension - '" + item + "' is greater than 1, can not create a unique numpy array."
                raise ValueError(msg)

        nrep = self.dims[SCAN_INDICES.index('repetition')]  # number of repetitions, may be 1
        nset = self.dims[SCAN_INDICES.index('set')]  # number of FIDs for spectroscopy
        ncha = self.dims[SCAN_INDICES.index('channel_id')]  # number of channels, may be 1
        npts = self.scans[0].samples_in_scan

        nparr = np.zeros([nrep, ncha, nset, npts], np.complex64)

        prep_arr = None
        if return_prep:
            header = _parse_protocol_data(self.evps[3][1])
            nprep_hdr = int(header['sSpecPara.lPreparingScans']) if 'sSpecPara.lPreparingScans' in header.keys() else 0
            nprep = int(len(self.scans) / ncha) - nrep * nset
            if nprep_hdr != nprep: print('Warning: get_data_numpy_rep_channel_set_order() - nprep_hdr != nprep')

            # I have seen in CMRR sLASER data that the scan.set value does not increment
            # while the prep scans are taken. We get around this empirically. We know
            # how many Averages are being taken, multiply by NCha and subtract from the
            # total number of scans and divide by NCha and we know how many prep averages
            # there will be.  I calc the iset index below using this sort of math. Just
            # in case someday it start to increment, I make allowances for that case too.

            prep_arr = None
            if nprep != 0:
                prep_arr = np.zeros([nrep, ncha, nprep, npts], np.complex64)
                for iscan, scan in enumerate(self.scans[0:int(nprep * ncha)]):
                    irep = scan.repetition
                    if scan.set == 0:
                        iset = iscan // ncha
                    else:
                        iset = scan.set
                    # iset = scan.set
                    icha = scan.channel_id
                    prep_arr[irep, icha, iset, :] = np.array(scan.data)
        else:
            nprep = 0

        for scan in self.scans[int(nprep * ncha):]:
            irep = scan.repetition
            iset = scan.set
            icha = scan.channel_id
            nparr[irep, icha, iset, :] = np.array(scan.data)

        return nparr, prep_arr

    def get_data_numpy_echo_coil_avg_npts_order(self, return_prep=False):
        return self.get_data_numpy_echo_channel_set_order(return_prep=return_prep)

    def get_data_numpy_echo_channel_set_order(self, return_prep=False):
        """
        Use the ICE loop indices to sort each scan into a numpy array with
        all 16 dimensions.  This option can not be used if we don't have a
        unique set of ICE indices. We use check_unique_indices() to determine
        if this is the case.  If everything is OK, we return a numpy array,
        if it is not, we return None.

        """
        items = ['ide', 'idd', 'idc', 'idb', 'ida', 'segment', 'phase', 'repetition', 'partition', 'slice',
                 'acquisition', 'line']
        for item in items:
            indx = SCAN_INDICES.index(item)
            if self.dims[indx] > 1:
                msg = "ICE dimension - '" + item + "' is greater than 1, can not create a unique numpy array."
                raise ValueError(msg)

        nrep = self.dims[SCAN_INDICES.index('echo')]  # number of repetitions, may be 1
        nset = self.dims[SCAN_INDICES.index('set')]  # number of FIDs for spectroscopy
        ncha = self.dims[SCAN_INDICES.index('channel_id')]  # number of channels, may be 1
        npts = self.scans[0].samples_in_scan

        nparr = np.zeros([nrep, ncha, nset, npts], np.complex64)

        prep_arr = None
        if return_prep:
            header = _parse_protocol_data(self.evps[3][1])
            nprep_hdr = int(header['sSpecPara.lPreparingScans']) if 'sSpecPara.lPreparingScans' in header.keys() else 0
            nprep = int(len(self.scans) / ncha) - nrep * nset
            if nprep_hdr != nprep: print('Warning: get_data_numpy_rep_channel_set_order() - nprep_hdr != nprep')

            prep_arr = None
            if nprep != 0:
                prep_arr = np.zeros([nrep, ncha, nprep, npts], np.complex64)
                for iscan, scan in enumerate(self.scans[0:int(nprep * ncha)]):
                    irep = scan.echo
                    if scan.set == 0:
                        iset = iscan // ncha
                    else:
                        iset = scan.set
                    # iset = scan.set
                    icha = scan.channel_id
                    prep_arr[irep, icha, iset, :] = np.array(scan.data)
        else:
            nprep = 0

        for scan in self.scans[int(nprep * ncha):]:
            irep = scan.echo
            iset = scan.set
            icha = scan.channel_id
            nparr[irep, icha, iset, :] = np.array(scan.data)

        return nparr, prep_arr











##################   Public functions start here   ##################   



def read(filename):
    """
    Moved this functionality into the TwixRaid class so I don't have
    duplicate code floating around.  But I leave this module level function
    in here for backwards compatibility.
    
    """
    measurement = TwixRaid()
    measurement.populate_from_file(filename)
    
    return measurement.scans, measurement.evps
    




# Python 2.5 lacks the bin() builtin which we use when printing a scan as
# a string. Here we cook up a reasonable fascimile. 
try:
    bin(42)
except NameError:
    # Apparently we are using Python 2.5
    
    # _BINARY_MAP maps hex characters to their equivalent in bits
    _BINARY_MAP = { '0' : '0000', '1' : '0001',
                    '2' : '0010', '3' : '0011',
                    '4' : '0100', '5' : '0101',
                    '6' : '0110', '7' : '0111',
                    '8' : '1000', '9' : '1001',
                    'a' : '1010', 'b' : '1011',
                    'c' : '1100', 'd' : '1101',
                    'e' : '1110', 'f' : '1111',
                  }

    def bin(value):
        """Given an int, returns a string representing the bits in the int.
        Mimics the bin() builtin of Python >= 2.6, except that this function
        only handles non-negative numbers.
        """
        # Swiped from here:
        # http://stackoverflow.com/questions/1993834/how-change-int-to-binary-on-python-2-5

        # Convert to a hex string
        value = '%x' % value

        # For each character in the string, substitute the binary value
        value = ''.join([_BINARY_MAP[c] for c in value])

        # Trim left 0s  
        value = value.lstrip('0')

        # Prefix with '0b' for compatibility with bin()
        return '0b' + value 


##################   "Private" functions start here   ##################   

def _bit_string(value, min_length=1):
    """
    Given an int value, returns a bitwise string representation that is 
    at least min_length long. 

    For example, read._bit_string(42, 8) returns '0b00101010'.
    """
    # Get value and trim leading '0b'
    value = bin(value)[2:]
    # Pad 
    value = value.rjust(min_length, '0')
    return '0b' + value


# The _read_xxx() functions are for reading specific types out of a file. All
# call _read_generic(). 

def _read_generic(source_file, type_, count=1):
    format = '<%d%s' % (count, type_)

    data = source_file.read(struct.calcsize(format))

    values = struct.unpack(format, data)
    if count == 1:
        return values[0]
    else:
        return values

def _read_byte(source_file, count=1):
    return _read_generic(source_file, 'B', count)

def _read_ushort(source_file, count=1):
    return _read_generic(source_file, 'H', count)

def _read_float(source_file, count=1):
    return _read_generic(source_file, 'f', count)

def _read_uint(source_file, count=1):
    return _read_generic(source_file, 'I', count)

def _read_int(source_file, count=1):
    return _read_generic(source_file, 'i', count)

def _read_cstr(source_file):
    # Reads until it hits a C NULL (0x00), returns the string (without the NULL)
    chars = [ ]
    char = source_file.read(1).decode()
    while char != '\x00':
        chars.append(char)
        char = source_file.read(1).decode()

    return ''.join(chars)


def _collapse_complexes(data):
    """Given a list or other iterable that's a series of (real, imaginary)
    pairs, returns a list of complex numbers. For instance, given this list --
       [a, b, c, d, e, f]
    this function returns --
       [complex(a, b), complex(c, d), complex(e, f)]

    The returned list is a new list; the original is unchanged.

    """
    # This code was chosen for speed and efficiency. It creates an iterator
    # over the original list which gets called by izip. (izip() is the same
    # as the builtin zip() except that it returns elements one by one instead
    # of creating the whole list in memory.)
    # It's the fastest method of the 5 or 6 I tried, and I think it is also
    # very memory-efficient. 
    # I stole it from here:
    # http://stackoverflow.com/questions/4628290/pairs-from-single-list
    data_iter = iter(data)

    return [complex(r, i) for r, i in zip(data_iter, data_iter)]


def _create_64bit_mask(vals):
    """ create a 64bit int value that corresponds to the bytes of the mask """
    if len(vals) != 8:
        # simple error check for number of vals to convert
        return 0
    
    bmask = 0
    for i, val in enumerate(vals):
        bmask += (val << i*8)
    return bmask


def _parse_protocol_data(prot):
    str, end = prot.find("### ASCCONV BEGIN"), prot.find("### ASCCONV END ###")
    prot  = prot[str+len("### ASCCONV BEGIN ###"):end]
    f = lambda pair: (pair[0].strip(), pair[1].strip())
    lines = prot.split('\n')
    return dict([f(line.split('=')) for line in lines[1:] if line])


#---------------------------------------------------------------------------------------------

if __name__ == "__main__":

    fname1 = 'D:\\bsoher\\code\\repository_svn\\sample_data\\siemens_twix_vd13\\meas_MID00099_FID06713_eja_svs_press.dat'
    fname2 = 'D:\\bsoher\\code\\repository_svn\\sample_data\\example_siemens_twix_wbnaa\\Sarah_WBNAA2.dat'
    
    twix = TwixRaid()
    twix.populate_from_file(fname1)

    # Uncomment the two lines below to print each scan to stdout.
    # for scan in twix.scans:
    #     print scan

    # Uncomment the two lines below to write the EVP data to files.
    # for name, data in twix.current.evps:
    #     open(name + ".txt", "wb").write(data)

    data11 = twix.get_data_numpy_ice_indexed()
    data12 = twix.get_data_numpy_rep_channel_scan_order()

