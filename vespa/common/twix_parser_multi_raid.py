"""
This is a Python parser for Siemens Multi-RAID Twix files created using version
VD11+ and VE software (see twix_parser.py for VB software twix files). As of
VERSION 0.2.0, I am only maintaining this for Python 3.6 and later. It contains
seven classes, that correspond roughly to the layout of data in the file, and
various helper functions.

  The class structure for Mult-RAID twix objects was derived from the
  Siemens presentation PDF entitled "ICE in IDEA for VD11A - News on RAID,
  MDH and Tools"

  Classes and inline functions and other info about the MDH header was taken
  from the Siemens IDEA environment for software VE11 from:

  n4/pkg/MrServers/MrMeasSrv/SeqIF/MDH/mdh.h
  n4/pkg/MrServers/MrMeasSrv/SeqIF/MDH/MdhProxy.h


(Original) Twix <= VB19 - Single RAID measurements per twix file
-------------------------------------------------------------------------------------------------

|  *  *  *  *  One Measurement  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  |
.                                                                                           .
.                                                                                           .
|  *  *  *  *   Header  *  *  *  *  |           Scan 1            | Scan 2 |  ...  | Scan N |
<---------    Header Size   --------> *                             *
.                                   .  *                              *
.                                   .   *                               *
| Header Size |   Header raw data   |   |    sMDH    | Measurement data |
<- uint32_t ->  (e.g. seqdata, evp)     <- 128 byte ->



MultiTwix >= VD11 - Multi-RAID measurements per twix file
-------------------------------------------------------------------------------------------------

SuperRaidFileHeader

<--- measurement offset (8+152*64 = 9736 bytes  [10240 bytes total with align] ---><-504b-><- measurement length ->
|       MrParcRaidFileHeader        |   MrParcRaidFileEntry 1   | ... | MrParcRaidFileEntry N <=64 | Align |     Measurement 1     | Align | ... | Measurement N <=64 |
.                                   .*                                                               .
.                                   .  *                                                                                          .
.                   count <=64      .    *
| ID (hdSize_) | # of meas (count_) |    |  MeasID   |  FileID  | Measurement offset| Measurement length | Patient name | Protocol name |
<-- uint32_t --><---- uint32_t ----->    <-uint32_t-><-uint32_t-><--- uint64_t ----><---- uint64_t -----><-- 64 byte --><--- 64 byte --->
<-------------  8 byte ------------->    <----------------------- 152 byte ------------------------------------------------------------->

 - The last measurement contains the 'real' image scan data, first one's are adjustments, pre-scans, noise-scans
 - The last measurement depends on the preceeding measurements
 - Use offset to determine the position of a measurement in file
 - There are always 64 MrParcRaidFileEntry structs, indepenend of the actual number of measurements
 - After the SuperRaidFileHeader and after every measurement there is a 512 byte alignment filled with '0'
   - except as noted above there are only 504 bytes after MrParcRaidFileHeader/Entries

Measurement sub-structure

|  *  *  *  *  One Measurement  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * |
.                                                                                                                                                                                                      .
.                                                                                                                                                                                                      .
|  *  *  *  *   Header  *  *  *  *  |                    Scan 1, Cha 0                    |      Scan1, Cha1         | ... | Scan1,ChaM (max 128) | ... | Scan2,Cha0 | Scan2,Cha1 | ... | ScanN,ChaM |
<---------    Header Size   -------> *                                                    *                          *
.                                   . *                                                   *                          *
.                                   .  *   scan header     chan header                    *                          *
| Header Size |   Header raw data   |   |    scanH 0    |     chH 0    | Measurement data | chH 1 | Measurement data | ...
<- uint32_t ->  (e.g. seqdata, evp)     <-- 192 byte --><-- 32 byte -->
<-- structures same as orig twix -->    <-- new style structures different from orig twix   --------------> ...


The classes are TwixSuperRaidFileHeader, TwixMdh and TwixScan. The first represents the
layout of all measurements in the complete file, the next one will holds a
single measurement data header of a single scan from a twix file. The last is
the two previous plus the scan data.


Brief sample usage is at the end of this file.

This file itself is executable; invoke like so:
   python twix.py  your_filename.dat

Enjoy,
Brian Soher on behalf of the Vespa project
http://scion.duhs.duke.edu/vespa/

==============================
 NEWS
==============================

2017-01-27  Released version 0.1.0, developed and tested only under VE11
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


_NDMA_FLAG_BYTES = 2            # There are 2 bytes of DMA flags
_NEVAL_INFO_FLAG_BYTES = 8      # THere are 8 bytes of eval info flags
_NICE_PARAMETERS = 24           # There are 4 ICE params (2 bytes each = 8 bytes)
_NFREE_PARAMETERS = 4           # 4 free params for pulseq use (2 bytes each = 8 bytes)
                                #  NB. in VD, these are named ReservedPara, but we keep
                                #  same name here as for VB

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


MDH_FLAGS = { 0 : "MDH_ACQEND",                # last scan
              1 : "MDH_RTFEEDBACK",            # Realtime feedback scan
              2 : "MDH_HPFEEDBACK",            # High perfomance feedback scan
              3 : "MDH_ONLINE",                # processing should be done online
              4 : "MDH_OFFLINE",               # processing should be done offline
              5 : "MDH_SYNCDATA",              # readout contains synchroneous data
              
              8 : "MDH_LASTSCANINCONCAT",      # Flag for last scan in concatination

             10 : "MDH_RAWDATACORRECTION",     # Correct the rawadata with the rawdata correction factor
             11 : "MDH_LASTSCANINMEAS",        # Flag for last scan in measurement
             12 : "MDH_SCANSCALEFACTOR",       # Flag for scan specific additional scale factor
             13 : "MDH_2NDHADAMARPULSE",       # 2nd RF exitation of HADAMAR
             14 : "MDH_REFPHASESTABSCAN",      # reference phase stabilization scan
             15 : "MDH_PHASESTABSCAN",         # phase stabilization scan
             16 : "MDH_D3FFT",                 # execute 3D FFT
             17 : "MDH_SIGNREV",               # sign reversal
             18 : "MDH_PHASEFFT",              # execute phase fft
             19 : "MDH_SWAPPED",               # swapped phase/readout direction
             20 : "MDH_POSTSHAREDLINE",        # shared line
             21 : "MDH_PHASCOR",               # phase correction data
             22 : "MDH_PATREFSCAN",            # additonal scan for PAT reference line/partition
             23 : "MDH_PATREFANDIMASCAN",      # additonal scan for PAT reference line/partition that is also used as image scan
             24 : "MDH_REFLECT",               # reflect line
             25 : "MDH_NOISEADJSCAN",          # noise adjust scan
             26 : "MDH_SHARENOW",              # all lines are acquired from the actual and previous e.g. phases
             27 : "MDH_LASTMEASUREDLINE",      # indicates that the current line is the last measured line of all succeeding e.g. phases
             28 : "MDH_FIRSTSCANINSLICE",      # indicates first scan in slice (needed for time stamps)
             29 : "MDH_LASTSCANINSLICE",       # indicates  last scan in slice (needed for time stamps)
             30 : "MDH_TREFFECTIVEBEGIN",      # indicates the begin time stamp for TReff (triggered measurement)
             31 : "MDH_TREFFECTIVEEND",        # indicates the   end time stamp for TReff (triggered measurement)
             32 : "MDH_MDS_REF_POSITION",      # indicates the reference position for move during scan images (must be set once per slice/partition in MDS mode)
             33 : "MDH_SLC_AVERAGED",          # indicates avveraged slice for slice partial averaging scheme
             34 : "MDH_TAGFLAG1",              # adjust scan
             35 : "MDH_CT_NORMALIZE",          # Marks scans used to calculate correction maps for TimCT-Prescan normalize
             36 : "MDH_SCAN_FIRST",            # Marks the first scan of a particular map
             37 : "MDH_SCAN_LAST",             # Marks the last scan of a particular map
             40 : "MDH_FIRST_SCAN_IN_BLADE",   # Marks the first line of a blade
             41 : "MDH_LAST_SCAN_IN_BLADE",    # Marks the last line of a blade
             42 : "MDH_LAST_BLADE_IN_TR",      # Set for all lines of the last BLADE in each TR interval

             44 : "MDH_PACE",                  # Distinguishes PACE scans from non PACE scans.
             45 : "MDH_RETRO_LASTPHASE",       # Marks the last phase in a heartbeat
             46 : "MDH_RETRO_ENDOFMEAS",             # Marks an ADC at the end of the measurement
             47 : "MDH_RETRO_REPEATTHISHEARTBEAT",   # Repeat the current heartbeat when this bit is found
             48 : "MDH_RETRO_REPEATPREVHEARTBEAT",   # Repeat the previous heartbeat when this bit is found
             49 : "MDH_RETRO_ABORTSCANNOW",          # Just abort everything
             50 : "MDH_RETRO_LASTHEARTBEAT",   # This adc is from the last heartbeat (a dummy)
             51 : "MDH_RETRO_DUMMYSCAN",       # This adc is just a dummy scan, throw it away
             52 : "MDH_RETRO_ARRDETDISABLED",  # Disable all arrhythmia detection when this bit is found
             53 : "MDH_B1_CONTROLLOOP",        # Marks the readout as to be used for B1 Control Loop
             54 : "MDH_SKIP_ONLINE_PHASCOR",   # Marks scans not to be online phase corrected, even if online phase correction is switched on
             55 : "MDH_SKIP_REGRIDDING"        # Marks scans not to be regridded, even if regridding is switched on
            }




class TwixScanHeader(object):
    """ Object for a Twix MDH (measurement data header), contains no data. """

    def __init__(self):
        # Inline comments are (I think) the names of the corresponding 
        # variables in Siemens' mdh.h. 

        # Some code treats DMA length as a single 32-bit int and some treats
        # it as a 16-bit in plus two bytes of flags. I'm guessing the former
        # is the easygoing interpretation and that the latter is more correct.
        self.dma_length         = 0             # ulDMALength
        self.dma_flags          = [0] * _NDMA_FLAG_BYTES
        self.measurement_userid = 0             # lMeasUID
        self.scan_counter       = 0             # ulScanCounter
        self.timestamp          = 0             # ulTimeStamp - 2.5 ms ticks since 00:00
        self.pmu_timestamp      = 0             # ulPMUTimeStamp - 2.5 ms ticks since last trigger

        self.system_type        = 0             # ushSystemType
        self.ptab_pos_delay     = 0             # ulPTABPosDelay - NB. CODE HAS uint16_t in here not uint32_t
        self.ptab_posx          = 0             # lPTABPosX - absolute PTAB position in [um]
        self.ptab_posy          = 0             # lPTABPosY
        self.ptab_posz          = 0             # lPTABPosZ
        self.reserved1          = 0             # ulReserved1 - reserved for future hardware signals
        
        self.eval_info_mask = [0] * _NEVAL_INFO_FLAG_BYTES  # aulEvalInfoMask

        self.samples_in_scan = 0        # Number of complex data points.

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

        self.ida = 0                            # ushIda - ICE dimension indices
        self.idb = 0                            # ushIdb
        self.idc = 0                            # ushIdc
        self.idd = 0                            # ushIdd
        self.ide = 0                            # ushIde

        self.pre = 0                            # ushPre - Cut-off data
        self.post = 0                           # ushPost

        self.kspace_center_column = 0           # ushKSpaceCentreColumn
        self.coil_select = 0                    # ushCoilSelect
        self.readout_off_center = 0.0           # fReadOutOffcentre
        self.time_since_last_rf = 0             # ulTimeSinceLastRF - timestamp since last RF pulse
        self.kspace_center_line = 0             # ushKSpaceCentreLineNo - number of K-space center line
        self.kspace_center_partition = 0        # ushKSpaceCentrePartitionNo

        # Slice data
        self.saggital   = 0.0                   # flSag
        self.coronal    = 0.0                   # flCor
        self.transverse = 0.0                   # flTra
        self.quaternion = [0.0] * 4             # aflQuaternion

        # 24 ushorts (2 bytes each = 48 total) are reserved for ICE program params
        # 4 ushorts  (2 bytes each = 8 total) are free for pulse sequence use
        self.ice_parameters = [0] * _NICE_PARAMETERS
        self.free_parameters = [0] * _NFREE_PARAMETERS

        self.application_counter = 0            # ushApplicationCounter
        self.application_mask    = 0            # ushApplicationMask

        self.crc  = 0                           # ulCRC - 32 checksum


    @property
    def is_last_acquisition(self):
        """
        Returns True if this scan is the last acquisition (which means it
        contains no FID), False otherwise.
        
        """
        return bool(self.test_eval_info_by_label("MDH_ACQEND"))


    @property
    def clock_time(self):
        """
        Returns the scan's timestamp as a string representing the clock
        time when the scan was taken. The string is in ISO format and 
        includes microseconds -- HH:MM:SS,uSEC

        This function assumes that the scanner's timestamp is in UTC and 
        converts it to local time. That assumption might be incorrect. 
        Caveat emptor.
        
        """
        # self.timestamp is the # of 2.5ms ticks since midnight, but when is
        # midnight? Is that in local time or UTC, since many Linux machines
        # use UTC for their system clock? 

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
        return self.__unicode__()


    def __unicode__(self):
        lines = [ ]
        lines.append("DMA length:      %d (0x%x)" % (self.dma_length, self.dma_length))

        lines.append("-- DMA Flags --")
        for i, byte_ in enumerate(self.dma_flags):
            bits = _bit_string(byte_, 8)
            lines.append("   byte[%d]:      0x%02x == %s" % (i, byte_, bits))

        lines.append("Meas userid:     %d" % self.measurement_userid)
        lines.append("Scan counter:    %d" % self.scan_counter)
        lines.append("Timestamp:       %d" % self.timestamp)
        lines.append("PMU timestamp:   %d" % self.pmu_timestamp)
        lines.append("Clock time:      %s" % self.clock_time)
        
        lines.append("System type:     %d" % self.system_type)
        lines.append("PTAB pos delay:  %d" % self.ptab_pos_delay)
        lines.append("PTAB position x: %d" % self.ptab_posx)
        lines.append("PTAB position y: %d" % self.ptab_posy)
        lines.append("PTAB position z: %d" % self.ptab_posz)
        
        lines.append("-- Eval info mask --")
        for i, mask in enumerate(self.eval_info_mask):
            bits = _bit_string(mask, 8)
            lines.append("   mask[%d]:      0x%02x == %s" % (i, mask, bits))

        lines.append("   mask 64 bit:    == %s" % _create_64bit_mask(self.eval_info_mask))

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

        lines.append("Saggital:      %f" % self.saggital)
        lines.append("Coronal:       %f" % self.coronal)
        lines.append("Transverse:    %f" % self.transverse)
        lines.append("Quaternion:    %f %f %f %f" % tuple(self.quaternion))

        lines.append("-- ICE parameters --")
        for i, parameter in enumerate(self.ice_parameters):
            lines.append("   param[%d]:  %d (0x%x)" % (i, parameter, parameter))

        lines.append("-- Free parameters --")
        for i, parameter in enumerate(self.free_parameters):
            lines.append("   param[%d]:  %d (0x%x)" % (i, parameter, parameter))


        lines.append("Application counter: %d " % self.application_counter)
        lines.append("Application mask:    %d " % self.application_mask)

        lines.append("Checksum crc:        %d " % self.crc)

        return '\n'.join(lines)


    def populate_from_file(self, infile):
        """
        Given an open file or file-like object (like a StringIO instance) 
        that's positioned at the first byte of a scan header, populates this  
        object from the file. The file pointer is advanced to the end of the 
        header.
        """
        self.dma_length             = _read_ushort(infile)

        self.dma_flags = [_read_byte(infile) for i in range(_NDMA_FLAG_BYTES)]

        self.measurement_userid     = _read_int(infile)
        self.scan_counter           = _read_uint(infile)
        self.timestamp              = _read_uint(infile)
        self.pmu_timestamp          = _read_uint(infile)

        self.system_type            = _read_ushort(infile)
        self.ptab_pos_delay         = _read_ushort(infile)
        self.ptab_posx              = _read_int(infile)
        self.ptab_posy              = _read_int(infile)
        self.ptab_posz              = _read_int(infile)
        self.reserved1              = _read_uint(infile)

        self.eval_info_mask = [_read_byte(infile) for i in range(_NEVAL_INFO_FLAG_BYTES)]

        self.samples_in_scan        = _read_ushort(infile) 
        self.used_channels          = _read_ushort(infile)
        self.line                   = _read_ushort(infile)
        self.acquisition            = _read_ushort(infile)
        self.slice                  = _read_ushort(infile)
        self.partition              = _read_ushort(infile)
        self.echo                   = _read_ushort(infile)
        self.phase                  = _read_ushort(infile)
        self.repetition             = _read_ushort(infile)
        self.set                    = _read_ushort(infile)
        self.segment                = _read_ushort(infile)
        # ICE dimension indices
        self.ida                    = _read_ushort(infile)
        self.idb                    = _read_ushort(infile)
        self.idc                    = _read_ushort(infile)
        self.idd                    = _read_ushort(infile)
        self.ide                    = _read_ushort(infile)
        # Cut-off data
        self.pre                    = _read_ushort(infile)
        self.post                   = _read_ushort(infile)

        self.kspace_center_column   = _read_ushort(infile)
        self.coil_select            = _read_ushort(infile)
        self.readout_off_center     = _read_float(infile)
        self.time_since_last_rf     = _read_uint(infile)
        self.kspace_center_line     = _read_ushort(infile)
        self.kspace_center_partition = _read_ushort(infile)

        # Slice data
        self.saggital               = _read_float(infile)
        self.coronal                = _read_float(infile)
        self.transverse             = _read_float(infile)
        self.quaternion             = [_read_float(infile) for i in range(4)]

        self.ice_parameters  = [_read_ushort(infile) for i in range(_NICE_PARAMETERS)]
        self.free_parameters = [_read_ushort(infile) for i in range(_NFREE_PARAMETERS)]

        self.application_counter    = _read_ushort(infile)
        self.application_mask       = _read_ushort(infile)

        self.crc                    = _read_uint(infile)


    def parse_evalinfomask(self):
        """ 
        Returns list of labels for all states that are set in the mask.
        
        I create a 64 bit bitmask from the 8 bytes of the eval info mask, then I
        run through all the flag values listed in the MDH_FLAGS dictionary. Note
        that the dict keys are the bit location, so I have to create a bit shifted
        value to compare to the bit mask. The & operator returns 0 if the bit is
        not set, and non-zero if it is set. I add the flag label to the return
        list if the bit is set.
        
        """
        
        set_flag_labels = []
        
        bmask = _create_64bit_mask(self.eval_info_mask)
        
        for item in list(MDH_FLAGS.keys()):
            if (1<<item) & bmask:
                set_flag_labels.append(MDH_FLAGS[item])
            
        return set_flag_labels
    
    
    def test_eval_info_by_bit(self, bit):
        """ test bit location in 64bit eval info mask to see if set """
        if bit < 0 or bit > 63:
            # very light error check, should raise exception
            return False
        
        bmask = _create_64bit_mask(self.eval_info_mask)
        if (1<<bit) & bmask:
            return True
        else:
            return False


    def test_eval_info_by_label(self, label):
        """ test bit location in 64bit eval info mask to see if set """
        
        bit = -1
        for val,item in MDH_FLAGS.items():
            if label == item:
                bit = val
                break

        if bit == -1:
            # very light error check, should raise exception
            return False

        return self.test_eval_info_by_bit(bit)
        


class TwixChannelHeader(object):
    """ header (32 bytes) preceding each channel's data in a scan """

    def __init__(self):
        self.type_and_channel_length =  0   # uint32_t
        self.meas_uid               =   0   # int32_t
        self.scan_counter           =   0   # uint32_t
        self.reserved1              =   0   # uint32_t
        self.sequence_time          =   0   # uint32_t
        self.unused2                =   0   # uint32_t
        self.channel_id             =   0   # uint16_t
        self.unused3                =   0   # uint16_t
        self.crc                    =   0   # uint32_t


    def __str__(self):
        return self.__unicode__()


    def __unicode__(self):
        lines = []
        lines.append("-- Channel Header --")
        lines.append("Type and channel Length: %d" % self.measurement_userid)
        lines.append("MeasUID:      %d" % self.meas_uid)
        lines.append("ScanCounter:  %d" % self.scan_counter)
        lines.append("Reserved1:    %d" % self.reserved1)
        lines.append("SequenceTime: %d" % self.sequence_time)
        
        lines.append("Unused2:      %d" % self.unused2)
        lines.append("ChannelId:    %d" % self.channel_id)
        lines.append("Unused3:      %d" % self.unused3)
        lines.append("Checksum crc: %d" % self.crc)
        
        return '\n'.join(lines)


    def populate_from_file(self, infile):
        """
        Given an open file or file-like object (like a StringIO instance) 
        that's positioned at the first byte of a channel header, populates 
        object from the file. The file pointer is advanced to the end of the 
        header.
        
        """
        self.type_and_channel_length = _read_uint(infile)
        self.meas_uid                = _read_int(infile)
        self.scan_counter            = _read_uint(infile)
        self.reserved1               = _read_uint(infile)
        self.sequence_time           = _read_uint(infile)

        self.unused2                 = _read_uint(infile)
        self.channel_id              = _read_ushort(infile)
        self.unused3                 = _read_ushort(infile)
        self.crc                     = _read_uint(infile)


class TwixScan(object):
    """
    This contains the single TwixScanHeader and a list of one or more tuples 
    each containing (channel header, channel data) scan data where the data is
    stored as an ordinary Python list.
    
    """
    def __init__(self):
        self.scan_header    = None      # the global scan header
        self.scan_index     = 0         # keeps track of scan order in twix file
        self.channels       = []        # list of tuples [(chan hdr, chan data), ...]


    def __unicode__(self):
        s = self.scan_header.__unicode__(self)
        c = self.channels[0][1].__unicode__(self)
        
        return ("Data points:     %d\n" % len(self.channels[0][1])) + s + c


    def populate_from_file(self, infile):
        """
        Given an open file or file-like object (like a StringIO instance) 
        that's positioned at the first byte of an MDH, populates this TwixScan 
        object's data and metadata from the file. The file pointer is advanced 
        to the end of the header.

        """
        start_loc = infile.tell()
        
        # Make my base class do its work
        
        self.scan_header = TwixScanHeader()
        self.scan_header.populate_from_file(infile)

        flags = self.scan_header.parse_evalinfomask()
        
        if "MDH_SYNCDATA" in flags:
            # do not know how to process this type of scan, so skip for now
            infile.seek(start_loc + self.scan_header.dma_length)
        else:
            for i in range(self.scan_header.used_channels):
    
                chan_hdr = TwixChannelHeader()
                chan_hdr.populate_from_file(infile)
    
                # Data is in complex #s that are stored as (real, imag) pairs.
                data = _read_float(infile, 2 * self.scan_header.samples_in_scan)
                data = _collapse_complexes(data)
                
                self.channels.append( (chan_hdr, data) )


class TwixMrParcRaidFileEntry(object):
    
    def __init__(self):
        """ implementation of MrParcRaidFileEntry """
        
        self.meas_id                = 0         # uint32_t
        self.file_id                = 0         # uint32_t
        self.measurement_offset     = 0         # uint64_t
        self.measurement_length     = 0         # uint64_t
        self.patient_name           = ''        # 64 byte
        self.protocol_name          = ''        # 64 byte


class TwixSuperRaidFileHeader(object):
    """
    Represents a Twix SuperRaidFileHeader. Describes the layout of the one or
    more measurements within the meas.dat file.
  
    By definition, there are no more than 64 entries
    Total SuperRaidFileHeader length = 8b + 152*64b = 9736b + 504b = 10240b
    
    """
    def __init__(self):

        self.super_raid_id      = 0
        self.super_raid_count   = 0
        self.entries            = []

    def __str__(self):
        lines = [ ]
        lines.append("Super raid id:        %d" % self.super_raid_id)
        lines.append("Super raid count:     %d" % self.super_raid_count)

        for i in range(self.super_raid_count):
            item = self.entries[i]
            lines.append(" ")
            lines.append("-- MrParcRaidFileEntry %d --" % i)
            lines.append("MeasID:               %d" % item.meas_id)
            lines.append("FileID:               %d" % item.file_id)
            lines.append("Measurement offset:   %d" % item.measurement_offset)
            lines.append("Measurement length:   %d" % item.measurement_length)
            lines.append("Patient name:         %s" % item.patient_name)
            lines.append("Protocol name:        %s" % item.protocol_name)

        return '\n'.join(lines)

    def populate_from_file(self, infile):
        """Given an open file or file-like object (like a StringIO instance) 
        that's positioned at the first byte of an MDH, populates this TwixMdh 
        object from the file. The file pointer is advanced to the end of the 
        header.
        """
        self.entries = []
        
        start_parc = infile.tell()

        self.super_raid_id      = _read_uint(infile)
        self.super_raid_count   = _read_uint(infile)

        for i in range(self.super_raid_count):

            start_parc = infile.tell()
        
            item = TwixMrParcRaidFileEntry()
            item.meas_id                = _read_uint(infile)
            item.file_id                = _read_uint(infile)
            item.measurement_offset     = _read_ulonglong(infile)
            item.measurement_length     = _read_ulonglong(infile)
            infile.seek(start_parc + 24)
            item.patient_name           = _read_cstr(infile)
            infile.seek(start_parc + 88)
            item.protocol_name          = _read_cstr(infile)
            infile.seek(start_parc + 152)

            self.entries.append(item)


class TwixMeasurement(object):

    def __init__(self):
    
        self.header_size     = 0
        self.evps            = None
        self.scans           = None
        self.free_parameters = []
        self.ice_parameters  = []

        # Attributes for Sort functionality

        self.indices_list   = []        # fill with self.create_ice_indices()
        self.indices_unique = []        # fill with self.check_unique_indices()
        self.dims           = []        # fill with self.get_dims()


    def get_free_parameters(self):
        free_params = []
        for scan in self.scans:
            free_params.append(scan.scan_header.free_parameters)
        return free_params


    def get_ice_parameters(self):
        ice_params = []
        for scan in self.scans:
            ice_params.append(scan.scan_header.ice_parameters)
        return ice_params


    def populate_from_file(self, infile, measurement_length):
        """
        Given a file pointer at the first byte of a Measurement header, this
        method reads the header and all the scans in the file. 

        Populates the evps and scans attributes, where scans is a list of 
        TwixScan objects (one for each scan in the file) and EVPs is a list 
        of 2-tuples. The EVP two-tuples are (name, data) and contain the 
        name and data associated with the ASCII EVP chunks in the measurement 
        header. They're returned in the order in which they appeared. 

        To summarize the attribute layouts are:
        
            self.scans = [scan1, scan2, ... scanN]
            self.evps  = [ (name1, data1), (name2, data2) (etc.) ]
        
        """

        # assumption is that the file points to the first byte of the header
        measurement_start = infile.tell()

#        print 'hdr_start = '+str(infile.tell())

        # Global header size is stashed in the first 4 bytes. This size includes
        # the size itself. In other words, the size can also be interpreted as 
        # an absolute offset from 0 that points to the the start of the first scan.
        # In fact, the doc says, "Skip data of this length for reading the real 
        # measurement data sorted as in former versions."
        self.header_size = _read_uint(infile)

#        print 'hdr_start + hdr_size = '+str(start_loc + self.header_size)

        # Next is # of EVPs (whatever they are)
        nevps = _read_uint(infile)

        evps = []
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
        infile.seek(measurement_start + self.header_size, os.SEEK_SET)

#        print 'after all evps = '+str(infile.tell())

        # Read the scans one by one until we hit the last (which should be flagged)
        # or run off the end of the file.
        scans = []
        more_scans = True
        index = 0
        while more_scans:
            scan = TwixScan()
            scan.populate_from_file(infile)
            scan.scan_index = index

            if scan.scan_header.is_last_acquisition:
                more_scans = False
            else:
                scans.append(scan)

                # PS - I'm not sure if this is necessary. In the samples I have,
                # the scan.scan_header.is_last_acquisition flag is set appropriately so I 
                # never run off the end of the file.
                if (infile.tell() - measurement_start) >= measurement_length:
                    more_scans = False

            index += 1

        self.scans = scans
        self.evps  = evps
        self.free_parameters = self.get_free_parameters()
        self.ice_parameters = self.get_ice_parameters()

        self.indices_list   = self.create_ice_indices()
        self.indices_unique = self.check_unique_indices()
        self.dims           = self.get_dims()

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
    #   [ide, idd, idc, idb, ... , partition, slice, line, chan, col]
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
        Return list of ICE loop values in MULTI_INDICES order given an index
        in the self.indices_list attribute. This method is typically used by
        the get_data_numpy() call, so we add an Ellipsis object to the end
        of the index list and convert the list into a tuple so it can be used
        to index a location in the numpy output array for the scan data.

        """
        if iscan > len(self.indices_list):
            msg = "Scan index outside range of indices_list"
            raise ValueError(msg)

        indx = list(self.indices_list[iscan])
        return indx

    def create_ice_indices(self):
        """ return list of the dims of each Scan """
        dims = []

        for scan in self.scans:
            vals = []
            for attr in MULTI_INDICES:
                vals.append(getattr(scan.scan_header, attr))
            dims.append(vals)

        return dims

    def get_dims(self):
        """ return list of the max value in each dim, plus samples per scan """

        # transpose the list of lists and check max val in each dimension
        dims = [max(item) + 1 for item in zip(*self.indices_list)]

        dims.append(self.scans[0].scan_header.used_channels)
        dims.append(self.scans[0].scan_header.samples_in_scan)

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

        NB. channel_id is not always 0 .. Nchan, sometimes it is an odd
            assorment of integers, maybe the elements turned on, so we have
            to use an enumeration to fill the icha index

        """
        if not self.indices_unique:
            msg = "ICE indices are not unique, can not return in a Loop Counter dimensioned numpy array."
            raise ValueError(msg)

        nparr = np.zeros(self.dims, np.complex)

        for i, scan in enumerate(self.scans):
            loops = self.get_ice_index(i)

            for icha, chan in enumerate(scan.channels):
                loops.append(icha)  # chan[0].channel_id
                loops.append(Ellipsis)
                indx = tuple(loops)

                print("List, index = " + str(self.indices_list[i]) + "  " + str(indx))
                nparr[indx] = np.array(chan[1])

        return nparr

    def get_data_numpy_scan_order(self):
        """
        Return numpy array of all FID data in scan order - no indexing

        NB. channel_id is not always 0 .. Nchan, sometimes it is an odd
            assorment of integers, maybe the elements turned on, so we have
            to use an enumeration to fill the icha index

        """

        ncha = self.scans[0].scan_header.used_channels
        npts = self.scans[0].scan_header.samples_in_scan

        nparr = np.zeros([len(self.scans) * ncha, npts], np.complex)

        for i, scan in enumerate(self.scans):
            for j, chan in enumerate(scan.channels):
                nparr[i * ncha + j, :] = np.array(chan[1])

        return nparr

    def get_data_numpy_channel_scan_order(self):
        """
        Sort all scans into numpy array in [channel_id, fid, spectral points] order.

        This does not require unique dimensions. It does require that total
        number of scans divided by number of channels is an integer. It is
        assumed that all channels for a given FID follow one after the other
        in the scan list.

        NB. channel_id is not always 0 .. Nchan, sometimes it is an odd
            assorment of integers, maybe the elements turned on, so we have
            to use an enumeration to fill the icha index

        """
        npts = self.scans[0].scan_header.samples_in_scan
        nscan = len(self.scans)
        ncha = self.scans[0].scan_header.used_channels

        nparr = np.zeros([ncha, nscan, npts], np.complex)

        for i, scan in enumerate(self.scans):
            for j, chan in enumerate(scan.channels):
                nparr[j, i, :] = np.array(chan[1])

        return nparr

    def get_data_numpy_scan_channel_order(self):
        """
        Sort all scans into numpy array in [fid, channel_id, spectral points] order.

        NB. Right now, this is the order funct_coil_combine module expects

        This does not require unique dimensions. It does require that total
        number of scans divided by number of channels is an integer. It is
        assumed that all channels for a given FID follow one after the other
        in the scan list.

        NB. channel_id is not always 0 .. Nchan, sometimes it is an odd
            assorment of integers, maybe the elements turned on, so we have
            to use an enumeration to fill the icha index

        """
        npts = self.scans[0].scan_header.samples_in_scan
        nscan = len(self.scans)
        ncha = self.scans[0].scan_header.used_channels

        nparr = np.zeros([nscan, ncha, npts], np.complex)

        for i, scan in enumerate(self.scans):
            for j, chan in enumerate(scan.channels):
                nparr[i, j, :] = np.array(chan[1])

        return nparr

    def get_data_numpy_rep_coil_avg_npts_order(self, return_prep=False):
        return self.get_data_numpy_rep_channel_set_order(return_prep=return_prep)

    def get_data_numpy_rep_channel_set_order(self, return_prep=False):
        """
        Use the ICE loop indices to sort each scan into a numpy array with
        all 16 dimensions.  This option can not be used if we don't have a
        unique set of ICE indices. We use check_unique_indices() to determine
        if this is the case.  If everything is OK, we return a numpy array,
        if it is not, we return None.

        NB. channel_id is not always 0 .. Nchan, sometimes it is an odd
            assorment of integers, maybe the elements turned on, so we have
            to use an enumeration to fill the icha index

        """
        items = ['ide', 'idd', 'idc', 'idb', 'ida', 'segment', 'phase', 'echo', 'partition', 'slice', 'acquisition',
                 'line']
        for item in items:
            indx = SCAN_INDICES.index(item)
            if self.dims[indx] > 1:
                msg = "ICE dimension - '" + item + "' is greater than 1, can not create a unique numpy array."
                raise ValueError(msg)

        nrep = self.dims[MULTI_INDICES.index('repetition')]  # number of repetitions, may be 1
        nset = self.dims[MULTI_INDICES.index('set')]  # number of FIDs for spectroscopy
        ncha = self.scans[0].scan_header.used_channels
        npts = self.scans[0].scan_header.samples_in_scan

        nparr = np.zeros([nrep, ncha, nset, npts], np.complex)

        prep_arr = None
        if return_prep:
            header = _parse_protocol_data(self.evps[3][1])
            nprep_hdr = int(header['sSpecPara.lPreparingScans'])
            nprep = len(self.scans) - nrep * nset
            if nprep_hdr != nprep: print('Warning: get_data_numpy_rep_channel_set_order() - nprep_hdr != nprep')

            # I have seen in CMRR sLASER (VE) data that the scan.set value does not
            # increment while the prep scans are taken. We get around this empirically.
            # We know how many Averages are being taken, subtract from the total number
            # of scans and know how many prep averages there will be.  I calc the iset
            # index below using this sort of math. Just in case someday it start to
            # increment, I make allowances for that case too.

            prep_arr = None
            if nprep != 0:
                prep_arr = np.zeros([nrep, ncha, nprep, npts], np.complex)
                for iscan, scan in enumerate(self.scans[0:nprep]):
                    irep = scan.scan_header.repetition
                    if scan.scan_header.set == 0:
                        iset = iscan
                    else:
                        iset = scan.scan_header.set
                    # iset = scan.scan_header.set
                    for icha, chan in enumerate(scan.channels):
                        prep_arr[irep, icha, iset, :] = np.array(chan[1])  # may need iset -> iscan here, need test data
        else:
            nprep = 0

        for scan in self.scans[nprep:]:
            irep = scan.scan_header.repetition
            iset = scan.scan_header.set
            for icha, chan in enumerate(scan.channels):
                nparr[irep, icha, iset, :] = np.array(chan[1])

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

        NB. channel_id is not always 0 .. Nchan, sometimes it is an odd
            assorment of integers, maybe the elements turned on, so we have
            to use an enumeration to fill the icha index

        """
        items = ['ide', 'idd', 'idc', 'idb', 'ida', 'segment', 'phase', 'partition', 'slice', 'acquisition',
                 'line', 'repetition']
        for item in items:
            indx = SCAN_INDICES.index(item)
            if self.dims[indx] > 1:
                msg = "ICE dimension - '" + item + "' is greater than 1, can not create a unique numpy array."
                raise ValueError(msg)

        nrep = self.dims[MULTI_INDICES.index('echo')]  # number of echos, may be 1 (NB. repurposed nrep here)
        nset = self.dims[MULTI_INDICES.index('set')]  # number of FIDs for spectroscopy

        ncha = self.scans[0].scan_header.used_channels
        npts = self.scans[0].scan_header.samples_in_scan

        nparr = np.zeros([nrep, ncha, nset, npts], np.complex)

        prep_arr = None
        if return_prep:
            header = _parse_protocol_data(self.evps[3][1])
            nprep_hdr = int(header['sSpecPara.lPreparingScans'])
            nprep = len(self.scans) - nrep * nset
            if nprep_hdr != nprep: print('Warning: get_data_numpy_rep_channel_set_order() - nprep_hdr != nprep')

            # I have seen in CMRR sLASER (VE) data that the scan.set value does not
            # increment while the prep scans are taken. We get around this empirically.
            # We know how many Averages are being taken, subtract from the total number
            # of scans and know how many prep averages there will be.  I calc the iset
            # index below using this sort of math. Just in case someday it start to
            # increment, I make allowances for that case too.

            prep_arr = None
            if nprep != 0:
                prep_arr = np.zeros([nrep, ncha, nprep, npts], np.complex)
                for iscan, scan in enumerate(self.scans[0:nprep]):
                    irep = scan.scan_header.echo
                    if scan.scan_header.set == 0:
                        iset = iscan
                    else:
                        iset = scan.scan_header.set
                    # iset = scan.scan_header.set
                    for icha, chan in enumerate(scan.channels):
                        prep_arr[irep, icha, iset, :] = np.array(chan[1])  # may need iset -> iscan here, need test data
        else:
            nprep = 0

        for scan in self.scans[nprep:]:
            irep = scan.scan_header.echo
            iset = scan.scan_header.set
            for icha, chan in enumerate(scan.channels):
                nparr[irep, icha, iset, :] = np.array(chan[1])

        return nparr, prep_arr









class TwixMultiRaid(object):
    """
    Before VD11A only a single measurement data can be copied by TWIX
    - Problem: If PreScan-Normalize or NoiseDecorrelation is used also
        these pre-measurement data have to copied to run ICE simulation.
    - Solution: New Multi-RAID file in VD11A (selectable in TWIX)
        collects automatically all dependent measUIDs needed for ICE
        simulation    
    
    TwixSuperRaidFileHeader 
    -------------------- 
    Describes how the measurements are laid out in the overall file. 

    Note. the last measurement contains the "real" image scan data, earlier 
    measurements are adjustments, pre-scans, noise-scans etc. 
    - The last measurement depends on the preceeding measurements
    - Use measurement_offset attribute to determine the position of a 
       measurement in file
    - There are always 64 MrParcRaidFileEntry structs, indepenent of the 
       actual number of measurements
    - After the SuperRaidFileHeader and after every measurement there is a 512 
       byte alignment filled with "0"
     
    The Measurement global header (one per Measurement) format has not changed
     
    The Mdh header format for each scan in a measurement has also been extended 
     (128 bytes to 192 bytes) to provide additional information about what is 
     happening during each scan

    
    """
    def __init__(self):
        
        self.super_raid_header = None
        self.measurements = []
        self.current_measurement_index = -1


    @property
    def current(self):
        return self.measurements[self.current_measurement_index]

    @property
    def evps(self):
        return self.current.evps

    @property
    def scans(self):
        return self.current.scans

    @property
    def free_parameters(self):
        return self.current.free_parameters

    @property
    def ice_parameters(self):
        return self.current.ice_parameters


    ########  General methods  ################################################

    def __str__(self):
        s = self.super_raid_header.__str__(self)
        return  s


    def populate_from_file(self, filename, current_measurement_index=-1):
        """
        Given the name of a Twix file, reads the global header and all the
        scans in the file. 

        """
        file_size = os.path.getsize(filename)
        infile = open(filename, 'rb')
        infile.seek(0)

        self.super_raid_header = TwixSuperRaidFileHeader()
        self.super_raid_header.populate_from_file(infile)
        
        for entry in self.super_raid_header.entries:
             
            meas = TwixMeasurement()
            infile.seek(entry.measurement_offset, os.SEEK_SET)    # from beginning
            meas.populate_from_file(infile, entry.measurement_length)
            self.measurements.append(meas)

        # checks if value given in consistent with measurements list
        self.current_measurement_index = self.get_indx(current_measurement_index)
            
            
    def get_measurement(self, index=-1):
        """ 
        Return a measurement from the multi-raid list, default is the last 
        one, since that contains the 'real' data, ie. all the earlier 
        measurements support the final one.
        
        """
        return self.measurements[index]


    ######## Sort functionality methods - convenience calls ###################
    #
    # Sort functionality for each measurement lies within each TwixMeasurement
    # object. Here we do pass through calls to either the default measurement
    # (via the current_measurement_index) or to the user set index in the call.

    def get_indx(self, index):
        indx = self.current_measurement_index
        if index is not None:
            if index >= 0 and index < len(self.measurements):
                indx = index
        return indx

    def get_ice_index(self, iscan, index=None):
        meas = self.measurements[self.get_indx(index)]
        return meas.get_ice_index(iscan)

    def create_ice_indices(self, index=None):
        meas = self.measurements[self.get_indx(index)]
        return meas.create_ice_indices()

    def get_dims(self, index=None):
        meas = self.measurements[self.get_indx(index)]
        return meas.get_dims()

    def check_unique_indices(self, index=None):
        meas = self.measurements[self.get_indx(index)]
        return meas.check_unique_indices()

    def get_data_numpy_ice_indexed(self, index=None):
        meas = self.measurements[self.get_indx(index)]
        return meas.get_data_numpy_ice_indexed()

    def get_data_numpy_scan_order(self, index=None):
        meas = self.measurements[self.get_indx(index)]
        return meas.get_data_numpy_scan_order()

    def get_data_numpy_channel_scan_order(self, index=None):
        meas = self.measurements[self.get_indx(index)]
        return meas.get_data_numpy_channel_scan_order()

    def get_data_numpy_scan_channel_order(self, index=None):
        meas = self.measurements[self.get_indx(index)]
        return meas.get_data_numpy_scan_channel_order()

    def get_data_numpy_rep_coil_avg_npts_order(self, return_prep=False, index=None):
        meas = self.measurements[self.get_indx(index)]
        return meas.get_data_numpy_rep_coil_avg_npts_order(return_prep=return_prep)

    def get_data_numpy_rep_channel_set_order(self, return_prep=False, index=None):
        meas = self.measurements[self.get_indx(index)]
        return meas.get_data_numpy_rep_channel_set_order(return_prep=return_prep)

    def get_data_numpy_echo_coil_avg_npts_order(self, return_prep=False, index=None):
        meas = self.measurements[self.get_indx(index)]
        return meas.get_data_numpy_echo_coil_avg_npts_order(return_prep=return_prep)

    def get_data_numpy_echo_channel_set_order(self, return_prep=False, index=None):
        meas = self.measurements[self.get_indx(index)]
        return meas.get_data_numpy_echo_channel_set_order(return_prep=return_prep)







##################   Public functions start here   ##################   


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


def _parse_protocol_data(prot):
    str, end = prot.find("### ASCCONV BEGIN"), prot.find("### ASCCONV END ###")
    prot = prot[str + len("### ASCCONV BEGIN ###"):end]
    f = lambda pair: (pair[0].strip(), pair[1].strip())
    lines = prot.split('\n')
    return dict([f(line.split('=')) for line in lines[1:] if line])


##################   "Private" functions start here   ##################   

def _bit_string(value, min_length=1):
    """Given an int value, returns a bitwise string representation that is 
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

def _read_double(source_file, count=1):
    return _read_generic(source_file, 'd', count)

def _read_uint(source_file, count=1):
    return _read_generic(source_file, 'I', count)

def _read_int(source_file, count=1):
    return _read_generic(source_file, 'i', count)

def _read_ulonglong(source_file, count=1):
    return _read_generic(source_file, 'Q', count)

def _read_longlong(source_file, count=1):
    return _read_generic(source_file, 'q', count)

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
    # The original source of this code is:
    # http://scion.duhs.duke.edu/vespa/project/browser/trunk/common/util/io_.py

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




if __name__ == "__main__":

    argv = sys.argv

    if len(argv) < 2:
        # this is a test default filename to make my debug life easier
        argv.append('c:\\bsoher\\\\VE11C_svs-se\\meas_MID00765_FID24095.dat')
    
    if len(argv) < 2:
        print("Please enter the name of the file you want to read.")
        sys.exit(0) 
    else:
    
        twix = TwixMultiRaid()
        twix.populate_from_file(argv[1])

        # at this point you've read in a twix file, getting at the data is another layer of complexity
        # see the (non-working) examples below
        bob = 10  

#        # Non-working Example 1. Not using twix_sort
#        # -------------------------------------------------------
#        meas = twix.measurements[-1]
#        scans, evps = meas.scans, meas.evps
#        ncoils   = scans[0].scan_header.used_channels
#        # alternative for getting number of coils
#        coil_ids = sorted(set([chan[0].channel_id for scan in scans for chan in scan.channels ]))
#        ncoils   = len(coil_ids)
#        nscans   = len(scans)
#        for scan in scans:
#            for item in scan.channels:
#                chan = item[1]
#                chan = np.array(chan[0:1024])     # note. spectral points hard set here, may vary in your data
#                # at this point, do coil combine or store in a numpy array or whatever
#
#
#        # Non-working Example 2. Using twix_sort
#        # -------------------------------------------------------
#        raw = twix.get_data_numpy_channel_scan_order()
#        nscans = int(raw.shape[-2])





