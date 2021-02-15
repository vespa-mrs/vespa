#
# Some of the code in this file was derived from the Python package 
# pfile-tools project, https://github.com/njvack/pfile-tools
# and as such we have included their BSD statement in this file.
#
# Copyright (c) 2012, Board of Regents of the University of Wisconsin
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
# 
# * Redistributions of source code must retain the above copyright notice, this list
#   of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice, this
#   list of conditions and the following disclaimer in the documentation and/or
#   other materials provided with the distribution.
# * Neither the name of the University of Wisconsin nor the names of its 
#   contributors may be used to endorse or promote products derived from this
#   software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Python modules
import os
import math
import sys
import csv
import ctypes

# Third party packages
import numpy as np

# Our Modules
#import vespa.common.ge_util as utilge
#import vespa.common.ge_pfile_mapper as pfile_mapper

from ctypes import *
from collections import namedtuple

StructInfo = namedtuple("StructInfo", ["label", "depth", "value", "field_type", "size", "offset"])



class UnknownPfile(RuntimeError):
    pass


class RevisionNumLittle(LittleEndianStructure):

    @property
    def major(self):
        return int(self.revision)

    _pack_   = 1
    _fields_ = [ ('revision', c_float) ]

class RevisionNumBig(BigEndianStructure):

    @property
    def major(self):
        return int(self.revision)

    _pack_   = 1
    _fields_ = [ ('revision', c_float) ]



class Pfile(object):
    """
    This class was based on the style of code from the pfile-tools 
    package written by Nathan Vack in that we use ctypes to organize 
    the reading of binary structures into a Python readable class
    instance. We have also incorporated code from their struct
    utilities modules as part of our class in order to dump out the
    header information to the stdout or as a list of strings.

    We use a subset of header variables that are sufficient to read the
    data from P-files in which we are interested. The structure for these
    variables was adapted from similar code found in the UCSF Sivic project.
    

    """
    def __init__(self, fname):
    
        self.file_name  = fname
        self.version    = 0
        self.hdr        = None
        self.map        = None
        self.endian     = 'little'  # def for version >= 11
        
        self.read_header()
        
        if not self.is_ge_file:
            raise UnknownPfile("Not a known GE Pfile - fname = %s" % fname)
        
        self.map_data()
        
    # Properties --------------------------------------------------

    @property
    def is_ge_file(self):
        if self.version < 12:
            if "GE" in self.hdr.rhr_rh_logo:
                return True
            else:
                return False
        else:
            offset = self.hdr.rhr_rdb_hdr_off_data
            if ( offset == 61464  or        # bjs from matlap script for ver = 9 
                 offset == 66072  or 
                 offset == 145908 or 
                 offset == 149788 or
                 offset == 150336 or
                 offset == 157276 or        # v24 empirical
                 offset == 213684 ):        # v26 empirical
                return True
            else:
                return False


    @property
    def is_svs(self):
        
        if self.map is None: 
            return False
        else:
            return self.map.is_svs     


    @property
    def get_mapper(self):
        
        if self.hdr is None: 
            return None
            
        psd = self.hdr.rhi_psdname.decode('utf-8').lower()
        
        if psd == 'probe-p':
            mapper = PfileMapper
        elif psd == 'oslaser':
            mapper = PfileMapperSlaser
        elif psd == 'presscsi':
            mapper = PfileMapper
        elif psd == 'fidcsi':
            # bjs - added for Pom's fidcsi 13C data
            mapper = PfileMapper
        elif psd == 'ia/stable/fidcsi':
            # bjs - added for Kearny's 13C data
            mapper = PfileMapper
        elif psd == 'presscsi_nfl':
            # bjs - added for Govind's SVS data off v25
            mapper = PfileMapper
        elif psd == 'epsi_3d_24':
            # bjs - added for soher check of MIDAS Browndyke data
            mapper = PfileMapper
        else:
            raise UnknownPfile("No Pfile mapper for pulse sequence = %s" % psd)

#        if psd == 'probe-p':
#            mapper = pfile_mapper.PfileMapper
#        elif psd == 'oslaser':
#            mapper = pfile_mapper.PfileMapperSlaser
#        elif psd == 'presscsi':
#            mapper = pfile_mapper.PfileMapper
#        elif psd == 'fidcsi':
#            # bjs - added for Pom's fidcsi 13C data
#            mapper = pfile_mapper.PfileMapper
#        elif psd == 'ia/stable/fidcsi':
#            # bjs - added for Kearny's 13C data
#            mapper = pfile_mapper.PfileMapper
#        elif psd == 'presscsi_nfl':
#            # bjs - added for Govind's SVS data off v25
#            mapper = pfile_mapper.PfileMapper
#        elif psd == 'epsi_3d_24':
#            # bjs - added for soher check of MIDAS Browndyke data
#            mapper = pfile_mapper.PfileMapper
#        else:
#            raise UnknownPfile("No Pfile mapper for pulse sequence = %s" % psd)
        
        return mapper


        
    def read_header(self):

        filelike = open(self.file_name, 'rb')

        # determine version number of this header from revision of rdbm.h
        version = self._major_version(filelike)
        if version == 0:
            raise UnknownPfile("Pfile not supported for version %s" % version)    
  
        # Here we dynamically configure the ctypes structures into which the
        # binary file will be read, based on the revision number
        #
        # Note. Determined empirically that I cannot declare the XxxHeader
        # class at the top level of the module with an attribute ._fields_ = [] 
        # and then append into it. I have to create a list and then assign 
        # _fields_ attribute to that list in one step.  Don't know why.
        #
        # Note 2. Had to move Class definition into this function so that the
        # class can be reconstituted more than once for multiple GE file reads.
        # At the top level of the module, the _fields_ attribute could be 
        # created once dynamically, but afterwards would stick around and
        # could not then be changed. 
        
        if version < 11:  # data taken on SGI - big endian
            class PfileHeaderBig(BigEndianStructure):
                """
                Contains the ctypes Structure for a GE P-file rdb header.
                Dynamically allocate the ctypes _fields_ list later depending on revision
                """
                _pack_   = 1            
                _fields_ = get_pfile_hdr_fields(version)
                # _fields_ = utilge.get_pfile_hdr_fields(version)
            hdr = PfileHeaderBig()
            self.endian = 'big'
        else:
            class PfileHeaderLittle(LittleEndianStructure):
                """
                Contains the ctypes Structure for a GE P-file rdb header.
                Dynamically allocate the ctypes _fields_ list later depending on revision
                """
                _pack_   = 1
                _fields_ = get_pfile_hdr_fields(version)
                # _fields_ = utilge.get_pfile_hdr_fields(version)
            hdr = PfileHeaderLittle()
            self.endian = 'little'
  
        try:
            # read  header information from start of file
            filelike.seek(0)
            filelike.readinto(hdr)
            filelike.close()
        except:
            filelike.close()
            raise UnknownPfile("Trouble reading file into header structure for version %s" % version)
        
        self.version = version
        self.hdr = hdr


    def map_data(self):
        """
        Select appropriate mapper class using the pulse sequence name string,
        instantiate and read the data from the file into the 'map' attribute
        
        """
        mapper = self.get_mapper
        self.map = mapper(self.file_name, self.hdr, self.version, self.endian)
        self.map.read_data()
        
        

    def _major_version(self, filelike):
        """
        Get the rdbm.h revision number from first 4 bytes. Then map the rdbm 
        revision to a platform number (e.g. 11.x, 12.x, etc.)
        
        """
        rev_little = RevisionNumLittle()
        rev_big    = RevisionNumBig()
        filelike.seek(0)
        filelike.readinto(rev_little)
        filelike.seek(0)
        filelike.readinto(rev_big)
        
        rev_little = rev_little.major
        rev_big    = rev_big.major
        
        version = 0

        if (rev_little > 6.95 and rev_little < 8.0) or (rev_big > 6.95 and rev_big < 8.0):
            version = 9.0; 
        elif ( rev_little == 9.0  or rev_big == 9.0  ): 
            version = 11.0;
        elif ( rev_little == 11.0 or rev_big == 11.0 ): 
            version = 12.0;
        elif ( rev_little == 14 or rev_big == 14 ): 
            version = 14.0;
        elif ( rev_little == 15 or rev_big == 15 ): 
            version = 15.0;
        elif ( rev_little == 16 or rev_big == 16 ): 
            version = 16.0;
        elif ( rev_little == 20 or rev_big == 20 ): 
            version = 20.0;
        elif ( rev_little == 21 or rev_big == 21 ): 
            version = 21.0;
        elif ( rev_little == 23 or rev_big == 23 ): 
            version = 21.0;
        elif ( rev_little == 24 or rev_big == 24 ): 
            version = 24.0;
        elif ( rev_little == 26 or rev_big == 26 ): 
            version = 26.0;
        else:
            raise UnknownPfile("Unknown header structure for revision %s" % rev_little)

        return version; 


    def dump_header_strarr(self):

        dumped = self._dump_struct(self.hdr)
        strarr = []
        for info in dumped:
            if (info.label.find("pad") == 0):
                continue
            # needed this because Gregor's UID values had some odd chars which
            # caused errors on the import of Probep data 
            try:
                val = info.value
                # needed this because Pom's UID values were causing errors
                # when VIFF was read back in
                if info.label == 'rhe_study_uid':    val = ' '
                if info.label == 'rhs_series_uid':   val = ' '
                if info.label == 'rhs_landmark_uid': val = ' '
                if info.label == 'rhi_image_uid':    val = ' '
                val = str(val)
            except:
                val = ' '
            strarr.append(str(info.label)+'        '+val)
            
        return strarr
        

    def dump_header(self):

        dumped = self._dump_struct(self.hdr)
        writer = csv.writer(sys.stdout, delimiter="\t")
        writer.writerow(["\n\nHeader All  -  field", "value"])
        for info in dumped:
            if (info.label.find("pad") == 0):
                continue
            writer.writerow([info.label, str(info.value)])


    def _dump_struct(self, struct, include_structs=False):
        """
        Recursively travels through a ctypes.Structure and returns a list of
        namedtuples, containing label, depth, value, size, and offset.
        If include_structs is true, output will include lines for individual
        structures and their sizes and offsets -- not just non-structure fields.
        """
        output = []
        self._dump_struct_rec(struct, output, include_structs)
        return output
    
    
    def _dump_struct_rec(self, struct, output, include_structs=False, prefix='', depth=0, base_offset=0):
        """
        Internal recursive method for dumping structures.
        Appends to the "output" parameter.
        
        """
        struct_class = type(struct)
        if include_structs:
            output.append(StructInfo(
                "%s (%s)" % (prefix, struct_class.__name__),
                depth, '', str(struct_class), ctypes.sizeof(struct_class), base_offset))
        for f in struct._fields_:
            name = f[0]
            field_type = f[1]
            field_meta = getattr(struct_class, name)
            field = getattr(struct, name)
            cur_prefix = "%s%s." % (prefix, name)
            field_offset = base_offset + field_meta.offset
            if isinstance(field, ctypes.Structure):
                self._dump_struct_rec(field, output, include_structs, cur_prefix,depth+1, field_offset)
            else:
                label = prefix+name
                output.append(StructInfo(label, depth, field, field_type.__name__, field_meta.size, field_offset))
    


#------------------------------------------------------------------------------
# originally ge_pfile_mapper.py
#------------------------------------------------------------------------------



class PfileMapper(object):

    def __init__(self, file_name, hdr, version, endian):
        """
        Given a file name, its header, version number and endianess, this
        class will parse the data section of the file for the suppressed and
        unsuppressed data.  
        
        All 'timePts' (aka. FID data arrays) are stored in the raw_data 
        attribute. It is a numpy ndarray with shape of:
        
        [cols, rows, slices, numTimePts, numCoils, numSpecPts], np.complex64
        
        For SVS data, cols, rows and slices are all equal to 1. 
        
            - raw_suppressed   is a view onto the  water suppressed fids data
            - raw_unsuppressed is a view onto the water unsuppressed fids data
            - avg_suppressed and avg_unsuppressed are numpy arrays where the
                  relevant raw_ views have been summed along the numTimePts
                  dimension. shape = [cols, rows, slices, numCoils, numSpecPts]
                  
        For non-SVS data, only the raw_data attribute has data in it.
        
        History:
        
        Derived from SIVIC file svkGEPFileMapper.cc which was used to map data
        from PROBE-P and PRESSCSI P-files.  SIVIC has other mapper classes for
        other types of P-file data. I will plan on using this model here, too.
        
        """
        
        self.file_name = file_name
        self.hdr       = hdr
        self.version   = version
        self.endian    = endian
        self.is_svs    = False
        
        self.raw_data         = None
        self.raw_suppressed   = None
        self.avg_suppressed   = None
        self.raw_unsuppressed = None
        self.avg_unsuppressed = None
    
        
        
    @property
    def get_select_box_center(self):
        """
        Center position is taken from user variables.  The Z "slice"
        position used to be taken from the image header "image.loc",
        but with the LX architecture, this held the table position only,
        so if Graphic RX was used to introduce an offset, it wouldn't
        be successfully extracted.
        
        """
        center0 = -1 * self.hdr.rhi_user11
        center1 = -1 * self.hdr.rhi_user12
        center2 =      self.hdr.rhi_user13
        
        return np.array([center0, center1, center2])
        

    @property
    def get_select_box_size(self):        
        
        boxsize = np.array([0.0, 0.0, 0.0])
        dcos = self.get_dcos
        
        if self.version > 9:

            lMax   = 0
            pMax   = 0
            sMax   = 0
            lIndex = 0
            pIndex = 0
            sIndex = 0
            for i in range(3):
                if abs( dcos[i][0] ) > lMax:
                    lIndex = i
                    lMax = abs( dcos[i][0] )
                if abs( dcos[i][1] ) > pMax:
                    pIndex = i
                    pMax = abs( dcos[i][1] ) 
                if abs( dcos[i][2] ) > sMax:
                    sIndex = i
                    sMax = abs( dcos[i][2] )

            boxsize[ lIndex ] = self.hdr.rhi_user8
            boxsize[ pIndex ] = self.hdr.rhi_user9
            boxsize[ sIndex ] = self.hdr.rhi_user10

        else:

            boxsize[0] = self.hdr.rhr_roilenx
            boxsize[1] = self.hdr.rhr_roileny
            boxsize[2] = self.hdr.rhr_roilenz

            if self.is_swap_on:
                ftemp = boxsize[0]
                boxsize[0] = boxsize[1]
                boxsize[1] = ftemp
        
        return boxsize


    @property
    def get_voxel_spacing(self):  
        """
        Get the voxel spacing in 3D. Note that the slice spacing may include
        a skip.
        Swaps the FOV if necessary based on freq_dir setting.

        """
        user19 = self.hdr.rhi_user19
        voxspace = np.array([0.0, 0.0, 0.0])
        
        if (user19 > 0)  and (self.version > 9):
            voxspace[0] = user19
            voxspace[1] = user19
            voxspace[2] = user19
        else:
            fov  = self.get_fov
            nvox = self.get_num_voxels
            voxspace[0] = fov[0]/nvox[0]
            voxspace[1] = fov[1]/nvox[1]
            voxspace[2] = fov[2]/nvox[2]       
    
        return voxspace
        

    @property
    def get_fov(self): 
        fov  = np.array([0.0, 0.0, 0.0])
        nvox = self.get_num_voxels
    
        dfov = self.hdr.rhi_dfov
        
        if self.version > 9:

            fov[0] = dfov
            fov[1] = dfov

            # 2D case vs 3D cases
            if self.is_2d:
                fov[2] = self.hdr.rhi_user10
            else:
                fov[2] = self.hdr.rhi_scanspacing * self.hdr.rhr_zcsi
        else:
            fov[0] = self.hdr.rhr_rh_user7
            fov[1] = self.hdr.rhr_rh_user8
            fov[2] = self.hdr.rhr_rh_user9

        #  Anisotropic voxels:
        if (self.version > 9) and (nvox[0] != nvox[1]):

            # CSI has already been reordered if needed - so fov  calculated 
            # with this CSI will not need reordering, need next power of 2:
            xdim = int(pow(2, math.ceil(math.log(nvox[0], 2))))
            ydim = int(pow(2, math.ceil(math.log(nvox[1], 2))))

            if( ydim > xdim ):
                fov_spatial_resolution = dfov/ydim
            else:
                fov_spatial_resolution = dfov/xdim

            fov[1] = fov_spatial_resolution * ydim
            fov[0] = fov_spatial_resolution * xdim

        elif self.is_swap_on:

            #  Swap the FOV if necessary based on freq dir:
            temp = fov[0]
            fov[0] = fov[1]
            fov[1] = temp

        return fov        
     
     
    @property
    def get_num_voxels(self):   
        """
        Get the 3D spatial dimensionality of the data set
        Returns an int array with 3 dimensions.  Swaps
        if necessary based on freq_dir setting.
        
        """
        nvox = np.array([0, 0, 0])
     
        if self.hdr.rhr_rh_file_contents == 0:
            nvox[0] = 1
            nvox[1] = 1
            nvox[2] = 1
        else:
            nvox[0] = int(self.hdr.rhr_xcsi)
            nvox[1] = int(self.hdr.rhr_ycsi)
            nvox[2] = int(self.hdr.rhr_zcsi)

        #  Swap dimensions if necessary:
        if self.is_swap_on:
            temp = nvox[0]
            nvox[0] = nvox[1]
            nvox[1] = temp

        return nvox
     

    @property
    def get_dcos(self):      
        dcos = np.zeros([3, 3], float)     
     
        dcos[0][0] = -( self.hdr.rhi_trhc_R - self.hdr.rhi_tlhc_R )
        dcos[0][1] = -( self.hdr.rhi_trhc_A - self.hdr.rhi_tlhc_A )
        dcos[0][2] =  ( self.hdr.rhi_trhc_S - self.hdr.rhi_tlhc_S )
        
        dcosLengthX = np.sqrt( dcos[0][0] * dcos[0][0]
                             + dcos[0][1] * dcos[0][1]
                             + dcos[0][2] * dcos[0][2] )

        dcos[0][0] /= dcosLengthX
        dcos[0][1] /= dcosLengthX
        dcos[0][2] /= dcosLengthX


        dcos[1][0] = -( self.hdr.rhi_brhc_R - self.hdr.rhi_trhc_R )
        dcos[1][1] = -( self.hdr.rhi_brhc_A - self.hdr.rhi_trhc_A )
        dcos[1][2] =  ( self.hdr.rhi_brhc_S - self.hdr.rhi_trhc_S )

        dcosLengthY = np.sqrt( dcos[1][0] * dcos[1][0]
                             + dcos[1][1] * dcos[1][1]
                             + dcos[1][2] * dcos[1][2] )

        dcos[1][0] /= dcosLengthY
        dcos[1][1] /= dcosLengthY
        dcos[1][2] /= dcosLengthY


        # third row is the vector product of the first two rows
        # actually, -1* vector product, at least for the axial and axial oblique
        # which is all that we support now
        dcos[2][0] = - dcos[0][1] * dcos[1][2] + dcos[0][2] * dcos[1][1]
        dcos[2][1] = - dcos[0][2] * dcos[1][0] + dcos[0][0] * dcos[1][2]
        dcos[2][2] = - dcos[0][0] * dcos[1][1] + dcos[0][1] * dcos[1][0]

        return dcos
     

    @property
    def is_swap_on(self):      
        """ Is frequency direction swapped? """     
        if self.hdr.rhi_freq_dir != 1:
            return True
        else:
            return False
   

    @property
    def is_2d(self):      
        """ Is this a 2D or 3D data set (spatial dimensions)? """        
        is2D = False
        ndims = self.hdr.rhr_csi_dims

        if ndims == 0:
            if self.hdr.rhr_xcsi >= 0:
                ndims += 1
            if self.hdr.rhr_ycsi >= 0:
                ndims += 1
            if self.hdr.rhr_zcsi >= 0:
                ndims += 1

        if ndims == 2:
            is2D = True

        return is2D
   

    @property
    def is_chop_on(self):      
        """ Is data chopped? """        
        chop  = False
        nex   = self.hdr.rhi_nex
        necho = self.hdr.rhi_numecho

        if ( math.ceil(nex) * necho ) <= 1:
            chop = True

        return chop
        
        
    @property
    def get_frequency_offset(self):      
        """ Returns the spectral frquency offset """        
        if self.version > 9:
            return 0.0
        else:
            return self.hdr.rhr_rh_user13

        
    @property
    def get_center_from_raw_file(self):      
        """ 
        Gets the center of the acquisition grid.  May vary between sequences.
        
        """
        center = np.array([0.0, 0.0, 0.0])
        if self.version < 11:
            center[0] = 0
            center[1] = 0
            center[2] = self.hdr.rhi_user13
        else:
            center[0] = -1 * self.hdr.rhi_user11
            center[1] = -1 * self.hdr.rhi_user12
            center[2] = self.hdr.rhi_user13
       
        return center

   
    @property
    def get_num_coils(self):      
        """ Determine number of coils of data in the PFile. """
        ncoils = 0
        for i in range(4):
            start_rcv = getattr(self.hdr, "rhr_rh_dab["+str(i)+"]_start_rcv")
            stop_rcv  = getattr(self.hdr, "rhr_rh_dab["+str(i)+"]_stop_rcv")

            if ( start_rcv != 0) or (stop_rcv != 0):
                ncoils += ( stop_rcv - start_rcv ) + 1

        #  Otherwise 1
        if ncoils == 0:
            ncoils = 1

        return int(ncoils)   
 
 
    @property
    def get_num_time_points(self):      
        """ 
        Determine number of time points in the PFile.
        Number of time points is determined from the file size,
        number of voxels and number of coils.
        """
        passSize      = float(self.hdr.rhr_rh_raw_pass_size)
        numCoils      = float(self.get_num_coils)
        numVoxels     = float(self.get_num_voxels_in_vol)
        dataWordSize  = float(self.hdr.rhr_rh_point_size)
        numFreqPoints = float(self.hdr.rhr_rh_frame_size)
        kSpacePoints  = float(self.get_num_kspace_points)

        numTimePoints = int( ( passSize ) / ( numCoils * 2 * dataWordSize * numFreqPoints ) - 1 ) / kSpacePoints

        # bjs - added this after Pom's fidcsi 13C data came up with 0 here
        if numTimePoints <= 0:
            numTimePoints = 1

        return int(numTimePoints)

       
    @property
    def get_num_dummy_scans(self):        
        """
        Determine number of dummy scans (FIDs) in the data block.
        This is the difference between the raw pass size and the
        expected size of the data based on numCoils, numTimePts, numKSpacePts
        and numFreqPts.        
        
        """
        passSize         = self.hdr.rhr_rh_raw_pass_size
        numCoils         = self.get_num_coils 
        numTimePoints    = self.get_num_time_points 
        numSampledVoxels = self.get_num_kspace_points
        numFreqPoints    = self.hdr.rhr_rh_frame_size
        dataWordSize     = self.hdr.rhr_rh_point_size

        dataRepresentation = "COMPLEX"  # this was hard set in DcmHeader code
        if ( dataRepresentation == "COMPLEX" ):
            numComponents = 2
        else:
            numComponents = 1

        #  Calc the diff between the size of the data buffer and the number of real data points
        #  then divide by the number of bytes in a single fid to get the number of dummy FIDs
        numDummyScans = passSize - ( numCoils * numTimePoints * numSampledVoxels * numFreqPoints * numComponents * dataWordSize )

        numDummyScans = numDummyScans / ( numFreqPoints * numComponents * dataWordSize)

        return int(numDummyScans)
 
   
    @property
    def get_num_frames(self):    
        """ Number of frames is number of slices * numCoils * numTimePoints """   

        nvox = self.get_num_voxels
        nframes = nvox[2] * self.get_num_coils * self.get_num_time_points

        return int(nframes)
   
   
    @property
    def get_num_voxels_in_vol(self):    

        nvox = self.get_num_voxels
        
        return int(nvox[0] * nvox[1] * nvox[2])   
   
   
    @property
    def get_num_kspace_points(self):    
        """
        Determine the number of sampled k-space points in the data set.
        This may differ from the number of voxels in the rectalinear grid,
        for example if elliptical or another non rectangular acquisition
        sampling strategy was employed.  GE product sequences pad the
        reduced k-space data with zeros so the number of k-space points
        is the same as the number of voxels, but that may not be true for
        custom sequences.        
        
        """
        return int(self.get_num_voxels_in_vol)


    @property
    def was_index_sampled(self):    
        """
        Determines whether a voxel (index) was sampled (or a zero padded
        point is present in the data set), or not, i.e. was it within
        the elliptical sampling volume if reduced k-space elliptical sampling
        was used. Could be extended to support other sparse sampling
        trajectories. Note that for product sequences this always returns true
        since GE  zero-pads reduced k-space data to a full rectilinear grid.
        
        """
        return True   


    @property
    def get_number_unsuppressed_acquisitions(self):    
        """
        For single voxel acquisitions, return the number of
        unsuppressed acquisitions.        
        
        """
        nex = self.hdr.rhi_nex
        return int(16 / nex)

        
    @property
    def get_number_suppressed_acquisitions(self):    
        """
        For single voxel acquisitions, return the number of
        suppressed acquisitions.        
        
        """
        nex   = self.hdr.rhi_nex
        user4 = self.hdr.rhi_user4
        return int( user4 / nex )


    def add_dummy(self, offset, coilNum, timePt):      
        """ 
        Determine whether to add a dummy scan. The assumption is that
        the number of dummy scans should be equal to the number of coils
        or numCoils * numTimePts (e.g. for a spectral editing sequence).
        If true, then the an FID worth of data should be skipped over when
        reading data (e.g. frame_size * numComponents, or numFreqPts * numComponents)        
        """
        numDummyScans    = self.get_num_dummy_scans
        numCoils         = self.get_num_coils
        numTimePoints    = self.get_num_time_points 
        numSampledVoxels = self.get_num_kspace_points
        numFreqPoints    = self.hdr.rhr_rh_frame_size
        numComponents    = 2
        numPointsPerFID  = numFreqPoints * numComponents

        #  subtract the number of dummy words from the current offset to see if another
        #  dummy scan should be skipped or not
        
        if numDummyScans == numCoils:
            numWordsBetweenDummies = numSampledVoxels * numPointsPerFID * numTimePoints
            offset = offset - (coilNum * numPointsPerFID)
            # additional time points have an offset that includes the per-coil dummy
            if timePt > 1:
                offset = offset - numPointsPerFID
        
        elif ( numDummyScans == (numCoils * numTimePoints) ):
            numWordsBetweenDummies = numSampledVoxels * numPointsPerFID
            offset = offset - (coilNum * numPointsPerFID) - ( ( coilNum + timePt ) * numPointsPerFID )
        
        elif numDummyScans == 0:    # bjs - added for fidcsi 13C data from Pom
            return False
        
        else:
            pass
            # "ERROR: Can not determine placement of dummy scans in raw file reader. \n"

        addDummy = False
        if ( ( offset % numWordsBetweenDummies ) == 0 ):
            addDummy = True

        return addDummy
   
   
    def get_xyz_indices(self, dataIndex): 
        """
        If swapping is turned on, the data will need to get mapped correctly
        from the input data buffer read from disk (specData) to the correct
        svkImageData arrays. If swap is true, then the data indices are swapped
        and ky is flipped.
        
        """

        numVoxels = self.get_num_voxels

        z = int( dataIndex/(numVoxels[0] * numVoxels[1]) )

        if self.is_swap_on:

            # If swap is on use numVoxels[1] for x dimension and numVoxels[0] for y dimension
            x = int((dataIndex - (z * numVoxels[0] * numVoxels[1]))/numVoxels[1])

            # In addition to swapping reverse the y direction
            y = numVoxels[1] - int( dataIndex % numVoxels[1] ) - 1

        else:
            x = int( dataIndex % numVoxels[0] )
            y = int((dataIndex - (z * numVoxels[0] * numVoxels[1]))/numVoxels[0])

        return x, y, z





    def get_center_from_origin(self, origin, numVoxels, voxelSpacing, dcos): 
        """
        Calculates the LPS center from the origin(toplc).
        
        """
        center = np.array([0.0, 0.0, 0.0])

        for i in range(3):
            center[i] = origin[i]
            for j in range(3):
                center[i] += dcos[j][i] * voxelSpacing[j] * ( numVoxels[j] / 2.0 - 0.5 )


    def get_origin_from_center(self, center, numVoxels, voxelSpacing, dcos): 
        """
        Calculates the LPS origin (toplc) from the center.
        
        """
        origin = np.array([0.0, 0.0, 0.0])

        for i in range(3):
            origin[i] = center[i]
            for j in range(3):
                origin[i] -= dcos[j][i] * voxelSpacing[j] * ( numVoxels[j] / 2.0 - 0.5 )
        
        


    def read_data(self):            
        """
        This method reads data from the pfile and puts the data into 
        the CellData arrays. If elliptical k-space sampling was used, 
        the data is zero-padded.  Other reduced k-space sampling 
        strategies aren't supported yet.
        
        """

        numCoils        = self.get_num_coils
        numTimePts      = self.get_num_time_points
        numSpecPts      = self.hdr.rhr_rh_frame_size
        numFreqPts      = numSpecPts
        numComponents   =  2
        dataWordSize    = self.hdr.rhr_rh_point_size

        numBytesInVol   = self.get_num_kspace_points * numSpecPts * numComponents * dataWordSize
        numBytesPerCoil = numBytesInVol * numTimePts

        numPtsPerSpectrum = numSpecPts * numComponents

        #  one dummy spectrum per volume/coil:
        numDummyBytes = self.get_num_dummy_scans * numPtsPerSpectrum * dataWordSize
        numDummyBytesPerCoil = int(numDummyBytes/numCoils)

        numBytesPerCoil += numDummyBytesPerCoil

        #  Only read in one coil at a time to reduce memory footprint 
        specData_size = int(np.round((numBytesPerCoil / dataWordSize),0))
        specData = np.zeros([specData_size, ],float) 

        try:
            readOffset = self.hdr.rhr_rdb_hdr_off_data
            filelike = open(self.file_name, 'rb')
            filelike.seek(0)
            filelike.seek(readOffset)
        except:
            pass
            # "ERROR: Exception opening/reading file " << pFileNames->GetValue(0) << " => " << e.what() << endl;

        numVoxels       = self.get_num_voxels
        cols            = numVoxels[0]
        rows            = numVoxels[1]
        slices          = numVoxels[2]
        arraysPerVolume = cols * rows * slices

        #  Preallocate data arrays. The API only permits dynamic assignment at end of CellData, so for
        #  swapped cases where we need to insert data out of order they need to be preallocated.
        
        data = np.zeros([cols, rows, slices, numTimePts, numCoils, numSpecPts], np.complex64)

        #  Blank scan prepended to data blocks.
        dummyOffset = numPtsPerSpectrum 

        #  If Chop On, then reinitialize chopVal:
        chopVal = 1
        if self.is_chop_on:
            chopVal = -1 

        #  pFileOOffset is the offset of the current data set within the 
        #  pFile (i.e. a global data offset).  #  a given coil. 
        pFileOffset = 0

        for coilNum in range(numCoils):
            #  Offset is the offset within the current block of loaded data (ie. within a given coil)
            offset = 0 

            if dataWordSize == 4:
                tempData = np.fromfile(filelike , dtype='i4' , count=int(numBytesPerCoil/dataWordSize) )
            elif dataWordSize == 2:
                tempData = np.fromfile(filelike , dtype='i2' , count=int(numBytesPerCoil/dataWordSize) )
            if self.endian != 'little':
                tempData.byteswap(True)     # swap in-place
            specData = tempData.astype(np.float32)     

            for timePt in range(numTimePts):

                #  Should a dummy scan be skipped over?
                if self.add_dummy( pFileOffset, coilNum, timePt ):
                    offset      += dummyOffset 
                    pFileOffset += dummyOffset

                for arrayNumber in range(arraysPerVolume):
                    x, y, z = self.get_xyz_indices(arrayNumber)

                    #  if k-space sampling was used check if the point was sampled, or if it needs
                    #  to be zero-padded in the grid.
                    #  if zero-padding don't increment the data pointer offset.

                    wasSampled = self.was_index_sampled

                    #  Default, chop is off, so multiply all values by 1
                    #  Only chop sampled data points.
                    if self.is_chop_on and wasSampled:
                        chopVal *= -1

                    if wasSampled:
                        tempFloat = specData[offset:offset+numFreqPts*numComponents]
                        data[x, y, z, timePt, coilNum, :] = chopVal * tempFloat.view(np.complex64)
                    else:
                        pass # do nothing because array is initialized to zeros 

                    if wasSampled:
                        offset      += numPtsPerSpectrum
                        pFileOffset += numPtsPerSpectrum


        filelike.close() 

        self.raw_data = data
        
        #  Modify the data loading behavior.  For single voxel multi-acq data
        #  this means return the averaged (suppresssed data, if applicable).

        numUnsuppressed = self.get_number_unsuppressed_acquisitions 
        numSuppressed   = self.get_number_suppressed_acquisitions

        # bjs - changed numTimePoints to >= 1 below due to Pom's fidcsi 13C data
        # bjs - changed numSuppressed to >= 1 below due to oslaser dv26 data

        # Check if single voxel
        self.is_svs = False
        if ( (numVoxels[0] * numVoxels[1] * numVoxels[2] == 1) and (numTimePts >= 1) and (numSuppressed >= 1)): 
            self.is_svs = True
        
        if self.is_svs:
            self.raw_unsuppressed = self.raw_data[:,:,:,0:numUnsuppressed,:,:]
            self.raw_suppressed   = self.raw_data[:,:,:,numUnsuppressed:,:,:]
            self.avg_unsuppressed = np.sum(self.raw_unsuppressed, axis=3) / float(numUnsuppressed)
            self.avg_suppressed   = np.sum(self.raw_suppressed,   axis=3) / float(numSuppressed)
            self.phase_first_point_deg = np.angle(self.avg_unsuppressed[0,0,0,:,0], True)
        else:
            self.raw_suppressed   = None
            self.avg_suppressed   = None
            self.raw_unsuppressed = None
            self.avg_unsuppressed = None
            self.phase_of_first_point_deg = None
            
        return
    


class PfileMapperSlaser(PfileMapper):

    def __init__(self, file_name, hdr, version, endian):
        """
        This info was provided by Ralph Noeske at GE who wrote the sLASER
        sequence in collaboration with Gulin Oz. Needed to get a few 
        parameters from different header locations, but otherwise the basic
        code seems to work OK.

        rdb_hdr_user0 = rhuser0 = 6024 (sampling frequency)
        rdb_hdr_user4 = rhuser4 = 64 (number of averages)
        rdb_hdr_user19 = rhuser19 = 8 (number of reference frames)
        
        image_user8 = opuser8 is used as starting for voxel dimension and location (opuser8,9,10,11,12,13)
        
        The sub-echo-times are TE1 = 8ms and TE2 = 12ms (opuser20 and 21).  

        
        """
        PfileMapper.__init__(self, file_name, hdr, version, endian)


    @property
    def get_number_unsuppressed_acquisitions(self):    
        """
        For single voxel acquisitions, return the number of
        unsuppressed acquisitions.        
        
        """
        user19 = self.hdr.rhr_rh_user19
        return int(user19)

        
    @property
    def get_number_suppressed_acquisitions(self):    
        """
        For single voxel acquisitions, return the number of
        suppressed acquisitions.        
        
        """
        user4 = self.hdr.rhr_rh_user4
        return int( user4 )
        
#------------------------------------------------------------------------------
# originally ge_util.py
#------------------------------------------------------------------------------



def get_pfile_hdr_fields(version):
    """
    This function returs a list of fields for mapping a GE P-file header
    to a ctypes class structure. We define different paddings and variable
    names depending on the software version that created the P-file.
    
    """
    plist = []

    if version == 9:
        plist.append( ('rhr_rh_rdbm_rev',          c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhr_rh_scan_date',         c_char * 10) )
        plist.append( ('rhr_rh_scan_time',         c_char * 8) )
        plist.append( ('rhr_rh_logo',              c_char * 10) )
        plist.append( ('rhr_rh_file_contents',     c_short) )
        plist.append( ('pad_xx',                   c_char * 10) )
        plist.append( ('rhr_rh_data_collect_type', c_short) )
        plist.append( ('pad_xx',                   c_char * 6) )
        plist.append( ('rhr_rh_npasses',           c_short) )
        plist.append( ('pad_xx',                   c_char * 2) )
        plist.append( ('rhr_rh_nslices',           c_short) )
        plist.append( ('pad_xx',                   c_char * 10) )
        plist.append( ('rhr_rh_frame_size',        c_short) )
        plist.append( ('rhr_rh_point_size',        c_short) )
        plist.append( ('pad_xx',                   c_char * 32) )
        plist.append( ('rhr_rh_raw_pass_size',     c_int) )
        plist.append( ('pad_xx',                   c_char * 80) )
        plist.append( ('rhr_rh_dab[0]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[0]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[1]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[1]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[2]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[2]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[3]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[3]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_user0',             c_float) )
        plist.append( ('rhr_rh_user1',             c_float) )
        plist.append( ('rhr_rh_user2',             c_float) )
        plist.append( ('rhr_rh_user3',             c_float) )
        plist.append( ('rhr_rh_user4',             c_float) )
        plist.append( ('rhr_rh_user5',             c_float) )
        plist.append( ('rhr_rh_user6',             c_float) )
        plist.append( ('rhr_rh_user7',             c_float) )
        plist.append( ('rhr_rh_user8',             c_float) )
        plist.append( ('rhr_rh_user9',             c_float) )
        plist.append( ('rhr_rh_user10',            c_float) )
        plist.append( ('rhr_rh_user11',            c_float) )
        plist.append( ('rhr_rh_user12',            c_float) )
        plist.append( ('rhr_rh_user13',            c_float) )
        plist.append( ('rhr_rh_user14',            c_float) )
        plist.append( ('rhr_rh_user15',            c_float) )
        plist.append( ('rhr_rh_user16',            c_float) )
        plist.append( ('rhr_rh_user17',            c_float) )
        plist.append( ('rhr_rh_user18',            c_float) )
        plist.append( ('rhr_rh_user19',            c_float) )
        plist.append( ('pad_xx',                   c_char * 72) )
        plist.append( ('rhr_spectral_width',       c_float) )
        plist.append( ('rhr_csi_dims',             c_short) )
        plist.append( ('rhr_xcsi',                 c_short) )
        plist.append( ('rhr_ycsi',                 c_short) )
        plist.append( ('rhr_zcsi',                 c_short) )
        plist.append( ('rhr_roilenx',              c_float) )
        plist.append( ('rhr_roileny',              c_float) )
        plist.append( ('rhr_roilenz',              c_float) )
        plist.append( ('pad_xx',                   c_char * 32) )
        plist.append( ('rhr_rh_ps_mps_freq',       c_int) )
        plist.append( ('pad_xx',                   c_char * 560) )
        plist.append( ('rhr_rh_user_usage_tag',    c_uint) )
        plist.append( ('pad_xx',                   c_char * 8) )
        plist.append( ('rhr_rh_user20',            c_float) )
        plist.append( ('rhr_rh_user21',            c_float) )
        plist.append( ('rhr_rh_user22',            c_float) )
        plist.append( ('rhr_rh_user23',            c_float) )
        plist.append( ('rhr_rh_user24',            c_float) )
        plist.append( ('rhr_rh_user25',            c_float) )
        plist.append( ('rhr_rh_user26',            c_float) )
        plist.append( ('rhr_rh_user27',            c_float) )
        plist.append( ('rhr_rh_user28',            c_float) )
        plist.append( ('rhr_rh_user29',            c_float) )
        plist.append( ('rhr_rh_user30',            c_float) )
        plist.append( ('rhr_rh_user31',            c_float) )
        plist.append( ('rhr_rh_user32',            c_float) )
        plist.append( ('rhr_rh_user33',            c_float) )
        plist.append( ('rhr_rh_user34',            c_float) )
        plist.append( ('rhr_rh_user35',            c_float) )
        plist.append( ('rhr_rh_user36',            c_float) )
        plist.append( ('rhr_rh_user37',            c_float) )
        plist.append( ('rhr_rh_user38',            c_float) )
        plist.append( ('rhr_rh_user39',            c_float) )
        plist.append( ('rhr_rh_user40',            c_float) )
        plist.append( ('rhr_rh_user41',            c_float) )
        plist.append( ('rhr_rh_user42',            c_float) )
        plist.append( ('rhr_rh_user43',            c_float) )
        plist.append( ('rhr_rh_user44',            c_float) )
        plist.append( ('rhr_rh_user45',            c_float) )
        plist.append( ('rhr_rh_user46',            c_float) )
        plist.append( ('rhr_rh_user47',            c_float) )
        plist.append( ('rhr_rh_user48',            c_float) )
        plist.append( ('pad_xx',                   c_char * 35764) )
        plist.append( ('rhe_ex_no',                c_ushort) )
        plist.append( ('rhe_hospname',             c_char * 33) )
        plist.append( ('pad_xx',                   c_char * 41) )
        plist.append( ('rhe_magstrength',          c_int) )
        plist.append( ('rhe_patid',                c_char * 13) )
        plist.append( ('rhe_patname',              c_char * 25) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhe_patsex',               c_short) )
        plist.append( ('pad_xx',                   c_char * 67) )
        plist.append( ('rhe_reqnum',               c_char * 13) )
        plist.append( ('rhe_ex_datetime',          c_int) )
        plist.append( ('rhe_refphy',               c_char * 33) )
        plist.append( ('pad_xx',                   c_char * 79) )
        plist.append( ('rhe_ex_sysid',             c_char * 9) )
        plist.append( ('pad_xx',                   c_char * 27) )
        plist.append( ('rhe_ex_verscre',           c_char * 2) )
        plist.append( ('pad_xx',                   c_char * 84) )
        plist.append( ('rhe_uniq_sys_id',          c_char * 16) )
        plist.append( ('pad_xx',                   c_char * 20) )
        plist.append( ('rhe_study_uid',            c_char * 32) )
        plist.append( ('pad_xx',                   c_char * 66) )
        plist.append( ('rhe_patnameff',            c_char * 65) )
        plist.append( ('rhe_patidff',              c_char * 65) )
        plist.append( ('rhe_reqnumff',             c_char * 17) )
        plist.append( ('rhe_dateofbirth',          c_char * 9) )
        plist.append( ('pad_xx',                   c_char * 310) )
        plist.append( ('rhs_se_no',                c_short) )
        plist.append( ('pad_xx',                   c_char * 8) )
        plist.append( ('rhs_se_desc',              c_char * 30) )
        plist.append( ('pad_xx',                   c_char * 26) )
        plist.append( ('rhs_position',             c_int) )
        plist.append( ('rhs_entry',                c_int) )
        plist.append( ('rhs_anref',                c_char * 3) )
        plist.append( ('pad_xx',                   c_char * 257) )
        plist.append( ('rhs_series_uid',           c_char * 32) )
        plist.append( ('rhs_landmark_uid',         c_char * 32) )
        plist.append( ('pad_xx',                   c_char * 656) )
        plist.append( ('rhi_dfov',                 c_float) )
        plist.append( ('pad_xx',                   c_char * 80) )
        plist.append( ('rhi_scanspacing',          c_float) )
        plist.append( ('pad_xx',                   c_char * 8) )
        plist.append( ('rhi_loc',                  c_float) )
        plist.append( ('rhi_ctr_R',                c_float) )
        plist.append( ('rhi_ctr_A',                c_float) )
        plist.append( ('rhi_ctr_S',                c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhi_tlhc_R',               c_float) )
        plist.append( ('rhi_tlhc_A',               c_float) )
        plist.append( ('rhi_tlhc_S',               c_float) )
        plist.append( ('rhi_trhc_R',               c_float) )
        plist.append( ('rhi_trhc_A',               c_float) )
        plist.append( ('rhi_trhc_S',               c_float) )
        plist.append( ('rhi_brhc_R',               c_float) )
        plist.append( ('rhi_brhc_A',               c_float) )
        plist.append( ('rhi_brhc_S',               c_float) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhi_tr',                   c_int) )
        plist.append( ('rhi_ti',                   c_int) )
        plist.append( ('rhi_te',                   c_int) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhi_numecho',              c_short) )
        plist.append( ('pad_xx',                   c_char * 6) )
        plist.append( ('rhi_nex',                  c_float) )
        plist.append( ('pad_xx',                   c_char * 32) )
        plist.append( ('rhi_mr_flip',              c_short) )
        plist.append( ('pad_xx',                   c_char * 58) )
        plist.append( ('rhi_psdname',              c_char * 33) )
        plist.append( ('pad_xx',                   c_char * 21) )
        plist.append( ('rhi_ctyp',                 c_short) )
        plist.append( ('rhi_cname',                c_char * 17) )
        plist.append( ('pad_xx',                   c_char * 31) )
        plist.append( ('rhi_user0',                c_float) )
        plist.append( ('rhi_user1',                c_float) )
        plist.append( ('rhi_user2',                c_float) )
        plist.append( ('rhi_user3',                c_float) )
        plist.append( ('rhi_user4',                c_float) )
        plist.append( ('rhi_user5',                c_float) )
        plist.append( ('rhi_user6',                c_float) )
        plist.append( ('rhi_user7',                c_float) )
        plist.append( ('rhi_user8',                c_float) )
        plist.append( ('rhi_user9',                c_float) )
        plist.append( ('rhi_user10',               c_float) )
        plist.append( ('rhi_user11',               c_float) )
        plist.append( ('rhi_user12',               c_float) )
        plist.append( ('rhi_user13',               c_float) )
        plist.append( ('rhi_user14',               c_float) )
        plist.append( ('rhi_user15',               c_float) )
        plist.append( ('rhi_user16',               c_float) )
        plist.append( ('rhi_user17',               c_float) )
        plist.append( ('rhi_user18',               c_float) )
        plist.append( ('rhi_user19',               c_float) )
        plist.append( ('rhi_user20',               c_float) )
        plist.append( ('rhi_user21',               c_float) )
        plist.append( ('rhi_user22',               c_float) )
        plist.append( ('rhi_user23',               c_float) )
        plist.append( ('rhi_user24',               c_float) )
        plist.append( ('pad_xx',                   c_char * 240) )
        plist.append( ('rhi_freq_dir',             c_short) )
        plist.append( ('pad_xx',                   c_char * 2) )
        plist.append( ('rhi_image_uid',            c_char * 32) )
        plist.append( ('pad_xx',                   c_char * 100) )
        plist.append( ('rhi_user25',               c_float) )
        plist.append( ('rhi_user26',               c_float) )
        plist.append( ('rhi_user27',               c_float) )
        plist.append( ('rhi_user28',               c_float) )
        plist.append( ('rhi_user29',               c_float) )
        plist.append( ('rhi_user30',               c_float) )
        plist.append( ('rhi_user31',               c_float) )
        plist.append( ('rhi_user32',               c_float) )
        plist.append( ('rhi_user33',               c_float) )
        plist.append( ('rhi_user34',               c_float) )
        plist.append( ('rhi_user35',               c_float) )
        plist.append( ('rhi_user36',               c_float) )
        plist.append( ('rhi_user37',               c_float) )
        plist.append( ('rhi_user38',               c_float) )
        plist.append( ('rhi_user39',               c_float) )
        plist.append( ('rhi_user40',               c_float) )
        plist.append( ('rhi_user41',               c_float) )
        plist.append( ('rhi_user42',               c_float) )
        plist.append( ('rhi_user43',               c_float) )
        plist.append( ('rhi_user44',               c_float) )
        plist.append( ('rhi_user45',               c_float) )
        plist.append( ('rhi_user46',               c_float) )
        plist.append( ('rhi_user47',               c_float) )
        plist.append( ('rhi_user48',               c_float) )

    elif version == 11:
        plist.append( ('rhr_rh_rdbm_rev',          c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhr_rh_scan_date',         c_char * 10) )
        plist.append( ('rhr_rh_scan_time',         c_char * 8) )
        plist.append( ('rhr_rh_logo',              c_char * 10) )
        plist.append( ('rhr_rh_file_contents',     c_short) )
        plist.append( ('pad_xx',                   c_char * 10) )
        plist.append( ('rhr_rh_data_collect_type', c_short) )
        plist.append( ('pad_xx',                   c_char * 6) )
        plist.append( ('rhr_rh_npasses',           c_short) )
        plist.append( ('pad_xx',                   c_char * 2) )
        plist.append( ('rhr_rh_nslices',           c_short) )
        plist.append( ('pad_xx',                   c_char * 10) )
        plist.append( ('rhr_rh_frame_size',        c_ushort) )
        plist.append( ('rhr_rh_point_size',        c_short) )
        plist.append( ('pad_xx',                   c_char * 32) )
        plist.append( ('rhr_rh_raw_pass_size',     c_int) )
        plist.append( ('pad_xx',                   c_char * 80) )
        plist.append( ('rhr_rh_dab[0]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[0]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[1]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[1]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[2]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[2]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[3]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[3]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_user0',             c_float) )
        plist.append( ('rhr_rh_user1',             c_float) )
        plist.append( ('rhr_rh_user2',             c_float) )
        plist.append( ('rhr_rh_user3',             c_float) )
        plist.append( ('rhr_rh_user4',             c_float) )
        plist.append( ('rhr_rh_user5',             c_float) )
        plist.append( ('rhr_rh_user6',             c_float) )
        plist.append( ('rhr_rh_user7',             c_float) )
        plist.append( ('rhr_rh_user8',             c_float) )
        plist.append( ('rhr_rh_user9',             c_float) )
        plist.append( ('rhr_rh_user10',            c_float) )
        plist.append( ('rhr_rh_user11',            c_float) )
        plist.append( ('rhr_rh_user12',            c_float) )
        plist.append( ('rhr_rh_user13',            c_float) )
        plist.append( ('rhr_rh_user14',            c_float) )
        plist.append( ('rhr_rh_user15',            c_float) )
        plist.append( ('rhr_rh_user16',            c_float) )
        plist.append( ('rhr_rh_user17',            c_float) )
        plist.append( ('rhr_rh_user18',            c_float) )
        plist.append( ('rhr_rh_user19',            c_float) )
        plist.append( ('pad_xx',                   c_char * 72) )
        plist.append( ('rhr_spectral_width',       c_float) )
        plist.append( ('rhr_csi_dims',             c_short) )
        plist.append( ('rhr_xcsi',                 c_short) )
        plist.append( ('rhr_ycsi',                 c_short) )
        plist.append( ('rhr_zcsi',                 c_short) )
        plist.append( ('rhr_roilenx',              c_float) )
        plist.append( ('rhr_roileny',              c_float) )
        plist.append( ('rhr_roilenz',              c_float) )
        plist.append( ('pad_xx',                   c_char * 32) )
        plist.append( ('rhr_rh_ps_mps_freq',       c_int) )
        plist.append( ('pad_xx',                   c_char * 560) )
        plist.append( ('rhr_rh_user_usage_tag',    c_uint) )
        plist.append( ('pad_xx',                   c_char * 8) )
        plist.append( ('rhr_rh_user20',            c_float) )
        plist.append( ('rhr_rh_user21',            c_float) )
        plist.append( ('rhr_rh_user22',            c_float) )
        plist.append( ('rhr_rh_user23',            c_float) )
        plist.append( ('rhr_rh_user24',            c_float) )
        plist.append( ('rhr_rh_user25',            c_float) )
        plist.append( ('rhr_rh_user26',            c_float) )
        plist.append( ('rhr_rh_user27',            c_float) )
        plist.append( ('rhr_rh_user28',            c_float) )
        plist.append( ('rhr_rh_user29',            c_float) )
        plist.append( ('rhr_rh_user30',            c_float) )
        plist.append( ('rhr_rh_user31',            c_float) )
        plist.append( ('rhr_rh_user32',            c_float) )
        plist.append( ('rhr_rh_user33',            c_float) )
        plist.append( ('rhr_rh_user34',            c_float) )
        plist.append( ('rhr_rh_user35',            c_float) )
        plist.append( ('rhr_rh_user36',            c_float) )
        plist.append( ('rhr_rh_user37',            c_float) )
        plist.append( ('rhr_rh_user38',            c_float) )
        plist.append( ('rhr_rh_user39',            c_float) )
        plist.append( ('rhr_rh_user40',            c_float) )
        plist.append( ('rhr_rh_user41',            c_float) )
        plist.append( ('rhr_rh_user42',            c_float) )
        plist.append( ('rhr_rh_user43',            c_float) )
        plist.append( ('rhr_rh_user44',            c_float) )
        plist.append( ('rhr_rh_user45',            c_float) )
        plist.append( ('rhr_rh_user46',            c_float) )
        plist.append( ('rhr_rh_user47',            c_float) )
        plist.append( ('rhr_rh_user48',            c_float) )
        plist.append( ('pad_xx',                   c_char * 56244) )
        plist.append( ('rhe_ex_no',                c_ushort) )
        plist.append( ('rhe_hospname',             c_char * 33) )
        plist.append( ('pad_xx',                   c_char * 41) )
        plist.append( ('rhe_magstrength',          c_int) )
        plist.append( ('rhe_patid',                c_char * 13) )
        plist.append( ('rhe_patname',              c_char * 25) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhe_patsex',               c_short) )
        plist.append( ('pad_xx',                   c_char * 67) )
        plist.append( ('rhe_reqnum',               c_char * 13) )
        plist.append( ('rhe_ex_datetime',          c_int) )
        plist.append( ('rhe_refphy',               c_char * 33) )
        plist.append( ('pad_xx',                   c_char * 79) )
        plist.append( ('rhe_ex_sysid',             c_char * 9) )
        plist.append( ('pad_xx',                   c_char * 27) )
        plist.append( ('rhe_ex_verscre',           c_char * 2) )
        plist.append( ('pad_xx',                   c_char * 84) )
        plist.append( ('rhe_uniq_sys_id',          c_char * 16) )
        plist.append( ('pad_xx',                   c_char * 20) )
        plist.append( ('rhe_study_uid',            c_char * 32) )
        plist.append( ('pad_xx',                   c_char * 66) )
        plist.append( ('rhe_patnameff',            c_char * 65) )
        plist.append( ('rhe_patidff',              c_char * 65) )
        plist.append( ('rhe_reqnumff',             c_char * 17) )
        plist.append( ('rhe_dateofbirth',          c_char * 9) )
        plist.append( ('pad_xx',                   c_char * 310) )
        plist.append( ('rhs_se_no',                c_short) )
        plist.append( ('pad_xx',                   c_char * 8) )
        plist.append( ('rhs_se_desc',              c_char * 30) )
        plist.append( ('pad_xx',                   c_char * 26) )
        plist.append( ('rhs_position',             c_int) )
        plist.append( ('rhs_entry',                c_int) )
        plist.append( ('rhs_anref',                c_char * 3) )
        plist.append( ('pad_xx',                   c_char * 257) )
        plist.append( ('rhs_series_uid',           c_char * 32) )
        plist.append( ('rhs_landmark_uid',         c_char * 32) )
        plist.append( ('pad_xx',                   c_char * 1164) )
        plist.append( ('rhi_dfov',                 c_float) )
        plist.append( ('pad_xx',                   c_char * 80) )
        plist.append( ('rhi_scanspacing',          c_float) )
        plist.append( ('pad_xx',                   c_char * 8) )
        plist.append( ('rhi_loc',                  c_float) )
        plist.append( ('rhi_ctr_R',                c_float) )
        plist.append( ('rhi_ctr_A',                c_float) )
        plist.append( ('rhi_ctr_S',                c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhi_tlhc_R',               c_float) )
        plist.append( ('rhi_tlhc_A',               c_float) )
        plist.append( ('rhi_tlhc_S',               c_float) )
        plist.append( ('rhi_trhc_R',               c_float) )
        plist.append( ('rhi_trhc_A',               c_float) )
        plist.append( ('rhi_trhc_S',               c_float) )
        plist.append( ('rhi_brhc_R',               c_float) )
        plist.append( ('rhi_brhc_A',               c_float) )
        plist.append( ('rhi_brhc_S',               c_float) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhi_tr',                   c_int) )
        plist.append( ('rhi_ti',                   c_int) )
        plist.append( ('rhi_te',                   c_int) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhi_numecho',              c_short) )
        plist.append( ('pad_xx',                   c_char * 6) )
        plist.append( ('rhi_nex',                  c_float) )
        plist.append( ('pad_xx',                   c_char * 32) )
        plist.append( ('rhi_mr_flip',              c_short) )
        plist.append( ('pad_xx',                   c_char * 58) )
        plist.append( ('rhi_psdname',              c_char * 33) )
        plist.append( ('pad_xx',                   c_char * 21) )
        plist.append( ('rhi_ctyp',                 c_short) )
        plist.append( ('rhi_cname',                c_char * 17) )
        plist.append( ('pad_xx',                   c_char * 31) )
        plist.append( ('rhi_user0',                c_float) )
        plist.append( ('rhi_user1',                c_float) )
        plist.append( ('rhi_user2',                c_float) )
        plist.append( ('rhi_user3',                c_float) )
        plist.append( ('rhi_user4',                c_float) )
        plist.append( ('rhi_user5',                c_float) )
        plist.append( ('rhi_user6',                c_float) )
        plist.append( ('rhi_user7',                c_float) )
        plist.append( ('rhi_user8',                c_float) )
        plist.append( ('rhi_user9',                c_float) )
        plist.append( ('rhi_user10',               c_float) )
        plist.append( ('rhi_user11',               c_float) )
        plist.append( ('rhi_user12',               c_float) )
        plist.append( ('rhi_user13',               c_float) )
        plist.append( ('rhi_user14',               c_float) )
        plist.append( ('rhi_user15',               c_float) )
        plist.append( ('rhi_user16',               c_float) )
        plist.append( ('rhi_user17',               c_float) )
        plist.append( ('rhi_user18',               c_float) )
        plist.append( ('rhi_user19',               c_float) )
        plist.append( ('rhi_user20',               c_float) )
        plist.append( ('rhi_user21',               c_float) )
        plist.append( ('rhi_user22',               c_float) )
        plist.append( ('pad_xx',                   c_char * 92) )
        plist.append( ('rhi_user23',               c_float) )
        plist.append( ('rhi_user24',               c_float) )
        plist.append( ('pad_xx',                   c_char * 148) )
        plist.append( ('rhi_freq_dir',             c_short) )
        plist.append( ('pad_xx',                   c_char * 2) )
        plist.append( ('rhi_image_uid',            c_char * 32) )
        plist.append( ('pad_xx',                   c_char * 100) )
        plist.append( ('rhi_user25',               c_float) )
        plist.append( ('rhi_user26',               c_float) )
        plist.append( ('rhi_user27',               c_float) )
        plist.append( ('rhi_user28',               c_float) )
        plist.append( ('rhi_user29',               c_float) )
        plist.append( ('rhi_user30',               c_float) )
        plist.append( ('rhi_user31',               c_float) )
        plist.append( ('rhi_user32',               c_float) )
        plist.append( ('rhi_user33',               c_float) )
        plist.append( ('rhi_user34',               c_float) )
        plist.append( ('rhi_user35',               c_float) )
        plist.append( ('rhi_user36',               c_float) )
        plist.append( ('rhi_user37',               c_float) )
        plist.append( ('rhi_user38',               c_float) )
        plist.append( ('rhi_user39',               c_float) )
        plist.append( ('rhi_user40',               c_float) )
        plist.append( ('rhi_user41',               c_float) )
        plist.append( ('rhi_user42',               c_float) )
        plist.append( ('rhi_user43',               c_float) )
        plist.append( ('rhi_user44',               c_float) )
        plist.append( ('rhi_user45',               c_float) )
        plist.append( ('rhi_user46',               c_float) )
        plist.append( ('rhi_user47',               c_float) )
        plist.append( ('rhi_user48',               c_float) )

    elif version == 12:
        plist.append( ('rhr_rh_rdbm_rev',          c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhr_rh_scan_date',         c_char * 10) )
        plist.append( ('rhr_rh_scan_time',         c_char * 8) )
        plist.append( ('rhr_rh_logo',              c_char * 10) )
        plist.append( ('rhr_rh_file_contents',     c_short) )
        plist.append( ('pad_xx',                   c_char * 10) )
        plist.append( ('rhr_rh_data_collect_type', c_short) )
        plist.append( ('pad_xx',                   c_char * 6) )
        plist.append( ('rhr_rh_npasses',           c_short) )
        plist.append( ('pad_xx',                   c_char * 2) )
        plist.append( ('rhr_rh_nslices',           c_short) )
        plist.append( ('pad_xx',                   c_char * 10) )
        plist.append( ('rhr_rh_frame_size',        c_ushort) )
        plist.append( ('rhr_rh_point_size',        c_short) )
        plist.append( ('pad_xx',                   c_char * 32) )
        plist.append( ('rhr_rh_raw_pass_size',     c_int) )
        plist.append( ('pad_xx',                   c_char * 80) )
        plist.append( ('rhr_rh_dab[0]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[0]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[1]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[1]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[2]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[2]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[3]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[3]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_user0',             c_float) )
        plist.append( ('rhr_rh_user1',             c_float) )
        plist.append( ('rhr_rh_user2',             c_float) )
        plist.append( ('rhr_rh_user3',             c_float) )
        plist.append( ('rhr_rh_user4',             c_float) )
        plist.append( ('rhr_rh_user5',             c_float) )
        plist.append( ('rhr_rh_user6',             c_float) )
        plist.append( ('rhr_rh_user7',             c_float) )
        plist.append( ('rhr_rh_user8',             c_float) )
        plist.append( ('rhr_rh_user9',             c_float) )
        plist.append( ('rhr_rh_user10',            c_float) )
        plist.append( ('rhr_rh_user11',            c_float) )
        plist.append( ('rhr_rh_user12',            c_float) )
        plist.append( ('rhr_rh_user13',            c_float) )
        plist.append( ('rhr_rh_user14',            c_float) )
        plist.append( ('rhr_rh_user15',            c_float) )
        plist.append( ('rhr_rh_user16',            c_float) )
        plist.append( ('rhr_rh_user17',            c_float) )
        plist.append( ('rhr_rh_user18',            c_float) )
        plist.append( ('rhr_rh_user19',            c_float) )
        plist.append( ('pad_xx',                   c_char * 72) )
        plist.append( ('rhr_spectral_width',       c_float) )
        plist.append( ('rhr_csi_dims',             c_short) )
        plist.append( ('rhr_xcsi',                 c_short) )
        plist.append( ('rhr_ycsi',                 c_short) )
        plist.append( ('rhr_zcsi',                 c_short) )
        plist.append( ('rhr_roilenx',              c_float) )
        plist.append( ('rhr_roileny',              c_float) )
        plist.append( ('rhr_roilenz',              c_float) )
        plist.append( ('pad_xx',                   c_char * 32) )
        plist.append( ('rhr_rh_ps_mps_freq',       c_int) )
        plist.append( ('pad_xx',                   c_char * 560) )
        plist.append( ('rhr_rh_user_usage_tag',    c_uint) )
        plist.append( ('pad_xx',                   c_char * 8) )
        plist.append( ('rhr_rh_user20',            c_float) )
        plist.append( ('rhr_rh_user21',            c_float) )
        plist.append( ('rhr_rh_user22',            c_float) )
        plist.append( ('rhr_rh_user23',            c_float) )
        plist.append( ('rhr_rh_user24',            c_float) )
        plist.append( ('rhr_rh_user25',            c_float) )
        plist.append( ('rhr_rh_user26',            c_float) )
        plist.append( ('rhr_rh_user27',            c_float) )
        plist.append( ('rhr_rh_user28',            c_float) )
        plist.append( ('rhr_rh_user29',            c_float) )
        plist.append( ('rhr_rh_user30',            c_float) )
        plist.append( ('rhr_rh_user31',            c_float) )
        plist.append( ('rhr_rh_user32',            c_float) )
        plist.append( ('rhr_rh_user33',            c_float) )
        plist.append( ('rhr_rh_user34',            c_float) )
        plist.append( ('rhr_rh_user35',            c_float) )
        plist.append( ('rhr_rh_user36',            c_float) )
        plist.append( ('rhr_rh_user37',            c_float) )
        plist.append( ('rhr_rh_user38',            c_float) )
        plist.append( ('rhr_rh_user39',            c_float) )
        plist.append( ('rhr_rh_user40',            c_float) )
        plist.append( ('rhr_rh_user41',            c_float) )
        plist.append( ('rhr_rh_user42',            c_float) )
        plist.append( ('rhr_rh_user43',            c_float) )
        plist.append( ('rhr_rh_user44',            c_float) )
        plist.append( ('rhr_rh_user45',            c_float) )
        plist.append( ('rhr_rh_user46',            c_float) )
        plist.append( ('rhr_rh_user47',            c_float) )
        plist.append( ('rhr_rh_user48',            c_float) )
        plist.append( ('pad_xx',                   c_char * 352) )
        plist.append( ('rhr_rdb_hdr_off_data',     c_int) )
        plist.append( ('pad_xx',                   c_char * 60088) )
        plist.append( ('rhe_magstrength',          c_int) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhe_ex_datetime',          c_int) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhe_ex_no',                c_ushort) )
        plist.append( ('pad_xx',                   c_char * 22) )
        plist.append( ('rhe_patsex',               c_short) )
        plist.append( ('pad_xx',                   c_char * 75) )
        plist.append( ('rhe_reqnum',               c_char * 13) )
        plist.append( ('rhe_refphy',               c_char * 33) )
        plist.append( ('pad_xx',                   c_char * 105) )
        plist.append( ('rhe_ex_sysid',             c_char * 9) )
        plist.append( ('pad_xx',                   c_char * 14) )
        plist.append( ('rhe_hospname',             c_char * 33) )
        plist.append( ('rhe_patid',                c_char * 13) )
        plist.append( ('rhe_patname',              c_char * 25) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhe_ex_verscre',           c_char * 2) )
        plist.append( ('pad_xx',                   c_char * 2) )
        plist.append( ('rhe_uniq_sys_id',          c_char * 16) )
        plist.append( ('pad_xx',                   c_char * 20) )
        plist.append( ('rhe_study_uid',            c_char * 32) )
        plist.append( ('pad_xx',                   c_char * 64) )
        plist.append( ('rhe_patnameff',            c_char * 65) )
        plist.append( ('rhe_patidff',              c_char * 65) )
        plist.append( ('rhe_reqnumff',             c_char * 17) )
        plist.append( ('rhe_dateofbirth',          c_char * 9) )
        plist.append( ('pad_xx',                   c_char * 422) )
        plist.append( ('rhs_position',             c_int) )
        plist.append( ('rhs_entry',                c_int) )
        plist.append( ('pad_xx',                   c_char * 62) )
        plist.append( ('rhs_se_no',                c_short) )
        plist.append( ('pad_xx',                   c_char * 74) )
        plist.append( ('rhs_se_desc',              c_char * 65) )
        plist.append( ('pad_xx',                   c_char * 18) )
        plist.append( ('rhs_anref',                c_char * 3) )
        plist.append( ('pad_xx',                   c_char * 27) )
        plist.append( ('rhs_series_uid',           c_char * 32) )
        plist.append( ('rhs_landmark_uid',         c_char * 32) )
        plist.append( ('pad_xx',                   c_char * 1573) )
        plist.append( ('rhi_dfov',                 c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhi_scanspacing',          c_float) )
        plist.append( ('rhi_loc',                  c_float) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhi_nex',                  c_float) )
        plist.append( ('pad_xx',                   c_char * 20) )
        plist.append( ('rhi_user0',                c_float) )
        plist.append( ('rhi_user1',                c_float) )
        plist.append( ('rhi_user2',                c_float) )
        plist.append( ('rhi_user3',                c_float) )
        plist.append( ('rhi_user4',                c_float) )
        plist.append( ('rhi_user5',                c_float) )
        plist.append( ('rhi_user6',                c_float) )
        plist.append( ('rhi_user7',                c_float) )
        plist.append( ('rhi_user8',                c_float) )
        plist.append( ('rhi_user9',                c_float) )
        plist.append( ('rhi_user10',               c_float) )
        plist.append( ('rhi_user11',               c_float) )
        plist.append( ('rhi_user12',               c_float) )
        plist.append( ('rhi_user13',               c_float) )
        plist.append( ('rhi_user14',               c_float) )
        plist.append( ('rhi_user15',               c_float) )
        plist.append( ('rhi_user16',               c_float) )
        plist.append( ('rhi_user17',               c_float) )
        plist.append( ('rhi_user18',               c_float) )
        plist.append( ('rhi_user19',               c_float) )
        plist.append( ('rhi_user20',               c_float) )
        plist.append( ('rhi_user21',               c_float) )
        plist.append( ('rhi_user22',               c_float) )
        plist.append( ('pad_xx',                   c_char * 8) )
        plist.append( ('rhi_user23',               c_float) )
        plist.append( ('rhi_user24',               c_float) )
        plist.append( ('pad_xx',                   c_char * 60) )
        plist.append( ('rhi_user25',               c_float) )
        plist.append( ('rhi_user26',               c_float) )
        plist.append( ('rhi_user27',               c_float) )
        plist.append( ('rhi_user28',               c_float) )
        plist.append( ('rhi_user29',               c_float) )
        plist.append( ('rhi_user30',               c_float) )
        plist.append( ('rhi_user31',               c_float) )
        plist.append( ('rhi_user32',               c_float) )
        plist.append( ('rhi_user33',               c_float) )
        plist.append( ('rhi_user34',               c_float) )
        plist.append( ('rhi_user35',               c_float) )
        plist.append( ('rhi_user36',               c_float) )
        plist.append( ('rhi_user37',               c_float) )
        plist.append( ('rhi_user38',               c_float) )
        plist.append( ('rhi_user39',               c_float) )
        plist.append( ('rhi_user40',               c_float) )
        plist.append( ('rhi_user41',               c_float) )
        plist.append( ('rhi_user42',               c_float) )
        plist.append( ('rhi_user43',               c_float) )
        plist.append( ('rhi_user44',               c_float) )
        plist.append( ('rhi_user45',               c_float) )
        plist.append( ('rhi_user46',               c_float) )
        plist.append( ('rhi_user47',               c_float) )
        plist.append( ('rhi_user48',               c_float) )
        plist.append( ('pad_xx',                   c_char * 76) )
        plist.append( ('rhi_ctr_R',                c_float) )
        plist.append( ('rhi_ctr_A',                c_float) )
        plist.append( ('rhi_ctr_S',                c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhi_tlhc_R',               c_float) )
        plist.append( ('rhi_tlhc_A',               c_float) )
        plist.append( ('rhi_tlhc_S',               c_float) )
        plist.append( ('rhi_trhc_R',               c_float) )
        plist.append( ('rhi_trhc_A',               c_float) )
        plist.append( ('rhi_trhc_S',               c_float) )
        plist.append( ('rhi_brhc_R',               c_float) )
        plist.append( ('rhi_brhc_A',               c_float) )
        plist.append( ('rhi_brhc_S',               c_float) )
        plist.append( ('pad_xx',                   c_char * 36) )
        plist.append( ('rhi_tr',                   c_int) )
        plist.append( ('rhi_ti',                   c_int) )
        plist.append( ('rhi_te',                   c_int) )
        plist.append( ('pad_xx',                   c_char * 170) )
        plist.append( ('rhi_numecho',              c_short) )
        plist.append( ('pad_xx',                   c_char * 36) )
        plist.append( ('rhi_mr_flip',              c_short) )
        plist.append( ('pad_xx',                   c_char * 22) )
        plist.append( ('rhi_ctyp',                 c_short) )
        plist.append( ('pad_xx',                   c_char * 64) )
        plist.append( ('rhi_freq_dir',             c_short) )
        plist.append( ('pad_xx',                   c_char * 38) )
        plist.append( ('rhi_psdname',              c_char * 33) )
        plist.append( ('pad_xx',                   c_char * 84) )
        plist.append( ('rhi_cname',                c_char * 17) )
        plist.append( ('pad_xx',                   c_char * 51) )
        plist.append( ('rhi_image_uid',            c_char * 32) )

    elif version == 14:
        plist.append( ('rhr_rh_rdbm_rev',          c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhr_rh_scan_date',         c_char * 10) )
        plist.append( ('rhr_rh_scan_time',         c_char * 8) )
        plist.append( ('rhr_rh_logo',              c_char * 10) )
        plist.append( ('rhr_rh_file_contents',     c_short) )
        plist.append( ('pad_xx',                   c_char * 10) )
        plist.append( ('rhr_rh_data_collect_type', c_short) )
        plist.append( ('pad_xx',                   c_char * 6) )
        plist.append( ('rhr_rh_npasses',           c_short) )
        plist.append( ('pad_xx',                   c_char * 2) )
        plist.append( ('rhr_rh_nslices',           c_short) )
        plist.append( ('pad_xx',                   c_char * 10) )
        plist.append( ('rhr_rh_frame_size',        c_ushort) )
        plist.append( ('rhr_rh_point_size',        c_short) )
        plist.append( ('pad_xx',                   c_char * 32) )
        plist.append( ('rhr_rh_raw_pass_size',     c_uint) )
        plist.append( ('pad_xx',                   c_char * 80) )
        plist.append( ('rhr_rh_dab[0]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[0]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[1]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[1]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[2]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[2]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[3]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[3]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_user0',             c_float) )
        plist.append( ('rhr_rh_user1',             c_float) )
        plist.append( ('rhr_rh_user2',             c_float) )
        plist.append( ('rhr_rh_user3',             c_float) )
        plist.append( ('rhr_rh_user4',             c_float) )
        plist.append( ('rhr_rh_user5',             c_float) )
        plist.append( ('rhr_rh_user6',             c_float) )
        plist.append( ('rhr_rh_user7',             c_float) )
        plist.append( ('rhr_rh_user8',             c_float) )
        plist.append( ('rhr_rh_user9',             c_float) )
        plist.append( ('rhr_rh_user10',            c_float) )
        plist.append( ('rhr_rh_user11',            c_float) )
        plist.append( ('rhr_rh_user12',            c_float) )
        plist.append( ('rhr_rh_user13',            c_float) )
        plist.append( ('rhr_rh_user14',            c_float) )
        plist.append( ('rhr_rh_user15',            c_float) )
        plist.append( ('rhr_rh_user16',            c_float) )
        plist.append( ('rhr_rh_user17',            c_float) )
        plist.append( ('rhr_rh_user18',            c_float) )
        plist.append( ('rhr_rh_user19',            c_float) )
        plist.append( ('pad_xx',                   c_char * 72) )
        plist.append( ('rhr_spectral_width',       c_float) )
        plist.append( ('rhr_csi_dims',             c_short) )
        plist.append( ('rhr_xcsi',                 c_short) )
        plist.append( ('rhr_ycsi',                 c_short) )
        plist.append( ('rhr_zcsi',                 c_short) )
        plist.append( ('rhr_roilenx',              c_float) )
        plist.append( ('rhr_roileny',              c_float) )
        plist.append( ('rhr_roilenz',              c_float) )
        plist.append( ('pad_xx',                   c_char * 32) )
        plist.append( ('rhr_rh_ps_mps_freq',       c_int) )
        plist.append( ('pad_xx',                   c_char * 560) )
        plist.append( ('rhr_rh_user_usage_tag',    c_uint) )
        plist.append( ('pad_xx',                   c_char * 8) )
        plist.append( ('rhr_rh_user20',            c_float) )
        plist.append( ('rhr_rh_user21',            c_float) )
        plist.append( ('rhr_rh_user22',            c_float) )
        plist.append( ('rhr_rh_user23',            c_float) )
        plist.append( ('rhr_rh_user24',            c_float) )
        plist.append( ('rhr_rh_user25',            c_float) )
        plist.append( ('rhr_rh_user26',            c_float) )
        plist.append( ('rhr_rh_user27',            c_float) )
        plist.append( ('rhr_rh_user28',            c_float) )
        plist.append( ('rhr_rh_user29',            c_float) )
        plist.append( ('rhr_rh_user30',            c_float) )
        plist.append( ('rhr_rh_user31',            c_float) )
        plist.append( ('rhr_rh_user32',            c_float) )
        plist.append( ('rhr_rh_user33',            c_float) )
        plist.append( ('rhr_rh_user34',            c_float) )
        plist.append( ('rhr_rh_user35',            c_float) )
        plist.append( ('rhr_rh_user36',            c_float) )
        plist.append( ('rhr_rh_user37',            c_float) )
        plist.append( ('rhr_rh_user38',            c_float) )
        plist.append( ('rhr_rh_user39',            c_float) )
        plist.append( ('rhr_rh_user40',            c_float) )
        plist.append( ('rhr_rh_user41',            c_float) )
        plist.append( ('rhr_rh_user42',            c_float) )
        plist.append( ('rhr_rh_user43',            c_float) )
        plist.append( ('rhr_rh_user44',            c_float) )
        plist.append( ('rhr_rh_user45',            c_float) )
        plist.append( ('rhr_rh_user46',            c_float) )
        plist.append( ('rhr_rh_user47',            c_float) )
        plist.append( ('rhr_rh_user48',            c_float) )
        plist.append( ('pad_xx',                   c_char * 352) )
        plist.append( ('rhr_rdb_hdr_off_data',     c_int) )
        plist.append( ('pad_xx',                   c_char * 139508) )
        plist.append( ('rhe_magstrength',          c_int) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhe_ex_datetime',          c_int) )
        plist.append( ('pad_xx',                   c_char * 52) )
        plist.append( ('rhe_ex_no',                c_ushort) )
        plist.append( ('pad_xx',                   c_char * 22) )
        plist.append( ('rhe_patsex',               c_short) )
        plist.append( ('pad_xx',                   c_char * 91) )
        plist.append( ('rhe_reqnum',               c_char * 13) )
        plist.append( ('rhe_refphy',               c_char * 33) )
        plist.append( ('pad_xx',                   c_char * 105) )
        plist.append( ('rhe_ex_sysid',             c_char * 9) )
        plist.append( ('pad_xx',                   c_char * 14) )
        plist.append( ('rhe_hospname',             c_char * 33) )
        plist.append( ('rhe_patid',                c_char * 13) )
        plist.append( ('rhe_patname',              c_char * 25) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhe_ex_verscre',           c_char * 2) )
        plist.append( ('pad_xx',                   c_char * 2) )
        plist.append( ('rhe_uniq_sys_id',          c_char * 16) )
        plist.append( ('pad_xx',                   c_char * 20) )
        plist.append( ('rhe_study_uid',            c_char * 32) )
        plist.append( ('pad_xx',                   c_char * 64) )
        plist.append( ('rhe_patnameff',            c_char * 65) )
        plist.append( ('rhe_patidff',              c_char * 65) )
        plist.append( ('rhe_reqnumff',             c_char * 17) )
        plist.append( ('rhe_dateofbirth',          c_char * 9) )
        plist.append( ('pad_xx',                   c_char * 358) )
        plist.append( ('rhs_position',             c_int) )
        plist.append( ('rhs_entry',                c_int) )
        plist.append( ('pad_xx',                   c_char * 126) )
        plist.append( ('rhs_se_no',                c_short) )
        plist.append( ('pad_xx',                   c_char * 122) )
        plist.append( ('rhs_se_desc',              c_char * 65) )
        plist.append( ('pad_xx',                   c_char * 18) )
        plist.append( ('rhs_anref',                c_char * 3) )
        plist.append( ('pad_xx',                   c_char * 27) )
        plist.append( ('rhs_series_uid',           c_char * 32) )
        plist.append( ('rhs_landmark_uid',         c_char * 32) )
        plist.append( ('pad_xx',                   c_char * 1429) )
        plist.append( ('rhi_dfov',                 c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhi_scanspacing',          c_float) )
        plist.append( ('rhi_loc',                  c_float) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhi_nex',                  c_float) )
        plist.append( ('pad_xx',                   c_char * 20) )
        plist.append( ('rhi_user0',                c_float) )
        plist.append( ('rhi_user1',                c_float) )
        plist.append( ('rhi_user2',                c_float) )
        plist.append( ('rhi_user3',                c_float) )
        plist.append( ('rhi_user4',                c_float) )
        plist.append( ('rhi_user5',                c_float) )
        plist.append( ('rhi_user6',                c_float) )
        plist.append( ('rhi_user7',                c_float) )
        plist.append( ('rhi_user8',                c_float) )
        plist.append( ('rhi_user9',                c_float) )
        plist.append( ('rhi_user10',               c_float) )
        plist.append( ('rhi_user11',               c_float) )
        plist.append( ('rhi_user12',               c_float) )
        plist.append( ('rhi_user13',               c_float) )
        plist.append( ('rhi_user14',               c_float) )
        plist.append( ('rhi_user15',               c_float) )
        plist.append( ('rhi_user16',               c_float) )
        plist.append( ('rhi_user17',               c_float) )
        plist.append( ('rhi_user18',               c_float) )
        plist.append( ('rhi_user19',               c_float) )
        plist.append( ('rhi_user20',               c_float) )
        plist.append( ('rhi_user21',               c_float) )
        plist.append( ('rhi_user22',               c_float) )
        plist.append( ('pad_xx',                   c_char * 8) )
        plist.append( ('rhi_user23',               c_float) )
        plist.append( ('rhi_user24',               c_float) )
        plist.append( ('pad_xx',                   c_char * 60) )
        plist.append( ('rhi_user25',               c_float) )
        plist.append( ('rhi_user26',               c_float) )
        plist.append( ('rhi_user27',               c_float) )
        plist.append( ('rhi_user28',               c_float) )
        plist.append( ('rhi_user29',               c_float) )
        plist.append( ('rhi_user30',               c_float) )
        plist.append( ('rhi_user31',               c_float) )
        plist.append( ('rhi_user32',               c_float) )
        plist.append( ('rhi_user33',               c_float) )
        plist.append( ('rhi_user34',               c_float) )
        plist.append( ('rhi_user35',               c_float) )
        plist.append( ('rhi_user36',               c_float) )
        plist.append( ('rhi_user37',               c_float) )
        plist.append( ('rhi_user38',               c_float) )
        plist.append( ('rhi_user39',               c_float) )
        plist.append( ('rhi_user40',               c_float) )
        plist.append( ('rhi_user41',               c_float) )
        plist.append( ('rhi_user42',               c_float) )
        plist.append( ('rhi_user43',               c_float) )
        plist.append( ('rhi_user44',               c_float) )
        plist.append( ('rhi_user45',               c_float) )
        plist.append( ('rhi_user46',               c_float) )
        plist.append( ('rhi_user47',               c_float) )
        plist.append( ('rhi_user48',               c_float) )
        plist.append( ('pad_xx',                   c_char * 76) )
        plist.append( ('rhi_ctr_R',                c_float) )
        plist.append( ('rhi_ctr_A',                c_float) )
        plist.append( ('rhi_ctr_S',                c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhi_tlhc_R',               c_float) )
        plist.append( ('rhi_tlhc_A',               c_float) )
        plist.append( ('rhi_tlhc_S',               c_float) )
        plist.append( ('rhi_trhc_R',               c_float) )
        plist.append( ('rhi_trhc_A',               c_float) )
        plist.append( ('rhi_trhc_S',               c_float) )
        plist.append( ('rhi_brhc_R',               c_float) )
        plist.append( ('rhi_brhc_A',               c_float) )
        plist.append( ('rhi_brhc_S',               c_float) )
        plist.append( ('pad_xx',                   c_char * 196) )
        plist.append( ('rhi_tr',                   c_int) )
        plist.append( ('rhi_ti',                   c_int) )
        plist.append( ('rhi_te',                   c_int) )
        plist.append( ('pad_xx',                   c_char * 306) )
        plist.append( ('rhi_numecho',              c_short) )
        plist.append( ('pad_xx',                   c_char * 36) )
        plist.append( ('rhi_mr_flip',              c_short) )
        plist.append( ('pad_xx',                   c_char * 22) )
        plist.append( ('rhi_ctyp',                 c_short) )
        plist.append( ('pad_xx',                   c_char * 64) )
        plist.append( ('rhi_freq_dir',             c_short) )
        plist.append( ('pad_xx',                   c_char * 112) )
        plist.append( ('rhi_psdname',              c_char * 33) )
        plist.append( ('pad_xx',                   c_char * 84) )
        plist.append( ('rhi_cname',                c_char * 17) )
        plist.append( ('pad_xx',                   c_char * 51) )
        plist.append( ('rhi_image_uid',            c_char * 32) )

    elif version == 15:
        plist.append( ('rhr_rh_rdbm_rev',          c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhr_rh_scan_date',         c_char * 10) )
        plist.append( ('rhr_rh_scan_time',         c_char * 8) )
        plist.append( ('rhr_rh_logo',              c_char * 10) )
        plist.append( ('rhr_rh_file_contents',     c_short) )
        plist.append( ('pad_xx',                   c_char * 10) )
        plist.append( ('rhr_rh_data_collect_type', c_short) )
        plist.append( ('pad_xx',                   c_char * 6) )
        plist.append( ('rhr_rh_npasses',           c_short) )
        plist.append( ('pad_xx',                   c_char * 2) )
        plist.append( ('rhr_rh_nslices',           c_short) )
        plist.append( ('pad_xx',                   c_char * 10) )
        plist.append( ('rhr_rh_frame_size',        c_ushort) )
        plist.append( ('rhr_rh_point_size',        c_short) )
        plist.append( ('pad_xx',                   c_char * 32) )
        plist.append( ('rhr_rh_raw_pass_size',     c_uint) )
        plist.append( ('pad_xx',                   c_char * 80) )
        plist.append( ('rhr_rh_dab[0]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[0]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[1]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[1]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[2]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[2]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[3]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[3]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_user0',             c_float) )
        plist.append( ('rhr_rh_user1',             c_float) )
        plist.append( ('rhr_rh_user2',             c_float) )
        plist.append( ('rhr_rh_user3',             c_float) )
        plist.append( ('rhr_rh_user4',             c_float) )
        plist.append( ('rhr_rh_user5',             c_float) )
        plist.append( ('rhr_rh_user6',             c_float) )
        plist.append( ('rhr_rh_user7',             c_float) )
        plist.append( ('rhr_rh_user8',             c_float) )
        plist.append( ('rhr_rh_user9',             c_float) )
        plist.append( ('rhr_rh_user10',            c_float) )
        plist.append( ('rhr_rh_user11',            c_float) )
        plist.append( ('rhr_rh_user12',            c_float) )
        plist.append( ('rhr_rh_user13',            c_float) )
        plist.append( ('rhr_rh_user14',            c_float) )
        plist.append( ('rhr_rh_user15',            c_float) )
        plist.append( ('rhr_rh_user16',            c_float) )
        plist.append( ('rhr_rh_user17',            c_float) )
        plist.append( ('rhr_rh_user18',            c_float) )
        plist.append( ('rhr_rh_user19',            c_float) )
        plist.append( ('pad_xx',                   c_char * 72) )
        plist.append( ('rhr_spectral_width',       c_float) )
        plist.append( ('rhr_csi_dims',             c_short) )
        plist.append( ('rhr_xcsi',                 c_short) )
        plist.append( ('rhr_ycsi',                 c_short) )
        plist.append( ('rhr_zcsi',                 c_short) )
        plist.append( ('rhr_roilenx',              c_float) )
        plist.append( ('rhr_roileny',              c_float) )
        plist.append( ('rhr_roilenz',              c_float) )
        plist.append( ('pad_xx',                   c_char * 32) )
        plist.append( ('rhr_rh_ps_mps_freq',       c_int) )
        plist.append( ('pad_xx',                   c_char * 560) )
        plist.append( ('rhr_rh_user_usage_tag',    c_uint) )
        plist.append( ('pad_xx',                   c_char * 8) )
        plist.append( ('rhr_rh_user20',            c_float) )
        plist.append( ('rhr_rh_user21',            c_float) )
        plist.append( ('rhr_rh_user22',            c_float) )
        plist.append( ('rhr_rh_user23',            c_float) )
        plist.append( ('rhr_rh_user24',            c_float) )
        plist.append( ('rhr_rh_user25',            c_float) )
        plist.append( ('rhr_rh_user26',            c_float) )
        plist.append( ('rhr_rh_user27',            c_float) )
        plist.append( ('rhr_rh_user28',            c_float) )
        plist.append( ('rhr_rh_user29',            c_float) )
        plist.append( ('rhr_rh_user30',            c_float) )
        plist.append( ('rhr_rh_user31',            c_float) )
        plist.append( ('rhr_rh_user32',            c_float) )
        plist.append( ('rhr_rh_user33',            c_float) )
        plist.append( ('rhr_rh_user34',            c_float) )
        plist.append( ('rhr_rh_user35',            c_float) )
        plist.append( ('rhr_rh_user36',            c_float) )
        plist.append( ('rhr_rh_user37',            c_float) )
        plist.append( ('rhr_rh_user38',            c_float) )
        plist.append( ('rhr_rh_user39',            c_float) )
        plist.append( ('rhr_rh_user40',            c_float) )
        plist.append( ('rhr_rh_user41',            c_float) )
        plist.append( ('rhr_rh_user42',            c_float) )
        plist.append( ('rhr_rh_user43',            c_float) )
        plist.append( ('rhr_rh_user44',            c_float) )
        plist.append( ('rhr_rh_user45',            c_float) )
        plist.append( ('rhr_rh_user46',            c_float) )
        plist.append( ('rhr_rh_user47',            c_float) )
        plist.append( ('rhr_rh_user48',            c_float) )
        plist.append( ('pad_xx',                   c_char * 352) )
        plist.append( ('rhr_rdb_hdr_off_data',     c_int) )
        plist.append( ('pad_xx',                   c_char * 139508) )
        plist.append( ('rhe_magstrength',          c_int) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhe_ex_datetime',          c_int) )
        plist.append( ('pad_xx',                   c_char * 52) )
        plist.append( ('rhe_ex_no',                c_ushort) )
        plist.append( ('pad_xx',                   c_char * 22) )
        plist.append( ('rhe_patsex',               c_short) )
        plist.append( ('pad_xx',                   c_char * 91) )
        plist.append( ('rhe_reqnum',               c_char * 13) )
        plist.append( ('rhe_refphy',               c_char * 33) )
        plist.append( ('pad_xx',                   c_char * 105) )
        plist.append( ('rhe_ex_sysid',             c_char * 9) )
        plist.append( ('pad_xx',                   c_char * 14) )
        plist.append( ('rhe_hospname',             c_char * 33) )
        plist.append( ('rhe_patid',                c_char * 13) )
        plist.append( ('rhe_patname',              c_char * 25) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhe_ex_verscre',           c_char * 2) )
        plist.append( ('pad_xx',                   c_char * 2) )
        plist.append( ('rhe_uniq_sys_id',          c_char * 16) )
        plist.append( ('pad_xx',                   c_char * 20) )
        plist.append( ('rhe_study_uid',            c_char * 32) )
        plist.append( ('pad_xx',                   c_char * 64) )
        plist.append( ('rhe_patnameff',            c_char * 65) )
        plist.append( ('rhe_patidff',              c_char * 65) )
        plist.append( ('rhe_reqnumff',             c_char * 17) )
        plist.append( ('rhe_dateofbirth',          c_char * 9) )
        plist.append( ('pad_xx',                   c_char * 358) )
        plist.append( ('rhs_position',             c_int) )
        plist.append( ('rhs_entry',                c_int) )
        plist.append( ('pad_xx',                   c_char * 126) )
        plist.append( ('rhs_se_no',                c_short) )
        plist.append( ('pad_xx',                   c_char * 122) )
        plist.append( ('rhs_se_desc',              c_char * 65) )
        plist.append( ('pad_xx',                   c_char * 18) )
        plist.append( ('rhs_anref',                c_char * 3) )
        plist.append( ('pad_xx',                   c_char * 27) )
        plist.append( ('rhs_series_uid',           c_char * 32) )
        plist.append( ('rhs_landmark_uid',         c_char * 32) )
        plist.append( ('pad_xx',                   c_char * 1429) )
        plist.append( ('rhi_dfov',                 c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhi_scanspacing',          c_float) )
        plist.append( ('rhi_loc',                  c_float) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhi_nex',                  c_float) )
        plist.append( ('pad_xx',                   c_char * 20) )
        plist.append( ('rhi_user0',                c_float) )
        plist.append( ('rhi_user1',                c_float) )
        plist.append( ('rhi_user2',                c_float) )
        plist.append( ('rhi_user3',                c_float) )
        plist.append( ('rhi_user4',                c_float) )
        plist.append( ('rhi_user5',                c_float) )
        plist.append( ('rhi_user6',                c_float) )
        plist.append( ('rhi_user7',                c_float) )
        plist.append( ('rhi_user8',                c_float) )
        plist.append( ('rhi_user9',                c_float) )
        plist.append( ('rhi_user10',               c_float) )
        plist.append( ('rhi_user11',               c_float) )
        plist.append( ('rhi_user12',               c_float) )
        plist.append( ('rhi_user13',               c_float) )
        plist.append( ('rhi_user14',               c_float) )
        plist.append( ('rhi_user15',               c_float) )
        plist.append( ('rhi_user16',               c_float) )
        plist.append( ('rhi_user17',               c_float) )
        plist.append( ('rhi_user18',               c_float) )
        plist.append( ('rhi_user19',               c_float) )
        plist.append( ('rhi_user20',               c_float) )
        plist.append( ('rhi_user21',               c_float) )
        plist.append( ('rhi_user22',               c_float) )
        plist.append( ('pad_xx',                   c_char * 8) )
        plist.append( ('rhi_user23',               c_float) )
        plist.append( ('rhi_user24',               c_float) )
        plist.append( ('pad_xx',                   c_char * 60) )
        plist.append( ('rhi_user25',               c_float) )
        plist.append( ('rhi_user26',               c_float) )
        plist.append( ('rhi_user27',               c_float) )
        plist.append( ('rhi_user28',               c_float) )
        plist.append( ('rhi_user29',               c_float) )
        plist.append( ('rhi_user30',               c_float) )
        plist.append( ('rhi_user31',               c_float) )
        plist.append( ('rhi_user32',               c_float) )
        plist.append( ('rhi_user33',               c_float) )
        plist.append( ('rhi_user34',               c_float) )
        plist.append( ('rhi_user35',               c_float) )
        plist.append( ('rhi_user36',               c_float) )
        plist.append( ('rhi_user37',               c_float) )
        plist.append( ('rhi_user38',               c_float) )
        plist.append( ('rhi_user39',               c_float) )
        plist.append( ('rhi_user40',               c_float) )
        plist.append( ('rhi_user41',               c_float) )
        plist.append( ('rhi_user42',               c_float) )
        plist.append( ('rhi_user43',               c_float) )
        plist.append( ('rhi_user44',               c_float) )
        plist.append( ('rhi_user45',               c_float) )
        plist.append( ('rhi_user46',               c_float) )
        plist.append( ('rhi_user47',               c_float) )
        plist.append( ('rhi_user48',               c_float) )
        plist.append( ('pad_xx',                   c_char * 76) )
        plist.append( ('rhi_ctr_R',                c_float) )
        plist.append( ('rhi_ctr_A',                c_float) )
        plist.append( ('rhi_ctr_S',                c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhi_tlhc_R',               c_float) )
        plist.append( ('rhi_tlhc_A',               c_float) )
        plist.append( ('rhi_tlhc_S',               c_float) )
        plist.append( ('rhi_trhc_R',               c_float) )
        plist.append( ('rhi_trhc_A',               c_float) )
        plist.append( ('rhi_trhc_S',               c_float) )
        plist.append( ('rhi_brhc_R',               c_float) )
        plist.append( ('rhi_brhc_A',               c_float) )
        plist.append( ('rhi_brhc_S',               c_float) )
        plist.append( ('pad_xx',                   c_char * 196) )
        plist.append( ('rhi_tr',                   c_int) )
        plist.append( ('rhi_ti',                   c_int) )
        plist.append( ('rhi_te',                   c_int) )
        plist.append( ('pad_xx',                   c_char * 306) )
        plist.append( ('rhi_numecho',              c_short) )
        plist.append( ('pad_xx',                   c_char * 32) )
        plist.append( ('rhi_mr_flip',              c_short) )
        plist.append( ('pad_xx',                   c_char * 22) )
        plist.append( ('rhi_ctyp',                 c_short) )
        plist.append( ('pad_xx',                   c_char * 64) )
        plist.append( ('rhi_freq_dir',             c_short) )
        plist.append( ('pad_xx',                   c_char * 116) )
        plist.append( ('rhi_psdname',              c_char * 33) )
        plist.append( ('pad_xx',                   c_char * 84) )
        plist.append( ('rhi_cname',                c_char * 17) )
        plist.append( ('pad_xx',                   c_char * 51) )
        plist.append( ('rhi_image_uid',            c_char * 32) )

    elif version == 16:
        plist.append( ('rhr_rh_rdbm_rev',          c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhr_rh_scan_date',         c_char * 10) )
        plist.append( ('rhr_rh_scan_time',         c_char * 8) )
        plist.append( ('rhr_rh_logo',              c_char * 10) )
        plist.append( ('rhr_rh_file_contents',     c_short) )
        plist.append( ('pad_xx',                   c_char * 10) )
        plist.append( ('rhr_rh_data_collect_type', c_short) )
        plist.append( ('pad_xx',                   c_char * 6) )
        plist.append( ('rhr_rh_npasses',           c_short) )
        plist.append( ('pad_xx',                   c_char * 2) )
        plist.append( ('rhr_rh_nslices',           c_short) )
        plist.append( ('pad_xx',                   c_char * 10) )
        plist.append( ('rhr_rh_frame_size',        c_ushort) )
        plist.append( ('rhr_rh_point_size',        c_short) )
        plist.append( ('pad_xx',                   c_char * 32) )
        plist.append( ('rhr_rh_raw_pass_size',     c_uint) )
        plist.append( ('pad_xx',                   c_char * 80) )
        plist.append( ('rhr_rh_dab[0]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[0]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[1]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[1]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[2]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[2]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[3]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[3]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_user0',             c_float) )
        plist.append( ('rhr_rh_user1',             c_float) )
        plist.append( ('rhr_rh_user2',             c_float) )
        plist.append( ('rhr_rh_user3',             c_float) )
        plist.append( ('rhr_rh_user4',             c_float) )
        plist.append( ('rhr_rh_user5',             c_float) )
        plist.append( ('rhr_rh_user6',             c_float) )
        plist.append( ('rhr_rh_user7',             c_float) )
        plist.append( ('rhr_rh_user8',             c_float) )
        plist.append( ('rhr_rh_user9',             c_float) )
        plist.append( ('rhr_rh_user10',            c_float) )
        plist.append( ('rhr_rh_user11',            c_float) )
        plist.append( ('rhr_rh_user12',            c_float) )
        plist.append( ('rhr_rh_user13',            c_float) )
        plist.append( ('rhr_rh_user14',            c_float) )
        plist.append( ('rhr_rh_user15',            c_float) )
        plist.append( ('rhr_rh_user16',            c_float) )
        plist.append( ('rhr_rh_user17',            c_float) )
        plist.append( ('rhr_rh_user18',            c_float) )
        plist.append( ('rhr_rh_user19',            c_float) )
        plist.append( ('pad_xx',                   c_char * 72) )
        plist.append( ('rhr_spectral_width',       c_float) )
        plist.append( ('rhr_csi_dims',             c_short) )
        plist.append( ('rhr_xcsi',                 c_short) )
        plist.append( ('rhr_ycsi',                 c_short) )
        plist.append( ('rhr_zcsi',                 c_short) )
        plist.append( ('rhr_roilenx',              c_float) )
        plist.append( ('rhr_roileny',              c_float) )
        plist.append( ('rhr_roilenz',              c_float) )
        plist.append( ('pad_xx',                   c_char * 32) )
        plist.append( ('rhr_rh_ps_mps_freq',       c_int) )
        plist.append( ('pad_xx',                   c_char * 560) )
        plist.append( ('rhr_rh_user_usage_tag',    c_uint) )
        plist.append( ('pad_xx',                   c_char * 8) )
        plist.append( ('rhr_rh_user20',            c_float) )
        plist.append( ('rhr_rh_user21',            c_float) )
        plist.append( ('rhr_rh_user22',            c_float) )
        plist.append( ('rhr_rh_user23',            c_float) )
        plist.append( ('rhr_rh_user24',            c_float) )
        plist.append( ('rhr_rh_user25',            c_float) )
        plist.append( ('rhr_rh_user26',            c_float) )
        plist.append( ('rhr_rh_user27',            c_float) )
        plist.append( ('rhr_rh_user28',            c_float) )
        plist.append( ('rhr_rh_user29',            c_float) )
        plist.append( ('rhr_rh_user30',            c_float) )
        plist.append( ('rhr_rh_user31',            c_float) )
        plist.append( ('rhr_rh_user32',            c_float) )
        plist.append( ('rhr_rh_user33',            c_float) )
        plist.append( ('rhr_rh_user34',            c_float) )
        plist.append( ('rhr_rh_user35',            c_float) )
        plist.append( ('rhr_rh_user36',            c_float) )
        plist.append( ('rhr_rh_user37',            c_float) )
        plist.append( ('rhr_rh_user38',            c_float) )
        plist.append( ('rhr_rh_user39',            c_float) )
        plist.append( ('rhr_rh_user40',            c_float) )
        plist.append( ('rhr_rh_user41',            c_float) )
        plist.append( ('rhr_rh_user42',            c_float) )
        plist.append( ('rhr_rh_user43',            c_float) )
        plist.append( ('rhr_rh_user44',            c_float) )
        plist.append( ('rhr_rh_user45',            c_float) )
        plist.append( ('rhr_rh_user46',            c_float) )
        plist.append( ('rhr_rh_user47',            c_float) )
        plist.append( ('rhr_rh_user48',            c_float) )
        plist.append( ('pad_xx',                   c_char * 352) )
        plist.append( ('rhr_rdb_hdr_off_data',     c_int) )
        plist.append( ('pad_xx',                   c_char * 139508) )
        plist.append( ('rhe_magstrength',          c_int) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhe_ex_datetime',          c_int) )
        plist.append( ('pad_xx',                   c_char * 52) )
        plist.append( ('rhe_ex_no',                c_ushort) )
        plist.append( ('pad_xx',                   c_char * 22) )
        plist.append( ('rhe_patsex',               c_short) )
        plist.append( ('pad_xx',                   c_char * 91) )
        plist.append( ('rhe_reqnum',               c_char * 13) )
        plist.append( ('rhe_refphy',               c_char * 33) )
        plist.append( ('pad_xx',                   c_char * 105) )
        plist.append( ('rhe_ex_sysid',             c_char * 9) )
        plist.append( ('pad_xx',                   c_char * 14) )
        plist.append( ('rhe_hospname',             c_char * 33) )
        plist.append( ('rhe_patid',                c_char * 13) )
        plist.append( ('rhe_patname',              c_char * 25) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhe_ex_verscre',           c_char * 2) )
        plist.append( ('pad_xx',                   c_char * 2) )
        plist.append( ('rhe_uniq_sys_id',          c_char * 16) )
        plist.append( ('pad_xx',                   c_char * 20) )
        plist.append( ('rhe_study_uid',            c_char * 32) )
        plist.append( ('pad_xx',                   c_char * 64) )
        plist.append( ('rhe_patnameff',            c_char * 65) )
        plist.append( ('rhe_patidff',              c_char * 65) )
        plist.append( ('rhe_reqnumff',             c_char * 17) )
        plist.append( ('rhe_dateofbirth',          c_char * 9) )
        plist.append( ('pad_xx',                   c_char * 358) )
        plist.append( ('rhs_position',             c_int) )
        plist.append( ('rhs_entry',                c_int) )
        plist.append( ('pad_xx',                   c_char * 126) )
        plist.append( ('rhs_se_no',                c_short) )
        plist.append( ('pad_xx',                   c_char * 122) )
        plist.append( ('rhs_se_desc',              c_char * 65) )
        plist.append( ('pad_xx',                   c_char * 18) )
        plist.append( ('rhs_anref',                c_char * 3) )
        plist.append( ('pad_xx',                   c_char * 27) )
        plist.append( ('rhs_series_uid',           c_char * 32) )
        plist.append( ('rhs_landmark_uid',         c_char * 32) )
        plist.append( ('pad_xx',                   c_char * 1477) )
        plist.append( ('rhi_dfov',                 c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhi_scanspacing',          c_float) )
        plist.append( ('rhi_loc',                  c_float) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhi_nex',                  c_float) )
        plist.append( ('pad_xx',                   c_char * 20) )
        plist.append( ('rhi_user0',                c_float) )
        plist.append( ('rhi_user1',                c_float) )
        plist.append( ('rhi_user2',                c_float) )
        plist.append( ('rhi_user3',                c_float) )
        plist.append( ('rhi_user4',                c_float) )
        plist.append( ('rhi_user5',                c_float) )
        plist.append( ('rhi_user6',                c_float) )
        plist.append( ('rhi_user7',                c_float) )
        plist.append( ('rhi_user8',                c_float) )
        plist.append( ('rhi_user9',                c_float) )
        plist.append( ('rhi_user10',               c_float) )
        plist.append( ('rhi_user11',               c_float) )
        plist.append( ('rhi_user12',               c_float) )
        plist.append( ('rhi_user13',               c_float) )
        plist.append( ('rhi_user14',               c_float) )
        plist.append( ('rhi_user15',               c_float) )
        plist.append( ('rhi_user16',               c_float) )
        plist.append( ('rhi_user17',               c_float) )
        plist.append( ('rhi_user18',               c_float) )
        plist.append( ('rhi_user19',               c_float) )
        plist.append( ('rhi_user20',               c_float) )
        plist.append( ('rhi_user21',               c_float) )
        plist.append( ('rhi_user22',               c_float) )
        plist.append( ('pad_xx',                   c_char * 8) )
        plist.append( ('rhi_user23',               c_float) )
        plist.append( ('rhi_user24',               c_float) )
        plist.append( ('pad_xx',                   c_char * 60) )
        plist.append( ('rhi_user25',               c_float) )
        plist.append( ('rhi_user26',               c_float) )
        plist.append( ('rhi_user27',               c_float) )
        plist.append( ('rhi_user28',               c_float) )
        plist.append( ('rhi_user29',               c_float) )
        plist.append( ('rhi_user30',               c_float) )
        plist.append( ('rhi_user31',               c_float) )
        plist.append( ('rhi_user32',               c_float) )
        plist.append( ('rhi_user33',               c_float) )
        plist.append( ('rhi_user34',               c_float) )
        plist.append( ('rhi_user35',               c_float) )
        plist.append( ('rhi_user36',               c_float) )
        plist.append( ('rhi_user37',               c_float) )
        plist.append( ('rhi_user38',               c_float) )
        plist.append( ('rhi_user39',               c_float) )
        plist.append( ('rhi_user40',               c_float) )
        plist.append( ('rhi_user41',               c_float) )
        plist.append( ('rhi_user42',               c_float) )
        plist.append( ('rhi_user43',               c_float) )
        plist.append( ('rhi_user44',               c_float) )
        plist.append( ('rhi_user45',               c_float) )
        plist.append( ('rhi_user46',               c_float) )
        plist.append( ('rhi_user47',               c_float) )
        plist.append( ('rhi_user48',               c_float) )
        plist.append( ('pad_xx',                   c_char * 76) )
        plist.append( ('rhi_ctr_R',                c_float) )
        plist.append( ('rhi_ctr_A',                c_float) )
        plist.append( ('rhi_ctr_S',                c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhi_tlhc_R',               c_float) )
        plist.append( ('rhi_tlhc_A',               c_float) )
        plist.append( ('rhi_tlhc_S',               c_float) )
        plist.append( ('rhi_trhc_R',               c_float) )
        plist.append( ('rhi_trhc_A',               c_float) )
        plist.append( ('rhi_trhc_S',               c_float) )
        plist.append( ('rhi_brhc_R',               c_float) )
        plist.append( ('rhi_brhc_A',               c_float) )
        plist.append( ('rhi_brhc_S',               c_float) )
        plist.append( ('pad_xx',                   c_char * 196) )
        plist.append( ('rhi_tr',                   c_int) )
        plist.append( ('rhi_ti',                   c_int) )
        plist.append( ('rhi_te',                   c_int) )
        plist.append( ('pad_xx',                   c_char * 306) )
        plist.append( ('rhi_numecho',              c_short) )
        plist.append( ('pad_xx',                   c_char * 32) )
        plist.append( ('rhi_mr_flip',              c_short) )
        plist.append( ('pad_xx',                   c_char * 22) )
        plist.append( ('rhi_ctyp',                 c_short) )
        plist.append( ('pad_xx',                   c_char * 64) )
        plist.append( ('rhi_freq_dir',             c_short) )
        plist.append( ('pad_xx',                   c_char * 116) )
        plist.append( ('rhi_psdname',              c_char * 33) )
        plist.append( ('pad_xx',                   c_char * 84) )
        plist.append( ('rhi_cname',                c_char * 17) )
        plist.append( ('pad_xx',                   c_char * 51) )
        plist.append( ('rhi_image_uid',            c_char * 32) )

    elif version == 20:
        plist.append( ('rhr_rh_rdbm_rev',          c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhr_rh_scan_date',         c_char * 10) )
        plist.append( ('rhr_rh_scan_time',         c_char * 8) )
        plist.append( ('rhr_rh_logo',              c_char * 10) )
        plist.append( ('rhr_rh_file_contents',     c_short) )
        plist.append( ('pad_xx',                   c_char * 10) )
        plist.append( ('rhr_rh_data_collect_type', c_short) )
        plist.append( ('pad_xx',                   c_char * 6) )
        plist.append( ('rhr_rh_npasses',           c_short) )
        plist.append( ('pad_xx',                   c_char * 2) )
        plist.append( ('rhr_rh_nslices',           c_short) )
        plist.append( ('pad_xx',                   c_char * 10) )
        plist.append( ('rhr_rh_frame_size',        c_ushort) )
        plist.append( ('rhr_rh_point_size',        c_short) )
        plist.append( ('pad_xx',                   c_char * 116) )
        plist.append( ('rhr_rh_dab[0]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[0]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[1]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[1]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[2]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[2]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[3]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[3]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_user0',             c_float) )
        plist.append( ('rhr_rh_user1',             c_float) )
        plist.append( ('rhr_rh_user2',             c_float) )
        plist.append( ('rhr_rh_user3',             c_float) )
        plist.append( ('rhr_rh_user4',             c_float) )
        plist.append( ('rhr_rh_user5',             c_float) )
        plist.append( ('rhr_rh_user6',             c_float) )
        plist.append( ('rhr_rh_user7',             c_float) )
        plist.append( ('rhr_rh_user8',             c_float) )
        plist.append( ('rhr_rh_user9',             c_float) )
        plist.append( ('rhr_rh_user10',            c_float) )
        plist.append( ('rhr_rh_user11',            c_float) )
        plist.append( ('rhr_rh_user12',            c_float) )
        plist.append( ('rhr_rh_user13',            c_float) )
        plist.append( ('rhr_rh_user14',            c_float) )
        plist.append( ('rhr_rh_user15',            c_float) )
        plist.append( ('rhr_rh_user16',            c_float) )
        plist.append( ('rhr_rh_user17',            c_float) )
        plist.append( ('rhr_rh_user18',            c_float) )
        plist.append( ('rhr_rh_user19',            c_float) )
        plist.append( ('pad_xx',                   c_char * 72) )
        plist.append( ('rhr_spectral_width',       c_float) )
        plist.append( ('rhr_csi_dims',             c_short) )
        plist.append( ('rhr_xcsi',                 c_short) )
        plist.append( ('rhr_ycsi',                 c_short) )
        plist.append( ('rhr_zcsi',                 c_short) )
        plist.append( ('rhr_roilenx',              c_float) )
        plist.append( ('rhr_roileny',              c_float) )
        plist.append( ('rhr_roilenz',              c_float) )
        plist.append( ('pad_xx',                   c_char * 32) )
        plist.append( ('rhr_rh_ps_mps_freq',       c_uint) )
        plist.append( ('pad_xx',                   c_char * 560) )
        plist.append( ('rhr_rh_user_usage_tag',    c_uint) )
        plist.append( ('pad_xx',                   c_char * 8) )
        plist.append( ('rhr_rh_user20',            c_float) )
        plist.append( ('rhr_rh_user21',            c_float) )
        plist.append( ('rhr_rh_user22',            c_float) )
        plist.append( ('rhr_rh_user23',            c_float) )
        plist.append( ('rhr_rh_user24',            c_float) )
        plist.append( ('rhr_rh_user25',            c_float) )
        plist.append( ('rhr_rh_user26',            c_float) )
        plist.append( ('rhr_rh_user27',            c_float) )
        plist.append( ('rhr_rh_user28',            c_float) )
        plist.append( ('rhr_rh_user29',            c_float) )
        plist.append( ('rhr_rh_user30',            c_float) )
        plist.append( ('rhr_rh_user31',            c_float) )
        plist.append( ('rhr_rh_user32',            c_float) )
        plist.append( ('rhr_rh_user33',            c_float) )
        plist.append( ('rhr_rh_user34',            c_float) )
        plist.append( ('rhr_rh_user35',            c_float) )
        plist.append( ('rhr_rh_user36',            c_float) )
        plist.append( ('rhr_rh_user37',            c_float) )
        plist.append( ('rhr_rh_user38',            c_float) )
        plist.append( ('rhr_rh_user39',            c_float) )
        plist.append( ('rhr_rh_user40',            c_float) )
        plist.append( ('rhr_rh_user41',            c_float) )
        plist.append( ('rhr_rh_user42',            c_float) )
        plist.append( ('rhr_rh_user43',            c_float) )
        plist.append( ('rhr_rh_user44',            c_float) )
        plist.append( ('rhr_rh_user45',            c_float) )
        plist.append( ('rhr_rh_user46',            c_float) )
        plist.append( ('rhr_rh_user47',            c_float) )
        plist.append( ('rhr_rh_user48',            c_float) )
        plist.append( ('pad_xx',                   c_char * 352) )
        plist.append( ('rhr_rdb_hdr_off_data',     c_int) )
        plist.append( ('pad_xx',                   c_char * 188) )
        plist.append( ('rhr_rh_raw_pass_size',     c_longlong) )
        plist.append( ('pad_xx',                   c_char * 141724) )
        plist.append( ('rhe_magstrength',          c_int) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhe_ex_datetime',          c_int) )
        plist.append( ('pad_xx',                   c_char * 112) )
        plist.append( ('rhe_ex_no',                c_ushort) )
        plist.append( ('pad_xx',                   c_char * 22) )
        plist.append( ('rhe_patsex',               c_short) )
        plist.append( ('pad_xx',                   c_char * 335) )
        plist.append( ('rhe_refphy',               c_char * 65) )
        plist.append( ('pad_xx',                   c_char * 198) )
        plist.append( ('rhe_ex_sysid',             c_char * 9) )
        plist.append( ('pad_xx',                   c_char * 14) )
        plist.append( ('rhe_hospname',             c_char * 33) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhe_ex_verscre',           c_char * 2) )
        plist.append( ('pad_xx',                   c_char * 2) )
        plist.append( ('rhe_uniq_sys_id',          c_char * 16) )
        plist.append( ('pad_xx',                   c_char * 20) )
        plist.append( ('rhe_study_uid',            c_char * 32) )
        plist.append( ('pad_xx',                   c_char * 64) )
        plist.append( ('rhe_patname',              c_char * 65) )
        plist.append( ('rhe_patid',                c_char * 65) )
        plist.append( ('rhe_reqnum',               c_char * 17) )
        plist.append( ('rhe_dateofbirth',          c_char * 9) )
        plist.append( ('pad_xx',                   c_char * 928) )
        plist.append( ('rhs_position',             c_int) )
        plist.append( ('rhs_entry',                c_int) )
        plist.append( ('pad_xx',                   c_char * 194) )
        plist.append( ('rhs_se_no',                c_short) )
        plist.append( ('pad_xx',                   c_char * 138) )
        plist.append( ('rhs_se_desc',              c_char * 65) )
        plist.append( ('pad_xx',                   c_char * 18) )
        plist.append( ('rhs_anref',                c_char * 3) )
        plist.append( ('pad_xx',                   c_char * 27) )
        plist.append( ('rhs_series_uid',           c_char * 32) )
        plist.append( ('rhs_landmark_uid',         c_char * 32) )
        plist.append( ('pad_xx',                   c_char * 1705) )
        plist.append( ('rhi_dfov',                 c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhi_scanspacing',          c_float) )
        plist.append( ('rhi_loc',                  c_float) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhi_nex',                  c_float) )
        plist.append( ('pad_xx',                   c_char * 20) )
        plist.append( ('rhi_user0',                c_float) )
        plist.append( ('rhi_user1',                c_float) )
        plist.append( ('rhi_user2',                c_float) )
        plist.append( ('rhi_user3',                c_float) )
        plist.append( ('rhi_user4',                c_float) )
        plist.append( ('rhi_user5',                c_float) )
        plist.append( ('rhi_user6',                c_float) )
        plist.append( ('rhi_user7',                c_float) )
        plist.append( ('rhi_user8',                c_float) )
        plist.append( ('rhi_user9',                c_float) )
        plist.append( ('rhi_user10',               c_float) )
        plist.append( ('rhi_user11',               c_float) )
        plist.append( ('rhi_user12',               c_float) )
        plist.append( ('rhi_user13',               c_float) )
        plist.append( ('rhi_user14',               c_float) )
        plist.append( ('rhi_user15',               c_float) )
        plist.append( ('rhi_user16',               c_float) )
        plist.append( ('rhi_user17',               c_float) )
        plist.append( ('rhi_user18',               c_float) )
        plist.append( ('rhi_user19',               c_float) )
        plist.append( ('rhi_user20',               c_float) )
        plist.append( ('rhi_user21',               c_float) )
        plist.append( ('rhi_user22',               c_float) )
        plist.append( ('pad_xx',                   c_char * 8) )
        plist.append( ('rhi_user23',               c_float) )
        plist.append( ('rhi_user24',               c_float) )
        plist.append( ('pad_xx',                   c_char * 60) )
        plist.append( ('rhi_user25',               c_float) )
        plist.append( ('rhi_user26',               c_float) )
        plist.append( ('rhi_user27',               c_float) )
        plist.append( ('rhi_user28',               c_float) )
        plist.append( ('rhi_user29',               c_float) )
        plist.append( ('rhi_user30',               c_float) )
        plist.append( ('rhi_user31',               c_float) )
        plist.append( ('rhi_user32',               c_float) )
        plist.append( ('rhi_user33',               c_float) )
        plist.append( ('rhi_user34',               c_float) )
        plist.append( ('rhi_user35',               c_float) )
        plist.append( ('rhi_user36',               c_float) )
        plist.append( ('rhi_user37',               c_float) )
        plist.append( ('rhi_user38',               c_float) )
        plist.append( ('rhi_user39',               c_float) )
        plist.append( ('rhi_user40',               c_float) )
        plist.append( ('rhi_user41',               c_float) )
        plist.append( ('rhi_user42',               c_float) )
        plist.append( ('rhi_user43',               c_float) )
        plist.append( ('rhi_user44',               c_float) )
        plist.append( ('rhi_user45',               c_float) )
        plist.append( ('rhi_user46',               c_float) )
        plist.append( ('rhi_user47',               c_float) )
        plist.append( ('rhi_user48',               c_float) )
        plist.append( ('pad_xx',                   c_char * 76) )
        plist.append( ('rhi_ctr_R',                c_float) )
        plist.append( ('rhi_ctr_A',                c_float) )
        plist.append( ('rhi_ctr_S',                c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhi_tlhc_R',               c_float) )
        plist.append( ('rhi_tlhc_A',               c_float) )
        plist.append( ('rhi_tlhc_S',               c_float) )
        plist.append( ('rhi_trhc_R',               c_float) )
        plist.append( ('rhi_trhc_A',               c_float) )
        plist.append( ('rhi_trhc_S',               c_float) )
        plist.append( ('rhi_brhc_R',               c_float) )
        plist.append( ('rhi_brhc_A',               c_float) )
        plist.append( ('rhi_brhc_S',               c_float) )
        plist.append( ('pad_xx',                   c_char * 300) )
        plist.append( ('rhi_tr',                   c_int) )
        plist.append( ('rhi_ti',                   c_int) )
        plist.append( ('rhi_te',                   c_int) )
        plist.append( ('pad_xx',                   c_char * 310) )
        plist.append( ('rhi_numecho',              c_short) )
        plist.append( ('pad_xx',                   c_char * 32) )
        plist.append( ('rhi_mr_flip',              c_short) )
        plist.append( ('pad_xx',                   c_char * 20) )
        plist.append( ('rhi_ctyp',                 c_short) )
        plist.append( ('pad_xx',                   c_char * 64) )
        plist.append( ('rhi_freq_dir',             c_short) )
        plist.append( ('pad_xx',                   c_char * 130) )
        plist.append( ('rhi_psdname',              c_char * 33) )
        plist.append( ('pad_xx',                   c_char * 84) )
        plist.append( ('rhi_cname',                c_char * 17) )
        plist.append( ('pad_xx',                   c_char * 51) )
        plist.append( ('rhi_image_uid',            c_char * 32) )

    elif version == 21:
        plist.append( ('rhr_rh_rdbm_rev',          c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhr_rh_scan_date',         c_char * 10) )
        plist.append( ('rhr_rh_scan_time',         c_char * 8) )
        plist.append( ('rhr_rh_logo',              c_char * 10) )
        plist.append( ('rhr_rh_file_contents',     c_short) )
        plist.append( ('pad_xx',                   c_char * 10) )
        plist.append( ('rhr_rh_data_collect_type', c_short) )
        plist.append( ('pad_xx',                   c_char * 6) )
        plist.append( ('rhr_rh_npasses',           c_short) )
        plist.append( ('pad_xx',                   c_char * 2) )
        plist.append( ('rhr_rh_nslices',           c_short) )
        plist.append( ('pad_xx',                   c_char * 10) )
        plist.append( ('rhr_rh_frame_size',        c_ushort) )
        plist.append( ('rhr_rh_point_size',        c_short) )
        plist.append( ('pad_xx',                   c_char * 116) )
        plist.append( ('rhr_rh_dab[0]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[0]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[1]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[1]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[2]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[2]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[3]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[3]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_user0',             c_float) )
        plist.append( ('rhr_rh_user1',             c_float) )
        plist.append( ('rhr_rh_user2',             c_float) )
        plist.append( ('rhr_rh_user3',             c_float) )
        plist.append( ('rhr_rh_user4',             c_float) )
        plist.append( ('rhr_rh_user5',             c_float) )
        plist.append( ('rhr_rh_user6',             c_float) )
        plist.append( ('rhr_rh_user7',             c_float) )
        plist.append( ('rhr_rh_user8',             c_float) )
        plist.append( ('rhr_rh_user9',             c_float) )
        plist.append( ('rhr_rh_user10',            c_float) )
        plist.append( ('rhr_rh_user11',            c_float) )
        plist.append( ('rhr_rh_user12',            c_float) )
        plist.append( ('rhr_rh_user13',            c_float) )
        plist.append( ('rhr_rh_user14',            c_float) )
        plist.append( ('rhr_rh_user15',            c_float) )
        plist.append( ('rhr_rh_user16',            c_float) )
        plist.append( ('rhr_rh_user17',            c_float) )
        plist.append( ('rhr_rh_user18',            c_float) )
        plist.append( ('rhr_rh_user19',            c_float) )
        plist.append( ('pad_xx',                   c_char * 72) )
        plist.append( ('rhr_spectral_width',       c_float) )
        plist.append( ('rhr_csi_dims',             c_short) )
        plist.append( ('rhr_xcsi',                 c_short) )
        plist.append( ('rhr_ycsi',                 c_short) )
        plist.append( ('rhr_zcsi',                 c_short) )
        plist.append( ('rhr_roilenx',              c_float) )
        plist.append( ('rhr_roileny',              c_float) )
        plist.append( ('rhr_roilenz',              c_float) )
        plist.append( ('pad_xx',                   c_char * 32) )
        plist.append( ('rhr_rh_ps_mps_freq',       c_uint) )
        plist.append( ('pad_xx',                   c_char * 560) )
        plist.append( ('rhr_rh_user_usage_tag',    c_uint) )
        plist.append( ('pad_xx',                   c_char * 8) )
        plist.append( ('rhr_rh_user20',            c_float) )
        plist.append( ('rhr_rh_user21',            c_float) )
        plist.append( ('rhr_rh_user22',            c_float) )
        plist.append( ('rhr_rh_user23',            c_float) )
        plist.append( ('rhr_rh_user24',            c_float) )
        plist.append( ('rhr_rh_user25',            c_float) )
        plist.append( ('rhr_rh_user26',            c_float) )
        plist.append( ('rhr_rh_user27',            c_float) )
        plist.append( ('rhr_rh_user28',            c_float) )
        plist.append( ('rhr_rh_user29',            c_float) )
        plist.append( ('rhr_rh_user30',            c_float) )
        plist.append( ('rhr_rh_user31',            c_float) )
        plist.append( ('rhr_rh_user32',            c_float) )
        plist.append( ('rhr_rh_user33',            c_float) )
        plist.append( ('rhr_rh_user34',            c_float) )
        plist.append( ('rhr_rh_user35',            c_float) )
        plist.append( ('rhr_rh_user36',            c_float) )
        plist.append( ('rhr_rh_user37',            c_float) )
        plist.append( ('rhr_rh_user38',            c_float) )
        plist.append( ('rhr_rh_user39',            c_float) )
        plist.append( ('rhr_rh_user40',            c_float) )
        plist.append( ('rhr_rh_user41',            c_float) )
        plist.append( ('rhr_rh_user42',            c_float) )
        plist.append( ('rhr_rh_user43',            c_float) )
        plist.append( ('rhr_rh_user44',            c_float) )
        plist.append( ('rhr_rh_user45',            c_float) )
        plist.append( ('rhr_rh_user46',            c_float) )
        plist.append( ('rhr_rh_user47',            c_float) )
        plist.append( ('rhr_rh_user48',            c_float) )
        plist.append( ('pad_xx',                   c_char * 352) )
        plist.append( ('rhr_rdb_hdr_off_data',     c_int) )
        plist.append( ('pad_xx',                   c_char * 188) )
        plist.append( ('rhr_rh_raw_pass_size',     c_longlong) )
        plist.append( ('pad_xx',                   c_char * 142272) )
        plist.append( ('rhe_magstrength',          c_int) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhe_ex_datetime',          c_int) )
        plist.append( ('pad_xx',                   c_char * 112) )
        plist.append( ('rhe_ex_no',                c_ushort) )
        plist.append( ('pad_xx',                   c_char * 22) )
        plist.append( ('rhe_patsex',               c_short) )
        plist.append( ('pad_xx',                   c_char * 335) )
        plist.append( ('rhe_refphy',               c_char * 65) )
        plist.append( ('pad_xx',                   c_char * 198) )
        plist.append( ('rhe_ex_sysid',             c_char * 9) )
        plist.append( ('pad_xx',                   c_char * 14) )
        plist.append( ('rhe_hospname',             c_char * 33) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhe_ex_verscre',           c_char * 2) )
        plist.append( ('pad_xx',                   c_char * 2) )
        plist.append( ('rhe_uniq_sys_id',          c_char * 16) )
        plist.append( ('pad_xx',                   c_char * 20) )
        plist.append( ('rhe_study_uid',            c_char * 32) )
        plist.append( ('pad_xx',                   c_char * 64) )
        plist.append( ('rhe_patname',              c_char * 65) )
        plist.append( ('rhe_patid',                c_char * 65) )
        plist.append( ('rhe_reqnum',               c_char * 17) )
        plist.append( ('rhe_dateofbirth',          c_char * 9) )
        plist.append( ('pad_xx',                   c_char * 928) )
        plist.append( ('rhs_position',             c_int) )
        plist.append( ('rhs_entry',                c_int) )
        plist.append( ('pad_xx',                   c_char * 194) )
        plist.append( ('rhs_se_no',                c_short) )
        plist.append( ('pad_xx',                   c_char * 138) )
        plist.append( ('rhs_se_desc',              c_char * 65) )
        plist.append( ('pad_xx',                   c_char * 18) )
        plist.append( ('rhs_anref',                c_char * 3) )
        plist.append( ('pad_xx',                   c_char * 27) )
        plist.append( ('rhs_series_uid',           c_char * 32) )
        plist.append( ('rhs_landmark_uid',         c_char * 32) )
        plist.append( ('pad_xx',                   c_char * 1705) )
        plist.append( ('rhi_dfov',                 c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhi_scanspacing',          c_float) )
        plist.append( ('rhi_loc',                  c_float) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhi_nex',                  c_float) )
        plist.append( ('pad_xx',                   c_char * 20) )
        plist.append( ('rhi_user0',                c_float) )
        plist.append( ('rhi_user1',                c_float) )
        plist.append( ('rhi_user2',                c_float) )
        plist.append( ('rhi_user3',                c_float) )
        plist.append( ('rhi_user4',                c_float) )
        plist.append( ('rhi_user5',                c_float) )
        plist.append( ('rhi_user6',                c_float) )
        plist.append( ('rhi_user7',                c_float) )
        plist.append( ('rhi_user8',                c_float) )
        plist.append( ('rhi_user9',                c_float) )
        plist.append( ('rhi_user10',               c_float) )
        plist.append( ('rhi_user11',               c_float) )
        plist.append( ('rhi_user12',               c_float) )
        plist.append( ('rhi_user13',               c_float) )
        plist.append( ('rhi_user14',               c_float) )
        plist.append( ('rhi_user15',               c_float) )
        plist.append( ('rhi_user16',               c_float) )
        plist.append( ('rhi_user17',               c_float) )
        plist.append( ('rhi_user18',               c_float) )
        plist.append( ('rhi_user19',               c_float) )
        plist.append( ('rhi_user20',               c_float) )
        plist.append( ('rhi_user21',               c_float) )
        plist.append( ('rhi_user22',               c_float) )
        plist.append( ('pad_xx',                   c_char * 8) )
        plist.append( ('rhi_user23',               c_float) )
        plist.append( ('rhi_user24',               c_float) )
        plist.append( ('pad_xx',                   c_char * 60) )
        plist.append( ('rhi_user25',               c_float) )
        plist.append( ('rhi_user26',               c_float) )
        plist.append( ('rhi_user27',               c_float) )
        plist.append( ('rhi_user28',               c_float) )
        plist.append( ('rhi_user29',               c_float) )
        plist.append( ('rhi_user30',               c_float) )
        plist.append( ('rhi_user31',               c_float) )
        plist.append( ('rhi_user32',               c_float) )
        plist.append( ('rhi_user33',               c_float) )
        plist.append( ('rhi_user34',               c_float) )
        plist.append( ('rhi_user35',               c_float) )
        plist.append( ('rhi_user36',               c_float) )
        plist.append( ('rhi_user37',               c_float) )
        plist.append( ('rhi_user38',               c_float) )
        plist.append( ('rhi_user39',               c_float) )
        plist.append( ('rhi_user40',               c_float) )
        plist.append( ('rhi_user41',               c_float) )
        plist.append( ('rhi_user42',               c_float) )
        plist.append( ('rhi_user43',               c_float) )
        plist.append( ('rhi_user44',               c_float) )
        plist.append( ('rhi_user45',               c_float) )
        plist.append( ('rhi_user46',               c_float) )
        plist.append( ('rhi_user47',               c_float) )
        plist.append( ('rhi_user48',               c_float) )
        plist.append( ('pad_xx',                   c_char * 76) )
        plist.append( ('rhi_ctr_R',                c_float) )
        plist.append( ('rhi_ctr_A',                c_float) )
        plist.append( ('rhi_ctr_S',                c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhi_tlhc_R',               c_float) )
        plist.append( ('rhi_tlhc_A',               c_float) )
        plist.append( ('rhi_tlhc_S',               c_float) )
        plist.append( ('rhi_trhc_R',               c_float) )
        plist.append( ('rhi_trhc_A',               c_float) )
        plist.append( ('rhi_trhc_S',               c_float) )
        plist.append( ('rhi_brhc_R',               c_float) )
        plist.append( ('rhi_brhc_A',               c_float) )
        plist.append( ('rhi_brhc_S',               c_float) )
        plist.append( ('pad_xx',                   c_char * 300) )
        plist.append( ('rhi_tr',                   c_int) )
        plist.append( ('rhi_ti',                   c_int) )
        plist.append( ('rhi_te',                   c_int) )
        plist.append( ('pad_xx',                   c_char * 310) )
        plist.append( ('rhi_numecho',              c_short) )
        plist.append( ('pad_xx',                   c_char * 32) )
        plist.append( ('rhi_mr_flip',              c_short) )
        plist.append( ('pad_xx',                   c_char * 20) )
        plist.append( ('rhi_ctyp',                 c_short) )
        plist.append( ('pad_xx',                   c_char * 64) )
        plist.append( ('rhi_freq_dir',             c_short) )
        plist.append( ('pad_xx',                   c_char * 130) )
        plist.append( ('rhi_psdname',              c_char * 33) )
        plist.append( ('pad_xx',                   c_char * 84) )
        plist.append( ('rhi_cname',                c_char * 17) )
        plist.append( ('pad_xx',                   c_char * 51) )
        plist.append( ('rhi_image_uid',            c_char * 32) )

    elif version == 23:
        plist.append( ('rhr_rh_rdbm_rev',          c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhr_rh_scan_date',         c_char * 10) )
        plist.append( ('rhr_rh_scan_time',         c_char * 8) )
        plist.append( ('rhr_rh_logo',              c_char * 10) )
        plist.append( ('rhr_rh_file_contents',     c_short) )
        plist.append( ('pad_xx',                   c_char * 10) )
        plist.append( ('rhr_rh_data_collect_type', c_short) )
        plist.append( ('pad_xx',                   c_char * 6) )
        plist.append( ('rhr_rh_npasses',           c_short) )
        plist.append( ('pad_xx',                   c_char * 2) )
        plist.append( ('rhr_rh_nslices',           c_short) )
        plist.append( ('pad_xx',                   c_char * 10) )
        plist.append( ('rhr_rh_frame_size',        c_ushort) )
        plist.append( ('rhr_rh_point_size',        c_short) )
        plist.append( ('pad_xx',                   c_char * 116) )
        plist.append( ('rhr_rh_dab[0]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[0]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[1]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[1]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[2]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[2]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[3]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[3]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_user0',             c_float) )
        plist.append( ('rhr_rh_user1',             c_float) )
        plist.append( ('rhr_rh_user2',             c_float) )
        plist.append( ('rhr_rh_user3',             c_float) )
        plist.append( ('rhr_rh_user4',             c_float) )
        plist.append( ('rhr_rh_user5',             c_float) )
        plist.append( ('rhr_rh_user6',             c_float) )
        plist.append( ('rhr_rh_user7',             c_float) )
        plist.append( ('rhr_rh_user8',             c_float) )
        plist.append( ('rhr_rh_user9',             c_float) )
        plist.append( ('rhr_rh_user10',            c_float) )
        plist.append( ('rhr_rh_user11',            c_float) )
        plist.append( ('rhr_rh_user12',            c_float) )
        plist.append( ('rhr_rh_user13',            c_float) )
        plist.append( ('rhr_rh_user14',            c_float) )
        plist.append( ('rhr_rh_user15',            c_float) )
        plist.append( ('rhr_rh_user16',            c_float) )
        plist.append( ('rhr_rh_user17',            c_float) )
        plist.append( ('rhr_rh_user18',            c_float) )
        plist.append( ('rhr_rh_user19',            c_float) )
        plist.append( ('pad_xx',                   c_char * 72) )
        plist.append( ('rhr_spectral_width',       c_float) )
        plist.append( ('rhr_csi_dims',             c_short) )
        plist.append( ('rhr_xcsi',                 c_short) )
        plist.append( ('rhr_ycsi',                 c_short) )
        plist.append( ('rhr_zcsi',                 c_short) )
        plist.append( ('rhr_roilenx',              c_float) )
        plist.append( ('rhr_roileny',              c_float) )
        plist.append( ('rhr_roilenz',              c_float) )
        plist.append( ('pad_xx',                   c_char * 32) )
        plist.append( ('rhr_rh_ps_mps_freq',       c_uint) )
        plist.append( ('pad_xx',                   c_char * 560) )
        plist.append( ('rhr_rh_user_usage_tag',    c_uint) )
        plist.append( ('pad_xx',                   c_char * 8) )
        plist.append( ('rhr_rh_user20',            c_float) )
        plist.append( ('rhr_rh_user21',            c_float) )
        plist.append( ('rhr_rh_user22',            c_float) )
        plist.append( ('rhr_rh_user23',            c_float) )
        plist.append( ('rhr_rh_user24',            c_float) )
        plist.append( ('rhr_rh_user25',            c_float) )
        plist.append( ('rhr_rh_user26',            c_float) )
        plist.append( ('rhr_rh_user27',            c_float) )
        plist.append( ('rhr_rh_user28',            c_float) )
        plist.append( ('rhr_rh_user29',            c_float) )
        plist.append( ('rhr_rh_user30',            c_float) )
        plist.append( ('rhr_rh_user31',            c_float) )
        plist.append( ('rhr_rh_user32',            c_float) )
        plist.append( ('rhr_rh_user33',            c_float) )
        plist.append( ('rhr_rh_user34',            c_float) )
        plist.append( ('rhr_rh_user35',            c_float) )
        plist.append( ('rhr_rh_user36',            c_float) )
        plist.append( ('rhr_rh_user37',            c_float) )
        plist.append( ('rhr_rh_user38',            c_float) )
        plist.append( ('rhr_rh_user39',            c_float) )
        plist.append( ('rhr_rh_user40',            c_float) )
        plist.append( ('rhr_rh_user41',            c_float) )
        plist.append( ('rhr_rh_user42',            c_float) )
        plist.append( ('rhr_rh_user43',            c_float) )
        plist.append( ('rhr_rh_user44',            c_float) )
        plist.append( ('rhr_rh_user45',            c_float) )
        plist.append( ('rhr_rh_user46',            c_float) )
        plist.append( ('rhr_rh_user47',            c_float) )
        plist.append( ('rhr_rh_user48',            c_float) )
        plist.append( ('pad_xx',                   c_char * 352) )
        plist.append( ('rhr_rdb_hdr_off_data',     c_int) )
        plist.append( ('pad_xx',                   c_char * 188) )
        plist.append( ('rhr_rh_raw_pass_size',     c_longlong) )
        plist.append( ('pad_xx',                   c_char * 141724) )
        plist.append( ('rhe_magstrength',          c_int) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhe_ex_datetime',          c_int) )
        plist.append( ('pad_xx',                   c_char * 112) )
        plist.append( ('rhe_ex_no',                c_ushort) )
        plist.append( ('pad_xx',                   c_char * 22) )
        plist.append( ('rhe_patsex',               c_short) )
        plist.append( ('pad_xx',                   c_char * 335) )
        plist.append( ('rhe_refphy',               c_char * 65) )
        plist.append( ('pad_xx',                   c_char * 198) )
        plist.append( ('rhe_ex_sysid',             c_char * 9) )
        plist.append( ('pad_xx',                   c_char * 22) )
        plist.append( ('rhe_hospname',             c_char * 33) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhe_ex_verscre',           c_char * 2) )
        plist.append( ('pad_xx',                   c_char * 2) )
        plist.append( ('rhe_uniq_sys_id',          c_char * 16) )
        plist.append( ('pad_xx',                   c_char * 20) )
        plist.append( ('rhe_study_uid',            c_char * 32) )
        plist.append( ('pad_xx',                   c_char * 64) )
        plist.append( ('rhe_patname',              c_char * 65) )
        plist.append( ('rhe_patid',                c_char * 65) )
        plist.append( ('rhe_reqnum',               c_char * 17) )
        plist.append( ('rhe_dateofbirth',          c_char * 9) )
        plist.append( ('pad_xx',                   c_char * 920) )
        plist.append( ('rhs_position',             c_int) )
        plist.append( ('rhs_entry',                c_int) )
        plist.append( ('pad_xx',                   c_char * 194) )
        plist.append( ('rhs_se_no',                c_short) )
        plist.append( ('pad_xx',                   c_char * 138) )
        plist.append( ('rhs_se_desc',              c_char * 65) )
        plist.append( ('pad_xx',                   c_char * 18) )
        plist.append( ('rhs_anref',                c_char * 3) )
        plist.append( ('pad_xx',                   c_char * 27) )
        plist.append( ('rhs_series_uid',           c_char * 32) )
        plist.append( ('rhs_landmark_uid',         c_char * 32) )
        plist.append( ('pad_xx',                   c_char * 1705) )
        plist.append( ('rhi_dfov',                 c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhi_scanspacing',          c_float) )
        plist.append( ('rhi_loc',                  c_float) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhi_nex',                  c_float) )
        plist.append( ('pad_xx',                   c_char * 20) )
        plist.append( ('rhi_user0',                c_float) )
        plist.append( ('rhi_user1',                c_float) )
        plist.append( ('rhi_user2',                c_float) )
        plist.append( ('rhi_user3',                c_float) )
        plist.append( ('rhi_user4',                c_float) )
        plist.append( ('rhi_user5',                c_float) )
        plist.append( ('rhi_user6',                c_float) )
        plist.append( ('rhi_user7',                c_float) )
        plist.append( ('rhi_user8',                c_float) )
        plist.append( ('rhi_user9',                c_float) )
        plist.append( ('rhi_user10',               c_float) )
        plist.append( ('rhi_user11',               c_float) )
        plist.append( ('rhi_user12',               c_float) )
        plist.append( ('rhi_user13',               c_float) )
        plist.append( ('rhi_user14',               c_float) )
        plist.append( ('rhi_user15',               c_float) )
        plist.append( ('rhi_user16',               c_float) )
        plist.append( ('rhi_user17',               c_float) )
        plist.append( ('rhi_user18',               c_float) )
        plist.append( ('rhi_user19',               c_float) )
        plist.append( ('rhi_user20',               c_float) )
        plist.append( ('rhi_user21',               c_float) )
        plist.append( ('rhi_user22',               c_float) )
        plist.append( ('pad_xx',                   c_char * 8) )
        plist.append( ('rhi_user23',               c_float) )
        plist.append( ('rhi_user24',               c_float) )
        plist.append( ('pad_xx',                   c_char * 60) )
        plist.append( ('rhi_user25',               c_float) )
        plist.append( ('rhi_user26',               c_float) )
        plist.append( ('rhi_user27',               c_float) )
        plist.append( ('rhi_user28',               c_float) )
        plist.append( ('rhi_user29',               c_float) )
        plist.append( ('rhi_user30',               c_float) )
        plist.append( ('rhi_user31',               c_float) )
        plist.append( ('rhi_user32',               c_float) )
        plist.append( ('rhi_user33',               c_float) )
        plist.append( ('rhi_user34',               c_float) )
        plist.append( ('rhi_user35',               c_float) )
        plist.append( ('rhi_user36',               c_float) )
        plist.append( ('rhi_user37',               c_float) )
        plist.append( ('rhi_user38',               c_float) )
        plist.append( ('rhi_user39',               c_float) )
        plist.append( ('rhi_user40',               c_float) )
        plist.append( ('rhi_user41',               c_float) )
        plist.append( ('rhi_user42',               c_float) )
        plist.append( ('rhi_user43',               c_float) )
        plist.append( ('rhi_user44',               c_float) )
        plist.append( ('rhi_user45',               c_float) )
        plist.append( ('rhi_user46',               c_float) )
        plist.append( ('rhi_user47',               c_float) )
        plist.append( ('rhi_user48',               c_float) )
        plist.append( ('pad_xx',                   c_char * 76) )
        plist.append( ('rhi_ctr_R',                c_float) )
        plist.append( ('rhi_ctr_A',                c_float) )
        plist.append( ('rhi_ctr_S',                c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhi_tlhc_R',               c_float) )
        plist.append( ('rhi_tlhc_A',               c_float) )
        plist.append( ('rhi_tlhc_S',               c_float) )
        plist.append( ('rhi_trhc_R',               c_float) )
        plist.append( ('rhi_trhc_A',               c_float) )
        plist.append( ('rhi_trhc_S',               c_float) )
        plist.append( ('rhi_brhc_R',               c_float) )
        plist.append( ('rhi_brhc_A',               c_float) )
        plist.append( ('rhi_brhc_S',               c_float) )
        plist.append( ('pad_xx',                   c_char * 300) )
        plist.append( ('rhi_tr',                   c_int) )
        plist.append( ('rhi_ti',                   c_int) )
        plist.append( ('rhi_te',                   c_int) )
        plist.append( ('pad_xx',                   c_char * 310) )
        plist.append( ('rhi_numecho',              c_short) )
        plist.append( ('pad_xx',                   c_char * 32) )
        plist.append( ('rhi_mr_flip',              c_short) )
        plist.append( ('pad_xx',                   c_char * 20) )
        plist.append( ('rhi_ctyp',                 c_short) )
        plist.append( ('pad_xx',                   c_char * 64) )
        plist.append( ('rhi_freq_dir',             c_short) )
        plist.append( ('pad_xx',                   c_char * 130) )
        plist.append( ('rhi_psdname',              c_char * 33) )
        plist.append( ('pad_xx',                   c_char * 84) )
        plist.append( ('rhi_cname',                c_char * 17) )
        plist.append( ('pad_xx',                   c_char * 51) )
        plist.append( ('rhi_image_uid',            c_char * 32) )

    elif version == 24:  # 25 and 25.1 same pfile
        plist.append( ('rhr_rh_rdbm_rev',          c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhr_rh_scan_date',         c_char * 10) )
        plist.append( ('rhr_rh_scan_time',         c_char * 8) )
        plist.append( ('rhr_rh_logo',              c_char * 10) )
        plist.append( ('rhr_rh_file_contents',     c_short) )
        plist.append( ('pad_xx',                   c_char * 10) )
        plist.append( ('rhr_rh_data_collect_type', c_short) )
        plist.append( ('pad_xx',                   c_char * 6) )
        plist.append( ('rhr_rh_npasses',           c_short) )
        plist.append( ('pad_xx',                   c_char * 2) )
        plist.append( ('rhr_rh_nslices',           c_short) )
        plist.append( ('pad_xx',                   c_char * 10) )
        plist.append( ('rhr_rh_frame_size',        c_ushort) )
        plist.append( ('rhr_rh_point_size',        c_short) )
        plist.append( ('pad_xx',                   c_char * 116) )
        plist.append( ('rhr_rh_dab[0]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[0]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[1]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[1]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[2]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[2]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[3]_start_rcv',  c_short) )
        plist.append( ('rhr_rh_dab[3]_stop_rcv',   c_short) )
        plist.append( ('rhr_rh_user0',             c_float) )
        plist.append( ('rhr_rh_user1',             c_float) )
        plist.append( ('rhr_rh_user2',             c_float) )
        plist.append( ('rhr_rh_user3',             c_float) )
        plist.append( ('rhr_rh_user4',             c_float) )
        plist.append( ('rhr_rh_user5',             c_float) )
        plist.append( ('rhr_rh_user6',             c_float) )
        plist.append( ('rhr_rh_user7',             c_float) )
        plist.append( ('rhr_rh_user8',             c_float) )
        plist.append( ('rhr_rh_user9',             c_float) )
        plist.append( ('rhr_rh_user10',            c_float) )
        plist.append( ('rhr_rh_user11',            c_float) )
        plist.append( ('rhr_rh_user12',            c_float) )
        plist.append( ('rhr_rh_user13',            c_float) )
        plist.append( ('rhr_rh_user14',            c_float) )
        plist.append( ('rhr_rh_user15',            c_float) )
        plist.append( ('rhr_rh_user16',            c_float) )
        plist.append( ('rhr_rh_user17',            c_float) )
        plist.append( ('rhr_rh_user18',            c_float) )
        plist.append( ('rhr_rh_user19',            c_float) )
        plist.append( ('pad_xx',                   c_char * 72) )
        plist.append( ('rhr_spectral_width',       c_float) )
        plist.append( ('rhr_csi_dims',             c_short) )
        plist.append( ('rhr_xcsi',                 c_short) )
        plist.append( ('rhr_ycsi',                 c_short) )
        plist.append( ('rhr_zcsi',                 c_short) )
        plist.append( ('rhr_roilenx',              c_float) )
        plist.append( ('rhr_roileny',              c_float) )
        plist.append( ('rhr_roilenz',              c_float) )
        plist.append( ('pad_xx',                   c_char * 32) )
        plist.append( ('rhr_rh_ps_mps_freq',       c_uint) )
        plist.append( ('pad_xx',                   c_char * 560) )
        plist.append( ('rhr_rh_user_usage_tag',    c_uint) )
        plist.append( ('pad_xx',                   c_char * 8) )
        plist.append( ('rhr_rh_user20',            c_float) )
        plist.append( ('rhr_rh_user21',            c_float) )
        plist.append( ('rhr_rh_user22',            c_float) )
        plist.append( ('rhr_rh_user23',            c_float) )
        plist.append( ('rhr_rh_user24',            c_float) )
        plist.append( ('rhr_rh_user25',            c_float) )
        plist.append( ('rhr_rh_user26',            c_float) )
        plist.append( ('rhr_rh_user27',            c_float) )
        plist.append( ('rhr_rh_user28',            c_float) )
        plist.append( ('rhr_rh_user29',            c_float) )
        plist.append( ('rhr_rh_user30',            c_float) )
        plist.append( ('rhr_rh_user31',            c_float) )
        plist.append( ('rhr_rh_user32',            c_float) )
        plist.append( ('rhr_rh_user33',            c_float) )
        plist.append( ('rhr_rh_user34',            c_float) )
        plist.append( ('rhr_rh_user35',            c_float) )
        plist.append( ('rhr_rh_user36',            c_float) )
        plist.append( ('rhr_rh_user37',            c_float) )
        plist.append( ('rhr_rh_user38',            c_float) )
        plist.append( ('rhr_rh_user39',            c_float) )
        plist.append( ('rhr_rh_user40',            c_float) )
        plist.append( ('rhr_rh_user41',            c_float) )
        plist.append( ('rhr_rh_user42',            c_float) )
        plist.append( ('rhr_rh_user43',            c_float) )
        plist.append( ('rhr_rh_user44',            c_float) )
        plist.append( ('rhr_rh_user45',            c_float) )
        plist.append( ('rhr_rh_user46',            c_float) )
        plist.append( ('rhr_rh_user47',            c_float) )
        plist.append( ('rhr_rh_user48',            c_float) )
        plist.append( ('pad_xx',                   c_char * 352) )
        plist.append( ('rhr_rdb_hdr_off_data',     c_int) )
        plist.append( ('pad_xx',                   c_char * 188) )
        plist.append( ('rhr_rh_raw_pass_size',     c_longlong) )
        plist.append( ('pad_xx',                   c_char * 141724) )
        plist.append( ('rhe_magstrength',          c_int) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhe_ex_datetime',          c_int) )
        plist.append( ('pad_xx',                   c_char * 112) )
        plist.append( ('rhe_ex_no',                c_ushort) )
        plist.append( ('pad_xx',                   c_char * 22) )
        plist.append( ('rhe_patsex',               c_short) )
        plist.append( ('pad_xx',                   c_char * 335) )
        plist.append( ('rhe_refphy',               c_char * 65) )
        plist.append( ('pad_xx',                   c_char * 198) )
        plist.append( ('rhe_ex_sysid',             c_char * 9) )
        plist.append( ('pad_xx',                   c_char * 22) )
        plist.append( ('rhe_hospname',             c_char * 33) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhe_ex_verscre',           c_char * 2) )
        plist.append( ('pad_xx',                   c_char * 2) )
        plist.append( ('rhe_uniq_sys_id',          c_char * 16) )
        plist.append( ('pad_xx',                   c_char * 20) )
        plist.append( ('rhe_study_uid',            c_char * 32) )
        plist.append( ('pad_xx',                   c_char * 64) )
        plist.append( ('rhe_patname',              c_char * 65) )
        plist.append( ('rhe_patid',                c_char * 65) )
        plist.append( ('rhe_reqnum',               c_char * 17) )
        plist.append( ('rhe_dateofbirth',          c_char * 9) )
        plist.append( ('pad_xx',                   c_char * 920) )
        plist.append( ('rhs_position',             c_int) )
        plist.append( ('rhs_entry',                c_int) )
        plist.append( ('pad_xx',                   c_char * 194) )
        plist.append( ('rhs_se_no',                c_short) )
        plist.append( ('pad_xx',                   c_char * 138) )
        plist.append( ('rhs_se_desc',              c_char * 65) )
        plist.append( ('pad_xx',                   c_char * 18) )
        plist.append( ('rhs_anref',                c_char * 3) )
        plist.append( ('pad_xx',                   c_char * 27) )
        plist.append( ('rhs_series_uid',           c_char * 32) )
        plist.append( ('rhs_landmark_uid',         c_char * 32) )
        plist.append( ('pad_xx',                   c_char * 1705) )
        plist.append( ('rhi_dfov',                 c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhi_scanspacing',          c_float) )
        plist.append( ('rhi_loc',                  c_float) )
        plist.append( ('pad_xx',                   c_char * 4) )
        plist.append( ('rhi_nex',                  c_float) )
        plist.append( ('pad_xx',                   c_char * 20) )
        plist.append( ('rhi_user0',                c_float) )
        plist.append( ('rhi_user1',                c_float) )
        plist.append( ('rhi_user2',                c_float) )
        plist.append( ('rhi_user3',                c_float) )
        plist.append( ('rhi_user4',                c_float) )
        plist.append( ('rhi_user5',                c_float) )
        plist.append( ('rhi_user6',                c_float) )
        plist.append( ('rhi_user7',                c_float) )
        plist.append( ('rhi_user8',                c_float) )
        plist.append( ('rhi_user9',                c_float) )
        plist.append( ('rhi_user10',               c_float) )
        plist.append( ('rhi_user11',               c_float) )
        plist.append( ('rhi_user12',               c_float) )
        plist.append( ('rhi_user13',               c_float) )
        plist.append( ('rhi_user14',               c_float) )
        plist.append( ('rhi_user15',               c_float) )
        plist.append( ('rhi_user16',               c_float) )
        plist.append( ('rhi_user17',               c_float) )
        plist.append( ('rhi_user18',               c_float) )
        plist.append( ('rhi_user19',               c_float) )
        plist.append( ('rhi_user20',               c_float) )
        plist.append( ('rhi_user21',               c_float) )
        plist.append( ('rhi_user22',               c_float) )
        plist.append( ('pad_xx',                   c_char * 8) )
        plist.append( ('rhi_user23',               c_float) )
        plist.append( ('rhi_user24',               c_float) )
        plist.append( ('pad_xx',                   c_char * 60) )
        plist.append( ('rhi_user25',               c_float) )
        plist.append( ('rhi_user26',               c_float) )
        plist.append( ('rhi_user27',               c_float) )
        plist.append( ('rhi_user28',               c_float) )
        plist.append( ('rhi_user29',               c_float) )
        plist.append( ('rhi_user30',               c_float) )
        plist.append( ('rhi_user31',               c_float) )
        plist.append( ('rhi_user32',               c_float) )
        plist.append( ('rhi_user33',               c_float) )
        plist.append( ('rhi_user34',               c_float) )
        plist.append( ('rhi_user35',               c_float) )
        plist.append( ('rhi_user36',               c_float) )
        plist.append( ('rhi_user37',               c_float) )
        plist.append( ('rhi_user38',               c_float) )
        plist.append( ('rhi_user39',               c_float) )
        plist.append( ('rhi_user40',               c_float) )
        plist.append( ('rhi_user41',               c_float) )
        plist.append( ('rhi_user42',               c_float) )
        plist.append( ('rhi_user43',               c_float) )
        plist.append( ('rhi_user44',               c_float) )
        plist.append( ('rhi_user45',               c_float) )
        plist.append( ('rhi_user46',               c_float) )
        plist.append( ('rhi_user47',               c_float) )
        plist.append( ('rhi_user48',               c_float) )
        plist.append( ('pad_xx',                   c_char * 76) )
        plist.append( ('rhi_ctr_R',                c_float) )
        plist.append( ('rhi_ctr_A',                c_float) )
        plist.append( ('rhi_ctr_S',                c_float) )
        plist.append( ('pad_xx',                   c_char * 12) )
        plist.append( ('rhi_tlhc_R',               c_float) )
        plist.append( ('rhi_tlhc_A',               c_float) )
        plist.append( ('rhi_tlhc_S',               c_float) )
        plist.append( ('rhi_trhc_R',               c_float) )
        plist.append( ('rhi_trhc_A',               c_float) )
        plist.append( ('rhi_trhc_S',               c_float) )
        plist.append( ('rhi_brhc_R',               c_float) )
        plist.append( ('rhi_brhc_A',               c_float) )
        plist.append( ('rhi_brhc_S',               c_float) )
        plist.append( ('pad_xx',                   c_char * 300) )
        plist.append( ('rhi_tr',                   c_int) )
        plist.append( ('rhi_ti',                   c_int) )
        plist.append( ('rhi_te',                   c_int) )
        plist.append( ('pad_xx',                   c_char * 310) )
        plist.append( ('rhi_numecho',              c_short) )
        plist.append( ('pad_xx',                   c_char * 32) )
        plist.append( ('rhi_mr_flip',              c_short) )
        plist.append( ('pad_xx',                   c_char * 20) )
        plist.append( ('rhi_ctyp',                 c_short) )
        plist.append( ('pad_xx',                   c_char * 64) )
        plist.append( ('rhi_freq_dir',             c_short) )
        plist.append( ('pad_xx',                   c_char * 130) )
        plist.append( ('rhi_psdname',              c_char * 33) )
        plist.append( ('pad_xx',                   c_char * 84) )
        plist.append( ('rhi_cname',                c_char * 17) )
        plist.append( ('pad_xx',                   c_char * 51) )
        plist.append( ('rhi_image_uid',            c_char * 32) )

    elif version == 26:
        plist.append( ('rhr_rh_rdbm_rev',           c_float) )
        plist.append( ('rhr_rdb_hdr_off_data',      c_int) )
        plist.append( ('pad_xx',                    c_char * 84) )
        plist.append( ('rhr_rh_scan_date',          c_char * 10) )
        plist.append( ('rhr_rh_scan_time',          c_char * 8) )
        plist.append( ('rhr_rh_logo',               c_char * 10) )
        plist.append( ('rhr_rh_file_contents',      c_short) )
        plist.append( ('pad_xx',                    c_char * 10) )
        plist.append( ('rhr_rh_data_collect_type',  c_short) )
        plist.append( ('pad_xx',                    c_char * 6) )
        plist.append( ('rhr_rh_npasses',            c_short) )
        plist.append( ('pad_xx',                    c_char * 2) )
        plist.append( ('rhr_rh_nslices',            c_short) )
        plist.append( ('pad_xx',                    c_char * 10) )
        plist.append( ('rhr_rh_frame_size',         c_short) )
        plist.append( ('rhr_rh_point_size',         c_short) )
        plist.append( ('pad_xx',                    c_char * 104) )
        plist.append( ('rhr_rh_dab[0]_start_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[0]_stop_rcv',    c_short) )
        plist.append( ('rhr_rh_dab[1]_start_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[1]_stop_rcv',    c_short) )
        plist.append( ('rhr_rh_dab[2]_start_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[2]_stop_rcv',    c_short) )
        plist.append( ('rhr_rh_dab[3]_start_rcv',   c_short) )
        plist.append( ('rhr_rh_dab[3]_stop_rcv',    c_short) )
        plist.append( ('rhr_rh_user0',              c_float) )
        plist.append( ('rhr_rh_user1',              c_float) )
        plist.append( ('rhr_rh_user2',              c_float) )
        plist.append( ('rhr_rh_user3',              c_float) )
        plist.append( ('rhr_rh_user4',              c_float) )
        plist.append( ('rhr_rh_user5',              c_float) )
        plist.append( ('rhr_rh_user6',              c_float) )
        plist.append( ('rhr_rh_user7',              c_float) )
        plist.append( ('rhr_rh_user8',              c_float) )
        plist.append( ('rhr_rh_user9',              c_float) )
        plist.append( ('rhr_rh_user10',             c_float) )
        plist.append( ('rhr_rh_user11',             c_float) )
        plist.append( ('rhr_rh_user12',             c_float) )
        plist.append( ('rhr_rh_user13',             c_float) )
        plist.append( ('rhr_rh_user14',             c_float) )
        plist.append( ('rhr_rh_user15',             c_float) )
        plist.append( ('rhr_rh_user16',             c_float) )
        plist.append( ('rhr_rh_user17',             c_float) )
        plist.append( ('rhr_rh_user18',             c_float) )
        plist.append( ('rhr_rh_user19',             c_float) )
        plist.append( ('pad_xx',                    c_char * 72) )
        plist.append( ('rhr_spectral_width',        c_float) )
        plist.append( ('rhr_csi_dims',              c_short) )
        plist.append( ('rhr_xcsi',                  c_short) )
        plist.append( ('rhr_ycsi',                  c_short) )
        plist.append( ('rhr_zcsi',                  c_short) )
        plist.append( ('rhr_roilenx',               c_float) )
        plist.append( ('rhr_roileny',               c_float) )
        plist.append( ('rhr_roilenz',               c_float) )
        plist.append( ('pad_xx',                    c_char * 32) )
        plist.append( ('rhr_rh_ps_mps_freq',        c_uint) )
        plist.append( ('pad_xx',                    c_char * 432) )
        plist.append( ('rhr_rh_user_usage_tag',     c_uint) )
        plist.append( ('pad_xx',                    c_char * 8) )
        plist.append( ('rhr_rh_user20',             c_float) )
        plist.append( ('rhr_rh_user21',             c_float) )
        plist.append( ('rhr_rh_user22',             c_float) )
        plist.append( ('rhr_rh_user23',             c_float) )
        plist.append( ('rhr_rh_user24',             c_float) )
        plist.append( ('rhr_rh_user25',             c_float) )
        plist.append( ('rhr_rh_user26',             c_float) )
        plist.append( ('rhr_rh_user27',             c_float) )
        plist.append( ('rhr_rh_user28',             c_float) )
        plist.append( ('rhr_rh_user29',             c_float) )
        plist.append( ('rhr_rh_user30',             c_float) )
        plist.append( ('rhr_rh_user31',             c_float) )
        plist.append( ('rhr_rh_user32',             c_float) )
        plist.append( ('rhr_rh_user33',             c_float) )
        plist.append( ('rhr_rh_user34',             c_float) )
        plist.append( ('rhr_rh_user35',             c_float) )
        plist.append( ('rhr_rh_user36',             c_float) )
        plist.append( ('rhr_rh_user37',             c_float) )
        plist.append( ('rhr_rh_user38',             c_float) )
        plist.append( ('rhr_rh_user39',             c_float) )
        plist.append( ('rhr_rh_user40',             c_float) )
        plist.append( ('rhr_rh_user41',             c_float) )
        plist.append( ('rhr_rh_user42',             c_float) )
        plist.append( ('rhr_rh_user43',             c_float) )
        plist.append( ('rhr_rh_user44',             c_float) )
        plist.append( ('rhr_rh_user45',             c_float) )
        plist.append( ('rhr_rh_user46',             c_float) )
        plist.append( ('rhr_rh_user47',             c_float) )
        plist.append( ('rhr_rh_user48',             c_float) )
        plist.append( ('pad_xx',                    c_char * 488) )
        plist.append( ('rhr_rh_raw_pass_size',      c_longlong) )
        plist.append( ('pad_xx',                    c_char * 192684) )
        plist.append( ('rhe_magstrength',           c_int) )
        plist.append( ('pad_xx',                    c_char * 4) )
        plist.append( ('rhe_ex_datetime',           c_int) )
        plist.append( ('pad_xx',                    c_char * 112) )
        plist.append( ('rhe_ex_no',                 c_short) )
        plist.append( ('pad_xx',                    c_char * 22) )
        plist.append( ('rhe_patsex',                c_short) )
        plist.append( ('pad_xx',                    c_char * 335) )
        plist.append( ('rhe_refphy',                c_char * 65) )
        plist.append( ('pad_xx',                    c_char * 198) )
        plist.append( ('rhe_ex_sysid',              c_char * 9) )
        plist.append( ('pad_xx',                    c_char * 22) )
        plist.append( ('rhe_hospname',              c_char * 33) )
        plist.append( ('pad_xx',                    c_char * 4) )
        plist.append( ('rhe_ex_verscre',            c_char * 2) )
        plist.append( ('pad_xx',                    c_char * 2) )
        plist.append( ('rhe_uniq_sys_id',           c_char * 16) )
        plist.append( ('pad_xx',                    c_char * 20) )
        plist.append( ('rhe_study_uid',             c_char * 32) )
        plist.append( ('pad_xx',                    c_char * 64) )
        plist.append( ('rhe_patname',               c_char * 65) )
        plist.append( ('rhe_patid',                 c_char * 65) )
        plist.append( ('rhe_reqnum',                c_char * 17) )
        plist.append( ('rhe_dateofbirth',           c_char * 9) )
        plist.append( ('pad_xx',                    c_char * 920) )
        plist.append( ('rhs_position',              c_int) )
        plist.append( ('rhs_entry',                 c_int) )
        plist.append( ('pad_xx',                    c_char * 88) )
        plist.append( ('rhs_se_no',                 c_int) )
        plist.append( ('pad_xx',                    c_char * 242) )
        plist.append( ('rhs_se_desc',               c_char * 65) )
        plist.append( ('pad_xx',                    c_char * 18) )
        plist.append( ('rhs_anref',                 c_char * 3) )
        plist.append( ('pad_xx',                    c_char * 77) )
        plist.append( ('rhs_series_uid',            c_char * 32) )
        plist.append( ('rhs_landmark_uid',          c_char * 32) )
        plist.append( ('pad_xx',                    c_char * 1655) )
        plist.append( ('rhi_dfov',                  c_float) )
        plist.append( ('pad_xx',                    c_char * 12) )
        plist.append( ('rhi_scanspacing',           c_float) )
        plist.append( ('rhi_loc',                   c_float) )
        plist.append( ('pad_xx',                    c_char * 4) )
        plist.append( ('rhi_nex',                   c_float) )
        plist.append( ('pad_xx',                    c_char * 20) )
        plist.append( ('rhi_user0',                 c_float) )
        plist.append( ('rhi_user1',                 c_float) )
        plist.append( ('rhi_user2',                 c_float) )
        plist.append( ('rhi_user3',                 c_float) )
        plist.append( ('rhi_user4',                 c_float) )
        plist.append( ('rhi_user5',                 c_float) )
        plist.append( ('rhi_user6',                 c_float) )
        plist.append( ('rhi_user7',                 c_float) )
        plist.append( ('rhi_user8',                 c_float) )
        plist.append( ('rhi_user9',                 c_float) )
        plist.append( ('rhi_user10',                c_float) )
        plist.append( ('rhi_user11',                c_float) )
        plist.append( ('rhi_user12',                c_float) )
        plist.append( ('rhi_user13',                c_float) )
        plist.append( ('rhi_user14',                c_float) )
        plist.append( ('rhi_user15',                c_float) )
        plist.append( ('rhi_user16',                c_float) )
        plist.append( ('rhi_user17',                c_float) )
        plist.append( ('rhi_user18',                c_float) )
        plist.append( ('rhi_user19',                c_float) )
        plist.append( ('rhi_user20',                c_float) )
        plist.append( ('rhi_user21',                c_float) )
        plist.append( ('rhi_user22',                c_float) )
        plist.append( ('pad_xx',                    c_char * 8) )
        plist.append( ('rhi_user23',                c_float) )
        plist.append( ('rhi_user24',                c_float) )
        plist.append( ('pad_xx',                    c_char * 60) )
        plist.append( ('rhi_user25',                c_float) )
        plist.append( ('rhi_user26',                c_float) )
        plist.append( ('rhi_user27',                c_float) )
        plist.append( ('rhi_user28',                c_float) )
        plist.append( ('rhi_user29',                c_float) )
        plist.append( ('rhi_user30',                c_float) )
        plist.append( ('rhi_user31',                c_float) )
        plist.append( ('rhi_user32',                c_float) )
        plist.append( ('rhi_user33',                c_float) )
        plist.append( ('rhi_user34',                c_float) )
        plist.append( ('rhi_user35',                c_float) )
        plist.append( ('rhi_user36',                c_float) )
        plist.append( ('rhi_user37',                c_float) )
        plist.append( ('rhi_user38',                c_float) )
        plist.append( ('rhi_user39',                c_float) )
        plist.append( ('rhi_user40',                c_float) )
        plist.append( ('rhi_user41',                c_float) )
        plist.append( ('rhi_user42',                c_float) )
        plist.append( ('rhi_user43',                c_float) )
        plist.append( ('rhi_user44',                c_float) )
        plist.append( ('rhi_user45',                c_float) )
        plist.append( ('rhi_user46',                c_float) )
        plist.append( ('rhi_user47',                c_float) )
        plist.append( ('rhi_user48',                c_float) )
        plist.append( ('pad_xx',                    c_char * 76) )
        plist.append( ('rhi_ctr_R',                 c_float) )
        plist.append( ('rhi_ctr_A',                 c_float) )
        plist.append( ('rhi_ctr_S',                 c_float) )
        plist.append( ('pad_xx',                    c_char * 12) )
        plist.append( ('rhi_tlhc_R',                c_float) )
        plist.append( ('rhi_tlhc_A',                c_float) )
        plist.append( ('rhi_tlhc_S',                c_float) )
        plist.append( ('rhi_trhc_R',                c_float) )
        plist.append( ('rhi_trhc_A',                c_float) )
        plist.append( ('rhi_trhc_S',                c_float) )
        plist.append( ('rhi_brhc_R',                c_float) )
        plist.append( ('rhi_brhc_A',                c_float) )
        plist.append( ('rhi_brhc_S',                c_float) )
        plist.append( ('pad_xx',                    c_char * 300) )
        plist.append( ('rhi_tr',                    c_int) )
        plist.append( ('rhi_ti',                    c_int) )
        plist.append( ('rhi_te',                    c_int) )
        plist.append( ('pad_xx',                    c_char * 40) )
        plist.append( ('rhi_rawrunnum',             c_int) )
        plist.append( ('pad_xx',                    c_char * 266) )
        plist.append( ('rhi_numecho',               c_short) )
        plist.append( ('pad_xx',                    c_char * 32) )
        plist.append( ('rhi_mr_flip',               c_short) )
        plist.append( ('pad_xx',                    c_char * 20) )
        plist.append( ('rhi_ctyp',                  c_short) )
        plist.append( ('pad_xx',                    c_char * 64) )
        plist.append( ('rhi_freq_dir',              c_short) )
        plist.append( ('pad_xx',                    c_char * 130) )
        plist.append( ('rhi_psdname',               c_char * 33) )
        plist.append( ('pad_xx',                    c_char * 84) )
        plist.append( ('rhi_cname',                 c_char * 17) )
        plist.append( ('pad_xx',                    c_char * 51) )
        plist.append( ('rhi_image_uid',             c_char * 32) )




    return plist
        



#        else if ( (int)(this->pfileVersion) == 26 ) {
#
#         offsets.assign (" \
#             rhr.rh_rdbm_rev                    , FLOAT_4, 1   , 0,\
#             rhr.rh_scan_date                   , CHAR   , 10  , 92,\
#             rhr.rh_scan_time                   , CHAR   , 8   , 102,\
#             rhr.rh_npasses                     , INT_2  , 1   , 140,\
#             rhr.rh_nslices                     , UINT_2 , 1   , 144,\
#             rhr.csi_dims                       , INT_2  , 1   , 436,\
#             rhr.rh_dab[0].start_rcv            , INT_2  , 1   , 264,\
#             rhr.rh_dab[0].stop_rcv             , INT_2  , 1   , 266,\
#             rhr.rh_dab[1].start_rcv            , INT_2  , 1   , 268,\
#             rhr.rh_dab[1].stop_rcv             , INT_2  , 1   , 270,\
#             rhr.rh_dab[2].start_rcv            , INT_2  , 1   , 272,\
#             rhr.rh_dab[2].stop_rcv             , INT_2  , 1   , 274,\
#             rhr.rh_dab[3].start_rcv            , INT_2  , 1   , 276,\
#             rhr.rh_dab[3].stop_rcv             , INT_2  , 1   , 278,\
#             rhr.rh_data_collect_type           , INT_2  , 1   , 132,\
#             rhr.rh_file_contents               , INT_2  , 1   , 120,\
#             rhr.rdb_hdr_off_data               , INT_4  , 1   , 4,\
#             rhr.rh_frame_size                  , UINT_2 , 1   , 156,\
#             rhr.rh_point_size                  , INT_2  , 1   , 158,\
#             rhr.rh_ps_mps_freq                 , UINT_4 , 1   , 488,\
#             rhr.rh_user_usage_tag              , UINT_4 , 1   , 924,\
#             rhr.roilenx                        , FLOAT_4, 1   , 444,\
#             rhr.roileny                        , FLOAT_4, 1   , 448,\
#             rhr.roilenz                        , FLOAT_4, 1   , 452,\
#             rhr.spectral_width                 , FLOAT_4, 1   , 432,\
#             rhr.xcsi                           , INT_2  , 1   , 438,\
#             rhr.ycsi                           , INT_2  , 1   , 440,\
#             rhr.zcsi                           , INT_2  , 1   , 442,\
#             rhr.rh_logo                        , CHAR   , 10  , 110,\
#             rhr.rh_raw_pass_size               , LINT_8 , 1   , 1540,\
#             rhr.rh_user0                       , FLOAT_4, 1   , 280,\
#             rhr.rh_user1                       , FLOAT_4, 1   , 284,\
#             rhr.rh_user2                       , FLOAT_4, 1   , 288,\
#             rhr.rh_user3                       , FLOAT_4, 1   , 292,\
#             rhr.rh_user4                       , FLOAT_4, 1   , 296,\
#             rhr.rh_user5                       , FLOAT_4, 1   , 300,\
#             rhr.rh_user6                       , FLOAT_4, 1   , 304,\
#             rhr.rh_user7                       , FLOAT_4, 1   , 308,\
#             rhr.rh_user8                       , FLOAT_4, 1   , 312,\
#             rhr.rh_user9                       , FLOAT_4, 1   , 316,\
#             rhr.rh_user10                      , FLOAT_4, 1   , 320,\
#             rhr.rh_user11                      , FLOAT_4, 1   , 324,\
#             rhr.rh_user12                      , FLOAT_4, 1   , 328,\
#             rhr.rh_user13                      , FLOAT_4, 1   , 332,\
#             rhr.rh_user14                      , FLOAT_4, 1   , 336,\
#             rhr.rh_user15                      , FLOAT_4, 1   , 340,\
#             rhr.rh_user16                      , FLOAT_4, 1   , 344,\
#             rhr.rh_user17                      , FLOAT_4, 1   , 348,\
#             rhr.rh_user18                      , FLOAT_4, 1   , 352,\
#             rhr.rh_user19                      , FLOAT_4, 1   , 356,\
#             rhr.rh_user20                      , FLOAT_4, 1   , 936,\
#             rhr.rh_user21                      , FLOAT_4, 1   , 940,\
#             rhr.rh_user22                      , FLOAT_4, 1   , 944,\
#             rhr.rh_user23                      , FLOAT_4, 1   , 948,\
#             rhr.rh_user24                      , FLOAT_4, 1   , 952,\
#             rhr.rh_user25                      , FLOAT_4, 1   , 956,\
#             rhr.rh_user26                      , FLOAT_4, 1   , 960,\
#             rhr.rh_user27                      , FLOAT_4, 1   , 964,\
#             rhr.rh_user28                      , FLOAT_4, 1   , 968,\
#             rhr.rh_user29                      , FLOAT_4, 1   , 972,\
#             rhr.rh_user30                      , FLOAT_4, 1   , 976,\
#             rhr.rh_user31                      , FLOAT_4, 1   , 980,\
#             rhr.rh_user32                      , FLOAT_4, 1   , 984,\
#             rhr.rh_user33                      , FLOAT_4, 1   , 988,\
#             rhr.rh_user34                      , FLOAT_4, 1   , 992,\
#             rhr.rh_user35                      , FLOAT_4, 1   , 996,\
#             rhr.rh_user36                      , FLOAT_4, 1   , 1000,\
#             rhr.rh_user37                      , FLOAT_4, 1   , 1004,\
#             rhr.rh_user38                      , FLOAT_4, 1   , 1008,\
#             rhr.rh_user39                      , FLOAT_4, 1   , 1012,\
#             rhr.rh_user40                      , FLOAT_4, 1   , 1016,\
#             rhr.rh_user41                      , FLOAT_4, 1   , 1020,\
#             rhr.rh_user42                      , FLOAT_4, 1   , 1024,\
#             rhr.rh_user43                      , FLOAT_4, 1   , 1028,\
#             rhr.rh_user44                      , FLOAT_4, 1   , 1032,\
#             rhr.rh_user45                      , FLOAT_4, 1   , 1036,\
#             rhr.rh_user46                      , FLOAT_4, 1   , 1040,\
#             rhr.rh_user47                      , FLOAT_4, 1   , 1044,\
#             rhr.rh_user48                      , FLOAT_4, 1   , 1048,\
#             rhi.psdname                        , CHAR   , 33  , 199812,\
#             rhi.scanspacing                    , FLOAT_4, 1   , 198500,\
#             rhi.te                             , INT_4  , 1   , 199244,\
#             rhi.ti                             , INT_4  , 1   , 199240,\
#             rhi.tr                             , INT_4  , 1   , 199236,\
#             rhi.tlhc_A                         , FLOAT_4, 1   , 198904,\
#             rhi.tlhc_R                         , FLOAT_4, 1   , 198900,\
#             rhi.tlhc_S                         , FLOAT_4, 1   , 198908,\
#             rhi.t                              , INT_4  , 1   , 199236,\
#             rhi.trhc_A                         , FLOAT_4, 1   , 198916,\
#             rhi.trhc_R                         , FLOAT_4, 1   , 198912,\
#             rhi.trhc_S                         , FLOAT_4, 1   , 198920,\
#             rhi.user0                          , FLOAT_4, 1   , 198536,\
#             rhi.user1                          , FLOAT_4, 1   , 198540,\
#             rhi.user2                          , FLOAT_4, 1   , 198544,\
#             rhi.user3                          , FLOAT_4, 1   , 198548,\
#             rhi.user4                          , FLOAT_4, 1   , 198552,\
#             rhi.user5                          , FLOAT_4, 1   , 198556,\
#             rhi.user6                          , FLOAT_4, 1   , 198560,\
#             rhi.user7                          , FLOAT_4, 1   , 198564,\
#             rhi.user8                          , FLOAT_4, 1   , 198568,\
#             rhi.user9                          , FLOAT_4, 1   , 198572,\
#             rhi.user10                         , FLOAT_4, 1   , 198576,\
#             rhi.user11                         , FLOAT_4, 1   , 198580,\
#             rhi.user12                         , FLOAT_4, 1   , 198584,\
#             rhi.user13                         , FLOAT_4, 1   , 198588,\
#             rhi.user14                         , FLOAT_4, 1   , 198592,\
#             rhi.user15                         , FLOAT_4, 1   , 198596,\
#             rhi.user16                         , FLOAT_4, 1   , 198600,\
#             rhi.user17                         , FLOAT_4, 1   , 198604,\
#             rhi.user18                         , FLOAT_4, 1   , 198608,\
#             rhi.user19                         , FLOAT_4, 1   , 198612,\
#             rhi.user20                         , FLOAT_4, 1   , 198616,\
#             rhi.user21                         , FLOAT_4, 1   , 198620,\
#             rhi.user22                         , FLOAT_4, 1   , 198624,\
#             rhi.user23                         , FLOAT_4, 1   , 198636,\
#             rhi.user24                         , FLOAT_4, 1   , 198640,\
#             rhi.user25                         , FLOAT_4, 1   , 198704,\
#             rhi.user26                         , FLOAT_4, 1   , 198708,\
#             rhi.user27                         , FLOAT_4, 1   , 198712,\
#             rhi.user28                         , FLOAT_4, 1   , 198716,\
#             rhi.user29                         , FLOAT_4, 1   , 198720,\
#             rhi.user30                         , FLOAT_4, 1   , 198724,\
#             rhi.user31                         , FLOAT_4, 1   , 198728,\
#             rhi.user32                         , FLOAT_4, 1   , 198732,\
#             rhi.user33                         , FLOAT_4, 1   , 198736,\
#             rhi.user34                         , FLOAT_4, 1   , 198740,\
#             rhi.user35                         , FLOAT_4, 1   , 198744,\
#             rhi.user36                         , FLOAT_4, 1   , 198748,\
#             rhi.user37                         , FLOAT_4, 1   , 198752,\
#             rhi.user38                         , FLOAT_4, 1   , 198756,\
#             rhi.user39                         , FLOAT_4, 1   , 198760,\
#             rhi.user40                         , FLOAT_4, 1   , 198764,\
#             rhi.user41                         , FLOAT_4, 1   , 198768,\
#             rhi.user42                         , FLOAT_4, 1   , 198772,\
#             rhi.user43                         , FLOAT_4, 1   , 198776,\
#             rhi.user44                         , FLOAT_4, 1   , 198780,\
#             rhi.user45                         , FLOAT_4, 1   , 198784,\
#             rhi.user46                         , FLOAT_4, 1   , 198788,\
#             rhi.user47                         , FLOAT_4, 1   , 198792,\
#             rhi.user48                         , FLOAT_4, 1   , 198796,\
#             rhi.cname                          , CHAR   , 17  , 199929,\
#             rhi.brhc_A                         , FLOAT_4, 1   , 198928,\
#             rhi.brhc_R                         , FLOAT_4, 1   , 198924,\
#             rhi.brhc_S                         , FLOAT_4, 1   , 198932,\
#             rhi.ctr_A                          , FLOAT_4, 1   , 198880,\
#             rhi.ctr_R                          , FLOAT_4, 1   , 198876,\
#             rhi.ctr_S                          , FLOAT_4, 1   , 198884,\
#             rhi.dfov                           , FLOAT_4, 1   , 198484,\
#             rhi.freq_dir                       , INT_2  , 1   , 199680,\
#             rhi.ctyp                           , INT_2  , 1   , 199614,\
#             rhi.loc                            , FLOAT_4, 1   , 198504,\
#             rhi.mr_flip                        , INT_2  , 1   , 199592,\
#             rhi.nex                            , FLOAT_4, 1   , 198512,\
#             rhi.numecho                        , INT_2  , 1   , 199558,\
#             rhi.image_uid                      , UID    , 32  , 199997,\
#             rhi.rawrunnum                      , INT_4  , 1   , 199288,\
#             rhe.ex_datetime                    , INT_4  , 1   , 194240,\
#             rhe.ex_no                          , UINT_2 , 1   , 194356,\
#             rhe.magstrength                    , INT_4  , 1   , 194232,\
#             rhe.patid                          , CHAR   , 65  , 195249,\
#             rhe.patname                        , CHAR   , 65  , 195184,\
#             rhe.refphy                         , CHAR   , 65  , 194717,\
#             rhe.reqnum                         , CHAR   , 17  , 195314,\
#             rhe.study_uid                      , UID    , 32  , 195088,\
#             rhe.dateofbirth                    , CHAR   , 9   , 195331,\
#             rhe.patsex                         , INT_2  , 1   , 194380,\
#             rhe.hospname                       , CHAR   , 33  , 195011,\
#             rhe.ex_sysid                       , CHAR   , 9   , 194980,\
#             rhe.uniq_sys_id                    , CHAR   , 16  , 195052,\
#             rhe.ex_verscre                     , CHAR   , 2   , 195048,\
#             rhs.se_no                          , INT_4  , 1   , 196356,\
#             rhs.se_desc                        , CHAR   , 65  , 196602,\
#             rhs.entry                          , INT_4  , 1   , 196264,\
#             rhs.position                       , INT_4  , 1   , 196260,\
#             rhs.series_uid                     , UID    , 32  , 196765,\
#             rhs.landmark_uid                   , UID    , 32  , 196797,\
#             rhs.anref                          , CHAR   , 3   , 196685,\


#----------------------------------------------------------
# Test routines
#----------------------------------------------------------



def main():

    # v16 data - 1 chan
    fname = 'C:\\Users\\bsoher\\code\\repository_svn\\sample_data\\example_ge_svs_1ch_pom\\P29184.7'
    
    # v11 data - 8 chan
#    fname = 'C:\\Users\\bsoher\\code\\repository_svn\\sample_data\\example_ge_pom_multi-channel\\P01024.7'

    # v20 data - 32 chan
#    fname = 'C:\\Users\\bsoher\\code\\repository_svn\\sample_data\\example_ge_svs_32ch_gregor\\P32256.7'

    hdr = Pfile(fname)
    hdr.dump_header()

    bob = 10
    bob += 1




if __name__ == '__main__':

    main()

