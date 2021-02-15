# PS - This is our version of the BSD-licensed configobj module.
# We have made one change to it.
# 1. We set the OPTION_DEFAULTS encoding to utf-8

import configobj

configobj.OPTION_DEFAULTS['encoding'] = 'utf-8'

class ConfigObj(configobj.ConfigObj):

    def __init__(self, infile=None, options=None, configspec=None, encoding=None,
                 interpolation=True, raise_errors=False, list_values=True,
                 create_empty=False, file_error=False, stringify=True,
                 indent_type=None, default_encoding=None, unrepr=False,
                 write_empty_values=False, _inspec=False):

        ''' call base class __init__ '''
        configobj.ConfigObj.__init__(self, infile=infile, options=options, configspec=configspec, encoding=encoding,
                                     interpolation=interpolation, raise_errors=raise_errors, list_values=list_values,
                                     create_empty=create_empty, file_error=file_error, stringify=stringify,
                                     indent_type=indent_type, default_encoding=default_encoding, unrepr=unrepr,
                                     write_empty_values=write_empty_values, _inspec=_inspec)



