# Python modules

import time
import datetime

# 3rd party modules

# Our modules

# DateTime objects have a .strftime() method that accepts a format string. 
# Below are some format strings for you to use. Please consider using one
# of these rather than creating a new format so that Vespa is at least 
# somewhat consistent about date formatting.
# To use, pass the timestamp format constant to .strftime(), e.g. --
#    print my_datetime_instance.strftime(util_time.DISPLAY_DATE_FORMAT)
#
# Note that DateTime objects have a method called .isoformat() which returns
# a string in ISO 8601 format. You can get similar results by using the
# ISO_TIMESTAMP_FORMAT constant below. e.g. if you already have a datetime
# object --
#    print my_datetime_instance.isoformat()
# is similar to --
#    print my_datetime_instance.strftime(util_time.ISO_TIMESTAMP_FORMAT)
# If you want the current time in ISO format -- 
#    print util_time.now().isoformat()
# or -- 
#    print util_time.now(util_time.ISO_TIMESTAMP_FORMAT)
#
# You can use whichever syntax you prefer. Be aware that .isoformat() includes
# microsecond & time zone info if available which may not be what you want.

ISO_TIMESTAMP_FORMAT = "%Y-%m-%dT%H:%M:%S"

ISO_DATE_FORMAT = "%Y-%m-%d"

# DISPLAY_DATE_FORMAT is a human-friendly date format. e.g. 21 April, 2010
# It spells out the month so there's no confusion over mm/dd/yy versus
# dd/mm/yy and so forth.
DISPLAY_DATE_FORMAT = "%d %B, %Y"

# DISPLAY_TIMESTAMP_FORMAT is the human-friendly date plus the "locale's
# appropriate time representation".
DISPLAY_TIMESTAMP_FORMAT = "%d %B, %Y %X"

# CLONE_TIMESTAMP_FORMAT is the human-friendly date plus the "locale's
# appropriate time representation" usable as the name of a cloned object
# (experiment, pulse project, metabolite, etc.)
CLONE_TIMESTAMP_FORMAT = "%Y-%m-%dT%H_%M_%S"


def now(format=None):
    """Returns the current local date & time.

       By default, the function returns a datetime object. If the format
       parameter is a string appropriate for time.strftime(), the function
       returns the current time as a string formatted per the format string.

       If a datetime object is returned, resolution is accurate only to one
       second. (That is, the object's .microsecond attribute is always 0.)
       This is by design. If you need sub-second accuracy, then this will
       return a datetime object with an accurate microsecond attribute:
           current = util_time.now().now()
       """
    if format:
        return time.strftime(format, time.localtime())
    else:
        return datetime.datetime(*time.localtime()[:6])


def datetime_from_iso(iso_timestamp):
    """Given an ISO timestamp string, returns an equivalent datetime object.
    The string must include a date and time. Microseconds are optional. If
    present, they're present in the datetime object.
    Examples --
        >>> util_time.datetime_from_iso("2010-06-11T16:16:24.387335")
        datetime.datetime(2010, 6, 11, 16, 16, 24, 387335)
        >>> util_time.datetime_from_iso("2010-06-11T16:16:24")
        datetime.datetime(2010, 6, 11, 16, 16, 24)
        >>> 
    """
    # Figure out whether or not usecs are present.
    i = iso_timestamp.rfind(".")
    if i == -1:
        microseconds = 0
    else:
        # Extract microsecond info and remove it from the string because 
        # microseconds confuse & upset time.strptime().
        microseconds = int(iso_timestamp[i + 1:])
        iso_timestamp = iso_timestamp[:i]
        
    params = time.strptime(iso_timestamp, ISO_TIMESTAMP_FORMAT)[:6]
    params = list(params) + [microseconds]
                        
    return datetime.datetime(*params)


def filename_timestamp():
    """Returns a timestamp appropriate for inclusion as part of a filename.
    The timestamp includes microseconds, and so subsequent calls to this
    function are guaranteed to return different filenames.
    """
    # FILENAME_TIMESTAMP_FORMAT is hidden inside this function because it
    # can't be used on its own to create unique filenames due to its lack of
    # sub-second information. It's easy enough to create multiple files in
    # one second, so without sub-second information appended the returned
    # filenames are not guaranteed to be different.
    FILENAME_TIMESTAMP_FORMAT = "%Y%m%d.%H%M%S"

    current = now().now()

    # Python's strftime relies on the underlying C library's strftime and
    # so can only guarantee the existence of the strftime formatting codes
    # that are guaranteed by C89. Unfortunately, this doesn't include
    # any sub-second values so I have to tack the microseconds on myself.
    s = current.strftime(FILENAME_TIMESTAMP_FORMAT)

    return s + ".%d" % current.microsecond
