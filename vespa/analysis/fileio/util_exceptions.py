# Python modules


# 3rd party modules


# Our modules



class MultifileAttributeMismatchError(Exception):
    """Raised when a class in this module opens multiple files and the
    attributes (dims and sweep width) don't match.
    """
    pass


class MultifileTypeMismatchError(Exception):
    """Raised when a class in this module opens multiple files and the
    the object types in the returned list are not all of the same type.
    """
    pass


class UnsupportedDimensionalityError(Exception):
    """Raised when a class in this module opens a file in which all
    dimensions are > 1.
    """
    pass


class IncorrectDimensionalityError(Exception):
    """Raised when a class in this module opens a file in which at least one
    dimensions that should be 1 is not 1.
    """
    pass


class SIDataError(UnsupportedDimensionalityError):
    """Raised when a class in this module opens a file containing SI data.
    This is a special case of UnsupportedDimensionalityError. We expect
    to support SI data eventually, and then this will go away.
    """
    pass


class OpenFileAttributeMismatchError(Exception):
    """Raised when a class in this module opens one or more files and the
    attributes (dims and sweep width) of the new files don't match the
    attributes of the currently open files.
    """
    pass


class OpenFileTypeMismatchError(Exception):
    """Raised when a class in this module opens one or more files and the
    object types in which the new files are stored don't match the object
    types of the currently open files.
    """
    pass


class FileNotFoundError(Exception):
    """
    Raised when a reader module can't find a matching data file for a params
    file or vice versa.
    """
    pass


class OpenFileUserReadRawError(Exception):
    """Raised when a user derived class want to indicate that there is a
    problem within their read_raw method.
    """
    pass

class IncompleteHeaderParametersError(Exception):
    """Raised when a user derived class want to indicate that there is a
    problem within their read_raw method.
    """
    pass

class UnsupportedPulseSequenceError(Exception):
    """Raised when a user selected file does not contain the desired pulse sequence """
    pass



