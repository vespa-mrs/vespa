# Python modules

#import exceptions


class TransformRunException(Exception):
    """
    Raised when am error occurs in one of the user defined algorithms in the
    Pulse app. Includes an error string, that should be displayed to the
    user.
    """
    
    def __init__(self, error_str, error_code):
        self._message = error_str
        self.code    = error_code 


    __doc = """Getter/setter for the message associated with this exception"""
    # exceptions.Exception has a .message attribute, but it's deprecated in
    # Python 2.6 so we re-implement it ourselves.
    def __get_message(self):
        return self._message
    def __set_message(self, message):
        self._message = message

    message = property(__get_message, __set_message, doc=__doc)
