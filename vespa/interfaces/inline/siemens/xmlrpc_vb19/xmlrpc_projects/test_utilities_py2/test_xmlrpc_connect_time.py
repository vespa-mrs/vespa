from __future__ import division

import xmlrpclib
import socket


class TestConnection(object):
    
    def __init__(self):
        
        self.address = "http://127.0.0.1:8000"
        self._server = xmlrpclib.ServerProxy(self.address)
        

    def call_server_with_startup(self):
        a = xmlrpclib.ServerProxy(self.address)
        try:
            a._()   # Call a fictive method.
        except xmlrpclib.Fault:
            # connected to the server and the method doesn't exist which is expected.
            pass
        except socket.error:
            # Not connected ; socket error mean that the service is unreachable.
            print "call_server_with_startup : connect failed"


    def call_server_already_running(self):
        try:
            self._server._()
        except xmlrpclib.Fault:
            pass
        except socket.error:
            print "call_server_already_running : connect failed"

    


def main(iterations=10000):
    """
    Blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah 
    blah blah blah blah  

    """

    item = TestConnection()

    for i in range(iterations):
        item.call_server_with_startup()
        
    for i in range(iterations):
        item.call_server_already_running()
        
    

if __name__ == "__main__":
    #main()
      
    import cProfile
    cProfile.run('main()')






