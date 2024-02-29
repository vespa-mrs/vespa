#!/usr/bin/env python

import xmlrpclib
from SimpleXMLRPCServer import SimpleXMLRPCServer


# Create server
server = SimpleXMLRPCServer(("127.0.0.1", 8080), logRequests=True, allow_none=True)
server.register_introspection_functions()

def adder_function(x,y):
    val = x+y
    return val

server.register_function(adder_function, 'add')

# Register an instance; all the methods of the instance are
# published as XML-RPC methods (in this case, just 'div').
class MyFuncs:
    
    def div(self, x, y):
        print "test_simple_server: div() --"
        return x // y
    
    def add(self, x, y):
        print "test_simple_server: add() --"
        return x + y

    def bump2(self, x, y):
        print "test_simple_server: bump2() --"
        return x + 1, y+1

    def is_available(self):
        """
        Currently a no-op.

        To be called when the scanner is initializing the ICE program and
        wants to check if there is an XML-RPD server listening for it.
        """
        print "test_simple_server: is_available() --"
        pass

server.register_instance(MyFuncs())

# Run the server's main loop
try:
    print 'Use Control-C to exit'
    server.serve_forever()

except KeyboardInterrupt:
    print 'Exiting'


