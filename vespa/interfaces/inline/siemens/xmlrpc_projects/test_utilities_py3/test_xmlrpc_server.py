#!/usr/bin/env python

# Python imports
from __future__ import division

import os
import sys
import xmlrpclib
from SimpleXMLRPCServer import SimpleXMLRPCServer


# _EXIT_TIMER_DELAY defines a recommended delay between when client calls 
#  "quit" and when this server will most likely actually have quitted. This 
#  is in case the client is trying to quit and re-start the sever. 
_EXIT_TIMER_DELAY = 0.1


#------------------------------------------------------------------------------
# Local Handlers - direct methods and instances

def adder_function(x,y):
    """A minimalist method for testing."""
    print "test_xmlrpc_server: adder_function() --"
    val = x+y
    return val


def kill_server():
    """
    Called by the ListenRequestHandler when it is time for this server to exit
    
    """
    global server
    server.quit = True
    

class ListenRequestHandler(object):
    """
    This class does most of the interesting work. It defines and provides 
    the API that's exposed to the world through XMLRPC.
    """
    def __init__(self, kill_server, verbose=False):
        
        self.verbose = verbose
        
        # kill_server is a function that this class can use to tell the listener to quit.
        self._kill_server = kill_server

    def div(self, x, y):
        """A minimalist method for testing - div()."""
        if self.verbose: print "test_xmlrpc_server: div() --"
        return x // y
    
    def add(self, x, y):
        """A minimalist method for testing - add()."""
        if self.verbose: print "test_xmlrpc_server: add() --"
        return x + y

    def bump2(self, x, y):
        ''' 
        Test for what is returned when more than 1 result returned.
        In C/C++ ICE, I was able to use xmlrpc_decompose_value() as:
           xmlrpc_decompose_value(&m_xmlrpc_env, result2, "(ii)", &res1, &res2);
        to parse a known return value set.
         
        '''
        if self.verbose: print "test_xmlrpc_server: bump2() --"
        return x+1, y+1

    def is_available(self):
        """
        Currently a no-op.

        To be called when the scanner is initializing the ICE program and
        wants to check if there is an XML-RPD server listening for it.
        """
        if self.verbose: print "test_simple_server: is_available() --"
        pass

    def quit(self):
        """
        Kills the listener. Should only be called from RtView.

        It returns a float that measures the number of seconds that will elapse
        before the app actually exits. RtView should not attempt to start
        a new listener until at least that amount of time has elapsed.
        
        """
        if self.verbose: print "listener: quit() --"
        #
        # I would like to call sys.exit() directly here, but an XMLRPC call 
        # has to return, otherwise xmlrpclib raises an error. So instead I set
        # a flag via the call to _kill_server() and ask the caller to wait for
        # _EXIT_TIMER_DELAY seconds.
        #
        self._kill_server()
        return _EXIT_TIMER_DELAY



#------------------------------------------------------------------------------
# Create server

class MyXmlRpcServer(SimpleXMLRPCServer):
    quit = False

    def serve_until_quit(self):
        self.quit = False
        while not self.quit:
            print "MyXmlRpcServer: handling next request"
            self.handle_request()
        print "MyXmlRpcServer: self.quit set, exiting..."



address = "127.0.0.1"
port    = 8080

server = MyXmlRpcServer((address, port), logRequests=True, allow_none=True)
server.register_introspection_functions()

rh = ListenRequestHandler(kill_server, verbose=True)

server.register_instance(rh)
server.register_function(kill_server)
server.register_function(adder_function, 'adder_function')

print "listener is listening..."
server.serve_until_quit()
print "listener is exiting."

   


# # Run the server's main loop
# try:
#     print 'Use Control-C to exit'
#     server.serve_forever()
# 
# except KeyboardInterrupt:
#     print 'Exiting'






#-------------------------------------
# Here starts another example from
#  https://pymotw.com/2/xmlrpclib/


# from SimpleXMLRPCServer import SimpleXMLRPCServer
# from xmlrpclib import Binary
# import datetime
# 
# server = SimpleXMLRPCServer(('localhost', 9000), logRequests=True, allow_none=True)
# server.register_introspection_functions()
# server.register_multicall_functions()
# 
# class ExampleService:
#     
#     def ping(self):
#         """Simple function to respond when called to demonstrate connectivity."""
#         return True
#         
#     def now(self):
#         """Returns the server current date and time."""
#         return datetime.datetime.now()
# 
#     def show_type(self, arg):
#         """Illustrates how types are passed in and out of server methods.
#         
#         Accepts one argument of any type.  
#         Returns a tuple with string representation of the value, 
#         the name of the type, and the value itself.
#         """
#         return (str(arg), str(type(arg)), arg)
# 
#     def raises_exception(self, msg):
#         "Always raises a RuntimeError with the message passed in"
#         raise RuntimeError(msg)
# 
#     def send_back_binary(self, bin):
#         "Accepts single Binary argument, unpacks and repacks it to return it"
#         data = bin.data
#         response = Binary(data)
#         return response
# 
# server.register_instance(ExampleService())
# 
# try:
#     print 'Use Control-C to exit'
#     server.serve_forever()
# except KeyboardInterrupt:
#     print 'Exiting'
