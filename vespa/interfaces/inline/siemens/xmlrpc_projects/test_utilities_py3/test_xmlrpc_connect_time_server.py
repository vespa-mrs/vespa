"""
Stole this example from ...

https://docs.python.org/2/library/simplexmlrpcserver.html


"""
from __future__ import division

from SimpleXMLRPCServer import SimpleXMLRPCServer
from SimpleXMLRPCServer import SimpleXMLRPCRequestHandler

# Restrict to a particular path.
class RequestHandler(SimpleXMLRPCRequestHandler):
    rpc_paths = ('/RPC2',)

# Create server
server = SimpleXMLRPCServer(("localhost", 8000), requestHandler=RequestHandler)
server.register_introspection_functions()

# Register an instance; all the methods of the instance are
# published as XML-RPC methods .
class MyFuncs:
    def div(self, x, y):
        return x // y

    def add(self, x, y):
        return x + y

server.register_instance(MyFuncs())

# Run the server's main loop
server.serve_forever()