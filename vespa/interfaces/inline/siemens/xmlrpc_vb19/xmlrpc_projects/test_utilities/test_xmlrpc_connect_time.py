

import xmlrpc.client
import socket


class TestConnection(object):
    
    def __init__(self):
        
        self.address = "http://127.0.0.1:8000"
        self._server = xmlrpc.client.ServerProxy(self.address)
        

    def call_server_with_startup(self):
        a = xmlrpc.client.ServerProxy(self.address)
        try:
            a._()   # Call a fictive method.
        except xmlrpc.client.Fault:
            # connected to the server and the method doesn't exist which is expected.
            pass
        except socket.error:
            # Not connected ; socket error mean that the service is unreachable.
            print("call_server_with_startup : connect failed")


    def call_server_already_running(self):
        try:
            self._server._()
        except xmlrpc.client.Fault:
            pass
        except socket.error:
            print("call_server_already_running : connect failed")

    


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






