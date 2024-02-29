#!/usr/bin/env python

import xmlrpc.client

s = xmlrpc.client.ServerProxy('http://127.0.0.1:8080')

print('bump2 - ', s.bump2(2,3))  # Returns [3,4]
print('add  - ', s.add(2,3))  # Returns 5
print('div   - ', s.div(5,2))  # Returns 5//2 = 2
#print('adder_function  - ', s.adder_function(7,9))  # Returns 16
# Print list of available methods
print(s.system.listMethods())

#img = s.get_fake_result(1)      # 0-return fake data, 1-return all zeros


bob = 10
bob += 1

s.quit()

bob += 1



#------------------------------------------
# Here starts another example, that matches the one in text_xmlrpc_server at the BOTTOM
#  from  https://pymotw.com/2/xmlrpclib/


# import xmlrpclib
# 
# server = xmlrpclib.ServerProxy('http://localhost:9000')
# 
# print 'Ping:', server.ping()
# 
# 
# server = xmlrpclib.ServerProxy('http://localhost:9000', allow_none=True)
# print 'Allowed:', server.show_type(None)
# 
# server = xmlrpclib.ServerProxy('http://localhost:9000', allow_none=False)
# print 'Not allowed:', server.show_type(None)
# 
# 
# 
# for t, v in [ ('boolean', True), 
#               ('integer', 1),
#               ('floating-point number', 2.5),
#               ('string', 'some text'), 
#               ('datetime', datetime.datetime.now()),
#               ('array', ['a', 'list']),
#               ('array', ('a', 'tuple')),
#               ('structure', {'a':'dictionary'}),
#             ]:
#     print '%-22s:' % t, server.show_type(v)
#     
#     
# 
# 
# s = 'This is a string with control characters' + '\0'
# print 'Local string:', s
# 
# data = xmlrpclib.Binary(s)
# print 'As binary:', server.send_back_binary(data)
# 
# print 'As string:', server.show_type(s)
# 
# 
# 
# try:
#     server.raises_exception('A message')
# except Exception, err:
#     print 'Fault code:', err.faultCode
#     print 'Message   :', err.faultString