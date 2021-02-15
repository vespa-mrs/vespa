"""
Routine for reading an <XProtocol> formatted header string from a Siemens Twix
file (only tested on VB at this point).

This version is only PARTLY WORKING. I don't know the actual rules for making
an XProtocol file, so it's all been trial and error. At this point the program
still 'miscounts' the right and left curlies '{ and }' and exits out too early
if I try to parse the whole string.

If I search for the start of a section (as shown in the test case below) as:

    pattern = '<ParamMap."MEAS">'
    twix = TwixRaid()
    twix.populate_from_file(fname)
    evps = twix.current.evps
    evp21 = evps[2][1].replace('\n', '')        # need to remove line breaks here
    evp21 = evp21[evp21.find(pattern):]
    p = XprotParserDict()
    out = p.parse(evp21)

Many of the sections (particularly '<ParamMap."sWiPMemBlock">') will parse out
correctly in their entirety.  But, some (like MEAS) will crap out early.  So.
Still early days.

2020-08-26 - bjs, I've stopped beating my head on this for now.

Here are some fun regex that I created to figure out how to parse the XProtocol
lines for useful info.  Not all work.

#        # bob1 = re.findall('<Param(\w+)\."(\w+)">\s*{([^}]*)', evps21)
#        # bob1 = re.findall('<Param(?:Bool|Long|String)\."(\w+)">\s*{([^}]*)', evps[0][1])
#        # bob2 = re.findall('<ParamDouble\."(\w+)">\s*{\s*(<Precision>\s*[0-9]*)?\s*([^}]*)', evps[0][1])
#        #
#        # tom1 = re.findall('<Param(?:Bool|Long|String)\."(\w+)">\s*{([^}]*)', evps[2][1])
#        # tom2 = re.findall('<ParamDouble\."(\w+)">\s*{\s*(<Precision>\s*[0-9]*)?\s*([^}]*)', evps[2][1])
#        # tom3 = re.findall('<ParamMap\."(\w+)">\s*{([^}]*)', evps[2][1])

"""


# Python modules
from __future__ import print_function
import re
from collections import OrderedDict

# 3rd party modules

# Our modules
from vespa.common.twix_parser import TwixRaid

import sre_compile, sre_parse


# --------------------------------------------------------------------
# experimental stuff (see python-dev discussions for details)

class Scanner:
    def __init__(self, lexicon, flags=0):
        from sre_constants import BRANCH, SUBPATTERN
        if isinstance(flags, re.RegexFlag):
            flags = flags.value
        self.lexicon = lexicon
        # combine phrases into a compound pattern
        p = []
        s = sre_parse.Pattern()
        s.flags = flags
        for phrase, action in lexicon:
            gid = s.opengroup()
            p.append(sre_parse.SubPattern(s, [
                (SUBPATTERN, (gid, 0, 0, sre_parse.parse(phrase, flags))),
                ]))
            s.closegroup(gid, p[-1])
        p = sre_parse.SubPattern(s, [(BRANCH, (None, p))])
        self.scanner = sre_compile.compile(p)
    def scan(self, string):
        result = []
        append = result.append
        match = self.scanner.scanner(string).match
        i = 0
        while True:
            m = match()
            if not m:
                break
            j = m.end()
            if i == j:
                break
            action = self.lexicon[m.lastindex-1][1]
            if callable(action):
                self.match = m
                action = action(self, m.group())
            if action is not None:
                append(action)
                break                   # bjs - one change from version in re.
            i = j
        return result, string[i:]



class Node(OrderedDict):
    def __init__(self, parent=None, ltype='', lkey='' ):
        self.parent = parent
        self.lastkey = ltype
        self.lasttype = lkey
        self['index'] = 0

class XprotParserDict(object):
    def __init__(self):
        lefttop = '(?:<XProtocol>\s*{)'
        leftlim = '(?:<Limit[^>]*>\s*{)'
        leftkey = '(?:<Param[^>]*>\s*{)'
        leftdef = '(?:<Default> <Param[^>]*>\s*{)'
        leftval = '{'
        right = '}'
        self.scanner = Scanner([
            ( lefttop, self.lefttop),
            ( leftlim, self.leftlim),
            ( leftkey, self.leftkey),
            ( leftdef, self.leftdef),
            ( leftval, self.leftval),
            ( right, self.right),
            ( r"\s+", None),
            ( ".+?(?=(%s|%s|%s|%s|%s|%s|$))" % (right, lefttop, leftlim, leftkey, leftdef, leftval), self.other),
        ])
        self.re1 = re.compile('<Param(.*)\.')
        self.re2 = re.compile('\."(.*)">')
        self.re3 = re.compile('<([^>]*)>')
        self.re4 = re.compile('<Default> <Param(.*)\.')
        self.re5 = re.compile('<Precision>\s*[0-9]*\s*(.*)')
        self.re6 = re.compile('<Default>\s*-?\d+\.?\d*\s*(.*)')

        self.result = Node()
        self.current = self.result

    def parse(self, content):
        self.done_parsing = False
        self.scanner.scan(content)
        return self.result

    def lefttop(self, scanner, token):
        self.current['XProtocol'] = Node(self.current)
        self.current = self.current['XProtocol']

    def leftlim(self, scanner, token):
        new_key = self.re3.search(token).group(1)
        self.current.lastkey = new_key
        self.current.lasttype = 'String'    # default bjs for now
        self.current[new_key] = ''

    def leftkey(self, scanner, token):
        dtype = self.re1.search(token).group(1)
        new_key = self.re2.search(token).group(1)
        if new_key == '': new_key = 'Top'

        if dtype in ['String','Long','Double']:
            self.current.lastkey  = new_key
            self.current.lasttype = dtype
            self.current[new_key] = ''
        else:
            self.current[new_key] = Node(self.current)
            self.current = self.current[new_key]
            self.current['type'] = dtype
            self.current.lastkey = ''
            self.current.lasttype = 'generic'

    def leftdef(self, scanner, token):
        dtype = self.re4.search(token).group(1)
        new_key = self.re2.search(token).group(1)
        if new_key == '': new_key = 'DefaultMap'
        self.current[new_key] = Node(self.current)
        self.current = self.current[new_key]
        self.current['type'] = dtype
        self.current.lastkey = ''
        self.current.lasttype = 'generic'

    def leftval(self, scanner, token):
        if self.current.lastkey == '':
            new_key = 'value'+str(self.current['index'])
            self.current[new_key] = ''
            self.current.lastkey = new_key
            self.current.lasttype = 'generic_array'
            self.current['index'] += 1
        else:
            self.current[self.current.lastkey] = Node(self.current)
            self.current = self.current[self.current.lastkey]
            new_key = 'value'+str(self.current['index'])
            self.current[new_key] = ''
            self.current.lastkey = new_key
            self.current.lasttype = 'generic_array'
            self.current['index'] += 1

    def right(self, scanner, token):
        if self.current.lastkey == '':
            if 'index' in self.current.keys():
                del self.current['index']
            self.current = self.current.parent
            if self.current.parent is None:
                return self.result
        self.current.lastkey = ''
        self.current.lasttype = ''

    def other(self, scanner, token):
        if 'index' in self.current.keys():
            if self.current.lastkey == '':
                index = self.current['index']
                self.current['value'+str(index)] = token.strip()
                self.current['index'] += 1
            else:
                try:
                    if self.current.lasttype == 'Long':
                        tmp = token
                        r6 = self.re6.search(tmp)
                        tmp = tmp if r6 is None else r6.group(1)       # for <Default>
                        val = [int(item) for item in tmp.split()]
                    elif self.current.lasttype == 'Double':
                        tmp = token
                        r5 = self.re5.search(tmp)
                        tmp = tmp if r5 is None else r5.group(1)       # for <Precision>
                        r6 = self.re6.search(tmp)
                        tmp = tmp if r6 is None else r6.group(1)       # for <Default>
                        val = [float(item) for item in tmp.split()]
                    else:
                        val = [token,]
                    val = val if val != [] else ['', ]
                    self.current[self.current.lastkey] = val if len(val) > 1 else val[0]
                except:
                    bob = 10        # for error checking





#------------------------------------------------------------------------------
# test code

def _test():

    fname = r"D:\Users\bsoher\code\repository_svn\sample_data\siemens_twix_vb19_svs_se\meas_MID00165_FID50843_eja_svs_press.dat"

    twix = TwixRaid()
    twix.populate_from_file(fname)
    evps = twix.current.evps

    evp21 = evps[2][1].replace('\n', '')

    # pattern = '<ParamMap."sWiPMemBlock">'
    # pattern = '<ParamMap."DICOM">'
    # pattern = '<ParamMap."">'
    pattern = '<ParamMap."MEAS">'
    indx = evp21.find(pattern)

    p = XprotParserDict()
    out = p.parse(evp21[indx:])

    from pprint import pprint
    logFile = open("D:\\users\\bsoher\\log_xprotocol.txt", 'w')
    pprint(out, logFile, indent=2, width=120)
    logFile.close()

    bob = 10
    bob += 1



if __name__ == '__main__':
    """ Just testing the read_raw method here """

    _test()


# class NestedParser(object):
#     def __init__(self, left='\{', right='\}'):
#         self.scanner = Scanner([
#             ( left, self.left),
#             ( right, self.right),
#             ( r"\s+", None),
#             ( ".+?(?=(%s|%s|$))" % (right, left), self.other),
#         ])
#         self.result = Node()
#         self.current = self.result
#
#     def parse(self, content):
#         self.scanner.scan(content)
#         return self.result
#
#     def left(self, scanner, token):
#         new = Node(self.current)
#         self.current.append(new)
#         self.current = new
#
#     def right(self, scanner, token):
#         self.current = self.current.parent
#
#     def other(self, scanner, token):
#         self.current.append(token.strip())


