# Python modules

import xml.etree.cElementTree as ElementTree
import gzip

# 3rd party modules

# Our modules
import vespa.common.constants as constants
import vespa.common.util.xml_ as util_xml
import vespa.common.util.time_ as util_time
import vespa.common.util.misc as util_misc

# The text below is added as an XML comment to each of our export files. 
# It's meant to help anyone who is unfamiliar with Vespa and somehow gets
# ahold of a VIFF file and wants to do something with it.
STANDARD_COMMENT = """
This XML file is in Vespa Interchange File Format (VIFF). You can download
applications that read and write VIFF files and learn more about VIFF here:
http://scion.duhs.duke.edu/vespa/

It was created with Vespa version %s.
""" % util_misc.get_vespa_version()

def export(filename, export_objects, db=None, comment=None, compress=False):
    """Given a list of objects to export, this exports them to filename.
    The objects must support the .deflate() method.
    
    The db parameter should be a valid database connection if the exported 
    objects exist in the database. They will be marked public.
    
    An optional comment can be included in the export. There's no need
    to sanitize the comment for XML meta-characters (e.g. it's OK if the 
    comment contains < and >).
    
    If compress is true, the output is gzipped.
    
    If writing the export fails, this code will raise IOError. Callers need
    to handle that exception.
    """
    # This code is simple, but not too smart. It builds the entire 
    # ElementTree in memory (and possibly turns it into a string) before
    # sending it to disk. Since each object's export is independent of the
    # others, I could export them one at a time and write each to disk as 
    # I export it, thus having just one object's tree in memory at any
    # given time rather than all of them.
    # That would be more complex, though, and I'd have to do things like 
    # writing the root element, timestamp and comment "by hand" (i.e. without
    # the support of ElementTree). 
    # Following Knuth's mantra of premature optimization being the root of 
    # evil, I'll stick with the simple and inefficient method until it 
    # becomes a problem.
    root = ElementTree.Element(constants.Export.ROOT_ELEMENT_NAME,
                               { "version" : constants.Export.VERSION })

    # We append an XML comment that we hope is informative. This is the same
    # for every exported file. Don't confuse it with the comment that the
    # user supplies.
    root.append(ElementTree.Comment(STANDARD_COMMENT))
                               
    util_xml.TextSubElement(root, "timestamp", util_time.now().isoformat())
    
    if comment:
        util_xml.TextSubElement(root, "comment", comment)
    else:
        # Add the element but leave it empty. 
        ElementTree.SubElement(root, "comment")
        
    for export_object in export_objects:
        element = export_object.deflate(constants.Deflate.ETREE)
        if not hasattr(element,'__iter__'):
            root.append(element)
        else:
            for item in element:
                try:
                    root.append(item)
                except:
                    bob = 10
    
    # Prettify the XML
    util_xml.indent(root)
    
    if compress:
        # Compression is done with gzip/zlib 
        xml = ElementTree.tostring(root, "utf-8")
        gzip_file = gzip.GzipFile(filename, "wb")
        gzip_file.write(xml)
        gzip_file.close()
    else:
        tree = ElementTree.ElementTree(root)
        tree.write(filename, "utf-8")

    if db:
        # Mark these as public now that they've been exported. 
        db.mark_public(export_objects)

