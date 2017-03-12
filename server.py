from pymol import cmd, util
from threading import Thread
from SimpleXMLRPCServer import SimpleXMLRPCServer
import cPickle
import os
import re
import string


def start_pymol_server():
    print "========================= Loading PyMOL XML-RPC extensions ============================="
    server = pymol_interface()
    print "========================================================================================"


class pymol_xmlrpc_server (SimpleXMLRPCServer):
    def __init__(self, addr, interface):
        self.interface = interface
        SimpleXMLRPCServer.__init__(self, addr, logRequests=0)

    def _dispatch(self, method, params):
        if not self.interface.enable_xmlrpc:
            return -1
        result = -1
        func = None
        if hasattr(self.interface, method):
            func = getattr(self.interface, method)
        elif hasattr(cmd, method):
            func = getattr(cmd, method)
        if not callable(func):
            print "%s is not a callable object" % method
        else:
            result = func(*params)
            if result is None:
                result = -1
        return result


class pymol_interface (object):
    def __init__(self):
        self.enable_xmlrpc = True
        # the port can be set via an environment variable - although it could
        # just as easily be passed to __init__
        port = string.atoi(os.environ.get("PYMOL_XMLRPC_PORT", "9123"))
        self._server = pymol_xmlrpc_server(("localhost", port), self)
        t = threading.Thread(target=self._server.serve_forever)
        t.setDaemon(1)
        t.start()
        print "Started XML-RPC server on port %d" % port
