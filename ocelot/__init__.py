from services import *
from converters import *
from features import *
from kernels import *
from go import *
from experiments import *

__version__ = "0.0-dontuse"

import configparser

_DEFAULT_INI = """\
[virtuoso]
ini = /home/virtuoso/virtuoso.ini
virtuoso = virtuoso-t
isql = isql
endpoint = http://localhost:8890/sparql
default_graph = http://ocelot
"""

# XXX add support for the xdg config path
_INI_PATHS = ["ocelot.ini"]

config = configparser.ConfigParser()
for path in _INI_PATHS:
    config.read(path)
    if len(config.sections()) > 0:
        print "Loaded config from '{}'".format(path)
        break

if len(config.sections()) == 0:
    print "Writing default config to '{}'".format(path)
    with open(_INI_PATHS[0], "wb") as fp:
        fp.writelines(_DEFAULT_INI)
    config.read(_INI_PATHS[0])
