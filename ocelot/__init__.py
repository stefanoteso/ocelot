from experiment import *
from services import *
from converters import *
from features import *
from kernels import *
from go import *
from utils import *

__version__ = "0.0-dontuse"

import os, configparser

_DEFAULT_INI = """\
[rdf]
endpoint = http://localhost:8890/sparql
default_graph = http://ocelot

[virtuoso]
ini = /home/virtuoso/virtuoso.ini
virtuoso = virtuoso-t
isql = isql
"""

_ini_paths = ["ocelot.ini"]
try:
    from appdirs import user_config_dir
    _ini_paths.append(os.path.join(user_config_dir(), "ocelot.ini"))
except ImportError:
    pass

config = configparser.ConfigParser()
for path in _ini_paths:
    config.read(path)
    if len(config.sections()) > 0:
        print "Loaded config from '{}'".format(path)
        break

if len(config.sections()) == 0:
    print "Writing default config to '{}'".format(_ini_paths[0])
    with open(_ini_paths[0], "wb") as fp:
        fp.writelines(_DEFAULT_INI)
    config.read(_ini_paths[0])
