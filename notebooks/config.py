import sys
import os

exeext = '.exe'
platform = sys.platform
if platform.lower() == 'darwin':
    exeext = '.mac'
    
mfexe = os.path.abspath(os.path.join('..', 'bin', 'mf2005{}'.format(exeext)))
mpexe = os.path.abspath(os.path.join('..', 'bin', 'mp6{}'.format(exeext)))
mtexe = os.path.abspath(os.path.join('..', 'bin', 'mt3dms{}'.format(exeext)))
swexe = os.path.abspath(os.path.join('..', 'bin', 'swt_v4{}'.format(exeext)))


