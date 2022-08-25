# Licensed under a 3-clause BSD style license - see LICENSE

import sys
import os.path
import numpy as np
import random

## Add rapyuta root to $PATH
rapyroot = os.path.dirname(os.path.abspath(__file__))
# try:
#     rapyroot = Path(__file__).parent.absolute()
# except NameError:
#     rapyroot = Path().absolute()
sys.path.insert(0, rapyroot)

## Local
from .version import __version__, release_date
import utbox as UT

## rapyuta logo
f = open(rapyroot+'/lib/logo/rapyuta.txt', 'r')
lines = UT.streplace(f.readlines(),[('\n','')])
lines.append("                          "
             f"Version {__version__} ({release_date})")
f.close()
width = 78
pad = int(np.floor(UT.term_width - width)/2)
colors = ['Black', 'Red', 'Green', 'Yellow',
          'Blue', 'Magenta', 'Cyan', 'White']
irand1, irand2 = random.sample(range(len(colors)), 2)
# irand1 = random.sample(range(1,len(colors)), 1)[0]
# irand2 = 0 # always black background
UT.print_text(lines,color=colors[irand1],back_color=colors[irand2],
           weight='Bold',bottom=1,pad=pad,width=width,blink=True)
print('')

# print("\n               \\  \\|/  /")
# # print("             \\  \\_ | _/  /")
# print("            \\ _\\_ ||| _/_ /")
# print("                  |||")
# print("              /\\/-/|\\-\\/\\")
# print("          _ | _M_ _|_ _M_ | _")
# print("        / .~.    _ ^ _    .~. \\")
# print("       |_|_0_|  |_&_&_|  |_0_|_|")
# print("      //    .~.    o    .~.    \\\\")
# print("     |_M_ _|_?_|__[_]__|_?_|_ _M_|")
# print("            |             |")
# print("             \\  rapyuta  /")
# print("               \\       /         _")
# print("    * ¨^  .~       =           ~^' _-")
# print("      ~ °\n")
# # print("\n            Author: D. HU")
# print(f"        Version {__version__} ({release_date})")
# print("\n")

def version_control():
    '''
	Version control

    '''

    version = __version__
    if 'dev' in version:
        version = 'latest'
    else:
        version = 'v' + version

    print(version)

def boost_rapyuta():
	'''
	Booster

	'''
	print("\nHouston, Tranquility Base here. The Eagle has landed.\n")

if __name__ == '__main__':
    
    boost_rapyuta()
