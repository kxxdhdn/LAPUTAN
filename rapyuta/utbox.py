#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Utility box

    Error:
        InputError,
    strike, streplace_scl, streplace, tcolor, print_text
    merge_aliases, is_different_value, codefold

    rapyroot, term_width

"""

import os, re, sys
from pathlib import Path
import numpy as np
import random
from colorama import Fore, Back, Style
from IPython.display import HTML, display

global rapyroot, term_width

## https://stackoverflow.com/questions/15411967/how-can-i-check-if-code-is-executed-in-the-ipython-notebook
def type_shell():
    '''
    Return interactive/terminal shell type
    
    '''
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            # Jupyter notebook or qtconsole
            return 'jupyter'
        elif shell == 'TerminalInteractiveShell':
            # Terminal running IPython
            return 'ipython'
        else:
            # Other type (?)
            return 'others'
    except NameError:
        # Probably standard Python interpreter
        return 'terminal'
        
## rapyroot
try:
    rapyroot = Path(__file__).parent.absolute()
except NameError:
    rapyroot = Path().absolute()
# sys.path.insert(0, rapyroot)

## term_width
if type_shell()=='terminal':
    term_width = int((os.popen('stty size', 'r').read().split())[1])
else:
    term_width = 80


##------------------------------------------------
##
##             <Error> based classes
##
##------------------------------------------------

class Error(Exception):
    '''
    Base class for exceptions in this module

    '''
    pass

class InputError(Error):
    '''
    Exception raised for errors in the input.

    ------ INPUT ------
    expression          input expression in which the error occurred
    message             explanation of the error
    ------ OUTPUT ------
    '''
    def __init__(self, expression, message):
        self.expression = expression
        self.message = message


##------------------------------------------------
##
##             French strike module
##
##------------------------------------------------

def strike(proc, why, cat=None):
    '''
    Raise error message

    Copyright: F. Galliano
    '''
    f = open(rapyroot / 'lib/logo/gameover.txt', 'r')
    
    lines = streplace(f.readlines(),[('\n','')])
    # text = f.readlines()
    # lines = text.copy()
    # subs = [('\n','')]
    # for i in range(len(lines)):
    #     pattern = '|'.join('(%s)' % re.escape(p) for p, s in subs)
    #     substs = [s for p, s in subs]
    #     replace = lambda m: substs[m.lastindex - 1]
    #     lines[i] = re.sub(pattern, replace, lines[i])

    lines.append("")
    lines.append("")
    lines.append("                           ***   INSERT COIN   ***")
    f.close()
    width = 78
    pad = int(np.floor(term_width - width)/2)

    colors = ['Black', 'Red', 'Green', 'Yellow',
              'Blue', 'Magenta', 'Cyan', 'White']
    irand1, irand2 = random.sample(range(len(colors)), 2)
    
    print_text(lines,color=colors[irand1],back_color=colors[irand2],
               weight='Bold',bottom=1,pad=pad,width=width,blink=True)
    print('')

    if cat=='ValueError':
        raise ValueError('<'+proc+'> went on strike: '+why)
    elif cat=='InputError':
        raise InputError('<'+proc+'>', why)
    else:
        raise Error('<'+proc+'> went on strike: '+why)

def streplace_scl(text,subs):
    '''
    String replacement (scalar)

    Copyright: F. Galliano
    '''
    pattern = '|'.join('(%s)' % re.escape(p) for p, s in subs)
    substs = [s for p, s in subs]
    replace = lambda m: substs[m.lastindex - 1]
    return(re.sub(pattern, replace, text))

def streplace(strarr0,subs):
    '''
    String replacement (1D)

    replaced = streplace(textlist,[('a','b'),('b','c'),...])
    Simultaneously perform all substitutions on the subject string.

    Copyright: F. Galliano
    '''
    if (np.isscalar(strarr0)):
        strarr = streplace_scl(strarr0,subs)
    else:
        strarr = strarr0.copy()
        for i in np.arange(len(strarr)):
            strarr[i] = streplace_scl(strarr[i],subs)
    return(strarr)

def tcolor(text,color=None,back_color=None,weight=None,
           reset=True,nocolor=False):
    '''
    Xterm colors

    Change the colors and weight of text printed in the terminal. The settings
    apply only to the text in parenthesis.
    color = 'Black', 'Red', 'Green', 'Yellow', 'Blue', 'Magenta', 'Cyan', 'White'
    back_color = 'Black', 'Red', 'Green', 'Yellow', 'Blue', 'Magenta', 'Cyan', 
                  'White'
    weight = 'Dim', 'Normal', 'Bold'

    Copyright: F. Galliano
    '''
    if (nocolor):
        return(text)
    colist = ['BLACK','RED','GREEN','YELLOW','BLUE','MAGENTA','CYAN','WHITE',
              'LIGHTBLACK_EX','LIGHTRED_EX','LIGHTGREEN_EX','LIGHTYELLOW_EX', \
              'LIGHTBLUE_EX','LIGHTMAGENTA_EX','LIGHTCYAN_EX','LIGHTWHITE_EX']
    weilist = ['DIM','NORMAL','BOLD','BRIGHT']
    if (not isinstance(color,str)): color = ""
    if (not isinstance(back_color,str)): back_color = ""
    if (not isinstance(weight,str)): weight = ""
    if (color.upper() in colist):
        line = eval("Fore."+color.upper())
    elif ("".join(color.upper().split(" "))+"_EX" in colist):
        line = eval("Fore."+"".join(color.upper().split(" "))+"_EX")
    else:
        line = ""
    if (back_color.upper() in colist):
        back = eval('Back.'+back_color.upper())
    elif ("".join(back_color.upper().split(" "))+"_EX" in colist):
        back = eval("Back."+"".join(back_color.upper().split(" "))+"_EX")
    else:
        back = ""
    if (weight.upper() in weilist):
        if (weight.upper() == "BOLD"): weight = "BRIGHT"
        wei = eval('Style.'+weight.upper())
    else:
        wei = ""
    text = wei + back + line + text
    if (reset): text += Style.RESET_ALL
    return(text)

def print_text(text,width=term_width,pad=0,margin=1,center=False,
               top=0,bottom=0,color=None,back_color=None,weight=None,
               blink=False,reset=True,fill=False,charpad=' ',nocolor=False,
               filobj=None,frame=False):
    '''
    Stylish printing
    General function to print text with a nice format

    Break-up long text over several lines. Print with colors and a margin.
    The text is aligned left, by default. If center=True, it is centered.
    If the text is a list of strings, each will be printed on a different line.
    If the background is not colored, margin and pad are degenerate. However,
    if back_color is not None, pad is the number of characters left blank on the
    left side of the box, while margin is the number of characters left blank 
    inside the colored box:

     pad                   width 
    <----><----------------------------------------------->
          _________________________________________________ 
          |                                               |  top
          |   Blabla blabla blabla blabla blabla blabla   |
          |   blabla blabla blabla blabla blabla.         |
          |                                               |  bottom 
          ________________________________________________
          <-->                                         <-->
         margin                                       margin
    
    Copyright: F. Galliano
    '''

    # Arguments
    if (center):
        lbox = width - 2*margin
        pad = int(np.floor((term_width-width)/2))
    else:
        lbox = width - 2*margin
    if (blink):
        pref = '\033[5m'
    else:
        pref = ''
    if (fill):
        fillsym = "·"
    else:
        fillsym = " "
    if (frame):
        horiframe = horisym
    else:
        horiframe = " "

    # Top
    try:
        for i in np.arange(top):
            if (i == 0):
                filobj.write(' '*pad+horiframe*width+'\n')
            else:
                filobj.write(' '*pad+" "*width+'\n')
    except:
        for i in np.arange(top):
            if (i == 0):
                print(' '*pad+ tcolor(horiframe*width,color=color, \
                                      back_color=back_color,weight=weight, \
                                      reset=reset,nocolor=nocolor))
            else:
                print(' '*pad+ tcolor(" "*width-2,color=color, \
                                      back_color=back_color,weight=weight, \
                                      reset=reset,nocolor=nocolor))

    # Chop the text
    if (np.isscalar(text)):
        text = [text]
    Nline = np.size(text)
    for i in np.arange(Nline):
        textlines = text[i]
        l = len(textlines)
        if (l > lbox):
            i0 = 0
            i1 = lbox
            while True:
                if (i1 < l):
                    ibreak = (textlines[i0:i1]).rfind(' ')
                else:
                    ibreak = i1 - i0
                line = textlines[i0:i0+ibreak]
                ll = len(line)
                if (center):
                    i2 = int(( lbox - ll ) / 2)
                else:
                    i2 = 0
                i3 = lbox - ll - i2
                padright = term_width - pad - width
                try:
                    filobj.write(charpad*pad+' '*(i2+margin) +pref+line \
                                + ' '*(i3+margin)+charpad*padright+'\n')
                except:
                    print(tcolor(charpad*pad,color=color,weight=weight, \
                                 nocolor=nocolor) \
                          + tcolor(' '*(i2+margin) +pref+line+ ' '*(i3+margin), \
                                   color=color,back_color=back_color, \
                                   weight=weight,reset=reset,nocolor=nocolor) \
                          + tcolor(charpad*padright,color=color,weight=weight, \
                                   nocolor=nocolor))
                if (i1 >= l): break
                i0 = i0 + ibreak + 1
                i1 = i0 + lbox
                if (i1 > l): i1 = l
        else:
            if (center):
                i2 = int(( lbox - l ) / 2)
            else:
                i2 = 0
            i3 = lbox - l - i2
            padright = term_width - pad - width
            if (i3 > 1):
                fills = " "+fillsym*(i3-1)
            else:
                fills = " "
            try:
                filobj.write(charpad*pad+' '*(i2+margin) + pref + textlines \
                            + fills + ' '*margin + charpad*padright + '\n')
            except:
                print(tcolor(charpad*pad,color=color,weight=weight, \
                             nocolor=nocolor) \
                      + tcolor(' '*(i2+margin) + pref+textlines \
                               + fills + ' '*margin, \
                               color=color,back_color=back_color,weight=weight, \
                               reset=reset,nocolor=nocolor) \
                      + tcolor(charpad*padright,color=color,weight=weight, \
                               nocolor=nocolor))

    # Bottom
    try:
        for i in np.arange(bottom):
            if (i == bottom-1):
                filobj.write(' '*pad+horiframe*width+'\n')
            else:
                filobj.write(' '*pad+" "*width+'\n')
    except:
        for i in np.arange(bottom):
            if (i == bottom-1):
                print(' '*pad \
                      + tcolor(horiframe*width,color=color, \
                               back_color=back_color,weight=weight, \
                               reset=reset,nocolor=nocolor))
            else:
                print(' '*pad \
                      + tcolor(" "*width,color=color,back_color=back_color, \
                               weight=weight,reset=reset,nocolor=nocolor))

##------------------------------------------------
##
##              Unsorted functions
##
##------------------------------------------------

## https://python-forum.io/thread-22874.html    
def merge_aliases(default, **kwargs):
    '''
    Merge aliases

    ------ INPUT ------
    default             default value of the keyword
                          should be the same with default value of each alias
    ------ OUTPUT ------
    return unique input value of target keyword or default value if no input
    '''
    d = {k: v for k, v in kwargs.items() if v is not default}
    if d:
        if len(d) > 1:
            raise InputError(
                'Only one of the following parameters can be set:',
                list(d.keys()))
        else:
            return d.popitem()[1]
    else:
        return default

## https://stackoverflow.com/questions/6346492/how-to-stop-a-for-loop
def is_different_value(iterable, i, j):
    '''
    The guru way to stop a for loop

    '''
    # value = iterable[i][j]
    # return  any(any((cell != value for cell in col)) for col in iterable)
    pass

## [Toggle] http://blog.nextgenetics.net/?e=102
## https://stackoverflow.com/questions/33159518/collapse-cell-in-jupyter-notebook
## [Not worked] https://github.com/Fmajor/grizli_extras/blob/main/grizli_extras/tools.py
def codefold(verbose=False, button=None):
    '''
    Toggle all input cells for the interactive display
    Counterpart of Jupyter menu : Cell > All Output > Toggle

    Copyright: Damian Kao
    '''
    script_html = '''
    <script>code_show=true; 
    function code_toggle() {
        if (code_show){
        $('div.input').hide();
        } else {
        $('div.input').show();
        }
        code_show = !code_show
    } 
    $( document ).ready(code_toggle);
    </script> '''
    
    if verbose:
        script_html += '''
        <p>The raw code for this IPython notebook is by default hidden for easier reading. 
        To toggle on/off the raw code, click the button below.
        <p></p> '''
        
    if button=='onclick':
        script_html += '''
        <button type="button" onclick="code_toggle(document.getElementById('hide'))">
        Show/hide inputs
        </button> '''
    else:
        script_html += '''
        <a href="javascript:code_toggle()">Show/hide inputs
        </a> '''

    display(HTML(script_html))
