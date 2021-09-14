#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Utilities

    Error
    InputError(Error)

"""

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class InputError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message
 
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
    
## Tests
if __name__ == '__main__':

    ## Test merge_aliases
    ##--------------------
    def putfoo(bar, position=None, pos=None):
        pos = merge_aliases(None, position=position, pos=pos)
        print('Position is:', repr(pos))
    
    putfoo('hello', position=7)
    putfoo('hello', pos=9)
    putfoo('hello', position= 11, pos=13)
