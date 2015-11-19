__author__ = 'Stephen Kan'

import os
import inspect
import string

def path(filename):
    frame = inspect.stack()
    (frame, filepath, line_number, function_name, lines, index) = frame[1]
    del frame
    return os.path.join(os.path.dirname(filepath), filename)

def only_abecedarian(self, str):
    all = string.maketrans('','')
    nodigs = all.translate(all, string.ascii_letters)
    return str.translate(all, nodigs)

def only_digits(str):
    all = string.maketrans('','')
    nodigs = all.translate(all, string.digits)
    return str.translate(all, nodigs)
