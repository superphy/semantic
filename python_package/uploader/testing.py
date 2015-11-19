__author__ = 'Stephen Kan'

import hashlib


m = hashlib.md5()
m.update("Nobody inspects")
m.update(" the spammish repetition")
print m.hexdigest()

n = hashlib.md5()
n.update("Nobody inspects the spammish repetition")
print n.hexdigest()
