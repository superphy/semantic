#!/usr/bin/env python

import argparse
import struct
from binascii import b2a_base64 as e64
from binascii import a2b_base64 as d64
import sys
import scrypt
import Crypto.Random
random = Crypto.Random.new().read
from passlib.utils import consteq

"""Use the scrypt library to generate a secure hash of the user password. \
   Return the hash of the password"""

parser = argparse.ArgumentParser()
parser.add_argument("password", help="The password to be hashed")
parser.add_argument("hash", help="The stored password")
args = parser.parse_args()

_PARAMS = struct.Struct("!BBBB")


def pack_verifier(logN, r, p, salt, hash):
    packed = _PARAMS.pack(logN, r, p, len(salt)) + salt + hash
    return packed


def unpack_verifier(verifier):
    logN, r, p, salt_bytes = _PARAMS.unpack_from(verifier)
    i = _PARAMS.size + salt_bytes
    salt = verifier[_PARAMS.size:i]
    hash = verifier[i:]
    return logN, r, p, salt, hash


def make_verifier( password, logN=14, r=8, p=1, salt_bytes=16, hash_bytes=16):
    salt = random(salt_bytes)
    hash_value = scrypt.hash(password, salt, 1<<logN, r, p, hash_bytes)
    return pack_verifier(logN,r,p,salt,hash_value)


def verify_password( password, verifier ):
    logN, r, p, salt, hash_value = unpack_verifier(verifier)
    newhash = scrypt.hash(password,salt,1<<logN,r,p,len(hash_value))
    return consteq(newhash, hash_value)


if __name__== "__main__":
    v = make_verifier(args.password)
    ev = e64(v).strip()

    # if not(verify_password(args.password, d64(args.hash))):
    #     print("Password not verified!!\n")
    #     sys.exit(1)
   
    print(ev)


   
