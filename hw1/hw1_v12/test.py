import sys
import struct

with open(sys.argv[1], 'rb') as f:
  byte = f.read(4)
  while byte != b'':
    print(byte, struct.unpack('f', byte)[0])	
    byte = f.read(4)