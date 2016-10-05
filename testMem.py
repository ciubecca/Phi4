import numpy
import sys

l1 = [0]*1000
l2 = [numpy.byte(0)]*1000

print(sys.getsizeof(l1))
print(sys.getsizeof(l2))
