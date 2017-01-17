from oscillators import *
from statefuncs import *

state = [(0,1),(1,2),(2,3),(4,4),(5,1)]
N = 4

print("occn", occn(state))
print("state", state)
ndPair = (2,3)
ndtotPair = (4,4)

print("1")
print(sorted(gendlistPairs(state, ndPair, ndtotPair, 10)))
print("2")
print(sorted(gendlistPairs2(state, ndPair, ndtotPair, 10)))
