import numpy as np
import matplotlib.pyplot as plt
import math

def primitive(y):

    tmp=-1.1547*math.atan(0.57735-1.1547*math.sqrt(y))-0.66*math.log(math.sqrt(y)+1)+0.33*math.log(y-math.sqrt(y)+1)
    return tmp


E0_N = 500000.
E0_P = 85000.
Ec = 117000.

print primitive(0)
print primitive(E0_N/Ec)

#power to ions
pi_N = Ec/E0_N*(primitive(E0_N/Ec)-primitive(0))
pi_P = Ec/E0_P*(primitive(E0_P/Ec)-primitive(0))

#power to electrons
pe_N = 1 - pi_N
pe_P = 1 - pi_P

print "NEGATIVE"
print "IONS:", pi_N, " EL:", pe_N

print "POSITIVE"
print "IONS:", pi_P, " EL:", pe_P
