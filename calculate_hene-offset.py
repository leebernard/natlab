'''
This is for calculating the offsets between the 632.8 nm alignment laser, and the 1.5 um chief ray of operating beam.
Pertinent surfaces are the camera focus, OAP 2 (step 5), and the final fold mirrors, Ch 1 0.5 in fold and Ch 2 0.5 in
fold (step 7).
'''
import matplotlib.pyplot as plt
import numpy as np

from spot0 import spot
ps = 18e-6  # pixel pitch, meters
f1, f2 = 0.1016, 0.27224  # OAP focal lengths, meters






