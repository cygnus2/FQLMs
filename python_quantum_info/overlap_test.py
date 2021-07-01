import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from copy import copy

import sys
sys.path.append('../python_data/')
from ed_result_wrapper import EDResult

# Read data.
param = {
    'L' : [2,2,2],
    'gauge_particles' : 'bosons',
    'lambda' : -100,
    'winding_sector' : [0,0,0]
}
result = EDResult(param, datadir="../parity_troubleshooting")
es = result.get_eigenstates(param['lambda'])
