import os
import sys
from pylal import SnglInspiralUtils


def readCoinc(CoincFile):
    coinc = SnglInspiralUtils.ReadSnglInspiralFromFiles(CoincFile)
    for coinc_index, coinc_row in enumerate(coinc):
        mass1 = coinc_row.mass1
        mass2 = coinc_row.mass2
        chi1 = coinc_row.spin1z
        time = coinc_row.end_time
    return [mass1, mass2, chi1, time]


