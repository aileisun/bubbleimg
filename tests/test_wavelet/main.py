# main.py
# ALS 2016 04 13

import noiselevel
reload(noiselevel)
import waveletdenoise
reload(waveletdenoise)
import measure
reload(measure)

def main():
    """
    auto do everything
    """


    # measure and plot noise level
    noiselevel.main()
    waveletdenoise.main()
    # measure.main()

