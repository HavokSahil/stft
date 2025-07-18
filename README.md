============================================================
   STFT: Short-Time Fourier Transform using PFFFT in C
============================================================

Author: Sahil Raj \
License: MIT \
Language: C (C99) \
FFT Backend: [PFFFT](https://github.com/marton78/pffft)

------------------------------------------------------------
CONTENTS
------------------------------------------------------------
1. Overview
2. Build Instructions
3. Run Instructions
4. File Descriptions
5. Plotting Output
6. Dependencies
7. License

------------------------------------------------------------
1. OVERVIEW
------------------------------------------------------------
This project implements a real-valued Short-Time Fourier Transform (STFT)
using the PFFFT library. It supports configurable window types and sizes,
and outputs interleaved complex FFT data for further processing or plotting.

------------------------------------------------------------
2. BUILD INSTRUCTIONS
------------------------------------------------------------
To build all targets:

    $ make

To clean:

    $ make clean

------------------------------------------------------------
3. RUN INSTRUCTIONS
------------------------------------------------------------
To run the test STFT computation:

    $ ./test_stft

To run the plotting test (requires Python):

    $ ./test_plot
    $ python3 stft_plot.py

Output files:
    - signal.txt     : Time-domain signal (text)
    - stft_out.bin   : STFT output (binary float32 complex matrix)

------------------------------------------------------------
4. FILE DESCRIPTIONS
------------------------------------------------------------
- stft.h              : STFT core header + implementation (single-header style)
- test_stft.c         : Unit test for STFT with synthetic signals
- stft_test_plot.c    : Test + signal generator + binary output dumper
- stft_plot.py        : Python script to plot spectrogram from STFT output
- Makefile            : Build system
- pffft/              : FFT backend (submodule or manually copied)

------------------------------------------------------------
5. PLOTTING OUTPUT
------------------------------------------------------------
To visualize STFT output:

1. Build and run the test_plot:
       $ make test_plot
       $ ./test_plot

2. Run the Python script:
       $ python3 stft_plot.py

The script uses `matplotlib` and shows:
    - The input time-domain signal
    - The STFT magnitude spectrogram in dB

------------------------------------------------------------
6. DEPENDENCIES
------------------------------------------------------------
C:
    - PFFFT (header + source in ./pffft/)
    - -lm (standard maths library)

Python (for plotting):
    - numpy
    - matplotlib

Install with:

    $ pip install numpy matplotlib

------------------------------------------------------------
7. LICENSE
------------------------------------------------------------
This project is released under the MIT License.  
See header of stft.h for details.

============================================================
