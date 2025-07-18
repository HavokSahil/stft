
#   STFT: Short-Time Fourier Transform using PFFFT in C

Author: Sahil Raj \
License: MIT \
Language: C (C99) \
FFT Backend: [PFFFT](https://github.com/marton78/pffft)

------------------------------------------------------------
## 1. OVERVIEW
------------------------------------------------------------
This project implements a real-valued Short-Time Fourier Transform (STFT)
using the PFFFT library. It supports configurable window types and sizes,
and outputs interleaved complex FFT data for further processing or plotting.

------------------------------------------------------------
## 2. BUILD INSTRUCTIONS
------------------------------------------------------------
To build all targets:

    $ make

To clean:

    $ make clean

------------------------------------------------------------
## 3. RUN INSTRUCTIONS
------------------------------------------------------------
To run the test STFT computation:

    $ ./test_stft

To run the plotting test (requires Python):

```shell
$ ./test_plot
$ python3 stft_plot.py
```

Output files:
    - signal.txt     : Time-domain signal (text)
    - stft_out.bin   : STFT output (binary float32 complex matrix)

------------------------------------------------------------
## 4. FILE DESCRIPTIONS
------------------------------------------------------------
| filename          | description                                             |
|-------------------|---------------------------------------------------------|
| stft.h            | STFT core header + implementation (single-header style) |
| test\_stft.c      | Unit test for STFT with synthetic signals               |
| stft\_test_plot.c | Test + signal generator + binary output dumper          |
| stft\_plot.py     | Python script to plot spectrogram from STFT output      |
| pffft/            | FFT backend (submodule or manually copied)              |

------------------------------------------------------------
## 5. PLOTTING OUTPUT
------------------------------------------------------------
To visualize STFT output:

1. Build and run the test_plot:\
```bash
$ make test_plot
$ ./test_plot
```

2. Run the Python script:
       $ python3 stft_plot.py

The script uses `matplotlib` and shows:
    - The input time-domain signal
    - The STFT magnitude spectrogram in dB

------------------------------------------------------------
## 6. DEPENDENCIES
------------------------------------------------------------
### C:
    - PFFFT (header + source in ./pffft/)
    - -lm (standard maths library)

### Python (for plotting):
    - numpy
    - matplotlib

Install with:

    $ pip install numpy matplotlib

------------------------------------------------------------
## 7. LICENSE
------------------------------------------------------------
This project is released under the MIT License. 
See header of stft.h for details.
