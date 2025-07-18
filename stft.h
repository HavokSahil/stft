/**
 * @file stft.h
 * @brief Short-Time Fourier Transform (STFT) implementation using PFFFT.
 *
 * This header defines functions for computing the STFT of a real-valued signal
 * using a configurable window and FFT size, optimized with the PFFFT library.
 *
 * @author Sahil
 * @license MIT
 */

/**
 * MIT License
 *
 * Copyright (c) 2025 Sahil
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef STFT_H_
#define STFT_H_

#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>

#define WTYPE float
#include "window-bank.h"

#ifndef CUSTOM_PFFFT_PATH
#include "pffft/pffft.h"
#else
#include CUSTOM_PFFFT_PATH
#endif

/** @brief Represents a complex number (packed as 2 floats) */
typedef struct __attribute__((__packed__)) {
    float re; /**< Real part */
    float im; /**< Imaginary part */
} Complex_t;

_Static_assert(sizeof(Complex_t) == sizeof(float) * 2, "Complex_t not packed");

#if defined(STFT_DEBUG)
#include <assert.h>
#include <stdio.h>
#define STFT_ERR(...) fprintf(stderr, __VA_ARGS__)
#else
#define STFT_ERR(...)
#endif

/** @brief STFT computation mode */
typedef enum sftft_mode {
    STFT_SLIDING, /**< Placeholder for sliding window STFT (not implemented) */
    STFT_FFT      /**< Simple windowed FFT mode */
} STFT_Mode_t;

/**
 * @brief Configuration structure for STFT processing.
 */
typedef struct stft_config {
    size_t hop;         /**< Hop size between frames */
    size_t win;         /**< Window size */
    size_t insize;      /**< Total input signal size */
    size_t outsize;     /**< Number of frames = (insize - win) / hop + 1 */
    size_t fftsize;     /**< FFT size (next power of 2 ≥ win) */
    PFFFT_Setup *setup; /**< Internal PFFFT setup object */
    float *input;       /**< Aligned input buffer for FFT */
    float *work;        /**< Aligned scratch buffer for FFT */
    float *output;      /**< Aligned FFT output buffer (interleaved) */
    Window *pwin;       /**< Window coefficients */
    STFT_Mode_t mode;   /**< Processing mode */
} STFT_Config_t;

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Return the next power of two ≥ n. */
unsigned long nextPow2(unsigned long n);

/** @brief Return the absolute value (magnitude) of a complex number. */
float cabs(Complex_t a);

/** @brief Return the sum of two complex numbers. */
Complex_t csum(Complex_t a, Complex_t b);

/** @brief Return the product of two complex numbers. */
Complex_t cprod(Complex_t a, Complex_t b);

/**
 * @brief Initialize the STFT configuration.
 *
 * @param hop Hop size between frames.
 * @param win Window size.
 * @param insize Length of the input signal.
 * @param wintype Type of window to apply (e.g., Hamming).
 * @param mode STFT mode to use.
 * @return Pointer to initialized STFT config, or NULL on failure.
 */
STFT_Config_t *stft_config_init(size_t hop, size_t win, size_t insize,
                                WinType wintype, STFT_Mode_t mode);

/**
 * @brief Deallocate the STFT configuration and associated buffers.
 *
 * @param config Pointer to STFT configuration.
 */
void stft_config_deinit(STFT_Config_t *config);

/**
 * @brief Compute STFT of the input signal.
 *
 * The output is written into a preallocated array of shape:
 * output[outsize][fftsize / 2], where each row contains complex frequency bins.
 *
 * @param config Initialized STFT configuration.
 * @param input Real input signal (length must be ≥ config->insize).
 * @param output 2D array of Complex_t to hold the results.
 * @return 0 on success, -1 on error.
 */
int stft_compute(STFT_Config_t *config, float *input, Complex_t **output);

#ifdef __cplusplus
}
#endif

/****************************/
/*  Implementation Section  */
/****************************/

STFT_Config_t *stft_config_init(size_t hop, size_t win, size_t insize,
                                WinType wintype, STFT_Mode_t mode) {
    if (win > insize) {
        STFT_ERR("input size too small.\n");
        return NULL;
    }

    if (hop <= 0) {
        STFT_ERR("invalid hop size.\n");
        return NULL;
    }

    size_t outsize = (insize - win) / hop + 1;
    if (outsize <= 0) {
        STFT_ERR("Invalid output size %zu.\n", outsize);
        return NULL;
    }

    STFT_Config_t *config = (STFT_Config_t *)malloc(sizeof(STFT_Config_t));
    if (!config) {
        STFT_ERR("Allocation failed for STFT_Config_t.\n");
        return NULL;
    }

    size_t fftsize = nextPow2(win);
    config->hop = hop;
    config->win = win;
    config->insize = insize;
    config->outsize = outsize;
    config->mode = mode;
    config->fftsize = fftsize;

    config->pwin = window_init(fftsize, wintype);
    if (!config->pwin) {
        STFT_ERR("Failed to allocate the window.\n");
        goto cleanup;
    }

    config->input = (float *)pffft_aligned_malloc(sizeof(float) * fftsize);
    config->work = (float *)pffft_aligned_malloc(sizeof(float) * fftsize);
    config->output = (float *)pffft_aligned_malloc(sizeof(float) * fftsize);
    if (!config->input || !config->work || !config->output) {
        STFT_ERR("Failed to allocate FFT buffers.\n");
        goto cleanup;
    }

    config->setup = pffft_new_setup(fftsize, PFFFT_REAL);
    if (!config->setup) {
        STFT_ERR("Failed to initialize PFFFT setup.\n");
        goto cleanup;
    }

    return config;

cleanup:
    stft_config_deinit(config);
    return NULL;
}

void stft_config_deinit(STFT_Config_t *config) {
    if (!config)
        return;
    if (config->setup)
        pffft_destroy_setup(config->setup);
    if (config->input)
        pffft_aligned_free(config->input);
    if (config->work)
        pffft_aligned_free(config->work);
    if (config->output)
        pffft_aligned_free(config->output);
    if (config->pwin)
        window_deinit(config->pwin);
    free(config);
}

int stft_compute(STFT_Config_t *config, float *input, Complex_t **output) {
    if (!config || !input || !output) {
        STFT_ERR("Null pointer passed to stft_compute.\n");
        return -1;
    }

#if defined(STFT_DEBUG)
    assert(config->pwin && config->input && config->work && config->output &&
           config->setup);
#endif

    if (config->mode == STFT_FFT) {
        for (size_t i = 0; i < config->outsize; ++i) {
            memset(config->input, 0, sizeof(float) * config->fftsize);
            memcpy(config->input, input + i * config->hop,
                   config->win * sizeof(float));

            for (size_t k = 0; k < config->fftsize; ++k) {
                config->input[k] *= config->pwin->values[k];
            }

            pffft_transform_ordered(config->setup, config->input,
                                    config->output, config->work,
                                    PFFFT_FORWARD);

            memcpy(output[i], config->output, sizeof(float) * config->fftsize);
        }
        return 0;
    }

    return -1;
}

unsigned long nextPow2(unsigned long n) {
    if (n == 0)
        return 1;
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
#if ULONG_MAX > 0xFFFFFFFF
    n |= n >> 32;
#endif
    return n + 1;
}

float cabs(Complex_t a) { return sqrtf(a.re * a.re + a.im * a.im); }

Complex_t csum(Complex_t a, Complex_t b) {
    return (Complex_t){a.re + b.re, a.im + b.im};
}

Complex_t cprod(Complex_t a, Complex_t b) {
    return (Complex_t){a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
}

#if defined(STFT_DEBUG)
/** @brief Debug print of complex number. */
void print_complex(Complex_t a) { printf("%.3f + %.3fi\n", a.re, a.im); }
#endif

#endif // STFT_H_
