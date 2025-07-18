/**
 * @file stft_test_plot.c
 * @brief Test program for STFT implementation with plotting capabilities
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

// Include your STFT header
#include "stft.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @brief Generate a test signal (chirp from f0 to f1)
 */
void generate_chirp(float *signal, size_t length, float sample_rate, float f0,
                    float f1, float duration) {
    for (size_t i = 0; i < length; i++) {
        float t = (float)i / sample_rate;
        float freq = f0 + (f1 - f0) * t / duration;
        signal[i] = sin(2.0 * M_PI * freq * t);
    }
}

/**
 * @brief Generate a multi-tone signal
 */
void generate_multitone(float *signal, size_t length, float sample_rate,
                        float *frequencies, float *amplitudes, int num_tones) {
    // Initialize to zero
    for (size_t i = 0; i < length; i++) {
        signal[i] = 0.0f;
    }

    // Add each tone
    for (int tone = 0; tone < num_tones; tone++) {
        for (size_t i = 0; i < length; i++) {
            float t = (float)i / sample_rate;
            signal[i] +=
                amplitudes[tone] * sin(2.0 * M_PI * frequencies[tone] * t);
        }
    }
}

/**
 * @brief Generate a noisy sine wave
 */
void generate_noisy_sine(float *signal, size_t length, float sample_rate,
                         float frequency, float noise_level) {
    for (size_t i = 0; i < length; i++) {
        float t = (float)i / sample_rate;
        float noise = noise_level * ((float)rand() / RAND_MAX - 0.5f) * 2.0f;
        signal[i] = sin(2.0 * M_PI * frequency * t) + noise;
    }
}

/**
 * @brief Save the original signal to a file
 */
void save_signal(const char *filename, float *signal, size_t length,
                 float sample_rate) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error opening file %s\n", filename);
        return;
    }

    fprintf(fp, "# Time Signal\n");
    for (size_t i = 0; i < length; i++) {
        float t = (float)i / sample_rate;
        fprintf(fp, "%.6f %.6f\n", t, signal[i]);
    }
    fclose(fp);
    printf("Signal saved to %s\n", filename);
}

void dump_stft_output(const char *filename, Complex_t **output, size_t rows,
                      size_t cols) {
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        fprintf(stderr, "Failed to open %s for writing.\n", filename);
        return;
    }

    for (size_t i = 0; i < rows; i++) {
        fwrite(output[i], sizeof(Complex_t), cols, fp);
    }

    fclose(fp);
}

/**
 * @brief Allocate 2D array for STFT output
 */
Complex_t **allocate_stft_output(size_t outsize, size_t fftsize) {
    Complex_t **output = (Complex_t **)malloc(outsize * sizeof(Complex_t *));
    if (!output)
        return NULL;

    for (size_t i = 0; i < outsize; i++) {
        output[i] = (Complex_t *)malloc((fftsize / 2) * sizeof(Complex_t));
        if (!output[i]) {
            // Clean up on failure
            for (size_t j = 0; j < i; j++) {
                free(output[j]);
            }
            free(output);
            return NULL;
        }
    }
    return output;
}

/**
 * @brief Free 2D array for STFT output
 */
void free_stft_output(Complex_t **output, size_t outsize) {
    if (!output)
        return;
    for (size_t i = 0; i < outsize; i++) {
        free(output[i]);
    }
    free(output);
}

int main() {
    const float sample_rate = 8000.0f; // 8 kHz
    const float duration = 1.0f;       // 1 second
    const size_t signal_len = (size_t)(sample_rate * duration);

    // Allocate signal
    float *signal = (float *)calloc(signal_len, sizeof(float));
    if (!signal) {
        fprintf(stderr, "Failed to allocate signal.\n");
        return EXIT_FAILURE;
    }

    // Simple chirp
    // generate_chirp(signal, signal_len, sample_rate, 100.0f, 3000.0f,
    // duration);

    // Multi-tone
    // float freqs[] = {440.0f, 880.0f, 1320.0f};
    // float amps[] = {1.0f, 0.5f, 0.3f};
    // generate_multitone(signal, signal_len, sample_rate, freqs, amps, 3);

    // Noisy sine
    generate_noisy_sine(signal, signal_len, sample_rate, 1000.0f, 0.2f);

    // Save raw signal
    save_signal("signal.txt", signal, signal_len, sample_rate);

    // === STFT Params ===
    size_t hop = 128;
    size_t win = 256;
    STFT_Config_t *config =
        stft_config_init(hop, win, signal_len, HAMMING, STFT_FFT);
    if (!config) {
        fprintf(stderr, "STFT config init failed.\n");
        free(signal);
        return EXIT_FAILURE;
    }

    // Allocate output buffer
    Complex_t **output = allocate_stft_output(config->outsize, config->fftsize);
    if (!output) {
        fprintf(stderr, "Failed to allocate STFT output.\n");
        stft_config_deinit(config);
        free(signal);
        return EXIT_FAILURE;
    }

    // Compute STFT
    if (stft_compute(config, signal, output) != 0) {
        fprintf(stderr, "STFT computation failed.\n");
        free_stft_output(output, config->outsize);
        stft_config_deinit(config);
        free(signal);
        return EXIT_FAILURE;
    }

    // Save output
    dump_stft_output("stft_out.bin", output, config->outsize,
                     config->fftsize / 2);

    printf("STFT output saved to stft_out.bin\n");

    // Cleanup
    free_stft_output(output, config->outsize);
    stft_config_deinit(config);
    free(signal);

    return EXIT_SUCCESS;
}
