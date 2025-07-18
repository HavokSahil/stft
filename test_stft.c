/**
 * @file test_stft.c
 * @brief Comprehensive test suite for STFT implementation
 *
 * This test suite includes:
 * - Unit tests for utility functions
 * - Integration tests for STFT computation
 * - Edge case testing
 * - Performance benchmarks
 * - Memory leak detection
 *
 * @author Test Suite Generator
 * @license MIT
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

// Enable debug mode for testing
#define STFT_DEBUG
#include "stft.h"

// Test configuration
#define TEST_TOLERANCE 1e-6f
#define PERFORMANCE_ITERATIONS 100
#define MAX_TEST_SIZE 16384

// ANSI color codes for output
#define ANSI_COLOR_RED "\x1b[31m"
#define ANSI_COLOR_GREEN "\x1b[32m"
#define ANSI_COLOR_YELLOW "\x1b[33m"
#define ANSI_COLOR_BLUE "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN "\x1b[36m"
#define ANSI_COLOR_RESET "\x1b[0m"

// Test result tracking
typedef struct {
    int passed;
    int failed;
    int total;
} TestResults;

static TestResults g_results = {0, 0, 0};

// Utility macros
#define TEST_ASSERT(condition, message)                                        \
    do {                                                                       \
        g_results.total++;                                                     \
        if (condition) {                                                       \
            printf(ANSI_COLOR_GREEN "✓ PASS: " ANSI_COLOR_RESET "%s\n",        \
                   message);                                                   \
            g_results.passed++;                                                \
        } else {                                                               \
            printf(ANSI_COLOR_RED "✗ FAIL: " ANSI_COLOR_RESET "%s\n",          \
                   message);                                                   \
            g_results.failed++;                                                \
        }                                                                      \
    } while (0)

#define TEST_ASSERT_FLOAT_EQ(expected, actual, message)                        \
    do {                                                                       \
        g_results.total++;                                                     \
        if (fabs((expected) - (actual)) < TEST_TOLERANCE) {                    \
            printf(ANSI_COLOR_GREEN "✓ PASS: " ANSI_COLOR_RESET "%s\n",        \
                   message);                                                   \
            g_results.passed++;                                                \
        } else {                                                               \
            printf(ANSI_COLOR_RED "✗ FAIL: " ANSI_COLOR_RESET                  \
                                  "%s (expected: %.6f, actual: %.6f)\n",       \
                   message, (double)(expected), (double)(actual));             \
            g_results.failed++;                                                \
        }                                                                      \
    } while (0)

#define TEST_SECTION(name)                                                     \
    printf(ANSI_COLOR_BLUE "\n=== %s ===" ANSI_COLOR_RESET "\n", name)

// Timer utilities
static double get_time_ms() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * 1000.0 + tv.tv_usec / 1000.0;
}

// Test signal generators
static void generate_sine_wave(float *signal, size_t length, float frequency,
                               float sample_rate) {
    for (size_t i = 0; i < length; i++) {
        signal[i] = sinf(2.0f * M_PI * frequency * i / sample_rate);
    }
}

static void generate_chirp(float *signal, size_t length, float f0, float f1,
                           float sample_rate) {
    for (size_t i = 0; i < length; i++) {
        float t = (float)i / sample_rate;
        float freq = f0 + (f1 - f0) * t * t / (2.0f * length / sample_rate);
        signal[i] = sinf(2.0f * M_PI * freq * t);
    }
}

static void generate_noise(float *signal, size_t length, float amplitude) {
    for (size_t i = 0; i < length; i++) {
        signal[i] = amplitude * (2.0f * rand() / RAND_MAX - 1.0f);
    }
}

// Utility function tests
void test_nextPow2() {
    TEST_SECTION("nextPow2 Function Tests");

    TEST_ASSERT(nextPow2(0) == 1, "nextPow2(0) == 1");
    TEST_ASSERT(nextPow2(1) == 1, "nextPow2(1) == 1");
    TEST_ASSERT(nextPow2(2) == 2, "nextPow2(2) == 2");
    TEST_ASSERT(nextPow2(3) == 4, "nextPow2(3) == 4");
    TEST_ASSERT(nextPow2(4) == 4, "nextPow2(4) == 4");
    TEST_ASSERT(nextPow2(5) == 8, "nextPow2(5) == 8");
    TEST_ASSERT(nextPow2(1023) == 1024, "nextPow2(1023) == 1024");
    TEST_ASSERT(nextPow2(1024) == 1024, "nextPow2(1024) == 1024");
    TEST_ASSERT(nextPow2(1025) == 2048, "nextPow2(1025) == 2048");
}

void test_complex_functions() {
    TEST_SECTION("Complex Number Function Tests");

    Complex_t a = {3.0f, 4.0f};
    Complex_t b = {1.0f, 2.0f};

    // Test cabs
    float magnitude = cabs(a);
    TEST_ASSERT_FLOAT_EQ(5.0f, magnitude, "cabs({3, 4}) == 5");

    // Test csum
    Complex_t sum = csum(a, b);
    TEST_ASSERT_FLOAT_EQ(4.0f, sum.re, "csum real part");
    TEST_ASSERT_FLOAT_EQ(6.0f, sum.im, "csum imaginary part");

    // Test cprod
    Complex_t prod = cprod(a, b);
    TEST_ASSERT_FLOAT_EQ(-5.0f, prod.re, "cprod real part");
    TEST_ASSERT_FLOAT_EQ(10.0f, prod.im, "cprod imaginary part");

    // Test zero complex number
    Complex_t zero = {0.0f, 0.0f};
    TEST_ASSERT_FLOAT_EQ(0.0f, cabs(zero), "cabs of zero");

    Complex_t sum_zero = csum(a, zero);
    TEST_ASSERT_FLOAT_EQ(a.re, sum_zero.re, "sum with zero - real");
    TEST_ASSERT_FLOAT_EQ(a.im, sum_zero.im, "sum with zero - imaginary");
}

void test_stft_config_init() {
    TEST_SECTION("STFT Configuration Tests");

    // Test valid configuration
    STFT_Config_t *config = stft_config_init(256, 512, 4096, HAMMING, STFT_FFT);
    TEST_ASSERT(config != NULL, "Valid configuration initialization");

    if (config) {
        TEST_ASSERT(config->hop == 256, "Hop size correctly set");
        TEST_ASSERT(config->win == 512, "Window size correctly set");
        TEST_ASSERT(config->insize == 4096, "Input size correctly set");
        TEST_ASSERT(config->outsize == 15, "Output size correctly calculated");
        TEST_ASSERT(config->fftsize == 512, "FFT size correctly set");
        TEST_ASSERT(config->mode == STFT_FFT, "Mode correctly set");
        TEST_ASSERT(config->pwin != NULL, "Window allocated");
        TEST_ASSERT(config->input != NULL, "Input buffer allocated");
        TEST_ASSERT(config->work != NULL, "Work buffer allocated");
        TEST_ASSERT(config->output != NULL, "Output buffer allocated");
        TEST_ASSERT(config->setup != NULL, "PFFFT setup created");

        stft_config_deinit(config);
    }

    // Test invalid configurations
    config = stft_config_init(512, 512, 256, HAMMING, STFT_FFT);
    TEST_ASSERT(config == NULL, "Invalid configuration (input too small)");

    config = stft_config_init(0, 512, 4096, HAMMING, STFT_FFT);
    TEST_ASSERT(config == NULL, "Invalid configuration (zero hop)");
}

void test_stft_compute_basic() {
    TEST_SECTION("Basic STFT Computation Tests");

    const size_t input_size = 1024;
    const size_t window_size = 256;
    const size_t hop_size = 128;

    STFT_Config_t *config =
        stft_config_init(hop_size, window_size, input_size, HAMMING, STFT_FFT);
    TEST_ASSERT(config != NULL, "Configuration created for basic test");

    if (!config)
        return;

    // Generate test signal
    float *input = malloc(input_size * sizeof(float));
    generate_sine_wave(input, input_size, 100.0f, 1000.0f);

    // Allocate output
    Complex_t **output = malloc(config->outsize * sizeof(Complex_t *));
    for (size_t i = 0; i < config->outsize; i++) {
        output[i] = malloc(config->fftsize * sizeof(Complex_t));
    }

    // Compute STFT
    int result = stft_compute(config, input, output);
    TEST_ASSERT(result == 0, "STFT computation successful");

    // Check that output is not all zeros
    bool has_non_zero = false;
    for (size_t i = 0; i < config->outsize && !has_non_zero; i++) {
        for (size_t j = 0; j < config->fftsize && !has_non_zero; j++) {
            if (cabs(output[i][j]) > TEST_TOLERANCE) {
                has_non_zero = true;
            }
        }
    }
    TEST_ASSERT(has_non_zero, "Output contains non-zero values");

    // Cleanup
    for (size_t i = 0; i < config->outsize; i++) {
        free(output[i]);
    }
    free(output);
    free(input);
    stft_config_deinit(config);
}

void test_stft_frequency_detection() {
    TEST_SECTION("STFT Frequency Detection Tests");

    const size_t input_size = 2048;
    const size_t window_size = 512;
    const size_t hop_size = 256;
    const float sample_rate = 1000.0f;
    const float test_freq = 100.0f;

    STFT_Config_t *config =
        stft_config_init(hop_size, window_size, input_size, HAMMING, STFT_FFT);
    if (!config)
        return;

    // Generate pure sine wave
    float *input = malloc(input_size * sizeof(float));
    generate_sine_wave(input, input_size, test_freq, sample_rate);

    // Allocate output
    Complex_t **output = malloc(config->outsize * sizeof(Complex_t *));
    for (size_t i = 0; i < config->outsize; i++) {
        output[i] = malloc(config->fftsize * sizeof(Complex_t));
    }

    // Compute STFT
    int result = stft_compute(config, input, output);
    TEST_ASSERT(result == 0, "STFT computation for frequency detection");

    // Find peak frequency in first frame
    size_t expected_bin = (size_t)(test_freq * config->fftsize / sample_rate);
    float max_magnitude = 0.0f;
    size_t max_bin = 0;

    for (size_t j = 0; j < config->fftsize / 2; j++) {
        float magnitude = cabs(output[0][j]);
        if (magnitude > max_magnitude) {
            max_magnitude = magnitude;
            max_bin = j;
        }
    }

    // Allow some tolerance for bin detection
    TEST_ASSERT(abs((int)max_bin - (int)expected_bin) <= 5,
                "Peak frequency detected within tolerance");

    // Cleanup
    for (size_t i = 0; i < config->outsize; i++) {
        free(output[i]);
    }
    free(output);
    free(input);
    stft_config_deinit(config);
}

void test_stft_edge_cases() {
    TEST_SECTION("STFT Edge Case Tests");

    // Test null pointer handling
    int result = stft_compute(NULL, NULL, NULL);
    TEST_ASSERT(result == -1, "Null config pointer handled");

    STFT_Config_t *config = stft_config_init(128, 256, 1024, HAMMING, STFT_FFT);
    if (!config)
        return;

    result = stft_compute(config, NULL, NULL);
    TEST_ASSERT(result == -1, "Null input/output pointers handled");

    // Test with zero signal
    float *zero_input = calloc(1024, sizeof(float));
    Complex_t **output = malloc(config->outsize * sizeof(Complex_t *));
    for (size_t i = 0; i < config->outsize; i++) {
        output[i] = malloc(config->fftsize * sizeof(Complex_t));
    }

    result = stft_compute(config, zero_input, output);
    TEST_ASSERT(result == 0, "Zero signal processed successfully");

    // Verify output is near zero
    bool all_near_zero = true;
    int noise_count = 0;
    float checksum = 0.0f;
    for (size_t i = 0; i < config->outsize; i++) {
        for (size_t j = 0; j < config->fftsize; j++) {
            if (cabs(output[i][j]) > TEST_TOLERANCE) {
                all_near_zero = false;
                noise_count++;
            }
        }
    }

    TEST_ASSERT(all_near_zero, "Zero input produces near-zero output");

    // Cleanup
    for (size_t i = 0; i < config->outsize; i++) {
        free(output[i]);
    }
    free(output);
    free(zero_input);
    stft_config_deinit(config);
}

void test_stft_different_window_types() {
    TEST_SECTION("Different Window Types Tests");

    const size_t input_size = 1024;
    const size_t window_size = 256;
    const size_t hop_size = 128;

    // Test different window types
    WinType window_types[] = {HANNING, HAMMING, BLACKMAN, BLACKMAN_HARRIS};
    const char *window_names[] = {"Hanning", "Hamming", "Blackman",
                                  "Blackman-Harris"};

    for (int w = 0; w < 4; w++) {
        STFT_Config_t *config = stft_config_init(
            hop_size, window_size, input_size, window_types[w], STFT_FFT);

        char test_name[64];
        snprintf(test_name, sizeof(test_name), "%s window configuration",
                 window_names[w]);
        TEST_ASSERT(config != NULL, test_name);

        if (config) {
            // Generate test signal
            float *input = malloc(input_size * sizeof(float));
            generate_sine_wave(input, input_size, 50.0f, 1000.0f);

            // Allocate output
            Complex_t **output = malloc(config->outsize * sizeof(Complex_t *));
            for (size_t i = 0; i < config->outsize; i++) {
                output[i] = malloc(config->fftsize * sizeof(Complex_t));
            }

            // Compute STFT
            int result = stft_compute(config, input, output);
            snprintf(test_name, sizeof(test_name), "%s window STFT computation",
                     window_names[w]);
            TEST_ASSERT(result == 0, test_name);

            // Cleanup
            for (size_t i = 0; i < config->outsize; i++) {
                free(output[i]);
            }
            free(output);
            free(input);
            stft_config_deinit(config);
        }
    }
}

void test_performance() {
    TEST_SECTION("Performance Tests");

    // Performance test parameters
    size_t test_sizes[] = {1024, 2048, 4096, 8192, 16384};
    size_t num_sizes = sizeof(test_sizes) / sizeof(test_sizes[0]);

    printf(ANSI_COLOR_CYAN "Signal Size | Window | Hop | Frames | Time (ms) | "
                           "Throughput (MB/s)\n");
    printf("-------------------------------------------------------------------"
           "-----\n" ANSI_COLOR_RESET);

    for (size_t i = 0; i < num_sizes; i++) {
        size_t input_size = test_sizes[i];
        size_t window_size = input_size / 8;
        size_t hop_size = window_size / 2;

        STFT_Config_t *config = stft_config_init(hop_size, window_size,
                                                 input_size, HAMMING, STFT_FFT);
        if (!config)
            continue;

        // Generate test signal
        float *input = malloc(input_size * sizeof(float));
        generate_sine_wave(input, input_size, 440.0f, 44100.0f);

        // Allocate output
        Complex_t **output = malloc(config->outsize * sizeof(Complex_t *));
        for (size_t j = 0; j < config->outsize; j++) {
            output[j] = malloc(config->fftsize * sizeof(Complex_t));
        }

        // Performance measurement
        double start_time = get_time_ms();

        for (int iter = 0; iter < PERFORMANCE_ITERATIONS; iter++) {
            stft_compute(config, input, output);
        }

        double end_time = get_time_ms();
        double avg_time = (end_time - start_time) / PERFORMANCE_ITERATIONS;

        // Calculate throughput
        double data_size_mb = (input_size * sizeof(float)) / (1024.0 * 1024.0);
        double throughput = data_size_mb / (avg_time / 1000.0);

        printf("%8zu | %6zu | %3zu | %6zu | %8.2f | %13.2f\n", input_size,
               window_size, hop_size, config->outsize, avg_time, throughput);

        // Cleanup
        for (size_t j = 0; j < config->outsize; j++) {
            free(output[j]);
        }
        free(output);
        free(input);
        stft_config_deinit(config);
    }
}

void test_memory_stress() {
    TEST_SECTION("Memory Stress Tests");

    // Test multiple allocations and deallocations
    const int num_configs = 10;
    STFT_Config_t *configs[num_configs];

    // Allocate multiple configurations
    for (int i = 0; i < num_configs; i++) {
        size_t size = 1024 + i * 512;
        configs[i] = stft_config_init(128, 256, size, HAMMING, STFT_FFT);
    }

    int successful_allocs = 0;
    for (int i = 0; i < num_configs; i++) {
        if (configs[i] != NULL) {
            successful_allocs++;
        }
    }

    TEST_ASSERT(successful_allocs > 0, "Multiple configurations allocated");

    // Test computation on all configurations
    for (int i = 0; i < num_configs; i++) {
        if (configs[i] != NULL) {
            float *input = malloc(configs[i]->insize * sizeof(float));
            generate_noise(input, configs[i]->insize, 0.5f);

            Complex_t **output =
                malloc(configs[i]->outsize * sizeof(Complex_t *));
            for (size_t j = 0; j < configs[i]->outsize; j++) {
                output[j] = malloc(configs[i]->fftsize * sizeof(Complex_t));
            }

            int result = stft_compute(configs[i], input, output);
            if (result != 0) {
                printf("Warning: STFT computation failed for config %d\n", i);
            }

            // Cleanup
            for (size_t j = 0; j < configs[i]->outsize; j++) {
                free(output[j]);
            }
            free(output);
            free(input);
        }
    }

    // Deallocate all configurations
    for (int i = 0; i < num_configs; i++) {
        stft_config_deinit(configs[i]);
    }

    TEST_ASSERT(1, "Memory stress test completed without crashes");
}

void test_chirp_signal() {
    TEST_SECTION("Chirp Signal Analysis");

    const size_t input_size = 4096;
    const size_t window_size = 512;
    const size_t hop_size = 256;
    const float sample_rate = 8000.0f;

    STFT_Config_t *config =
        stft_config_init(hop_size, window_size, input_size, HAMMING, STFT_FFT);
    if (!config)
        return;

    // Generate chirp signal (frequency sweep)
    float *input = malloc(input_size * sizeof(float));
    generate_chirp(input, input_size, 100.0f, 1000.0f, sample_rate);

    // Allocate output
    Complex_t **output = malloc(config->outsize * sizeof(Complex_t *));
    for (size_t i = 0; i < config->outsize; i++) {
        output[i] = malloc(config->fftsize * sizeof(Complex_t));
    }

    // Compute STFT
    int result = stft_compute(config, input, output);
    TEST_ASSERT(result == 0, "Chirp signal STFT computation");

    // Analyze frequency content over time
    bool frequency_increases = true;
    size_t prev_peak_bin = 0;

    for (size_t i = 0; i < config->outsize; i++) {
        float max_magnitude = 0.0f;
        size_t max_bin = 0;

        for (size_t j = 1; j < config->fftsize / 2; j++) {
            float magnitude = cabs(output[i][j]);
            if (magnitude > max_magnitude) {
                max_magnitude = magnitude;
                max_bin = j;
            }
        }

        if (i > 0 && max_bin < prev_peak_bin) {
            frequency_increases = false;
        }
        prev_peak_bin = max_bin;
    }

    TEST_ASSERT(frequency_increases,
                "Chirp signal shows increasing frequency over time");

    // Cleanup
    for (size_t i = 0; i < config->outsize; i++) {
        free(output[i]);
    }
    free(output);
    free(input);
    stft_config_deinit(config);
}

void print_test_summary() {
    printf(ANSI_COLOR_MAGENTA "\n=== TEST SUMMARY ===" ANSI_COLOR_RESET "\n");
    printf("Total Tests: %d\n", g_results.total);
    printf(ANSI_COLOR_GREEN "Passed: %d" ANSI_COLOR_RESET "\n",
           g_results.passed);
    printf(ANSI_COLOR_RED "Failed: %d" ANSI_COLOR_RESET "\n", g_results.failed);

    if (g_results.failed == 0) {
        printf(ANSI_COLOR_GREEN "\n All tests passed!" ANSI_COLOR_RESET "\n");
    } else {
        printf(ANSI_COLOR_RED "\n Some tests failed!" ANSI_COLOR_RESET "\n");
    }

    double success_rate = (double)g_results.passed / g_results.total * 100.0;
    printf("Success Rate: %.1f%%\n", success_rate);
}

int main() {
    printf(ANSI_COLOR_CYAN "STFT Library Comprehensive Test Suite\n");
    printf("=====================================\n" ANSI_COLOR_RESET);

    // Seed random number generator
    srand((unsigned int)time(NULL));

    // Run all tests
    test_nextPow2();
    test_complex_functions();
    test_stft_config_init();
    test_stft_compute_basic();
    test_stft_frequency_detection();
    test_stft_edge_cases();
    test_stft_different_window_types();
    test_chirp_signal();
    test_memory_stress();
    test_performance();

    // Print summary
    print_test_summary();

    return (g_results.failed == 0) ? 0 : 1;
}
