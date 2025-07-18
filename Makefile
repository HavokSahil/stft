CC := gcc
CFLAGS := -g -O0 -Wall -lm -I. -I./pffft

%.o: pffft/%.c
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test_stft: test_stft.o pffft.o
	$(CC) $(CFLAGS) -o $@ $^

test_plot: test_plot.o pffft.o
	$(CC) $(CFLAGS) -o $@ $^

.PHONY: test
test: test_stft
	@./test_stft

.PHONY: clean
clean:
	@rm -f *.o test_stft test_plot

.DEFAULT_TARGET: test
