#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <signal.h>
#include <string.h>

#define BLOCK_SIZE 10000000ULL // 10 million numbers per segment
#define PRIME_BATCH 1000000ULL   // batch writes to file
#define U64FMT "%" PRIu64

static volatile sig_atomic_t stop = 0;
static uint64_t biggest = 2;
static uint64_t start_from = 2;

const uint32_t wheel_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23};
const int wheel_pcount = sizeof(wheel_primes) / sizeof(wheel_primes[0]);

static uint64_t wheel_lcm = 0;
static uint64_t *wheel_residues = NULL;
static uint64_t *wheel_deltas = NULL;
static int wheel_size = 0;

void handle_signal(int sig) { stop = 1; }

uint64_t gcd64(uint64_t a, uint64_t b) {
    while (b) { uint64_t t = b; b = a % b; a = t; }
    return a;
}

uint64_t lcm64(uint64_t a, uint64_t b) {
    return a / gcd64(a, b) * b;
}

void compute_wheel(void) {
    wheel_lcm = 1;
    for (int i = 0; i < wheel_pcount; i++)
        wheel_lcm = lcm64(wheel_lcm, wheel_primes[i]);

    wheel_residues = malloc(wheel_lcm * sizeof(uint64_t));
    if (!wheel_residues) { perror("malloc wheel_residues"); exit(1); }

    int count = 0;
    for (uint64_t n = 1; n <= wheel_lcm; n++) {
        int coprime = 1;
        for (int i = 0; i < wheel_pcount; i++)
            if (n % wheel_primes[i] == 0) { coprime = 0; break; }
        if (coprime) wheel_residues[count++] = n;
    }
    wheel_size = count;

    wheel_deltas = malloc(wheel_size * sizeof(uint64_t));
    if (!wheel_deltas) { perror("malloc wheel_deltas"); exit(1); }

    for (int i = 0; i < wheel_size - 1; i++)
        wheel_deltas[i] = wheel_residues[i + 1] - wheel_residues[i];
    wheel_deltas[wheel_size - 1] = wheel_lcm - wheel_residues[wheel_size - 1] + wheel_residues[0];
}

#define GET_BIT64(arr,i) (((arr[(i)>>6] >> ((i)&63)) & 1ULL))
#define SET_BIT64(arr,i) (arr[(i)>>6] |= (1ULL << ((i)&63)))
#define CLEAR_BLOCK(arr,words) memset(arr,0,(words)*sizeof(uint64_t))

void simple_sieve(uint64_t limit, uint64_t **primes, size_t *count) {
    char *mark = calloc(limit + 1, 1);
    if (!mark) { perror("simple_sieve calloc"); exit(1); }

    for (uint64_t i = 2; i * i <= limit; i++)
        if (!mark[i])
            for (uint64_t j = i * i; j <= limit; j += i)
                mark[j] = 1;

    size_t c = 0;
    for (uint64_t i = 2; i <= limit; i++)
        if (!mark[i]) c++;

    *primes = malloc(c * sizeof(uint64_t));
    if (!*primes) { perror("simple_sieve malloc"); exit(1); }

    *count = 0;
    for (uint64_t i = 2; i <= limit; i++)
        if (!mark[i])
            (*primes)[(*count)++] = i;

    free(mark);
}

void save_prime_buffer(uint64_t *buffer, size_t count, FILE *f) {
    if (fwrite(buffer, sizeof(uint64_t), count, f) != count)
        perror("fwrite prime buffer");
}

void load_last_prime_from_file(const char *filename) {
    FILE *f = fopen(filename, "rb");
    if (!f) { biggest = 2; start_from = 2; return; }
    if (fseek(f, -sizeof(uint64_t), SEEK_END) == 0) {
        uint64_t last_prime;
        if (fread(&last_prime, sizeof(uint64_t), 1, f) == 1) {
            biggest = last_prime;
            start_from = last_prime + 1;
        }
    }
    fclose(f);
}

int main(void) {
    signal(SIGINT, handle_signal);
    signal(SIGTERM, handle_signal);

    const char *filename = "primes64.bin";
    load_last_prime_from_file(filename);
    compute_wheel();

    uint64_t *base_primes = NULL;
    size_t base_count = 0;
    simple_sieve(1000000ULL, &base_primes, &base_count);

    size_t block_words = (BLOCK_SIZE + 63) / 64;
    uint64_t *is_prime = calloc(block_words, sizeof(uint64_t));
    if (!is_prime) { perror("calloc"); return 1; }

    FILE *outfile = fopen(filename, "ab");
    if (!outfile) { perror("fopen"); return 1; }

    uint64_t prime_buffer[PRIME_BATCH];
    size_t prime_buffer_count = 0;

    uint64_t base = (start_from / wheel_lcm) * wheel_lcm;

    printf("Resuming from " U64FMT " using wheel of " U64FMT " (%d residues)\n",
           start_from, wheel_lcm, wheel_size);

    while (!stop) {
        uint64_t segment_start = base;
        uint64_t segment_end = base + BLOCK_SIZE;
        CLEAR_BLOCK(is_prime, block_words);

        for (size_t i = 0; i < base_count; i++) {
            uint64_t p = base_primes[i];
            uint64_t p2 = p * p;
            if (p2 >= segment_end) break;

            uint64_t first = (segment_start + p - 1) / p * p;
            if (first < p2) first = p2;

            for (uint64_t j = first - segment_start; j < BLOCK_SIZE; j += p)
                SET_BIT64(is_prime, j);
        }

        uint64_t n = segment_start + wheel_residues[0];
        int w = 0;

        while (n < segment_end && !stop) {
            if (!GET_BIT64(is_prime, n - segment_start)) {
                prime_buffer[prime_buffer_count++] = n;
                if (prime_buffer_count == PRIME_BATCH) {
                    save_prime_buffer(prime_buffer, PRIME_BATCH, outfile);
                    prime_buffer_count = 0;
                }
                if (n > biggest) biggest = n;
            }

            n += wheel_deltas[w];
            w = (w + 1) % wheel_size;
        }

        if (prime_buffer_count > 0) {
            save_prime_buffer(prime_buffer, prime_buffer_count, outfile);
            prime_buffer_count = 0;
        }

        // Accurate total primes based on file size
        long long file_size = ftell(outfile);
        uint64_t total_primes_found = file_size / sizeof(uint64_t);

        printf("Total prime numbers found: %" PRIu64 " | Current candidate: " U64FMT " | Largest prime: " U64FMT "\r",
               total_primes_found, n, biggest);
        fflush(stdout);

        base += BLOCK_SIZE;
    }

    fclose(outfile);

    printf("\nExiting. Largest prime found: " U64FMT "\n", biggest);

    free(base_primes);
    free(is_prime);
    free(wheel_residues);
    free(wheel_deltas);
    return 0;
}
