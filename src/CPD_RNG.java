import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

/**
 * CPD_RNG — Compound Poisson Distribution Random Number Generator
 *          with R2.3 same-seed comparison experiments.
 *
 * Per cycle:
 *   1. Draw n ~ Poisson(mu) from the LCG state.
 *   2. Run n Fisher-Yates permutations on a length-4 array.
 *   3. Measure elapsed time T (in 100-ns ticks).
 *   4. Emit N_p = n mod R and advance state via seed = (seed << nbits) + (T mod R).
 */
public class CPD_RNG {

    // ------------------------------------------------------------------
    // State
    // ------------------------------------------------------------------
    private static Random random = new Random(347563409L);
    private static long   seed   = 347563409L;

    private static final long INIT_SEED = 347563409L;

    // ------------------------------------------------------------------
    // 48-bit LCG
    // ------------------------------------------------------------------
    protected static synchronized int next(int bits) {
        seed = (seed * 0x5DEECE66DL + 0xBL) & ((1L << 48) - 1);
        return (int) (seed >>> (48 - bits));
    }

    public static int nextInt(int n) {
        if (n <= 0) throw new IllegalArgumentException("n must be positive");
        int bits, val;
        do {
            bits = next(31);
            val  = bits % n;
        } while (bits - val + (n - 1) < 0);
        return val;
    }

    public static double nextDouble() {
        long bits = ((long) next(26) << 27) + next(27);
        return bits / (double) (1L << 53);
    }

    // ------------------------------------------------------------------
    // Poisson sampler
    // ------------------------------------------------------------------
    private static int samplePoisson(double lambda) {
        if (lambda <= 200) {
            double L = Math.exp(-lambda);
            int k = 0;
            double p = 1.0;
            do {
                k++;
                p *= nextDouble();
            } while (p > L);
            return k - 1;
        } else {
            return (int) Math.round(random.nextGaussian() * Math.sqrt(lambda) + lambda);
        }
    }

    // ------------------------------------------------------------------
    // Fisher-Yates
    // ------------------------------------------------------------------
    private static void permutArray(int[] array) {
        for (int i = array.length - 1; i > 0; i--) {
            int j = nextInt(i + 1);
            int temp = array[j];
            array[j] = array[i];
            array[i] = temp;
        }
    }

    // ------------------------------------------------------------------
    // One Poisson cycle: returns { T (ticks), n_p }
    // ------------------------------------------------------------------
    public static long[] QPP_RNG_Poi(int[] array2) {
        int n = samplePoisson(100);

        long start = System.nanoTime();
        long[] ret = new long[2];
        int[] array = Arrays.copyOf(array2, array2.length);

        start = System.nanoTime();
        for (int k = 0; k < n; k++) {
            permutArray(array);
        }
        long end = System.nanoTime();

        ret[0] = (end - start) / 100;   // ns -> 100-ns ticks
        ret[1] = n;
        return ret;
    }

    // ------------------------------------------------------------------
    // Byte builder: pack 8/n symbols per byte, advance state each sub-cycle
    //
    // n = bits per symbol; R = 1 << n; 8/n sub-cycles per byte.
    // ------------------------------------------------------------------
    public static long[] QPP_RNG_Byte(int n, int[] array1) {
        long[] ret = new long[2];
        int R = 1 << n;
        ret[0] = 0;
        ret[1] = 0;
        for (int i = 0; i < 8 / n; i++) {
            long[] u = QPP_RNG_Poi(array1);
            seed = ((seed << n) + (u[0] % R)) & ((1L << 48) - 1);
            ret[0] += ((u[0] % R) << (i * n));
            ret[1] += ((u[1] % R) << (i * n));
        }
        ret[0] = ret[0] % 256;
        ret[1] = ret[1] % 256;
        return ret;
    }

    // ------------------------------------------------------------------
    // Main uniformity experiment
    // ------------------------------------------------------------------
    public static void iidPermutation() {
        int iterations = 1_000_000;

        long[] deltas  = new long[iterations];
        long[] delta1  = new long[iterations];
        int[]  delta2  = new int[iterations];
        int[]  deltas2 = new int[iterations];
        byte[] delta3  = new byte[iterations];   // N_p byte stream
        byte[] delta4  = new byte[iterations];   // T   byte stream

        int[] array01 = { 3, 1, 2, 4 };

        // Warm-up: discard first 80 cycles
        seed = INIT_SEED;
        for (int i = 0; i < 80; i++) {
            QPP_RNG_Byte(8, array01);
        }

        long startTime = System.nanoTime();

        seed = INIT_SEED;

        for (int i = 0; i < iterations; i++) {
            long[] q = QPP_RNG_Byte(4, array01);

            delta1[i]  = q[0];                 // T (ticks)      — physical
            deltas[i]  = q[1] % 256;           // N_p mod 256    — combinatorial

            delta2[i]  = (int) delta1[i];
            deltas2[i] = (int) deltas[i];
            delta3[i]  = (byte) deltas[i];     // N_p byte stream
            delta4[i]  = (byte) delta1[i];     // T   byte stream
        }

        long endTime = System.nanoTime();

        System.out.println("\nShannon Entropy N_p=" + calculateEntropy(deltas)
                + " : T = " + calculateEntropy(delta1)
                + " bits, time delta= " + (endTime - startTime) / 1_000_000.0 + " ms");
        System.out.println("NIST 90B, min-entropy N_p= " + estimateMinEntropy(delta3)
                + "; T = " + estimateMinEntropy(delta4));
        System.out.println("Chi-Square Test Score: N_p=" + chiSquareLarge(deltas2)
                + " : T=" + chiSquareLarge(delta2));

        // Write N_p byte stream to file (for NIST SP 800-90B)
        try {
            writeByteArrayAsBytes(delta3, "np_output.bytes");
            System.out.println("N_p byte stream written to np_output.bytes");
        } catch (IOException e) {
            System.err.println("Error writing np_output.bytes: " + e.getMessage());
        }

        // Write T byte stream to file (for NIST SP 800-90B)
        try {
            writeByteArrayAsBytes(delta4, "t_output.bytes");
            System.out.println("T byte stream written to t_output.bytes");
        } catch (IOException e) {
            System.err.println("Error writing t_output.bytes: " + e.getMessage());
        }
    }

    // ==================================================================
    // R2.3 experiments: same initial condition, multiple runs
    // ==================================================================

    /**
     * One cycle for the R2.3 experiment.
     * Returns { T, n mod R }.
     * If useFeedback, advances the class-level seed by (T mod R).
     */
    private static long[] cycleForR23(int[] array, int R, boolean useFeedback) {
        int n = samplePoisson(100.0);

        long start = System.nanoTime();
        for (int k = 0; k < n; k++) permutArray(array);
        long end = System.nanoTime();

        long T = (end - start) / 100;      // ns -> 100-ns ticks
        long y = n % R;

        if (useFeedback) {
            int nbits = Integer.numberOfTrailingZeros(R);   // log2(R)
            seed = ((seed << nbits) + (T % R)) & ((1L << 48) - 1);
        }

        return new long[]{ T, y };
    }

    /**
     * Same-seed comparison with optional warm-up.
     *
     * @param seedArray   initial array
     * @param warmupCycles number of discarded cycles before recording
     * @param numCycles   number of recorded cycles
     * @param R           output modulus
     * @param useFeedback true -> state advance by T mod R
     * @return long[numCycles][2] with [i][0]=T, [i][1]=n mod R
     */
    public static long[][] R23_TimeSeries(int[] seedArray,
                                          int warmupCycles,
                                          int numCycles,
                                          int R,
                                          boolean useFeedback) {

        int[] array = Arrays.copyOf(seedArray, seedArray.length);

        // Warm-up (discarded)
        for (int w = 0; w < warmupCycles; w++) {
            cycleForR23(array, R, useFeedback);
        }

        // Recording
        long[][] series = new long[numCycles][2];
        for (int i = 0; i < numCycles; i++) {
            long[] c = cycleForR23(array, R, useFeedback);
            series[i][0] = c[0];
            series[i][1] = c[1];
        }
        return series;
    }

    private static void printR23Table(String title,
                                      long[][][] runs,
                                      int numRuns,
                                      int numCycles) {
        System.out.println("\n=== " + title + " ===");
        System.out.print("Run  |");
        for (int i = 1; i <= numCycles; i++) System.out.printf(" c%-2d", i);
        System.out.println(" | div@");

        for (int r = 0; r < numRuns; r++) {
            System.out.printf("%4d |", r + 1);
            for (int i = 0; i < numCycles; i++) {
                System.out.printf(" %2d ", runs[r][i][1]);
            }
            int div = -1;
            for (int i = 0; i < numCycles; i++) {
                if (runs[r][i][1] != runs[0][i][1]) { div = i + 1; break; }
            }
            System.out.println(" | " + (div < 0 ? "none" : div));
        }
    }

    /**
     * Run the four R2.3 experiments:
     *   - no warm-up, feedback off
     *   - no warm-up, feedback on
     *   - warm-up 8, feedback off
     *   - warm-up 8, feedback on
     */
    public static void runR23Experiments() {
        int[] seedArray = { 3, 1, 2, 4 };
        int R = 16;
        int numCycles = 10;
        int numRuns = 10;

        for (int warmup : new int[]{ 0, 8 }) {
            for (boolean feedback : new boolean[]{ false, true }) {
                String title = "Warm-up " + warmup
                        + ", feedback " + (feedback ? "ON" : "OFF");
                long[][][] runs = new long[numRuns][][];
                for (int r = 0; r < numRuns; r++) {
                    seed = INIT_SEED;
                    runs[r] = R23_TimeSeries(seedArray, warmup, numCycles, R, feedback);
                }
                printR23Table(title, runs, numRuns, numCycles);
            }
        }
    }

    // ------------------------------------------------------------------
    // Statistics helpers
    // ------------------------------------------------------------------
    public static double calculateEntropy(long[] deltas) {
        Map<Long, Long> freqMap = new HashMap<>();
        for (long d : deltas) freqMap.put(d, freqMap.getOrDefault(d, 0L) + 1);
        long total = deltas.length;
        double entropy = 0.0;
        for (long f : freqMap.values()) {
            double p = (double) f / total;
            if (p > 0) entropy -= p * (Math.log(p) / Math.log(2));
        }
        return entropy;
    }

    public static double estimateMinEntropy(byte[] data) {
        if (data == null || data.length == 0)
            throw new IllegalArgumentException("Data must be non-empty");
        int[] counts = new int[256];
        for (byte b : data) counts[b & 0xFF]++;
        int maxCount = Arrays.stream(counts).max().orElse(0);
        double pMax = (double) maxCount / data.length;
        double hMin = -Math.log(pMax) / Math.log(2);
        return Math.min(hMin, 8.0);
    }

    public static double chiSquareTest(int[] frequencies, int expectedFrequency) {
        double chi = 0.0;
        double expFreq = 1.0 * expectedFrequency;
        for (int freq : frequencies) {
            double diff = freq - expFreq;
            chi += (diff * diff) / expFreq;
        }
        return chi;
    }

    public static double chiSquareLarge(int[] data) {
        final int BLOCK_SIZE = 1_000_000;
        final int NUM_BINS   = 256;
        int length   = data.length;
        int segments = (length + BLOCK_SIZE - 1) / BLOCK_SIZE;
        double sumChi = 0.0;
        for (int s = 0; s < segments; s++) {
            int start = s * BLOCK_SIZE;
            int end   = Math.min(start + BLOCK_SIZE, length);
            int size  = end - start;
            int[] freq = new int[NUM_BINS];
            for (int i = start; i < end; i++) {
                int v = data[i];
                if (v >= 0 && v < NUM_BINS) freq[v]++;
                else throw new IllegalArgumentException("Value out of range: " + v);
            }
            int expectedFreq = size / NUM_BINS;
            sumChi += chiSquareTest(freq, expectedFreq);
        }
        return sumChi / segments;
    }

    /**
     * Writes a byte array directly to a file.
     * Used for the NIST SP 800-90B byte streams (N_p and T).
     */
    public static void writeByteArrayAsBytes(byte[] data, String filePath)
            throws IOException {
        try (FileOutputStream fos = new FileOutputStream(filePath)) {
            fos.write(data);
        }
    }

    // ------------------------------------------------------------------
    // main — mode selection
    // ------------------------------------------------------------------
    public static void main(String[] args) {
        int test_mode = 2;

        if (test_mode == 2) {
            for (int i = 0; i < 5; i++) {
                random = new Random(INIT_SEED);
                iidPermutation();
            }
        } else if (test_mode == 4) {
            runR23Experiments();
        }
    }
}