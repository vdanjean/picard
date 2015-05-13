package picard.analysis;

import java.util.*;

/**
 * Created by David Benjamin on 5/13/15.
 */
public class TheoreticalSensitivity {
    /*
    @param the probability of depth n is depthDistribution[n] for n = 0, 1. . . N - 1
    @param the probability of quality q is qualityDistribution[q] for q = 0, 1. . . Q
    @param sample size is the number of random sums of quality scores for each m
    @param logOddsThreshold is the log_10 of the likelihood ratio required to call a SNP,
    for example 5 if the variant likelihood must be 10^5 times greater
     */
    public double hetSNPSensitivity(final List<Double> depthDistribution, final List<Double> qualityDistribution,
                                  final int sampleSize, final double logOddsThreshold) {
        int N = depthDistribution.size();
        RouletteWheel qualitySampler = new RouletteWheel(qualityDistribution);

        //qualitySums[m] is a random sample of sums of m quality scores, for m = 0, 1, N - 1
        List<ArrayList<Integer>> qualitySums = qualitySampler.sampleCumulativeSums(N, sampleSize);

        //if a quality sum of m qualities exceeds the quality sum threshold for n total reads, a SNP is called
        ArrayList<Double> qualitySumThresholds = new ArrayList<Double>(N);
        for (int n = 0; n < N; n++) qualitySumThresholds.add(10*(n*Math.log10(2) + logOddsThreshold));


        //probabilityToExceedThreshold[m][n] is the probability that the sum of m quality score
        //exceeds the nth quality sum threshold
        List<ArrayList<Double>> probabilityToExceedThreshold = proportionsAboveThresholds(qualitySums, qualitySumThresholds);
        List<ArrayList<Double>> altDepthDistribution = hetAltDepthDistribution(N);
        double result = 0.0;
        for (int n = 0; n < N; n++) {
            for (int m = 0; m <= n; m++) {
                result += depthDistribution.get(n) * altDepthDistribution.get(n).get(m)* probabilityToExceedThreshold.get(m).get(n);
            }
        }
        return result;
    }

    //given L lists of lists and N thresholds, count the proportion of each list above each threshold
    public static List<ArrayList<Double>> proportionsAboveThresholds(final List<ArrayList<Integer>> sums, List<Double> thresholds) {
        int N = thresholds.size();
        int L = sums.size();
        ArrayList<ArrayList<Double>> result = new ArrayList<ArrayList<Double>>(L);

        for (int m = 0; m < L; m++) {
            ArrayList<Double> mthRow = new ArrayList<Double>(Collections.nCopies(N,0.0));
            Collections.sort(sums.get(m));
            int sampleSize = sums.get(m).size();
            int n = 0;
            int j = 0;  //index within the ordered set of qualitySum samples
            while (n < N && j < sampleSize) {
                if (thresholds.get(n) > sums.get(m).get(j)) {
                    j++;
                }
                else {
                    n++;
                    result.get(m).set(n, (double) (sampleSize - j)/sampleSize);
                }
            }
            result.add(mthRow);
        }
        return result;

    }


    //Utility function for making table of binomial distribution probabilities nCm * (0.5)^n
    //for n = 0, 1 . . . N - 1 and m = 0, 1. . . n
    public static List<ArrayList<Double>> hetAltDepthDistribution(final int N) {
        List<ArrayList<Double>> table = new ArrayList<ArrayList<Double>>();
        for (int n = 0; n < N; n++) {
            ArrayList<Double> nthRow = new ArrayList<Double>();

            ArrayList<Double> prevRow = n == 0 ? null : table.get(n - 1);
            //Add the m = 0 element
            nthRow.add(n == 0 ? 1 : 0.5*prevRow.get(0));

            //add elements 1 through n - 1.  Note that nCm = (n-1)C(m-1) * (n/m)
            for (int m = 1; m < n; m++) nthRow.add(0.5*n*prevRow.get(m - 1)/m);

            //add element m = n, which equals the m = 0 element
            nthRow.add(nthRow.get(0));

            table.add(nthRow);
        }

        return table;

    }


    /*
    Perform random draws from {0, 1. . . N - 1} according to a list of relative probabilities.

    We use an O(1) stochastic acceptance algorithm -- see Physica A, Volume 391, Page 2193 (2012) --
    which works well when the ratio of maximum weight to average weight is not large.
     */
    public static class RouletteWheel {
        private List<Double> probabilities;
        private int N;

        RouletteWheel(final List<Double> weights) {
            N = weights.size();

            probabilities = new ArrayList<Double>();
            double wMax = Collections.max(weights);
            for (final double w : weights) {
                probabilities.add(w/wMax);
            }
        }

        public int draw() {
            while (true) {
                int n = (int) (N * Math.random());
                if (Math.random() < probabilities.get(n)) return n;
            }
        }

        //get samples of sums of 0, 1, 2,. . .  N - 1 draws
        public List<ArrayList<Integer>> sampleCumulativeSums(final int maxNumberOfSummands, final int sampleSize) {
            List<ArrayList<Integer>> result = new ArrayList<ArrayList<Integer>>();
            for (int m = 0; m < maxNumberOfSummands; m++) result.add(new ArrayList<Integer>());

            for (int iteration = 0; iteration < sampleSize; iteration++) {
                int cumulativeSum = 0;
                for (int m = 0; m < maxNumberOfSummands; m++) {
                    result.get(m).add(cumulativeSum);
                    cumulativeSum += draw();
                }
            }
            return result;
        }

    }

}
