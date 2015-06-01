package picard.analysis;

import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;

/**
 * Created by davidben on 5/18/15.
 */
public class TheoreticalSensitivityTest {

    @Test
    public void testRouletteWheel() throws Exception {

        //test that a deterministic roulette wheel only gives one value
        List<Double> deterministicWeights = Arrays.asList(0.0, 1.0, 0.0);
        final TheoreticalSensitivity.RouletteWheel deterministicWheel = new TheoreticalSensitivity.RouletteWheel(deterministicWeights);
        for (int n = 0; n < 10; n++) Assert.assertEquals(deterministicWheel.draw(), 1);

        //test the sums of this deterministic wheel: a sum of n 1's equals n
        List<ArrayList<Integer>> deterministicSums = deterministicWheel.sampleCumulativeSums(10, 1);
        for (int n = 0; n < 10; n++) Assert.assertEquals(deterministicSums.get(n).get(0), (Integer) n);

        //test that an unfair coin w/ probability 1/4 of heads gives roughly 1/4 heads
        double p = 0.25;
        int total_heads = 0;
        int N = 10000;
        double stdDev = Math.sqrt(N*p*(1-p));   //the stddev of the sample total heads
        List<Double> unfairCoinWeights = Arrays.asList(1-p, p);
        final TheoreticalSensitivity.RouletteWheel coinFlipWheel = new TheoreticalSensitivity.RouletteWheel(unfairCoinWeights);
        for (int n = 0; n < N; n++) total_heads += coinFlipWheel.draw();
        Assert.assertEquals(total_heads, N*p, 10*stdDev);
    }

    @Test
    public void testProportionsAboveThresholds() throws Exception {
        List<ArrayList<Integer>> sums = new ArrayList<ArrayList<Integer>>();
        sums.add(new ArrayList<Integer>(Arrays.asList(0,0,0)));
        sums.add(new ArrayList<Integer>(Arrays.asList(10, 10)));
        sums.add(new ArrayList<Integer>(Arrays.asList(5, 11, -2, 4)));
        List<Double> thresholds = Arrays.asList(-1.0, 1.0, 6.0);
        Assert.assertEquals(sums.size(), 3);
        Assert.assertEquals(thresholds.size(), 3);

        List<ArrayList<Double>> proportions = TheoreticalSensitivity.proportionsAboveThresholds(sums, thresholds);
        Assert.assertEquals(proportions.size(), 3);

        Assert.assertEquals(proportions.get(0).get(0), (double) 3/3);
        Assert.assertEquals(proportions.get(0).get(1), (double) 0/3);
        Assert.assertEquals(proportions.get(0).get(2), (double) 0/3);
        Assert.assertEquals(proportions.get(1).get(0), (double) 2/2);
        Assert.assertEquals(proportions.get(1).get(1), (double) 2/2);
        Assert.assertEquals(proportions.get(1).get(2), (double) 2/2);
        Assert.assertEquals(proportions.get(2).get(0), (double) 3/4);
        Assert.assertEquals(proportions.get(2).get(1), (double) 3/4);
        Assert.assertEquals(proportions.get(2).get(2), (double) 1/4);
    }

    @Test
    public void testHetAltDepthDistribution() throws Exception {
        int N = 6;
        double p = 0.5;
        List<ArrayList<Double>> distribution = TheoreticalSensitivity.hetAltDepthDistribution(N);

        for (int n = 0; n < N-1; n++) {
            for (int m = 0; m <= n; m++) {
                //java has no built-in binomial coefficient
                //when this is in hellbender, use apache commons
                int binomialCoefficient = 1;
                for (int i = n; i > (n - m); i--) binomialCoefficient *= i;
                for (int i = m; i > 0; i--) binomialCoefficient /= i;

                Assert.assertEquals(distribution.get(n).get(m), binomialCoefficient*Math.pow(p,n));
            }
        }
    }

    //test that large-sample sums from the RouletteWheel converge to a normal distribution
    //using the empirical CDF as measured by proportionsAboveThresholds
    @Test
    public void testCentralLimitTheorem() throws Exception {
        //use a RouletteWheel that gives 0, 1, 2 with equal probability
        List<Double> weights = Arrays.asList(1.0, 1.0, 1.0);
        final TheoreticalSensitivity.RouletteWheel wheel = new TheoreticalSensitivity.RouletteWheel(weights);

        int sampleSize = 1000;
        int numSummands = 100;

        //the mean and standard deviation of a single roulette draw and of many draws
        double muSingleDraw = 1.0;
        double sigmaSingleDraw = Math.sqrt(2.0 / 3.0);
        double mu = numSummands * muSingleDraw;
        double sigma = Math.sqrt(numSummands) * sigmaSingleDraw;

        //test the sums of this deterministic wheel: a sum of n 1's equals n
        List<ArrayList<Integer>> sums = wheel.sampleCumulativeSums(numSummands, sampleSize);
        //we only want the last set of sums, those with numSummands summands
        sums.subList(0, sums.size() - 1).clear();

        Assert.assertEquals(sums.size(), 1);

        //test whether the number of elements within one standard deviation agrees with the normal distribution
        List<Double> thresholds = Arrays.asList(mu - sigma, mu + sigma);

        //sums is 1 x sampleSize, thresholds is a 2-vector, so proportions is 1 x 2
        List<ArrayList<Double>> proportions = TheoreticalSensitivity.proportionsAboveThresholds(sums, thresholds);
        double empiricalProportionWithinOneSigma = proportions.get(0).get(0) - proportions.get(0).get(1);

        //the proportion within one sigma for the normal distribution
        //hence whether any element falls within one sigma is a Bernoulli variable
        double theoreticalProportionWithinOneSigma = 0.682689492;
        double samplingStandardDeviationOfProportion = Math.sqrt(theoreticalProportionWithinOneSigma*(1-theoreticalProportionWithinOneSigma) /  sampleSize);

        Assert.assertEquals(empiricalProportionWithinOneSigma, theoreticalProportionWithinOneSigma, 5*samplingStandardDeviationOfProportion);

    }
}
