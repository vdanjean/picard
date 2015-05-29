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

}
