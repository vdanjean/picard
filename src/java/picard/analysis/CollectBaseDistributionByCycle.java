package picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.util.Arrays;
import java.util.List;

import htsjdk.samtools.util.StringUtil;
import picard.PicardException;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;
import picard.util.RExecutor;



@CommandLineProgramProperties(
        usage = CollectBaseDistributionByCycle.USAGE_SUMMARY + CollectBaseDistributionByCycle.USAGE_DETAILS,
        usageShort = CollectBaseDistributionByCycle.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class CollectBaseDistributionByCycle extends SinglePassSamProgram {
        static final String USAGE_SUMMARY = "Program to chart the nucleotide distribution per cycle in a SAM or BAM file";
        static final String USAGE_DETAILS = "This program charts the nucleotide distribution per cycle in a SAM or BAM file " +
                "to enable assessments of systematic errors at specific positions in the reads.<br /><br />" +
                "" +
                "Although the the use of this tool is described for the Illumina sequencing platform, it can be used for" +
                " any high-throughput sequencing technology including SOLiD/Ion Torrent (Life Technologies), 454 (Roche), " +
                "Nanopore (Oxford Nanopore Technolgies), and PACBIO RS II (Pacific Biosciences).<br /><br /> " +

                "The Illumina sequencing platform carries out a sequence by synthesis process, whereby each cycle represents " +
                "the addition of a new nucleotide to a growing oligonucleotide chain.  " +
                "The growing synthetic oligonucleotide chains are represented as the observed reads in a sequencing run.  " +
                "However, as the run progresses, the accuracy of nucleotide incorporation by the DNA polymerase can degrade, such that " +
                "base calls become ambiguous towards the end of a run and are indicated as unassigned nucleotides \"N\"s in" +
                " a read.  Other sources of errors include the inadequate flushing of the flow cell, which lead to " +
                "inaccurate base calls that accumulate at each subsequent cycle.  However, each sequencing instrument platform can" +
                " have unique sources of systematic errors e.g. A-T bias (SOLiD), deletions (Oxford Nanopore), and " +
                "homopolymers (Roche) etc.  Sequencing platform evaluations and comparisons can be found at: " +
                "http://www.molecularecologist.com/next-gen-table-3c-2014/ <br /><br />" +
                "" +
                "Using the CollectBaseDistributionByCycle tool, increased numbers of miscalled bases will be reflected " +
                "in both base distribution changes and increases in the number of \"N\"s.  " +
                "In general, we expect that for any given cycle, or position within reads, the relative proportions of " +
                "A, T, C and G should reflect the AT:GC content of the organism's genome.  Thus, for all four nucleotides, flattish lines " +
                "would be expected.  Deviations from this expectation, for example a spike of A at a particular " +
                "cycle (position within reads), would suggest a systematic sequencing error."+
                "" +
                "Although previous workflows involved the discarding of low-quality tails through quality trimming, " +
                "GATK's current Best Practices includes the Base Quality Score Recalibrator (BQSR) tool.  The BQSR tool is " +
                "quality-aware and will handle systemic biases that covary with the reads, thus reducing systematic errors " +
                "in a platform-independent manner.  For more information on the GATK Best Practices workflow, visit:" +
                " http://www.broadinstitute.org/gatk/guide/best-practices/"+
                "" +
                "<br /><h4>Usage example:</h4>" +
                "<pre>" +
                "java -jar picard.jar CollectBaseDistributionByCycle \\<br />" +
                "      CHART=collectbasedistbycycle.pdf \\<br />" +
                "      I=input.bam \\<br />" +
                "      O=output.txt" +
                "</pre>" +
                "<hr />";
    @Option(shortName = "CHART", doc = "A file (with .pdf extension) to write the chart to.")
    public File CHART_OUTPUT;

    @Option(doc = "If set to true, calculate the base distribution over aligned reads only.")
    public boolean ALIGNED_READS_ONLY = false;

    @Option(doc = "If set to true calculate the base distribution over PF reads only.")
    public boolean PF_READS_ONLY = false;

    private HistogramGenerator hist;
    private String plotSubtitle = "";
    private final Log log = Log.getInstance(CollectBaseDistributionByCycle.class);

    public static void main(String[] args) {
        System.exit(new CollectBaseDistributionByCycle().instanceMain(args));
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(CHART_OUTPUT);
        final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
        if (readGroups.size() == 1) {
            plotSubtitle = StringUtil.asEmptyIfNull(readGroups.get(0).getLibrary());
        }
        hist = new HistogramGenerator();
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        if ((PF_READS_ONLY) && (rec.getReadFailsVendorQualityCheckFlag())) {
            return;
        }
        if ((ALIGNED_READS_ONLY) && (rec.getReadUnmappedFlag())) {
            return;
        }
        if (rec.isSecondaryOrSupplementary()) {
            return;
        }
        hist.addRecord(rec);
    }

    @Override
    protected void finish() {
        final MetricsFile<BaseDistributionByCycleMetrics, ?> metrics = getMetricsFile();
        hist.addToMetricsFile(metrics);
        metrics.write(OUTPUT);
        if (hist.isEmpty()) {
            log.warn("No valid bases found in input file. No plot will be produced.");
        } else {
            final int rResult = RExecutor.executeFromClasspath("picard/analysis/baseDistributionByCycle.R",
                    OUTPUT.getAbsolutePath(),
                    CHART_OUTPUT.getAbsolutePath(),
                    INPUT.getName(),
                    plotSubtitle);
            if (rResult != 0) {
                throw new PicardException("R script nucleotideDistributionByCycle.R failed with return code " + rResult);
            }
        }
    }

    private class HistogramGenerator {
        private int maxLengthSoFar = 0;
        final private long[][] firstReadTotalsByCycle = new long[5][maxLengthSoFar];
        private long[] firstReadCountsByCycle = new long[maxLengthSoFar];
        final private long[][] secondReadTotalsByCycle = new long[5][maxLengthSoFar];
        private long[] secondReadCountsByCycle = new long[maxLengthSoFar];
        private boolean seenSecondEnd = false;

        private int baseToInt(final byte base) {
            switch (base) {
                case 'A':
                case 'a':
                    return 0;
                case 'C':
                case 'c':
                    return 1;
                case 'G':
                case 'g':
                    return 2;
                case 'T':
                case 't':
                    return 3;
            }
            return 4;
        }

        void addRecord(final SAMRecord rec) {
            final byte[] bases = rec.getReadBases();
            if (bases == null) {
                return;
            }
            final int length = bases.length;
            final boolean rc = rec.getReadNegativeStrandFlag();
            ensureArraysBigEnough(length + 1);
            if ((rec.getReadPairedFlag()) && (rec.getSecondOfPairFlag())) {
                seenSecondEnd = true;
                for (int i = 0; i < length; i++) {
                    final int cycle = rc ? length - i : i + 1;
                    secondReadTotalsByCycle[baseToInt(bases[i])][cycle] += 1;
                    secondReadCountsByCycle[cycle] += 1;
                }
            } else {
                for (int i = 0; i < length; i++) {
                    final int cycle = rc ? length - i : i + 1;
                    firstReadTotalsByCycle[baseToInt(bases[i])][cycle] += 1;
                    firstReadCountsByCycle[cycle] += 1;
                }
            }
        }

        private void ensureArraysBigEnough(final int length) {
            if (length > maxLengthSoFar) {
                for (int i = 0; i < 5; i++) {
                    firstReadTotalsByCycle[i] = Arrays.copyOf(firstReadTotalsByCycle[i], length);
                    secondReadTotalsByCycle[i] = Arrays.copyOf(secondReadTotalsByCycle[i], length);
                }
                firstReadCountsByCycle = Arrays.copyOf(firstReadCountsByCycle, length);
                secondReadCountsByCycle = Arrays.copyOf(secondReadCountsByCycle, length);
                maxLengthSoFar = length;
            }
        }

        boolean isEmpty() {
            return maxLengthSoFar == 0;
        }

        public void addToMetricsFile(final MetricsFile<BaseDistributionByCycleMetrics, ?> metrics) {
            int firstReadLength = 0;
            for (int i = 0; i < maxLengthSoFar; i++) {
                if (0 != firstReadCountsByCycle[i]) {
                    final BaseDistributionByCycleMetrics metric = new BaseDistributionByCycleMetrics();
                    metric.READ_END = 1;
                    metric.CYCLE = i;
                    metric.PCT_A = (100.0 * firstReadTotalsByCycle[0][i] / firstReadCountsByCycle[i]);
                    metric.PCT_C = (100.0 * firstReadTotalsByCycle[1][i] / firstReadCountsByCycle[i]);
                    metric.PCT_G = (100.0 * firstReadTotalsByCycle[2][i] / firstReadCountsByCycle[i]);
                    metric.PCT_T = (100.0 * firstReadTotalsByCycle[3][i] / firstReadCountsByCycle[i]);
                    metric.PCT_N = (100.0 * firstReadTotalsByCycle[4][i] / firstReadCountsByCycle[i]);
                    metrics.addMetric(metric);
                    firstReadLength = i;
                }
            }
            if (seenSecondEnd) {
                for (int i = 0; i < maxLengthSoFar; i++) {
                    if (0 != secondReadCountsByCycle[i]) {
                        final BaseDistributionByCycleMetrics metric = new BaseDistributionByCycleMetrics();
                        metric.READ_END = 2;
                        metric.CYCLE = (i + firstReadLength);
                        metric.PCT_A = (100.0 * secondReadTotalsByCycle[0][i] / secondReadCountsByCycle[i]);
                        metric.PCT_C = (100.0 * secondReadTotalsByCycle[1][i] / secondReadCountsByCycle[i]);
                        metric.PCT_G = (100.0 * secondReadTotalsByCycle[2][i] / secondReadCountsByCycle[i]);
                        metric.PCT_T = (100.0 * secondReadTotalsByCycle[3][i] / secondReadCountsByCycle[i]);
                        metric.PCT_N = (100.0 * secondReadTotalsByCycle[4][i] / secondReadCountsByCycle[i]);
                        metrics.addMetric(metric);
                    }
                }
            }
        }
    }
}
