package picard.sam;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.DownsamplingIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

import java.io.File;
import java.util.Random;

/**
 * Class to randomly downsample a BAM file while respecting that we should either get rid
 * of both ends of a pair or neither end of the pair!
 */
@CommandLineProgramProperties(
        usage = "Randomly down-sample a SAM or BAM file to retain " +
                "a random subset of the reads. Mate-pairs are either both kept or both discarded. Reads marked as not primary " +
                "alignments are all discarded. Each read is given a probability P of being retained - results with the exact " +
                "same input in the same order and with the same value for RANDOM_SEED will produce the same results.",
        usageShort = "Down-sample a SAM or BAM file to retain a random subset of the reads",
        programGroup = SamOrBam.class
)
public class DownsampleSam extends CommandLineProgram {

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to downsample.")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output, downsampled, SAM or BAM file to write.")
    public File OUTPUT;

    @Option(shortName = "R", doc = "Random seed to use if reproducibilty is desired.  " +
            "Setting to null will cause multiple invocations to produce different results.")
    public Long RANDOM_SEED = 1L;

    @Option(shortName = "P", doc = "The probability of keeping any individual read, between 0 and 1.")
    public double PROBABILITY = 1;

    private final Log log = Log.getInstance(DownsampleSam.class);

    public static void main(final String[] args) {
        new DownsampleSam().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        final Random r = RANDOM_SEED == null ? new Random() : new Random(RANDOM_SEED);
        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        final SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(in.getFileHeader(), true, OUTPUT);
        final ProgressLogger progress = new ProgressLogger(log, (int) 1e7, "Wrote");
        final DownsamplingIterator iterator = new DownsamplingIterator(in.iterator(), r, PROBABILITY);

        for (final SAMRecord rec : iterator) {
            out.addAlignment(rec);
            progress.record(rec);
        }

        out.close();
        CloserUtil.close(in);
        log.info("Finished! Kept " + iterator.getKeptReads() + " out of " + iterator.getTotalReads() + " reads.");

        return 0;
    }
}
