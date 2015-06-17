package picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.PositionalArguments;
import picard.cmdline.programgroups.Metrics;

import java.io.File;
import java.io.FileReader;
import java.util.List;

/**
 * Compare two metrics files.
 */
@CommandLineProgramProperties(
        usage = CompareMetrics.USAGE_SUMMARY + CompareMetrics.USAGE_DETAIL,
        usageShort = CompareMetrics.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class CompareMetrics extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Compares two metrics files";
    static final String USAGE_DETAIL = "This tool compares the metrics and histograms generated from metric tools to determine " +
            "if the generated results are identical.  This tool is useful to test and compare outputs when code changes are implemented.<br /><br />  " +
            "Output indicates that two metrics are either equal or not equal. <br /> "  +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar CompareMetrics \\<br />" +
            "      mymetricfile1.txt \\<br />" +
            "      mymetricfile2.txt" +
            "</pre>" +
            "<hr />";

    @PositionalArguments(minElements = 2, maxElements = 2)
    public List<File> metricsFiles;

    private static final Log log = Log.getInstance(CompareMetrics.class);

    @Override
    protected int doWork() {
        IOUtil.assertFilesAreReadable(metricsFiles);
        final MetricsFile<?, ?> metricsA = new MetricsFile();
        final MetricsFile<?, ?> metricsB = new MetricsFile();
        try {
            metricsA.read(new FileReader(metricsFiles.get(0)));
            metricsB.read(new FileReader(metricsFiles.get(1)));
            final boolean areEqual = metricsA.areMetricsEqual(metricsB) && metricsA.areHistogramsEqual(metricsB);
            final String status = areEqual ? "EQUAL" : "NOT EQUAL";
            log.info("Files " + metricsFiles.get(0) + " and " + metricsFiles.get(1) + "are " + status);
        } catch (final Exception e) {
            throw new PicardException(e.getMessage());
        }
        return 0;
    }

    public static void main(String[] argv) {
        new CompareMetrics().instanceMainWithExit(argv);
    }
}
