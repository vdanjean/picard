package picard.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Iso8601Date;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

import java.io.File;
import java.util.Arrays;

/**
 * Replaces read groups in a BAM file
 *
 * @author mdepristo
 */
@CommandLineProgramProperties(
        usage = AddOrReplaceReadGroups.USAGE_SUMMARY + AddOrReplaceReadGroups.USAGE_DETAILS,
        usageShort = AddOrReplaceReadGroups.USAGE_SUMMARY,
        programGroup = SamOrBam.class
)
public class AddOrReplaceReadGroups extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Replaces read groups in a BAM file.  ";
    static final String USAGE_DETAILS = "This tool enables the user to merge all the read group in the INPUT files into a single " +
            "new read group.  The definition of a read group can vary depending on the sequencing platform used.  " +
            "For example, a read group refers to an instrument lane for the Illumina platform and for SOLiD, a read group indicates a slide.  " +
            "Typically, an Illumina sequencer will produce up to eight read groups of data, corresponding to the 8 lanes of the instrument.  " +
            "Merging read groups enables the user to combine data from multiple lanes or runs." +
            "<br /><br /> " +
            "This tool accepts INPUT \".bam\" files or URLs from the Global Alliance for Genomics and Health (GA4GH)*.  " +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar AddOrReplaceReadGroups \\<br />" +
            "      I=bamfile_1.bam \\<br />" +
            "      I=bamfile_2.bam \\<br />" +
            "      O=bamfile_1_2.bam \\<br />" +
            "      RGID=4 \\<br />" +
            "      RGLB=lib1 \\<br />" +
            "      RGPL=illumina \\<br />" +
            "      RGPU=unit1 \\<br />" +
            "      RGSM=20" +
            "</pre>" +
            "*For information on GA4GH, please see: http://ga4gh.org/#/documentation. " +
            "<hr />" ;
    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input file (bam or sam or a GA4GH url).")
    public String INPUT = null;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (bam or sam).")
    public File OUTPUT = null;

    @Option(shortName = StandardOptionDefinitions.SORT_ORDER_SHORT_NAME, optional = true,
            doc = "Optional sort order to output in. If not supplied OUTPUT is in the same order as INPUT.")
    public SortOrder SORT_ORDER;

    @Option(shortName = "ID", doc = "Read Group ID")
    public String RGID = "1";

    @Option(shortName = "LB", doc = "Read Group Library")
    public String RGLB;

    @Option(shortName = "PL", doc = "Read Group platform (e.g. illumina, solid)")
    public String RGPL;

    @Option(shortName = "PU", doc = "Read Group platform unit (eg. run barcode)")
    public String RGPU;

    @Option(shortName = "SM", doc = "Read Group sample name")
    public String RGSM;

    @Option(shortName = "CN", doc = "Read Group sequencing center name", optional = true)
    public String RGCN;

    @Option(shortName = "DS", doc = "Read Group description", optional = true)
    public String RGDS;

    @Option(shortName = "DT", doc = "Read Group run date", optional = true)
    public Iso8601Date RGDT;

    @Option(shortName = "PI", doc = "Read Group predicted insert size", optional = true)
    public Integer RGPI;
    
    @Option(shortName = "PG", doc = "Read Group program group", optional = true)
    public String RGPG;
    
    @Option(shortName = "PM", doc = "Read Group platform model", optional = true)
    public String RGPM;

    private final Log log = Log.getInstance(AddOrReplaceReadGroups.class);

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        new AddOrReplaceReadGroups().instanceMainWithExit(argv);
    }

    protected int doWork() {
        IOUtil.assertInputIsValid(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        final SamReader in = SamReaderFactory.makeDefault()
            .referenceSequence(REFERENCE_SEQUENCE)
            .open(SamInputResource.of(INPUT));

        // create the read group we'll be using
        final SAMReadGroupRecord rg = new SAMReadGroupRecord(RGID);
        rg.setLibrary(RGLB);
        rg.setPlatform(RGPL);
        rg.setSample(RGSM);
        rg.setPlatformUnit(RGPU);
        if (RGCN != null) rg.setSequencingCenter(RGCN);
        if (RGDS != null) rg.setDescription(RGDS);
        if (RGDT != null) rg.setRunDate(RGDT);
        if (RGPI != null) rg.setPredictedMedianInsertSize(RGPI);
        if (RGPG != null) rg.setProgramGroup(RGPG);
        if (RGPM != null) rg.setPlatformModel(RGPM);

        log.info(String.format("Created read group ID=%s PL=%s LB=%s SM=%s%n", rg.getId(), rg.getPlatform(), rg.getLibrary(), rg.getSample()));

        // create the new header and output file
        final SAMFileHeader inHeader = in.getFileHeader();
        final SAMFileHeader outHeader = inHeader.clone();
        outHeader.setReadGroups(Arrays.asList(rg));
        if (SORT_ORDER != null) outHeader.setSortOrder(SORT_ORDER);

        final SAMFileWriter outWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(outHeader,
                outHeader.getSortOrder() == inHeader.getSortOrder(),
                OUTPUT);

        final ProgressLogger progress = new ProgressLogger(log);
        for (final SAMRecord read : in) {
            read.setAttribute(SAMTag.RG.name(), RGID);
            outWriter.addAlignment(read);
            progress.record(read);
        }

        // cleanup
        CloserUtil.close(in);
        outWriter.close();
        return 0;
    }
}
