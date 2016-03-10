package guttmanlab.core.test;

import static org.junit.Assert.*;
import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.BEDFileRecord;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.PopulatedWindow;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMPairedFragmentCollection;

import java.io.File;
import java.io.IOException;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.CloseableIterator;

import org.apache.log4j.Logger;
import org.junit.Before;
import org.junit.Test;

public class PopulatedWindowsTest {
	
	private BAMPairedFragmentCollection data;
	private String fname;
	private BEDFileIO io;
	private AnnotationCollection<BEDFileRecord> features;	
	private BEDFileRecord gene;
	private static Logger logger = Logger.getLogger(PopulatedWindowsTest.class.getName());
	
	@Before
	public void setUp() throws IOException
	{
		data = new BAMPairedFragmentCollection(new File("/Users/prussell/Documents/lncRNA/Peak_caller_testing/input.bam"));
		fname = "/Users/prussell/Documents/lncRNA/Peak_caller_testing/genes.bed";  
		io =  new BEDFileIO("/Users/prussell/Documents/lncRNA/Peak_caller_testing/sizes");
		features = io.loadFromFile(fname);
		CloseableIterator<BEDFileRecord> iter = features.sortedIterator();
		while(iter.hasNext()) {
			BEDFileRecord g = iter.next();
			if(g.getName().equals("NM_029432")) {
				gene = g;
				break;
			}
		}
		iter.close();
	}
	
	//@Test
	public void testNumFragmentsInFragmentBam() {
		logger.info("");
		File fragmentFile = new File("/Users/prussell/Documents/lncRNA/Peak_caller_testing/input.pe.bam");
		SAMFileReader reader = new SAMFileReader(fragmentFile);
		SAMRecordIterator iter = reader.iterator();
		int numFound = 0;
		while(iter.hasNext()) {
			SAMRecord record = iter.next();
			if(record.getReadName().equals("HISEQ:665:H1110ADXX:2:2108:13851:11667")) {
				logger.info("Fragment bam file contains read " + record.getSAMString());
				numFound++;
			}
		}
		reader.close();
		assertEquals(1.0, numFound, 0.01);
	}
	
	//@Test
	public void testNumFragmentsInOriginalBam() {
		logger.info("");
		File bamFile = new File("/Users/prussell/Documents/lncRNA/Peak_caller_testing/input.bam");
		SAMFileReader reader = new SAMFileReader(bamFile);
		SAMRecordIterator iter = reader.iterator();
		int numFound = 0;
		while(iter.hasNext()) {
			SAMRecord record = iter.next();
			if(record.getReadName().equals("HISEQ:665:H1110ADXX:2:2108:13851:11667")) {
				logger.info("Original bam file contains read " + record.getSAMString());
				numFound++;
			}
		}
		reader.close();
		assertEquals(2.0, numFound, 0.01);
	}
	
	@Test
	public void getPopulatedWindowsTest() {
		CloseableIterator<? extends PopulatedWindow<? extends Annotation>> windowIter = data.getPopulatedWindows(gene, 50, 10);
		PopulatedWindow<? extends Annotation> window = null;
		while(windowIter.hasNext()) {
			PopulatedWindow<? extends Annotation> w = windowIter.next();
			if(w.toUCSC().equals("chr2:130539080-130539130")) {
				window = w;
				break;
			}
		}
		windowIter.close();
		double n = window.getNumberOfAnnotationsInWindow();
		logger.info("Window " + window.toUCSC() + " contains " + n + " fragments.");
		assertEquals(1.0, n, 0.01);
	}
	
}
