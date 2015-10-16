package guttmanlab.core.test;

import static org.junit.Assert.*;
import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import guttmanlab.core.coordinatespace.CoordinateSpace;

import java.io.File;
import java.io.IOException;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.util.CloseableIterator;

import org.junit.Before;
import org.junit.Test;

public class BAMWriteTest {

	private BAMSingleReadCollection bam;
	private SAMFileHeader fhead;
	private CoordinateSpace refSpace;
	private String fname;
	private BEDFileIO io;
	private AnnotationCollection<Gene> features;
	
	@Before
	public void setUp() throws IOException
	{
		this.bam = new BAMSingleReadCollection(new File("/Users/cburghard/Downloads/chr19.clean.sorted.bam"));
		this.fhead = bam.getFileHeader(); 
		this.refSpace = new CoordinateSpace(fhead);  
		this.fname = "/Users/cburghard/Downloads/RefSeqStrandTest.bed";
		this.io =  new BEDFileIO("/Users/cburghard/Downloads/refspace.txt"); 
		this.features = io.loadFromFile(fname);
	}
	

	@Test
	public void BAMReadWritetest() {
		String fname = "/Users/cburghard/Downloads/newGeneTest.bam";
		File f = new File(fname);
		CloseableIterator<Gene> iter = features.sortedIterator();
		Annotation a = null;
		
		while(iter.hasNext()) 
		{
			a = iter.next();
			if(a.getName().equals("NM_025741"))
				break;
		}
		iter.close();
		
		bam.writeToBAM(fname, a, false);
		BAMSingleReadCollection bamRead = new BAMSingleReadCollection(f);
		CloseableIterator<SAMFragment> b_iter = bamRead.sortedIterator();
		
		while(b_iter.hasNext())
		{
			System.out.println(b_iter.next().toString());
		}
		
		fail("Not yet implemented");
	}

}
