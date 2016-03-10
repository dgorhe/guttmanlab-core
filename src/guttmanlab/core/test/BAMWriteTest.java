package guttmanlab.core.test;

import static org.junit.Assert.*;
import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.BEDFileRecord;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;

import java.io.File;
import java.io.IOException;

import net.sf.samtools.util.CloseableIterator;

import org.junit.Before;
import org.junit.Test;

public class BAMWriteTest {

	private BAMSingleReadCollection bam;
	private String fname;
	private BEDFileIO io;
	private AnnotationCollection<BEDFileRecord> features;
	
	@Before
	public void setUp() throws IOException
	{
		this.bam = new BAMSingleReadCollection(new File("/storage/shared/CoreTestData/chr19.clean.sorted.bam"));
		this.fname = "/storage/shared/CoreTestData/RefSeqStrandTest.bed";
		this.io =  new BEDFileIO("/storage/shared/CoreTestData/refspace.txt"); 
		this.features = io.loadFromFile(fname);
	}
	

	@Test
	public void BAMReadWritetest() {
		String fname = "/storage/shared/CoreTestData/newGeneTest.bam";
		File f = new File(fname);
		CloseableIterator<BEDFileRecord> iter = features.sortedIterator();
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
