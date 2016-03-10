package guttmanlab.core.test;

import static org.junit.Assert.*;
import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.BEDFileRecord;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotation.predicate.MaximumLengthFilter;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;

import java.io.File;
import java.io.IOException;

import net.sf.samtools.util.CloseableIterator;

import org.junit.Before;
import org.junit.Test;

public class BAMSingleReadTest {

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

	//@Test //Pass
	//Tests that the sortedIterator returns reads that only overlap blocks, and not introns
	public void SortedIteratorBlockTest() {
		
		BlockedAnnotation multi = new BlockedAnnotation();
		BlockedAnnotation single = new  BlockedAnnotation();
		SingleInterval block1 = new SingleInterval("chr19", 30090800, 30090948);
		SingleInterval block2 = new SingleInterval("chr19", 30091691, 30091891);
		SingleInterval block = new SingleInterval("chr19", 30090800, 30091891);
		
		multi.addBlocks(block1);
		multi.addBlocks(block2);
		single.addBlocks(block);
		
		multi.setOrientation(Strand.BOTH);
		single.setOrientation(Strand.BOTH);
		
		CloseableIterator<SAMFragment> multi_iter = bam.sortedIterator(multi, false);
		
		int mcount = 0;
		while(multi_iter.hasNext())
		{
			multi_iter.next();
			mcount++;
		}
		multi_iter.close();
		
		CloseableIterator<SAMFragment> singl_iter = bam.sortedIterator(single, false);

		int scount = 0;
		while(singl_iter.hasNext())
		{
			singl_iter.next();
			scount++;
		}
		singl_iter.close();
		
		System.out.println("mcount: " + mcount + "\nscount: " + scount);
		assertEquals("mcount = 0.", 0,mcount);
		assertEquals("scount = 2.", 2,scount);
	}
	

	//@Test
	public void IteratorStrandMatchingTest() throws IOException{
		//System.out.println("\n\nCcdc87 Mapped Reads:");
		CloseableIterator<BEDFileRecord> iter = features.sortedIterator();

		Annotation a = null;
		while(iter.hasNext()) 
		{
			a = iter.next();
			//b.setOrientation(Strand.BOTH);
			if(a.getName().equals("NM_025741")) //Trpd52l3
				break;
		}
		iter.close();
		
		int count =0;
		CloseableIterator<SAMFragment> f_iter = bam.sortedIterator(a, false);
		while(f_iter.hasNext())
		{
			f_iter.next();
			count++;
		}
		
		f_iter.close();
		assertEquals("3 positive read should overlap Trpd52l3.",3,count); 
		
	}
	
	@Test
	public void MultipleSortedIteratorsOnSameCollection() throws IOException{

		Annotation a = new SingleInterval("chr19", 30267000, 30272000, Strand.POSITIVE);
		
		int count =0;
		CloseableIterator<SAMFragment> f_iter = bam.sortedIterator(a, false);
		while(f_iter.hasNext())
		{
			SAMFragment f = f_iter.next();
			System.out.println(f.toBED());
			count++;
		}
		
		f_iter.close();
		assertEquals("2 positive reads should overlap region.",2,count);
		
		Annotation a2 = new SingleInterval("chr19", 30267000, 30272000, Strand.NEGATIVE);
		System.out.println(a.overlaps(a2));
		count =0;
		CloseableIterator<SAMFragment> f_iter2 = bam.sortedIterator(a2, false);
		while(f_iter2.hasNext())
		{
			f_iter2.next();
			count++;
		}
		
		f_iter2.close();
		assertEquals("6 negative reads should overlap region.",6,count); 
		
	}
	
	
	//@Test
	public void AnnotationCollectionGetCount() {
		int count = bam.getNumAnnotations();
		System.out.println(count);
		int count2 = bam.getNumAnnotations();
		bam.addFilter(new MaximumLengthFilter<SAMFragment>(1000000));
		System.out.println(count2);
		int count3 = bam.getNumAnnotations();
		System.out.println(count3);
		
		assertEquals(count,2267045);
		assertEquals(count2,2267045);
		assertEquals(count3,2267045);
	}
	

}
