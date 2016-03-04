package guttmanlab.core.test;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.PopulatedWindow;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import guttmanlab.core.coordinatespace.CoordinateSpace;
import junit.framework.TestCase;
import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.util.CloseableIterator;

import org.junit.Before;
import org.junit.Test;

public class GetWindowsTest {
	
	private BAMSingleReadCollection bam;
	private SAMFileHeader fhead;
	private CoordinateSpace refSpace;
	private String fname;
	private BEDFileIO io;
	private AnnotationCollection<Gene> features;	
	@Before
	public void setUp() throws IOException
	{
		this.bam = new BAMSingleReadCollection(new File("/storage/shared/CoreTestData/chr19.clean.sorted.bam"));
		this.fhead = bam.getFileHeader(); 
		this.refSpace = new CoordinateSpace(fhead);  
		this.fname = "/storage/shared/CoreTestData/RefSeq.bed";  
		this.io =  new BEDFileIO("/storage/shared/CoreTestData/refspace.txt"); //use existing coordinateSpace?
		this.features = io.loadFromFile(fname);
	}
	//test
	
	//Test a single exon positive gene with only three reads
	@Test //Pass
	public void Ccdc87GetWindowReadCounts() throws IOException{
		CloseableIterator<Gene> iter = features.sortedIterator();
		Map<Integer,Integer> map1 = new HashMap<Integer,Integer>();
		
		//expect 164 windows w/1 read, 188 with 2, and 12 with 3.
		Map<Integer,Integer> expectedMap = new HashMap<Integer,Integer>();
		expectedMap.put(1, 164);
		expectedMap.put(2, 188);
		expectedMap.put(3, 12);
		
		Annotation a = null;
		while(iter.hasNext())
		{
			a = iter.next();
			if(a.getName().equals("NM_025741"))
				break;
		}
		
		iter.close();
		
		CloseableIterator<? extends PopulatedWindow<SAMFragment>> windows = bam.getPopulatedWindows(a, 50);
		
		int count = 0;
		while(windows.hasNext())
		{
			PopulatedWindow<SAMFragment> win = windows.next();
			int reads = win.getNumberOfAnnotationsInWindow();

			if(map1.containsKey(reads))
				{
					count = map1.get(reads);
					map1.remove(reads);
					map1.put(reads, count+1);
				}
			else
				map1.put(reads, 1);
		}
		
		System.out.println(map1.toString());
		assertEquals("Window read counts for Ccdc87",true,map1.equals(expectedMap));
	}
	
	//Test a single exon negative gene 
	@Test //Pass
	public void Neat1GetWindowReadCount() throws IOException{
		CloseableIterator<Gene> iter = features.sortedIterator();
		Map<Integer,Integer> map1 = new HashMap<Integer,Integer>();
		Annotation a = null;
		while(iter.hasNext())
		{
			a = iter.next();
			if(a.getName().equals("NR_003513"))
				break;
		}
		
		iter.close();
		CloseableIterator<? extends PopulatedWindow<SAMFragment>> windows = bam.getPopulatedWindows(a, 50);
		
		int winCount = 0;
		while(windows.hasNext())
		{
			PopulatedWindow<SAMFragment> win = windows.next();
			winCount++;
			int reads = win.getNumberOfAnnotationsInWindow();
			int count = 0;
			if(map1.containsKey(reads))
				{
					count = map1.get(reads);
					map1.remove(reads);
					map1.put(reads, count+1);
				}
			else
				map1.put(reads, 1);
		}
		
		System.out.println(map1.toString());
		assertEquals("shoud return 3226 nonempty windows of size 50",3226,winCount);
	}
	
	@Test
	//Verifies the correct reads are returned for CCdc87, a multi exon negative strand gene.
	public void Neat1ConvertCoordinatesGetWindows() throws IOException{
		CloseableIterator<Gene> iter = features.sortedIterator();
		Map<Integer,Integer> map1 = new HashMap<Integer,Integer>();
		Annotation a = null;
		while(iter.hasNext()) 
		{
			a = iter.next();
			if(a.getName().equals("NR_003513"))
				break;
		}
		iter.close();
		
		Annotation b = new SingleInterval(a.getName(), 0, a.size()-1);
		b.setOrientation(Strand.BOTH);
		AnnotationCollection<DerivedAnnotation<SAMFragment>> converted = features.convertCoordinates(bam, refSpace, true);
		CloseableIterator<? extends PopulatedWindow<DerivedAnnotation<SAMFragment>>> windows = converted.getPopulatedWindows(b, 50);

		int winCount = 0;
		while(windows.hasNext())
		{
			PopulatedWindow<DerivedAnnotation<SAMFragment>> win = windows.next();
			winCount++;
			int reads = win.getNumberOfAnnotationsInWindow();
			int count = 0;
			if(map1.containsKey(reads))
				{
					count = map1.get(reads);
					map1.remove(reads);
					map1.put(reads, count+1);
				}
			else
				map1.put(reads, 1);
		}
		
		System.out.println(map1.toString());
		assertEquals("converted coordinates windows",3177,winCount);
	}
	
	@Test
	public void SimpleWIndowsTest() {
		SingleInterval fifty = new SingleInterval("chr19",5843330,5843380,Strand.BOTH);
		int count = 0;
		CloseableIterator<? extends PopulatedWindow<SAMFragment>> windows = bam.getPopulatedWindows(fifty, 1);
		while(windows.hasNext())
		{
			PopulatedWindow<SAMFragment> win = windows.next();
			System.out.println(win.toBED());
			count++;
		}
		windows.close();
		System.out.println(count);
		assertEquals("there should be 50 windows.",50,count);
	}
	
	@Test
	public void MissingWindowsTest() {
		Annotation a = new SingleInterval("chr19",5845800,5847200,Strand.BOTH);
		CloseableIterator<? extends PopulatedWindow<SAMFragment>> windows = bam.getPopulatedWindows(a, 1);
		int count = 0;
		while(windows.hasNext())
		{
			PopulatedWindow<SAMFragment> win = windows.next();
			count++;
		}
		windows.close();
		assertEquals("nonzero windows = ?",435,count);
		
		CloseableIterator<? extends PopulatedWindow<SAMFragment>> all_windows = bam.getPopulatedWindows(a, 1, 1, true);
		count = 0;
		while(all_windows.hasNext())
		{
			PopulatedWindow<SAMFragment> win = all_windows.next();
			System.out.println(win.toBED());
			count++;
		}
		
		assertEquals("all windows = 1400",1400,count);
	}
	
}
