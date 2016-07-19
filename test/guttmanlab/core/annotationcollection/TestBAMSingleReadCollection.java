package guttmanlab.core.annotationcollection;

import static org.junit.Assert.assertEquals;
import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.BEDFileRecord;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.MappedFragment;
import guttmanlab.core.annotation.PairedMappedFragment;
import guttmanlab.core.annotation.PopulatedWindow;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotation.predicate.MaximumLengthFilter;
import guttmanlab.core.coordinatespace.CoordinateSpace;
import guttmanlab.core.coordinatespace.GenomeSize;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.util.CloseableIterator;

import org.junit.BeforeClass;
import org.junit.Test;

public class TestBAMSingleReadCollection {
	
	private static SingleInterval malat1 = new SingleInterval("chr19", 5795689, 5802671, Strand.BOTH, "Malat1");
	private static File singleBam = new File("resources/SingleCollectionTest.bam");
	private static CoordinateSpace refSpace;
	private static AnnotationCollection<BEDFileRecord> refSeqFeatures;	
	private static File chr19bam = new File("resources/chr19.clean.sorted.bam");

	@BeforeClass
	public static void setUp() throws IOException
	{
		BAMSingleReadCollection bam = new BAMSingleReadCollection(chr19bam);
		SAMFileHeader fhead = bam.getFileHeader(); 
		refSpace = new CoordinateSpace(fhead);  
		BEDFileIO io =  new BEDFileIO(new CoordinateSpace(GenomeSize.MM9));
		refSeqFeatures = io.loadFromFile(new File("resources/RefSeq_mm9.bed"));
	}

	@Test
	public void testSingleNumOverlappersPartial() throws IOException {
		AbstractAnnotationCollection<? extends MappedFragment> data = BAMFragmentCollectionFactory.createFromBam(singleBam);
		assertEquals(4, data.numOverlappers(malat1, false));
	}
	
	@Test
	public void testSingleNumOverlappersFull() throws IOException {
		AbstractAnnotationCollection<? extends MappedFragment> data = BAMFragmentCollectionFactory.createFromBam(singleBam);
		assertEquals(2, data.numOverlappers(malat1, true));
	}

	//Test a single exon positive gene with only three reads
	@Test //Pass
	public void ccdc87GetWindowReadCounts() throws IOException{
		CloseableIterator<BEDFileRecord> iter = refSeqFeatures.sortedIterator();
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
		
		BAMSingleReadCollection bam = new BAMSingleReadCollection(chr19bam);
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
		
		assertEquals("Window read counts for Ccdc87",true,map1.equals(expectedMap));
	}
	
	//Test a single exon negative gene 
	@Test //Pass
	public void neat1GetWindowReadCount() throws IOException{
		CloseableIterator<BEDFileRecord> iter = refSeqFeatures.sortedIterator();
		Map<Integer,Integer> map1 = new HashMap<Integer,Integer>();
		Annotation a = null;
		while(iter.hasNext())
		{
			a = iter.next();
			if(a.getName().equals("NR_003513"))
				break;
		}
		
		iter.close();
		BAMSingleReadCollection bam = new BAMSingleReadCollection(chr19bam);
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
		
		assertEquals("shoud return 3226 nonempty windows of size 50",3226,winCount);
	}
	
	@Test
	//Verifies the correct reads are returned for CCdc87, a multi exon negative strand gene.
	public void neat1ConvertCoordinatesGetWindows() throws IOException{
		CloseableIterator<BEDFileRecord> iter = refSeqFeatures.sortedIterator();
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
		BAMSingleReadCollection bam = new BAMSingleReadCollection(chr19bam);
		AnnotationCollection<DerivedAnnotation<SAMFragment>> converted = refSeqFeatures.convertCoordinates(bam, refSpace, true);
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
		
		assertEquals("converted coordinates windows",3177,winCount);
	}
	
	@Test
	public void simpleWIndowsTest() {
		SingleInterval fifty = new SingleInterval("chr19",5843330,5843380,Strand.BOTH);
		int count = 0;
		BAMSingleReadCollection bam = new BAMSingleReadCollection(chr19bam);
		CloseableIterator<? extends PopulatedWindow<SAMFragment>> windows = bam.getPopulatedWindows(fifty, 1);
		while(windows.hasNext())
		{
			@SuppressWarnings("unused")
			PopulatedWindow<SAMFragment> win = windows.next();
			count++;
		}
		windows.close();
		assertEquals("there should be 50 windows.",50,count);
	}
	
	@Test
	public void missingWindowsTest1() {
		Annotation a = new SingleInterval("chr19",5845800,5847200,Strand.BOTH);
		BAMSingleReadCollection bam = new BAMSingleReadCollection(chr19bam);
		CloseableIterator<? extends PopulatedWindow<SAMFragment>> windows = bam.getPopulatedWindows(a, 1);
		int count = 0;
		while(windows.hasNext())
		{
			windows.next();
			count++;
		}
		windows.close();
		assertEquals("nonzero windows = ?",435,count);
		
	}
	
	@Test
	public void missingWindowsTest2() {
		Annotation a = new SingleInterval("chr19",5845800,5847200,Strand.BOTH);
		BAMSingleReadCollection bam = new BAMSingleReadCollection(chr19bam);
		
		CloseableIterator<? extends PopulatedWindow<SAMFragment>> all_windows = bam.getPopulatedWindows(a, 1, 1, true);
		int count = 0;
		while(all_windows.hasNext())
		{
			@SuppressWarnings("unused")
			PopulatedWindow<SAMFragment> win = all_windows.next();
			count++;
		}
		
		assertEquals("all windows = 1400",1400,count);
	}

	//Tests that the sortedIterator returns reads that only overlap blocks, and not introns
	@Test
	public void sortedIteratorBlockTest1() {
		
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
		
		BAMSingleReadCollection bam = new BAMSingleReadCollection(chr19bam);

		CloseableIterator<SAMFragment> multi_iter = bam.sortedIterator(multi, false);
		
		int mcount = 0;
		while(multi_iter.hasNext())
		{
			multi_iter.next();
			mcount++;
		}
		multi_iter.close();
		
		assertEquals("mcount = 0.", 0,mcount);
	}
	
	//Tests that the sortedIterator returns reads that only overlap blocks, and not introns
	@Test
	public void sortedIteratorBlockTest2() {
		
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
		
		BAMSingleReadCollection bam = new BAMSingleReadCollection(chr19bam);
		
		CloseableIterator<SAMFragment> singl_iter = bam.sortedIterator(single, false);

		int scount = 0;
		while(singl_iter.hasNext())
		{
			singl_iter.next();
			scount++;
		}
		singl_iter.close();
		
		assertEquals("scount = 2.", 2,scount);
	}
	
	@Test
	public void iteratorStrandMatchingTest() throws IOException{
		CloseableIterator<BEDFileRecord> iter = refSeqFeatures.sortedIterator();

		Annotation a = null;
		while(iter.hasNext()) 
		{
			a = iter.next();
			//b.setOrientation(Strand.BOTH);
			if(a.getName().equals("NM_025741")) //Trpd52l3
				break;
		}
		iter.close();
		
		BAMSingleReadCollection bam = new BAMSingleReadCollection(chr19bam);

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
	public void multipleSortedIteratorsOnSameCollection1() throws IOException{

		
		BAMSingleReadCollection bam = new BAMSingleReadCollection(chr19bam);

		int count =0;		
		Annotation a2 = new SingleInterval("chr19", 30267000, 30272000, Strand.NEGATIVE);
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
	
	@Test
	public void multipleSortedIteratorsOnSameCollection2() throws IOException{

		Annotation a = new SingleInterval("chr19", 30267000, 30272000, Strand.POSITIVE);
		
		BAMSingleReadCollection bam = new BAMSingleReadCollection(chr19bam);

		int count =0;
		CloseableIterator<SAMFragment> f_iter = bam.sortedIterator(a, false);
		while(f_iter.hasNext())
		{
			@SuppressWarnings("unused")
			SAMFragment f = f_iter.next();
			count++;
		}
		
		f_iter.close();
		assertEquals("2 positive reads should overlap region.",2,count);
		
		
	}
	
	
	@Test
	public void annotationCollectionGetCount() {
		BAMSingleReadCollection bam = new BAMSingleReadCollection(chr19bam);
		int count = bam.getNumAnnotations();
		int count2 = bam.getNumAnnotations();
		bam.addFilter(new MaximumLengthFilter<SAMFragment>(1000000));
		int count3 = bam.getNumAnnotations();
		
		assertEquals(count,2267045);
		assertEquals(count2,2267045);
		assertEquals(count3,2267045);
	}
	
	@Test //Pass 
	//Verifies that the correct number of reads are returned for Ccdc87, a single exon gene on the positive strand.
	public void ccdc87OverlapReadCount() throws IOException{
		System.out.println("\n\nCcdc87 Mapped Reads:");
		CloseableIterator<BEDFileRecord> iter = refSeqFeatures.sortedIterator();

		//locate the feature of interest from the feature iterator
		Annotation a = null;
		while(iter.hasNext()) 
		{
			a = iter.next();
			if(a.getName().equals("NM_207268"))
				break;
		}
		iter.close();
		
		int count =0;
		BAMSingleReadCollection bam = new BAMSingleReadCollection(chr19bam);
		CloseableIterator<SAMFragment> f_iter = bam.sortedIterator(a, false);
		while(f_iter.hasNext())
		{
			SAMFragment f = f_iter.next();
			System.out.println(f.toString());
			count++;
		}
		
		f_iter.close();
		assertEquals("43 unconverted reads should overlap Ccdc87.",43,count); 
		
	}
	
	
	@Test //Pass
	//Verifies the correct number of overlapping reads are returned after conversion to feature space for Cd248.
	public void cd248ConvertCoordinates() throws IOException{
		System.out.println("\n\nCd248 Converted Reads:");
		CloseableIterator<BEDFileRecord> iter = refSeqFeatures.sortedIterator();

		Annotation a = null;
		while(iter.hasNext()) 
		{
			a = iter.next();
			if(a.getName().equals("NM_054042")) //Cd248
				break;
		}
		iter.close();
		
		Annotation b = new SingleInterval(a.getName(), 0, a.size()-1);
		b.setOrientation(Strand.NEGATIVE);
		BAMSingleReadCollection bam = new BAMSingleReadCollection(chr19bam);
		AnnotationCollection<DerivedAnnotation<SAMFragment>> converted = refSeqFeatures.convertCoordinates(bam, refSpace, false);	
		CloseableIterator<DerivedAnnotation<SAMFragment>> c_iter = converted.sortedIterator(b,false);
		
		int count = 0;
		while(c_iter.hasNext())
		{
			DerivedAnnotation<SAMFragment> c = c_iter.next();
			System.out.println(c.toString());
			count++;
		}
		
		assertEquals("6 converted annotations should overlap Cd248",6,count);
	}
	
	
	@Test //Pass
	//Should fail; region is treated as start and end, should not include intron overlaps
	public void kcnk4Overlaps() throws IOException{
		System.out.println("\n\nKcnk4 Mapped Reads:");
		CloseableIterator<BEDFileRecord> iter = refSeqFeatures.sortedIterator();

		Annotation a = null;
		while(iter.hasNext()) 
		{
			a = iter.next();
			if(a.getName().equals("NM_008431"))
				break;
		}
		iter.close();
		BAMSingleReadCollection bam = new BAMSingleReadCollection(chr19bam);
		CloseableIterator<SAMFragment> c_iter = bam.sortedIterator(a,true);
		
		int count = 0;
		while(c_iter.hasNext())
		{
			SAMFragment c = c_iter.next();
			//System.out.println(c.getParentAnnotation().toString());
			System.out.println(c.toString());
			count++;
		}
		
		assertEquals("10 converted annotations should overlap Kcnk4",10,count);
	}
	
	
	@Test //Pass
	//Verifies the correct reads are returned for CCdc87, a multi exon negative strand gene.
	public void kcnk4ConvertCoordinates() throws IOException{
		System.out.println("\n\nKcnk4 Converted Reads:");
		CloseableIterator<BEDFileRecord> iter = refSeqFeatures.sortedIterator();

		Annotation a = null;
		while(iter.hasNext()) 
		{
			a = iter.next();
			if(a.getName().equals("NM_008431"))
				break;
		}
		iter.close();
		
		Annotation b = new SingleInterval(a.getName(), 0, a.size()-1);
		b.setOrientation(Strand.BOTH);
		BAMSingleReadCollection bam = new BAMSingleReadCollection(chr19bam);
		AnnotationCollection<DerivedAnnotation<SAMFragment>> converted = refSeqFeatures.convertCoordinates(bam, refSpace, true);	
		CloseableIterator<DerivedAnnotation<SAMFragment>> c_iter = converted.sortedIterator(b,true);
		
		int count = 0;
		while(c_iter.hasNext())
		{
			DerivedAnnotation<SAMFragment> c = c_iter.next();
			//System.out.println(c.getParentAnnotation().toString());
			System.out.println(c.toString());
			count++;
		}
		
		assertEquals("6 converted annotations should overlap Kcnk4",2,count);
	}
	
	
	@Test 
	//Test a multi-exon negative gene for correct number of unconverted overlaps
	public void dkk1OverlapReadCount()
	{
		System.out.println("\n\nDkk1 Mapped Reads:");
		CloseableIterator<BEDFileRecord> iter = refSeqFeatures.sortedIterator();

		Annotation a = null;
		while(iter.hasNext()) 
		{
			a = iter.next();
			if(a.getName().equals("NM_010051")) //Dkk1
				break;
		}
		iter.close();
		
		int count =0;
		BAMSingleReadCollection bam = new BAMSingleReadCollection(chr19bam);
		bam.addFilter(new MaximumLengthFilter<SAMFragment>(50000));
		//bam.writeToBAM("/home/burghard/Desktop/dkk1splicederror", a, false);
		CloseableIterator<SAMFragment> f_iter = bam.sortedIterator(a, false);
		while(f_iter.hasNext())
		{
			SAMFragment f = f_iter.next();
			System.out.println(f.toString());
			count++;
		}
		
		f_iter.close();
		assertEquals("77 unconverted reads should overlap Dkk1.",77,count); 
	}
	
	
	@Test //Fail no reads
	//Test a multi-exon negative gene for correct number of converted, fully contained overlaps
	public void dkk1ConvertCoordinatesFullyContained() throws IOException{
		System.out.println("\n\nDkk1 Converted Reads (Fully Contained):");
		CloseableIterator<BEDFileRecord> iter = refSeqFeatures.sortedIterator();

		Annotation a = null;
		while(iter.hasNext()) 
		{
			a = iter.next();
			if(a.getName().equals("NM_010051")) //Dkk1
				break;
		}
		iter.close();
		
		Annotation b = new SingleInterval(a.getName(), 0, a.size()-1);
		b.setOrientation(Strand.NEGATIVE);
		BAMSingleReadCollection bam = new BAMSingleReadCollection(chr19bam);
		AnnotationCollection<DerivedAnnotation<SAMFragment>> converted = refSeqFeatures.convertCoordinates(bam, refSpace, true);		
		CloseableIterator<DerivedAnnotation<SAMFragment>> c_iter = converted.sortedIterator(b,true);

		int count = 0;
		while(c_iter.hasNext())
		{
			DerivedAnnotation<SAMFragment> c = c_iter.next();
			//System.out.println(c.getParentAnnotation().toString());
			System.out.println(c.toString());
			count++;
		}
		
		assertEquals("72 converted annotations should overlap Dkk1",72,count);
	}
	
	@Test //Pass
	//Test a multi-exon negative gene for correct number of converted overlaps
	public void dkk1ConvertCoordinates() throws IOException{
		System.out.println("\n\nDkk1 Converted Reads:");
		CloseableIterator<BEDFileRecord> iter = refSeqFeatures.sortedIterator();

		Annotation a = null;
		while(iter.hasNext()) 
		{
			a = iter.next();
			if(a.getName().equals("NM_010051"))
				break;
		}
		iter.close();
		
		Annotation b = new SingleInterval(a.getName(), 0, a.size()-1);
		b.setOrientation(Strand.BOTH);
		System.out.println("Orientation: "+b.getOrientation());
		BAMSingleReadCollection bam = new BAMSingleReadCollection(chr19bam);
		AnnotationCollection<DerivedAnnotation<SAMFragment>> converted = refSeqFeatures.convertCoordinates(bam, refSpace, false);	
		//converted.writeToBAM("/home/burghard/Desktop/Dkk1Converted.bam");
		CloseableIterator<DerivedAnnotation<SAMFragment>> c_iter = converted.sortedIterator(b,false);
		
		int count = 0;
		while(c_iter.hasNext())
		{
			DerivedAnnotation<SAMFragment> c = c_iter.next();
			System.out.println(c.toString());
			count++;
		}
		
		assertEquals("75 converted annotations should overlap Dkk1",75,count);
	}
	
	
	
	//Utility Methods
	public ArrayList<String> iterToNameList(CloseableIterator<? extends Annotation> c_iter)
	{
		ArrayList<String> c_list = new ArrayList<String>();

		while(c_iter.hasNext())
			c_list.add(c_iter.next().getName());
		c_iter.close();
		
		return c_list;
	}
	
	public void compareOutputs(ArrayList<String> a_list, ArrayList<String> b_list)
	{ 
		int a_count = a_list.size(), b_count = b_list.size();
		//System.out.println("\n\nTotal in A: " + a_count + "\tTotal in B: " + b_count);
		
		//System.out.print("\nMismatched Annotations only in B: \n");
		for(String name : b_list)
		{
			if(!a_list.contains(name))
			{
				//System.out.println(name);
			}
			else
				a_list.remove(name);
		}
		
		//System.out.print("\n\nMismatched Annotations only in A: \n");
		for(String name : a_list)
			{
			//System.out.println(name);
			}	
	}


}
