package guttmanlab.core.test;

import static org.junit.Assert.*;
import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.PairedMappedFragment;
import guttmanlab.core.annotation.PopulatedWindow;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotation.predicate.MaximumLengthFilter;
import guttmanlab.core.annotation.predicate.MinimumLengthFilter;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMPairedFragmentCollection;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import guttmanlab.core.coordinatespace.CoordinateSpace;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.util.CloseableIterator;

import org.junit.Before;
import org.junit.Test;

public class ConvertCoordinatesTest {
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
	
	@Test //Pass 
	//Verifies that the correct number of reads are returned for Ccdc87, a single exon gene on the positive strand.
	public void Ccdc87OverlapReadCount() throws IOException{
		System.out.println("\n\nCcdc87 Mapped Reads:");
		CloseableIterator<Gene> iter = features.sortedIterator();

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
	public void Cd248ConvertCoordinates() throws IOException{
		System.out.println("\n\nCd248 Converted Reads:");
		CloseableIterator<Gene> iter = features.sortedIterator();

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
		AnnotationCollection<DerivedAnnotation<SAMFragment>> converted = features.convertCoordinates(bam, refSpace, false);	
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
	public void Kcnk4Overlaps() throws IOException{
		System.out.println("\n\nKcnk4 Mapped Reads:");
		CloseableIterator<Gene> iter = features.sortedIterator();

		Annotation a = null;
		while(iter.hasNext()) 
		{
			a = iter.next();
			if(a.getName().equals("NM_008431"))
				break;
		}
		iter.close();
		
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
	public void Kcnk4ConvertCoordinates() throws IOException{
		System.out.println("\n\nKcnk4 Converted Reads:");
		CloseableIterator<Gene> iter = features.sortedIterator();

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
		AnnotationCollection<DerivedAnnotation<SAMFragment>> converted = features.convertCoordinates(bam, refSpace, true);	
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
	public void Dkk1OverlapReadCount()
	{
		System.out.println("\n\nDkk1 Mapped Reads:");
		CloseableIterator<Gene> iter = features.sortedIterator();

		Annotation a = null;
		while(iter.hasNext()) 
		{
			a = iter.next();
			if(a.getName().equals("NM_010051")) //Dkk1
				break;
		}
		iter.close();
		
		int count =0;
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
	public void Dkk1ConvertCoordinatesFullyContained() throws IOException{
		System.out.println("\n\nDkk1 Converted Reads (Fully Contained):");
		CloseableIterator<Gene> iter = features.sortedIterator();

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
		AnnotationCollection<DerivedAnnotation<SAMFragment>> converted = features.convertCoordinates(bam, refSpace, true);		
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
	public void Dkk1ConvertCoordinates() throws IOException{
		System.out.println("\n\nDkk1 Converted Reads:");
		CloseableIterator<Gene> iter = features.sortedIterator();

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
		AnnotationCollection<DerivedAnnotation<SAMFragment>> converted = features.convertCoordinates(bam, refSpace, false);	
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
	
	
	@Test  //Pass, paired ends take a long time
	//Verifies that the correct number of reads are returned for Ccdc87, a single exon gene on the positive strand, using paired end reads
	public void Ccdc87SPairedEndOverlapReadCount() throws IOException{
		
		BAMPairedFragmentCollection bamPair = new BAMPairedFragmentCollection(new File("/storage/shared/CoreTestData/chr19.clean.sorted.bam"));
		System.out.println("\n\nCcdc87 Mapped Paired Reads:");
		CloseableIterator<Gene> iter = features.sortedIterator();

		Annotation a = null;
		while(iter.hasNext()) 
		{
			a = iter.next();
			if(a.getName().equals("NM_207268"))
				break;
		}
		iter.close();
		
		int count =0;
		bamPair.addFilter(new MaximumLengthFilter<PairedMappedFragment<SAMFragment>>(10000));
		CloseableIterator<PairedMappedFragment<SAMFragment>> f_iter = bamPair.sortedIterator(a, true);
		while(f_iter.hasNext())
		{
			PairedMappedFragment<SAMFragment> f = f_iter.next();
			System.out.println(f.toString());
			count++;
		}
		
		f_iter.close();
		assertEquals("20 unconverted reads should overlap Ccdc87.",20,count);
		
	}

	
	@Test 
	//Verifies that the correct number of reads are returned for Ccdc87, a single exon gene on the positive strand, using paired end reads
	public void Ccdc87SPairedEndConvertCoordinates() throws IOException{
		
		BAMPairedFragmentCollection bamPair = new BAMPairedFragmentCollection(new File("/storage/shared/CoreTestData/chr19.clean.sorted.bam"));
		
		System.out.println("\n\nCcdc87 Converted Paired Reads:");
		CloseableIterator<Gene> iter = features.sortedIterator();

		Annotation a = null;
		while(iter.hasNext()) 
		{
			a = iter.next();
			if(a.getName().equals("NM_207268"))
				break;
		}
		iter.close();
		
		int count =0;
		Annotation b = new SingleInterval(a.getName(), 0, a.size()-1);
		b.setOrientation(Strand.BOTH);
		AnnotationCollection<DerivedAnnotation<PairedMappedFragment<SAMFragment>>> converted = features.convertCoordinates(bamPair, refSpace, true);	
		CloseableIterator<DerivedAnnotation<PairedMappedFragment<SAMFragment>>> c_iter = converted.sortedIterator(b,true);
		
		while(c_iter.hasNext())
		{
			DerivedAnnotation<PairedMappedFragment<SAMFragment>> c = c_iter.next();
			System.out.println(c.toString());
			count++;
		}
		
		assertEquals("20 unconverted reads should overlap Ccdc87.",20,count);
		
	}
	
	
	//Utility Methods
	public ArrayList<String> IterToNameList(CloseableIterator<? extends Annotation> c_iter)
	{
		ArrayList<String> c_list = new ArrayList<String>();

		while(c_iter.hasNext())
			c_list.add(c_iter.next().getName());
		c_iter.close();
		
		return c_list;
	}
	
	public void CompareOutputs(ArrayList<String> a_list, ArrayList<String> b_list)
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
