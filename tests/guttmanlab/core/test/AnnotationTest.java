package guttmanlab.core.test;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.PairedMappedFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotation.predicate.PairedFilterWrapper;
import net.sf.samtools.util.CloseableIterator;

import org.junit.Before;
import org.junit.Test;

public class AnnotationTest {
		private Annotation pos;
		private Annotation pos2;
		private Annotation neg;
		private Annotation both;
		
	@Before
	public void setUp() throws Exception {
		 pos = new BlockedAnnotation();
		 pos2 = new BlockedAnnotation();
		 neg = new BlockedAnnotation();
		 both = new BlockedAnnotation();
		
		pos.setOrientation(Strand.POSITIVE);
		pos2.setOrientation(Strand.POSITIVE);
		neg.setOrientation(Strand.NEGATIVE);
		both.setOrientation(Strand.BOTH);
	}

	@Test
	public void EqualsShouldCompareEnumValues() {
		Strand strand1 = pos.getOrientation();
		Strand strand2 = pos2.getOrientation();
		Strand strand3 = neg.getOrientation();
		
		assertEquals("two positive strands should be equal",strand1.equals(strand2),true);
		assertEquals("two different orientations should be not equal",strand2.equals(strand3),false);
	}
	
	@Test 
	public void trimAnnotationTest()
	{

		BlockedAnnotation blocked = new BlockedAnnotation();
		//DerivedAnnotation converted = new DerivedAnnotation(both, both);
		//PairedMappedFragment paired = new PairedMappedFragment<>(pair));
		SingleInterval block1 = new SingleInterval("a1",100,300);
		SingleInterval block2 = new SingleInterval("a1",350,500);
		SingleInterval block3 = new SingleInterval("a1",600,700);
		
		blocked.addBlocks(block1);
		blocked.addBlocks(block2);
		blocked.addBlocks(block3);
		
		blocked.setOrientation(Strand.POSITIVE);
		Annotation blocked1 = blocked.trim(0,800);    //no trimming
		Annotation blocked2 = blocked.trim(150,650);  //trim first and last exon
		Annotation blocked3 = blocked.trim(330,550);  //between exons
		
		System.out.println(blocked.toString());
		System.out.println(blocked1.toString());
		System.out.println(blocked2.toString());
		System.out.println(blocked3.toString()+"\n");
		
		blocked.setOrientation(Strand.NEGATIVE);
		blocked1 = blocked.trim(0,800);    //no trimming
	    blocked2 = blocked.trim(150,650);  //trim first and last exon
		blocked3 = blocked.trim(330,550);  //between exons
		
		System.out.println(blocked.toString());
		System.out.println(blocked1.toString());
		System.out.println(blocked2.toString());
		System.out.println(blocked3.toString());
		blocked.setOrientation(Strand.BOTH);
		
		/*
		PairedMappedFragment paired1 = paired.trim(0,800);
		PairedMappedFragment paired2 = paired.trim(0,800);
		PairedMappedFragment paired3 = paired.trim(0,800);
		
		DerivedAnnotation converted1 = converted.trim(0,800);
		DerivedAnnotation converted2 = converted.trim(0,800);
		DerivedAnnotation converted3 = converted.trim(0,800);
		*/
		
		
	}
	
	@Test 
	public void testHashCodeEquals()
	{
		BlockedAnnotation blocked = new BlockedAnnotation();
		BlockedAnnotation blocked2 = new BlockedAnnotation();
		BlockedAnnotation blocked3 = new BlockedAnnotation();

		//DerivedAnnotation converted = new DerivedAnnotation(both, both);
		//PairedMappedFragment paired = new PairedMappedFragment<>(pair));
		SingleInterval block1 = new SingleInterval("a1",100,300);
		SingleInterval block2 = new SingleInterval("a1",350,500);
		SingleInterval block3 = new SingleInterval("a1",600,700);
	
		blocked.addBlocks(block1);
		blocked.addBlocks(block2);
		blocked.addBlocks(block3);
		
		blocked2.addBlocks(block1);
		blocked2.addBlocks(block2);
		blocked2.addBlocks(block3);
		
		blocked3.addBlocks(block3);
		
		System.out.println(blocked.hashCode());
		System.out.println(blocked2.hashCode());
		System.out.println(blocked3.hashCode());
		System.out.println(block3.hashCode());
		
		assertEquals("blocked1 and 3 not equal.",false,blocked.equals(blocked3));
		assertEquals("blocked1 and 2 equal.",true,blocked.equals(blocked2));
		assertEquals("blocked3 and single interval 3 equal.",true,block3.equals(blocked3));
		assertEquals("blocked3 and single interval 3 have same hash code.",true,blocked3.hashCode()==block3.hashCode());
		assertEquals("blocked2 and 3 have different hash codes.",false,blocked3.hashCode()==blocked2.hashCode());
	}
	
	@Test 
	public void testAnnotationSize()
	{
		//single exon positive strand
		BlockedAnnotation a1 = new BlockedAnnotation();
		SingleInterval b1 = new SingleInterval("chr1",100,500,Strand.POSITIVE);
		a1.addBlocks(b1);
		assertEquals("single exon + strand size is 400bp",400,a1.size());
		
		//single exon negative strand
		BlockedAnnotation a2 = new BlockedAnnotation();
		SingleInterval b2 = new SingleInterval("chr1",100,500,Strand.NEGATIVE);
		a2.addBlocks(b2);
		assertEquals("single exon + strand size is 400bp",400,a2.size());
		
		//multi exon positive strand
		BlockedAnnotation a3 = new BlockedAnnotation();
		SingleInterval b3_1 = new SingleInterval("chr1",100,500,Strand.POSITIVE);
		SingleInterval b3_2 = new SingleInterval("chr1",900,1500,Strand.POSITIVE);
		SingleInterval b3_3 = new SingleInterval("chr1",2000,2500,Strand.POSITIVE);
		a3.addBlocks(b3_1);
		a3.addBlocks(b3_2);
		a3.addBlocks(b3_3);
		assertEquals("multi exon + strand size is 1500bp",1500,a3.size());
		//multi-exon negative strand
	}
	@Test
	public void testOverlapsOnUnknownStrand(){
		//get an interval with unknown strand
		//get reads that should overlap
		
	}
	
	@Test
	public void testRefSeq() throws NumberFormatException, IOException{
			
		
		//sample refseq genes
		BEDFileIO io =  new BEDFileIO("/Users/cburghard/Downloads/sizes");
		String FeatureFile = "/Users/cburghard/Downloads/RefSeq.bed";
		CloseableIterator<? extends Annotation> features = io.loadFromFile(FeatureFile).sortedIterator();
		PrintWriter writer = new PrintWriter("/Users/cburghard/Downloads/geneLengthsNew.txt");
		
		while(features.hasNext())
		{
			Annotation feature = features.next();
			writer.println(feature.getName()+"\t"+feature.size());
		}
		
		features.close();
		writer.close();
		
		
		//read in bed file as normal IO
		PrintWriter writer2 = new PrintWriter("/Users/cburghard/Downloads/geneLengthsBed.txt");
		File refseq = new File("/Users/cburghard/Downloads/RefSeq.bed");
		BufferedReader br = new BufferedReader(new FileReader(refseq));
	    String line;
	    while ((line = br.readLine()) != null) {
		    int length = 0;
	    	String name = line.split("\t")[3];
	    	String[] exons = line.split("\t")[10].split(",");
	    	for(String num : exons)
	    	{
	    		length+=Integer.valueOf(num);
	    	}
	    	writer2.println(name+"\t"+length);
	    }
	    br.close();
		writer2.close();
	}
	

}
