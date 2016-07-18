package guttmanlab.core.annotationcollection;

import static org.junit.Assert.assertEquals;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.coordinatespace.CoordinateSpace;

import java.util.Map;
import java.util.TreeMap;

import org.junit.Test;

public class TestFeatureCollection {
	
	@Test
	public void featureCollectionMerge2to1() {
		Map<String,Integer> mapping = new TreeMap<String,Integer>();
		mapping.put("a1",1);
		mapping.put("a2",1000);
		CoordinateSpace fcspace = new CoordinateSpace(mapping);
		FeatureCollection<BlockedAnnotation> fc = new FeatureCollection<BlockedAnnotation>(fcspace);
		
		BlockedAnnotation a1 = new BlockedAnnotation();
		SingleInterval block1 = new SingleInterval("a1",100,300);
		a1.addBlocks(block1);
		
		BlockedAnnotation a2 = new BlockedAnnotation();
		SingleInterval block2 = new SingleInterval("a1",200,500);		
		a2.addBlocks(block2);
		
		fc.add(a1);
		fc.add(a2);
		
		FeatureCollection<BlockedAnnotation> fc_merged = fc.merge();
		assertEquals("merged fc contains one annotation",1,fc_merged.getNumAnnotations());
	}
	
	@Test
	public void featureCollectionMerge3to1() {
		Map<String,Integer> mapping = new TreeMap<String,Integer>();
		mapping.put("a1",1);
		mapping.put("a2",1000);
		CoordinateSpace fcspace = new CoordinateSpace(mapping);
		FeatureCollection<BlockedAnnotation> fc = new FeatureCollection<BlockedAnnotation>(fcspace);
		
		BlockedAnnotation a1 = new BlockedAnnotation();
		SingleInterval block1 = new SingleInterval("a1",100,300);
		a1.addBlocks(block1);
		
		BlockedAnnotation a2 = new BlockedAnnotation();
		SingleInterval block2 = new SingleInterval("a1",200,500);		
		a2.addBlocks(block2);
		
		BlockedAnnotation a3 = new BlockedAnnotation();
		SingleInterval block3 = new SingleInterval("a1",350,600);		
		a3.addBlocks(block3);
		
		fc.add(a1);
		fc.add(a2);
		fc.add(a3);
		
		FeatureCollection<BlockedAnnotation> fc_merged = fc.merge();
		assertEquals("merged fc contains one annotation",1,fc_merged.getNumAnnotations());
	}

	@Test
	public void featureCollectionMerge4to2() {
		Map<String,Integer> mapping = new TreeMap<String,Integer>();
		mapping.put("a1",1);
		mapping.put("a2",1000);
		CoordinateSpace fcspace = new CoordinateSpace(mapping);
		FeatureCollection<BlockedAnnotation> fc = new FeatureCollection<BlockedAnnotation>(fcspace);
		
		BlockedAnnotation a1 = new BlockedAnnotation();
		SingleInterval block1 = new SingleInterval("a1",100,300);
		a1.addBlocks(block1);
		
		BlockedAnnotation a2 = new BlockedAnnotation();
		SingleInterval block2 = new SingleInterval("a1",200,500);		
		a2.addBlocks(block2);
		
		BlockedAnnotation a3 = new BlockedAnnotation();
		SingleInterval block3 = new SingleInterval("a1",600,700);		
		a3.addBlocks(block3);
		
		BlockedAnnotation a4 = new BlockedAnnotation();
		SingleInterval block4 = new SingleInterval("a1",650,750);		
		a4.addBlocks(block4);
		
		fc.add(a1);
		fc.add(a2);
		fc.add(a3);
		fc.add(a4);
		
		FeatureCollection<BlockedAnnotation> fc_merged = fc.merge();
		assertEquals("merged fc contains two annotations",2,fc_merged.getNumAnnotations());
	}

	
}
