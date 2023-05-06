package guttmanlab.core.annotation.io;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.FeatureCollection;
import guttmanlab.core.coordinatespace.CoordinateSpace;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.splicing.speckle.GTFToJunctions;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * This implementation will parse an entire file and read it into memory
 * @author mguttman
 *
 */
public class BEDFileIO implements AnnotationFileIO<Gene> {

	private CoordinateSpace referenceSpace;
	
	/**
	 * @param referenceSizes The reference coordinate information containing names and sizes
	 */
	public BEDFileIO(String referenceSizes){
		this.referenceSpace=new CoordinateSpace(referenceSizes);
	}
	
	public BEDFileIO(CoordinateSpace referenceSizes){
		this.referenceSpace=referenceSizes;
	}
	
	/**
	 * Write collection of annotations to a file
	 * @param regions The annotations to write
	 * @param outputBed Output bed file
	 * @throws IOException
	 */
	public static void writeToFile(AnnotationCollection<? extends Annotation> regions, String outputBed) throws IOException {
		FileWriter w = new FileWriter(outputBed);
		Iterator<? extends Annotation> iter = regions.sortedIterator();
		while(iter.hasNext()) {
			w.write(iter.next().toBED() + "\n");
		}
		w.close();
	}
	
	/**
	 * Static method to get the annotation collection represented in a bed file
	 * @param geneFile Bed file name 
	 * @param coordinateSpace The reference coordinate information containing names and sizes
	 * @return The collection of genes described in the bed file
	 * @throws IOException
	 */
	public static AnnotationCollection<Gene> loadFromFile(File geneFile, CoordinateSpace coordinateSpace) throws IOException {
		BEDFileIO bfio = new BEDFileIO(coordinateSpace);
		return bfio.loadFromFile(geneFile.getAbsolutePath());
	}
	
	/**
	 * Static method to get the annotation collection represented in a bed file
	 * @param geneFile Bed file name 
	 * @param chrSizeFile Table of chromosome names and sizes
	 * @return The collection of genes described in the bed file
	 * @throws IOException
	 */
	public static AnnotationCollection<Gene> loadFromFile(String fileName, String chrSizeFile) throws IOException {
		return loadFromFile(fileName, new CoordinateSpace(chrSizeFile));
	}
	
	public static AnnotationCollection<Gene> loadFromFile(String fileName, CoordinateSpace space) throws IOException {
		BEDFileIO bfio = new BEDFileIO(space);
		return bfio.loadFromFile(fileName);
	}
	
	/**
	 * Static method to get the annotation collection represented in a bed file, organized by reference name
	 * @param fileName Bed file name 
	 * @param referenceSizes The reference coordinate information containing names and sizes
	 * @return Map of reference name to the collection of genes on that reference described in the bed file
	 * @throws IOException
	 */
	public static Map<String, FeatureCollection<Gene>> loadFromFileByReferenceName(String fileName, String referenceSizes) throws IOException {
		CoordinateSpace refSpace = new CoordinateSpace(referenceSizes);
		Map<String, FeatureCollection<Gene>> rtrn = new TreeMap<String, FeatureCollection<Gene>>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			Gene annotation=parse(nextLine);
			String reference = annotation.getReferenceName();
			if(!rtrn.containsKey(reference)) {
				rtrn.put(reference, new FeatureCollection<Gene>(refSpace));
			}
			rtrn.get(reference).addAnnotation(annotation);
		}
		reader.close();
		return rtrn;
	}
	
	public static AnnotationCollection<Gene> loadFromFile(String fileName) throws IOException {
		FeatureCollection<Gene> collection=new FeatureCollection<Gene>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			Gene annotation=parse(nextLine);
			collection.addAnnotation(annotation);
		}
		reader.close();
		return collection;
	}
	
	public static Collection<SingleInterval> loadCustom(String fileName) throws IOException {
		return loadCustom(fileName, 1.0, 1.0, 0.0);
	}
		
	
	public static Collection<SingleInterval> loadCustom(String fileName, double windowP, double localP, double enrichmentCutoff) throws IOException {
		Collection<SingleInterval> collection=new TreeSet<SingleInterval>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
		String nextLine=reader.readLine();
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String chr=tokens[0].split(":")[0];
			String strand=tokens[0].substring(tokens[0].length()-1, tokens[0].length());
			String noStrand=tokens[0].substring(0, tokens[0].length()-1);
			int start=new Integer(noStrand.split(":")[1].split("-")[0]);
			int end=new Integer(noStrand.split(":")[1].split("-")[1]);
			String name=tokens[1];
			SingleInterval annotation=new SingleInterval(chr, start, end, Strand.fromString(strand), name);
			
			double pVal=new Double(tokens[8]);
			double localPVal=new Double(tokens[9]);
			double enrichment=new Double(tokens[7]);
			
			if(pVal<=windowP && localPVal<=localP && enrichment>enrichmentCutoff){
				collection.add(annotation);
			}
		}
		reader.close();
		return collection;
	}
	
	public static List<String> loadLines(String fileName) throws IOException {
		List<String> rtrn=new ArrayList<String>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			rtrn.add(nextLine);
		}
		reader.close();
		return rtrn;
	}
	
	
	public static List<String> loadLines(String fileName, int skip) throws IOException {
		List<String> rtrn=new ArrayList<String>();
		
		int counter=0;
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			if(counter>=skip) {rtrn.add(nextLine);}
			counter++;
		}
		reader.close();
		return rtrn;
	}
	
	
	public static List<SingleInterval> loadSingleIntervalsInOrder(String fileName) throws IOException {
		List<SingleInterval> collection=new ArrayList<SingleInterval>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			SingleInterval annotation=parseSingleInterval(nextLine);
			collection.add(annotation);
		}
		reader.close();
		return collection;
	}
	
	public static Collection<SingleInterval> loadSingleIntervalFromFile(String fileName) throws IOException {
		Collection<SingleInterval> collection=new TreeSet<SingleInterval>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			SingleInterval annotation=parseSingleInterval(nextLine);
			collection.add(annotation);
		}
		reader.close();
		return collection;
	}
	
	public static Collection<SingleInterval> loadSingleIntervalFromFileFromGTF(String fileName, boolean onlyproteins) throws IOException {
		Collection<SingleInterval> collection=new TreeSet<SingleInterval>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
		String nextLine;
		
		while ((nextLine = reader.readLine()) != null) {
			if(!nextLine.startsWith("#")) {
				String[] tokens=nextLine.split("\t");
				if(tokens[2].equals("gene")) {
					SingleInterval exon=new SingleInterval(tokens[0], Integer.parseInt(tokens[3])-1, Integer.parseInt(tokens[4]));
					exon.setOrientation(Strand.fromString(tokens[6]));
					String geneName=getGeneName(tokens[8]);
					exon.setName(geneName);
					if(onlyproteins) {
						String geneType=getTag(tokens[8], "gene_type");
						if(geneType.equalsIgnoreCase("protein_coding")) {collection.add(exon);}
					}
					else {
						collection.add(exon);
					}
				}
			}
		}
		
		
		
	
		reader.close();
		return collection;
	}

	
	private static String getGeneName(String string) {
		return getTag(string, "gene_name");
	}
	
	private static String getTag(String string, String tag) {
		String[] tokens=string.split(";");
		
		for(int i=0; i<tokens.length; i++) {
			String token=tokens[i].trim();
			String key=token.split(" ")[0];
			String val=token.split(" ")[1];
			val=val.replaceAll("\"","");
			if(key.equals(tag)) {return val;}	
		}
		return null;
	}
	

	public static Collection<SingleInterval> loadSingleIntervalFromFileStrand(String fileName) throws IOException {
		Collection<SingleInterval> collection=new TreeSet<SingleInterval>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			SingleInterval annotation=parseSingleIntervalStrand(nextLine);
			collection.add(annotation);
		}
		reader.close();
		return collection;
	}
	
	/**
	 * For each line, parse it into an Annotation with all information
	 * @param rawData
	 * @return An Annotation with all the information from the BED line
	 */
	public static Gene parse(String rawData) {
		String[] tokens = rawData.split("\t");
		String chr=(tokens[0]);
		int start=new Integer(tokens[1]);
		int end=new Integer(tokens[2]);
		Strand orientation=Strand.UNKNOWN;
		if(tokens.length > 3) {
			String name=tokens[3];
			if(tokens.length > 4) {
				@SuppressWarnings("unused")
				double score = new Double(tokens[4]); //TODO Currently unused
				if(tokens.length > 5){
					orientation= Annotation.Strand.fromString(tokens[5]);
					if(tokens.length > 11) { //Needs to have ALL block info to be parsed
						int cdsStart=Integer.parseInt(tokens[6]);
						int cdsEnd=Integer.parseInt(tokens[7]);
						String[] blockSizes=tokens[10].split(",");
						String[] blockStarts=tokens[11].split(",");
						Pair<List<Integer>> exonStartEnd=getBlockStartsAndEnds(blockStarts, blockSizes, new Integer(tokens[9]), start);
						Collection<Annotation> exons=new ArrayList<Annotation>();
						for(int i = 0; i < blockSizes.length; i++ ) {
							Annotation exon=new SingleInterval(chr, exonStartEnd.getValue1().get(i), exonStartEnd.getValue2().get(i), orientation, name);
							exons.add(exon);
						}
						Gene blockedAnnotation = new Gene(exons, cdsStart, cdsEnd, name);
						return blockedAnnotation;
					}
					else{ //size=6: Has orientation
						Gene g=new Gene(new SingleInterval(chr, start, end, orientation, name));
						return g;
					}
				}
				else{ //size=5: Has score
					Gene g=new Gene(new SingleInterval(chr, start, end, orientation, name));
					return g;
				}
			}
			else{ //size=4: Has name
				Gene g = new Gene(new SingleInterval(chr, start, end, orientation, name));
				return g;
			}
		}
		else{ // size=3: Has positions only
			String name=chr+":"+start+"-"+end;
			Gene g = new Gene(new SingleInterval(chr, start, end, orientation, name));
			return g;
		}
		
	}
	
	private static SingleInterval parseSingleInterval(String rawData) {
		String[] tokens = rawData.split("\t");
		String chr=(tokens[0]);
		int start=new Integer(tokens[1]);
		int end=new Integer(tokens[2]);
		Strand orientation=Strand.UNKNOWN;
		SingleInterval rtrn=new SingleInterval(chr, start, end);
		if(tokens.length > 3) {
			String name=tokens[3];
			rtrn.setName(name);
		}
		if(tokens.length>5) {
			Strand o=Strand.fromString(tokens[5]);
			rtrn.setOrientation(o);
		}
		
		return rtrn;
	}
	
	private static SingleInterval parseSingleIntervalStrand(String rawData) {
		String[] tokens = rawData.split("\t");
		String chr=(tokens[0]);
		int start=new Integer(tokens[1]);
		int end=new Integer(tokens[2]);
		Strand orientation=Strand.UNKNOWN;
		SingleInterval rtrn=new SingleInterval(chr, start, end);
		if(tokens.length > 3) {
			orientation=Strand.fromString(tokens[3]);
			rtrn.setOrientation(orientation);
		}
		
		return rtrn;
	}
	
	private static Pair<List<Integer>> getBlockStartsAndEnds(String[] blockStarts, String[] blockSizes, int size, int start){
		List<Integer> starts=new ArrayList<Integer> ();
		List<Integer>  end=new ArrayList<Integer> ();
		for(int i=0; i<size; i++){
			starts.add(start+new Integer(blockStarts[i].replaceAll("\"", "").trim()));
			end.add((Integer)starts.get(i)+new Integer(blockSizes[i].replaceAll("\"", "").trim()));
		}
		
		Pair<List<Integer>> rtrn=new Pair<List<Integer>>(starts, end);
		return rtrn;
	}

	public static Collection<Gene> loadRegionsFromFile(String file) throws IOException {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		AnnotationCollection<Gene> genes=loadFromFile(file);
		Iterator<Gene> iter=genes.sortedIterator();
		while(iter.hasNext()){
			Gene gene=iter.next();
			rtrn.add(gene);
		}
		return rtrn;
	}
	
	public static Collection<Gene> loadRegionsFromFile(String file, boolean reverse) throws IOException {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		AnnotationCollection<Gene> genes=loadFromFile(file);
		Iterator<Gene> iter=genes.sortedIterator();
		while(iter.hasNext()){
			Gene gene=iter.next();
			if(reverse) {
				Strand s=Strand.antisense(gene.getOrientation());
				gene.setOrientation(s);
			}
			rtrn.add(gene);
		}
		return rtrn;
	}
	
	public static Collection<SingleInterval> loadSingleIntervalsFromFile(String file) throws IOException {
		Collection<SingleInterval> rtrn=loadSingleIntervalFromFile(file);
		return rtrn;
	}

	public static Map<String, Collection<Gene>> loadRegionsFromFileByChr(String bEDFile) throws IOException {
		Collection<Gene> genes=loadRegionsFromFile(bEDFile);
		Map<String, Collection<Gene>> rtrn=new TreeMap<String, Collection<Gene>>();
		
		for(Gene gene: genes){
			String chr=gene.getReferenceName();
			Collection<Gene> list=new TreeSet<Gene>();
			if(rtrn.containsKey(chr)){
				list=rtrn.get(chr);
			}
			list.add(gene);
			rtrn.put(chr, list);
		}
		
		return rtrn;
	}
	
	public static Map<String, Collection<Gene>> loadRegionsFromFileByChr(String bEDFile, boolean reverseStrand) throws IOException {
		Collection<Gene> genes=loadRegionsFromFile(bEDFile, reverseStrand);
		Map<String, Collection<Gene>> rtrn=new TreeMap<String, Collection<Gene>>();
		
		for(Gene gene: genes){
			String chr=gene.getReferenceName();
			Collection<Gene> list=new TreeSet<Gene>();
			if(rtrn.containsKey(chr)){
				list=rtrn.get(chr);
			}
			list.add(gene);
			rtrn.put(chr, list);
		}
		
		return rtrn;
	}
	
	
	public static Map<String, Collection<Gene>> loadIntronsByChr(String bEDFile) throws IOException {
		Collection<Gene> genes=loadRegionsFromFile(bEDFile);
		Map<String, Collection<Gene>> rtrn=new TreeMap<String, Collection<Gene>>();
		
		for(Gene gene: genes){
			String chr=gene.getReferenceName();
			Collection<Gene> list=new TreeSet<Gene>();
			if(rtrn.containsKey(chr)){
				list=rtrn.get(chr);
			}
			for(Annotation intron:gene.getIntrons()){
				list.add(new Gene(intron));
			}
			
			rtrn.put(chr, list);
		}
		
		return rtrn;
	}

	public static Collection<String> loadNames(String fileName) throws IOException {
		Collection<String> rtrn=new TreeSet<String>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String name=nextLine.trim();
			rtrn.add(name.toUpperCase());
		}
		reader.close();
		return rtrn;
	}

	public static Map<String, String> loadIDToName(String fileName) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String id=nextLine.split("\t")[0];
			String name=nextLine.split("\t")[1];
			rtrn.put(id, name);
		}
		reader.close();
		return rtrn;
	}

	
	public static Map<String, IntervalTree<Annotation>> loadIntergenicRegions(String file) throws IOException{
		Map<String, IntervalTree<Annotation>> rtrn=new TreeMap<String, IntervalTree<Annotation>>();
		Map<String, IntervalTree<Gene>> trees=BEDFileIO.loadTree(file);
		for(String chr: trees.keySet()) {
			IntervalTree<Annotation> temp=new IntervalTree<Annotation>();
			IntervalTree<Gene> tree=trees.get(chr);
			Iterator<Gene> iter=tree.valueIterator();
			int start=0;
			int end=0;
			while(iter.hasNext()) {
				Gene g=iter.next();
				end=g.getReferenceStartPosition();
				//Annotation r=new SingleInterval(chr, start, end);
				if(end>start) {
					Annotation pos=new SingleInterval(chr, start, end);
					pos.setOrientation(Strand.POSITIVE);
					Annotation neg=new SingleInterval(chr, start, end);
					neg.setOrientation(Strand.NEGATIVE);
					temp.put(start, end, pos);
					temp.put(start, end, neg);
				}
				start=Math.max(start, g.getReferenceEndPosition());
			}
			rtrn.put(chr, temp);
		}
		return rtrn;
	}
	
	public static Map<String, IntervalTree<Gene>> loadTree(String fileName) throws IOException {
		Map<String, Collection<Gene>> genes=loadRegionsFromFileByChr(fileName);
		Map<String, IntervalTree<Gene>> rtrn=new TreeMap<String, IntervalTree<Gene>>();
		
		for(String chr: genes.keySet()){
			Collection<Gene> geneList=genes.get(chr);
			IntervalTree<Gene> tree=new IntervalTree<Gene>();
			for(Gene gene: geneList){tree.put(gene.getReferenceStartPosition(), gene.getReferenceEndPosition(), gene);}
			rtrn.put(chr, tree);
		}
		return rtrn;
	}
	
	public static Map<String, IntervalTree<Double>> loadTreePlusExpression(String fileName) throws IOException {
		Map<String, IntervalTree<Double>> rtrn=new TreeMap<String, IntervalTree<Double>>();
		
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String chr=tokens[0];
			int start=Integer.parseInt(tokens[1]);
			int end=Integer.parseInt(tokens[2]);
			
			Double val=Double.parseDouble(tokens[5]);
			if(!rtrn.containsKey(chr)) {rtrn.put(chr, new IntervalTree<Double>());}
			
			IntervalTree<Double> tree=rtrn.get(chr);
			tree.put(start, end, val);
			
		}
		reader.close();
		return rtrn;
	}
	
	
	public static Map<String, IntervalTree<Annotation>> loadTreeAnnotation(String fileName) throws IOException {
		Map<String, Collection<Gene>> genes=loadRegionsFromFileByChr(fileName);
		Map<String, IntervalTree<Annotation>> rtrn=new TreeMap<String, IntervalTree<Annotation>>();
		
		for(String chr: genes.keySet()){
			Collection<Gene> geneList=genes.get(chr);
			IntervalTree<Annotation> tree=new IntervalTree<Annotation>();
			for(Gene gene: geneList){tree.put(gene.getReferenceStartPosition(), gene.getReferenceEndPosition(), gene);}
			rtrn.put(chr, tree);
		}
		return rtrn;
	}
	
	
	public static Map<String, IntervalTree<String>> loadGeneNamesFromRefFlat(String fileName) throws IOException {
		List<String> lines=loadLines(fileName);
		
		Map<String, IntervalTree<String>> rtrn=new TreeMap<String, IntervalTree<String>>();
		
		for(String line: lines){
			String chr=line.split("\t")[2];
			int start=Integer.parseInt(line.split("\t")[4]);
			int end=Integer.parseInt(line.split("\t")[5]);
			
			
			if(!rtrn.containsKey(chr)) {
				IntervalTree<String> tree=new IntervalTree<String>();
				rtrn.put(chr, tree);
			}
			IntervalTree<String> tree=rtrn.get(chr);
			tree.put(start, end, line.split("\t")[0]);
			
			rtrn.put(chr, tree);
		}
		return rtrn;
	}
	
	
	public static Map<String, IntervalTree<Gene>> loadTreePlusIntrons(String fileName) throws IOException {
		Map<String, Collection<Gene>> genes=loadRegionsFromFileByChr(fileName);
		Map<String, IntervalTree<Gene>> rtrn=new TreeMap<String, IntervalTree<Gene>>();
		
		for(String chr: genes.keySet()){
			Collection<Gene> geneList=genes.get(chr);
			IntervalTree<Gene> tree=new IntervalTree<Gene>();
			for(Gene gene: geneList){
				tree.put(gene.getReferenceStartPosition(), gene.getReferenceEndPosition(), gene);
			
				Collection<Annotation> introns=gene.getIntrons();
				if(!introns.isEmpty()) {
					Gene newIntron=new Gene(introns, gene.getName()+"_intron");
					tree.put(newIntron.getReferenceStartPosition(), newIntron.getReferenceEndPosition(), newIntron);
				}
				
				
			}
			rtrn.put(chr, tree);
		}
		return rtrn;
	}
	
	public static Map<String, IntervalTree<Gene>> loadTreeGenomeBin(String fileName) throws IOException {
		Map<String, Collection<Gene>> genes=loadRegionsFromFileByChr(fileName);
		Map<String, IntervalTree<Gene>> rtrn=new TreeMap<String, IntervalTree<Gene>>();
		
		for(String chr: genes.keySet()){
			Collection<Gene> geneList=genes.get(chr);
			IntervalTree<Gene> tree=new IntervalTree<Gene>();
			for(Gene g: geneList){
				Gene gene=new Gene(g.getGenomicRegion());
				tree.put(gene.getReferenceStartPosition(), gene.getReferenceEndPosition(), gene);
			}
			rtrn.put(chr, tree);
		}
		return rtrn;
	}
	
	public static Map<String, IntervalTree<Annotation>> loadTreePlusIntronsAnnotation(String fileName) throws IOException {
		Map<String, Collection<Gene>> genes=loadRegionsFromFileByChr(fileName);
		Map<String, IntervalTree<Annotation>> rtrn=new TreeMap<String, IntervalTree<Annotation>>();
		
		for(String chr: genes.keySet()){
			Collection<Gene> geneList=genes.get(chr);
			IntervalTree<Annotation> tree=new IntervalTree<Annotation>();
			for(Gene gene: geneList){
				Gene newGene=new Gene(gene.getBlockSet(), gene.getName()+"_exon");
				tree.put(gene.getReferenceStartPosition(), gene.getReferenceEndPosition(), newGene);
			
				
				Collection<Annotation> introns=gene.getIntrons();
				if(!introns.isEmpty()) {
					Gene newIntron=new Gene(introns, gene.getName()+"_intron");
					tree.put(newIntron.getReferenceStartPosition(), newIntron.getReferenceEndPosition(), newIntron);
				}
				
				
			}
			rtrn.put(chr, tree);
		}
		return rtrn;
	}
	
	
	public static Map<String, IntervalTree<Annotation>> loadTreePlusIntronsAnnotation(String fileName, boolean reverseStrand) throws IOException {
		Map<String, Collection<Gene>> genes=loadRegionsFromFileByChr(fileName, reverseStrand);
		Map<String, IntervalTree<Annotation>> rtrn=new TreeMap<String, IntervalTree<Annotation>>();
		
		for(String chr: genes.keySet()){
			Collection<Gene> geneList=genes.get(chr);
			IntervalTree<Annotation> tree=new IntervalTree<Annotation>();
			for(Gene gene: geneList){
				Gene newGene=new Gene(gene.getBlockSet(), gene.getName()+"_exon");
				//if(reverseStrand) {newGene.setOrientation(Strand.antisense(gene.getOrientation()));}
				
				tree.put(gene.getReferenceStartPosition(), gene.getReferenceEndPosition(), newGene);
			
				
				Collection<Annotation> introns=gene.getIntrons();
				if(!introns.isEmpty()) {
					Gene newIntron=new Gene(introns, gene.getName()+"_intron");
					//if(reverseStrand) {newIntron.setOrientation(Strand.antisense(gene.getOrientation()));}
					tree.put(newIntron.getReferenceStartPosition(), newIntron.getReferenceEndPosition(), newIntron);
				}
				
				//System.err.println(gene.getName()+" "+gene.getOrientation()+" "+newGene.getOrientation());
			}
			rtrn.put(chr, tree);
		}
		return rtrn;
	}
	
	
	public static Map<String, IntervalTree<Gene>> loadIntronTree(String fileName) throws IOException {
		Map<String, Collection<Gene>> genes=loadIntronsByChr(fileName);
		Map<String, IntervalTree<Gene>> rtrn=new TreeMap<String, IntervalTree<Gene>>();
		
		for(String chr: genes.keySet()){
			Collection<Gene> geneList=genes.get(chr);
			IntervalTree<Gene> tree=new IntervalTree<Gene>();
			for(Gene gene: geneList){tree.put(gene.getReferenceStartPosition(), gene.getReferenceEndPosition(), gene);}
			rtrn.put(chr, tree);
		}
		return rtrn;
	}

	
	public static Map<String, IntervalTree<SingleInterval>> loadRepeatTree(String fileName) throws IOException {
		Collection<SingleInterval> genes=loadSingleIntervalFromFile(fileName);
				
		Map<String, IntervalTree<SingleInterval>> rtrn=new TreeMap<String, IntervalTree<SingleInterval>>();
		
		for(SingleInterval region: genes){
			String chr=region.getReferenceName();
			IntervalTree<SingleInterval> tree=new IntervalTree<SingleInterval>();
			if(rtrn.containsKey(chr)){tree=rtrn.get(chr);}
			tree.put(region.getReferenceStartPosition(), region.getReferenceEndPosition(), region);
			rtrn.put(chr, tree);
		}
		return rtrn;
	}

	public static TreeMap<SingleInterval, Double> loadbedgraph(File fileName) throws IOException {
		TreeMap<SingleInterval, Double> collection=new TreeMap<SingleInterval, Double>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens = nextLine.split("\t");
			String chr=(tokens[0]);
			int start=new Integer(tokens[1]);
			int end=new Integer(tokens[2]);
			SingleInterval annotation=new SingleInterval(chr, start, end);
			collection.put(annotation, new Double(tokens[3]));
		}
		reader.close();
		return collection;
	}

	public static List<String> getGeneNames(String genes) throws IOException {
		Collection<Gene> list=BEDFileIO.loadRegionsFromFile(genes);
		List<String> rtrn=new ArrayList<String>();
		
		for(Gene gene: list){
			if(!rtrn.contains(gene.getName())){
				rtrn.add(gene.getName());
			}
		}
		return rtrn;
	}

	public static Map<String, IntervalTree<SingleInterval>> loadSingleIntervalTree(String fileName) throws IOException {
			Collection<SingleInterval> regions=loadSingleIntervalFromFile(fileName);
					
			Map<String, IntervalTree<SingleInterval>> rtrn=new TreeMap<String, IntervalTree<SingleInterval>>();
			
			for(SingleInterval region: regions){
				IntervalTree<SingleInterval> tree=new IntervalTree<SingleInterval>();
				if(rtrn.containsKey(region.getReferenceName())) {tree=rtrn.get(region.getReferenceName());}
				tree.put(region.getReferenceStartPosition(), region.getReferenceEndPosition(), region);
				rtrn.put(region.getReferenceName(), tree);
			
			}
			return rtrn;
		}
	
	
	
	public static Map<String, IntervalTree<SingleInterval>> loadSingleIntervalTree(Collection<SingleInterval> regions) {
		Map<String, IntervalTree<SingleInterval>> rtrn=new TreeMap<String, IntervalTree<SingleInterval>>();
		
		for(SingleInterval region: regions){
			IntervalTree<SingleInterval> tree=new IntervalTree<SingleInterval>();
			if(rtrn.containsKey(region.getReferenceName())) {tree=rtrn.get(region.getReferenceName());}
			tree.put(region.getReferenceStartPosition(), region.getReferenceEndPosition(), region);
			rtrn.put(region.getReferenceName(), tree);
		
		}
		return rtrn;
	}

	
	public static Map<String, IntervalTree<SingleInterval>> loadSingleIntervalTreeFromGTF(String fileName, boolean onlyproteins) throws IOException {
		Collection<SingleInterval> regions=loadSingleIntervalFromFileFromGTF(fileName, onlyproteins);
				
		Map<String, IntervalTree<SingleInterval>> rtrn=new TreeMap<String, IntervalTree<SingleInterval>>();
		
		for(SingleInterval region: regions){
			IntervalTree<SingleInterval> tree=new IntervalTree<SingleInterval>();
			if(rtrn.containsKey(region.getReferenceName())) {tree=rtrn.get(region.getReferenceName());}
			tree.put(region.getReferenceStartPosition(), region.getReferenceEndPosition(), region);
			rtrn.put(region.getReferenceName(), tree);
		
		}
		return rtrn;
	}

	
	public static Map<String, IntervalTree<String>> parseSNPTree(File snpFile, Collection<String> chrs) throws NumberFormatException, IOException {
		Map<String, IntervalTree<String>> rtrn=new TreeMap<String, IntervalTree<String>>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(snpFile)));
		String nextLine;
		int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			if(chrs.contains(tokens[0])) {
			SingleInterval region=new SingleInterval(tokens[0], new Integer(tokens[1]), new Integer(tokens[2]));
			String genotype=tokens[3];
			IntervalTree<String> tree;
			if(rtrn.containsKey(region.getReferenceName())){tree=rtrn.get(region.getReferenceName());}
			else{tree=new IntervalTree<String>();}
			tree.put(region.getReferenceStartPosition(), region.getReferenceEndPosition(), genotype);
			rtrn.put(region.getReferenceName(), tree);
			}
			counter++;
			if(counter%1000000==0){System.err.println(counter);}
		}
		reader.close();
		return rtrn;
}
	
	public static Map<String, IntervalTree<String>> parseSNPTree(File snpFile) throws NumberFormatException, IOException {
			Map<String, IntervalTree<String>> rtrn=new TreeMap<String, IntervalTree<String>>();
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(snpFile)));
			String nextLine;
			int counter=0;
			while ((nextLine = reader.readLine()) != null) {
				String[] tokens=nextLine.split("\t");
				SingleInterval region=new SingleInterval(tokens[0], new Integer(tokens[1]), new Integer(tokens[2]));
				String genotype=tokens[3];
				IntervalTree<String> tree;
				if(rtrn.containsKey(region.getReferenceName())){tree=rtrn.get(region.getReferenceName());}
				else{tree=new IntervalTree<String>();}
				tree.put(region.getReferenceStartPosition(), region.getReferenceEndPosition(), genotype);
				rtrn.put(region.getReferenceName(), tree);
				counter++;
				if(counter%1000000==0){System.err.println(counter);}
			}
			reader.close();
			return rtrn;
	}

	public static void writeBEDGraph(Map<SingleInterval, Double> map, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: map.keySet()) {
			double val=map.get(region);
			writer.write(region.toBedgraph(val)+"\n");
		}
		
		writer.close();
		
	}
	
	
	public static void writeBEDGraph(Map<SingleInterval, Double> map, String save, String chr) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: map.keySet()) {
			if(region.getReferenceName().equals(chr)) {
			double val=map.get(region);
			writer.write(region.toBedgraph(val)+"\n");
			}
		}
		
		writer.close();
		
	}

	public static Set<String> loadLineSet(String readFile) throws IOException {
		Set<String> rtrn=new HashSet<String>();
		
		int counter=0;
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(readFile)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			rtrn.add(nextLine);
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		reader.close();
		return rtrn;
	}

	public static void writeBED(Collection<? extends Annotation> randomRegions, String string) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		for(Annotation r: randomRegions) {
			writer.write(r.toBED()+"\n");
		}
		
		writer.close();
	}
	
	public static void writeShortBED(Collection<SingleInterval> randomRegions, String string) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		for(Annotation r: randomRegions) {
			writer.write(r.toBED()+"\n");
		}
		
		writer.close();
	}

	public static void writeUCSCScore(Map<SingleInterval, Double> scores, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval r: scores.keySet()) {
			double score=scores.get(r);
			writer.write(r.toUCSC()+"\t"+score+"\n");
		}
		
		writer.close();	
	}

	public static void write(Map<String, ? extends Number> counts, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String chr: counts.keySet()) {writer.write(chr+"\t"+counts.get(chr)+"\n");}
		
		writer.close();
	}

	public static Map<String, String> parseRefFlat(String fileName) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		List<String> lines=loadLines(fileName);
		for(String line: lines) {
			String[] tokens=line.split("\t");
			rtrn.put(tokens[1], tokens[0]);
		}
		return rtrn;	
	}

	public static Map<String, Gene> loadGenesByName(String string) throws IOException {
		Collection<Gene> genes=loadRegionsFromFile(string);
		Map<String, Gene> rtrn=new TreeMap<String, Gene>();
		
		for(Gene g: genes) {
			rtrn.put(g.getName(), g);
		}
		
		return rtrn;
	}
	
	
	public static double writeBEDGraphInteger(Map<SingleInterval, Integer> map, String save, Collection<String> chromosomesToUse) throws IOException {
		Map<String, Double> maxByChr=new TreeMap<String, Double>();
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: map.keySet()) {
			if(chromosomesToUse.contains(region.getReferenceName())){
				double val=map.get(region);
				double max=val;
				if(maxByChr.containsKey(region.getReferenceName())) {
					max=maxByChr.get(region.getReferenceName());
				}
				max=Math.max(max, val);
				maxByChr.put(region.getReferenceName(), max);
				writer.write(region.toBedgraph(val)+"\n");
			}
		}
		
		double rtrn=getMin(maxByChr);
		
		writer.close();
		return rtrn;
	}
	
	
	public static void writeBEDGraphInteger(Map<SingleInterval, Integer> map, String save, String chrToUse) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: map.keySet()) {
			if(chrToUse.equals(region.getReferenceName())){
				double val=map.get(region);
				writer.write(region.toBedgraph(val)+"\n");
			}
		}
		
		
		writer.close();
	}

	public static void writeBEDGraphInteger(Map<SingleInterval, Integer> map, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: map.keySet()) {
			double val=map.get(region);
			writer.write(region.toBedgraph(val)+"\n");
		}
		
		
		writer.close();
	}

	private static double getMin(Map<String, Double> maxByChr) {
		if(maxByChr.isEmpty()) {return 0;}
		double min=maxByChr.values().iterator().next();
		
		for(String chr: maxByChr.keySet()) {
			min=Math.min(min, maxByChr.get(chr));
		}
		
		return min;
	}

	public static Map<String, String> loadNameToChromosome(String fileName) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String chr=tokens[0];
			String name=tokens[3].split("\\.")[0];
			rtrn.put(name, chr);
		}
		reader.close();
		return rtrn;
	}

	public static Collection<SingleInterval> collapse(Collection<SingleInterval> regions) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		Iterator<SingleInterval> iter=regions.iterator();
		
		SingleInterval current=null;
		while(iter.hasNext()) {
			SingleInterval next=iter.next();
			if(overlaps(current, next)) {current=merge(current, next);}
			else {
				if(current!=null) {rtrn.add(current);}
				current=next;
			}
		}
		if(current!=null) {
			rtrn.add(current);
		}
		return rtrn;
		
		/*Map<String, IntervalTree<SingleInterval>> tree=makeTree(regions);
		
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		for(SingleInterval interval: regions) {
			SingleInterval collapsed=getOverlappers(interval, tree);
			rtrn.add(collapsed);
		}
		
		return rtrn;*/
	}

	private static SingleInterval merge(SingleInterval current, SingleInterval next) {
		SingleInterval merged=new SingleInterval(current.getReferenceName(), Math.min(current.getReferenceStartPosition(), next.getReferenceStartPosition()), Math.max(current.getReferenceEndPosition(), next.getReferenceEndPosition()));
		return merged;
	}
	

	private static boolean overlaps(SingleInterval current, SingleInterval next) {
		if(current==null) {return false;}
		return current.overlaps(next);
		
	}

	private static Map<String, IntervalTree<SingleInterval>> makeTree(Collection<SingleInterval> regions) {
		Map<String, IntervalTree<SingleInterval>> rtrn=new TreeMap<String, IntervalTree<SingleInterval>>();
		
		for(SingleInterval region: regions) {
			String chr=region.getReferenceName();
			if(!rtrn.containsKey(chr)) {rtrn.put(chr, new IntervalTree<SingleInterval>());}
			rtrn.get(chr).put(region.getReferenceStartPosition(), region.getReferenceEndPosition(), region);
		}
		
		return rtrn;
	}

	private static SingleInterval getOverlappers(SingleInterval interval, Map<String, IntervalTree<SingleInterval>> tree) {
		String chr=interval.getReferenceName();
		int start=interval.getReferenceStartPosition();
		int end=interval.getReferenceEndPosition();
		
		if(tree.containsKey(chr)){
			Iterator<SingleInterval> iter=tree.get(chr).overlappingValueIterator(start, end);
			
			while(iter.hasNext()) {
				SingleInterval next=iter.next();
				start=Math.min(next.getReferenceStartPosition(), start);
				end=Math.max(next.getReferenceEndPosition(), end);
			}
			
		}
		
		return new SingleInterval(chr, start, end);
	}

	
	public static void main(String[] args) throws IOException {
		
		FileWriter writer=new FileWriter(args[1]);
		Map<String, IntervalTree<Annotation>> tree=BEDFileIO.loadIntergenicRegions(args[0]);
		
		for(String chr: tree.keySet()) {
			System.err.println(chr);
			IntervalTree<Annotation> itree=tree.get(chr);
			Iterator<Annotation> iter=itree.valueIterator();
			while(iter.hasNext()) {
				Annotation g=iter.next();
				writer.write(g.toBED()+"\n");
			}
		}
		writer.close();
		
	}

	public static Map<String, Integer> loadChrSizes(String sizes2) throws IOException {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		Collection<String> lines=loadLines(sizes2);
		
		for(String line: lines) {
			String[] tokens=line.split("\t");
			rtrn.put(tokens[0], Integer.parseInt(tokens[1]));
		}
		return rtrn;
	}

	public static Map<String, IntervalTree<String>> loadRefFlatTree(String fileName) throws IOException {
		
		List<String> lines=loadLines(fileName);
		
		Map<String, Collection<SingleInterval>> map=new TreeMap<String, Collection<SingleInterval>>();
		
		for(String line: lines){
			String chr=line.split("\t")[2];
			int start=Integer.parseInt(line.split("\t")[4]);
			int end=Integer.parseInt(line.split("\t")[5]);
			
			SingleInterval region=new SingleInterval(chr, start, end);
			
			String name=line.split("\t")[0];
			
			
			if(!map.containsKey(name)) {map.put(name, new TreeSet<SingleInterval>());}
			
			Collection<SingleInterval> set=map.get(name);
			set.add(region);
		}
		
		return makeTree(map);
		
	}

	private static Map<String, SingleInterval> remove(Map<String, SingleInterval> map, Collection<String> list) {
		for(String name: list) {
			System.err.println("Removed: "+name);
			map.remove(name);
		}
		return map;
	}

	private static Map<String, IntervalTree<String>> makeTree(Map<String, Collection<SingleInterval>> map) {
		Map<String, IntervalTree<String>> rtrn=new TreeMap<String, IntervalTree<String>>();
		
		for(String name: map.keySet()) {
			Collection<SingleInterval> regions=map.get(name);
			for(SingleInterval region: regions) {
				if(!rtrn.containsKey(region.getReferenceName())) {rtrn.put(region.getReferenceName(), new IntervalTree<String>());}
				IntervalTree<String> tree=rtrn.get(region.getReferenceName());
				tree.put(region.getReferenceStartPosition(), region.getReferenceEndPosition(), name);
			}
		}
		return rtrn;
	}

	public static Collection<Gene> getJunctions(Collection<Gene> genes) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		
		for(Gene gene: genes) {
			rtrn.addAll(gene.getJunctions());
		}
		
		return rtrn;
	}

	public static void writeGeneBED(Collection<Gene> junctions, String string) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		for(Gene r: junctions) {
			writer.write(r.toBED()+"\n");
		}
		
		writer.close();
		
	}

	
	public static Collection<SingleInterval> getIntrons(String file) throws IOException {
		Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(file);
		Collection<SingleInterval> exons=new TreeSet<SingleInterval>();
		
		for(Gene g: genes) {
			Collection<Annotation> blocks=g.getIntrons();
			for(Annotation intron: blocks) {exons.add(intron.getSingleInterval());}
		}
		
		return exons;
	}
	
	public static Collection<SingleInterval> getExons(String file) throws IOException {
		Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(file);
		Collection<SingleInterval> exons=new TreeSet<SingleInterval>();
		
		for(Gene g: genes) {
			Collection<Annotation> blocks=g.getExonSet();
			for(Annotation exon: blocks) {exons.add(exon.getSingleInterval());}
		}
		
		return exons;
	}
	
	public static Collection<SingleInterval> getGeneInterval(String file) throws IOException {
		Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(file);
		Collection<SingleInterval> exons=new TreeSet<SingleInterval>();
		
		for(Gene g: genes) {
			SingleInterval region=g.getSingleInterval();
			region.setOrientation(Strand.POSITIVE);
			exons.add(region);
			//Collection<Annotation> blocks=g.getExonSet();
			//for(Annotation exon: blocks) {exons.add(exon.getSingleInterval());}
		}
		
		return exons;
	}

	public static Collection<Gene> loadTranscriptsFromGTF(File file) throws IOException{
		return new GTFToJunctions(file).getGenes();
	}

	


}
