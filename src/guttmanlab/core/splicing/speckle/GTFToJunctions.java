package guttmanlab.core.splicing.speckle;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.math.Statistics;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class GTFToJunctions {
	
	private Collection<Gene> overlappingJunctions;
	private Collection<Gene> uniqueJunctions;
	private Collection<Gene> allJunctions;
	private Map<String, IntervalTree<Gene>> junctionTree;
	private Map<Gene, Integer> junctionCounts;
	private Map<String, Collection<Gene>> geneToJunctions;
	private Collection<Gene> allGenes;
	//private Map<Gene, String> transcriptToGeneName;
	private File gtfFile;
	private Map<String, Collection<String>> geneType;
	private Map<String, IntervalTree<Gene>> geneCoordinates;
	private Map<String, SingleInterval> geneCoordinateMap;
	
	private Map<SingleInterval, Collection<Integer>> junctionDistance;
	
	

	public GTFToJunctions(File gtfFile) throws IOException {
		this.gtfFile=gtfFile;
		junctionDistance=new TreeMap<SingleInterval, Collection<Integer>>();
		junctionCounts=null;
		allJunctions=getJunctions(gtfFile);
		overlappingJunctions=overlappingJunctions(allJunctions);
		uniqueJunctions=uniqueJunctions(allJunctions, overlappingJunctions);
		makeJunctionTree();
	}
	
	public GTFToJunctions(File gtfFile, File bamFile) throws IOException {
		this.gtfFile=gtfFile;
		Map<SingleInterval, Integer> intronReads=pullSplicedReads(bamFile);
		allJunctions=getJunctions(gtfFile);
		junctionCounts=filter(allJunctions, intronReads);
		allJunctions=junctionCounts.keySet();
		overlappingJunctions=overlappingJunctions(allJunctions);
		uniqueJunctions=uniqueJunctions(allJunctions, overlappingJunctions);
		makeJunctionTree();
	}
	
	public GTFToJunctions(File gtfFile, File bamFile, double percent) throws IOException {
		this.gtfFile=gtfFile;
		Map<SingleInterval, Integer> intronReads=pullSplicedReads(bamFile);
		allJunctions=getJunctions(gtfFile);
		junctionCounts=filter(allJunctions, intronReads);
		allJunctions=this.getJunctions(percent);
		//allJunctions=junctionCounts.keySet();
		overlappingJunctions=overlappingJunctions(allJunctions);
		uniqueJunctions=uniqueJunctions(allJunctions, overlappingJunctions);
		makeJunctionTree();
	}
	
	private void makeJunctionTree() {
		Map<String, IntervalTree<Gene>> tree=new TreeMap<String, IntervalTree<Gene>>();
			
		for(Gene gene: allJunctions){
			IntervalTree<Gene> temp=new IntervalTree<Gene>();
			if(tree.containsKey(gene.getReferenceName())){
				temp=tree.get(gene.getReferenceName());
			}
			temp.put(gene.getReferenceStartPosition(), gene.getReferenceEndPosition(), gene);
			tree.put(gene.getReferenceName(), temp);
		}
		this.junctionTree=tree;
	}

	
	private Map<Gene, Integer> filter(Collection<Gene> allJunctions2, Map<SingleInterval, Integer> intronReads) {
		Map<Gene, Integer> rtrn=new TreeMap<Gene, Integer>();
		for(Gene junction: allJunctions2) {
			SingleInterval intron=junction.getIntrons().iterator().next().getSingleInterval();
			if(intronReads.containsKey(intron)) {
				int counts=intronReads.get(intron);
				rtrn.put(junction, counts);
			}
		}
		return rtrn;
	}

	private static Map<SingleInterval, Integer> pullSplicedReads(File bamFile) {
		SAMFileReader reader=new SAMFileReader(bamFile);
		SAMRecordIterator reads=reader.iterator();
		
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			
			if(SAMFragment.isSpliced(record)) {
				
				SAMFragment frag=new SAMFragment(SAMFragment.compressCIGAR(record));
				Collection<Annotation> introns=frag.getIntrons();
				for(Annotation intron: introns) {
					int count=0;
					SingleInterval intronSI=intron.getSingleInterval();
					if(rtrn.containsKey(intronSI)) {
						count=rtrn.get(intronSI);
					}
					count++;
					rtrn.put(intronSI, count);
				}
			}
				
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		reader.close();
		reads.close();
		
		
		return rtrn;
		
	}
	
	private Collection<Gene> getJunctions(File gtfFile) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(gtfFile)));
		
		Map<String, Map<String, Collection<SingleInterval>>> map=new TreeMap<String, Map<String, Collection<SingleInterval>>>();
		
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			if(!nextLine.startsWith("#")) {
				String[] tokens=nextLine.split("\t");
				if(tokens[2].equals("exon")) {
					SingleInterval exon=new SingleInterval(tokens[0], Integer.parseInt(tokens[3])-1, Integer.parseInt(tokens[4]));
					exon.setOrientation(Strand.fromString(tokens[6]));
					String geneName=getGeneName(tokens[8]);
					String transcript=getTranscriptID(tokens[8]);
					if(transcript!=null) {
						set(map, exon, geneName, transcript);
					}
				}
			}
		}
		reader.close();
		return getJunctions(map);
	}
	
	private Collection<Gene> uniqueJunctions(Collection<Gene> allJunctions2, Collection<Gene> overlappingJunctions2) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		
		for(Gene j: allJunctions2) {
			if(!overlappingJunctions2.contains(j)) {rtrn.add(j);}
		}
		
		return rtrn;
	}

	public Collection<Gene> getUniqueJunctions(){return this.uniqueJunctions;}
	public Collection<Gene> getOverlappingJunctions(){return this.overlappingJunctions;}
	public Collection<Gene> getAllJunctions(){return allJunctions;}
	

	private Collection<Gene> overlappingJunctions(Collection<Gene> junctions) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		
		Map<String, IntervalTree<SingleInterval>> intronTree=makeTree(junctions);
		
		for(Gene junction: junctions) {
			SingleInterval intron=junction.getIntrons().iterator().next().getSingleInterval();
			if(hasOverlaps(intronTree, intron)) {
				//System.out.println(junction.toBED());
				rtrn.add(junction);
			}
			
			
		}
		return rtrn;
		
	}




	private Map<String, IntervalTree<SingleInterval>> makeTree(Collection<Gene> junctions) {
		Map<String, IntervalTree<SingleInterval>> rtrn=new TreeMap<String, IntervalTree<SingleInterval>>();
		
		for(Gene junction: junctions) {
			SingleInterval intron=junction.getIntrons().iterator().next().getSingleInterval();
			if(!rtrn.containsKey(intron.getReferenceName())) {rtrn.put(intron.getReferenceName(), new IntervalTree<SingleInterval>());}
			IntervalTree<SingleInterval> tree=rtrn.get(intron.getReferenceName());
			tree.put(intron.getReferenceStartPosition(), intron.getReferenceEndPosition(), intron);
		}
		
		return rtrn;
	}



	private boolean hasOverlaps(Map<String, IntervalTree<SingleInterval>> intronTree, SingleInterval intron) {
		if(intronTree.containsKey(intron.getReferenceName())) {
			Iterator<SingleInterval> iter=intronTree.get(intron.getReferenceName()).overlappingValueIterator(intron.getReferenceStartPosition(), intron.getReferenceEndPosition());
			while(iter.hasNext()) {
				SingleInterval overlapper=iter.next();
				if(!overlapper.equals(intron)) {return true;}
			}
		}
		return false;
	}



	private Collection<Gene> getJunctions(Map<String, Map<String, Collection<SingleInterval>>> map) {
		this.geneToJunctions=new TreeMap<String, Collection<Gene>>();
		this.allGenes=new TreeSet<Gene>();
		//this.transcriptToGeneName=new TreeMap<Gene, String>();
		Collection<Gene> set=new TreeSet<Gene>();
		Map<SingleInterval, Collection<Gene>> rtrn=new TreeMap<SingleInterval, Collection<Gene>>();
		for(String gene: map.keySet()) {
			for(String transcript: map.get(gene).keySet()) {
				Collection<SingleInterval> exons=map.get(gene).get(transcript);
				Gene g=new Gene(exons, gene);
				allGenes.add(g);
				//this.transcriptToGeneName.put(g, gene);
				Collection<Gene> junctions=g.getJunctions();
				for(Gene j: junctions) {
					j.setName(gene);
					SingleInterval intron=j.getIntrons().iterator().next().getSingleInterval();
					int distance=getDistance(g, j);
					addDistance(intron, distance);
					set(rtrn, intron, j);
				}
				//System.out.println(g.toBED());
			}
		}
		
		for(SingleInterval intron: rtrn.keySet()) {
			Collection<Gene> junctions=rtrn.get(intron);
			Gene collapsed=collapse(junctions);
			Collection<String> names=getAllNames(junctions);
			set.add(collapsed);
			add(geneToJunctions, names, collapsed);
			//System.out.println(collapsed.toBED());
		}
		
		return set;
	}

	


	

	private int getDistance(Gene g, Gene j) {
		if(g.getOrientation().equals(Strand.POSITIVE)) {
			int ss3=j.getIntrons().iterator().next().getReferenceEndPosition();
			int tss=g.getReferenceStartPosition();
			return ss3-tss;
		}
		
		int tss=g.getReferenceEndPosition();
		int ss3=j.getIntrons().iterator().next().getReferenceEndPosition();
		return tss-ss3;
	}

	private void addDistance(SingleInterval j, int distance) {
		if(!this.junctionDistance.containsKey(j)) {
			this.junctionDistance.put(j, new ArrayList<Integer>());
		}
		
		Collection<Integer> list=this.junctionDistance.get(j);
		list.add(distance);
	}

	private Collection<String> getAllNames(Collection<Gene> junctions) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(Gene g: junctions) {rtrn.add(g.getName());}
		
		return rtrn;
	}

	private void add(Map<String, Collection<Gene>> geneToJunctions, Collection<String> names, Gene collapsed) {
		for(String name: names) {
			if(!geneToJunctions.containsKey(name)) {
				geneToJunctions.put(name, new TreeSet<Gene>());
			}
			Collection<Gene> list=geneToJunctions.get(name);
			list.add(collapsed);
		}
	}


	private Gene collapse(Collection<Gene> junctions) {
		if(junctions.size()==1) {return junctions.iterator().next();}
		
		Collection<Annotation> blocks=new TreeSet<Annotation>();
		int start=-Integer.MAX_VALUE;
		int end=Integer.MAX_VALUE;
		String name="";
		for(Gene g: junctions) {
			blocks.addAll(g.getExonSet());
			start=Math.max(start, g.getReferenceStartPosition());
			end=Math.min(end, g.getReferenceEndPosition());
			name=g.getName();
		}
		
		Gene rtrn=new Gene(blocks, name);
		rtrn=new Gene(rtrn.trim(start, end));
		rtrn.setName(name);
		return rtrn;
	}

	private void set(Map<SingleInterval, Collection<Gene>> rtrn, SingleInterval intron, Gene j) {
		if(!rtrn.containsKey(intron)) {rtrn.put(intron, new TreeSet<Gene>());}
		Collection<Gene> list=rtrn.get(intron);
		list.add(j);
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

	private static String getTranscriptID(String string) {
		return getTag(string, "transcript_id");
	}

	private void set(Map<String, Map<String, Collection<SingleInterval>>> map, SingleInterval exon, String geneName, String transcript) {
		if(!map.containsKey(geneName)) {
			map.put(geneName, new TreeMap<String, Collection<SingleInterval>>());
		}
		
		Map<String, Collection<SingleInterval>> rtrn=map.get(geneName);
		if(!rtrn.containsKey(transcript)) {
			rtrn.put(transcript, new TreeSet<SingleInterval>());
		}
		
		Collection<SingleInterval> set=rtrn.get(transcript);
		set.add(exon);
	}
	
	

	public boolean isUniqueJunction(Gene gene) {
		if(this.getUniqueJunctions().contains(gene)) {return true;}
		return false;
	}

	public Map<String, IntervalTree<Gene>> getJunctionTree() {
		return junctionTree;
	}


	public double getJunctionCount(Gene gene) {
		if(junctionCounts==null || !junctionCounts.containsKey(gene)) {return -1;}
		return junctionCounts.get(gene);
	}
	
	public Map<String, Collection<Gene>> getGeneToJunctions() {
		return this.geneToJunctions;
		
	}
	
	private Collection<Gene> getJunctions(double percent){
		Collection<Gene> rtrn=new TreeSet<Gene>();
		for(String gene: geneToJunctions.keySet()) {
			Collection<Gene> junctions=geneToJunctions.get(gene);
			ArrayList<Double> vals=new ArrayList<Double>();
			for(Gene jun: junctions) {
				double count=getJunctionCount(jun);
				if(count>-1) {
					vals.add(count);
					//System.out.println(gene+"\t"+jun.toUCSC()+"\t"+count);
				}
			}
			if(!vals.isEmpty()) { 
				//double median=Statistics.quantile(vals, 0.5);
				double max=Statistics.max(vals);
				//System.out.println(gene+"\t"+vals.size()+"\t"+Statistics.mean(vals)+"\t"+Statistics.max(vals));
				for(Gene jun: junctions) {
					double count=getJunctionCount(jun);
					double threshold=percent*max;
					if(count>threshold) {
						rtrn.add(jun);
						//System.out.println(gene+"\t"+jun.toUCSC()+"\t"+count);
					}
				}
			}
		}
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		File file=new File(args[0]);
		GTFToJunctions gtf=	new GTFToJunctions(file, new File(args[1]));
		Collection<Gene> list=gtf.getJunctions(0.05);
		for(Gene jun: list) {System.out.println(jun.toBED());}
		
	}

	public Collection<Gene> getGenes() {
		return this.allGenes;
	}
	
	
	private void makeGeneCoordinates() throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(gtfFile)));
			
		Map<String, SingleInterval> rtrn=new TreeMap<String, SingleInterval>();
			
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			if(!nextLine.startsWith("#")) {
				String[] tokens=nextLine.split("\t");
				if(tokens[2].equals("gene")) {
					SingleInterval gene=new SingleInterval(tokens[0], Integer.parseInt(tokens[3])-1, Integer.parseInt(tokens[4]));
					gene.setOrientation(Strand.fromString(tokens[6]));
					String geneName=getGeneName(tokens[8]);
					gene.setName(geneName);
					rtrn.put(geneName, gene);
					//System.out.println(gene.toBED());
				}
			}
		}
		reader.close();
		this.geneCoordinateMap=rtrn;
	}
	
	
	public SingleInterval getGeneCoordinates(String name) {
		if(geneCoordinateMap==null || geneCoordinateMap.isEmpty()) {
			try {makeGeneCoordinates();} catch (IOException e) {e.printStackTrace();}
		}
		
		return geneCoordinateMap.get(name);
	}

	public Map<String, IntervalTree<Gene>> getGeneTree() {
		Map<String, IntervalTree<Gene>> rtrn=new TreeMap<String, IntervalTree<Gene>>();
		
		for(Gene gene: allGenes){
			String chr=gene.getReferenceName();
			if(!rtrn.containsKey(chr)) {rtrn.put(chr, new IntervalTree<Gene>());}
			IntervalTree<Gene> tree=rtrn.get(chr);
			tree.put(gene.getReferenceStartPosition(), gene.getReferenceEndPosition(), gene);
		}
		
		return rtrn;
	}
	
	public Map<String, IntervalTree<SingleInterval>> getGeneCoordinateTree() throws IOException {
		Map<String, IntervalTree<SingleInterval>> rtrn=new TreeMap<String, IntervalTree<SingleInterval>>();
		
		Collection<SingleInterval> allGenes=getGeneCoordinates();
		for(SingleInterval gene: allGenes){
			String chr=gene.getReferenceName();
			if(!rtrn.containsKey(chr)) {rtrn.put(chr, new IntervalTree<SingleInterval>());}
			IntervalTree<SingleInterval> tree=rtrn.get(chr);
			tree.put(gene.getReferenceStartPosition(), gene.getReferenceEndPosition(), gene);
		}
		
		return rtrn;
	}

	private Collection<SingleInterval> getGeneCoordinates() throws NumberFormatException, IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(gtfFile)));
		
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
			
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			if(!nextLine.startsWith("#")) {
				String[] tokens=nextLine.split("\t");
				if(tokens[2].equals("gene")) {
					SingleInterval gene=new SingleInterval(tokens[0], Integer.parseInt(tokens[3])-1, Integer.parseInt(tokens[4]));
					gene.setOrientation(Strand.fromString(tokens[6]));
					String geneName=getGeneName(tokens[8]);
					gene.setName(geneName);
					rtrn.add(gene);
					//System.out.println(gene.toBED());
				}
			}
		}
		reader.close();
		return rtrn;
	}

	public String getGeneType(String name) {
		if(this.geneType==null || this.geneType.isEmpty()) {
			try {makeGeneType();} catch (IOException e) {e.printStackTrace();}
		}
		return consensusType(name);
	}

	private void makeGeneType() throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(gtfFile)));
		
		this.geneType=new TreeMap<String, Collection<String>>();
			
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			if(!nextLine.startsWith("#")) {
				String[] tokens=nextLine.split("\t");
				if(tokens[2].equals("gene")) {
					//SingleInterval gene=new SingleInterval(tokens[0], Integer.parseInt(tokens[3])-1, Integer.parseInt(tokens[4]));
					//gene.setOrientation(Strand.fromString(tokens[6]));
					String geneName=getGeneName(tokens[8]);
					String type=getTag(tokens[8], "gene_type");
					
					//System.out.println(geneName+"\t"+type);
					
					if(!this.geneType.containsKey(geneName)) {this.geneType.put(geneName, new TreeSet<String>());}
					Collection<String> set=geneType.get(geneName);
					set.add(type);
				}
			}
		}
		reader.close();
	}

	private String consensusType(String name) {
		Collection<String> types=this.geneType.get(name);
		//System.out.println(name);
		if(types.size()>1) {System.out.println(name+"\t"+types.size());}
		
		return types.iterator().next();
	}

	public Collection<String> getGene(SingleInterval intron) {
		Collection<String> rtrn=new TreeSet<String>();
		if(this.junctionTree.containsKey(intron.getReferenceName())) {
			Iterator<Gene> iter=this.junctionTree.get(intron.getReferenceName()).overlappingValueIterator(intron.getReferenceStartPosition(), intron.getReferenceEndPosition());
			while(iter.hasNext()) {
				rtrn.add(iter.next().getName());
			}
		}
		return rtrn;
	}
	
	public Map<SingleInterval, Collection<Integer>> getJunctionDistance(){return this.junctionDistance;}

	public Collection<SingleInterval> getIntrons() {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		Collection<Gene> junctions= getAllJunctions();
		for(Gene g: junctions) {
			SingleInterval intron=g.getIntrons().iterator().next().getSingleInterval();
			rtrn.add(intron);
		}
		return rtrn;
	}

	public static Collection<Gene> getTranscripts(File gtfFile) throws IOException{
		return new GTFToJunctions(gtfFile).getGenes();
	}
	
	
	
	

}
