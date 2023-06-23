package guttmanlab.core.snps;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.math3.stat.inference.ChiSquareTest;

import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.IntervalTree.Node;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.sequence.Sequence;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class QuantifyBySNP {

	static String ambiguous="amb";
	Map<SingleInterval, Pair<Integer>> snpCounts;
	
	
	
	
	
	
	public QuantifyBySNP(File bamFile, Map<String, IntervalTree<String>> snpTree, String fileName) throws IOException{
		
		
		SAMFileReader reader=new SAMFileReader(bamFile);
		
		
		SAMFileHeader fileHeader=reader.getFileHeader();
		
		SAMFileWriter writer1=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(fileHeader, false, new File(fileName+".SNP1.bam"));
		SAMFileWriter writer2=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(fileHeader, false, new File(fileName+".SNP2.bam"));
		SAMFileWriter writer3=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(fileHeader, false, new File(fileName+".amb.bam"));
		
		
		
		SAMRecordIterator reads=reader.iterator();
		
		Map<String, Pair<SAMRecord>> set1=new TreeMap<String, Pair<SAMRecord>>();
		Map<String, Pair<SAMRecord>> set2=new TreeMap<String, Pair<SAMRecord>>();
		Map<String, Pair<SAMRecord>> set3=new TreeMap<String, Pair<SAMRecord>>();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			String genotype=getOverlappingRegion(record, snpTree);
			
			
			if(genotype.equals("SNP1")){
				add(set1, set3, record, writer1);
			}
			else if(genotype.equals("SNP2")){
				add(set2, set3, record, writer2);
			}
			else{
				addAmbiguous(record, set1, set2, set3, writer1, writer2, writer3);
			}
			
			
			
			counter++;
			if(counter%1000000 ==0){System.err.println(counter+" "+set1.size()+" "+set2.size()+" "+set3.size());}
		}
		
		System.err.println(set1.size()+" "+set2.size()+" "+set3.size());
		//mergeSets(set1, set2, set3, writer1, writer2, writer3);
		
		write(set1, writer1);
		write(set2, writer2);
		write(set3, writer3);
		
		
		reader.close();
		reads.close();
		writer1.close();
		writer2.close();
		writer3.close();
		
		
		this.snpCounts=getSNPCounts(new File(fileName+".SNP1.bam"), new File(fileName+".SNP2.bam"), snpTree);
	}
	
	
	
	public QuantifyBySNP(File bamFile1, File bamFile2, Map<String, IntervalTree<String>> snpTree, String save) throws IOException{
		
		List<String> rows=getRows(snpTree);
		List<String> columns=new ArrayList<String>();
		columns.add("snp1");
		columns.add("snp2");
		MatrixWithHeaders data=new MatrixWithHeaders(rows, columns);
		
		
		
		
		snpCount(bamFile1, snpTree, data, "snp1");
		snpCount(bamFile2, snpTree, data, "snp2");
		
		
		data.write(save);
	}
	
	
	
	private List<String> getRows(Map<String, IntervalTree<String>> snpTree) {
		List<String> rtrn=new ArrayList<String>();
		
		for(String chr: snpTree.keySet()) {
			Iterator<Node<String>> iter=snpTree.get(chr).iterator();
			while(iter.hasNext()) {
				Node<String> n=iter.next();
				String name=chr+":"+n.getStart()+"-"+n.getEnd()+"_"+n.getValue();
				rtrn.add(name);
			}
		}
		
		return rtrn;
	}



	private void snpCount(File bamFile1, Map<String, IntervalTree<String>> snpTree, MatrixWithHeaders data, String col) {
		SAMFileReader reader=new SAMFileReader(bamFile1);
		
		SAMRecordIterator reads=reader.iterator();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			Collection<String> genotype=getOverlappingSNP(record, snpTree);
			
			for(String row: genotype) {
				data.incrementCount(row, col);
			}
			
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		reader.close();
		reads.close();
		
	}



	public QuantifyBySNP(File bamFile, Map<String, IntervalTree<String>> snpTree, String fileName, boolean rewrite) throws IOException{	
		if(rewrite){new QuantifyBySNP(bamFile, snpTree, fileName);}
		else{this.snpCounts=getSNPCounts(new File(fileName+".SNP1.bam"), new File(fileName+".SNP2.bam"), snpTree);}
	}
	
	
	
	
	private static void write(String string, Map<SingleInterval, Pair<Integer>> snpCounts, int minCount) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		for(SingleInterval region: snpCounts.keySet()){
			double total=(snpCounts.get(region).getValue1()+snpCounts.get(region).getValue2());
			if(total>minCount){
				double ratio=(double)snpCounts.get(region).getValue1()/total;
				writer.write(region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+region.getReferenceEndPosition()+"\t"+ratio+"\n");
			}
		}
		
		writer.close();
	}


	private Map<SingleInterval, Pair<Integer>> getSNPCounts(File file1, File file2, Map<String, IntervalTree<String>> snpTree) {
		Map<SingleInterval, Pair<Integer>> rtrn=new TreeMap<SingleInterval, Pair<Integer>>();
		
		
		updateValues(file1, snpTree, rtrn, 1);
		updateValues(file2, snpTree, rtrn, 2);
		
		
		return rtrn;
	}


	private void updateValues(File file1, Map<String, IntervalTree<String>> snpTree, Map<SingleInterval, Pair<Integer>> rtrn, int i) {
		SAMFileReader reader=new SAMFileReader(file1);
		SAMRecordIterator reads=reader.iterator();
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			Collection<SingleInterval> regions=getSNPs(record, snpTree);
			for(SingleInterval region: regions){
				Pair<Integer> temp=new Pair<Integer>(0,0);
				if(rtrn.containsKey(region)){temp=rtrn.get(region);}
				int score=temp.getValue1()+1;
				if(i==2){
					score=temp.getValue2()+1;
					temp.setValue2(score);
				}
				else{temp.setValue1(score);}
				rtrn.put(region, temp);
			}
		}
		
		reader.close();
		reads.close();
		
	}


	private Collection<SingleInterval> getSNPs(SAMRecord record, Map<String, IntervalTree<String>> snpTree) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		String chr=getChromosome(record.getReferenceName());
		if(!snpTree.containsKey(chr)){return rtrn;}
		
		Iterator<Node<String>> snps=snpTree.get(chr).overlappers(record.getAlignmentStart(), record.getAlignmentEnd());
				
		while(snps.hasNext()){
			Node<String> g=snps.next();
			SingleInterval region=new SingleInterval(chr, g.getStart(), g.getEnd());
			rtrn.add(region);
		}
		
		return rtrn;
	}


	private void write(Map<String, Pair<SAMRecord>> set1, SAMFileWriter writer1) {
		for(String readName: set1.keySet()){
			Pair<SAMRecord> records=set1.get(readName);
			if(records.hasValue1()){writer1.addAlignment(records.getValue1());}
			if(records.hasValue2()){writer1.addAlignment(records.getValue2());}
		}
		
	}


	private void addAmbiguous(SAMRecord record, Map<String, Pair<SAMRecord>> set1, Map<String, Pair<SAMRecord>> set2, Map<String, Pair<SAMRecord>> set3, SAMFileWriter writer1, SAMFileWriter writer2, SAMFileWriter writer3) {
		//Pair<SAMRecord> pair=new Pair<SAMRecord>(null, null);
		
		if(set1.containsKey(record.getReadName())){
			add(set1, set3, record, writer1);
		}
		
		if(set2.containsKey(record.getReadName())){
			add(set2, set3, record, writer2);
		}
	
		else{
			add(set3, set3, record, writer3);
		}
		
	}


	/*public SplitBySNP(File bamFile, File snpFile, String fileName) throws IOException{
		this.files=new TreeMap<String, String>();
		SAMFileReader reader=new SAMFileReader(bamFile);
		Map<String, IntervalTree<String>> snpTree=parseTree(snpFile);
		
		
		SAMRecordIterator reads=reader.iterator();
		
		Map<String, Pair<SAMRecord>> set1=new TreeMap<String, Pair<SAMRecord>>();
		Map<String, Pair<SAMRecord>> set2=new TreeMap<String, Pair<SAMRecord>>();
		Map<String, Pair<SAMRecord>> set3=new TreeMap<String, Pair<SAMRecord>>();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			//System.err.println(record.getReadName());
			String genotype=getOverlappingRegion(record, snpTree);
			
			
			if(genotype.equals("SNP1")){
				add(set1, record);
			}
			else if(genotype.equals("SNP2")){
				add(set2, record);
			}
			else{
				add(set3, record);
			}
			
			
			
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		mergeSets(set1, set2, set3, fileName, reader.getFileHeader(), "SNP1", "SNP2");
		
		
		reader.close();
		reads.close();
	}*/
	
	public static Map<String, IntervalTree<String>> parseTree(File snpFile) throws NumberFormatException, IOException {
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
	
	
	
	
	private void add(Map<String, Pair<SAMRecord>> set1, Map<String, Pair<SAMRecord>> set3, SAMRecord record, SAMFileWriter writer) {
		Pair<SAMRecord> pair=new Pair<SAMRecord>(null, null);
		if(set1.containsKey(record.getReadName())){pair=set1.remove(record.getReadName());}
		else if(set3.containsKey(record.getReadName())){pair=set3.remove(record.getReadName());}
		boolean first=record.getFirstOfPairFlag();
		if(first){pair.setValue1(record);}
		else{pair.setValue2(record);}
		
		if(pair.hasValue1() && pair.hasValue2()){
			writer.addAlignment(pair.getValue1());
			writer.addAlignment(pair.getValue2());
			//System.err.println("Write to file");
		}
		
		else{
			set1.put(record.getReadName(), pair);
		}
	}


	private void mergeSets(Map<String, Pair<SAMRecord>> set1, Map<String, Pair<SAMRecord>> set2,Map<String, Pair<SAMRecord>> set3, SAMFileWriter writer1, SAMFileWriter writer2, SAMFileWriter writer3) {
		Map<String, Pair<SAMRecord>> temp=new TreeMap<String, Pair<SAMRecord>>();
		for(String readName: set3.keySet()){
			Pair<SAMRecord> pair=set3.get(readName);
			if(set1.containsKey(readName)){
				Pair<SAMRecord> merged=merge(pair, set1.get(readName));
				set1.put(readName, merged);
			}
			else if(set2.containsKey(readName)){
				Pair<SAMRecord> merged=merge(pair, set2.get(readName));
				set2.put(readName, merged);
			}
			else{
				temp.put(readName, pair);
			}
		}
		
		set3=temp;
		
		for(String readName: set1.keySet()){
			Pair<SAMRecord> pair=set1.get(readName);
			if(pair.getValue1()!=null){writer1.addAlignment(pair.getValue1());}
			if(pair.getValue2()!=null){writer1.addAlignment(pair.getValue2());}
		}
		
		for(String readName: set2.keySet()){
			Pair<SAMRecord> pair=set2.get(readName);
			if(pair.getValue1()!=null){writer2.addAlignment(pair.getValue1());}
			if(pair.getValue2()!=null){writer2.addAlignment(pair.getValue2());}
		}
		
		for(String readName: set3.keySet()){
			Pair<SAMRecord> pair=set3.get(readName);
			if(pair.getValue1()!=null){writer3.addAlignment(pair.getValue1());}
			if(pair.getValue2()!=null){writer3.addAlignment(pair.getValue2());}
		}
		
	}


	private Pair<SAMRecord> merge(Pair<SAMRecord> pair, Pair<SAMRecord> pair2) {
		Pair<SAMRecord> rtrn=new Pair<SAMRecord>();
		SAMRecord record1=getVal(pair.getValue1(), pair2.getValue1());
		SAMRecord record2=getVal(pair.getValue2(), pair2.getValue2());
		rtrn.setValue1(record1);
		rtrn.setValue2(record2);
		return rtrn;
	}


	private SAMRecord getVal(SAMRecord val1, SAMRecord val2) {
		if(val1==null){return val2;}
		return val1;
	}


	private String getOverlappingRegion(SAMRecord record, Map<String, IntervalTree<String>> tree) {
		//CloseableIterator<VariantContext> iter=snps.query(record.getReferenceName().replace("chr", ""), record.getAlignmentStart(), record.getAlignmentEnd());
		
		String chr=getChromosome(record.getReferenceName());
		if(!tree.containsKey(chr)){return ambiguous;}
		Iterator<Node<String>> snps=tree.get(chr).overlappers(record.getAlignmentStart(), record.getAlignmentEnd());
		
		Collection<String> rtrn=new TreeSet<String>();
		
		while(snps.hasNext()){
			Node<String> g=snps.next();
			String g1=g.getValue().split(",")[0];
			String g2=g.getValue().split(",")[1];
			String base=getSequence(record, g, record.getReferenceName());
			if(base.equals(g1)){rtrn.add("SNP1");}
			if(base.equals(g2)){rtrn.add("SNP2");}
		}
		
		
		if(rtrn.size()==1){return rtrn.iterator().next();}
		//else if(rtrn.size()>1){System.err.println(rtrn.size());}
		return ambiguous;
	}
	
	
	private Collection<String> getOverlappingSNP(SAMRecord record, Map<String, IntervalTree<String>> tree) {
		Collection<String> rtrn=new ArrayList<String>();
		String chr=getChromosome(record.getReferenceName());
		if(!tree.containsKey(chr)){return rtrn;}
		Iterator<Node<String>> snps=tree.get(chr).overlappers(record.getAlignmentStart(), record.getAlignmentEnd());
		
		
		while(snps.hasNext()){
			Node<String> n=snps.next();
			String name=chr+":"+n.getStart()+"-"+n.getEnd()+"_"+n.getValue();
			rtrn.add(name);
		}
		
		return rtrn;
	}
	
	private String getSequence(SAMRecord record, Node<String> genotype, String chr) {
		SAMFragment fragment=new SAMFragment(record);
		char seq=fragment.getSequenceAtPosition(chr, genotype.getStart()-1);
		String val=seq+"/"+seq;
		//if(seq!=' '){System.err.println(variant.getStart()+" "+val);}
		return val;
		
	}
	
	private String getChromosome(String referenceName) {
		if(referenceName.contains("chr")){return referenceName;}
		return "chr"+referenceName;
	}
	
	private String getOverlappingRegion(SAMRecord record, VCFFileReader snps, String sp1, String sp2) {
		CloseableIterator<VariantContext> iter=snps.query(record.getReferenceName().replace("chr", ""), record.getAlignmentStart(), record.getAlignmentEnd());
		
		Collection<String> rtrn=new TreeSet<String>();
		
		while(iter.hasNext()){
			VariantContext variant=iter.next();
			//System.err.println(variant);
			boolean isHet=isHetero(variant, sp1, sp2);
			/*if(isHet){
				System.err.println(variant.getStart()+" "+variant.getGenotype(sp1) +" "+variant.getGenotype(sp2));
			}*/
			if(isHet){
				for(Genotype g: variant.getGenotypes()){
					String name=g.getSampleName();
					if(name.equals(sp1) || name.equals(sp2)){
						String genotype=g.getGenotypeString();
						boolean matches=checkRead(record, variant, genotype);
						if(matches){rtrn.add(name);}
					}
				}
			}
		}
		
		iter.close();
		if(rtrn.size()==1){return rtrn.iterator().next();}
		//else if(rtrn.size()>1){System.err.println(rtrn.size());}
		return ambiguous;
	}


	private boolean isHetero(VariantContext variant, String sp1, String sp2) {
		return !variant.getGenotype(sp1).sameGenotype(variant.getGenotype(sp2));
	}


	private boolean checkRead(SAMRecord record, VariantContext variant, String genotype) {
		SAMFragment fragment=new SAMFragment(record);
		char seq=fragment.getSequenceAtPosition("chr"+variant.getChr(), variant.getStart()-1);
		String val=seq+"/"+seq;
		//if(seq!=' '){System.err.println(variant.getStart()+" "+val);}
		return val.equals(genotype);
		
	}


	/*public SplitBySNP(File bamFile, Map<SingleInterval, Pair<String>> snps, String fileName) throws IOException{
		SAMFileReader reader=new SAMFileReader(bamFile);
		
		SAMRecordIterator reads=reader.iterator();
		
		SAMFileWriter writer1=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(reader.getFileHeader(), false, new File(fileName+".SNP1.bam"));
		SAMFileWriter writer2=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(reader.getFileHeader(), false, new File(fileName+".SNP2.bam"));
		SAMFileWriter writer3=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(reader.getFileHeader(), false, new File(fileName+".amb.bam"));
		
		
		Map<SingleInterval, Pair<Integer>> counts=new TreeMap<SingleInterval, Pair<Integer>>();
		
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			//SAMFragment fragment=new SAMFragment(record);
			//fragment.getRelativePositionFrom5PrimeOfFeature(referenceStart)
			//String sequence=record.getReadString();
			SingleInterval overlappingRegion=getOverlappingRegion(record, snps);
			String SNP=getStringAtPosition(record, snps);
			Pair<Integer> pairCount=new Pair<Integer>(0,0);
			if(overlappingRegion!=null &&counts.containsKey(overlappingRegion)){pairCount=counts.get(overlappingRegion);}
			
			if(SNP.equals("SNP1")){
				writer1.addAlignment(record);
				pairCount.setValue1(pairCount.getValue1()+1);
			}
			else if(SNP.equals("SNP2")){
				writer2.addAlignment(record);
				pairCount.setValue2(pairCount.getValue2()+1);
			}
			else{writer3.addAlignment(record);}
			if(overlappingRegion!=null){counts.put(overlappingRegion, pairCount);}
		}
		
		
		FileWriter writer=new FileWriter(fileName);
		
		for(SingleInterval region: counts.keySet()){
			writer.write(region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+counts.get(region).getValue1()+"\t"+counts.get(region).getValue2()+"\n");
		}
		
		
		writer.close();
		reader.close();
		reads.close();
		writer1.close();
		writer2.close();
		writer3.close();
		
	}*/
	
	
	
	/*private SingleInterval getOverlappingRegion(SAMRecord record, Map<SingleInterval, Pair<String>> snps) {
		SAMFragment fragment=new SAMFragment(record);
		
		for(SingleInterval region: snps.keySet()){
			if(fragment.overlaps(region)){
				return region;
			}
		}
		return null;
	}*/



	private String getStringAtPosition(SAMRecord record, Map<SingleInterval, Pair<String>> snps) {
		SAMFragment fragment=new SAMFragment(record);
		
		for(SingleInterval region: snps.keySet()){
			//SingleInterval region2=new SingleInterval(region.getReferenceName(), region.getReferenceStartPosition()-1, region.getReferenceEndPosition());
			if(fragment.overlaps(region)){
				Pair<String> snp=snps.get(region);
				char seq=getPosition(region, fragment);
				if(seq!=' '){
					System.err.println(seq+" "+snp.getValue1()+" "+snp.getValue2());
					if(seq == (snp.getValue1().toCharArray()[0])){return "SNP1";}
					if(seq == snp.getValue2().toCharArray()[0]){return "SNP2";}
				}
			}
		}
		return "";
	}

	/*private String getStringAtPosition(SAMRecord record, Map<SingleInterval, Pair<String>> snps) {
		for(SingleInterval region: snps.keySet()){
			if(overlaps(record, region)){
				Pair<String> snp=snps.get(region);
				char seq=getPosition(region, record);
				if(seq!=' '){
					System.err.println(seq+" "+snp.getValue1()+" "+snp.getValue2());
					if(seq == (snp.getValue1().toCharArray()[0])){return "SNP1";}
					if(seq == snp.getValue2().toCharArray()[0]){return "SNP2";}
				}
			}
		}
		return "";
	}*/

	private boolean overlaps(SAMRecord record, SingleInterval region) {
		SingleInterval temp=new SingleInterval("chr"+record.getReferenceName(), record.getAlignmentStart()-1, record.getAlignmentEnd());
		//System.err.println(temp.toUCSC()+" "+region.toUCSC());
		return temp.overlaps(region);
	}
	

	private char getPosition(SingleInterval region, SAMRecord record) {
		int offset=region.getReferenceStartPosition()-record.getAlignmentStart();
		//System.err.println(region.getReferenceStartPosition()+" "+record.getAlignmentStart()+" "+record.getAlignmentEnd()+" "+record.getCigarString()+" "+offset+" "+record.getReadString().length());
		if(offset<record.getReadLength()){
		char[] seq=record.getReadString().toCharArray();
		return seq[offset];
		}
		return ' ';
	}
	
	private char getPosition(SingleInterval region, SAMFragment record) {
		int offset=record.getRelativePositionFrom5PrimeOfFeature(region.getReferenceStartPosition())-1;
		//System.err.println(region.getReferenceStartPosition()+" "+offset);
		char[] seq=record.getSamRecord().getReadString().toCharArray();
		return seq[offset];
	}
	
	private char getPosition(VariantContext region, SAMFragment record) {
		
		
		char[] seq=record.getSamRecord().getReadString().toCharArray();
		for(int i=0; i<seq.length; i++){
			int pos=record.getSamRecord().getReferencePositionAtReadPosition(i);
			if(pos==region.getStart()){return seq[i];}
		}
		
		return ' ';
	}
	
	private static Map<SingleInterval, Pair<String>> parse(String string) throws IOException {
		List<String> lines=BEDFileIO.loadLines(string);
		Map<SingleInterval, Pair<String>> rtrn=new TreeMap<SingleInterval, Pair<String>>();
		
		int counter=0;
		for(String line: lines){
			if(counter>0){
				String[] tokens=line.split("\t");
				String chr=tokens[1];
				int start=new Integer(tokens[2]);
				int end=start+1;
				Strand st=Strand.fromString(tokens[3]);
				Pair<String> snps=new Pair<String>(tokens[4], tokens[5]);
				SingleInterval region=new SingleInterval(chr, start, end);
				//if(st.equals(Strand.NEGATIVE)){snps=reverseComplement(snps);}
				rtrn.put(region, snps);
			}
			counter++;
		}
		
		return rtrn;
	}


	
	
	public Map<SingleInterval, Pair<Integer>> getSNPCounts() {
		return this.snpCounts;
	}

	
	
	private static void ratio(String save, Map<SingleInterval, Pair<Integer>> sampleSnpCounts, Map<SingleInterval, Pair<Integer>> controlSnpCounts, int minCount) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: sampleSnpCounts.keySet()){
			Pair<Integer> sampleVals=sampleSnpCounts.get(region);
			if(controlSnpCounts.containsKey(region)){
				Pair<Integer> controlVals=controlSnpCounts.get(region);
				if(sum(sampleVals)>minCount && sum(controlVals)>minCount){
					double ratio=ratio(sampleVals, controlVals);
					double p=hypergeometric(sampleVals, controlVals, minCount);
					if(p<0.05){
					writer.write(region.toBedgraph(ratio)+"\n");
					}
				}
			}
		}
		
		writer.close();
	}
	
	private static double hypergeometric(Pair<Integer> sampleVal, Pair<Integer> controlVal, int minCount) {
		ChiSquareTest test=new ChiSquareTest();
		
		long[] observed1=new long[2];
		long[] observed2=new long[2];
		
		observed1[0]=sampleVal.getValue1();
		observed1[1]=sampleVal.getValue2();
		
		observed2[0]=controlVal.getValue1();
		observed2[1]=controlVal.getValue2();
		
		double ratio1=((double)sampleVal.getValue1())/((double)sum(sampleVal));
		double ratio2=((double)controlVal.getValue1())/((double)sum(controlVal));
		
		double p=1.0;
		if(ratio1==0 && ratio2==0){
			p=1.0;
		}
		else if(ratio1==1 && ratio2==1){
			p=1.0;
		}
		else if(sum(sampleVal)>minCount && sum(controlVal)>minCount){
			p=test.chiSquareTestDataSetsComparison(observed1, observed2);
		}
		
		
		return p;
	}
	
	private static double ratio(Pair<Integer> sampleVals, Pair<Integer> controlVals) {
		double ratio1=ratio(sampleVals);
		double ratio2=ratio(controlVals);
		
		return Math.log(ratio1/ratio2)/Math.log(2);
		
		/*double expected=((double)controlVals.getValue1()/(sum(controlVals)))*sum(sampleVals);
		double enrich=(double)sampleVals.getValue1()/expected;
		return enrich;*/
	}


	private static double ratio(Pair<Integer> sampleVals) {
		return ((double)sampleVals.getValue1()+1.0)/((double)sum(sampleVals)+2.0);
	}


	private static int sum(Pair<Integer> controlVals) {
		return controlVals.getValue1()+controlVals.getValue2();
	}

	

	private static void mergeGenes(AnnotationCollection<Gene> genes, Map<SingleInterval, Pair<Integer>> sampleSnpCounts, Map<SingleInterval, Pair<Integer>> controlSnpCounts, String save, int minCount) throws IOException {
		Map<Gene, Pair<Integer>> sampleGeneCounts=getCountsPerGene(genes, sampleSnpCounts);
		Map<Gene, Pair<Integer>> controlGeneCounts=getCountsPerGene(genes, controlSnpCounts);
		
		FileWriter writer=new FileWriter(save);
		
		for(Gene gene: sampleGeneCounts.keySet()){
			Pair<Integer> sampleVals=sampleGeneCounts.get(gene);
			if(controlGeneCounts.containsKey(gene)){
				Pair<Integer> controlVals=controlGeneCounts.get(gene);
				if(sum(sampleVals)>minCount && sum(controlVals)>minCount){
					//System.err.println(gene.getName()+" "+sum(sampleVals)+" "+sum(controlVals));
					double ratio=ratio(sampleVals, controlVals);
					double p=hypergeometric(sampleVals, controlVals, minCount);
					if(p<0.05){
						writer.write(gene.getReferenceName()+"\t"+gene.getReferenceStartPosition()+"\t"+gene.getReferenceEndPosition()+"\t"+ratio+"\n");
					}
				}
			}
		}
		
		writer.close();
		
		
	}

	private static Map<Gene, Pair<Integer>> getCountsPerGene(AnnotationCollection<Gene> genes, Map<SingleInterval, Pair<Integer>> sampleSnpCounts) {
		Map<Gene, Pair<Integer>> rtrn=new TreeMap<Gene, Pair<Integer>>();
		for(SingleInterval region: sampleSnpCounts.keySet()){
			net.sf.samtools.util.CloseableIterator<Gene> iter=genes.sortedIterator(region, true);
			Pair<Integer> sampleVals=sampleSnpCounts.get(region);
			while(iter.hasNext()){
				Gene gene=iter.next();
				Pair<Integer> s=new Pair<Integer>(sampleVals.getValue1(), sampleVals.getValue2());
				if(rtrn.containsKey(gene)){
					s=merge2(rtrn.get(gene), s);
				}
				
				rtrn.put(gene,  s);
			}
			iter.close();
		}
		return rtrn;
	}


	private static Pair<Integer> merge2(Pair<Integer> pair, Pair<Integer> sampleVals) {
		Pair<Integer> rtrn=new Pair<Integer>();
		rtrn.setValue1(pair.getValue1()+sampleVals.getValue1());
		rtrn.setValue2(pair.getValue2()+sampleVals.getValue2());
		
		return rtrn;
	}


	public static void main(String[] args) throws IOException{
		
		
		
		
		
		if(args.length>4){
			File file1=new File(args[0]);
			File file2=new File(args[1]);
			File snpFile=new File(args[2]);
			String save=args[3];
			
			
			Map<String, IntervalTree<String>> snpTree=parseTree(snpFile);
			
			QuantifyBySNP snp2=new QuantifyBySNP(file1, file2, snpTree, save);
			
			
			
		}
		else{System.err.println(usage);}
		
		
	}
	
	

	










	static String usage=" args[0]=bam file1 \n args[1]=bam file2 \n args[2]=snp file \n args[3]=save file name";


	
	
}