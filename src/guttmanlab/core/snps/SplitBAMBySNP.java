package guttmanlab.core.snps;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.datastructures.IntervalTree.Node;
import guttmanlab.core.math.Statistics;
import guttmanlab.core.datastructures.Pair;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;


public class SplitBAMBySNP {

	static String ambiguous="amb";
	int minCount=5;
	
	
	
	public SplitBAMBySNP(File bamFile, File snpFile, String save) throws IOException{
		
		
		SAMFileReader reader=new SAMFileReader(bamFile);
		Map<String, IntervalTree<String>> snpTree=parseTree(snpFile);
		
		SAMFileHeader fileHeader=reader.getFileHeader();
		SAMFileWriter writer1=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(fileHeader, false, new File(save+".SNP1.bam"));
		SAMFileWriter writer2=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(fileHeader, false, new File(save+".SNP2.bam"));
		SAMFileWriter writer3=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(fileHeader, false, new File(save+".amb.bam"));
		
		
		
		
		List<String> rows=getRows(snpTree);
		List<String> columns=new ArrayList<String>();
		columns.add("snp1");
		columns.add("snp2");
		MatrixWithHeaders data=new MatrixWithHeaders(rows, columns);
		
		
		
		SAMRecordIterator reads=reader.iterator();
		
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			
			String genotype=getOverlappingRegion(record, snpTree);
			
			String row=genotype.split("_")[0];
			if(genotype.endsWith("SNP1")){
				data.incrementCount(row, "snp1");
				writer1.addAlignment(record);
			}
			else if(genotype.endsWith("SNP2")){
				data.incrementCount(row, "snp2");
				writer2.addAlignment(record);
			}
			else {writer3.addAlignment(record);}
			
			
			
			
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		reads.close();
		reader.close();
		data=filter(data, minCount);
		data.write(save);
		
		writer1.close();
		writer2.close();
		writer3.close();
	}
	
	
	private MatrixWithHeaders filter(MatrixWithHeaders data, int minCount2) {
		List<String> newRows=new ArrayList<String>();
		
		for(String row: data.getRowNames()) {
			double[] vals=data.getRow(row);
			if(Statistics.sum(vals)>minCount2) {newRows.add(row);}
		}
		
		return data.submatrixByRowNames(newRows);
	}


	private List<String> getRows(Map<String, IntervalTree<String>> snpTree) {
		List<String> rtrn=new ArrayList<String>();
		
		for(String chr: snpTree.keySet()) {
			Iterator<Node<String>> iter=snpTree.get(chr).iterator();
			while(iter.hasNext()) {
				Node<String> n=iter.next();
				String name=chr+":"+n.getStart()+"-"+n.getEnd();
				rtrn.add(name);
			}
		}
		
		return rtrn;
	}
	
	private boolean overlaps(SAMRecord record, Map<String, IntervalTree<String>> snpTree) {
		boolean overlap=false;
		if(snpTree.containsKey(record.getReferenceName())) {
			overlap=snpTree.get(record.getReferenceName()).hasOverlappers(record.getAlignmentStart(), record.getAlignmentEnd());
		}
		return overlap;
	}


	public SplitBAMBySNP(File bamFile, File snpFile, Map<String, IntervalTree<Gene>> genes, String fileName) throws IOException{
		
		
		SAMFileReader reader=new SAMFileReader(bamFile);
		Map<String, IntervalTree<String>> snpTree=parseTree(snpFile);
		
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
			
			if(overlapsGene(record, genes)){
			
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
	}
	
	
	private boolean overlapsGene(SAMRecord record, Map<String, IntervalTree<Gene>> genes) {
		if(genes.containsKey(record.getReferenceName())){
			return genes.get(record.getReferenceName()).hasOverlappers(record.getAlignmentStart(), record.getAlignmentEnd());
		}
		return false;
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


	
	
	private Map<String, IntervalTree<String>> parseTree(File snpFile) throws NumberFormatException, IOException {
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
			if(base.equals(g1)){rtrn.add(chr+":"+g.getStart()+"-"+g.getEnd()+"_SNP1");}
			if(base.equals(g2)){rtrn.add(chr+":"+g.getStart()+"-"+g.getEnd()+"_SNP2");}
		}
		
		
		if(rtrn.size()==1){return rtrn.iterator().next();}
		//else if(rtrn.size()>1){System.err.println(rtrn.size());}
		return ambiguous;
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
	
	
	public static void main(String[] args) throws IOException{
		
		
		if(args.length>2){
			File file=new File(args[0]);
			File snp=new File(args[1]);
			String save=args[2];
			new SplitBAMBySNP(file, snp, save);
			
		}
		else{System.err.println(usage);}
		
		
	}
	
	static String usage=" args[0]=bam file \n args[1]=snp file \n args[2]=save file name";


	
	
}
