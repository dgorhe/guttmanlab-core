package guttmanlab.core.test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class CountReadsPerGene {

	public CountReadsPerGene(File bamFile, Map<String, IntervalTree<Gene>> geneTree, String save) throws IOException{
		SAMFileReader reader=new SAMFileReader(bamFile);
		SAMRecordIterator reads=reader.iterator();
		
		Map<Gene, Integer> geneCounts=new TreeMap<Gene, Integer>();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			Collection<Gene> overlappingGenes=getOverlappingGenes(record, geneTree);
			for(Gene gene: overlappingGenes){
				int count=0;
				if(geneCounts.containsKey(gene)){count=geneCounts.get(gene);}
				count++;
				geneCounts.put(gene, count);
			}
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		reads.close();
		reader.close();
		
		write(save, geneCounts);
	}
	
	private void write(String save, Map<Gene, Integer> geneCounts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Gene gene: geneCounts.keySet()){
			writer.write(gene.getName()+"\t"+gene.toUCSC()+"\t"+geneCounts.get(gene)+"\n");
		}
		
		writer.close();
	}

	private Collection<Gene> getOverlappingGenes(SAMRecord record, Map<String, IntervalTree<Gene>> geneTree) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		
		if(!geneTree.containsKey(record.getReferenceName())){return rtrn;}
		Iterator<Gene> iter=geneTree.get(record.getReferenceName()).overlappingValueIterator(record.getAlignmentStart(), record.getAlignmentEnd());
		while(iter.hasNext()){
			rtrn.add(iter.next());
		}
		
		return rtrn;
	}

	public static void main(String[] args) throws IOException{
		File bamFile=new File(args[0]);
		Map<String, IntervalTree<Gene>> tree=BEDFileIO.loadTree(args[1]);
		String save=args[2];
		new CountReadsPerGene(bamFile, tree, save);
	}
	
}
