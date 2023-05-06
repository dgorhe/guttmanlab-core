package guttmanlab.core.splicing.speckle;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.IntervalTree;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class ComputeGeneExpression {

	private double totalReads;
	private Map<String, Double> rpkm;
	private Map<SingleInterval, Integer> scores;
	
	public ComputeGeneExpression(File bam, GTFToJunctions gtf) throws IOException {
		Map<String, IntervalTree<SingleInterval>> trees=gtf.getGeneCoordinateTree();
		scores=score(bam, trees);
		this.rpkm=computeRPKM(scores);
	}
	
	public Map<String, Double> getRPKM(){return this.rpkm;}
	

	private Map<String, Double> computeRPKM(Map<SingleInterval, Integer> scores) {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		for(SingleInterval g: scores.keySet()) {
			double score=scores.get(g);
			double gLength=(double)g.getGenomicLength()/1000.0;
			double total=totalReads/1000000.0;
			double rpkm=score/(gLength*total);
			//numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
			rtrn.put(g.getName(), rpkm);
			
		}
		
		return rtrn;
	}

	private Map<SingleInterval, Integer> score(File bam, Map<String, IntervalTree<SingleInterval>> trees) {
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		SAMFileReader reader=new SAMFileReader(bam);
		SAMRecordIterator reads=reader.iterator();
		
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			if(!record.getDuplicateReadFlag()) {
			
				SAMFragment read=new SAMFragment(record);
				
				Collection<SingleInterval> genes=findGene(read, trees);
				
				//Collection<Gene> overlappingGenes=overlapsGene(read, genes);
				
				for(SingleInterval g: genes) {
					int score=1;
					if(rtrn.containsKey(g)) {
						score+=rtrn.get(g);
					}
					rtrn.put(g, score);
				}
				counter++;
			}
				
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		this.totalReads=counter;
		reader.close();
		reads.close();
		return rtrn;
	}

	private Collection<SingleInterval> findGene(Annotation fragment, Map<String, IntervalTree<SingleInterval>> geneTree) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		IntervalTree<SingleInterval> tree=geneTree.get(fragment.getReferenceName());
		if(tree!=null){
			Iterator<SingleInterval> genes=tree.overlappingValueIterator(fragment.getReferenceStartPosition(), fragment.getReferenceEndPosition());
			while(genes.hasNext()){
				SingleInterval gene=genes.next();
				if(fragment.getOrientation().equals(gene.getOrientation())){
					rtrn.add(gene);
				}
			}
		}
		return rtrn;
		
	}
	
	public static void main(String[] args) throws IOException {
		File bam=new File(args[0]);
		GTFToJunctions gtf=new GTFToJunctions(new File(args[1]));
		new ComputeGeneExpression(bam, gtf);
	}

	public Map<String, Double> getRPKM(double rpkmFilter) {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		for(String g: this.rpkm.keySet()) {
			double score=this.rpkm.get(g);
			if(score>=rpkmFilter) {rtrn.put(g, score);}
		}
		
		return rtrn;
	}

	public void writeRPKM(String string) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		for(SingleInterval g: scores.keySet()) {
			double score=scores.get(g);
			double gLength=(double)g.getGenomicLength()/1000.0;
			double total=totalReads/1000000.0;
			double rpkm=score/(gLength*total);
			//numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
			writer.write(g.toShortBED()+"\t"+g.getName()+"\t"+rpkm+"\n");
			
		}
		
		writer.close();
	}
	
}
