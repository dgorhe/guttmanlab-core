package guttmanlab.core.orfcloning;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.IntervalTree.Node;
import net.sf.samtools.util.CloseableIterator;

public class ExpressionLevels {
	
	Map<Gene, Double> geneExpression;

	public ExpressionLevels(BAMSingleReadCollection bam1, Collection<Gene> genes) throws IOException{
		this.geneExpression=getExpression(bam1, genes);
		
		//write(save, geneExpression);
		
	}
	
	
	public Map<Gene, Double> getGeneExpression(){
		return this.geneExpression;
	}
	
	/*private void write(String save, Map<Gene, Double> geneExpression) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Gene gene: geneExpression.keySet()){
			if(gene.hasCodingRegion()){
				writer.write(gene.getName()+"\t"+geneExpression.get(gene)+"\n");
			}
		}
		
		writer.close();
	}*/


	private Annotation get3UTRProbe(Gene gene, int probeLength) {
		return gene.get3UTR().trimByRelativePositions(0, probeLength);	
	}


	private Annotation get5UTRProbe(Gene gene, int probeLength) {
		return gene.get5UTR().trimByRelativePositions(gene.get5UTR().size()-probeLength, gene.get5UTR().size());
	}
	
	
	private Map<Gene, Double> getExpression(BAMSingleReadCollection bam1, Collection<Gene> genes) {
		int totalNumReads=0;
		Map<String, IntervalTree<Gene>> geneTree=makeTree(genes);
		
		Map<Gene, Integer> geneCount=new TreeMap<Gene, Integer>();
		CloseableIterator<SAMFragment> iter=bam1.sortedIterator();
		
		int counter=0;
		while(iter.hasNext()){
			SAMFragment fragment=iter.next();
			totalNumReads++;
			Collection<Gene> overlapping=getGenes(fragment, geneTree);
			for(Gene gene: overlapping){
				int count=0;
				if(geneCount.containsKey(gene)){
					count=geneCount.get(gene);
				}
				count++;
				geneCount.put(gene, count);
			}
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		return normalize(geneCount);
		
	}


	private Map<String, IntervalTree<Gene>> makeTree(Collection<Gene> genes) {
		Map<String, IntervalTree<Gene>> rtrn=new TreeMap<String, IntervalTree<Gene>>();
		
		for(Gene gene: genes){
			if(gene.hasCodingRegion()){
				IntervalTree<Gene> tree=new IntervalTree<Gene>();
				if(rtrn.containsKey(gene.getReferenceName())){
					tree=rtrn.get(gene.getReferenceName());
				}
				tree.put(gene.getReferenceStartPosition(), gene.getReferenceEndPosition(), gene);
				rtrn.put(gene.getReferenceName(), tree);
			}
		}
		
		return rtrn;
	}


	private Collection<Gene> getGenes(SAMFragment fragment, Map<String, IntervalTree<Gene>> geneTree) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		
		if(geneTree.containsKey(fragment.getReferenceName())){
			IntervalTree<Gene> tree=geneTree.get(fragment.getReferenceName());
			
			//Iterator<Node<Gene>> overlappers=tree.overlappers(fragment.getReferenceStartPosition(), fragment.getReferenceEndPosition());
			Iterator<Gene> overlappers=tree.overlappingValueIterator(fragment.getReferenceStartPosition(), fragment.getReferenceEndPosition());
			while(overlappers.hasNext()){
				Gene gene=overlappers.next();
				if(fragment.overlaps(gene)){rtrn.add(gene);}
			}
		
		}
		
		return rtrn;
	}


	private Map<Gene, Double> normalize(Map<Gene, Integer> geneCount) {
		Map<Gene, Double> rtrn=new TreeMap<Gene, Double>();
		
		double sumTotal=getTotal(geneCount);
		
		for(Gene gene: geneCount.keySet()){
			double count=geneCount.get(gene);
			double length=gene.size();
			double ratio=count/length;
			
			//Divide the read counts by the length of each gene in kilobases
			//Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
			//Divide the RPK values by the “per million” scaling factor. This gives you TPM.
			double tpm=Math.pow(10, 6)*(ratio/sumTotal);
			
			rtrn.put(gene, tpm);
		}
		
		return rtrn;
	}


	private double getTotal(Map<Gene, Integer> geneCount) {
		double sumTotal=0;
		
		for(Gene gene: geneCount.keySet()){
			double count=geneCount.get(gene);
			double length=gene.size();
			double ratio=count/length;
			sumTotal+=ratio;
		}
		
		return sumTotal;
	}


	


	private static void write(String save, Map<String, Map<Gene, Double>> temp) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<Gene> allGenes=getAllGenes(temp);
		
		writer.write("geneName");
		for(String name: temp.keySet()){
			writer.write("\t"+name);
		}
		writer.write("\n");
		
		for(Gene gene: allGenes){
			writer.write(gene.getName());
			for(String name: temp.keySet()){
				double value=0.0;
				if(temp.get(name).containsKey(gene)){
					value=temp.get(name).get(gene);
				}
				writer.write("\t"+value);
			}
			writer.write("\n");
		}
		
		writer.close();
	}


	private static Collection<Gene> getAllGenes(Map<String, Map<Gene, Double>> temp) {
		Collection<Gene> allGenes=new TreeSet<Gene>();
		
		for(String name: temp.keySet()){
			allGenes.addAll(temp.get(name).keySet());
		}
		
		return allGenes;
	}
	
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
			/*String dir=args[0];
			Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[1]);
			String save=args[2];
			
			Map<String, Map<Gene, Double>> temp=new TreeMap<String, Map<Gene, Double>>();
			
			File[] files=new File(dir).listFiles();
			for(int i=0; i<files.length; i++){
				if(files[i].isDirectory()){
					String name=files[i].getName();
					File bamFile=new File(files[i].getAbsolutePath()+"/aligned.MAPQ255STAR.combo.bam");
					System.err.println(bamFile.getAbsolutePath()+" "+name);
					BAMSingleReadCollection bam=new BAMSingleReadCollection(bamFile);
					ExpressionLevels expression=new ExpressionLevels(bam, genes);
					Map<Gene, Double> geneExpression=expression.getGeneExpression();
					temp.put(name, geneExpression);
				}
			}
			
			write(save, temp);*/
			
			System.err.println("updated");
			File[] bamFiles=new File(args[0]).listFiles();
			Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[1]);
			String save=args[2];
			Map<String, Map<Gene, Double>> temp=new TreeMap<String, Map<Gene, Double>>();
			
			
			for(int i=0; i<bamFiles.length; i++){
				String name=bamFiles[i].getName();
				System.err.println(name);
				BAMSingleReadCollection bam=new BAMSingleReadCollection(bamFiles[i]);
				ExpressionLevels expression=new ExpressionLevels(bam, genes);
				Map<Gene, Double> geneExpression=expression.getGeneExpression();
				
				temp.put(name, geneExpression);
			}
			write(save, temp);
			
		}
		
		else{
			System.err.println(usage);
		}
		
	}
	
	static String usage=" args[0]=bamFiles \n args[1]=genes \n args[2]=save";
	
}
