package guttmanlab.core.clap.old;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.codehaus.jackson.map.DeserializationConfig.Feature;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import guttmanlab.core.annotationcollection.FeatureCollection;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.simulation.CoordinateSpace;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.util.CloseableIterator;

public class PercentOfReadsByFeature {
	
	static String intergenic="intergenic";
	static String intronic="intron";
	static String exonic="exon";
	static String utr3="utr3";
	static String utr5="utr5";
	static String cds="cds";
	static String amb="ambiguous";
	static String ncRNA="ncRNA";
	static String repeat="repeat";
	static String rnaRepeat="RNARepeat";
	
	
	/*static String first="first";
	static String second="second";
	static String third="third";
	static String last="last";
	static String secondLast="2last";
	static String thirdLast="3last";
	static String none="none";
	static String after="after";*/

	public PercentOfReadsByFeature(File bamFile, Map<String, Collection<Gene>> regions, Map<String, IntervalTree<SingleInterval>> repeats, String save) throws IOException{
		BAMSingleReadCollection bam1=new BAMSingleReadCollection(bamFile); //TODO Get fragments
		Map<String, Integer> quant=quantification(bam1, regions, repeats, save);
		//Map<String, Integer> featureLength=computeFeatureLengths(regions);
		write(save, quant);
	}
	
	public PercentOfReadsByFeature(Collection<SingleInterval> reads, Map<String, Collection<Gene>> regions, String save) throws IOException{
		Map<String, Integer> quant=quantification(reads, regions);
		//Map<String, Integer> featureLength=computeFeatureLengths(regions);
		write(save, quant);
	}
	
	private Map<String, Integer> computeFeatureLengths(Map<String, Collection<Gene>> regions) {
		int intronSize=0;
		int cdsSize=0;
		int utr5Size=0;
		int utr3Size=0;
		int ncRNASize=0;
		
		for(String chr: regions.keySet()){
			for(Gene gene: regions.get(chr)){
				for(Annotation intron: gene.getIntrons()){intronSize+=intron.size();}
				
				if(gene.hasCodingRegion()){
					cdsSize+=gene.getCodingRegion().size();
					if(gene.get5UTR()!=null){utr5Size+=gene.get5UTR().size();}
					if(gene.get3UTR()!=null){utr3Size+=gene.get3UTR().size();}
				}
				else{
					ncRNASize+=gene.size();
				}
			}
		}
		
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		rtrn.put(intronic, intronSize);
		rtrn.put(cds, cdsSize);
		rtrn.put(ncRNA, ncRNASize);
		rtrn.put(utr5, utr5Size);
		rtrn.put(utr3, utr3Size);
		
		return rtrn;
	}

	private Map<String, Integer> quantification(BAMSingleReadCollection reads, Map<String, Collection<Gene>> regions, Map<String, IntervalTree<SingleInterval>> repeatTree, String save){
		Map<String, IntervalTree<Gene>> geneBodyTree=makeGeneTree(regions);
		Map<String, Integer> featureCount=new HashMap<String, Integer>();
		//Map<String, IntervalTree<Gene>> repeatTree=makeGeneTree(repeats);
		
		CloseableIterator<SAMFragment> iter=reads.sortedIterator();
		
		SAMFileHeader header=CoordinateSpace.MM9.getSAMFileHeaderOldVersion();
		header.setSortOrder(SortOrder.coordinate);
		//SAMFileWriter alignmentWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, new File(save+".intergenic.bam"));
		
		
		int counter=0;
		while(iter.hasNext()){
			Set<String> features=new TreeSet<String>();
			SAMFragment read=iter.next();
			//Does read overlap gene?
			String consensus=intergenic;
			
			if(geneBodyTree.containsKey(read.getReferenceName())){
				Iterator<Gene> geneOverlaps=geneBodyTree.get(read.getReferenceName()).overlappingValueIterator(read.getReferenceStartPosition(), read.getReferenceEndPosition());
				
				if(geneOverlaps.hasNext()){
					while(geneOverlaps.hasNext()){
						Gene gene=geneOverlaps.next();
						if(gene.getOrientation().equals(read.getOrientation())){
							String f=overlappingRegion(read, gene);
							features.add(f);
							//System.err.println(gene.getName()+" "+read.getName()+" "+read.getOrientation());
						}
					}
					consensus=getConsensusFeature(features, read);
				}
			}
			
			if(consensus.equalsIgnoreCase(intergenic)){
				if(repeatTree.containsKey(read.getReferenceName())){
					Iterator<SingleInterval> repeatOverlaps=repeatTree.get(read.getReferenceName()).overlappingValueIterator(read.getReferenceStartPosition(), read.getReferenceEndPosition());
					if(repeatOverlaps.hasNext()){
						consensus=ncRNA;	
					}
				}
			}
			updateFeatureCount(featureCount, consensus);
			
			/*if(consensus==intergenic){
				alignmentWriter.addAlignment(read.getSamRecord());
			}*/
			
			counter++;
			if(counter%1000000 ==0){System.err.println("V2 "+counter);}
		}
		
		iter.close();
		System.out.flush();
		//alignmentWriter.close();
		
		return featureCount;
	
	}

	

	private Map<String, Integer> quantification(Collection<SingleInterval> reads, Map<String, Collection<Gene>> regions){
		Map<String, IntervalTree<Gene>> geneBodyTree=makeGeneTree(regions);
		
		System.err.println("Made gene tree");
		
		Map<String, Integer> featureCount=new HashMap<String, Integer>();
		
		Iterator<SingleInterval> iter=reads.iterator();
		
		int counter=0;
		while(iter.hasNext()){
			Set<String> features=new TreeSet<String>();
			Annotation read=iter.next();
			//Does read overlap gene?
			String consensus=intergenic;
			
			if(geneBodyTree.containsKey(read.getReferenceName())){
				Iterator<Gene> geneOverlaps=geneBodyTree.get(read.getReferenceName()).overlappingValueIterator(read.getReferenceStartPosition(), read.getReferenceEndPosition());
				
				if(geneOverlaps.hasNext()){
					while(geneOverlaps.hasNext()){
						Gene gene=geneOverlaps.next();
						if(gene.getOrientation().equals(read.getOrientation())){
							String f=overlappingRegion(read, gene);
							features.add(f);
							//System.err.println(gene.getName()+" "+read.getName()+" "+read.getOrientation());
						}
					}
					consensus=getConsensusFeature(features, read);
				}
			}
			
			/*if(consensus.equalsIgnoreCase(intergenic)){
				if(repeatTree.containsKey(read.getReferenceName())){
					Iterator<Gene> repeatOverlaps=repeatTree.get(read.getReferenceName()).overlappingValueIterator(read.getReferenceStartPosition(), read.getReferenceEndPosition());
					if(repeatOverlaps.hasNext()){
						consensus=repeat;
						while(repeatOverlaps.hasNext()){
							Annotation o=repeatOverlaps.next();
							if(o.getName().contains("RNA") || o.getName().startsWith("U")){consensus=rnaRepeat;}
						}
					}
				}
			}*/
			updateFeatureCount(featureCount, consensus);
			
			counter++;
			//if(counter%10 ==0){System.err.println("V2 "+counter);}
		}
		
		System.out.flush();
		
		
		return featureCount;
	
	}
	
	

	private Map<String, IntervalTree<Annotation>> makeGeneTree(Collection<? extends Annotation> repeats) {
		// TODO Auto-generated method stub
		return null;
	}

	private String overlappingRegion(SAMFragment read, Gene gene) {
		boolean inIntron=gene.inIntron(read.getReferenceStartPosition()) || gene.inIntron(read.getReferenceEndPosition());
		if(inIntron){return intronic;}
		
		//TODO NEED TO USE STRAND INFO
		
		if(gene.hasCodingRegion()){
			//else both within CDS?
			boolean inCDS=gene.getCodingRegion().isWithinBlock(read.getReferenceStartPosition()) && gene.getCodingRegion().isWithinBlock(read.getReferenceEndPosition());
			if(inCDS){return cds;}
			
			//else in 5' UTR?
			Annotation utr5Annotation=gene.get5UTR();
			if(utr5Annotation!=null){
				boolean inUTR5=utr5Annotation.isWithinBlock(read.getReferenceStartPosition())|| utr5Annotation.isWithinBlock(read.getReferenceEndPosition());
				if(inUTR5){return utr5;}
			}
			
			//else in 3' UTR
			Annotation utr3Annotation=gene.get3UTR();
			if(utr3Annotation!=null){
				boolean inUTR3=utr3Annotation.isWithinBlock(read.getReferenceStartPosition())|| utr3Annotation.isWithinBlock(read.getReferenceEndPosition());
				if(inUTR3){return utr3;}
			}
		}
		
		
		
		return ncRNA;
	}
	
	private String overlappingRegion(Annotation read, Gene gene) {
		boolean inIntron=gene.inIntron(read.getReferenceStartPosition()) || gene.inIntron(read.getReferenceEndPosition());
		if(inIntron){return intronic;}
		
		//TODO NEED TO USE STRAND INFO
		
		if(gene.hasCodingRegion()){
			//else both within CDS?
			boolean inCDS=gene.getCodingRegion().isWithinBlock(read.getReferenceStartPosition()) && gene.getCodingRegion().isWithinBlock(read.getReferenceEndPosition());
			if(inCDS){return cds;}
			
			//else in 5' UTR?
			Annotation utr5Annotation=gene.get5UTR();
			if(utr5Annotation!=null){
				boolean inUTR5=utr5Annotation.isWithinBlock(read.getReferenceStartPosition())|| utr5Annotation.isWithinBlock(read.getReferenceEndPosition());
				if(inUTR5){return utr5;}
			}
			
			//else in 3' UTR
			Annotation utr3Annotation=gene.get3UTR();
			if(utr3Annotation!=null){
				boolean inUTR3=utr3Annotation.isWithinBlock(read.getReferenceStartPosition())|| utr3Annotation.isWithinBlock(read.getReferenceEndPosition());
				if(inUTR3){return utr3;}
			}
		}
		
		
		
		return ncRNA;
	}

	private void updateFeatureCount(Map<String, Integer> featureCount, String consensus) {
		int count=0;
		if(featureCount.containsKey(consensus)){count=featureCount.get(consensus);}
		count++;
		featureCount.put(consensus, count);
	}

	private String getConsensusFeature(Set<String> features, SAMFragment read) {
		if(features.isEmpty()){return intergenic;}
		if(features.size()==1){return features.iterator().next();}
		
		//TODO If intron and other, choose other
		if(features.size()==2){
			if(features.contains(intronic)){
				for(String f: features){
					if(!f.equals(intronic)){return f;}
				}
			}
			if(features.contains(ncRNA)){
				for(String f: features){
					if(!f.equals(ncRNA)){return f;}
				}
			}
		}
		
		//else{System.err.println(read.toBED()+" "+features.toString());}
		
		return amb; //TODO Fix this
	}

	private String getConsensusFeature(Set<String> features, Annotation read) {
		if(features.isEmpty()){return intergenic;}
		if(features.size()==1){return features.iterator().next();}
		
		//TODO If intron and other, choose other
		if(features.size()==2){
			if(features.contains(intronic)){
				for(String f: features){
					if(!f.equals(intronic)){return f;}
				}
			}
			if(features.contains(ncRNA)){
				for(String f: features){
					if(!f.equals(ncRNA)){return f;}
				}
			}
		}
		
		//else{System.err.println(read.toBED()+" "+features.toString());}
		
		return amb; //TODO Fix this
	}
	
	private Map<String, IntervalTree<Annotation>> makeExonTree(Map<String, Collection<Annotation>> regions) {
		Map<String, IntervalTree<Annotation>> rtrn=new TreeMap<String, IntervalTree<Annotation>>();
		double geneLength=0;
		double exonLength=0;
		
		for(String chr: regions.keySet()){
			IntervalTree<Annotation> tree=new IntervalTree<Annotation>();
			for(Annotation gene: regions.get(chr)){
				geneLength+=(gene.getReferenceEndPosition()-gene.getReferenceStartPosition());
				for(Annotation exon:gene.getBlockSet()){
					exonLength+=exon.size();
					tree.put(exon.getReferenceStartPosition(), exon.getReferenceEndPosition(), exon);
				}
			}
			rtrn.put(chr, tree);
		}
		
		System.err.println(exonLength+" "+geneLength+" "+(exonLength/geneLength));
		
		return rtrn;
	}
	
	
	private Map<String, IntervalTree<Gene>> makeGeneTree(Map<String, Collection<Gene>> regions) {
		Map<String, IntervalTree<Gene>> rtrn=new TreeMap<String, IntervalTree<Gene>>();
		
		for(String chr: regions.keySet()){
			IntervalTree<Gene> tree=new IntervalTree<Gene>();
			for(Gene gene: regions.get(chr)){
				tree.put(gene.getReferenceStartPosition(), gene.getReferenceEndPosition(), gene);
			}
			rtrn.put(chr, tree);
		}
		
		return rtrn;
	}
	
	

	private Map<String, Integer> featureCounts(Gene gene, CloseableIterator<SAMFragment> readIter) {
		int ambCount=0;
		int intronCount=0;
		int exonCount=0;
		
		IntervalTree<SingleInterval> exons=new IntervalTree<SingleInterval>();
		IntervalTree<SingleInterval> introns=new IntervalTree<SingleInterval>();
		for(Annotation exon:gene.getExonSet()){
			exons.put(exon.getReferenceStartPosition(), exon.getReferenceEndPosition(), exon.getSingleInterval());
		}
		for(SingleInterval intron: gene.getIntronSet()){
			introns.put(intron.getReferenceStartPosition(), intron.getReferenceEndPosition(), intron);
		}
		
		while(readIter.hasNext()){
			SAMFragment read=readIter.next();
			boolean overlapsExon=exons.hasOverlappers(read.getReferenceStartPosition(), read.getReferenceEndPosition());
			boolean overlapsIntron=introns.hasOverlappers(read.getReferenceStartPosition(), read.getReferenceEndPosition());
			
			if(overlapsExon && overlapsIntron){ambCount++;}
			else if(overlapsExon){exonCount++;}
			else if(overlapsIntron){intronCount++;}
		}
		
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		rtrn.put("ambiguous", ambCount);
		rtrn.put("exon", exonCount);
		rtrn.put("intron", intronCount);
		return rtrn;
	}

	private String getFeature(SAMFragment read, CloseableIterator<Gene> overlappingGenes) {
		IntervalTree<SingleInterval> exons=new IntervalTree<SingleInterval>();
		IntervalTree<SingleInterval> introns=new IntervalTree<SingleInterval>();
		
		while(overlappingGenes.hasNext()){
			Gene gene=overlappingGenes.next();
			for(Annotation exon:gene.getExonSet()){
				exons.put(exon.getReferenceStartPosition(), exon.getReferenceEndPosition(), exon.getSingleInterval());
			}
			for(SingleInterval intron: gene.getIntronSet()){
				introns.put(intron.getReferenceStartPosition(), intron.getReferenceEndPosition(), intron);
			}
		}
		
		boolean overlapsExon=exons.hasOverlappers(read.getReferenceStartPosition(), read.getReferenceEndPosition());
		boolean overlapsIntron=introns.hasOverlappers(read.getReferenceStartPosition(), read.getReferenceEndPosition());
		
		if(overlapsExon && overlapsIntron){return "ambigous";}
		if(overlapsExon){return "exon";}
		if(overlapsIntron){return "intron";}
		
		return "intergenic";
	}

	private void write(String save, Map<String, Integer> featureCount) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String feature: featureCount.keySet()){
			int count=featureCount.get(feature);
			writer.write(feature+"\t"+count+"\n");
		}
		
		writer.close();
	}
	
	/*public static void main (String[] args) throws IOException{
		if(args.length>2){
			Collection<SingleInterval> sigWindows=BEDFileIO.loadCustom(args[0]);
			//File file=new File(args[0]);
			Map<String, Collection<Gene>> regions=BEDFileIO.loadRegionsFromFileByChr(args[1]);
			
			String save=args[2];
			System.err.println("Started "+sigWindows.size());
			new PercentOfReadsByFeature(sigWindows, regions, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=.sig file \n args[1]=Annotation BED \n args[2]=save";*/
	
	public static void main(String[] args) throws IOException{
		if(args.length>2){
			File file=new File(args[0]);
			Map<String, Collection<Gene>> regions=BEDFileIO.loadRegionsFromFileByChr(args[1]);
			//Map<String, Collection<Gene>> regions=new TreeMap();
			Map<String, IntervalTree<SingleInterval>> repeats=BEDFileIO.loadSingleIntervalTreeFromGTF(args[2], false);
					
			String save=args[3];
			
			new PercentOfReadsByFeature(file, regions, repeats, save);
		}
		else{System.err.println(usage);}
	}
	
	

	static String usage=" args[0]=BAM file \n args[1]=Gene BED \n args[2]=repeats (ncRNA gtf) \n args[3]=save";
}
