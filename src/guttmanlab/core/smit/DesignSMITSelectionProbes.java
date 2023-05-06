package guttmanlab.core.smit;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.smit.AssignReadsToSplicingStates.SpliceSite;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.BlockedAnnotation;

public class DesignSMITSelectionProbes {

	int probeLength=120;
	
	public DesignSMITSelectionProbes(Collection<Gene> genes, String save) throws IOException{
		Collection<Gene> regions=new TreeSet<Gene>();
		//probe directly at 3'SS
		for(Gene gene: genes){
			Collection<Gene> unsplicedProbes=getUnsplicedProbesCenteredOn3SS(gene);
			Collection<Gene> splicedProbes=getSplicedProbesCenteredOn3SS(gene);
			regions.addAll(unsplicedProbes);
			regions.addAll(splicedProbes);
		}
		write(save, regions);
	}
	
	private Collection<Gene> getUnsplicedProbesCenteredOn3SS(Gene gene) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		List<Integer> spliceSites=get3PrimeSS(gene);
		for(Integer ss: spliceSites){
			int extension=probeLength/2;
			//Make unspliced probe
			SingleInterval probe=new SingleInterval(gene.getReferenceName(), ss-extension, ss+extension, gene.getOrientation());
			Gene g=new Gene(probe);
			g.setName(gene.getName());
			rtrn.add(g);
		}
		return rtrn;
	}
	
	private Collection<Gene> getSplicedProbesCenteredOn3SS(Gene gene){
		Collection<Gene> rtrn=new TreeSet<Gene>();
		
		Collection<Annotation> pairs=gene.getExonIntronPairs();
		
		for(Annotation junction: pairs){
			BlockedAnnotation a=new BlockedAnnotation();
			
			//SingleInterval exon5 =get5PrimeExon(junction);
			SingleInterval exon5 = get5PrimeExon(junction);
			Annotation intron=getIntron(junction);
			
			int ss=intron.get5PrimePosition();
			int extension=probeLength/2;
			
			SingleInterval probe1=new SingleInterval(gene.getReferenceName(), ss, ss+extension, gene.getOrientation());
			if(gene.getOrientation().equals(Strand.NEGATIVE)){
				probe1=new SingleInterval(gene.getReferenceName(), ss-extension, ss, gene.getOrientation());
			}
			
			SingleInterval probe2=new SingleInterval(gene.getReferenceName(), exon5.get3PrimePosition()-probeLength, exon5.get3PrimePosition(), gene.getOrientation());
			if(gene.getOrientation().equals(Strand.NEGATIVE)){
				probe2=new SingleInterval(gene.getReferenceName(), exon5.get3PrimePosition(), exon5.get3PrimePosition()+probeLength, gene.getOrientation());
			}
			
			a.addBlocks(probe1);
			a.addBlocks(probe2);
			a.setOrientation(gene.getOrientation());
			a.setName(gene.getName());
			rtrn.add(new Gene(a));
		}
		
		return rtrn;
	}

	private Annotation getIntron(Annotation junction) {
		Annotation intron=junction.getIntrons().iterator().next();
		return intron;
	}

	
	
	private SingleInterval get5PrimeExon(Annotation junction) {
		Iterator<SingleInterval> exons=junction.getBlocks();
		SingleInterval exon1=exons.next();
		SingleInterval exon2=exons.next();
		
		if(exon2.getReferenceStartPosition()>exon1.getReferenceStartPosition()){
			if(junction.getOrientation().equals(Strand.NEGATIVE)){return exon2;}
			else{return exon1;}
		}
		
		else{
			if(junction.getOrientation().equals(Strand.NEGATIVE)){return exon1;}
			else{return exon2;}
		}
		
		
	}

	private void write(String save, Collection<Gene> genes) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Gene gene: genes){
			writer.write(gene.toBED()+"\n");
		}
		
		writer.close();
	}
	
	

	private List<Integer> get3PrimeSS(Gene gene){
		List<Integer> rtrn=new ArrayList<Integer>();
		
		if(gene.getNumberOfBlocks()>1){
			SingleInterval firstExon=getFirstExon(gene);
			//iterate through blocks and get 3'SS
			Iterator<SingleInterval> exons=gene.getBlocks();
			//Test if fragment contains it
			while(exons.hasNext()){
				SingleInterval exon=exons.next();
				if(exon!=firstExon){
					int ss=exon.get5PrimePosition();
					rtrn.add(ss);
				}
			}
		}
		return rtrn;
	}
	
	
	
	private SingleInterval get5PrimeRegion(SingleInterval exon) {
		SingleInterval region=new SingleInterval(exon.getReferenceName(), exon.getReferenceStartPosition(), exon.getReferenceStartPosition()+probeLength, exon.getOrientation());
		if(exon.getOrientation().equals(Strand.NEGATIVE)){
			region=new SingleInterval(exon.getReferenceName(),exon.getReferenceEndPosition()-probeLength, exon.getReferenceEndPosition(), exon.getOrientation());
		}
		return region;
	}

	private SingleInterval getFirstExon(Gene gene) {
		SingleInterval firstExon=gene.getFirstBlock();
		if(gene.getOrientation().equals(Strand.NEGATIVE)){
			firstExon=gene.getLastBlock();
		}
		return firstExon;
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>1){
			Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[0]);
			String save=args[1];
			new DesignSMITSelectionProbes(genes, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=BED \n args[1]=save";
	
}
