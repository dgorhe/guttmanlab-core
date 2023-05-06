package guttmanlab.core.orfcloning;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;

//Probes that anneal to 3' and 5' UTRs to cleave by RNaseH
public class DesignORFCaptureProbes {
	
	int probeLength=90;
	Collection<Annotation> overlappers;

	public DesignORFCaptureProbes(Map<String, Collection<Gene>> genesByChr, String save) throws IOException{
		Collection<Annotation> utr5Probes=new TreeSet<Annotation>();
		Collection<Annotation> utr3Probes=new TreeSet<Annotation>();
		this.overlappers=new TreeSet<Annotation>();
		
		for(String chr: genesByChr.keySet()){
			System.err.println(chr);
			Collection<Gene> genes=genesByChr.get(chr);
			utr5Probes.addAll(get5UTRProbes(genes));
			utr3Probes.addAll(get3UTRProbes(genes));
		}
		
		write(save+".5Utr.bed", utr5Probes);
		write(save+".3Utr.bed", utr3Probes);
		
		//TODO For overlappers we need to also resynthesize orthogonal probe to make sure selection works
		write(save+".overlappers.bed", overlappers);
		
	}

	
	private Collection<Annotation> get3UTRProbes(Collection<Gene> genes) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		for(Gene gene: genes){
			if(gene.hasCodingRegion()){
				if(gene.get3UTR()!=null){
					Annotation utr3Probe=gene.get3UTR().trimByRelativePositions(0, probeLength);
					boolean overlaps=overlapsCDS(utr3Probe, genes);
					if(!overlaps){rtrn.add(utr3Probe);}
					else{overlappers.add(utr3Probe);}
				}
			}
		}
		return rtrn;
	}


	private Collection<Annotation> get5UTRProbes(Collection<Gene> genes) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		
		for(Gene gene: genes){
			if(gene.hasCodingRegion()){
				if(gene.get5UTR()!=null){ 
					Annotation utr5Probe=gene.get5UTR().trimByRelativePositions(gene.get5UTR().size()-probeLength, gene.get5UTR().size());
					boolean overlaps=overlapsCDS(utr5Probe, genes);
					if(!overlaps){rtrn.add(utr5Probe);}
					else{overlappers.add(utr5Probe);}
				}
			}
		}
		
		return rtrn;
	}



	private void write(String save, Collection<Annotation> probes) throws IOException {
		FileWriter writer=new FileWriter(save);
		for(Annotation probe: probes){
			writer.write(probe.toBED()+"\n");
		}
		writer.close();
	}


	private boolean overlapsCDS(Annotation probe, Collection<Gene> genes) {
		for(Gene gene: genes){
			boolean overlaps=gene.getCodingRegion().overlaps(probe);
			if(overlaps){return true;}
		}
		return false;
	}
	

	private void overlapsCDS(Collection<Annotation> probes, Collection<Gene> genes) {
		for(Annotation probe: probes){
			for(Gene gene: genes){
				boolean overlaps=gene.getCodingRegion().overlaps(probe);
				if(overlaps){System.out.println(probe.toBED());}
			}
		}
		
	}



	public static void main(String[] args) throws IOException{
		Map<String, Collection<Gene>> genesByChr=BEDFileIO.loadRegionsFromFileByChr(args[0]);
		//Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[0]);
		String save=args[1];
		new DesignORFCaptureProbes(genesByChr, save);
	}
	
}
