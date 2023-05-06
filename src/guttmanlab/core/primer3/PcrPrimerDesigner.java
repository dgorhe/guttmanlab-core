package guttmanlab.core.primer3;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import org.apache.log4j.Logger;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.IntervalTree.Node;
import guttmanlab.core.primer3.Primer3SequenceInputTags.SequenceRegionCoordinates;
import guttmanlab.core.sequence.Sequence;



public class PcrPrimerDesigner  {
	public static int n=1;
	public static final int MIN_PRIMER_DES_SPACE = 150; 
	public static final int MIN_PROD_SIZE = 100;
	public static final String USAGE = "Usage: PrimerDesigner TASK=<task_num> <task_args>\n" +
	"\tTasks:\n" +
	"\t\t1. prime gaps: OUTDIR=output directory IN=<sequence AGP> file MAXPRODSIZE=maximum product size BUFFER=size of buffer between primers and target\n" +
	"\t\t2. Generate short amplicon set. This is useful to say desing qPCR exeperiments \n"+ 
		"\t\t\tIN=<sequence AGP>\n"+
		"\t\t\tOUTDIR=<output directory>\n"+
		"\t\t\tDIST=<Distance between amplicons> OR NUM<number of amplicons>\n"+
		"\t\t\tAMPSIZE=<Desired amplicon size>\n"+
		"\t\t\tTARGETSTART=<If supplied, the start of region to cover>\n"+
		"\t\t\tTARGETEND=<If supplied, the end of region to cover>\n" +
		"\t\t\tREPEATOUT=<RepeatMasker output file for the sequence>\n" +
	"\t\t3. Generate gene set amplification set -in <gene list, a gene symbol per line> -minsize <minimum amplicon size> -cdsonly <include this flag if only interested in cds only amplification>\n" +
		"\n\t\t-maxsize <max amplicon size> [-numOfDesigns <number of different designs to ouput (default 1)> -buffer <intronic bases to include (default 0)>] " +
	"\n\t\t4. Generate flanking primers -in <Multifasta file default is standard in> -outdir <directory where to write output files> -optimalDistFromEnds <Optimal distance from sequence begining and end" +
	"\n\t\t\t -maxDistFromEnds <Maximum distance from ends, used when no primer could be found within optimal distance> [-outprefix <a prefix to prepend to output, if not specified no prefix will be prepended to file names>]";
	
	private ArrayList<SequenceRegionCoordinates> exludedAnnotations = new ArrayList<SequenceRegionCoordinates>();
	private static Logger logger = Logger.getLogger(PcrPrimerDesigner.class.getName());

	public PcrPrimerDesigner() {}
	
	
	
	

	/*public static class PrimedRegion {
		ArrayList<PrimerPair> primers = new ArrayList<PrimerPair>();
		GenomicAnnotation targetAnnotation;
		
		public PrimedRegion(GenomicAnnotation annotation) {
			this.targetAnnotation = annotation;
		}
		
		public Collection<? extends Sequence> getPrimerSequences() {
			ArrayList<Sequence> primerSequences = new ArrayList<Sequence>(); 
			for(int i = 0; i < primers.size(); i++) {
				primerSequences.addAll(primers.get(i).getPrimerSequences());
			}
			return primerSequences;
		}

		public void addPrimerPair(PrimerPair pp) {
			primers.add(pp);
		}
		
		public PrimerPair getLastPrimerPair() {
			return primers.size() > 0 ? primers.get(primers.size() - 1) : null;
		}
		
		public List<PrimerPair> getPrimers() { return primers; }
		
		public void write(BufferedWriter goodPrimerWriter,  BufferedWriter goodPrimerConfig, Primer3IO p3io, boolean printExpectedProduct) throws IOException {
			Iterator<PrimerPair> it = primers.iterator();
			PrimerPair pp = null;
			while(it.hasNext()) {
				pp = it.next();
				logger.info(pp);
				if(pp.hasPrimers()) {
					goodPrimerWriter.write(pp.toString(printExpectedProduct));
					goodPrimerWriter.newLine();
					pp.writePrimer3Record(goodPrimerConfig, p3io);
				} 
			}
		}
		
		public void write(BufferedWriter goodPrimerWriter,  BufferedWriter goodPrimerConfig, Primer3IO p3io) throws IOException {
			write(goodPrimerWriter, goodPrimerConfig, p3io, false);
		}
		
	}*/


	/*public static Collection<PrimerPair> designIntronPrimers(Sequence chr, Gene gene, boolean repeatMask, String sequencePrimer, String sequencePrimerRevComp, int numDesigns, int min3, int min5, String pathPrimer3core) throws Exception {
		Primer3Configuration best = Primer3ConfigurationFactory.getJenRTPCRConfiguration();
		
		if(numDesigns>0){best.setPrimerNumReturn(numDesigns);}
		if(min3>0){best.setPrimerMin3PrimeOverlapOfJunction(min3);}
		if(min5>0){best.setPrimerMin5PrimeOverlapOfJunction(min5);}
		
		Primer3IO p3io = new Primer3IO(pathPrimer3core);
		p3io.startPrimer3Communications();
				
		Sequence mRNA=new Sequence(gene.getName());
		mRNA.setSequenceBases(gene.getSequence(chr, repeatMask, false));
		ArrayList<Integer> splicePositions=gene.getSpliceJunctionCoordinates();
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(gene.getName());
		p3sit.addJunctions(splicePositions);
		if(sequencePrimer!=null && !sequencePrimer.isEmpty()){p3sit.setPrimerLeftInput(sequencePrimer);}
		if(sequencePrimerRevComp!=null && !sequencePrimerRevComp.isEmpty()){p3sit.setPrimerRightInput(sequencePrimerRevComp);}
		
		Collection<PrimerPair> pp = p3io.findPrimerPair(p3sit, best);
		

		p3io.endPrimer3Communications();
		
		TreeSet<PrimerPair> rtrn=new TreeSet<PrimerPair>();
		for(PrimerPair pair: pp){
			if(pair.getLeftPrimer()!=null && pair.getRightPrimer()!=null){rtrn.add(pair);}
		}
		
		return rtrn;
	}*/
	
	
	/*public static Collection<PrimerPair> designRACEPrimers(Sequence mRNA, String leftPrimer, String rightPrimer, String pathPrimer3core) throws Exception {
		Primer3Configuration best = Primer3ConfigurationFactory.getJenRTPCRConfiguration();
				
		Primer3IO p3io = new Primer3IO(pathPrimer3core);
		p3io.startPrimer3Communications();
			
		//ArrayList<Integer> splicePositions=gene.getSpliceJunctionCoordinates();
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(mRNA.getId());
		//p3sit.addJunctions(splicePositions);
		if(leftPrimer!=null && !leftPrimer.isEmpty()){p3sit.setPrimerLeftInput(leftPrimer);}
		if(rightPrimer!=null && !rightPrimer.isEmpty()){p3sit.setPrimerRightInput(rightPrimer);}
		
		Collection<PrimerPair> pp = p3io.findPrimerPair(p3sit, best);
		

		p3io.endPrimer3Communications();
		
		
		
		return pp;
	}*/

	/*public static Collection<PrimerPair> designPCRPrimers(Sequence chr, Gene gene, boolean repeatMask,	int numDesigns, boolean crossJunction, String pathPrimer3core) throws Exception {
		Primer3Configuration best = Primer3ConfigurationFactory.getJenRTPCRConfiguration();
		
		if(numDesigns>0){best.setPrimerNumReturn(numDesigns);}
		if(!crossJunction){repeatMask=true;}		
		
		Primer3IO p3io = new Primer3IO(pathPrimer3core);
		p3io.startPrimer3Communications();
			
		
				
		Sequence mRNA=new Sequence(gene.getName());
		mRNA.setSequenceBases(gene.getSequence(chr, repeatMask, false));
		ArrayList<Integer> splicePositions=gene.getSpliceJunctionCoordinates();
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(gene.getName());
		
		if(crossJunction){
			ArrayList<SequenceRegionCoordinates> targetList=new ArrayList<SequenceRegionCoordinates>();
			for(Integer spliceJunction: splicePositions){
				SequenceRegionCoordinates region=new SequenceRegionCoordinates(spliceJunction.intValue()-1, spliceJunction.intValue()+1);
				targetList.add(region);
			}
			p3sit.addTargets(targetList);
		}
		
		
		Collection<PrimerPair> pp = p3io.findPrimerPair(p3sit, best);
		

		p3io.endPrimer3Communications();
		
		TreeSet<PrimerPair> rtrn=new TreeSet<PrimerPair>();
		for(PrimerPair pair: pp){
			if(pair.getLeftPrimer()!=null && pair.getRightPrimer()!=null){rtrn.add(pair);}
		}
		
		return rtrn;
	}*/
	
	/*public static Collection<PrimerPair> designPCRPrimers(Sequence chr, Gene gene, boolean repeatMask,	int numDesigns,  SequenceRegionCoordinates target, String pathPrimer3core) throws Exception {
		Primer3Configuration best = Primer3ConfigurationFactory.getJenRTPCRConfiguration();
		
		if(numDesigns>0){best.setPrimerNumReturn(numDesigns);}	
		
		Primer3IO p3io = new Primer3IO(pathPrimer3core);
		p3io.startPrimer3Communications();		
				
		Sequence mRNA=new Sequence(gene.getName());
		mRNA.setSequenceBases(gene.getSequence(chr, repeatMask, false));
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(gene.getName());
		p3sit.addTarget(target);
		
		
		Collection<PrimerPair> pp = p3io.findPrimerPair(p3sit, best);
		

		p3io.endPrimer3Communications();
		
		TreeSet<PrimerPair> rtrn=new TreeSet<PrimerPair>();
		for(PrimerPair pair: pp){
			if(pair.getLeftPrimer()!=null && pair.getRightPrimer()!=null){rtrn.add(pair);}
		}
		
		return rtrn;
	}*/
	
	
	public static Collection<PrimerPair> designTilingPrimers(Sequence geneSequence, String pathPrimer3core) throws Exception {
		Primer3Configuration best=Primer3ConfigurationFactory.getOptimalPCRConfiguration();
		
		
		Primer3IO p3io = new Primer3IO(pathPrimer3core);
		p3io.startPrimer3Communications();
			
				
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(geneSequence);	
		p3sit.setPrimerSequenceId(geneSequence.getName());
		
		Collection<PrimerPair> primers=new TreeSet<PrimerPair>();
		Collection<PrimerPair> pp = p3io.findPrimerPair(p3sit, best);
			for(PrimerPair primer: pp){
			if(primer!=null && primer.getLeftPrimer()!=null && primer.getRightPrimer()!=null){
				//p3sit.addExcludedRegion(new SequenceRegionCoordinates(primer.getLeftPrimerPosition(), primer.getLeftPrimerPosition()+primer.getLeftPrimer().toCharArray().length));
				//p3sit.addExcludedRegion(new SequenceRegionCoordinates(primer.getRightPrimerPosition()-primer.getRightPrimer().toCharArray().length, primer.getRightPrimerPosition()));
				
				primers.add(primer);
			}
			}
			
		
		
		p3io.endPrimer3Communications();
		
		Collection<PrimerPair> rtrn=new TreeSet<PrimerPair>();
		for(PrimerPair pair: primers){
			if(pair.getLeftPrimer()!=null && pair.getRightPrimer()!=null){rtrn.add(pair);}
		}
		
		
		rtrn=chooseOptimalLayout(rtrn);
		
		
		return rtrn;
	}
	
	private static Collection<PrimerPair> chooseOptimalLayout(Collection<PrimerPair> list) {
		Collection<PrimerPair> rtrn=new ArrayList<PrimerPair>();
		
		//Make Interval tree
		IntervalTree<PrimerPair> tree=makeTree(list);
		
		//Start from the beginning
		Node<PrimerPair> first=tree.min();
		rtrn.add(first.getValue());
		
		Node<PrimerPair> current=first;
		while(current!=null && tree.hasNodeAfterInterval(current.getStart(), current.getEnd())){
			//System.err.println(current.getStart()+" "+current.getEnd());
			//then find all overlappers and take latest
			Iterator<Node<PrimerPair>> overlappers=tree.overlappers(current.getStart(), current.getEnd());
			Node<PrimerPair> largest=getLargest(overlappers);
			if(largest.equals(current)){
				Node<PrimerPair> node=tree.getNodeAfterInterval(largest.getStart(), largest.getEnd());
				//System.err.println("Node "+node.getStart()+" "+node.getEnd());
				largest=node;
			}
			rtrn.add(largest.getValue());
			current=largest;
		}
		
		return rtrn;
	}





	private static Node<PrimerPair> getLargest(Iterator<Node<PrimerPair>> overlappers) {
		Node<PrimerPair> rtrn=null;
		
		while(overlappers.hasNext()){
			rtrn=overlappers.next();
		}
		
		return rtrn;
	}





	private static IntervalTree<PrimerPair> makeTree(Collection<PrimerPair> list) {
		IntervalTree<PrimerPair> rtrn=new IntervalTree<PrimerPair>();
		for(PrimerPair pp:list){
			rtrn.put(pp.getLeftPrimerPosition(), pp.getRightPrimerPosition(), pp);
		}
		return rtrn;
	}





	public static Collection<PrimerPair> designCloningPrimers(Sequence geneSequence, int numDesigns, SequenceRegionCoordinates target, String pathPrimer3core, Primer3Configuration best) throws Exception {
		//if(numDesigns>0){best.setPrimerNumReturn(numDesigns);}
		
		
		Primer3IO p3io = new Primer3IO(pathPrimer3core);
		p3io.startPrimer3Communications();
			
				
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(geneSequence);	
		p3sit.setPrimerSequenceId(geneSequence.getName());
		if(target!=null){p3sit.addTarget(target);}
		//Collection<PrimerPair> pp = p3io.findPrimerPair(p3sit, best);
		
		Collection<PrimerPair> primers=new TreeSet<PrimerPair>();
		for(int i=0; i<numDesigns; i++){
			Collection<PrimerPair> pp = p3io.findPrimerPair(p3sit, best);
			if(pp.iterator().hasNext()){
			PrimerPair primer=pp.iterator().next();
			if(primer!=null && primer.getLeftPrimer()!=null && primer.getRightPrimer()!=null){
				p3sit.addExcludedRegion(new SequenceRegionCoordinates(primer.getLeftPrimerPosition(), primer.getLeftPrimerPosition()+primer.getLeftPrimer().toCharArray().length));
				p3sit.addExcludedRegion(new SequenceRegionCoordinates(primer.getRightPrimerPosition()-primer.getRightPrimer().toCharArray().length, primer.getRightPrimerPosition()));
				//p3sit.addExcludedRegion(new SequenceRegionCoordinates(primer.getLeftPrimerPosition()+100, primer.getRightPrimerPosition()-100));
				//logger.info(mRNA.getId());
				primers.add(primer);
			}
			}
			//return pp;
		}
		
		p3io.endPrimer3Communications();
		
		Collection<PrimerPair> rtrn=new TreeSet<PrimerPair>();
		for(PrimerPair pair: primers){
			if(pair.getLeftPrimer()!=null && pair.getRightPrimer()!=null){rtrn.add(pair);}
		}
		
		return rtrn;
	}
	
	public static Collection<PrimerPair> designCloningPrimers(Sequence geneSequence, int numDesigns, String pathPrimer3core, Primer3Configuration best) throws Exception {
		return designCloningPrimers(geneSequence, numDesigns, null, pathPrimer3core, best);
	}
	
	public static Collection<PrimerPair> designCloningPrimers(Sequence geneSequence, int numDesigns, SequenceRegionCoordinates target, String pathPrimer3core) throws Exception {
		Primer3Configuration best = Primer3ConfigurationFactory.getLongRangePCRConfiguration();
		return designCloningPrimers(geneSequence, numDesigns, target, pathPrimer3core, best);
	}
		
	
	

	public static Collection<PrimerPair> designSyntheticPrimers(String seq, int numDesigns, String pathPrimer3core, double optimalMeltingTemp) throws IOException {
		Primer3Configuration best = Primer3ConfigurationFactory.getSyntheticConfiguration(optimalMeltingTemp);
		
		if(numDesigns>0){best.setPrimerNumReturn(numDesigns);}
				
		Primer3IO p3io = new Primer3IO(pathPrimer3core);
		p3io.startPrimer3Communications();
			
		
				
		Sequence mRNA=new Sequence("Gene", seq);
		
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(mRNA.getName());
		
				
		Collection<PrimerPair> pp = p3io.findPrimerPair(p3sit, best);
		

		p3io.endPrimer3Communications();
		
		TreeSet<PrimerPair> rtrn=new TreeSet<PrimerPair>();
		for(PrimerPair pair: pp){
			if(pair.getLeftPrimer()!=null && pair.getRightPrimer()!=null){rtrn.add(pair);}
		}
		
		return rtrn;
	}
	
	/**
	 * Design a primer pair that flanks a region of interest
	 * One primer is in each flanking region
	 * Neither primer overlaps the actual region
	 * @param config primer3 configuration
	 * @param sequence The sequence containing the region of interest
	 * @param pathPrimer3core primer3core executable
	 * @param regionStart Start position of region of interest
	 * @param regionEnd Position after last position of region of interest
	 * @param flankingRegionSize Size of flanking regions to search for primers
	 * @return Primer pair with one primer in each flank or null if none exist
	 * @throws IOException
	 */
	public static PrimerPair designPrimerPairFlankingWindow(Primer3Configuration config, Sequence sequence, String pathPrimer3core, int regionStart, int regionEnd, int flankingRegionSize) throws IOException {

		Sequence leftFlank = sequence.getSubSequence("", regionStart - flankingRegionSize, regionStart);
		Sequence rightFlank = sequence.getSubSequence("", regionEnd, regionEnd + flankingRegionSize);
		
		// Make a string of Ns to put between the two flanking regions
		int innerLength = Math.max(regionEnd - regionStart, flankingRegionSize);
		char[] ns = new char[innerLength];
		for(int i = 0; i < ns.length; i++) ns[i] = 'N';
		String nStr = new String(ns);
		
		String modifiedSequence = leftFlank.getSequenceBases() + nStr + rightFlank.getSequenceBases(); // Sequence to design primers against
		
		// Change the product size in the primer3 config
		config.minProductSize = innerLength;
		config.maxProductSize = modifiedSequence.length();
		
		// Change number of primers to return
		config.maxNumPrimersToReturn = 1;
		
		// Get the primer
		return designBestPrimer(config, modifiedSequence, pathPrimer3core);
		
	}
	
	
	public static PrimerPair designPrimerPairAroundSNP(Primer3Configuration config, String sequence, String pathPrimer3core, int regionStart, int regionEnd, int flankingRegionSize) throws IOException {

		String leftFlank = sequence.substring(0, regionStart);
		String rightFlank = sequence.substring(regionEnd, sequence.length());
		
		// Make a string of Ns to put between the two flanking regions
		int innerLength = Math.max(regionEnd - regionStart, flankingRegionSize);
		char[] ns = new char[innerLength];
		for(int i = 0; i < ns.length; i++) ns[i] = 'N';
		String nStr = new String(ns);
		
		String modifiedSequence = leftFlank + nStr + rightFlank; // Sequence to design primers against
		
		// Change the product size in the primer3 config
		config.minProductSize = innerLength;
		config.maxProductSize = modifiedSequence.length();
		
		// Change number of primers to return
		config.maxNumPrimersToReturn = 1;
		
		// Get the primer
		return designBestPrimer(config, modifiedSequence, pathPrimer3core);
		
	}
	
	/**
	 * Design primers using primer3
	 * @param config Primer3 configuration
	 * @param seq Sequence
	 * @param pathPrimer3core primer3core executable
	 * @return Collection of primers
	 * @throws IOException
	 */
	public static Collection<PrimerPair> designPrimers(Primer3Configuration config, String seq, String pathPrimer3core) throws IOException {
		Sequence sequence = new Sequence("", seq);
		return designPrimers(config, sequence, pathPrimer3core);
	}
	
	/**
	 * Design primers that amplify the entire sequence using primer3 and get the one with the lowest penalty
	 * @param seq Sequence
	 * @param pathPrimer3core primer3core executable
	 * @return Primer with lowest penalty whose product is the entire sequence, or null if there are no valid primers
	 * @throws IOException
	 */
	public static PrimerPair designBestPrimersFullSequence(Sequence seq, String pathPrimer3core) throws IOException {
		Primer3Configuration config = Primer3ConfigurationFactory.getSpecificLengthPCRConfiguration(seq.getLength());
		return designBestPrimer(config, seq, pathPrimer3core);
	}
	
	/**
	 * Design primers using primer3
	 * @param config Primer3 configuration
	 * @param seq Sequence
	 * @param pathPrimer3core primer3core executable
	 * @return Collection of primers
	 * @throws IOException
	 */
	public static Collection<PrimerPair> designPrimers(Primer3Configuration config, Sequence seq, String pathPrimer3core) throws IOException {
		Primer3IO p3io = new Primer3IO(pathPrimer3core);
		p3io.startPrimer3Communications();
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(seq);	
		p3sit.setPrimerSequenceId(seq.getName());
		Collection<PrimerPair> pp = p3io.findPrimerPair(p3sit, config);
		p3io.endPrimer3Communications();
		TreeSet<PrimerPair> rtrn=new TreeSet<PrimerPair>();
		for(PrimerPair pair: pp){
			if(pair.getLeftPrimer()!=null && pair.getRightPrimer()!=null){rtrn.add(pair);}
		}
		return rtrn;
	}
	/**
	 * Design primers using primer3 and get the one with the lowest penalty
	 * @param config Primer3 configuration
	 * @param seq Sequence
	 * @param pathPrimer3core primer3core executable
	 * @return Primer with lowest penalty, or null if there are no primers
	 * @throws IOException
	 */
	public static PrimerPair designBestPrimer(Primer3Configuration config, String seq, String pathPrimer3core) throws IOException {
		Sequence sequence = new Sequence("", seq);
		return designBestPrimer(config, sequence, pathPrimer3core);
	}
	
	/**
	 * Design primers using primer3 and get the one with the lowest penalty
	 * @param config Primer3 configuration
	 * @param seq Sequence
	 * @param pathPrimer3core primer3core executable
	 * @return Primer with lowest penalty, or null if there are no primers
	 * @throws IOException
	 */
	public static PrimerPair designBestPrimer(Primer3Configuration config, Sequence seq, String pathPrimer3core) throws IOException {
		Collection<PrimerPair> allPrimers = designPrimers(config, seq, pathPrimer3core);
		if(allPrimers.isEmpty()) {
			return null;
		}
		float minPenalty = Float.MAX_VALUE;
		PrimerPair bestPrimer = null;
		for(PrimerPair primer : allPrimers) {
			if(primer.getPrimerPairPenalty() < minPenalty) {
				minPenalty = primer.getPrimerPairPenalty();
				bestPrimer = primer;
			}
		}
		return bestPrimer;
	}
	
	


	public static Collection<PrimerPair> designSyntheticPrimers(String seq, String sequencePrimer, String sequencePrimerRevComp, int numDesigns, String pathPrimer3core, double optimalMeltingTemp) throws Exception {
		Primer3Configuration best = Primer3ConfigurationFactory.getSyntheticConfiguration2(optimalMeltingTemp);
		
		if(numDesigns>0){best.setPrimerNumReturn(numDesigns);}
				
		Primer3IO p3io = new Primer3IO(pathPrimer3core);
		p3io.startPrimer3Communications();
			
		Sequence mRNA=new Sequence("Primer", seq);
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(mRNA.getName());
		if(sequencePrimer!=null && !sequencePrimer.isEmpty()){p3sit.setPrimerLeftInput(sequencePrimer);}
		if(sequencePrimerRevComp!=null && !sequencePrimerRevComp.isEmpty()){p3sit.setPrimerRightInput(sequencePrimerRevComp);}
		
		Collection<PrimerPair> pp = p3io.findPrimerPair(p3sit, best);
		

		p3io.endPrimer3Communications();
		
		TreeSet<PrimerPair> rtrn=new TreeSet<PrimerPair>();
		for(PrimerPair pair: pp){
			if(pair.getLeftPrimer()!=null && pair.getRightPrimer()!=null){rtrn.add(pair);}
		}
		
		return rtrn;
	}
	
	private static final String SYNTHETIC_CONFIG_NAME = "synthetic";
	private static final String QPCR_CONFIG_NAME = "qPCR";
	private static final String RAP_QPCR_CONFIG_NAME = "RAPqPCR";
	private static final String DELETION_PLASMID_CONFIG_NAME = "deletion_plasmid";
	
	private static final String[] CONFIG_NAMES = {SYNTHETIC_CONFIG_NAME,QPCR_CONFIG_NAME,RAP_QPCR_CONFIG_NAME,DELETION_PLASMID_CONFIG_NAME};
	
	private static String configNamesList() {
		StringBuilder sb = new StringBuilder(CONFIG_NAMES[0]);
		for(int i = 1; i < CONFIG_NAMES.length; i++) {
			sb.append(", " + CONFIG_NAMES[i]);
		}
		return sb.toString();
	}
	
	
	public static void main(String [] args) throws IOException{
		Primer3Configuration primer3config = Primer3ConfigurationFactory.getRAPqPCRConfiguration();
		String seq=args[0];
		String primer3core=args[1];
		PrimerPair primer=PcrPrimerDesigner.designBestPrimer(primer3config, seq, primer3core);
		System.err.println(primer.getLeftPrimer()+" "+primer.getRightPrimer());
	}
	
	/*public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-s","Fasta file of sequences to design primers against",true);
		p.addStringArg("-r","Optional file of region coordinates to design primers against (format: sequence_name start_pos end_pos)",false,null);
		p.addStringArg("-c","Primer3 configuration name (" + configNamesList() + ")",true);
		p.addBooleanArg("-rc", "Design primers against antisense strand", true);
		p.addStringArg("-o", "Outfile", true);
		p.addStringArg("-p3c", "Primer3core executable", true);
		
		p.parse(args);
		
		String inputfile = p.getStringArg("-s");
		String regions = p.getStringArg("-r");
		String config = p.getStringArg("-c");
		String outfile = p.getStringArg("-o");
		boolean rc = p.getBooleanArg("-rc");
		String primer3core = p.getStringArg("-p3c");
		
		
		boolean configOk = false;
		for(int i=0; i<CONFIG_NAMES.length; i++) {
			if(config.equals(CONFIG_NAMES[i])) configOk = true;
		}
		
		// validate configuration name
		if(!configOk) {
			logger.info("\nValid Primer3 configuration names:");
			for(int i=0; i<CONFIG_NAMES.length; i++) {
				logger.info(CONFIG_NAMES[i]);
			}
			logger.info("");
			throw new IllegalArgumentException();
		}
		
		// the sequences to design primers for
		ArrayList<Sequence> seqs = new ArrayList<Sequence>();
		
		// the input sequences
		FastaSequenceIO fsio = new FastaSequenceIO(inputfile);
		List<Sequence> inputseqs = fsio.loadAll();
		
		// if regions are not specified, just keep the whole sequences
		if(regions == null) {
			for(Sequence seq : inputseqs) {
				if(rc) seqs.add(seq.getAntisense());
				else seqs.add(seq);
			}
		}
		
		// if regions are specified, extract the regions
		if(regions != null) {
			
			FileReader r = new FileReader(regions);
			BufferedReader b = new BufferedReader(r);
			
			StringParser s = new StringParser();
			
			while(b.ready()) {
				
				String line = b.readLine();
				s.parse(line);
				if(s.getFieldCount() != 3 && s.getFieldCount() != 4) {
					b.close();
					throw new IllegalArgumentException("Region line is not valid:\n" + line + "\nFormat: sequence_name start_pos end_pos <optional ID>");
				}
				String seqname = s.asString(0);
				int start = s.asInt(1);
				int end = s.asInt(2);
				
				Collection<String> seqnames = new ArrayList<String>();
				seqnames.add(seqname);
				
				Collection<Sequence> theseqs = fsio.extractRecords(seqnames);
				Iterator<Sequence> iter = theseqs.iterator();
				Sequence theseq = iter.next();
				
				Sequence subseq = theseq.getSubSequence(theseq.getId() + ":" + start + "-" + end, start, end);
				if(s.getFieldCount() == 4) {
					subseq.setId(s.asString(3));
				}
				if(rc) seqs.add(subseq.getAntisense());
				else seqs.add(subseq);
				
			}
			
			b.close();
			
		}
		
		// set up primer3
		Primer3Configuration primer3config = new Primer3Configuration();
		
		if(config.equals(SYNTHETIC_CONFIG_NAME)) primer3config = Primer3ConfigurationFactory.getSyntheticConfiguration(60);
		if(config.equals(RAP_QPCR_CONFIG_NAME)) primer3config = Primer3ConfigurationFactory.getRAPqPCRConfiguration();
		if(config.equals(QPCR_CONFIG_NAME)) primer3config = Primer3ConfigurationFactory.getQpcrConfiguration();
		// If plasmids, will get separate configuration for each plasmid
		if(config.equals(DELETION_PLASMID_CONFIG_NAME)) primer3config = null;
		
		// output file
		FileWriter writer = new FileWriter(outfile);
		
		String header = "primer_ID\t";
		header += "left_primer\t";
		header += "right_primer\t";
		header += "left_primer_TM\t";
		header += "right_primer_TM\t";
		header += "primer_pair_penalty";
		writer.write(header + "\n");

		for(Sequence seq : seqs) {
			
			PrimerPair primer = PcrPrimerDesigner.designBestPrimer(primer3config, seq, primer3core);
			String lineToWrite = seq.getId() + "\t";
			if(primer == null) {
				lineToWrite += "NO_PRIMERS";
				writer.write(lineToWrite + "\n");
				continue;
			}
			lineToWrite += primer.getLeftPrimer().toUpperCase() + "\t";
			lineToWrite += primer.getRightPrimer().toUpperCase() + "\t";
			lineToWrite += primer.getLeftPrimerTM() + "\t";
			lineToWrite += primer.getRightPrimerTM() + "\t";
			lineToWrite += primer.getPrimerPairPenalty();
			writer.write(lineToWrite + "\n");

			
		}
		
		writer.close();
		
		logger.info("");
		logger.info("All done.");
		
	}*/
	
	
	
}
