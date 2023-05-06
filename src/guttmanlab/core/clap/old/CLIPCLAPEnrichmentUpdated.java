package guttmanlab.core.clap.old;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.PopulatedWindow;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.PairedMappedFragment;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotation.predicate.UniqueMapperFilter;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMPairedFragmentCollection;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import guttmanlab.core.annotationcollection.FeatureCollection;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.math.ScanStat;
import jsc.distributions.Binomial;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.util.CloseableIterator;

public class CLIPCLAPEnrichmentUpdated {
	
	static String totalNumberReadsTag="totalreads";
	static String numWindowsTag="numberofwindows";
	
	static String totalNumberReadsTagMouse="mousereads";
	static String numWindowsTagMouse="mousewindows";
	static String totalNumberReadsTagHuman="humanreads";
	static String numWindowsTagHuman="humanwindows";
	
	private int clipTotalCounts;
	private int inputTotalCounts;
	private int clipWindows;
	private int inputWindows;
	
	private int minReadCounts=0;
	private double percentile;
	private int stepSize;
	private int windowSize;
	
	private boolean writeBED;
	private FileWriter writerTotal;
	private FileWriter writerBEDPos;
	private FileWriter writerBEDNeg;
	
	private boolean isMixedSample;
	
	private int mouseClipTotalCounts;
	private int mouseInputTotalCounts;;
	private int mouseInputWindows;
	private int mouseClipWindows;
	
	private int humanClipTotalCounts;
	private int humanInputTotalCounts;;
	private int humanInputWindows;
	private int humanClipWindows;
	
	private double sigCutoff=0.0001;
	private Collection<String> written;
	
	

	/*public CLIPCLAPEnrichmentUpdated(File clipBAM, File inputBAM, int windowSize, String save, Collection<? extends Annotation> regions, double percentile, int stepSize, boolean skipIntrons, boolean writeBED, boolean isMixedSample) throws IOException{
		this.writeBED=writeBED;
		this.isMixedSample=isMixedSample;
		this.windowSize=windowSize;
		this.stepSize=stepSize;
		this.percentile=percentile;
		if(stepSize==1){this.minReadCounts=-1;}
		if(isMixedSample){setTags(clipBAM, inputBAM);}
		else{getCounts(clipBAM, inputBAM);}
		this.written=new TreeSet<String>();
		
		//computeEnrichment(clipBAM, inputBAM, regions, skipIntrons, save);
		System.err.println("Compute Enrichment");
		computeEnrichmentSpeedup(clipBAM, inputBAM, regions, skipIntrons, save);
		
	}*/
	
	public CLIPCLAPEnrichmentUpdated(File clipBAM, File inputBAM, int windowSize, String save, Collection<? extends Annotation> regions, double percentile, int stepSize, boolean skipIntrons, boolean writeBED, boolean isMixedSample) throws IOException{
		this.writeBED=writeBED;
		this.isMixedSample=isMixedSample;
		this.windowSize=windowSize;
		this.stepSize=stepSize;
		this.percentile=percentile;
		if(stepSize==1){this.minReadCounts=-1;}
		if(isMixedSample){setTags(clipBAM, inputBAM);}
		else{getCounts(clipBAM, inputBAM);}
		this.written=new TreeSet<String>();
		
		System.err.println("Compute Enrichment");
		computeEnrichmentSpeedupBuffered(clipBAM, inputBAM, regions, skipIntrons, save, writeBED);
		
	}
	
	
	private void getCounts(File clipBAM, File inputBAM) {
		WindowCounts clipCount=new WindowCounts(clipBAM, this.windowSize);
		this.clipWindows=clipCount.getNumberOfWindows();
		this.clipTotalCounts=clipCount.getTotalNumberOfReads();
		
		
		WindowCounts inputCount=new WindowCounts(inputBAM, this.windowSize);
		this.inputWindows=inputCount.getNumberOfWindows();
		this.inputTotalCounts=inputCount.getTotalNumberOfReads();
	}

	private void setTags(File clipBAM, File inputBAM) {
		this.clipTotalCounts=getTotalNumberOfReads(clipBAM, totalNumberReadsTag);
		this.inputTotalCounts=getTotalNumberOfReads(inputBAM, totalNumberReadsTag);
		
		this.inputWindows=getTotalNumberOfReads(inputBAM, numWindowsTag);
		this.clipWindows=getTotalNumberOfReads(clipBAM, numWindowsTag);
		
		if(isMixedSample){
			this.mouseClipTotalCounts=getTotalNumberOfReads(clipBAM, totalNumberReadsTagMouse);
			this.mouseInputTotalCounts=getTotalNumberOfReads(inputBAM, totalNumberReadsTagMouse);
			this.mouseInputWindows=getTotalNumberOfReads(inputBAM, numWindowsTagMouse);
			this.mouseClipWindows=getTotalNumberOfReads(clipBAM, numWindowsTagMouse);
			
			
			this.humanClipTotalCounts=getTotalNumberOfReads(clipBAM, totalNumberReadsTagHuman);
			this.humanInputTotalCounts=getTotalNumberOfReads(inputBAM, totalNumberReadsTagHuman);
			this.humanInputWindows=getTotalNumberOfReads(inputBAM, numWindowsTagHuman);
			this.humanClipWindows=getTotalNumberOfReads(clipBAM, numWindowsTagHuman);
		}
	}

	
	private int getTotalNumberOfReads(File bamFile, String tag) {
		BAMSingleReadCollection bam=new BAMSingleReadCollection(bamFile);
		SAMFileHeader header=bam.getFileHeader();
		List<String> comments=header.getComments();
		for(String co: comments){
			System.err.println(co);
			String line=co.split("\t")[1];
			String parsedTag=line.split(":")[0];
			if(parsedTag.equalsIgnoreCase(tag)){
				return new Integer(line.split(":")[1]);
			}
		}
		throw new IllegalArgumentException("Could not find "+ tag);
		//TODO compute the tags and add to file
	}
	
	
	private void computeEnrichment(File sampleBamFile, File inputBamFile, Collection<? extends Annotation> regions, boolean skipIntrons, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write("window\tgene\tnormalized score\tglobal enrichment\tsample count\tinput count\t max input count\twindow normalized enrichment\twindow p-value\tlocal p-value\tlocal scan p-value\n");	
		
		
		Pair<FileWriter> bedWriter=new Pair<FileWriter>();
		if(writeBED){bedWriter=new Pair<FileWriter>(new FileWriter(save+".pos.bedgraph"), new FileWriter(save+".neg.bedgraph"));}
		
		
		int count=0;
		
		for(Annotation a: regions){
			Map<Annotation, PairedScore> scores=computeEnrichment(sampleBamFile, inputBamFile, a.getOrientation(), a);
			write(writer, scores);
			if(writeBED){writeBEDGraph(bedWriter, scores);}
			if(count%1000 ==0){System.err.println(count+" "+regions.size()+" "+a.getName()+" "+a.toUCSC());}
			count++;
		}
		
		count=0;
		if(!skipIntrons){
			for(Annotation a: regions){
				Collection<Annotation> introns=a.getIntrons();
				for(Annotation intron: introns){
					Annotation newIntron=new SingleInterval(intron.getReferenceName(), intron.getReferenceStartPosition()+1, intron.getReferenceEndPosition(), intron.getOrientation(), intron.getName());
					Map<Annotation, PairedScore> scores=computeEnrichment(sampleBamFile, inputBamFile, newIntron.getOrientation(), newIntron);
					write(writer, scores);
					if(writeBED){writeBEDGraph(bedWriter, scores);}
				}
				
				if(count %1000 ==0){System.err.println(count+" "+regions.size() +" introns "+a.getName());}
				count++;
			}
			
		}
		
		if(writeBED){bedWriter.getValue1().close(); bedWriter.getValue2().close();}
		writer.close();
	}
	
	
	private void computeEnrichmentSpeedupBuffered(File sampleBamFile, File inputBamFile, Collection<? extends Annotation> regions, boolean skipIntrons, String save, boolean writeBedgraph) throws IOException {
		//Add introns
		Collection<Annotation> regionsToUse=addIntrons(regions, skipIntrons);
		
		//Score all regions
		FeatureCollection<SingleInterval> geneCounts=scoreAllWindowsMap(inputBamFile, regionsToUse);
		System.err.println("Scored Genes");
		
		int numWindows=numWindows(regionsToUse);
		
		
		FeatureCollection<SingleInterval> sampleCollection;
		FeatureCollection<SingleInterval> inputCollection;
		FeatureCollection<SingleInterval> pairedInputMap;
		
		System.err.println("Num windows in genes "+numWindows+" total windows "+this.clipWindows);
		if(numWindows<this.clipWindows){
			System.err.println("Scoring by gene windows");
			//TODO Get all windows in regions
			FeatureCollection<Annotation> windows=getWindows(regionsToUse);
			System.err.println("Got windows "+windows.size());
			
			//Score all windows in sample
			sampleCollection=scoreAllWindowsMap(sampleBamFile, windows.getCollection());
			System.err.println("Scored samples");
			//Score all windows in input
			pairedInputMap=scoreAllWindowsMap(inputBamFile, sampleCollection.getCollection());
			inputCollection=pairedInputMap;
			System.err.println("Scored Input");
		}
		
		else{
			System.err.println("Scoring ALL windows");
			sampleCollection=makeWindowCollectionMap(sampleBamFile);
			pairedInputMap=scoreAllWindowsMap(inputBamFile, sampleCollection.getCollection());
			inputCollection=makeWindowCollectionMap(inputBamFile);
			System.err.println("Scored Input");
		}
		
		
		
		
		/*System.err.println("Make window collection");
		FeatureCollection<SingleInterval> sampleCollection=makeWindowCollectionMap(sampleBamFile);
		
		System.err.println("Score input");	
		//score all windows in input
		FeatureCollection<SingleInterval> pairedInputMap=scoreAllWindowsMap(inputBamFile, sampleCollection.getCollection());
		FeatureCollection<SingleInterval> inputCollection=makeWindowCollectionMap(inputBamFile);
		
		
		
		Collection<Annotation> regionsToUse=addIntrons(regions, skipIntrons);
		FeatureCollection<SingleInterval> geneCounts=scoreAllWindowsMap(inputBamFile, regionsToUse);*/
		
		FileWriter writer=new FileWriter(save);
		FileWriter writerBedPos=new FileWriter(save+".pos.bedgraph");
		FileWriter writerBedNeg=new FileWriter(save+".neg.bedgraph");
		writer.write("window\tgene\tnormalized score\tglobal enrichment\tsample count\tinput count\t max input count\twindow normalized enrichment\twindow p-value\tlocal p-value\tlocal scan p-value\n");	
		
		System.err.println("Scoring regions");
		
		System.err.println("Iterate over regions");
		
		int counter=0;
		//iterate through genes and get windows
		for(Annotation a: regionsToUse){
			CloseableIterator<SingleInterval> sampleWindows=sampleCollection.sortedIterator(a, true);
			Map<Annotation, PairedScore> scores=computeEnrichmentMap(sampleWindows, pairedInputMap, inputCollection, geneCounts, a.getOrientation(), a);
			write(writer, scores);
			writeBedGraph(writerBedPos, writerBedNeg, scores, writeBedgraph);
			counter++;
			if(counter%1000==0){System.err.println("processing "+counter);}
		}
		
		writer.close();
		writerBedPos.close();
		writerBedNeg.close();
	}
	
	
	
	/*private void computeEnrichmentSpeedup(File sampleBamFile, File inputBamFile, Collection<? extends Annotation> regions, boolean skipIntrons, String save) throws IOException {
		System.err.println("Make window collection");
		FeatureCollection<PopulatedWindow<SAMFragment>> sampleCollection=makeWindowCollection(sampleBamFile);
		
		System.err.println("Score input");	
		//score all windows in input
		FeatureCollection<PopulatedWindow<SAMFragment>> pairedInputMap=scoreAllWindows(inputBamFile, sampleCollection);
		FeatureCollection<PopulatedWindow<SAMFragment>> inputCollection=makeWindowCollection(inputBamFile);
		
		Collection<Annotation> regionsToUse=addIntrons(regions, skipIntrons);
		
		FileWriter writer=new FileWriter(save);
		writer.write("window\tgene\tnormalized score\tglobal enrichment\tsample count\tinput count\t max input count\twindow normalized enrichment\twindow p-value\tlocal p-value\tlocal scan p-value\n");	
		
		System.err.println("Scoring regions");
		
		System.err.println("Iterate over regions");
		
		int counter=0;
		//iterate through genes and get windows
		for(Annotation a: regionsToUse){
			CloseableIterator<PopulatedWindow<SAMFragment>> sampleWindows=sampleCollection.sortedIterator(a, true);
			Map<Annotation, PairedScore> scores=computeEnrichment(sampleWindows, pairedInputMap, inputCollection, a.getOrientation(), a);
			write(writer, scores);
			counter++;
			if(counter%1000==0){System.err.println("processing "+counter);}
		}
		
		writer.close();
	}*/
	
	
	private int numWindows(Collection<Annotation> regionsToUse) {
		int size=0;
		for(Annotation region: regionsToUse){
			size+=region.size();
		}
		return size/this.stepSize;
	}


	private FeatureCollection<Annotation> getWindows(Collection<Annotation> regionsToUse) {
		FeatureCollection<Annotation> rtrn=new FeatureCollection<Annotation>();
		
		for(Annotation region: regionsToUse){
			CloseableIterator<DerivedAnnotation<? extends Annotation>> windows=region.getWindows(windowSize, stepSize).sortedIterator();
			while(windows.hasNext()){
				DerivedAnnotation<?extends Annotation> window=windows.next();
				rtrn.add(window);
			}
		}
		
		return rtrn;
	}


	private FeatureCollection<SingleInterval> scoreAllWindowsMap(File inputBamFile, Collection<? extends Annotation> sampleCollection) {
		BAMSingleReadCollection bam=new BAMSingleReadCollection(inputBamFile);
		bam.addFilter(new UniqueMapperFilter());
		
		FeatureCollection<SingleInterval> windows=bam.scoreAllWindows(sampleCollection);
		
		
		bam.close();
		return windows;
	}
	
	
	/*private FeatureCollection<PopulatedWindow<SAMFragment>> scoreAllWindows(File inputBamFile, FeatureCollection<PopulatedWindow<SAMFragment>> sampleCollection) {
		BAMSingleReadCollection bam=new BAMSingleReadCollection(inputBamFile);
		bam.addFilter(new UniqueMapperFilter());
		
		Collection<Annotation> regions=sampleCollection.getCollection();
		FeatureCollection<PopulatedWindow<SAMFragment>> windows=bam.scoreAllWindows(regions);
		
		
		bam.close();
		return windows;
	}*/


	private Collection<Annotation> addIntrons(Collection<? extends Annotation> regions, boolean skipIntrons) {
		Collection<Annotation> regionsToUse=new TreeSet<Annotation>();
		regionsToUse.addAll(regions);
		if(!skipIntrons){
			for(Annotation region: regions){
				regionsToUse.addAll(region.getIntrons());
			}
		}
		return regionsToUse;
	}

	
	private Map<Annotation, PairedScore> computeEnrichmentMap(CloseableIterator<SingleInterval> sampleWindows, FeatureCollection<SingleInterval> pairedInputMap, FeatureCollection<SingleInterval> inputCollection, FeatureCollection<SingleInterval> geneCounts, Strand orientation, Annotation region) {
		Map<Annotation, PairedScore> rtrn=new TreeMap<Annotation, PairedScore>();
		
		int geneSampleTotal=0; //sum across all sample windows
		//For each window get scores in sample and input
		while(sampleWindows.hasNext()){
			SingleInterval window=sampleWindows.next();
			window.setOrientation(orientation);
			int sampleCount=window.getCount();
			int inputCount=pairedInputMap.get(window).getCount(); 
			
			PairedScore clipScore=new PairedScore(window, sampleCount, inputCount, this.clipTotalCounts, this.inputTotalCounts);
			clipScore.setParentAnnotation(region);
			clipScore.setNumberInputWindows(this.inputWindows);
			clipScore.setNumberSampleWindows(this.clipWindows);
			rtrn.put(window, clipScore);
			geneSampleTotal+=sampleCount;
		}
		sampleWindows.close();
		
		
		//System.err.println(region);
		//compute global scores from entire iterator
		int inputPercentile=percentile2(inputCollection.sortedIterator(region, true), orientation);
		int geneInputTotal=geneCounts.get(region.getSingleInterval()).getCount();
				
		//Update each window with constants		
		for(Annotation window: rtrn.keySet()){
			PairedScore score=rtrn.get(window);
			score.setInputPercentile(inputPercentile);
			score.setGeneSampleTotal(geneSampleTotal);
			score.setGeneInputTotal(geneInputTotal);
		}		
		
		return rtrn;
	}

	private Map<Annotation, PairedScore> computeEnrichment(CloseableIterator<PopulatedWindow<SAMFragment>> sampleWindows, FeatureCollection<PopulatedWindow<SAMFragment>> pairedInputMap, FeatureCollection<PopulatedWindow<SAMFragment>> inputCollection, Strand orientation, Annotation region) {
		Map<Annotation, PairedScore> rtrn=new TreeMap<Annotation, PairedScore>();
		
		int geneSampleTotal=0; //sum across all sample windows
		//For each window get scores in sample and input
		while(sampleWindows.hasNext()){
			PopulatedWindow<SAMFragment> window=sampleWindows.next();
			window.setOrientation(orientation);
			int sampleCount=window.getNumberOfAnnotationsInWindow(orientation);
			int inputCount=this.sum(pairedInputMap.sortedIterator(window, true), orientation); //TODO Score all windows
			//System.err.println(window.toUCSC(window.getOrientation())+" "+inputCount);
			
			PairedScore clipScore=new PairedScore(window, sampleCount, inputCount, this.clipTotalCounts, this.inputTotalCounts);
			clipScore.setParentAnnotation(region);
			clipScore.setNumberInputWindows(this.inputWindows);
			clipScore.setNumberSampleWindows(this.clipWindows);
			rtrn.put(window, clipScore);
			geneSampleTotal+=sampleCount;
		}
		sampleWindows.close();
		
		
		
		//compute global scores from entire iterator
		int inputPercentile=percentile(inputCollection.sortedIterator(region, true), orientation);
		int geneInputTotal=sum(inputCollection.sortedIterator(region, true), orientation); //sum across all input windows;
				
		//Update each window with constants		
		for(Annotation window: rtrn.keySet()){
			PairedScore score=rtrn.get(window);
			score.setInputPercentile(inputPercentile);
			score.setGeneSampleTotal(geneSampleTotal);
			score.setGeneInputTotal(geneInputTotal);
		}		
		
		return rtrn;
	}
	
	


	private int sum(CloseableIterator<PopulatedWindow<SAMFragment>> sortedIterator, Strand orientation) {
		int count=0;
		
		Collection<String> visited=new TreeSet<String>();
		while(sortedIterator.hasNext()){
			PopulatedWindow<SAMFragment> window=sortedIterator.next();
			if(!visited.contains(window.toUCSC())){
				count+=window.getNumberOfAnnotationsInWindow(orientation);
			}
			visited.add(window.toUCSC());
		}
		sortedIterator.close();
		
		return count;
	}
	
	private int sum2(CloseableIterator<SingleInterval> sortedIterator, Strand orientation) {
		int count=0;
		
		Collection<String> visited=new TreeSet<String>();
		while(sortedIterator.hasNext()){
			SingleInterval window=sortedIterator.next();
			if(!visited.contains(window.toUCSC())){
				count+=window.getCount();
			}
			visited.add(window.toUCSC());
		}
		sortedIterator.close();
		
		return count;
	}

	
	private FeatureCollection<SingleInterval> makeWindowCollectionMap(File sampleBamFile) throws IOException {
		BAMSingleReadCollection bam=new BAMSingleReadCollection(sampleBamFile);
		bam.addFilter(new UniqueMapperFilter());
		CloseableIterator<? extends PopulatedWindow<SAMFragment>> windows1=bam.getPopulatedWindows(windowSize, stepSize);
		
		
		FeatureCollection<SingleInterval> rtrn=new FeatureCollection<SingleInterval>();
		//Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		
		int index=0;
		while(windows1.hasNext()){
			PopulatedWindow<SAMFragment> window=windows1.next();
			SingleInterval newInterval=new SingleInterval(window.getReferenceName(), window.getReferenceStartPosition(), window.getReferenceEndPosition(), window.getOrientation());
			newInterval.setCount(window.getNumberOfAnnotationsInWindow());
			rtrn.add(newInterval);
			
			index++;
			if(index%1000000 ==0){System.err.println("processing "+index);}
		}
		windows1.close();
		bam.close();
		
		return rtrn;
		//return null;
	}
	

	private FeatureCollection<PopulatedWindow<SAMFragment>> makeWindowCollection(File sampleBamFile) throws IOException {
		FileWriter temp=new FileWriter("temp");
		BAMSingleReadCollection bam=new BAMSingleReadCollection(sampleBamFile);
		bam.addFilter(new UniqueMapperFilter());
		CloseableIterator<? extends PopulatedWindow<SAMFragment>> windows1=bam.getPopulatedWindows(windowSize, stepSize);
		FeatureCollection<PopulatedWindow<SAMFragment>> sampleCollection=new FeatureCollection<PopulatedWindow<SAMFragment>>();
		
		int index=0;
		while(windows1.hasNext()){
			PopulatedWindow<SAMFragment> window=windows1.next();
			temp.write(window.toUCSC(window.getOrientation())+"\t"+window.getNumberOfAnnotationsInWindow()+"\n");
			sampleCollection.addAnnotation(window);
			index++;
			if(index%1000000 ==0){System.err.println("processing "+index);}
		}
		windows1.close();
		bam.close();
		temp.close();
		return sampleCollection;
		//return null;
	}
	
	private FeatureCollection<PopulatedWindow<SAMFragment>> makeWindowCollection(File sampleBamFile, FeatureCollection<PopulatedWindow<SAMFragment>> inputMap) throws IOException {
		FeatureCollection<PopulatedWindow<SAMFragment>> rtrn=makeWindowCollection(sampleBamFile);
		
		CloseableIterator<PopulatedWindow<SAMFragment>> iter=inputMap.sortedIterator();
		while(iter.hasNext()){
			rtrn.addAnnotation(iter.next());
		}
		return rtrn;
	}


	private void writeBEDGraph(Pair<FileWriter> bedWriter, Map<Annotation, PairedScore> scores) throws IOException {
		for(Annotation window: scores.keySet()){
			PairedScore score=scores.get(window);
			FileWriter writer=bedWriter.getValue1();
			if(window.getOrientation().equals(Strand.NEGATIVE)){
				writer=bedWriter.getValue2();
			}
			writer.write(window.getReferenceName()+"\t"+window.getReferenceStartPosition()+"\t"+(window.getReferenceStartPosition()+this.stepSize)+"\t"+score.getWindowNormEnrichment()+"\n");
			
		}
	}
		
		
	
	private void writeBedGraph(FileWriter writerPos, FileWriter writerNeg, Map<Annotation, PairedScore> scores, boolean writeBedGraph) throws IOException{
		if(!writeBedGraph){return;}
		
		FileWriter writer=writerPos;
		for(Annotation window: scores.keySet()){
			PairedScore score=scores.get(window);
			if(window.getOrientation().equals(Strand.NEGATIVE)){writer=writerNeg;}
			int start=window.getReferenceStartPosition();
			int end=start+this.stepSize;
			writer.write(window.getReferenceName()+"\t"+start+"\t"+end+"\t"+score.getWindowNormEnrichment()+"\n");
		}
		writer.flush();
	}	
	
	private void write(FileWriter writer, Map<Annotation, PairedScore> scores) throws IOException{
		for(Annotation window: scores.keySet()){
			PairedScore score=scores.get(window);
			double maxP=Math.max(score.getWindowPVal(), score.getLocalPVal());
			//if(maxP<0.05){
				String line=window.toUCSC(window.getOrientation())+"\t"+score.toString();
				if(!written.contains(window.toUCSC(window.getOrientation()))){
					writer.write(line+"\n");
				}
				//else{System.err.println("Already written "+window.toUCSC(window.getOrientation()));}
				written.add(window.toUCSC(window.getOrientation()));
			//}
		}
		writer.flush();
	}	
		
	
	private Map<Annotation, PairedScore> computeEnrichment(File sampleBamFile, File inputBamFile, Strand strand, Annotation region) throws IOException {
		BAMSingleReadCollection sampleBam=new BAMSingleReadCollection(sampleBamFile);
		BAMSingleReadCollection inputBam=new BAMSingleReadCollection(inputBamFile);
		
		int geneSampleTotal=countOverlappers(sampleBam, region);//TODO Need to ensure this doesn't count reads where skipped regions overlaps region
		int geneInputTotal=countOverlappers(inputBam, region);
		
		
		//System.err.println(region.getName()+" "+region.toUCSC()+" "+geneSampleTotal+" "+geneInputTotal+" "+inputBam.numOverlappers(region, false));
		
		Map<Annotation, PairedScore> rtrn=new TreeMap<Annotation, PairedScore>();
		
		
		
		Map<Annotation, Integer> sampleMap=new TreeMap<Annotation, Integer>();
		Map<Annotation, Integer> inputMap=new TreeMap<Annotation, Integer>();
				
		//Iterate over windows
		CloseableIterator<? extends PopulatedWindow<SAMFragment>> windows1=sampleBam.getPopulatedWindows(region, windowSize, stepSize);
		
		while(windows1.hasNext()){
			PopulatedWindow<SAMFragment> window=windows1.next();
			window.setOrientation(strand);
			int count=window.getNumberOfAnnotationsInWindow();
			//TODO Make sure the window is fully within region
			sampleMap.put(window, count);
		}
		windows1.close();
		
		//Iterate over windows
		CloseableIterator<? extends PopulatedWindow<SAMFragment>> windows2=inputBam.getPopulatedWindows(region, windowSize, stepSize);
		while(windows2.hasNext()){
			PopulatedWindow<SAMFragment> window=windows2.next();
			window.setOrientation(strand);
			int count=window.getNumberOfAnnotationsInWindow();
			//System.err.println(window.toUCSC()+" "+count);
			if(window.getReferenceStartPosition()>=0){
				inputMap.put(window, count);
			}
		}
		windows2.close();
		sampleBam.close();
		inputBam.close();
		
		
		int inputPercentile=0; //Floor at 1 count
		if(!inputMap.isEmpty()){
			inputPercentile=percentile(inputMap);
		}
		
		for(Annotation r: sampleMap.keySet()){
			int sampleCount=sampleMap.get(r);
			int inputCount=0;
			if(inputMap.containsKey(r)){
				inputCount=inputMap.get(r);
			}
			else{
				inputCount=inputBam.numOverlappers(r, false);
			}
			
			PairedScore clipScore=new PairedScore(r, sampleCount, inputCount, this.clipTotalCounts, this.inputTotalCounts);
			clipScore.setInputPercentile(inputPercentile);
			clipScore.setParentAnnotation(region);
			clipScore.setNumberInputWindows(inputWindows);
			clipScore.setNumberSampleWindows(clipWindows);
			clipScore.setGeneSampleTotal(geneSampleTotal);
			clipScore.setGeneInputTotal(geneInputTotal);
			rtrn.put(r, clipScore);
		}
		
		return rtrn;
	}

	
	private int countOverlappers(BAMSingleReadCollection sampleBam, Annotation region) {
		CloseableIterator<SAMFragment> iter=sampleBam.sortedIterator(region, false);
		
		int count=0;
		while(iter.hasNext()){
			SAMFragment fragment=iter.next();
			//Does fragment overlap region?
			if(overlap(fragment, region)){count++;}
		}
		
		iter.close();
		return count;
	}


	private boolean overlap(SAMFragment fragment, Annotation region) {
		return fragment.overlaps(region);
	}


	private int percentile(Map<Annotation, Integer> inputMap) {
		List<Integer> values=new ArrayList<Integer>(inputMap.values());
		Collections.sort(values);
		
		int position=new Double(values.size()*this.percentile).intValue();
		return values.get(position);
	}
	
	private int percentile(CloseableIterator<PopulatedWindow<SAMFragment>> inputMap, Strand orientation) {
		List<Integer> values=new ArrayList<Integer>();
		if(!inputMap.hasNext()){return 0;}
		while(inputMap.hasNext()){
			values.add(inputMap.next().getNumberOfAnnotationsInWindow(orientation));
		}
		inputMap.close();
		
		Collections.sort(values);
		
		int position=new Double(values.size()*this.percentile).intValue();
		return values.get(position);
	}
	
	private int percentile2(CloseableIterator<SingleInterval> inputMap, Strand orientation) {
		List<Integer> values=new ArrayList<Integer>();
		if(!inputMap.hasNext()){return 0;}
		while(inputMap.hasNext()){
			values.add(inputMap.next().getCount());
		}
		inputMap.close();
		
		Collections.sort(values);
		
		int position=new Double(values.size()*this.percentile).intValue();
		return values.get(position);
	}

	
	public class PairedScore{
		private int sampleCount;
		private int inputCount;
		private int sampleTotal;
		private int inputTotal;
		Annotation parentAnnotation;
		private int inputPercentile;
		boolean hasPercentile;
		private int numberInputWindows;
		private int numberSampleWindows;
		private int geneSampleTotal;
		private int geneInputTotal;
		private int geneLength;
		private Annotation window;
		

		public PairedScore(Annotation window,int sampleCount, int inputCount, int sampleTotal, int inputTotal){
			this.sampleCount=sampleCount;
			this.sampleTotal=sampleTotal;
			this.inputCount=inputCount;
			this.inputTotal=inputTotal;
			hasPercentile=false;
			this.window=window;
		}

		public void setGeneInputTotal(int geneInputTotal) {
			this.geneInputTotal=geneInputTotal;
		}

		public void setGeneSampleTotal(int geneSampleTotal) {
			this.geneSampleTotal=geneSampleTotal;
		}

		public String toString(){
			String rtrn=getParentAnnotation().getName()+"\t"+getNormalizedSampleScore()+"\t"+getGlobalEnrichment()+"\t"+sampleCount+"\t"+inputCount+"\t"+getMaxInputCount()+"\t"+getWindowNormEnrichment()+"\t"+getWindowPVal()+"\t"+this.getLocalPVal()+"\t"+getLocalScanPValue(window.size());	
			return rtrn;
		}
		
		

		public void setNumberSampleWindows(int clipWindows) {
			numberSampleWindows=clipWindows;
		}

		public void setNumberInputWindows(int inputWindows) {
			numberInputWindows=inputWindows;
		}
		
		public int getNumberSampleWindows(){
			return this.numberSampleWindows;
		}
		
		public int getNumberInputWindows(){
			return this.numberInputWindows;
		}
		
		public double getWindowNormEnrichment(){
			double numerator=(double)getSampleReadCount()/(double)getMaxInputCount();
			double elution=(double)this.getSampleTotalCount()/(double)this.getNumberSampleWindows();
			double input=(double)this.getInputTotalCount()/(double)this.getNumberInputWindows();
			double denominator=elution/input;
			return numerator/denominator;
		}
		

		public double getScanPValue(int windowSize, int totalSize) {
			return ScanStat.getPValue(getMaxInputCount(), getSampleReadCount(), inputTotal, sampleTotal, windowSize, totalSize);
		}
		
		
		public double getLocalScanPValue(int windowSize) {
			return ScanStat.getPValue(getMaxInputCount(), getSampleReadCount(), this.geneInputTotal, this.geneSampleTotal, windowSize, this.geneLength);
		}

		/**
		 * P values based on local "gene" level
		 * @return window within gene p-value
		 */
		public double getLocalPVal() {
			int k=this.getSampleReadCount();
			int n=getSampleReadCount()+getMaxInputCount();
			
			if(n>0){
				double elution=(double)this.geneSampleTotal;
				double input=(double) this.geneInputTotal;
				if(elution==0 && input==0){return 1.0;}
				double p=elution/(elution+input);
				if(p==1.0){return 1.0;}
				if(p==0.0){return 1.0;}
				Binomial b=new Binomial(n, p);
				return 1-b.cdf(k);
			}
			return 1.0;
			
		}
		
		
		
		/**
		 * P values based on window normalization
		 * @return window normalized p-value
		 */
		public double getWindowPVal() {
			int k=this.getSampleReadCount();
			int n=getSampleReadCount()+getMaxInputCount();
			
			if(n>0){
				double elution=(double)this.getSampleTotalCount()/(double)this.getNumberSampleWindows();
				double input=(double) this.getInputTotalCount()/(double)this.getNumberInputWindows();
				double p=elution/(elution+input);
				Binomial b=new Binomial(n, p);
				return 1-b.cdf(k);
			}
			return 1.0;
			
		}

		/**
		 * 
		 * @return Return the maximum of the input percentile of the value in this window
		 */
		public int getMaxInputCount() {
			int max=Math.max(this.inputPercentile, getInputReadCount());
			return max;
		}

		public double getGlobalEnrichment() {
			return getNormalizedSampleScore()/getPercentileNormalizedInputScore();
		}

		
		private double getPercentileNormalizedInputScore() {
			double num=(double)getMaxInputCount();
			return num/(double)inputTotal;
		}

		public double getNormalizedSampleScore() {
			int numerator=getSampleReadCount();
			return (double)numerator/(double)sampleTotal;
		}

		public void setInputPercentile(int inputPercentile) {
			this.inputPercentile=inputPercentile;
			this.hasPercentile=true;
		}

		public double getPValue() {
			//Return binomial p
			
			int n=getSampleReadCount()+getInputReadCount();
			if(n>0){
				double p=(double)getSampleTotalCount()/((double)getInputTotalCount()+(double)getSampleTotalCount());
				Binomial b=new Binomial(n, p);
				return 1-b.cdf(getSampleReadCount());
			}
			return 1.0;
		}

		private int getInputTotalCount() {
			return this.inputTotal;
		}

		private int getSampleTotalCount() {
			return this.sampleTotal;
		}

		

		public int getSampleReadCount(){
			int count=this.sampleCount;
			if(count==0){count=1;}
			return count;
		}
		
		public int getInputReadCount(){
			int count= this.inputCount;
			if(count==0){count=1;}
			return count;
		}	

		public Annotation getParentAnnotation() {
			return this.parentAnnotation;
		}

		public void setParentAnnotation(Annotation region) {
			this.parentAnnotation=region;
		}
		
	}
	
	private static Collection<? extends Annotation> getRegions(Map<String, Collection<Gene>> regions, String chr, String toSkip){
		Collection<Gene> list=new TreeSet<Gene>();
		if(chr.equalsIgnoreCase("ALL")){
			for(String c: regions.keySet()){
				list.addAll(regions.get(c));
			}
		}
		else if(chr.equalsIgnoreCase("human")){
			for(String c: regions.keySet()){
				if(c.contains("human")){
					list.addAll(regions.get(c));
				}
			}
		}
		else if(chr.equalsIgnoreCase("mouse")){
			for(String c: regions.keySet()){
				if(c.contains("mouse")){
					list.addAll(regions.get(c));
				}
			}
		}
		else{
			list=regions.get(chr);
		}
		
		if(toSkip.isEmpty()){return list;}
		
		Collection<Annotation> list2=new TreeSet<Annotation>();
		
		if(!toSkip.isEmpty()){
			for(Annotation region: list){
				if(!region.getName().equals(toSkip)){
					list2.add(region);
				}
				else{System.err.println("Skipped "+region.getName()+" "+region);}
			}
		}
		
		return list2;
	}
	
	
	
	public static void main(String[] args) throws IOException{
		if(args.length>4){
			System.err.println("Step 1");
			File clipBam=new File(args[0]);
			File inputBam=new File(args[1]);
			String save=args[2];
			
			int windowSize=100;
			
			
			String BEDFile=(args[3]);
			
			String chr="ALL";
			if(args.length>4){
				chr=args[4];
			}
			
			int stepSize=windowSize;
			if(args.length>5){
				stepSize=new Integer(args[5]);
			}
			
			boolean skipIntrons=true;
			if(args.length>6){
				skipIntrons=new Boolean(args[6]);
			}
			
			boolean writeBED=false;
			if(args.length>7){
				writeBED=new Boolean(args[7]);
			}
			
			boolean isMixed=false;
			if(args.length>8){
				isMixed=new Boolean(args[8]);
			}
			
			String toSkip="";
			
			if(args.length>9){
				windowSize=new Integer(args[9]);
			}
			
			double percentile=0.5;
			if(args.length>10){
				percentile=new Double(args[10]);
			}
			
			Map<String, Collection<Gene>> regions=BEDFileIO.loadRegionsFromFileByChr(BEDFile);
			
			System.err.println("loaded regions");
			
			Collection<? extends Annotation> list=getRegions(regions, chr, toSkip);
			
			System.err.println("Get regions");
			
			new CLIPCLAPEnrichmentUpdated(clipBam, inputBam, windowSize, save, list, percentile, stepSize, skipIntrons, writeBED, isMixed);
		}
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=clip BAM file \n arg[1]=input BAM file \n  args[2]=save \n args[3]=BED file of regions \n args[4]=chr to use (optional, mouse, human, or all (all default))\n args[5]=step size (optional, default=windowSize) \n args[6]=skip introns (optional, default=true) \n args[7]=write bedgraph (optional, default=false) \n args[8]= is a mixed human/mouse sample (optional, default=false) \n args[9]=window size (optional, default=100) \n args[10]=percentile (optional, default=0.5)";
	
}
