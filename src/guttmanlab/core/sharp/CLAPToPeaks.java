package guttmanlab.core.sharp;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import net.sf.samtools.util.CloseableIterator;

public class CLAPToPeaks {

	static String intra="intra";
	static String inter="inter";
	
	public static void convertToPeaks(String file, int minCount, double maxP, double minEnrich, String type) throws IOException {
		List<String> lines=BEDFileIO.loadLines(file,1);
		
		for(String line: lines) {
			filter(line, minCount, maxP, minEnrich, type);
		}
		
	}
	
	
	public static void toBed(String file, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		List<String> lines=BEDFileIO.loadLines(file,1);
		for(String line: lines) {
			Annotation annotation=parse(line);
			writer.write(annotation+"\n");		
		}
		
		writer.close();
	}
	
	public static void toIntraBedGraphs(String file, String save, int minCount, double minPVal, double enrichCutoff, double windowPValCutoff) throws IOException {
		FileWriter pos=new FileWriter(save+".intra.pos.bedgraph");
		FileWriter neg=new FileWriter(save+".intra.neg.bedgraph");
		
		List<String> lines=BEDFileIO.loadLines(file,1);
		for(String line: lines) {
			Annotation annotation=parse(line);
			double intraEnrich=intraEnrichment(line);
			double intraP=intraP(line);
			int intraCount=count(line);
			//window-normalized enrichment and window norm p-val? (column 10 for p-val, column 11 for enrichment; these are from 0 start count)
			double windowEnrichment=windowEnrichment(line);
			double windowP=windowPValue(line);
			
			if(intraCount>minCount && intraP<minPVal && windowEnrichment>enrichCutoff && windowP<windowPValCutoff) {
				if(annotation.getOrientation().equals(Strand.POSITIVE)) {pos.write(annotation.tobedgraph(intraEnrich)+"\n");}
				else {neg.write(annotation.tobedgraph(intraEnrich)+"\n");}
			}
			
		}
		
		pos.close();
		neg.close();
	}
	
	
	
	public static void toIntraBed(String file, String save, int minCount, double minPVal, double enrichCutoff) throws IOException {
		FileWriter writer=new FileWriter(save+".intra.bed");
		
		List<String> lines=BEDFileIO.loadLines(file,1);
		for(String line: lines) {
			Annotation annotation=parse(line);
			double intraEnrich=intraEnrichment(line);
			double intraP=intraP(line);
			int intraCount=count(line);
			
			if(intraCount>minCount && intraP<minPVal && intraEnrich>enrichCutoff) {
				writer.write(annotation.getSingleInterval().toShortBED()+"\n");
			}
			
		}
		
		writer.close();
	}
	
	public static void toInterBedGraphs(String file, String save, int minCount, double minPVal, double enrichCutoff, double windowPValCutoff) throws IOException {
		FileWriter pos=new FileWriter(save+".inter.pos.bedgraph");
		FileWriter neg=new FileWriter(save+".inter.neg.bedgraph");
		
		List<String> lines=BEDFileIO.loadLines(file,1);
		for(String line: lines) {
			Annotation annotation=parse(line);
			double interEnrich=interEnrichment(line);
			double intraP=interP(line);
			int interCount=count(line);
			double windowEnrichment=windowEnrichment(line);
			double windowP=windowPValue(line);
			
			if(interCount>minCount && intraP<minPVal  && windowEnrichment>enrichCutoff && windowP<windowPValCutoff) {
				if(annotation.getOrientation().equals(Strand.POSITIVE)) {pos.write(annotation.tobedgraph(interEnrich)+"\n");}
				else {neg.write(annotation.tobedgraph(interEnrich)+"\n");}
			}
		}
		
		pos.close();
		neg.close();
	}
	
	
	public static void toInterBed(String file, String save, int minCount, double minPVal, double enrichCutoff) throws IOException {
		FileWriter writer=new FileWriter(save+".inter.bed");
		
		
		List<String> lines=BEDFileIO.loadLines(file,1);
		for(String line: lines) {
			Annotation annotation=parse(line);
			double interEnrich=interEnrichment(line);
			double interP=interP(line);
			int interCount=count(line);
			//double windowEnrichment=windowEnrichment(line);
			//double windowP=windowPValue(line);
			
			if(interCount>minCount && interP<minPVal  && interEnrich>enrichCutoff) {
				writer.write(annotation.getSingleInterval().toShortBED()+"\n");
			}
		}
		
		writer.close();
	}
	
	

	public static void writeBedgraphs(String file, String save, int minCount, double minPVal, double enrichCutoff, double windowPValCutoff) throws IOException {
		toInterBedGraphs(file, save, minCount, minPVal, enrichCutoff, windowPValCutoff);
		toIntraBedGraphs(file, save, minCount, minPVal, enrichCutoff, windowPValCutoff);
	}
	
	
	public static void writeBed(String file, String save, int minCount, double minPVal, double enrichCutoff) throws IOException {
		toInterBed(file, save, minCount, minPVal, enrichCutoff);
		toIntraBed(file, save, minCount, minPVal, enrichCutoff);
	}

	//window-normalized enrichment and window norm p-val? (column 10 for p-val, column 11 for enrichment; these are from 0 start count)
	private static double windowPValue(String line) {
		String[] tokens=line.split("\t");
		double count=Double.parseDouble(tokens[10]);
		return count;
	}

	
	private static double windowEnrichment(String line) {
		String[] tokens=line.split("\t");
		double count=Double.parseDouble(tokens[11]);
		return count;
	}
	
	private static int count(String line) {
		String[] tokens=line.split("\t");
		int count=Integer.parseInt(tokens[4]);
		return count;
	}


	private static double intraP(String line) {
		String[] tokens=line.split("\t");
		double p=Double.parseDouble(tokens[8]);
		return p;
	}


	private static double intraEnrichment(String line) {
		String[] tokens=line.split("\t");
		return Double.parseDouble(tokens[6]);
	}
	
	private static double interEnrichment(String line) {
		String[] tokens=line.split("\t");
		return Double.parseDouble(tokens[7]);
	}
	
	private static double sampleCount(String line) {
		String[] tokens=line.split("\t");
		return Double.parseDouble(tokens[5]);
	}
	
	private static double interP(String line) {
		String[] tokens=line.split("\t");
		return Double.parseDouble(tokens[9]);
	}
	
	private static String getType(String line) {
		String[] tokens=line.split("\t");
		return tokens[4];
	}


	private static Annotation parse(String line) {
		String[] tokens=line.split("\t");
		
		Annotation region;
		if(Integer.parseInt(tokens[1])==1){
			region=new SingleInterval(tokens[0]);
			region.setOrientation(Strand.fromString(tokens[2]));
		}
		else {
			region=new BlockedAnnotation(parse(tokens[0], tokens[2]));
		}
		
		return region;
	}


	private static String toIntraBedgraph(String line) {
		String[] tokens=line.split("\t");
		
		Annotation region;
		if(Integer.parseInt(tokens[1])==1){
			region=new SingleInterval(tokens[0]);
			region.setOrientation(Strand.fromString(tokens[2]));
		}
		else {
			region=new BlockedAnnotation(parse(tokens[0], tokens[2]));
		}
		
		return region.tobedgraph(Double.parseDouble(tokens[6]));
	}


	private static void filter(String line, int minCount, double maxP, double minEnrich, String type) {
		String[] tokens=line.split("\t");
		
		Annotation region;
		if(Integer.parseInt(tokens[1])==1){
			region=new SingleInterval(tokens[0]);
			region.setOrientation(Strand.fromString(tokens[2]));
		}
		else {
			region=new BlockedAnnotation(parse(tokens[0], tokens[2]));
		}
		
		region.setName(tokens[3]);
		int count=Integer.parseInt(tokens[4]);
		
		
		double enrich=Double.parseDouble(tokens[6]);
		double p=Double.parseDouble(tokens[8]);
		
		if(type.equals(inter)) {
			enrich=Double.parseDouble(tokens[7]);
			p=Double.parseDouble(tokens[9]);
		}
		
		if(count>minCount && enrich>minEnrich && p<maxP) {
			System.out.println(region.toBED());
		}
		
		
	}

	private static Collection<SingleInterval> parse(String r, String string2) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		String[] tokens=r.split(",");
		for(int i=0; i<tokens.length; i++) {
			SingleInterval interval=new SingleInterval(tokens[i]);
			interval.setOrientation(Strand.fromString(string2));
			rtrn.add(interval);
		}
		return rtrn;
	}
	
	
	public static void toIntronBedgraph(String file, String save) throws IOException {
		FileWriter pos=new FileWriter(save+".pos.bedgraph");
		FileWriter neg=new FileWriter(save+".neg.bedgraph");
		
		List<String> lines=BEDFileIO.loadLines(file);
		System.err.println(lines.size());
		for(String line: lines) {
			Annotation annotation=parse(line);
			double interEnrich=interEnrichment(line);
			String type=getType(line);
			if(type.equals(AssignReads.intron)) {
				if(annotation.getOrientation().equals(Strand.POSITIVE)) {pos.write(annotation.tobedgraph(interEnrich)+"\n");}
				else {neg.write(annotation.tobedgraph(interEnrich)+"\n");}
			}
		}
		
		pos.close();
		neg.close();
	}
	
	public static void toExonBedgraph(String file, String save) throws IOException {
		FileWriter pos=new FileWriter(save+".pos.bedgraph");
		FileWriter neg=new FileWriter(save+".neg.bedgraph");
		
		List<String> lines=BEDFileIO.loadLines(file);
		System.err.println(lines.size());
		for(String line: lines) {
			Annotation annotation=parse(line);
			double interEnrich=interEnrichment(line);
			double p=interP(line);
			String type=getType(line);
			if(type.equals(AssignReads.exon)) {
				if(p<0.0001) {
				Iterator<SingleInterval> iter=annotation.getBlocks();
				while(iter.hasNext()) {
					SingleInterval block=iter.next();
					if(annotation.getOrientation().equals(Strand.POSITIVE)) {pos.write(block.tobedgraph(interEnrich)+"\n");}
					else {neg.write(block.tobedgraph(interEnrich)+"\n");}
				}
			}
			}
		}
		
		pos.close();
		neg.close();
	}
	
	
	public static void toExonCount(String file, String save) throws IOException {
		FileWriter pos=new FileWriter(save+".pos.bedgraph");
		FileWriter neg=new FileWriter(save+".neg.bedgraph");
		
		List<String> lines=BEDFileIO.loadLines(file);
		System.err.println(lines.size());
		for(String line: lines) {
			Annotation annotation=parse(line);
			double count=sampleCount(line);
			String type=getType(line);
			if(type.equals(AssignReads.exon)) {
				Iterator<SingleInterval> iter=annotation.getBlocks();
				while(iter.hasNext()) {
					SingleInterval block=iter.next();
					if(annotation.getOrientation().equals(Strand.POSITIVE)) {pos.write(block.tobedgraph(count)+"\n");}
					else {neg.write(block.tobedgraph(count)+"\n");}
				}
			}
		}
		
		pos.close();
		neg.close();
	}
	
	public static void toExonBed(String file, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		List<String> lines=BEDFileIO.loadLines(file);
		System.err.println(lines.size());
		for(String line: lines) {
			Annotation annotation=parse(line);
			String type=getType(line);
			if(type.equals(AssignReads.exon)) {
				writer.write(annotation.toBED()+"\n");
			}
		}
		
		writer.close();
	}
	
	public static void main(String[] args) throws IOException {
		
		//Can we also add it to do window-normalized enrichment and window norm p-val? (column 10 for p-val, column 11 for enrichment; these are from 0 start count)
		
		if(args.length>4) {
			String file=args[0];
			String save=args[1];
			int minCount=Integer.parseInt(args[2]);
			double minPVal=Double.parseDouble(args[3]);
			double enrichCutoff=Double.parseDouble(args[4]);
			//double normPValCutoff=Double.parseDouble(args[5]);
			
			//writeBedgraphs(file, save, minCount, minPVal, enrichCutoff, normPValCutoff);
			writeBed(file, save, minCount, minPVal, enrichCutoff);
			
		}
		else {System.err.println(usage);}	
	}
	
	

	static String usage=" args[0]=file \n args[1]=save (base name) \n args[2]=min sample count \n args[3]=min p val \n args[4]=enrichment cutoff";
	
	
}
