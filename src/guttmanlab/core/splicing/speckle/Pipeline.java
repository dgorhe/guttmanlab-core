package guttmanlab.core.splicing.speckle;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;

public class Pipeline {

	public static void main(String[] args) throws IOException {
		if(args.length>5) {
			List<File> bamFiles=parseBAM(args[0]);
			GTFToJunctions gtf=new GTFToJunctions(new File(args[1]));
			String saveDir=args[2];
			int numPerm=Integer.parseInt(args[3]);
			
		
			
			ComputeGeneExpression geneExpression=new ComputeGeneExpression(new File(args[4]), gtf);
			
			double rpkmFilter=Double.parseDouble(args[5]);
			
			geneExpression.writeRPKM(saveDir+"/genes.expression");
			Map<String, Double> rpkm=geneExpression.getRPKM(rpkmFilter);
			
			writeAnalyzedGenes(saveDir+"/analyzedGenes.bed", rpkm, gtf);
			
			List<String> splicingFiles=new ArrayList<String>();
			for(File bam: bamFiles) {
				System.err.println(bam.getName());
				String save=saveDir+"/"+bam.getName()+".splicing";
				new SplicingRatiosPerIsoform(bam, gtf.getGenes(), save, false);
				splicingFiles.add(save);
			}
			
			/*String save=saveDir+"/junctionCounts.matrix";
			Compute4PLFit.writeJunctionCounts(splicingFiles, save, rpkm);
			*/
					
			String save=saveDir+"/junctionFraction.matrix";
			Compute4PLFit.writeJunctionFraction(splicingFiles, save, rpkm); //TODO Also filter junctions overlapping that come from other gene
			
			String observed=save;
			
			Map<SingleInterval, Double> uniqueJunctions=computeUniqueJunctions(observed);
			
			List<String> permFiles=new ArrayList<String>();
			for(int i=0; i<numPerm; i++) {
				System.err.println("processing perm"+i);
				save=saveDir+"/perm"+i+".matrix";
				Compute4PLFit.processPerm(splicingFiles, save, i, rpkm);
				permFiles.add(save);
			}
			
			writeSummary(saveDir+"/results.matrix", observed, permFiles, uniqueJunctions, gtf);
			
			System.exit(0);
		}
		else {System.err.println(usage);}
	}

	

	private static void writeAnalyzedGenes(String save, Map<String, Double> rpkm, GTFToJunctions gtf) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String gene: rpkm.keySet()) {
			SingleInterval region=gtf.getGeneCoordinates(gene);
			writer.write(region.toBED()+"\n");
		}
		
		writer.close();
	}



	private static Map<SingleInterval, Double> computeUniqueJunctions(String observed) throws IOException {
		Map<SingleInterval, Integer> junctionCounts=parseJunctionCounts(observed);
		Map<String, Collection<SingleInterval>> junctionsByGene=parseJunctions(observed);
		
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		for(String gene: junctionsByGene.keySet()) {
			Collection<SingleInterval> junctions=junctionsByGene.get(gene);
			for(SingleInterval junction: junctions) {
				double overlaps=fractionOfOverlappers(junction, junctions, junctionCounts);
				rtrn.put(junction, overlaps);
			}
		}
		
		return rtrn;
	}



	private static Map<String, Collection<SingleInterval>> parseJunctions(String file) throws IOException {
		Map<String, Collection<SingleInterval>> rtrn=new TreeMap<String, Collection<SingleInterval>>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String name=tokens[0];
			SingleInterval junction=new SingleInterval(tokens[1]);
			if(!rtrn.containsKey(name)) {rtrn.put(name, new TreeSet<SingleInterval>());}
			Collection<SingleInterval> list=rtrn.get(name);
			list.add(junction);
		}
		
		reader.close();
		return rtrn;
	}
	
	private static Map<SingleInterval, Integer> parseJunctionCounts(String file) throws IOException {
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			//String name=tokens[0];
			SingleInterval junction=new SingleInterval(tokens[1]);
			Integer count=Integer.parseInt(tokens[3]);
			rtrn.put(junction, count);
		}
		
		reader.close();
		return rtrn;
	}


	private static double fractionOfOverlappers(SingleInterval junction1, Collection<SingleInterval> junctions, Map<SingleInterval, Integer> junctionCounts) {
		double o=junctionCounts.get(junction1);
		
		//Collection<SingleInterval> overlappers=gtf.getOverlappingJunctions(junction1);
		
		
		Collection<SingleInterval> overlappers=new TreeSet<SingleInterval>();
		
		for(SingleInterval junction2: junctions) {
			if(!junction1.equals(junction2)) {
				if(junction1.overlaps(junction2) || junction2.overlaps(junction1)) {
					overlappers.add(junction2);
				}
			}
		}
		
		double sum=o;
		for(SingleInterval val: overlappers) {
			sum+=junctionCounts.get(val);
		}
		
		
		return o/sum;
	}
	

	private static boolean overlaps(SingleInterval junction1, Collection<SingleInterval> junctions) {
		boolean rtrn=false;
		for(SingleInterval junction2: junctions) {
			if(!junction1.equals(junction2)) {
				if(junction1.overlaps(junction2) || junction2.overlaps(junction1)) {rtrn=true;}
			}
		}
		return rtrn;
	}



	private static void writeSummary(String save, String observed, List<String> permFiles, Map<SingleInterval, Double> uniqueJunctions, GTFToJunctions gtf) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		//+rpkm+"\t"+junctionCount+"\t"+junctionFraction+"\t"+fit[2]
		writer.write("gene name\tgene type\tjunction coordinates\tfraction of overlapping\tRPKM (0min)\tnumber of junction reads\tfraction of junction reads relative to max of gene\tobserved C50");
		for(int i=0; i<permFiles.size(); i++) {writer.write("\tperm"+i+" C50");}
		writer.write("\n");
		
		Map<String, List<Double>> scores=new TreeMap<String, List<Double>>();
		
		for(String file: permFiles) {
			Map<String, Double> vals=parse(file);
			for(String gene: vals.keySet()) {
				if(!scores.containsKey(gene)) {scores.put(gene, new ArrayList<Double>());}
				List<Double> list=scores.get(gene);
				if(vals.containsKey(gene)) {
					list.add(vals.get(gene));
				}
			}
		}
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(observed)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String name=tokens[0]+"_"+tokens[1];
			SingleInterval region=new SingleInterval(tokens[1]);
			double unique=uniqueJunctions.get(region);
			String geneType=gtf.getGeneType(tokens[0]);
			
			writer.write(tokens[0]+"\t"+geneType+"\t"+tokens[1]+"\t"+unique+"\t"+tokens[2]+"\t"+tokens[3]+"\t"+tokens[4]+"\t"+tokens[5]);
			
			List<Double> perms=scores.get(name);
			for(Double perm: perms) {writer.write("\t"+perm);}
			writer.write("\n");
		}
		
		reader.close();
		
		writer.close();
	}



	private static Map<String, Double> parse(String file) throws NumberFormatException, IOException {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String name=tokens[0]+"_"+tokens[1];
			double val=Double.parseDouble(tokens[5]);
			rtrn.put(name, val);
		}
		
		reader.close();
		return rtrn;
	}



	private static List<File> parseBAM(String string) {
		File[] files=new File(string).listFiles();
		
		List<File> rtrn=new ArrayList<File>();
		
		for(int i=0; i<files.length; i++) {
			if(files[i].getName().endsWith(".bam")) {rtrn.add(files[i]);}
		}
		
		return rtrn;
	}
	
	static String usage="-Djava.awt.headless=true\n args[0]=files \n args[1]=gtf \n args[2]=save dir \n args[3]=num perm \n args[4]=0min bam \n args[5]=min rpkm";
}
