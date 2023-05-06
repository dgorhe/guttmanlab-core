package guttmanlab.core.sharp;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import guttmanlab.core.datastructures.Pair;

public class TieredAlignmentPipeline {
	
	static String picard="java -jar /groups/guttman/software/picard.2.18.7/picard.jar";
	List<String> outputSamFiles;
	static String index45S="/central/groups/guttman/mguttman/bowtieIndex/mouse/bowtieIndex/45S.mouse";
	static String indexNC="/central/groups/guttman/mguttman/bowtieIndex/mouse/bowtieIndex/smallRNA.mouse";
	static String indexRepeats="/central/groups/guttman/mguttman/bowtieIndex/mouse/bowtieIndex/allRepeats.mouse";
	static String indexGenome="/central/groups/guttman/genomes/mus_musculus/GRCm38/bowtie2/GRCm38.primary_assembly.genome";
	static String indexMRNA="/central/groups/guttman/mguttman/bowtieIndex/mouse/bowtieIndex/mm10.refGenes";
	static String indexIntron="/central/groups/guttman/mguttman/bowtieIndex/mouse/bowtieIndex/introns.mouse";
	
	
	public TieredAlignmentPipeline(File fastq1, File fastq2, String save) throws IOException, InterruptedException {
		//loadBowtie();
		
		this.outputSamFiles=new ArrayList<String>();
		String name=fastq1.getName();
		String base="/central/groups/guttman/mguttman/temp/"+name+"/";
		new File(base).mkdir();
		
		
		Pair<String> unalignedFastq=bowtie(index45S, fastq1.getAbsolutePath(), fastq2.getAbsolutePath(), save, base, "45S", false);
		unalignedFastq=bowtie(indexRepeats, unalignedFastq.getValue1(), unalignedFastq.getValue2(), save, base, "repeats", true);
		
		
		unalignedFastq=bowtie(indexMRNA, unalignedFastq.getValue1(), unalignedFastq.getValue2(), save, base, "mRNA", true);
		unalignedFastq=bowtie(indexGenome, unalignedFastq.getValue1(), unalignedFastq.getValue2(), save, base, "genome", true);
		
		
		//Align to 45S
		/*Pair<String> unalignedFastq=bowtie(index45S, fastq1.getAbsolutePath(), fastq2.getAbsolutePath(), save, base, "45S");
		
		//Align to ncRNAs
		unalignedFastq=bowtie(indexNC, unalignedFastq.getValue1(), unalignedFastq.getValue2(), save, base, "classicalRNA");
				
		//Align to repeats
		unalignedFastq=bowtie(indexRepeats, unalignedFastq.getValue1(), unalignedFastq.getValue2(), save, base, "repeats");
				
		//Align to mRNA
		unalignedFastq=bowtie(indexMRNA, unalignedFastq.getValue1(), unalignedFastq.getValue2(), save, base, "mRNA");
			
		//Align to intron
		//unalignedFastq=bowtie(indexIntron, unalignedFastq.getValue1(), unalignedFastq.getValue2(), save, base, "intron");
		
		//Align to genome
		unalignedFastq=bowtie(indexGenome, unalignedFastq.getValue1(), unalignedFastq.getValue2(), save, base, "genome");*/
				
		
		//mergeSAM(save);
	}
	
	private void loadBowtie() throws IOException, InterruptedException {
		String cmd="module add bowtie/2.3.4.1";
		System.err.println(cmd);
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		System.err.print(p.getErrorStream());
		
	}

	private void mergeSAM(String save) throws IOException, InterruptedException {
		String cmd=picard+" MergeSamFiles";
		
		for(String file: this.outputSamFiles){
			cmd+=" I="+file;
		}
		
		cmd+=" O="+save+".merged.bam";
		cmd+=" SO=coordinate MSD=true";
		System.err.println(cmd);
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		System.err.print(p.getErrorStream());
	}
	
	private Pair<String> bowtie(String index, String fq1, String fq2, String output, String baseFolder, String name, boolean sort) throws IOException, InterruptedException{
		File tempFolder=new File(baseFolder+"/"+name+"/");
		tempFolder.mkdir();
		tempFolder.deleteOnExit();
		
		String unalignedPath=tempFolder.getAbsolutePath();
		
		String save=output+"."+name+".sam";
		String cmd="bowtie2 --no-unal --no-mixed --un-conc "+ unalignedPath+" -x "+index +" -1 "+fq1+" -2 "+fq2+" -S "+save;
		System.err.println(cmd);
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		
		if(sort) {sortAndIndex(save);}
		
		this.outputSamFiles.add(save);
		String out1=unalignedPath+"/un-conc-mate.1";
		String out2=unalignedPath+"/un-conc-mate.2";
		return new Pair<String>(out1, out2);
	}
	

	private void sortAndIndex(String save) throws IOException, InterruptedException {
		String out=sort(save);
		out=removeDuplicates(out);
		index(out);
	}

	private String removeDuplicates(String input) throws IOException, InterruptedException {
		String out=input+".duplicates.bam";
		// I= O= M= REMOVE_DUPLICATES=true
		String cmd=picard+" MarkDuplicates REMOVE_DUPLICATES=true";
		cmd+=" I="+input;
		cmd+=" O="+out;
		cmd+=" M="+input+".metrics";
		System.err.println(cmd);
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		System.err.print(p.getErrorStream());
		return out;
	}

	private void index(String input) throws IOException, InterruptedException {
		String cmd=picard+" BuildBamIndex";
		cmd+=" I="+input;
		System.err.println(cmd);
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		System.err.print(p.getErrorStream());	
	}

	private String sort(String input) throws IOException, InterruptedException {
		String cmd=picard+" SortSam";
		
		cmd+=" I="+input;
		
		String output=input+".sorted.bam";
		
		cmd+=" O="+output;
		cmd+=" SO=coordinate";
		System.err.println(cmd);
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		System.err.print(p.getErrorStream());
		return output;
	}

	public static void main(String[] args) throws IOException, InterruptedException {
		File fastq1=new File(args[0]);
		File fastq2=new File(args[1]);
		String save=args[2];
		new TieredAlignmentPipeline(fastq1, fastq2, save);
		
	}
	
	
}
