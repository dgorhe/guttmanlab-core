package guttmanlab.core.sars;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import guttmanlab.core.datastructures.Pair;

public class AlignmentPipeline {
	
	static String classicalIndex="/central/groups/guttman/data/SARSCov2/bowtieIndex/bowtie2/classicalRNAs";
	static String snoRNAIndex="/central/groups/guttman/data/SARSCov2/bowtieIndex/bowtie2/snoRNAs";
	static String tRNAIndex="/central/groups/guttman/data/SARSCov2/bowtieIndex/bowtie2/tRNAs";
	static String mRNAIndex="/central/groups/guttman/data/SARSCov2/bowtieIndex/bowtie2/mRNAs";

	static String picard="java -jar /groups/guttman/software/picard.2.18.7/picard.jar";
	
	List<String> outputSamFiles;
	
	public AlignmentPipeline(File fastq1, File fastq2, String save) throws IOException, InterruptedException{
		
		this.outputSamFiles=new ArrayList<String>();
		
		String name=fastq1.getName();
		
		String base="/central/groups/guttman/mguttman/temp/"+name+"/";
		new File(base).mkdir();
		
		
		
		
		//Filter low complexity sequences
		
		
		//Align to mRNAs
		Pair<String> unalignedFastq=bowtie(mRNAIndex, fastq1.getAbsolutePath(), fastq2.getAbsolutePath(), save, base, "mRNA");
		
		
		//Align to classical ncRNAs
		unalignedFastq=bowtie(classicalIndex, unalignedFastq.getValue1(), unalignedFastq.getValue2(), save, base, "classicalRNA");
		
		//Align to snoRNAs
		unalignedFastq=bowtie(snoRNAIndex, unalignedFastq.getValue1(), unalignedFastq.getValue2(), save, base, "snoRNA");
		
		//Align to tRNAs
		unalignedFastq=bowtieLocal(tRNAIndex, unalignedFastq.getValue1(), unalignedFastq.getValue2(), save, base, "tRNA");
		
		//Align to introns
		
		//Get remaining unaligned
		bowtieWithUnaligned(classicalIndex, unalignedFastq.getValue1(), unalignedFastq.getValue2(), save+".unaligned.sam");
		
		//Merge BAM file
		mergeSAM(save);
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


	private Pair<String> bowtie(String index, String fq1, String fq2, String output, String baseFolder, String name) throws IOException, InterruptedException{
		File tempFolder=new File(baseFolder+"/"+name+"/");
		tempFolder.mkdir();
		tempFolder.deleteOnExit();
		
		String unalignedPath=tempFolder.getAbsolutePath();
		
		String save=output+"."+name+".sam";
		String cmd="bowtie2 --no-unal --un-conc "+ unalignedPath+" -x "+index +" -1 "+fq1+" -2 "+fq2+" -S "+save;
		System.err.println(cmd);
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		
		this.outputSamFiles.add(save);
		String out1=unalignedPath+"/un-conc-mate.1";
		String out2=unalignedPath+"/un-conc-mate.2";
		return new Pair<String>(out1, out2);
	}
	
	private Pair<String> bowtieLocal(String index, String fq1, String fq2, String output, String baseFolder, String name) throws IOException, InterruptedException{
		File tempFolder=new File(baseFolder+"/"+name+"/");
		tempFolder.mkdir();
		tempFolder.deleteOnExit();
		
		String unalignedPath=tempFolder.getAbsolutePath();
		
		String save=output+"."+name+".sam";
		String cmd="bowtie2 --local --no-unal --un-conc "+ unalignedPath+" -x "+index +" -1 "+fq1+" -2 "+fq2+" -S "+save;
		System.err.println(cmd);
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		
		this.outputSamFiles.add(save);
		String out1=unalignedPath+"/un-conc-mate.1";
		String out2=unalignedPath+"/un-conc-mate.2";
		return new Pair<String>(out1, out2);
	}
	
	private void bowtieWithUnaligned(String index, String fq1, String fq2, String output) throws IOException, InterruptedException{
		String cmd="bowtie2 -x "+index +" -1 "+fq1+" -2 "+fq2+" -S "+output;
		System.err.println(cmd);
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		this.outputSamFiles.add(output);
	}
	
	public static void main(String[] args) throws IOException, InterruptedException{
		File fq1=new File(args[0]);
		File fq2=new File(args[1]);
		String save=args[2];
		new AlignmentPipeline(fq1, fq2, save);
	}
	
}
