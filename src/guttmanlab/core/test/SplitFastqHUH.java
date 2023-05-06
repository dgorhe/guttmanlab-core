package guttmanlab.core.test;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.io.BEDFileIO;
import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;

public class SplitFastqHUH {

	public SplitFastqHUH(File fastq, Map<String, String> dpmToName, Map<String, String> antibodyNames, String save){
		FastqReader reader=new FastqReader(fastq);
		FastqWriter writerAb=new BasicFastqWriter(new File(save+".HUH.fq"));
		FastqWriter writerDpm=new BasicFastqWriter(new File(save+".DPM.fq"));
		while(reader.hasNext()){
			FastqRecord read=reader.next();
			String name=parse(read, dpmToName, antibodyNames);
			if(name!=null){
				FastqWriter writer=writerAb;
				if(name.contains("DPM")){writer=writerDpm;}
				String newName=read.getReadHeader()+"["+name+"]";
				FastqRecord newRead=new FastqRecord(newName, read.getReadString(), read.getBaseQualityHeader(), read.getBaseQualityString());
				writer.write(newRead);
			}
			
			//Append to name
		}
		reader.close();
		writerDpm.close();
		writerAb.close();
		
	}

	private String parse(FastqRecord read, Map<String, String> dpmToName, Map<String, String> antibodyNames) {
		//DPM
		//Check if read has DPM
		String dpm=hasDPM(read, dpmToName);
		
		if(dpm!=null){return dpm;}
		
		//UMI (8)+ HUH (35)
		return parseHUH(read, antibodyNames);
		
		
	}
	
	private String parseHUH(FastqRecord read, Map<String, String> antibodyNames) {
		String odd=read.getReadString().substring(52, 69);
		//System.err.println("HUH: "+odd+" "+read.getReadString());
		
		if(antibodyNames.containsKey(odd)){return antibodyNames.get(odd);}
		return null;
	}

	private String hasDPM(FastqRecord read, Map<String, String> dpmToName) {
		String start=read.getReadString().substring(0, 10);
		//System.err.println("DPM: "+start+" "+read.getReadString());
		
		if(dpmToName.containsKey(start)){return dpmToName.get(start);}
		return null;
	}

	private static Map<String, String> parseSeqToName(String string) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		List<String> lines=BEDFileIO.loadLines(string);
		
		for(String line: lines){
			String[] tokens=line.split("\t");
			//System.err.println(line+"\t"+tokens.length+" "+tokens[0]);
			rtrn.put(tokens[1], tokens[0]);
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>3){
			File fastq=new File(args[0]);
			Map<String, String> dpm=parseSeqToName(args[1]);
			Map<String, String> antibody=parseSeqToName(args[2]);
			String save=args[3];
			new SplitFastqHUH(fastq, dpm, antibody, save);
		}
		else{System.err.println(usage);}
	}
	
	

	static String usage=" args[0]=fastq \n args[1]=dpm \n args[2]=antibody \n args[3]=save";
}
