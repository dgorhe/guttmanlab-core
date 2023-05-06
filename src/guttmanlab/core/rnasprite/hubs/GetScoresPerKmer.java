package guttmanlab.core.rnasprite.hubs;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.rnasprite.Kmer;

public class GetScoresPerKmer {

	public GetScoresPerKmer(File[] files, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		for(File file: files) {
			String score=parseMaxScore(file);
			if(!score.isEmpty()) {
				writer.write(score+"\n");
			}
		}
		writer.close();
	}

	private String parseMaxScore(File file) throws IOException {
		List<String> lines=BEDFileIO.loadLines(file.getAbsolutePath());
		
		double max=-1;
		String maxLine="";
		for(String line: lines) {
			String[] tokens=line.split("\t");
			double randomVal=Double.parseDouble(tokens[3]);
			if(randomVal>max) {
				max=randomVal;
				maxLine=line;
			}
		}
		if(max>-1) {
			maxLine=fix(maxLine);
		}
		return maxLine;
		
	}

	private String fix(String maxLine) {
		String[] tokens=maxLine.split("\t");
		Kmer kmer=new Kmer(tokens[0]);
		
		String rtrn=maxLine+"\t"+kmer.getSize();
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		if(args.length>1) {
			File[] files=new File(args[0]).listFiles();
			String save=args[1];
			new GetScoresPerKmer(files, save);
		}
		else {System.err.println(usage);}
	}
	static String usage=" args[0]=files \n args[1]=save";
}
