package guttmanlab.core.clap.old;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;

public class GetReadPAirs {

	public GetReadPAirs(File samFile1, File samFile2, File fastq1, File fastq2, String save) throws IOException{
		//Get reads
		Collection<String> names=getNames(samFile1);
		names.addAll(getNames(samFile2));
		
		FileWriter writer=new FileWriter(save);
		//Get fastq
		getLines(names, fastq1, writer);
		getLines(names, fastq2, writer);
		writer.close();
		
	}
	
	private void getLines(Collection<String> names, File fastq1, FileWriter writer) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fastq1)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String header=nextLine;
			String line2=reader.readLine();
			String line3=reader.readLine();
			String line4=reader.readLine();
			
			String name=header.split(" ")[0];
			if(names.contains(name)){
				writer.write(header+"\n"+line2+"\n"+line3+"\n"+line4+"\n");
			}
			
		}
		reader.close();
		
	}

	private Collection<String> getNames(File samFile) throws IOException{
		Collection<String> rtrn=new TreeSet<String>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(samFile)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String name="@"+nextLine.split("\t")[0];
			rtrn.add(name);
		}
		reader.close();
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException{
		File sam1=new File(args[0]);
		File sam2=new File(args[1]);
		File f1=new File(args[2]);
		File f2=new File(args[3]);
		String save=args[4];
		new GetReadPAirs(sam1, sam2, f1, f2, save);
	}
}
