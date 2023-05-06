package guttmanlab.core.barcoding.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.SingleInterval;

public class KmersToGeneBED {
	
	int resolution=1000000;

	public KmersToGeneBED(File file, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		int counter=0;
		while((nextLine=reader.readLine())!=null && (nextLine.trim().length()>0)){
			String[] tokens=nextLine.split(",");
			String[] regions=tokens[0].split(" ");
			Collection<SingleInterval> blocks=new TreeSet<SingleInterval>();
			Collection<String> chromosomes=new TreeSet<String>();
			for(int i=0; i<regions.length; i++){
				String chr=regions[i].split("_")[0];
				int start=new Integer(regions[i].split("_")[1]);
				int end=start+resolution;
				SingleInterval region=new SingleInterval(chr, start, end);
				blocks.add(region);
				chromosomes.add(chr);
			}
			
			if(chromosomes.size()==1){
				Annotation spliced=new BlockedAnnotation(blocks, "kmer"+counter);
				writer.write(spliced.toBED()+"\n");
			}
			
			counter++;
			if(counter %10000000 ==0 ){System.err.println(counter);}
		}
		
		writer.close();
		reader.close();
		
	}

	private SingleInterval parse(String line) {
		String chr=line.split(":")[0];
		int start=new Integer(line.split(":")[1].split("-")[0]);
		int end=new Integer(line.split(":")[1].split("-")[1]);
		return new SingleInterval(chr, start, end);
	}
	
	public static void main(String[] args)throws IOException{
		File file=new File(args[0]);
		String save=args[1];
		new KmersToGeneBED(file, save);
	}
	
}
