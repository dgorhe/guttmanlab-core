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
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;

public class FilterSameChromosome {

	public FilterSameChromosome(File barcodeFile) throws IOException{
		//go through all clusters and split by homogeneity
		FileWriter writer=new FileWriter("C:/Data/Barcoding/nonHomo.txt");
		FileWriter writerHomo=new FileWriter("C:/Data/Barcoding/Homo.txt");

		
		int countHomo=0;
		int homoNonZero=0;
		int nonHomo=0;
		int total=0;
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(barcodeFile)));
		String nextLine;
		int counter=0;
		while((nextLine=reader.readLine())!=null && (nextLine.trim().length()>0)){
			String[] tokens=nextLine.split("\t");
			
			
			//if(tokens.length>2 && tokens.length<11){
			
			String barcode=tokens[0];
			
			
			boolean homo=isHomogeneous(nextLine, tokens);
			if(homo){	
				countHomo++;
				if(tokens.length>2){
					homoNonZero++;
					writerHomo.write(nextLine+"\n");
				}
			}
			else{
				nonHomo++;
				writer.write(nextLine+"\n");
			}
			total++;
	
			if(counter%1000000 ==0){System.err.println("Read in "+counter +" clusters");}
			counter++;
			
		//}
		
			
		}
		reader.close();
		writer.close();
		writerHomo.close();
		System.err.println(countHomo +" "+total+" "+((double)countHomo/(double)total)+" "+((double)homoNonZero)/((double)homoNonZero+nonHomo));

	}

	private boolean isHomogeneous(String nextLine, String[] tokens) {
		//if regions in cluster are on different chromosomes
		Collection<String> chr=new TreeSet<String>();
		
		for(int i=1; i<tokens.length; i++){
			chr.add(tokens[i].split(":")[0]);
			
		}
		
		if(chr.size() == 1){return true;}
		return false;
	}
	
	public static void main(String[] args) throws IOException{
		new FilterSameChromosome(new File(args[0]));
	}
}
