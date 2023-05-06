package guttmanlab.core.scSPRITE;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class HiCToSPRITE {

	public HiCToSPRITE(File file, String save, double proportion) throws IOException{
		FileWriter writer=new FileWriter(save);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			double random=Math.random();
			if(random<proportion){
				String[] tokens=nextLine.split("\t");
				String barcode=tokens[0];
				String region1=tokens[1]+":"+tokens[2];
				String region2=tokens[4]+":"+tokens[5];
				writer.write(barcode+"\t"+region1+"\t"+region2+"\n");
			}
			counter++;
			if(counter%10000000 ==0){System.err.println(counter);}
		}
		reader.close();
		writer.close();
	}
	
	public static void main(String[] args) throws IOException{
		File file=new File(args[0]);
		String save=args[1];
		double proportion=new Double(args[2]);
		new HiCToSPRITE(file, save, proportion);
	}
	
}
