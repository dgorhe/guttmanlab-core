package guttmanlab.core.nanobody;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;

public class ProcessFiles {

	public static void main(String[] args) throws IOException {
		File file=new File(args[0]);
		int index1=Integer.parseInt(args[1]);
		int index2=Integer.parseInt(args[2]);
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		int counter=0;
		
		while ((nextLine = reader.readLine()) != null) {
			if(counter!=0) {
				String[] tokens=nextLine.split(",");
				String sequence=tokens[1];
				
				Double val1=Double.parseDouble(tokens[index1]);
				
				
				Double val2=Double.parseDouble(tokens[index2]);
				
				
				System.out.println(sequence+"\t"+val1.intValue()+"\t"+val2.intValue());
				
			}
			counter++;
			if(counter%1000000==0) {System.err.println(counter);}
		}
		reader.close();
	}
	
	
}
