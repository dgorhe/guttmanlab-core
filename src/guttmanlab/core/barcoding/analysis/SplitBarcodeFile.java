package guttmanlab.core.barcoding.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class SplitBarcodeFile {

	public SplitBarcodeFile(File barcodeFile, int number, String saveDir) throws IOException{
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(barcodeFile)));
		String nextLine;
		int counter=0;
		FileWriter writer=new FileWriter(saveDir+".0.barcode");
		while((nextLine=reader.readLine())!=null && (nextLine.trim().length()>0)){
			if(counter%number ==0){
				writer.close();
				int fileNum=counter/number;
				writer=new FileWriter(saveDir+"."+(fileNum)+".barcode");
				System.err.println(counter+" "+fileNum);
			}
			
			writer.write(nextLine+"\n");
			
			counter++;
		}
		writer.close();
		reader.close();
		
	}
	
	public static void main(String[] args) throws IOException{
		File barcodeFile=new File(args[0]);
		int number=new Integer(args[1]);
		String saveDir=args[2];
		new SplitBarcodeFile(barcodeFile, number, saveDir);
	}
	
}
