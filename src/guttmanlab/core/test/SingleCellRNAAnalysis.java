package guttmanlab.core.test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class SingleCellRNAAnalysis {

	public SingleCellRNAAnalysis(File file) throws IOException{
		
		SamReader reader=SamReaderFactory.makeDefault().open((file));
		SAMRecordIterator iter=reader.iterator();
		
		//Map<String, Collection<SAMRecord>> cellBarcodes=new TreeMap<String, Collection<SAMRecord>>();
		
		Map<String, Integer> cellBarcodes=new TreeMap<String, Integer>();
		while(iter.hasNext()){
			SAMRecord record=iter.next();
			String cellBarcode=(String)record.getAttribute("CB");
			if(cellBarcode!=null){
				//Collection<SAMRecord> list=new ArrayList<SAMRecord>();
				int count=0;
				if(cellBarcodes.containsKey(cellBarcode)){
					count=cellBarcodes.get(cellBarcode);
				}
				count++;
				//list.add(record);
				cellBarcodes.put(cellBarcode, count);
				//System.err.println(cellBarcode);
			}
		}
		
		iter.close();
		
		int count=0;
		Collection<String> barcodes=new TreeSet<String>();
		for(String barcode: cellBarcodes.keySet()){
			System.out.println(barcode+"\t"+cellBarcodes.get(barcode));
			if(cellBarcodes.get(barcode)>10000){
				count++;
				System.err.println(barcode);
				barcodes.add(barcode);
			}
		}
		
		System.err.println(count);
		
		/*iter=reader.iterator();
		
		
		while(iter.hasNext()){
			SAMRecord record=iter.next();
			String cellBarcode=(String)record.getAttribute("CB");
			if(barcodes.contains(cellBarcode)){
				
			}
		}*/
		
		reader.close();
		
		
	}
	
	public static void main (String[] args) throws IOException{
		new SingleCellRNAAnalysis(new File(args[0]));
	}
	
	
}
