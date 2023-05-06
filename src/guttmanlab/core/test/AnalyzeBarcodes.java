package guttmanlab.core.test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import com.sleepycat.je.rep.util.ldiff.Record;

import guttmanlab.core.sequence.Sequence;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class AnalyzeBarcodes {

	private Sequence primer1=new Sequence("cctgcaaaggccatgctata");
	private Sequence primer2=new Sequence("atatccaaaaccgctccttc");
	
	public AnalyzeBarcodes(File bamFile, String save) throws IOException{
		SamReader reader=SamReaderFactory.makeDefault().open((bamFile));
		SAMRecordIterator iter=reader.iterator();
		
		Map<String, Map<String, Collection<String>>> splitCollection=new TreeMap<String, Map<String, Collection<String>>>();
		
		Map<String, Integer> recordCount=new TreeMap<String, Integer>();
		
		//System.err.println(primer1.reverseComplement().toString().toUpperCase());
		
		int counter=0;
		while(iter.hasNext()){
			SAMRecord record=iter.next();
			
			boolean hasPrimer1=record.getReadString().contains(primer1.reverseComplement().toString().toUpperCase());
			boolean hasPrimer2=(record.getReadString().contains(primer2.reverseComplement().toString().toUpperCase()));
			
			if(hasPrimer1 || hasPrimer2){
			
				String name=record.getReadName();
				String barcode=name.split("::")[1];
				ArrayList<String> barcodes=split(barcode);
				String sampleID=barcodes.remove(barcodes.size()-1);
				String remainder=makeBarcode(barcodes);
				
				int count=0;
				
				String recordString=record.getReferenceName()+"_"+record.getAlignmentStart();
				if(hasPrimer1){ recordString="primer1";}
				if(hasPrimer2){recordString="primer2";}
				if(hasPrimer1 && hasPrimer2){recordString="AHHH";}
				
				if(recordCount.containsKey(recordString)){count=recordCount.get(recordString);}
				count++;
				recordCount.put(recordString, count);
				
				add(sampleID, remainder, recordString, splitCollection);
			}
			
			
			counter++;
			
			if(counter%100000 ==0){System.err.println(counter);}
		}
		iter.close();
		
		write(save, splitCollection);
		collisions(splitCollection);
		
		
	}
	
	private void collisions(Map<String, Map<String, Collection<String>>> splitCollection) {
		Map<String, Collection<String>> allBarcodes=new TreeMap<String, Collection<String>>();
		
		for(String sample: splitCollection.keySet()){
			for(String barcode: splitCollection.get(sample).keySet()){
				Collection<String> list=new TreeSet<String>();
				if(allBarcodes.containsKey(barcode)){
					list=allBarcodes.get(barcode);
				}
				list.add(sample);
				allBarcodes.put(barcode, list);
			}
		}
		
		int count=0;
		for(String barcode: allBarcodes.keySet()){
			if(allBarcodes.get(barcode).size()>1){
				//System.err.println(barcode+" "+allBarcodes.get(barcode));
				count++;
			}
		}
		
		System.err.println(count+" "+allBarcodes.size());
	}

	private void write(String save, Map<String, Map<String, Collection<String>>> splitCollection) throws IOException {
		for(String sample: splitCollection.keySet()){
			writeSub(save+"."+sample, splitCollection.get(sample));
		}
		
	}

	private void writeSub(String string, Map<String, Collection<String>> map) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		int countPaired=0;
		int countUnpaired=0;
		
		for(String barcode: map.keySet()){
			Collection<String> records=map.get(barcode);
				if(records.size()==1){countUnpaired++;}
				if(records.size()==2){countPaired++;}
				writer.write(barcode);
				for(String record: records){
					writer.write("\t"+record);
				}
				writer.write("\n");
			
		}
		
		double percent=(double)countPaired/(double)(countPaired+countUnpaired);
		System.err.println(string+" "+countPaired+" "+countUnpaired+" "+percent);
		
		writer.close();
	}

	private void add(String sampleID, String remainder, String record, Map<String, Map<String, Collection<String>>> splitCollection) {
		Map<String, Collection<String>> map=new TreeMap<String, Collection<String>>();
		if(splitCollection.containsKey(sampleID)){
			map=splitCollection.get(sampleID);
		}
		
		Collection<String> records=new TreeSet<String>();
		if(map.containsKey(remainder)){
			records=map.get(remainder);
		}
		
		records.add(record);
		
		map.put(remainder, records);
		splitCollection.put(sampleID, map);
		
	}

	private ArrayList<String> split(String barcode) {
		ArrayList<String> rtrn=new ArrayList<String>();
		String[] tokens=barcode.split("\\[");
		//System.err.println("\n"+barcode);
		for(int i=1; i<tokens.length; i++){
			//System.err.println(i+" "+tokens[i]);
			String newToken=tokens[i].replaceAll("\\]", "");
			rtrn.add(newToken);
		}
		return rtrn;
	}

	private String makeBarcode(ArrayList<String> barcodes) {
		String rtrn="";
		for(String barcode: barcodes){
			rtrn+="["+barcode+"]";
		}
		return rtrn;
	}

	public static void main(String[] args) throws IOException{
		File bamFile=new File(args[0]);
		String save=args[1];
		new AnalyzeBarcodes(bamFile, save);
	}
	
}
