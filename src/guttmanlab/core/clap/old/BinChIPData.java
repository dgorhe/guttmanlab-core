package guttmanlab.core.clap.old;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.PairedMappedFragment;
import guttmanlab.core.annotation.PopulatedWindow;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotationcollection.BAMPairedFragmentCollection;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.CloseableIterator;

public class BinChIPData {

	public BinChIPData(File file, int windowSize, String save) throws IOException{
		SAMFileReader inputReader= new SAMFileReader(file);
		
		SAMRecordIterator reads=inputReader.iterator();
		
		Map<SingleInterval, Integer> map=new TreeMap<SingleInterval, Integer>();
		
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			SAMFragment f=new SAMFragment(read);
			SingleInterval interval=f.bin(windowSize);
			int counter=0;
			if(map.containsKey(interval)){counter=map.get(interval);}
			counter++;
			map.put(interval, counter);
		}
		
		reads.close();
		inputReader.close();
		
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: map.keySet()){
			writer.write(region.toBedgraph(map.get(region))+"\n");
		}
		
		writer.close();	
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>2){
			File bamFile=new File(args[0]);
			int windowSize=new Integer(args[1]);
			String save=args[2];
			new BinChIPData(bamFile, windowSize, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam file \n args[1]=window size \n args[2]=save";
}
