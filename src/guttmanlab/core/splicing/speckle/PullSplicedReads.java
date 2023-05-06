package guttmanlab.core.splicing.speckle;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class PullSplicedReads {

	private static void pullSplicedReads(File[] files, String save) {
		SAMFileReader reader=new SAMFileReader(files[0]);
		SAMFileWriter writer=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(reader.getFileHeader(), false, new File(save));
		reader.close();
		
		
		for(int i=0; i<files.length; i++) {
			System.err.println(files[i].getName());
			pullSplicedReads(files[i], writer);
		}
		
		writer.close();
		
	}
	
	private static Map<SingleInterval, Integer> pullSplicedReads(File bamFile) {
		SAMFileReader reader=new SAMFileReader(bamFile);
		SAMRecordIterator reads=reader.iterator();
		
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			
			if(SAMFragment.isSpliced(record)) {
				
				SAMFragment frag=new SAMFragment(SAMFragment.compressCIGAR(record));
				Collection<Annotation> introns=frag.getIntrons();
				for(Annotation intron: introns) {
					int count=0;
					SingleInterval intronSI=intron.getSingleInterval();
					if(rtrn.containsKey(intronSI)) {
						count=rtrn.get(intronSI);
					}
					count++;
					rtrn.put(intronSI, count);
				}
			}
				
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		reader.close();
		reads.close();
		
		
		return rtrn;
		
	}
	
	
	
	private static void pullSplicedReads(File bamFile, SAMFileWriter writer) {
		SAMFileReader reader=new SAMFileReader(bamFile);
		SAMRecordIterator reads=reader.iterator();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			
			if(SAMFragment.isSpliced(record)) {writer.addAlignment(record);}
				
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		reader.close();
		reads.close();
		
	}

	private static void write(String save, Set<SingleInterval> introns) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval intron: introns) {
			writer.write(intron.toBED()+"\n");
		}
		
		writer.close();
	}
	
	
	private static void write(String save, Map<SingleInterval, Integer> intronReads) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval intron: intronReads.keySet()) {
			writer.write(intron.toBedgraph(intronReads.get(intron))+"\n");
		}
		
		writer.close();
		
	}


	public static void main(String[] args) throws IOException {
		if(args.length>1) {
		File file=new File(args[0]);
		String save=args[1];
		
		Collection<Gene> junctions=BEDFileIO.loadRegionsFromFile(args[2]);
		
		Map<SingleInterval, Integer> intronReads=pullSplicedReads(file);
		
		write(save+".all.reads.bed", intronReads.keySet());
		
		FileWriter writer=new FileWriter(save);
		for(Gene junction: junctions) {
			SingleInterval intron=junction.getIntrons().iterator().next().getSingleInterval();
			if(intronReads.containsKey(intron)) {
				int counts=intronReads.get(intron);
				writer.write(junction.toBED("c="+counts)+"\n");
			}
		}
		writer.close();
		
		//write(save+".counts.bedgraph", intronReads);
		
		} else {System.err.println(usage);}
	}

	

	



	static String usage=" java -jar -Dsnappy.disable=true\n args[0]=bam files \n args[1]=save \n args[2]=junctions";
}
