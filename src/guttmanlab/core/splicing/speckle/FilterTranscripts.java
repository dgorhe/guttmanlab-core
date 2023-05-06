package guttmanlab.core.splicing.speckle;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class FilterTranscripts {

	public FilterTranscripts(File bamFile, File gtfFile) throws IOException {
		Map<SingleInterval, Integer> intronReads=pullSplicedReads(bamFile);
		Collection<Gene> transcripts=BEDFileIO.loadTranscriptsFromGTF(gtfFile);
		
		getTranscriptsWithSupport(transcripts, intronReads);
		
	}
	
	private void getTranscriptsWithSupport(Collection<Gene> transcripts, Map<SingleInterval, Integer> intronReads) {
		//for each gene, get junctions, check if junction has reads, if not, toss junction
		
		for(Gene transcript: transcripts) {
			boolean retain=true;
			Collection<Annotation> junctions=transcript.getIntrons();
			for(Annotation jun: junctions) {
				SingleInterval intron=jun.getSingleInterval();
				if(!intronReads.containsKey(intron)) {retain=false;}
			}
			if(retain) {System.out.println(transcript.toBED());}
		}
		
	}
	
	/*private Map<Gene, Integer> filter(Collection<Gene> allJunctions2, Map<SingleInterval, Integer> intronReads) {
		Map<Gene, Integer> rtrn=new TreeMap<Gene, Integer>();
		for(Gene junction: allJunctions2) {
			SingleInterval intron=junction.getIntrons().iterator().next().getSingleInterval();
			if(intronReads.containsKey(intron)) {
				int counts=intronReads.get(intron);
				rtrn.put(junction, counts);
			}
		}
		return rtrn;
	}*/

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
	
	public static void main(String[] args) throws IOException {
		File bam=new File(args[0]);
		File gtf=new File(args[1]);
		new FilterTranscripts(bam, gtf);
	}
	
}
