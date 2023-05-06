package guttmanlab.core.util;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.SAMFragment;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class JunctionCounts {

	public JunctionCounts(File bam, String save) throws IOException {
		Map<Annotation, Integer> intronCounts=new TreeMap<Annotation, Integer>();
		
		SAMFileReader reader=new SAMFileReader(bam);
		SAMRecordIterator reads=reader.iterator();
		
		int counter=0;
		while(reads.hasNext()) {
			SAMRecord read=reads.next();
			SAMFragment frag=new SAMFragment(read);
			Collection<Annotation> introns=frag.getIntrons();
			for(Annotation intron: introns) {
				int count=0;
				if(intronCounts.containsKey(intron)) {count=intronCounts.get(intron);}
				count++;
				intronCounts.put(intron, count);
			}
			counter++;
			if(counter%100000 ==0) {System.err.println(counter);}
		}
		
		reads.close();
		reader.close();
		
		write(save, intronCounts);
		
	}

	private void write(String save, Map<Annotation, Integer> intronCounts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Annotation intron: intronCounts.keySet()) {
			int score=intronCounts.get(intron);
			writer.write(intron.getReferenceName()+"\t"+intron.getReferenceStartPosition()+"\t"+intron.getReferenceEndPosition()+"\t"+score+"\n");
		}
		
		writer.close();
	}
	
	public static void main(String[] args) throws IOException {
		if(args.length>1) {
		File bam=new File(args[0]);
		String save=args[1];
		new JunctionCounts(bam,save);
		}
		else {System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam file \n args[1]=output";
	
}
