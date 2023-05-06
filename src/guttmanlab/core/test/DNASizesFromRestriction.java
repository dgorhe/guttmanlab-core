package guttmanlab.core.test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;

public class DNASizesFromRestriction {

	public DNASizesFromRestriction(String fastaDir, Map<String, Sequence> restrictionEnzymes, String save) throws IOException{
		
		File[] fastaFiles=new File(fastaDir).listFiles();
		
		Collection<SingleInterval> fragments=new TreeSet<SingleInterval>();
		
		for(int i=0; i<1; i++){
		//for(int i=0; i<fastaFiles.length; i++){
			System.err.println(fastaFiles[i].getName());
			
			Sequence seq=FastaFileIOImpl.readFromFile(fastaFiles[i].getAbsolutePath()).iterator().next();

			System.err.println(seq.getName());
			String chr=seq.getName();
			
			fragments.addAll(getFragmentSizes(seq, restrictionEnzymes, chr));
		}
		
		write(save, fragments);
		
	}

	private void write(String save, Collection<SingleInterval> allFragments) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval frag: allFragments){
			writer.write(frag.getReferenceName()+"\t"+frag.getReferenceStartPosition()+"\t"+frag.getReferenceEndPosition()+"\t"+frag.getGenomicLength()+"\n");
		}
		
		writer.close();
	}

	private Collection<SingleInterval> getFragmentSizes(Sequence seq, Map<String, Sequence> restrictionEnzymes, String chr) {
		Map<String, Collection<SingleInterval>> regions=find(seq, restrictionEnzymes);
		return makeIntervals(regions);
		
	}

	private Collection<SingleInterval> makeIntervals(Map<String, Collection<SingleInterval>> regions) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		Collection<SingleInterval> temp=new TreeSet<SingleInterval>();
		for(String name: regions.keySet()){
			temp.addAll(regions.get(name));
		}
		
		Iterator<SingleInterval> iter=temp.iterator();
		if(!iter.hasNext()){return rtrn;}
		SingleInterval current=iter.next();
		while(iter.hasNext()){
			SingleInterval next=iter.next();

			//TODO Check to make sure chromosomes are the same
			if(current.getReferenceName().equals(next.getReferenceName())){
				SingleInterval interval=new SingleInterval(current.getReferenceName(), current.getReferenceEndPosition(), next.getReferenceStartPosition());
				//System.err.println(interval.toUCSC());
				rtrn.add(interval);
			}
			current=next;
			
		}
		
		return rtrn;
	}

	private Map<String, Collection<SingleInterval>> find(Sequence seq, Map<String, Sequence> restrictionEnzymes) {
		Map<String, Collection<SingleInterval>> rtrn=new TreeMap<String, Collection<SingleInterval>>();
		for(String re: restrictionEnzymes.keySet()){
			Sequence reSeq=restrictionEnzymes.get(re);
			Collection<SingleInterval> hits=seq.find(reSeq);
			System.err.println(reSeq.getName()+" "+reSeq.getSequenceBases()+" "+hits.size());
			rtrn.put(re, hits);
		}
		return rtrn;
	}
	
	private static Map<String, Sequence> parse(String file) throws IOException {
		Map<String, Sequence> rtrn=new TreeMap<String, Sequence>();
		Collection<String> lines=BEDFileIO.loadLines(file);
		for(String line: lines){
			String[] tokens=line.split("\t");
			Sequence seq=new Sequence(tokens[0], tokens[1]);
			rtrn.put(seq.getName(), seq);
		}
		return rtrn;
	}

	public static void main(String[] args) throws IOException{
		String fastaDir=args[0];
		Map<String, Sequence> re=parse(args[1]);
		String save=args[2];
		
		//for(String enzyme1: re.keySet()){
			//for(String enzyme2: re.keySet()){
			//	if(!enzyme1.equals(enzyme2)){
					//System.err.println(enzyme1+" "+enzyme2);
					/*Map<String, Sequence> temp=new TreeMap<String, Sequence>();
					temp.put(enzyme1, re.get(enzyme1));
					temp.put(enzyme2, re.get(enzyme2));*/
					new DNASizesFromRestriction(fastaDir, re, save);
				//}
			//}
		//}
		
		
	}

	
	
	
	
}
