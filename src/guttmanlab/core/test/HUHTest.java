package guttmanlab.core.test;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;
import guttmanlab.core.simulation.CoordinateSpace;

public class HUHTest {

	private static Map<Integer, Map<SingleInterval, Integer>> parse(String string, int binResolution) throws IOException {
		Collection<String> lines=BEDFileIO.loadLines(string);
		
		
		Map<Integer, Map<SingleInterval, Integer>> mapByHUHSize=new TreeMap<Integer, Map<SingleInterval, Integer>>();
		
		
		Map<Integer, Integer> huhCounts=new TreeMap<Integer, Integer>();
		
		for(String line: lines){
			//Collection<String> huhTags=parseHUH(line);
			int huhSize=getHUHSize(line);
			Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
			if(mapByHUHSize.containsKey(huhSize)){rtrn=mapByHUHSize.get(huhSize);}
			
			
			Collection<SingleInterval> regions=parse(line);
			//System.out.println(huhTags.size()+"\t"+regions.size());
			int counter=0;
			if(huhCounts.containsKey(huhSize)){counter=huhCounts.get(huhSize);}
			counter++;
			huhCounts.put(huhSize, counter);
			for(SingleInterval region: regions){
					SingleInterval binned=region.bin(binResolution);
					int count=0;
					if(rtrn.containsKey(binned)){
						count=rtrn.get(binned);
					}
					count++;
					rtrn.put(binned, count);
				}
				
				mapByHUHSize.put(huhSize, rtrn);
				
			//}
		}
		
		
		
		return mapByHUHSize;
	}
	

	private static List<Integer> get(String file, String chr) throws IOException {
		Collection<String> lines=BEDFileIO.loadLines(file);
		List<Integer> list=new ArrayList<Integer>();
		
		
		for(String line: lines){
			//Collection<String> huhTags=parseHUH(line);
			int huhSize=getHUHSize(line);
			Collection<SingleInterval> regions=parse(line);
			if(contains(chr, regions)){list.add(huhSize);}
			
		}
		
		
		
		return list;
	}
	
	
	private static boolean contains(String chr, Collection<SingleInterval> regions) {
		boolean containsChr=false;
		boolean oneChr=true;
		Collection<String> chrs=new TreeSet<String>();
		for(SingleInterval region: regions){
			chrs.add(region.getReferenceName());
			if(region.getReferenceName().equalsIgnoreCase(chr)){containsChr=true;}
		}
		oneChr=(chrs.size()==1);
		return containsChr&&oneChr;
		//return false;
	}


	private static MatrixWithHeaders makePairwiseContacts(String string, int binResolution, String chr) throws IOException {
		Collection<String> lines=BEDFileIO.loadLines(string);
		
		List<String> list=CoordinateSpace.HG19.getBins(binResolution, chr);
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(list, list);
		
		for(String line: lines){
			int huhSize=getHUHSize(line);
			
			Collection<SingleInterval> regions=parse(line, chr);
			for(SingleInterval region1: regions){
				
				SingleInterval binned1=region1.bin(binResolution);
				for(SingleInterval region2: regions){
					
					if(!region1.equals(region2)){
						SingleInterval binned2=region2.bin(binResolution);
						//System.err.println(binned1.toUCSC());
						//System.err.println(binned2.toUCSC());
						double score=mwh.get(binned1.toUCSC(), binned2.toUCSC());
						score+=1;
						mwh.set(binned1.toUCSC(), binned2.toUCSC(), score);
					}
				}
			}
		}
		
		
		
		
		return mwh;
	}
	
	
	
	
	
	private static MatrixWithHeaders makePairwiseContacts(String string, int binResolution, SingleInterval includeRegion) throws IOException {
		Collection<String> lines=BEDFileIO.loadLines(string);
		
		List<String> list=CoordinateSpace.HG19.getBins(binResolution, includeRegion);
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(list, list);
		
		for(String line: lines){
			int huhSize=getHUHSize(line);
			
			Collection<SingleInterval> regions=parse(line, includeRegion);
			for(SingleInterval region1: regions){
				
				SingleInterval binned1=region1.bin(binResolution);
				for(SingleInterval region2: regions){
					
					if(!region1.equals(region2)){
						SingleInterval binned2=region2.bin(binResolution);
						//System.err.println(binned1.toUCSC());
						//System.err.println(binned2.toUCSC());
						if(mwh.containsColumn(binned2.toUCSC()) && mwh.containsRow(binned1.toUCSC())){
							double score=mwh.get(binned1.toUCSC(), binned2.toUCSC());
							score+=1;
							mwh.set(binned1.toUCSC(), binned2.toUCSC(), score);
						}
						else{System.err.println(binned1.toUCSC()+ " "+binned2.toUCSC());}
					}
				}
			}
		}
		
		
		
		
		return mwh;
	}
	
	
	private static int getHUHSize(String line) {
		Collection<String> rtrn=new TreeSet<String>();
		String[] tokens=line.split("\t");
			
		for(int i=1; i<tokens.length; i++){
			if(tokens[i].startsWith("HUHNOSPLINT")){rtrn.add(tokens[i]);}
		}
		
		if(rtrn.isEmpty()){return 0;}
		return rtrn.size();
	}


	private static Collection<SingleInterval> parse(String line) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		//NYBot10_Stg.Odd2Bo51.Even2Bo19.Odd2Bo29.Even2Bo90.Odd2Bo17.Even2Bo9.Odd2Bo86	HUHNOSPLINT	HUHNOSPLINT	HUHNOSPLINT	HUHNOSPLINT	HUHNOSPLINT	HUHNOSPLINT	HEK-EZH2-16-chr1:90057393	HEK-EZH2-16-chr6:78716936
		
		String[] tokens=line.split("\t");
		boolean huh=false;
		for(int i=1; i<tokens.length; i++){
			if(tokens[i].startsWith("HUHNOSPLINT")){huh=true;}
			else{
				String[] region=tokens[i].split("-");
				String pos=region[region.length-1];
				int start=new Integer(pos.split(":")[1]);
				int end=start+1;
				SingleInterval r=new SingleInterval(pos.split(":")[0], start, end);
				rtrn.add(r);
			}
		}
		//if(huh){return rtrn;}
		return rtrn;
	}
	
	
	private static Collection<SingleInterval> parse(String line, String chr) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		//NYBot10_Stg.Odd2Bo51.Even2Bo19.Odd2Bo29.Even2Bo90.Odd2Bo17.Even2Bo9.Odd2Bo86	HUHNOSPLINT	HUHNOSPLINT	HUHNOSPLINT	HUHNOSPLINT	HUHNOSPLINT	HUHNOSPLINT	HEK-EZH2-16-chr1:90057393	HEK-EZH2-16-chr6:78716936
		
		String[] tokens=line.split("\t");
		boolean huh=false;
		for(int i=1; i<tokens.length; i++){
			if(tokens[i].startsWith("HUHNOSPLINT")){huh=true;}
			else{
				String[] region=tokens[i].split("-");
				String pos=region[region.length-1];
				int start=new Integer(pos.split(":")[1]);
				int end=start+1;
				SingleInterval r=new SingleInterval(pos.split(":")[0], start, end);
				if(r.getReferenceName().equalsIgnoreCase(chr)){rtrn.add(r);}
			}
		}
		return rtrn;
	}
	
	
	private static Collection<SingleInterval> parse(String line, SingleInterval includeRegion) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		//NYBot10_Stg.Odd2Bo51.Even2Bo19.Odd2Bo29.Even2Bo90.Odd2Bo17.Even2Bo9.Odd2Bo86	HUHNOSPLINT	HUHNOSPLINT	HUHNOSPLINT	HUHNOSPLINT	HUHNOSPLINT	HUHNOSPLINT	HEK-EZH2-16-chr1:90057393	HEK-EZH2-16-chr6:78716936
		
		String[] tokens=line.split("\t");
		boolean huh=false;
		for(int i=1; i<tokens.length; i++){
			if(tokens[i].startsWith("HUHNOSPLINT")){huh=true;}
			else{
				String[] region=tokens[i].split("-");
				String pos=region[region.length-1];
				int start=new Integer(pos.split(":")[1]);
				int end=start+1;
				SingleInterval r=new SingleInterval(pos.split(":")[0], start, end);
				if(r.overlaps(includeRegion)){rtrn.add(r);}
			}
		}
		return rtrn;
	}

	
	private static Collection<String> parseHUH(String line) {
		Collection<String> rtrn=new TreeSet<String>();
		//NYBot10_Stg.Odd2Bo51.Even2Bo19.Odd2Bo29.Even2Bo90.Odd2Bo17.Even2Bo9.Odd2Bo86	HUHNOSPLINT	HUHNOSPLINT	HUHNOSPLINT	HUHNOSPLINT	HUHNOSPLINT	HUHNOSPLINT	HEK-EZH2-16-chr1:90057393	HEK-EZH2-16-chr6:78716936
		
		String[] tokens=line.split("\t");
		
		for(int i=1; i<tokens.length; i++){
			if(tokens[i].startsWith("HUHNOSPLINT")){rtrn.add(tokens[i]);}
			
		}
		//if(huh){return rtrn;}
		return rtrn;
	}
	
	
	

	private static void write(Map<SingleInterval, Integer> data, String string) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		for(SingleInterval region: data.keySet()){
			writer.write(region.toBedgraph(data.get(region))+"\n");
		}
		
		writer.close();
	}
	
	private static Collection<SingleInterval> getNeighbors(SingleInterval region, int binResolution, int n) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		//upstream
		for(int i=0; i<n; i++){
			int start=region.getReferenceStartPosition()-((i+1)*binResolution);
			int end=start+binResolution;
			SingleInterval newRegion=new SingleInterval(region.getReferenceName(), start, end);
			rtrn.add(newRegion);
		}
		
		
		//downstream
		for(int i=0; i<n; i++){
			int start=region.getReferenceStartPosition()+((i+1)*binResolution);
			int end=start+binResolution;
			SingleInterval newRegion=new SingleInterval(region.getReferenceName(), start, end);
			rtrn.add(newRegion);
		}
		
		return rtrn;
	}
	
	private static double average(Collection<SingleInterval> regions, Map<SingleInterval, Integer> map) {
		double sum=0;
		double count=0;
		
		for(SingleInterval region: regions){
			if(map.containsKey(region)){
				sum+=map.get(region);
				count++;
			}
		}
		
		return sum/count;
	}
	
	
	private static void writeSmoothed(String save, Map<SingleInterval, Integer> map, int binResolution, int n) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: map.keySet()){
			Collection<SingleInterval> regions=getNeighbors(region, binResolution, n);
			double val=average(regions, map);
			writer.write(region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+region.getReferenceEndPosition()+"\t"+val+"\n");
		}
		writer.close();
	}
	
	
	private static void writeAll(Map<Integer, Map<SingleInterval, Integer>> mapByHUHSize, String string) throws IOException {
		Collection<SingleInterval> regions=new TreeSet<SingleInterval>();
		for(Integer huhSize: mapByHUHSize.keySet()){
			String save=string+".huh"+huhSize+".bedgraph";
			regions.addAll(mapByHUHSize.get(huhSize).keySet());
			//write(mapByHUHSize.get(huhSize), save);
		}
		
		
		FileWriter writer=new FileWriter(string+".0_10.bedgraph");
		
		for(SingleInterval region: regions){
			int sum=sum(mapByHUHSize, region, 0, 10);
			writer.write(region.toBedgraph(sum)+"\n");
		}
		
		writer.close();
		
		
		 writer=new FileWriter(string+".11_100.bedgraph");
		
		for(SingleInterval region: regions){
			int sum=sum(mapByHUHSize, region, 11, 100);
			writer.write(region.toBedgraph(sum)+"\n");
		}
		
		writer.close();
		
		 writer=new FileWriter(string+".101_1000.bedgraph");
			
			for(SingleInterval region: regions){
				int sum=sum(mapByHUHSize, region, 101, 1000);
				writer.write(region.toBedgraph(sum)+"\n");
			}
			
			writer.close();
			
			
			writer=new FileWriter(string+".10_10K.bedgraph");
			
			for(SingleInterval region: regions){
				int sum=sum(mapByHUHSize, region, 10, 10000);
				writer.write(region.toBedgraph(sum)+"\n");
			}
			
			writer.close();
			
			
			writer=new FileWriter(string+".100_10K.bedgraph");
			
			for(SingleInterval region: regions){
				int sum=sum(mapByHUHSize, region, 100, 10000);
				writer.write(region.toBedgraph(sum)+"\n");
			}
			
			writer.close();
			
			
			writer=new FileWriter(string+".10_1K.bedgraph");
			
			for(SingleInterval region: regions){
				int sum=sum(mapByHUHSize, region, 10, 1000);
				writer.write(region.toBedgraph(sum)+"\n");
			}
			
			writer.close();
		
		writer=new FileWriter(string+".sum.bedgraph");
		
		for(SingleInterval region: regions){
			int sum=sum(mapByHUHSize, region, 0, 10000);
			writer.write(region.toBedgraph(sum)+"\n");
		}
		
		writer.close();
	}
	
	private static int sum(Map<Integer, Map<SingleInterval, Integer>> mapByHUHSize, SingleInterval region, int minSize, int maxSize) {
		int sum=0;
		for(Integer key: mapByHUHSize.keySet()){
			if(key>= minSize && key<=maxSize){
				if(mapByHUHSize.get(key).containsKey(region)){sum+=mapByHUHSize.get(key).get(region);}
			}
		}
		return sum;
	}
	
	private static MatrixWithHeaders makePairwiseContacts(Map<Integer, Map<SingleInterval, Integer>> mapByHUHSize) {
		// TODO Auto-generated method stub
		return null;
	}

	
	private static MatrixWithHeaders remove(MatrixWithHeaders mwh, SingleInterval exclude) {
		List<String> rowsToUse=new ArrayList<String>();
		for(String region: mwh.getRowNames()){
			SingleInterval temp=new SingleInterval(region);
			if(!temp.overlaps(exclude)){rowsToUse.add(region);}
		}
		return mwh.submatrixByRowNames(rowsToUse).submatrixByColumnNames(rowsToUse);
		
	}
	
	private static void writeList(List<Integer> list, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Integer val: list){
			writer.write(val+"\n");
		}
		
		writer.close();
	}

	

	public static void main(String[] args) throws IOException{
		if(args.length>2){
			
			//int huhSize=new Integer(args[3]);	
			//Map<Integer, Map<SingleInterval, Integer>> mapByHUHSize=parse(args[0], new Integer(args[2]));
			//writeAll(mapByHUHSize, args[1]);
		
			
			
			//SingleInterval includeRegion=new SingleInterval(args[3]);
			
			//MatrixWithHeaders mwh= makePairwiseContacts(args[0], new Integer(args[2]), includeRegion);
			
			
			for(String chr:CoordinateSpace.HG19.getRefSizes().keySet()){
				List<Integer> list=get(args[0], chr);
				if(list.size()>20){
					double val=Statistics.quantileInt(list, 0.99);
					double max=Statistics.max(list);
					double average=Statistics.mean(list);
					System.err.println(chr+" "+list.size()+" "+val+" "+max+" "+average);
				}
			}
			//writeList(list, args[1]+"."+args[3]+".sizes");
			
		//	mwh=remove(mwh, exclude);
			
			
			//mwh.write(args[1]+".pairwise.matrix");
		}
		else{System.err.println(usage);}
	}

	
	static String usage=" args[0]=clusters \n args[1]=save \n args[2]=bin resolution \n args[3]=chr";


	

	
}
