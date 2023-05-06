package guttmanlab.core.scSPRITE;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.rnasprite.Cluster;

public class ContactsByBin {

	public ContactsByBin(BarcodingDataStreaming data, SingleInterval region, int binResolution, String save) throws IOException{
		MatrixWithHeaders mwh=getDNADNAContactMatrix(data, binResolution, false, region);
		mwh.write(save);
	}
	
	
	public MatrixWithHeaders getDNADNAContactMatrix(BarcodingDataStreaming data, int binResolution, boolean weight, SingleInterval region){
		List<String> columns=getGenomePositions(data, binResolution, region.getReferenceName());
		List<String> rows=add(region, binResolution);
		
		MatrixWithHeaders counts=new MatrixWithHeaders(rows, columns);
		
		scoreDNADNA(data, counts, binResolution, weight);
		
		return counts;
	}
	
	
	private List<String> getGenomePositions(BarcodingDataStreaming data, int binResolution, String referenceName) {
		TreeSet<SingleInterval> set=new TreeSet<SingleInterval>();
		while(data.hasNext()){	
			Cluster c=data.next();
			Collection<SingleInterval> regions=c.getAllDNAIntervals();
			for(SingleInterval region: regions){
				if(region.getReferenceName().equals(referenceName)){
					SingleInterval binned=region.bin(binResolution);
					set.add(binned);
				}
			}
		}
		
		List<String> rtrn=new ArrayList<String>();
		for(SingleInterval r: set){rtrn.add(r.toUCSC());}
		
		data.close();
		return rtrn;
	}


	private List<String> add(SingleInterval region, int binResolution) {
		List<String> rtrn=new ArrayList<String>();
		
		SingleInterval binned=region.bin(binResolution);
		rtrn.add(binned.toUCSC());
		
		return rtrn;
	}


	private void scoreDNADNA(BarcodingDataStreaming data, MatrixWithHeaders counts, int binResolution, boolean weight) {
		
		while(data.hasNext()){
			Cluster c=data.next();
			Cluster binned=c.bin(binResolution);
			for(SingleInterval region1: binned.getAllDNAIntervals()){
				for(SingleInterval region2: binned.getAllDNAIntervals()){
					String row=region1.toUCSC();
					String column=region2.toUCSC();
					if(counts.containsColumn(column) && counts.containsRow(row)){
						//System.err.println(row+" "+column);
						double count=counts.get(row, column);
						double score=2.0/c.getAllDNAIntervals().size();
						//count++;
						
						if(weight){
							count+=score;
						}
						else{
							count++;
						}
						//System.err.println(count);
						counts.set(row, column, count);
					}
				}
			}
		
		}
		data.close();
		
		
	}
	
	
	public static void main(String[] args) throws IOException{
		if(args.length>3){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			SingleInterval region=new SingleInterval(args[1]);
			int binResolution=new Integer(args[2]);
			String save=args[3];
			new ContactsByBin(data, region, binResolution, save);
		}
		else{System.err.println(usage);}
		
	}
	
	static String usage=" args[0]=cluster file \n args[1]=region \n args[2]=bin resolution \n args[3]=save";
	
}
