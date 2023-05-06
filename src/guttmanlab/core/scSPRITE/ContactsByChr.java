package guttmanlab.core.scSPRITE;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.rnasprite.Cluster;

public class ContactsByChr {
	
	int maxClusterSize=10000;

	public ContactsByChr(BarcodingDataStreaming data, String chr, int binResolution, String save) throws IOException, InterruptedException{
		MatrixWithHeaders mwh=getDNADNAContactMatrix(data, binResolution, false, chr);
		mwh.write(save);
		/*writeAndNorm(mwh, save);
		
		FileWriter writer=new FileWriter(save+".binned.cluster");
		while(data.hasNext()){
			Cluster c=data.next();
			Cluster binned=c.bin(binResolution);
			Cluster temp=new Cluster(binned.getBarcode());
			for(SingleInterval region: binned.getAllDNAIntervals()){
				if(region.getReferenceName().equals(chr)){temp.addDNARead(region);}
			}
			if(temp.getClusterSize()>0){
				writer.write(temp.toString()+"\n");
			}
			
		}
		writer.close();*/
	}
	
	
	private File writeAndNorm(MatrixWithHeaders observed, String input) throws IOException, InterruptedException {
		observed.write(input);
		
		String cmd="python /groups/guttman/SPRITE2/ice_matrix.py "+input;
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		
		File rtrn=new File(input+".iced");
		return rtrn;
	}
	
	
	public MatrixWithHeaders getDNADNAContactMatrix(BarcodingDataStreaming data, int binResolution, boolean weight, String chr){
		List<String> columns=getGenomePositions(data, binResolution, chr);
		List<String> rows=columns;
		
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
				if(region.getReferenceName().equals(referenceName) || referenceName.equalsIgnoreCase("all")){
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
			if(c.getClusterSize()<this.maxClusterSize){
				Cluster binned=c.bin(binResolution);
				for(SingleInterval region1: binned.getAllDNAIntervals()){
					for(SingleInterval region2: binned.getAllDNAIntervals()){
						String row=region1.toUCSC();
						String column=region2.toUCSC();
						if(counts.containsColumn(column) && counts.containsRow(row)){
							//System.err.println(row+" "+column);
							double count=counts.get(row, column);
							double score=1.0/c.getAllDNAIntervals().size();
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
		}
		data.close();
		
		
	}
	
	
	public static void main(String[] args) throws IOException, InterruptedException{
		if(args.length>3){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			String chr=args[1];
			int binResolution=new Integer(args[2]);
			String save=args[3];
			new ContactsByChr(data, chr, binResolution, save);
		}
		else{System.err.println(usage);}
		
	}
	
	static String usage=" args[0]=cluster file \n args[1]=chr \n args[2]=bin resolution \n args[3]=save";
	
}
