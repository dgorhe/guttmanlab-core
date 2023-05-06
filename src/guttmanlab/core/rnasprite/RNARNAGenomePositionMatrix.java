package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;
import guttmanlab.core.simulation.CoordinateSpace;

public class RNARNAGenomePositionMatrix {
	
	int minCount=1000;
	int binResolution=1000000;

	public RNARNAGenomePositionMatrix(BarcodingDataStreaming data, String save, CoordinateSpace space) throws IOException, InterruptedException{
		MatrixWithHeaders mwh=data.getRNARNAContactMatrix(binResolution, false, space);
		
		mwh=filter(mwh, minCount);
		
		//Get gene intron matrix
		/*Map<String, SingleInterval> names=getGeneIntronNames(data);
		
		//Initialize matrix
		MatrixWithHeaders mwh=initialize(names);
		
		
		//populate matrix
		scoreMatrix(data, mwh);*/
		
		
		//write, ICE normalize, output annotation order
		normalize(mwh, save);
		
		
		//write matrix + annotation order
		
	}
	
	


	private MatrixWithHeaders filter(MatrixWithHeaders mwh, int minCount2) {
		List<String> list=new ArrayList<String>();
		for(String row: mwh.getRowNames()){
			double[] vals=mwh.getRow(row);
			double sum=Statistics.sum(vals);
			if(sum>minCount2){list.add(row);}
		}
		mwh=mwh.submatrixByRowNames(list);
		mwh=mwh.submatrixByColumnNames(list);
		return mwh;
	}




	private void normalize(MatrixWithHeaders mwh, String save) throws IOException, InterruptedException {
		mwh.write(save);
		
		Process p=Runtime.getRuntime().exec("python /groups/guttman/SPRITE2/ice_matrix.py "+save);
		p.waitFor();
	}
	
	private void normalize(MatrixWithHeaders mwh, Map<String, SingleInterval> names, String save) throws IOException, InterruptedException {
		mwh.write(save);
		
		Process p=Runtime.getRuntime().exec("python /groups/guttman/SPRITE2/ice_matrix.py "+save);
		p.waitFor();
		
		//write annotations
		FileWriter writer=new FileWriter(save+".annotation");
		
		writer.write("ID\tchr\tstart\tposition\n");
		for(String name: names.keySet()){
			SingleInterval interval=names.get(name);
			writer.write(name+"\t"+ interval.getReferenceName()+"\t"+interval.getReferenceStartPosition()+"\t"+interval.toUCSC()+"\n");
		}
		writer.close();
	}


	private void scoreMatrix(BarcodingDataStreaming data, MatrixWithHeaders mwh) {
		
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			Collection<String> intronList=new TreeSet<String>();
			
			for(RNAInterval r: c.getAllRNARegions()){
				if(r.isIntron() && mwh.hasColumn(r.getName())){
					intronList.add(r.getName());
				}
			}
			
			for(String name1: intronList){
				for(String name2: intronList){
					mwh.incrementCount(name1, name2);
				}
			}
			counter++;
			if(counter %1000000 ==0){System.err.println(counter);}
		}
		data.close();
	}


	


	private MatrixWithHeaders initialize(Map<String, SingleInterval> names) {
		List<String> rows=new ArrayList<String>();
		rows.addAll(names.keySet());
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, rows);
		return rtrn;
	}


	private Map<String, SingleInterval> getGeneIntronNames(BarcodingDataStreaming data){
		Map<String, SingleInterval> names=new TreeMap<String, SingleInterval>();
		Map<String, Integer> counts=new TreeMap<String, Integer>();
		
		
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			for(RNAInterval r: c.getAllRNARegions()){
				if(r.isIntron()){
					String name=r.getName();
					names.put(name, r.getSingleInterval());
					
					int count=0;
					if(counts.containsKey(name)){count=counts.get(name);}
					count++;
					counts.put(name, count);
				}
			}
			counter++;
			if(counter %1000000 ==0){System.err.println(counter);}
		}
		data.close();
		
		
		Map<String, SingleInterval> rtrn=new TreeMap<String, SingleInterval>();
		
		for(String name: names.keySet()){
			SingleInterval region=names.get(name);
			int score=counts.get(name);
			if(score>this.minCount){rtrn.put(name, region);}
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException, InterruptedException{
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String save=args[1];
		CoordinateSpace space=CoordinateSpace.HG38;
		new RNARNAGenomePositionMatrix(data, save, space);
		
	}
	
}
