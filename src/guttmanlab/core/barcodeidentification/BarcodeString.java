package guttmanlab.core.barcodeidentification;

import java.util.ArrayList;
import java.util.Collection;

public class BarcodeString {

	Collection<Barcode>[] order;
	
	public BarcodeString(int num){
		order=new ArrayList[num];
		for(int i=0; i<order.length; i++){
			order[i]=new ArrayList<Barcode>();
		}
	}
	
	public void add(String barcodeName, String sequence, int position, int index){
		Barcode b=new Barcode(barcodeName, sequence, index);
		order[position].add(b);	
	}
	
	public String toString(){
		int count=0;
		String rtrn="";
		for(int i=0; i<order.length; i++){
			String name=resolve(order[i], i);
			rtrn+="\t"+i+":"+order[i].size()+":"+name;
			if(i>0){count=Math.max(order[i].size(), count);}
		}
		
		if(count>1){System.err.println(rtrn);}
		
		return rtrn;
	}
	
	
	public boolean isComplete(){
		for(int i=0; i<order.length; i++){
			if(order[i].isEmpty() || order[i].size()==0){return false;}
		}
		return true;
	}
	
	private String resolve(Collection<Barcode> collection, int position) {
		if(collection.isEmpty() || collection.size()==0){return "NF";}
		else if(collection.size()==1){return collection.iterator().next().name;}
		
		else if(position ==0){
			for(Barcode b: collection){
				if(b.position==0){return b.name;}
			}
			
		}
		
		String rtrn="";
		for(Barcode b: collection){rtrn+=":"+b.position;}
		return rtrn;
	}


	private class Barcode{
		String name;
		String sequence;
		int position;
		
		Barcode(String name, String sequence, int position){
			this.name=name;
			this.sequence=sequence;
			this.position=position;
		}
		
	}
	
}


