package guttmanlab.core.rnasprite;

import java.util.Collection;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;

public class RNAInterval extends SingleInterval{

	private Collection<String> introns;
	private Collection<String> exons;
	private Collection<String> other;
	private boolean hasName;
	private String name;
	
	
	
	public RNAInterval(String refName, int start, int end, Strand orientation) {
		super(refName, start, end, orientation);
		this.exons=new TreeSet<String>();
		this.introns=new TreeSet<String>();
		this.other=new TreeSet<String>();
		this.hasName=false;
	}

	
	public RNAInterval(SingleInterval interval) {
		super(interval.getReferenceName(), interval.getReferenceStartPosition(), interval.getReferenceEndPosition(), interval.getOrientation(), interval.getName());
		this.exons=new TreeSet<String>();
		this.introns=new TreeSet<String>();
		this.other=new TreeSet<String>();
		this.hasName=false;
	}

	public Collection<String> getGeneNames(){
		Collection<String> rtrn=new TreeSet<String>();
		rtrn.addAll(exons);
		rtrn.addAll(introns);
		rtrn.addAll(other);
		return rtrn;
	}
	
	public Collection<String> getIntronNames(){return introns;}
	
	public Collection<String> getExonNames(){return exons;}
	
	public Collection<String> getOtherNames(){return other;}
	
	public Collection<String> getExonOnlyNames(){
		Collection<String> rtrn=new TreeSet<String>();
		for(String exon: exons){
			if(!introns.contains(exon)){rtrn.add(exon);}
		}	
		return rtrn;
	}
	
	public Collection<String> getIntronOnlyNames(){
		Collection<String> rtrn=new TreeSet<String>();
		for(String intron: introns){
			if(!exons.contains(intron)){rtrn.add(intron);}
		}	
		return rtrn;
	}
	
	
	public boolean isExon(){
		if(!exons.isEmpty() && introns.isEmpty()){return true;}
		return false;
	}
	
	public boolean hasExon() {
		return !exons.isEmpty();
	}
	
	public boolean isIntron(){
		if(!introns.isEmpty() && exons.isEmpty()){return true;}
		return false;
	}

	public String getType(){
		if(isExon()){return "exon";}
		if(isIntron()){return "intron";}
		return "ambiguous";
		
	}
	
	
	@Override
	public void setName(String name){
		hasName=true;
		this.name=name;
	}
	
	@Override
	public String getName(){
		if(hasName) {return this.name;}
		if(!exons.isEmpty()){return exons.iterator().next();}
		if(!introns.isEmpty()){return introns.iterator().next();}
		if(!other.isEmpty()) {return other.iterator().next();}
		return "";
	}
	

	public void addNameType(String name, String type) {
		if(type.equalsIgnoreCase("exon")){exons.add(name);}
		if(type.equalsIgnoreCase("intron")){introns.add(name);}
		else {other.add(name);}
		
	}


	public Collection<String> getAllRNANames() {
		Collection<String> rtrn=new TreeSet<String>();
		
		rtrn.addAll(getExonNames());
		rtrn.addAll(getIntronNames());
		rtrn.addAll(getOtherNames());
		
		return rtrn;
	}
	
	
}
