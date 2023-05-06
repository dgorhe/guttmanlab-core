package guttmanlab.core.rnasprite;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.barcoding.analysis.Cluster;


public class Kmer implements Comparable<Kmer>{

	Collection<String> regions;
	Collection<SingleInterval> coordinates;
	boolean isInput=false;
	String name;
	
	public Kmer(){
		this.regions=new TreeSet<String>();
		this.coordinates=new TreeSet<SingleInterval>();
		name="";
	}
	
	public Kmer(String str){
		initializeFromString(str);
	}
	
	
	public Kmer(Kmer cluster) {
		this.regions=new TreeSet<String>();
		this.coordinates=new TreeSet<SingleInterval>();
		
		this.regions.addAll(cluster.getRegions());
		this.coordinates.addAll(cluster.getIntervals());
		this.name=cluster.getName();
	}

	public boolean isInput(){return this.isInput;}
	
	private void initializeFromString(String str) {
		this.name="";
		if(str.equalsIgnoreCase("input")){isInput=true;}
		this.regions=new TreeSet<String>();
		String[] tokens=str.split("_");
		for(int i=0; i<tokens.length; i++){
			addRegion(tokens[i]);
		}
	}

	public void addRegion(String region){
		regions.add(region);
	}
	
	public void addInterval(SingleInterval interval){
		coordinates.add(interval);
	}
	
	public Collection<SingleInterval> getIntervals(){
		return coordinates;
	}
	
	@Override
	public boolean equals(Object o){
		Kmer other=(Kmer)o;
		
		
		if(this.getSize()!=other.getSize()){return false;}
		
		for(String read: this.regions){
			if(!other.regions.contains(read)){return false;}
		}
		
		return true;
	}
	
	public int getSize() {
		return regions.size();
	}

	public String toString(){
		String rtrn="";
		int counter=0;
		for(String region: regions){
			rtrn+=region;
			counter++;
			if(counter<regions.size()) {rtrn+="_";}
		}
		
		return rtrn;
	}

	@Override
	public int compareTo(Kmer cluster2) {
		if(this.getSize()!=cluster2.getSize()){return getSize()-cluster2.getSize();}
		
		Iterator<String> iter1=regions.iterator();
		Iterator<String> iter2=cluster2.regions.iterator();
		
		while(iter1.hasNext()){
			String interval1=iter1.next();
			String interval2=iter2.next();
			if(!interval1.equals(interval2)){return interval1.compareTo(interval2);}
		}
		
		if(!this.getName().equals(cluster2.getName())) {return this.getName().compareTo(cluster2.getName());}
		
		return 0;
	}

	/**
	 * Test if kmer2 is a subset of this kmer
	 * @param kmer2
	 * @return
	 */
	public boolean isSubset(Kmer kmer2) {
		//Are all of the regions in kmer2 within this kmer
		for(String region: getRegions()){
			if(!kmer2.hasRegion(region)){return false;}
		}
		return true;
	}

	private boolean hasRegion(String region) {
		return this.regions.contains(region);
	}

	public Collection<String> getRegions() {
		return this.regions;
	}

	public String toFileName() {
		String rtrn="";
		int counter=0;
		for(String region: regions){
			rtrn+=region;
			if(counter<regions.size()) {rtrn+="_";}
			counter++;
		}
		return rtrn;
	}

	public Kmer remove(String region) {
		Kmer rtrn=new Kmer();
		for(String r: this.regions){
			if(!r.equalsIgnoreCase(region)){
				rtrn.addRegion(r);
			}
		}
		return rtrn;
	}

	public void addRegions(Collection<String> names) {
		for(String name: names){
			addRegion(name);
		}
		
	}

	public void addIntervals(Collection<SingleInterval> set) {
		coordinates.addAll(set);
		
	}
	
	
	public boolean removeInterval(SingleInterval r) {
		return coordinates.remove(r);
	}
	
	public void setName(String name) {this.name=name;}

	public String getName() {
		return name;
	}

	public Collection<Kmer> enumerateSubK(int k) {
		Collection<String> regions=this.getRegions();		
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		
		for(String r: regions){
				Kmer cluster=new Kmer();
				cluster.setName(getName());
				cluster.addRegion(r);
				rtrn.add(cluster);
			}
				
			for(int i=2; i<=k; i++){
				rtrn=add(regions, rtrn, i);
			}
				
			return rtrn;
		}
	
	private Collection<Kmer> add(Collection<String> regions, Collection<Kmer> clusters, int k) {
		//Iterate through each cluster and add each interval
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		
		for(Kmer cluster: clusters){
			for(String interval: regions){
				Kmer newCluster=new Kmer(cluster);
				newCluster.addRegion(interval);
				if(newCluster.getSize()==k){
					rtrn.add(newCluster);
				}
			}
		}
		
		return rtrn;
	}

	public List<String> getRegionList() {
		List<String> rtrn=new ArrayList<String>();
		
		rtrn.addAll(getRegions());
		
		return rtrn;
	}

	public boolean containsRegion(String col) {
		return regions.contains(col);
	}

	public Collection<Kmer> subtractOne() {
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		
		for(String region: this.getRegions()) {
			Kmer k=this.remove(region);
			rtrn.add(k);
		}
		
		return rtrn;
	}

	public String toString(String sep) {
		String rtrn="";
		int counter=0;
		for(String region: regions){
			rtrn+=region;
			counter++;
			if(counter<regions.size()) {rtrn+=sep;}
		}
		
		return rtrn;
	}

	public void addRegions(String[] regions) {
		for(String r: regions) {addRegion(r);}
		
	}

	public int distance(Kmer k2) {
		if(this.getSize()!=k2.getSize()) {return -1;}
		
		Iterator<String> iter1=this.getRegions().iterator();
		Iterator<String> iter2=k2.getRegions().iterator();
		
		int d=0;
		
		while(iter1.hasNext()) {
			String s1=iter1.next();
			String s2=iter2.next();
			if(!s1.equals(s2)) {d++;}
		}
		
		return d;
	}

	public Collection<Kmer> enumerateSubK() {
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		rtrn.add(this);
		
		for(int i=1; i<this.getSize(); i++) {
			rtrn.addAll(enumerateSubK(i));
		}
		return rtrn;
	}

	public int getCoordinateSize() {
		return this.coordinates.size();
	}
	
	
}
