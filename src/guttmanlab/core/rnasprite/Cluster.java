package guttmanlab.core.rnasprite;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.Pair;

public class Cluster {

	private Collection<RNAInterval> rnaRegions;
	private Collection<SingleInterval> dnaRegions;
	private Collection<String> rnaNames;
	private Collection<String> proteins;
	private List<String> rnaList;
	String barcode;
	String line;
	
	public Cluster(String barcode){
		this.barcode=barcode;
		this.rnaRegions=new TreeSet<RNAInterval>();
		this.dnaRegions=new TreeSet<SingleInterval>();
		this.rnaNames=new TreeSet<String>();
		this.proteins=new ArrayList<String>();
		this.rnaList=new ArrayList<String>();
	}

	

	public void addDNARead(SingleInterval interval) {
		this.dnaRegions.add(interval);
	}
	
	public void addRNARead(RNAInterval interval) {
		this.rnaRegions.add(interval);
		this.rnaNames.add(interval.getName());
	}

	public Collection<RNAInterval> getAllRNARegions() {
		return rnaRegions;
	}
	
	@Override
	public boolean equals(Object o){
		Cluster other=(Cluster)o;
		
		//if(!this.barcode.equalsIgnoreCase(other.getBarcode())){return false;} //TODO remove
		
		if(this.rnaRegions.size()!=other.rnaRegions.size()){return false;}
		if(this.dnaRegions.size()!=other.dnaRegions.size()){return false;}
		
		for(SingleInterval read: this.dnaRegions){
			if(!other.dnaRegions.contains(read)){return false;}
		}
		for(SingleInterval read: this.rnaRegions){
			if(!other.rnaRegions.contains(read)){return false;}
		}
		
		return true;
	}

	public Collection<SingleInterval> getAllDNAIntervals() {
		return dnaRegions;
	}

	public int size() {
		return this.dnaRegions.size()+this.rnaRegions.size()+this.proteins.size();
	}

	public int getClusterSize() {
		return size();
	}
	
	public int getRNAClusterSize() {
		return this.rnaRegions.size();
	}

	public String getBarcode() {
		return this.barcode;
	}
	
	//TODO This needs to be fixed
	public Cluster bin(int resolution) {
		Cluster rtrn=new Cluster(getBarcode());
		for(SingleInterval interval: getAllDNAIntervals()){
			int startIndex=interval.getReferenceStartPosition()/resolution;
			int newStart=startIndex*resolution;
			int newEnd=newStart+Math.max(interval.getLength(), resolution);
			SingleInterval newInterval=new SingleInterval(interval.getReferenceName(), newStart, newEnd);
			rtrn.addDNARead(newInterval);	
		}
		
		for(RNAInterval interval: getAllRNARegions()){
			int startIndex=interval.getReferenceStartPosition()/resolution;
			int newStart=startIndex*resolution;
			int newEnd=newStart+Math.max(interval.getLength(), resolution);
			SingleInterval newInterval=new SingleInterval(interval.getReferenceName(), newStart, newEnd);
			RNAInterval bin=new RNAInterval(newInterval);
			rtrn.addRNARead(bin);	
		}
		
		return rtrn;
	}

	/*public void addRNAName(String rnaName) {
		this.rnaNames.add(rnaName);
	}*/
	
	public Collection<String> getRNANames(){return this.rnaNames;}
	
	public Collection<String> getProteins(){return this.proteins;}
	
	public String toString(){
		return toString(this.barcode);
	}
	
	public String toString(String newBarcode){
		String rtrn=newBarcode;
		
		for(SingleInterval dna: this.dnaRegions){
			rtrn+="\tDPM["+dna.getOrientation()+"]_"+dna.toUCSC();
		}
		
		for(RNAInterval rna: this.getAllRNARegions()) {
			rtrn+="\tRPM["+rna.getName()+"]_"+rna.toUCSC();
		}
		
		for(String protein: this.proteins) {
			rtrn+="\tBPM[]_"+protein;
		}
		
		return rtrn;
		//return this.line;
	}
	
	
	public String toStringByChr(String chr){
		String rtrn=this.barcode;
		
		for(SingleInterval dna: this.dnaRegions){
			if(dna.getReferenceName().equals(chr)) {
				rtrn+="\tDPM["+dna.getOrientation()+"]_"+dna.toUCSC();
			}
		}
		
		for(RNAInterval rna: this.getAllRNARegions()) {
			rtrn+="\tRPM["+rna.getName()+"]_"+rna.toUCSC();
		}
		
		for(String protein: this.getProteinSet()) {
			rtrn+="\tBPM[]_"+protein;
		}
		
		return rtrn;
		//return this.line;
	}
	
	public String toCompressedString(){
		String rtrn=this.barcode;
		
		for(SingleInterval dna: this.dnaRegions){
			rtrn+="\tDPM["+dna.getOrientation()+"]_"+dna.toUCSC();
		}
		
		for(RNAInterval rna: this.getAllRNARegions()) {
			rtrn+="\tRPM["+rna.getName()+"]_"+rna.toUCSC();
		}
		
		for(String protein: this.getProteinSet()) {
			rtrn+="\tBPM[]_"+protein;
		}
		
		return rtrn;
		//return this.line;
	}
	
	
	public String toRNAString(){
		String rtrn=this.barcode;
		
		for(RNAInterval rna: this.getAllRNARegions()) {
			rtrn+="\tRPM("+rna.getName()+")_"+rna.toUCSC();
		}
		
		
		return rtrn;
		//return this.line;
	}
	
	public String toUniqueRNAString() {
		String rtrn=this.barcode;
		
		for(String rna: this.getRNANames()) {
			rtrn+="\tRPM("+rna+")_"+rna+":0-1";
		}
		
		
		return rtrn;
		//return this.line;
	}
	
	public String toDNAString(){
		String rtrn=this.barcode;
		
		for(SingleInterval dna: this.dnaRegions){
			rtrn+="\tDNA_"+dna.toUCSC();
		}
		
		return rtrn;
		//return this.line;
	}
	
	public void setLine(String line){this.line=line;}

	public boolean containsOverlappingDNA(SingleInterval region) {
		for(SingleInterval dna: this.dnaRegions){
			if(region.overlaps(dna)){return true;}
		}
		return false;
	}
	
	public boolean containsChromosome(String chr) {
		for(SingleInterval dna: this.dnaRegions){
			if(dna.getReferenceName().equals(chr)){return true;}
		}
		return false;
	}

	/*public Cluster permuteRNA(ArrayList<SingleInterval> allRNAs) {
		int clusterSize=this.getRNANames().size();
		Cluster c=new Cluster("Random");
		
		for(int i=0; i<clusterSize; i++){
			int index=new Double(Math.random()*allRNAs.size()).intValue();
			c.addRNARead(allRNAs.get(index));
		}
		
		return c;
		
	}*/

	public boolean containsRNA(String gene) {
		return rnaNames.contains(gene);
	}
	
	
	public boolean containsRNA(SingleInterval geneRegion) {
		for(RNAInterval region: rnaRegions) {
			if(region.overlaps(geneRegion)) {return true;}
		}
		
		return false;
	}
	
	public boolean containsRNA(String gene, String type) {
		for(RNAInterval rna: getAllRNARegions()) {
			//System.err.println(rna.getName() +" "+rna.getType());
			if(rna.getName().equals(gene) && rna.getType().equals(type)) {return true;}
		}
		
		return false;
	}
	
	public boolean containsRNA(Kmer kmer, boolean requireAll) {
		if(requireAll){
			return rnaNames.containsAll(kmer.getRegions());
		}
		
		for(String r: kmer.getRegions()){
			if(rnaNames.contains(r)){return true;}
		}
		
		return false;
	}
	
	public boolean containsRNA(Collection<String> genes) {
		
		for(String r: genes){
			if(rnaNames.contains(r)){return true;}
		}
		
		return false;
	}
	
	
	public boolean containsRNA(Kmer kmer, boolean requireAll, Map<String, Kmer> collapseSets) {
		if(requireAll){
			//boolean hasAll=true;
			for(String region: kmer.getRegions()){
				Kmer kmer1=get(collapseSets,region);
				boolean has1=false;
				for(String r: kmer1.getRegions()){
					if(rnaNames.contains(r)){has1=true;}
				}
				if(!has1){return false;}
			}
			return true;
			
		}
		
		else{
			Collection<String> regions=new TreeSet<String>();
			for(String region: kmer.getRegions()){
				regions.addAll(get(collapseSets,region).getRegions());
			}
			
			for(String r: regions){
				if(rnaNames.contains(r)){return true;}
			}
		}
		
		return false;
	}
	
	private Kmer get(Map<String, Kmer> collapseSets, String region) {
		if(collapseSets.containsKey(region)){return collapseSets.get(region);}
		
		Kmer rtrn=new Kmer();
		rtrn.addRegion(region);
		return rtrn;
	}

	public boolean containsRNA(Kmer kmer, Map<String, Kmer> collapseSets, int n) {
		//boolean hasAll=true;
		int counter=0;
		for(String region: kmer.getRegions()){
			Kmer kmer1=collapseSets.get(region);
			boolean has1=false;
			for(String r: kmer1.getRegions()){
				if(rnaNames.contains(r)){has1=true;}
			}
			if(has1){counter++;}
		}
		
		return counter>=n;
	}

	public void addRNAReads(Collection<RNAInterval> rnaRegions2) {
		for(RNAInterval rna: rnaRegions2){
			addRNARead(rna);
		}
		
	}

	/*public Cluster renameRNA(Map<String, String> alias) {
		Cluster rtrn=new Cluster(this.barcode);
		
		rtrn.addDNAReads(this.getAllDNAIntervals());
		
		for(SingleInterval rna: this.getAllRNARegions()){
			SingleInterval newRegion=rna;
			if(alias.containsKey(rna.getName())){
				String newName=alias.get(rna.getName());
				newRegion=new SingleInterval(rna.getReferenceName(), rna.getReferenceStartPosition(), rna.getReferenceEndPosition());
				newRegion.setName(newName);
			}
			rtrn.addRNARead(newRegion);
		}
		
		String newLine=getLine(rtrn);
		rtrn.setLine(newLine);
		return rtrn;
	}*/

	private String getLine(Cluster rtrn) {
		String string=rtrn.barcode;
		
		for(SingleInterval dna: rtrn.dnaRegions){
			string+="\tDPM("+dna.getOrientation().toString()+")_"+dna.toUCSC();
		}
		
		for(SingleInterval rna: rtrn.rnaRegions){
			string+="\tRPM("+rna.getName()+")_"+rna.toUCSC();
		}
		
		return string;
	}

	public void addDNAReads(Collection<SingleInterval> allDNAIntervals) {
		for(SingleInterval region: allDNAIntervals){
			addDNARead(region);
		}
	}

	//TODO New parsing
	public static Cluster parseCluster(String line){
		String[] tokens=line.split("\t");
		if(tokens.length==1) {tokens=line.split(" ");}
		
		String barcode=tokens[0];
		Cluster cluster=new Cluster(barcode);
		cluster.setLine(line);
		
		for(int i=1; i<tokens.length; i++){
			parse(tokens[i], cluster);
		}
		return cluster;
	}
	
	
		//TODO New parsing
		public static Cluster parseDNACluster(String line){
			String[] tokens=line.split("\t");
			String barcode=tokens[0];
			Cluster cluster=new Cluster(barcode);
			cluster.setLine(line);
			
			for(int i=1; i<tokens.length; i++){
				String chr=tokens[i].split(":")[0];
				int start=new Integer(tokens[i].split(":")[1]);
				int end= start+1;
				SingleInterval region=new SingleInterval(chr, start, end);
				cluster.addDNARead(region);
			}
			return cluster;
		}
	
	private static void parse(String string, Cluster cluster) {
		boolean isDNA=isDNA(string);
		SingleInterval interval=getCoordinates(string); //get strand and name
		Strand orientation=getOrientation(string);
		interval.setOrientation(orientation);
		
		//System.err.println(string+" "+interval.toUCSC());
		
		if(isDNA){cluster.addDNARead(interval);}
		else{
			RNAInterval rna=new RNAInterval(interval);
			getTypes(string, rna);
			cluster.addRNARead(rna);
			//System.err.println(rna.getName());
		}
	}

	private static void getTypes(String string, RNAInterval rna) {
		String delim1="\\(";
		String delim2="\\)";
		if(string.contains("]")) {
			delim1="\\[";
			delim2="\\]";
		}
		
		String internal=string.split(delim1)[1].split(delim2)[0];
		
		String[] tokens=internal.split(";");
		
		
		
			for(int i=0; i<tokens.length; i++){
				//System.err.println(tokens[i]);
				if(tokens[i].contains(".")){
					String[] split=tokens[i].split("\\.");
					String name=split[0];
					String type=split[1];
					//System.err.println(tokens[i]+ " "+name+" "+type+" "+rna.getReferenceName());
					rna.addNameType(name, type);
				}
				else if(tokens[i].equals("+") || tokens[i].equals("-")) {}
				else {
					String name=tokens[i];
					String type="exon";
					rna.addNameType(name, type);
				}
				
			}
		
		
	}

	private static Strand getOrientation(String token) {
		String delim1="\\(";
		String delim2="\\)";
		if(token.contains("]")) {
			delim1="\\[";
			delim2="\\]";
		}
		
		String internal=token.split(delim1)[1].split(delim2)[0];
		String[] tokens=internal.split(";");
		for(int i=0; i<tokens.length; i++){
			if(tokens[i].equals("+") || tokens[i].equals("-")){
				return Strand.fromString(tokens[i]);
			}
		}
		return Strand.UNKNOWN;
	}

	//TODO Original parsing method
	/*public static Cluster parseCluster(String line){
		String[] tokens=line.split("\t");
		String barcode=tokens[0];
		Cluster cluster=new Cluster(barcode);
		cluster.setLine(line);
		
		//TODO SPLIT ON PARENTHESES
		for(int i=1; i<tokens.length; i++){
			//System.err.println(tokens[i]);
			if(tokens[i].contains(":")){
				String name=splitOnParentheses(tokens[i]);
				boolean dna=isDNA(tokens[i]);
				tokens[i]=tokens[i].split("\\)_")[1];
				String chr=tokens[i].split(":")[0];
				int start=new Integer(tokens[i].split(":")[1].split("-")[0]);
				int end=new Integer(tokens[i].split(":")[1].split("-")[1]);
				SingleInterval interval=new SingleInterval(chr, start, end);
				if(dna){cluster.addDNARead(interval);}
				else{
					Pair<String> rnaName=getName(name); //TODO this needs to be fixed to get multiple names
					//String type=getType(name);
					RNAInterval rnaInterval=new RNAInterval(interval);
					rnaInterval.addNameType(rnaName.getValue1(), rnaName.getValue2());
					cluster.addRNARead(rnaInterval);
				}
			}
		}
		return cluster;
	}*/
	
	private static String getType(String token) {
		String internal=token.split("\\[")[1].split("\\]")[0];
		String type=internal.split(";")[2].split("\\.")[1];
		return type;
	}

	private static SingleInterval getCoordinates(String token) {
		
		
		//System.err.println(token);
		//NYBot21_Stg.Odd2Bo82.Even2Bo26.Odd2Bo64.HL522DSXX       DPM(129S1_SvImJ;+;Unassigned_NoFeatures.none;Pde6a.intron;B4A,B4A_dup103463,B4,SINE.repeat)_chr18:61228228-61228368     DPM(CAST_EiJ;-;Clstn2.intron;Unassigned_NoFeatures.none)_chr9:97688696-97688835 DPM(CAST_EiJ;+;U

		
		
		/*String internal=token.split("\\[")[1].split("\\]")[0];
		System.err.println("token "+token);
		System.err.println("internal "+internal);
		String strand=internal.split(";")[1];
		String name=internal.split(";")[2].split("\\.")[0];*/
		
		//System.err.println(token);
		
		String delim="\\)";
		if(token.contains("]")) {delim="]";}
		
		//String delim="_";
		
		//String external=token.split("\\)_")[1];
		String external=token.split(delim+"_")[1];
		
		String[] tokens=external.split(":");
		String startEnd=tokens[tokens.length-1];
		String[] split=startEnd.split("-");
		int start=Integer.parseInt(split[0]);
		int end=start+1;
		if(split.length>1) {
			end=Integer.parseInt(split[1]);
		}
		String chr=external.replace(":"+startEnd, "");
		SingleInterval rtrn=new SingleInterval(chr, start, end);
		
		//SingleInterval rtrn=new SingleInterval(external);
		/*rtrn.setName(name);
		rtrn.setOrientation(Strand.fromString(strand));*/
		return rtrn;
	}

	private static String splitOnParentheses(String string) {
		return string.split("\\(")[1].split("\\)")[0];
	}

	private static boolean isDNA(String string) {
		if(string.startsWith("DPM") || string.startsWith("DNA")){return true;}
		return false;
	}
	
	private static Pair<String> getName(String string) {
		String name=string;
		//String name=string.split("\\(")[1].replaceAll("\\)", "");
		//System.err.println(name);
		String[] tokens=name.split(";");
		Collection<Pair<String>> list=new ArrayList<Pair<String>>();
		for(int i=0; i<tokens.length; i++){
			if(tokens[i].contains(".")){
				//System.err.println(tokens[i]);
				String feature=tokens[i].split("\\.")[1];
				//System.err.println(feature);
				if(!feature.equalsIgnoreCase("none")){
					//System.err.println(feature);
					String name1=tokens[i].split("\\.")[0];
					String type=tokens[i].split("\\.")[1];
					Pair<String> nameType=new Pair<String>(name1, type);
					list.add(nameType);
				}
			}
			else{
				//System.err.println(tokens[i]);
				Pair<String> pair=new Pair<String>(tokens[i], "none");
				list.add(pair);
			}
		}
		
		if(list.size()==0){return new Pair<String>("unassigned", "none");}
		if(list.size()==1){return list.iterator().next();}
		
		//System.err.println(name+" "+list.size());
		Pair<String> rtrn=collapse(list);
		return rtrn;
		//return list.iterator().next();
		//return name;
	}
	
	private static Pair<String> collapse(Collection<Pair<String>> list) {
		Pair<String> ambigous= new Pair<String>("unassigned", "none");
		
		Pair<String> current=list.iterator().next();
		String name=current.getValue1();
		String type=current.getValue2();
		
		for(Pair<String> pair: list){
			if(!name.equalsIgnoreCase(pair.getValue1())){return ambigous;}
			if(!type.equalsIgnoreCase(pair.getValue2())){type="ambigous";}
		}
		
		return new Pair<String>(name, type);
	}

	public String getLine() {
		return this.line;
	}

	public boolean hasDNA() {
		return !this.dnaRegions.isEmpty();
	}
	
	public boolean hasRNA() {
		return !this.rnaRegions.isEmpty();
	}

	public Cluster subsetRNA(Kmer kmer) {
		Cluster newCluster=new Cluster(this.getBarcode());
		
		for(RNAInterval rna: this.getAllRNARegions()){
			if(kmer.getRegions().contains(rna.getName())){newCluster.addRNARead(rna);}
		}
		
		
		return newCluster;
	}
	
	public Cluster subsetRNA(Kmer kmer, Collection<String> mRNA) {
		Cluster newCluster=new Cluster(this.getBarcode());
		
		for(RNAInterval rna: this.getAllRNARegions()){
			if(kmer.getRegions().contains(rna.getName())){newCluster.addRNARead(rna);}
			else if(mRNA.contains(rna.getName())){
				//System.err.println(rna.getName());
				if(rna.isExon()){
					RNAInterval exon=new RNAInterval(rna);
					exon.addNameType(rna.getName()+".exon", "exon");
					//exon.setName(rna.getName()+".exon");
					newCluster.addRNARead(exon);
				}
				else if(rna.isIntron()){
					RNAInterval intron=new RNAInterval(rna);
					intron.addNameType(rna.getName()+".intron", "intron");
					//intron.setName(rna.getName()+".intron");
					newCluster.addRNARead(intron);
				}
			}
		}
		
		
		return newCluster;
	}
	
	public Cluster subsetRNA(Kmer kmer, Collection<String> intron, Collection<String> exon) {
		Cluster newCluster=new Cluster(this.getBarcode());
		
		for(RNAInterval rna: this.getAllRNARegions()){
			if(kmer.getRegions().contains(rna.getName())){newCluster.addRNARead(rna);}
			if(exon.contains(rna.getName())){
				//System.err.println(rna.getName());
				if(rna.isExon()){
					RNAInterval e=new RNAInterval(rna);
					e.addNameType(rna.getName()+".exon", "exon");
					//exon.setName(rna.getName()+".exon");
					newCluster.addRNARead(e);
				}
			}
			if(intron.contains(rna.getName())){
				if(rna.isIntron()){
					RNAInterval i=new RNAInterval(rna);
					i.addNameType(rna.getName()+".intron", "intron");
					//intron.setName(rna.getName()+".intron");
					newCluster.addRNARead(i);
				}
			}
		}
		
		
		return newCluster;
	}

	public static Map<SingleInterval, Integer> generateCounts(Collection<Cluster> clusters) {
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		
		for(Cluster cluster: clusters){
			Collection<SingleInterval> regions=cluster.getAllDNAIntervals();
			for(SingleInterval region: regions){
				for(int start=region.getReferenceStartPosition(); start<region.getReferenceEndPosition(); start++){
					SingleInterval temp=new SingleInterval(region.getReferenceName(), start, start+1);
					int count=0;
					if(rtrn.containsKey(temp)){count=rtrn.get(temp);}
					count++;
					rtrn.put(temp, count);
				}
			}
		}
		
		return rtrn;
	}

	public void addName(String name) {
		this.rnaList.add(name);
		
	}

	public List<String> getRNANameList() {
		return this.rnaList;
	}

	public Cluster renameRNA(Map<String, String> names) {
		Cluster rtrn=new Cluster(this.barcode);
		rtrn.addDNAReads(getAllDNAIntervals());
		
		for(RNAInterval rna: this.getAllRNARegions()) {
			RNAInterval newRNA=rna;
			if(names.containsKey(rna.getName())) {
				newRNA=new RNAInterval(rna);
				newRNA.setName(names.get(rna.getName()));
				//System.err.println(rna.getName()+"\t"+newRNA.getName()+"\t"+names.get(rna.getName()));
			}
			rtrn.addRNARead(newRNA);
		}
		return rtrn;
	}

	public void addProtein(String protein) {
		this.proteins.add(protein);
		
	}

	public void addProteins(Collection<String> proteins2) {
		this.proteins.addAll(proteins2);
		
	}

	public Collection<String> getProteinSet() {
		Collection<String> rtrn=new TreeSet<String>();
		
		rtrn.addAll(this.getProteins());
		
		return rtrn;
	}

	public Map<String, Double> getProteinWeights() {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		double increment=1.0/(double)this.getProteins().size();
		
		for(String protein: getProteins()) {
			double count=0;
			if(rtrn.containsKey(protein)) {
				count=rtrn.get(protein);
			}
			count+=increment;
			rtrn.put(protein, count);
		}
		
		return rtrn;
	}

	

	
	

	
	
}
