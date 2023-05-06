package guttmanlab.core.nanobody;

import java.util.Collection;
import java.util.Iterator;
import java.util.TreeSet;


public class Cluster implements Comparable<Cluster>{
		
		Collection<String> proteins;
		String name;
		int score;
		
		public Cluster(Collection<String> relatedProteins) {
			this.proteins=new TreeSet<String>();
			this.proteins.addAll(relatedProteins);
		}
		
		public Cluster(String relatedProtein) {
			this.proteins=new TreeSet<String>();
			this.proteins.add(relatedProtein);
		}
		
		public void setScore(int score) {
			this.score=score;
		}
		
		
		@Override
		public int compareTo(Cluster o) {
			if(this.score!=o.score) {
				return o.score-this.score;
			}
			
			if(this.proteins.size()!=o.proteins.size()) {
				return o.proteins.size()-this.proteins.size();
			}
			
			Iterator<String> iter1=this.proteins.iterator();
			Iterator<String> iter2=o.proteins.iterator();
				
			while(iter1.hasNext()) {
				String p1=iter1.next();
				String p2=iter2.next();
				if(!p1.equals(p2)) {return p1.compareTo(p2);}
			}
			return 0;
		}
		
		@Override
		public String toString() {
			String rtrn=""+proteins.size();
			for(String protein: proteins) {
				rtrn+="_"+protein;
			}
			return rtrn;
		}

		public void setName(String string) {
			this.name=string;
		}

		public String getName() {
			return this.name;
		}

		public int getScore() {
			return this.score;
		}

		public boolean isSubset(Cluster c) {
			if(this.proteins.size()<c.proteins.size()) {
				for(String p: this.proteins) {
					if(!c.proteins.contains(p)) {return false;}
				}
				return true;
			}
			return false;
		}

		public void addProtein(String str2) {
			this.proteins.add(str2);
		}

		public int getSize() {
			return this.proteins.size();
		}

		public Collection<String> getProteins() {
			return this.proteins;
		}

		public void addProteins(Collection<String> proteins2) {
			for(String p: proteins2) {addProtein(p);}
			
		}

	
}
