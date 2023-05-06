package guttmanlab.core.rnasprite.hubs;


public class Score implements Comparable<Score>{
	double score;
	String name;
		
	public Score(String name, double score) {
		this.score=score;
		this.name=name;
	}
		
	@Override
	public int compareTo(Score o) {
		if(score!=o.score) {
			return Double.valueOf(o.score).compareTo(Double.valueOf(score));
		}
			
		return o.name.compareTo(name);
	}
		
	

}
