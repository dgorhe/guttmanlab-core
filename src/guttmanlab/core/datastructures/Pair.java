package guttmanlab.core.datastructures;

public class Pair<T1>{

	T1 value1;
	T1 value2;
	String name;
	
	/**
	 * An empty constructor
	 */
	public Pair(){}
	
	public Pair(T1 v1, T1 v2){
		this.value1=v1;
		this.value2=v2;
	}
	
	public void setValue1(T1 v1){this.value1=v1;}
	public void setValue2(T1 v2){this.value2=v2;}
	
	public T1 getValue1(){return value1;}
	public T1 getValue2(){return value2;}
	
	public boolean hasValue2(){
		if(value2==null){return false;}
		return true;
	}

	public boolean hasValue1(){
		if(value1==null){return false;}
		return true;
	}
	
	public boolean isEmpty() {
		if(value1==null && value2==null){return true;}
		return false;
	}
	
	public boolean isComplete() {
		if(value1!=null && value2!=null){return true;}
		return false;
	}
	
	public boolean equals(Pair<T1> other){
		if(other.value1.equals(value1) && other.value2.equals(value2)){return true;}
		return false;
	}

	
	public String toString(){
		return getValue1().toString()+","+getValue2().toString();
	}

	public boolean equals(Object other){
		if(!other.getClass().equals(Pair.class)) return false;
		@SuppressWarnings("unchecked")
		Pair<T1> t=(Pair<T1>)other;
		return equals(t);
	}
	
	public int hashCode() {
		if(value1==null){return value2.hashCode();}
		if(value2==null){return value1.hashCode();}
		String h = Integer.valueOf(value1.hashCode()).toString() + "_" + Integer.valueOf(value2.hashCode()).toString();
		return h.hashCode();
	}

	public void setName(String name) {
		this.name=name;
	}
	
	public String getName(){return name;}
}
