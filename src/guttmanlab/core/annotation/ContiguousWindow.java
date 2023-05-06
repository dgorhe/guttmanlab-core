package guttmanlab.core.annotation;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;

/**
 * A window that is a single contiguous interval
 * @author mguttman
 *
 */
public class ContiguousWindow<T extends Annotation> extends SingleInterval implements PopulatedWindow<T>{

	Collection<T> annotations;
	
	public ContiguousWindow(String refName, int start, int end) {
		this(refName, start, end, Strand.BOTH);
	}
	
	public ContiguousWindow(String refName, int start, int end, Strand orientation){
		super(refName, start, end, orientation);
		this.annotations=new ArrayList<T>();
	}

	@Override
	public Annotation getParentAnnotation() {
		// FIXME Auto-generated method stub
		throw new UnsupportedOperationException("TODO");
	}

	@Override
	public Score getScore() {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException("TODO");
	}

	@Override
	public void addAnnotation(T annotation) {
		this.annotations.add(annotation);
	}

	@Override
	public int getNumberOfAnnotationsInWindow() {
		return this.annotations.size();
	}

	@Override
	public Iterator<T> getAnnotationsInWindow() {
		return this.annotations.iterator();
	}

	@Override
	public int getNumberOfAnnotationsInWindow(Strand orientation) {
		int count=0;
		for(T annotation: annotations){
			if(annotation.getOrientation().equals(orientation)){count++;}
		}
		return count;
	}

}
