package guttmanlab.core.annotation;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;

import guttmanlab.core.annotation.Annotation.Strand;

public class BlockedWindow<T extends Annotation> extends BlockedAnnotation implements PopulatedWindow<T>{

	Collection<T> annotations;
	Annotation parent;
	
	public BlockedWindow(DerivedAnnotation<? extends Annotation> a) {
		super(a);
		this.parent=a.getParentAnnotation();
		this.annotations=new ArrayList<T>();
	}
	
	@Override
	public Annotation getParentAnnotation() {
		return this.parent;
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


	@Override
	public SingleInterval getMidPoint() {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException("TODO");
	}

	
}
