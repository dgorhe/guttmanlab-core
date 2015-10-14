package guttmanlab.core.annotation;

/**
 * Represents an annotation that was derived from another Annotation
 * @author mguttman
 *
 * @param <T> The type of the Annotation it was derived from
 */
public class DerivedAnnotation<T extends Annotation> extends BlockedAnnotation{

	private T parentDerivedFrom;
	
	public DerivedAnnotation(Annotation annotation, T parent){
		super(annotation);
		this.parentDerivedFrom=parent;
	}
	
	public DerivedAnnotation(Annotation annotation, String newName, T parent){
		super(annotation, newName);
		this.parentDerivedFrom=parent;
	}
	
	/**
	 * The parent annotation this was derived from
	 * @return A pointer to the parent this was derived from
	 */
	public T getParentAnnotation(){return parentDerivedFrom;}
	
	public String getParentName()
	{
		return parentDerivedFrom.getName();
	}
	
	public Strand getOrientation()
	{
		return parentDerivedFrom.getOrientation();
	}
	
	public DerivedAnnotation<T> merge(Annotation other)
	{
		Annotation merged = super.merge(other);
		return new DerivedAnnotation<T>(merged,this.parentDerivedFrom);
	}
}
