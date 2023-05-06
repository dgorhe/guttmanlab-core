package guttmanlab.core.annotation;

import org.apache.commons.lang.builder.HashCodeBuilder;

/**
 * A class for annotations that need to be distinguishable by name even when they have the same structure
 * @author prussell
 *
 */
public abstract class AbstractNamedAnnotation extends BlockedAnnotation {

	public boolean equals(Object o) {
		if(o instanceof Annotation) {
			String otherName = ((Annotation)o).getName();
			if(!otherName.equals(getName())) {
				return false;
			}
			return equals(o);
		}
		return false;
	}
	
	public int hashCode() {
		HashCodeBuilder h = hashCodeBuilder();
		h.append(getName());
		return h.toHashCode();
	}
	
}
