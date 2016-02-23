package guttmanlab.core.datastructures;

/**
 * A class which represents a pair of objects with the same type. Either or both of the elements
 * of this Pair may be null.
 *
 * @param <T> is the type of the elements contained in the Pair
 */
public class Pair<T> {

	private T value1;
	private T value2;
	
	public Pair(){}
	
	public Pair(T v1, T v2) {
		value1 = v1;
		value2 = v2;
	}
	
	/**
	 * Creates a new Pair from two objects, inferring the generic types.
	 * @param the first element of this Pair
	 * @param the second element of this Pair
	 * @return a Pair containing v1 and v2
	 */
	public static <T> Pair<T> of(T v1, T v2) {
		return new Pair<T>(v1, v2);
	}
	
	public void setValue1(T v1) {
		value1 = v1;
	}
	
	public void setValue2(T v2){
		value2 = v2;
	}
	
	public T getValue1() {
		return value1;
	}
	
	public T getValue2() {
		return value2;
	}

	/**
	 * Returns true if this Pair's first value is not null; otherwise, returns false.
	 * @return if this Pair has a first value
	 */	
	public boolean hasValue1() {
		return value1 != null;
	}

	/**
	 * Returns true if this Pair's second value is not null; otherwise, returns false.
	 * @return if this Pair has a second value
	 */	
	public boolean hasValue2() {
		return value2 != null;
	}
	
	/**
	 * Returns whether this Pair is empty, i.e., whether this Pair has null values
	 * for both of its elements.
	 * @return if this Pair is empty
	 */
	public boolean isEmpty() {
		return value1 == null && value2 == null;
	}

	/**
	 * Returns whether this Pair is complete, i.e., whether this Pair has non-null
	 * values for both of its elements.
	 * @return if this Pair is complete
	 */
	public boolean isComplete() {
		return value1 != null && value2 != null;
	}
	
	/**
	 * {@inheritDoc}
	 * 
	 * Two Pairs are considered equal if their first elements are equal and if their
	 * second elements are equal. Since null values within the Pair may be meaningful,
	 * two elements are considered equal if both are null. Note that if the Pair itself
	 * is null (as opposed to one of its elements), this method will return false.
	 */
	@Override
	public boolean equals(Object other) {

		// No need for null check. The instanceof operator returns false if (other == null).
		if (!(other instanceof Pair)) {
			return false;
		}

		// OK to cast this. Class was explicitly checked above
		@SuppressWarnings("unchecked")
		Pair<T> o = (Pair<T>)other;
		
		boolean cond1 = (value1 == null ? o.value1 == null : value1.equals(o.value1));
		boolean cond2 = (value2 == null ? o.value2 == null : value2.equals(o.value2));
		return cond1 && cond2;
	}
	
	@Override
	public int hashCode() {
		int hashCode = 17;
		hashCode = value1 == null ? 31 * hashCode : 31 * hashCode + value1.hashCode();
		hashCode = value2 == null ? 31 * hashCode : 31 * hashCode + value2.hashCode();
		return hashCode;
	}
	
	/**
	 *  Returns a String representation of this Pair by calling toString() on its members. The
	 *  exact details of this representation are unspecified and subject to change, but the
	 *  following may be regarded as typical:
	 *  
	 *  "(value1, value2)"
	 */
	@Override
	public String toString() {
		if (isComplete()) {
			return "(" + value1.toString() + ", " + value2.toString() + ")";
		}
		if (isEmpty()) {
			return "(null, null)";
		}
		if (!hasValue1()) {
			return "(null, " + value2.toString() + ")";
		}
		// else value2 == null
		return "(" + value1.toString() + ", null)";
	}
}