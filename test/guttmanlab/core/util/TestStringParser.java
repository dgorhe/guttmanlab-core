package guttmanlab.core.util;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.Test;

public class TestStringParser {
	
	private StringParser sp;
	private static final String empty = "";
	private static final String ints = "1 2 3";
	private static final String fruits = "apple orange banana";
	private static final String diffSpaces = "a b  c   d";
	private static final String tabs = "a\tb\tc";
	private static final String tabSpace = "a b\tc";
	
	
	@Before
	public void initialize() {
		sp = new StringParser();
	}
	
	@Test
	public void testInt() {
		sp.parse(ints);
		assertEquals("Field 0 should be an int: 1", sp.asInt(0), 1);
	}
	
}
