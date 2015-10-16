package guttmanlab.core.annotationcollection;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.ContiguousWindow;
import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.PairedMappedFragment;
import guttmanlab.core.annotation.PopulatedWindow;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.coordinatespace.CoordinateSpace;
import guttmanlab.core.datastructures.IntervalTree;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.collections15.Predicate;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.util.CloseableIterator;

public abstract class AbstractAnnotationCollection<T extends Annotation> implements AnnotationCollection<T>{

	private Collection<Predicate<T>> filters;
	private int numAnnotations;

	public AbstractAnnotationCollection(){
		filters=new ArrayList<Predicate<T>>();
	}

	@Override
	public void addFilter(Predicate<T> filter) {
		filters.add(filter);
		numAnnotations = 0;
	}
	
	public void addFilter(Collection<Predicate<PairedMappedFragment<SAMFragment>>> filters) {
		Iterator<Predicate<PairedMappedFragment<SAMFragment>>> iter = filters.iterator();
		while(iter.hasNext())
			filters.add(iter.next());
	}

	@Override
	public Collection<Predicate<T>> getFilters(){
		return filters;
	}

	@Override
	public FeatureCollection<BlockedAnnotation> merge() {
		//if(stranded)
			//Graph plusExons = new Graph();
			//Graph minusExons = new Graph();
		//else
			//Graph exons = new Graph();
		Graph<SingleInterval> graph = new Graph<SingleInterval>();
		//TODO: should this be in blocked annotation?
		CloseableIterator<T> iter = this.sortedIterator();
		while ( iter.hasNext())
		{
			T next = iter.next();
			Iterator<SingleInterval> blocks = next.getBlocks();
			
			SingleInterval block = blocks.next();
			String name = block.getReferenceName()+":"+block.getReferenceStartPosition()+"-"+block.getReferenceEndPosition();

			Vertex<SingleInterval> from = null;
			if(block != null)
			{
				from = graph.addVertex(new Vertex<SingleInterval>(name,block));
			}
			while ( blocks.hasNext())
			{
				block = blocks.next();
				name = block.getReferenceName()+":"+block.getReferenceStartPosition()+"-"+block.getReferenceEndPosition();
				Vertex<SingleInterval> to = graph.addVertex(new Vertex<SingleInterval>(name,block));
				if(!from.hasEdge(to))
					from.addOutgoingEdge(to,0);
				from = to;
			}
		}
		
		/*
		List<Vertex<SingleInterval>> vertices = graph.getVerticies();
		for (Vertex<SingleInterval> v : vertices)
		{
			System.out.println(v);

			List<Edge> edges = v.getOutgoingEdges();
			for(Edge e : edges)
			{
				System.out.println(e);
			}
		}*/
		
		FeatureCollection<BlockedAnnotation> merged = new FeatureCollection<BlockedAnnotation>(getReferenceCoordinateSpace());
		List<Vertex<SingleInterval>> vertices = graph.getVerticies();
		for (Vertex<SingleInterval> v : vertices)
		{
			for( BlockedAnnotation b : graph.allPaths(v))
			{
				//System.out.println(b.toBED());
				merged.addAnnotation(b);
			}
		}
		
		return merged;
	}

	 public static class Graph<T> {
		 
		 public List<BlockedAnnotation> allPaths(Vertex<SingleInterval> v)
		 {
			 ArrayList<BlockedAnnotation> result = new ArrayList<BlockedAnnotation>();
			 List<Edge> edges = v.getOutgoingEdges();
			 boolean added = false;

			 for (Edge e : edges)
			 {
				 Vertex<SingleInterval> child = e.getTo();
				 List<BlockedAnnotation> subpaths = allPaths(child);
				 
				 for (BlockedAnnotation path : subpaths)
				 {
					 path.addBlocks(v.getData());
					 result.add(path);
					 added = true;
				 }
				 
			 }
			 if(!added)
			 {
				 BlockedAnnotation b = new BlockedAnnotation();
				 b.addBlocks(v.getData());
				 result.add(b);
			 }
			 return result;
		 }
		 
		  /** Color used to mark unvisited nodes */
		  public static final int VISIT_COLOR_WHITE = 1;

		  /** Color used to mark nodes as they are first visited in DFS order */
		  public static final int VISIT_COLOR_GREY = 2;

		  /** Color used to mark nodes after descendants are completely visited */
		  public static final int VISIT_COLOR_BLACK = 3;

		  /** Vector<Vertex> of graph verticies */
		  private List<Vertex<T>> verticies;

		  /** Vector<Edge> of edges in the graph */
		  private List<Edge<T>> edges;

		  /** The vertex identified as the root of the graph */
		  private Vertex<T> rootVertex;

		  /**
		   * Construct a new graph without any vertices or edges
		   */
		  public Graph() {
		    verticies = new ArrayList<Vertex<T>>();
		    edges = new ArrayList<Edge<T>>();
		  }

		  /**
		   * Are there any verticies in the graph
		   * 
		   * @return true if there are no verticies in the graph
		   */
		  public boolean isEmpty() {
		    return verticies.size() == 0;
		  }

		  /**
		   * Add a vertex to the graph
		   * 
		   * @param v
		   *          the Vertex to add
		   * @return true if the vertex was added, false if it was already in the graph.
		   */
		  public Vertex<T> addVertex(Vertex<T> v) {
		    Vertex<T> added = hasCompatible(v);
		    if (hasCompatible(v) == null) {
		      verticies.add(v);
		      return v;
		    }
		    return added;
		  }
		  
		  //TODO: this ought to actually use isCompatible()
		  private Vertex<T> hasCompatible(Vertex<T> v){
			  //fix this casting
			  SingleInterval e1 = (SingleInterval)v.getData();
			  for ( Vertex<T> next : verticies)
			  {
				  SingleInterval e2 = (SingleInterval)next.getData();
				  if (e1.getReferenceStartPosition()==e2.getReferenceStartPosition() &&
						  e1.getReferenceEndPosition() == e2.getReferenceEndPosition())
					  return next;
			  }
			  return null;
		  }

		  /**
		   * Get the vertex count.
		   * 
		   * @return the number of verticies in the graph.
		   */
		  public int size() {
		    return verticies.size();
		  }

		  /**
		   * Get the root vertex
		   * 
		   * @return the root vertex if one is set, null if no vertex has been set as
		   *         the root.
		   */
		  public Vertex<T> getRootVertex() {
		    return rootVertex;
		  }

		  /**
		   * Set a root vertex. If root does no exist in the graph it is added.
		   * 
		   * @param root -
		   *          the vertex to set as the root and optionally add if it does not
		   *          exist in the graph.
		   */
		  public void setRootVertex(Vertex<T> root) {
		    this.rootVertex = root;
		    if (verticies.contains(root) == false)
		      this.addVertex(root);
		  }

		  /**
		   * Get the given Vertex.
		   * 
		   * @param n
		   *          the index [0, size()-1] of the Vertex to access
		   * @return the nth Vertex
		   */
		  public Vertex<T> getVertex(int n) {
		    return verticies.get(n);
		  }

		  /**
		   * Get the graph verticies
		   * 
		   * @return the graph verticies
		   */
		  public List<Vertex<T>> getVerticies() {
		    return this.verticies;
		  }

		  /**
		   * Insert a directed, weighted Edge<T> into the graph.
		   * 
		   * @param from -
		   *          the Edge<T> starting vertex
		   * @param to -
		   *          the Edge<T> ending vertex
		   * @param cost -
		   *          the Edge<T> weight/cost
		   * @return true if the Edge<T> was added, false if from already has this Edge<T>
		   * @throws IllegalArgumentException
		   *           if from/to are not verticies in the graph
		   */
		  public boolean addEdge(Vertex<T> from, Vertex<T> to, int cost) throws IllegalArgumentException {
		    if (verticies.contains(from) == false)
		      throw new IllegalArgumentException("from is not in graph");
		    if (verticies.contains(to) == false)
		      throw new IllegalArgumentException("to is not in graph");

		    Edge<T> e = new Edge<T>(from, to, cost);
		    if (from.findEdge(to) != null)
		      return false;
		    else {
		      from.addEdge(e);
		      to.addEdge(e);
		      edges.add(e);
		      return true;
		    }
		  }

		  /**
		   * Insert a bidirectional Edge<T> in the graph
		   * 
		   * @param from -
		   *          the Edge<T> starting vertex
		   * @param to -
		   *          the Edge<T> ending vertex
		   * @param cost -
		   *          the Edge<T> weight/cost
		   * @return true if edges between both nodes were added, false otherwise
		   * @throws IllegalArgumentException
		   *           if from/to are not verticies in the graph
		   */
		  public boolean insertBiEdge(Vertex<T> from, Vertex<T> to, int cost)
		      throws IllegalArgumentException {
		    return addEdge(from, to, cost) && addEdge(to, from, cost);
		  }

		  /**
		   * Get the graph edges
		   * 
		   * @return the graph edges
		   */
		  public List<Edge<T>> getEdges() {
		    return this.edges;
		  }

		  /**
		   * Remove a vertex from the graph
		   * 
		   * @param v
		   *          the Vertex to remove
		   * @return true if the Vertex was removed
		   */
		  public boolean removeVertex(Vertex<T> v) {
		    if (!verticies.contains(v))
		      return false;

		    verticies.remove(v);
		    if (v == rootVertex)
		      rootVertex = null;

		    // Remove the edges associated with v
		    for (int n = 0; n < v.getOutgoingEdgeCount(); n++) {
		      Edge<T> e = v.getOutgoingEdge(n);
		      v.remove(e);
		      Vertex<T> to = e.getTo();
		      to.remove(e);
		      edges.remove(e);
		    }
		    for (int n = 0; n < v.getIncomingEdgeCount(); n++) {
		      Edge<T> e = v.getIncomingEdge(n);
		      v.remove(e);
		      Vertex<T> predecessor = e.getFrom();
		      predecessor.remove(e);
		    }
		    return true;
		  }

		  /**
		   * Remove an Edge<T> from the graph
		   * 
		   * @param from -
		   *          the Edge<T> starting vertex
		   * @param to -
		   *          the Edge<T> ending vertex
		   * @return true if the Edge<T> exists, false otherwise
		   */
		  public boolean removeEdge(Vertex<T> from, Vertex<T> to) {
		    Edge<T> e = from.findEdge(to);
		    if (e == null)
		      return false;
		    else {
		      from.remove(e);
		      to.remove(e);
		      edges.remove(e);
		      return true;
		    }
		  }

		  /**
		   * Clear the mark state of all verticies in the graph by calling clearMark()
		   * on all verticies.
		   * 
		   * @see Vertex#clearMark()
		   */
		  public void clearMark() {
		    for (Vertex<T> w : verticies)
		      w.clearMark();
		  }

		  /**
		   * Clear the mark state of all edges in the graph by calling clearMark() on
		   * all edges.
		   */
		  public void clearEdges() {
		    for (Edge<T> e : edges)
		      e.clearMark();
		  }

		  /**
		   * Perform a depth first serach using recursion.
		   * 
		   * @param v -
		   *          the Vertex to start the search from
		   * @param visitor -
		   *          the vistor to inform prior to
		   * @see Visitor#visit(Graph, Vertex)
		   */
		  public void depthFirstSearch(Vertex<T> v, final Visitor<T> visitor) {
		    VisitorEX<T, RuntimeException> wrapper = new VisitorEX<T, RuntimeException>() {
		      public void visit(Graph<T> g, Vertex<T> v) throws RuntimeException {
		        if (visitor != null)
		          visitor.visit(g, v);
		      }
		    };
		    this.depthFirstSearch(v, wrapper);
		  }

		  /**
		   * Perform a depth first serach using recursion. The search may be cut short
		   * if the visitor throws an exception.
		   * 
		   * @param <E>
		   * 
		   * @param v -
		   *          the Vertex to start the search from
		   * @param visitor -
		   *          the vistor to inform prior to
		   * @see Visitor#visit(Graph, Vertex)
		   * @throws E
		   *           if visitor.visit throws an exception
		   */
		  public <E extends Exception> void depthFirstSearch(Vertex<T> v, VisitorEX<T, E> visitor) throws E {
		    if (visitor != null)
		      visitor.visit(this, v);
		    v.visit();
		    for (int i = 0; i < v.getOutgoingEdgeCount(); i++) {
		      Edge<T> e = v.getOutgoingEdge(i);
		      if (!e.getTo().visited()) {
		        depthFirstSearch(e.getTo(), visitor);
		      }
		    }
		  }

		  /**
		   * Perform a breadth first search of this graph, starting at v.
		   * 
		   * @param v -
		   *          the search starting point
		   * @param visitor -
		   *          the vistor whose vist method is called prior to visting a vertex.
		   */
		  public void breadthFirstSearch(Vertex<T> v, final Visitor<T> visitor) {
		    VisitorEX<T, RuntimeException> wrapper = new VisitorEX<T, RuntimeException>() {
		      public void visit(Graph<T> g, Vertex<T> v) throws RuntimeException {
		        if (visitor != null)
		          visitor.visit(g, v);
		      }
		    };
		    this.breadthFirstSearch(v, wrapper);
		  }

		  /**
		   * Perform a breadth first search of this graph, starting at v. The vist may
		   * be cut short if visitor throws an exception during a vist callback.
		   * 
		   * @param <E>
		   * 
		   * @param v -
		   *          the search starting point
		   * @param visitor -
		   *          the vistor whose vist method is called prior to visting a vertex.
		   * @throws E
		   *           if vistor.visit throws an exception
		   */
		  public <E extends Exception> void breadthFirstSearch(Vertex<T> v, VisitorEX<T, E> visitor)
		      throws E {
		    LinkedList<Vertex<T>> q = new LinkedList<Vertex<T>>();

		    q.add(v);
		    if (visitor != null)
		      visitor.visit(this, v);
		    v.visit();
		    while (q.isEmpty() == false) {
		      v = q.removeFirst();
		      for (int i = 0; i < v.getOutgoingEdgeCount(); i++) {
		        Edge<T> e = v.getOutgoingEdge(i);
		        Vertex<T> to = e.getTo();
		        if (!to.visited()) {
		          q.add(to);
		          if (visitor != null)
		            visitor.visit(this, to);
		          to.visit();
		        }
		      }
		    }
		  }

		  /**
		   * Find the spanning tree using a DFS starting from v.
		   * 
		   * @param v -
		   *          the vertex to start the search from
		   * @param visitor -
		   *          visitor invoked after each vertex is visited and an edge is added
		   *          to the tree.
		   */
		  public void dfsSpanningTree(Vertex<T> v, DFSVisitor<T> visitor) {
		    v.visit();
		    if (visitor != null)
		      visitor.visit(this, v);

		    for (int i = 0; i < v.getOutgoingEdgeCount(); i++) {
		      Edge<T> e = v.getOutgoingEdge(i);
		      if (!e.getTo().visited()) {
		        if (visitor != null)
		          visitor.visit(this, v, e);
		        e.mark();
		        dfsSpanningTree(e.getTo(), visitor);
		      }
		    }
		  }

		  /**
		   * Search the verticies for one with name.
		   * 
		   * @param name -
		   *          the vertex name
		   * @return the first vertex with a matching name, null if no matches are found
		   */
		  public Vertex<T> findVertexByName(String name) {
		    Vertex<T> match = null;
		    for (Vertex<T> v : verticies) {
		      if (name.equals(v.getName())) {
		        match = v;
		        break;
		      }
		    }
		    return match;
		  }

		  /**
		   * Search the verticies for one with data.
		   * 
		   * @param data -
		   *          the vertex data to match
		   * @param compare -
		   *          the comparator to perform the match
		   * @return the first vertex with a matching data, null if no matches are found
		   */
		  public Vertex<T> findVertexByData(T data, Comparator<T> compare) {
		    Vertex<T> match = null;
		    for (Vertex<T> v : verticies) {
		      if (compare.compare(data, v.getData()) == 0) {
		        match = v;
		        break;
		      }
		    }
		    return match;
		  }

		  /**
		   * Search the graph for cycles. In order to detect cycles, we use a modified
		   * depth first search called a colored DFS. All nodes are initially marked
		   * white. When a node is encountered, it is marked grey, and when its
		   * descendants are completely visited, it is marked black. If a grey node is
		   * ever encountered, then there is a cycle.
		   * 
		   * @return the edges that form cycles in the graph. The array will be empty if
		   *         there are no cycles.
		   */
		  public Edge<T>[] findCycles() {
		    ArrayList<Edge<T>> cycleEdges = new ArrayList<Edge<T>>();
		    // Mark all verticies as white
		    for (int n = 0; n < verticies.size(); n++) {
		      Vertex<T> v = getVertex(n);
		      v.setMarkState(VISIT_COLOR_WHITE);
		    }
		    for (int n = 0; n < verticies.size(); n++) {
		      Vertex<T> v = getVertex(n);
		      visit(v, cycleEdges);
		    }

		    Edge<T>[] cycles = new Edge[cycleEdges.size()];
		    cycleEdges.toArray(cycles);
		    return cycles;
		  }

		  private void visit(Vertex<T> v, ArrayList<Edge<T>> cycleEdges) {
		    v.setMarkState(VISIT_COLOR_GREY);
		    int count = v.getOutgoingEdgeCount();
		    for (int n = 0; n < count; n++) {
		      Edge<T> e = v.getOutgoingEdge(n);
		      Vertex<T> u = e.getTo();
		      if (u.getMarkState() == VISIT_COLOR_GREY) {
		        // A cycle Edge<T>
		        cycleEdges.add(e);
		      } else if (u.getMarkState() == VISIT_COLOR_WHITE) {
		        visit(u, cycleEdges);
		      }
		    }
		    v.setMarkState(VISIT_COLOR_BLACK);
		  }

		  public String toString() {
		    StringBuffer tmp = new StringBuffer("Graph[");
		    for (Vertex<T> v : verticies)
		      tmp.append(v);
		    tmp.append(']');
		    return tmp.toString();
		  }

		}

		/*
		 * JBoss, Home of Professional Open Source Copyright 2006, Red Hat Middleware
		 * LLC, and individual contributors by the @authors tag. See the copyright.txt
		 * in the distribution for a full listing of individual contributors.
		 * 
		 * This is free software; you can redistribute it and/or modify it under the
		 * terms of the GNU Lesser General Public License as published by the Free
		 * Software Foundation; either version 2.1 of the License, or (at your option)
		 * any later version.
		 * 
		 * This software is distributed in the hope that it will be useful, but WITHOUT
		 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
		 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
		 * details.
		 * 
		 * You should have received a copy of the GNU Lesser General Public License
		 * along with this software; if not, write to the Free Software Foundation,
		 * Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA, or see the FSF
		 * site: http://www.fsf.org.
		 */

		/**
		 * A directed, weighted edge in a graph
		 * 
		 * @author Scott.Stark@jboss.org
		 * @version $Revision$
		 * @param <T>
		 */
		static class Edge<T> {
		  private Vertex<T> from;

		  private Vertex<T> to;

		  private int cost;

		  private boolean mark;

		  /**
		   * Create a zero cost edge between from and to
		   * 
		   * @param from
		   *          the starting vertex
		   * @param to
		   *          the ending vertex
		   */
		  public Edge(Vertex<T> from, Vertex<T> to) {
		    this(from, to, 0);
		  }

		  /**
		   * Create an edge between from and to with the given cost.
		   * 
		   * @param from
		   *          the starting vertex
		   * @param to
		   *          the ending vertex
		   * @param cost
		   *          the cost of the edge
		   */
		  public Edge(Vertex<T> from, Vertex<T> to, int cost) {
		    this.from = from;
		    this.to = to;
		    this.cost = cost;
		    mark = false;
		  }

		  /**
		   * Get the ending vertex
		   * 
		   * @return ending vertex
		   */
		  public Vertex<T> getTo() {
		    return to;
		  }

		  /**
		   * Get the starting vertex
		   * 
		   * @return starting vertex
		   */
		  public Vertex<T> getFrom() {
		    return from;
		  }

		  /**
		   * Get the cost of the edge
		   * 
		   * @return cost of the edge
		   */
		  public int getCost() {
		    return cost;
		  }

		  /**
		   * Set the mark flag of the edge
		   * 
		   */
		  public void mark() {
		    mark = true;
		  }

		  /**
		   * Clear the edge mark flag
		   * 
		   */
		  public void clearMark() {
		    mark = false;
		  }

		  /**
		   * Get the edge mark flag
		   * 
		   * @return edge mark flag
		   */
		  public boolean isMarked() {
		    return mark;
		  }

		  /**
		   * String rep of edge
		   * 
		   * @return string rep with from/to vertex names and cost
		   */
		  public String toString() {
		    StringBuffer tmp = new StringBuffer("Edge[from: ");
		    tmp.append(from.getName());
		    tmp.append(",to: ");
		    tmp.append(to.getName());
		    tmp.append(", cost: ");
		    tmp.append(cost);
		    tmp.append("]");
		    return tmp.toString();
		  }
		}

		/*
		 * JBoss, Home of Professional Open Source Copyright 2006, Red Hat Middleware
		 * LLC, and individual contributors by the @authors tag. See the copyright.txt
		 * in the distribution for a full listing of individual contributors.
		 * 
		 * This is free software; you can redistribute it and/or modify it under the
		 * terms of the GNU Lesser General Public License as published by the Free
		 * Software Foundation; either version 2.1 of the License, or (at your option)
		 * any later version.
		 * 
		 * This software is distributed in the hope that it will be useful, but WITHOUT
		 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
		 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
		 * details.
		 * 
		 * You should have received a copy of the GNU Lesser General Public License
		 * along with this software; if not, write to the Free Software Foundation,
		 * Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA, or see the FSF
		 * site: http://www.fsf.org.
		 */

		/**
		 * A named graph vertex with optional data.
		 * 
		 * @author Scott.Stark@jboss.org
		 * @version $Revision$
		 * @param <T>
		 */
		@SuppressWarnings("unchecked")
		static class Vertex<T> {
		  private List<Edge<T>> incomingEdges;

		  private List<Edge<T>> outgoingEdges;

		  private String name;

		  private boolean mark;

		  private int markState;

		  private T data;
		  
		  /**
		   * Calls this(null, null).
		   */
		  public Vertex() {
		    this(null, null);
		  }

		  /**
		   * Create a vertex with the given name and no data
		   * 
		   * @param n
		   */
		  public Vertex(String n) {
		    this(n, null);
		  }

		  /**
		   * Create a Vertex with name n and given data
		   * 
		   * @param n -
		   *          name of vertex
		   * @param data -
		   *          data associated with vertex
		   */
		  public Vertex(String n, T data) {
		    incomingEdges = new ArrayList<Edge<T>>();
		    outgoingEdges = new ArrayList<Edge<T>>();
		    name = n;
		    mark = false;
		    this.data = data;
		  }

		  /**
		   * @return the possibly null name of the vertex
		   */
		  public String getName() {
		    return name;
		  }

		  /**
		   * @return the possibly null data of the vertex
		   */
		  public T getData() {
		    return this.data;
		  }

		  /**
		   * @param data
		   *          The data to set.
		   */
		  public void setData(T data) {
		    this.data = data;
		  }

		  /**
		   * Add an edge to the vertex. If edge.from is this vertex, its an outgoing
		   * edge. If edge.to is this vertex, its an incoming edge. If neither from or
		   * to is this vertex, the edge is not added.
		   * 
		   * @param e -
		   *          the edge to add
		   * @return true if the edge was added, false otherwise
		   */
		  public boolean addEdge(Edge<T> e) {
		    if (e.getFrom() == this)
		      outgoingEdges.add(e);
		    else if (e.getTo() == this)
		      incomingEdges.add(e);
		    else
		      return false;
		    return true;
		  }

		  /**
		   * Add an outgoing edge ending at to.
		   * 
		   * @param to -
		   *          the destination vertex
		   * @param cost
		   *          the edge cost
		   */
		  public void addOutgoingEdge(Vertex<T> to, int cost) {
		    Edge<T> out = new Edge<T>(this, to, cost);
		    outgoingEdges.add(out);
		  }

		  /**
		   * Add an incoming edge starting at from
		   * 
		   * @param from -
		   *          the starting vertex
		   * @param cost
		   *          the edge cost
		   */
		  public void addIncomingEdge(Vertex<T> from, int cost) {
		    Edge<T> out = new Edge<T>(this, from, cost);
		    incomingEdges.add(out);
		  }

		  /**
		   * Check the vertex for either an incoming or outgoing edge matching e.
		   * 
		   * @param e
		   *          the edge to check
		   * @return true it has an edge
		   */
		  public boolean hasEdge(Edge<T> e) {
		    if (e.getFrom() == this)
		      return incomingEdges.contains(e);
		    else if (e.getTo() == this)
		      return outgoingEdges.contains(e);
		    else
		      return false;
		  }

		  /**
		   * Remove an edge from this vertex
		   * 
		   * @param e -
		   *          the edge to remove
		   * @return true if the edge was removed, false if the edge was not connected
		   *         to this vertex
		   */
		  public boolean remove(Edge<T> e) {
		    if (e.getFrom() == this)
		      incomingEdges.remove(e);
		    else if (e.getTo() == this)
		      outgoingEdges.remove(e);
		    else
		      return false;
		    return true;
		  }

		  /**
		   * 
		   * @return the count of incoming edges
		   */
		  public int getIncomingEdgeCount() {
		    return incomingEdges.size();
		  }

		  /**
		   * Get the ith incoming edge
		   * 
		   * @param i
		   *          the index into incoming edges
		   * @return ith incoming edge
		   */
		  public Edge<T> getIncomingEdge(int i) {
		    return incomingEdges.get(i);
		  }

		  /**
		   * Get the incoming edges
		   * 
		   * @return incoming edge list
		   */
		  public List getIncomingEdges() {
		    return this.incomingEdges;
		  }

		  /**
		   * 
		   * @return the count of incoming edges
		   */
		  public int getOutgoingEdgeCount() {
		    return outgoingEdges.size();
		  }

		  /**
		   * Get the ith outgoing edge
		   * 
		   * @param i
		   *          the index into outgoing edges
		   * @return ith outgoing edge
		   */
		  public Edge<T> getOutgoingEdge(int i) {
		    return outgoingEdges.get(i);
		  }

		  /**
		   * Get the outgoing edges
		   * 
		   * @return outgoing edge list
		   */
		  public List getOutgoingEdges() {
		    return this.outgoingEdges;
		  }

		  /**
		   * Search the outgoing edges looking for an edge whose's edge.to == dest.
		   * 
		   * @param dest
		   *          the destination
		   * @return the outgoing edge going to dest if one exists, null otherwise.
		   */
		  public Edge<T> findEdge(Vertex<T> dest) {
		    for (Edge<T> e : outgoingEdges) {
		      if (e.getTo() == dest)
		        return e;
		    }
		    return null;
		  }

		  /**
		   * Search the outgoing edges for a match to e.
		   * 
		   * @param e -
		   *          the edge to check
		   * @return e if its a member of the outgoing edges, null otherwise.
		   */
		  public Edge<T> findEdge(Edge<T> e) {
		    if (outgoingEdges.contains(e))
		      return e;
		    else
		      return null;
		  }

		  /**
		   * What is the cost from this vertex to the destination vertex.
		   * 
		   * @param dest -
		   *          the destination vertex.
		   * @return Return Integer.MAX_VALUE if we have no edge to dest, 0 if dest is
		   *         this vertex, the cost of the outgoing edge otherwise.
		   */
		  public int cost(Vertex<T> dest) {
		    if (dest == this)
		      return 0;

		    Edge<T> e = findEdge(dest);
		    int cost = Integer.MAX_VALUE;
		    if (e != null)
		      cost = e.getCost();
		    return cost;
		  }

		  /**
		   * Is there an outgoing edge ending at dest.
		   * 
		   * @param dest -
		   *          the vertex to check
		   * @return true if there is an outgoing edge ending at vertex, false
		   *         otherwise.
		   */
		  public boolean hasEdge(Vertex<T> dest) {
		    return (findEdge(dest) != null);
		  }

		  /**
		   * Has this vertex been marked during a visit
		   * 
		   * @return true is visit has been called
		   */
		  public boolean visited() {
		    return mark;
		  }

		  /**
		   * Set the vertex mark flag.
		   * 
		   */
		  public void mark() {
		    mark = true;
		  }

		  /**
		   * Set the mark state to state.
		   * 
		   * @param state
		   *          the state
		   */
		  public void setMarkState(int state) {
		    markState = state;
		  }

		  /**
		   * Get the mark state value.
		   * 
		   * @return the mark state
		   */
		  public int getMarkState() {
		    return markState;
		  }

		  /**
		   * Visit the vertex and set the mark flag to true.
		   * 
		   */
		  public void visit() {
		    mark();
		  }

		  /**
		   * Clear the visited mark flag.
		   * 
		   */
		  public void clearMark() {
		    mark = false;
		  }

		  /**
		   * @return a string form of the vertex with in and out edges.
		   */
		  public String toString() {
		    StringBuffer tmp = new StringBuffer("Vertex(");
		    tmp.append(name);
		    tmp.append(", data=");
		    tmp.append(data);
		    tmp.append("), in:[");
		    for (int i = 0; i < incomingEdges.size(); i++) {
		      Edge<T> e = incomingEdges.get(i);
		      if (i > 0)
		        tmp.append(',');
		      tmp.append('{');
		      tmp.append(e.getFrom().name);
		      tmp.append(',');
		      tmp.append(e.getCost());
		      tmp.append('}');
		    }
		    tmp.append("], out:[");
		    for (int i = 0; i < outgoingEdges.size(); i++) {
		      Edge<T> e = outgoingEdges.get(i);
		      if (i > 0)
		        tmp.append(',');
		      tmp.append('{');
		      tmp.append(e.getTo().name);
		      tmp.append(',');
		      tmp.append(e.getCost());
		      tmp.append('}');
		    }
		    tmp.append(']');
		    return tmp.toString();
		  }
		}

		/*
		 * JBoss, Home of Professional Open Source Copyright 2006, Red Hat Middleware
		 * LLC, and individual contributors by the @authors tag. See the copyright.txt
		 * in the distribution for a full listing of individual contributors.
		 * 
		 * This is free software; you can redistribute it and/or modify it under the
		 * terms of the GNU Lesser General Public License as published by the Free
		 * Software Foundation; either version 2.1 of the License, or (at your option)
		 * any later version.
		 * 
		 * This software is distributed in the hope that it will be useful, but WITHOUT
		 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
		 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
		 * details.
		 * 
		 * You should have received a copy of the GNU Lesser General Public License
		 * along with this software; if not, write to the Free Software Foundation,
		 * Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA, or see the FSF
		 * site: http://www.fsf.org.
		 */

		/**
		 * A graph visitor interface.
		 * 
		 * @author Scott.Stark@jboss.org
		 * @version $Revision$
		 * @param <T>
		 */
		interface Visitor<T> {
		  /**
		   * Called by the graph traversal methods when a vertex is first visited.
		   * 
		   * @param g -
		   *          the graph
		   * @param v -
		   *          the vertex being visited.
		   */
		  public void visit(Graph<T> g, Vertex<T> v);
		}

		/*
		 * JBoss, Home of Professional Open Source Copyright 2006, Red Hat Middleware
		 * LLC, and individual contributors by the @authors tag. See the copyright.txt
		 * in the distribution for a full listing of individual contributors.
		 * 
		 * This is free software; you can redistribute it and/or modify it under the
		 * terms of the GNU Lesser General Public License as published by the Free
		 * Software Foundation; either version 2.1 of the License, or (at your option)
		 * any later version.
		 * 
		 * This software is distributed in the hope that it will be useful, but WITHOUT
		 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
		 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
		 * details.
		 * 
		 * You should have received a copy of the GNU Lesser General Public License
		 * along with this software; if not, write to the Free Software Foundation,
		 * Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA, or see the FSF
		 * site: http://www.fsf.org.
		 */

		/**
		 * A graph visitor interface that can throw an exception during a visit
		 * callback.
		 * 
		 * @author Scott.Stark@jboss.org
		 * @version $Revision$
		 * @param <T>
		 * @param <E>
		 */
		interface VisitorEX<T, E extends Exception> {
		  /**
		   * Called by the graph traversal methods when a vertex is first visited.
		   * 
		   * @param g -
		   *          the graph
		   * @param v -
		   *          the vertex being visited.
		   * @throws E
		   *           exception for any error
		   */
		  public void visit(Graph<T> g, Vertex<T> v) throws E;
		}

		/*
		 * JBoss, Home of Professional Open Source Copyright 2006, Red Hat Middleware
		 * LLC, and individual contributors by the @authors tag. See the copyright.txt
		 * in the distribution for a full listing of individual contributors.
		 * 
		 * This is free software; you can redistribute it and/or modify it under the
		 * terms of the GNU Lesser General Public License as published by the Free
		 * Software Foundation; either version 2.1 of the License, or (at your option)
		 * any later version.
		 * 
		 * This software is distributed in the hope that it will be useful, but WITHOUT
		 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
		 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
		 * details.
		 * 
		 * You should have received a copy of the GNU Lesser General Public License
		 * along with this software; if not, write to the Free Software Foundation,
		 * Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA, or see the FSF
		 * site: http://www.fsf.org.
		 */

		/**
		 * A spanning tree visitor callback interface
		 * 
		 * @see Graph#dfsSpanningTree(Vertex, DFSVisitor)
		 * 
		 * @author Scott.Stark@jboss.org
		 * @version $Revision$
		 * @param <T>
		 */
		interface DFSVisitor<T> {
		  /**
		   * Called by the graph traversal methods when a vertex is first visited.
		   * 
		   * @param g -
		   *          the graph
		   * @param v -
		   *          the vertex being visited.
		   */
		  public void visit(Graph<T> g, Vertex<T> v);

		  /**
		   * Used dfsSpanningTree to notify the visitor of each outgoing edge to an
		   * unvisited vertex.
		   * 
		   * @param g -
		   *          the graph
		   * @param v -
		   *          the vertex being visited
		   * @param e -
		   *          the outgoing edge from v
		   */
		  public void visit(Graph<T> g, Vertex<T> v, Edge<T> e);
		}

		
		
		
	@Override
	public boolean overlaps(Annotation other) {
		// FIXME Auto-generated method stub
		throw new UnsupportedOperationException("TODO");
	}

	@Override
	public <X extends Annotation> AnnotationCollection<DerivedAnnotation<X>> convertCoordinates(AnnotationCollection<X> readCollection, CoordinateSpace referenceSpaceForAnnotations, boolean fullyContained){
		//TODO Should check the coordinate space and convert appropriately
		return convertFromReferenceSpace(readCollection, referenceSpaceForAnnotations, fullyContained);
	}

	private <X extends Annotation> AnnotationCollection<DerivedAnnotation<X>> convertFromReferenceSpace(AnnotationCollection<X> readCollection, CoordinateSpace referenceSpace, boolean fullyContained){
		return new ConvertedSpace<X>(readCollection, this, referenceSpace, fullyContained);
	}

	
	
	@Override
	public CloseableIterator<? extends PopulatedWindow<T>> getPopulatedWindows(Annotation region, int windowLength){
		CloseableIterator<T> iter=sortedIterator(region, false);
		return new WindowIterator<T>(iter, windowLength, region);
	}

	@Override
	public CloseableIterator<? extends PopulatedWindow<T>> getPopulatedWindows(Annotation region, int winSize, int stepSize) {
		CloseableIterator<T> iter=sortedIterator(region, false);
		return new WindowIterator<T>(iter,winSize,region,stepSize);
	}	
	
	@Override
	public CloseableIterator<? extends PopulatedWindow<T>> getPopulatedWindows(Annotation region, int winSize, int stepSize, boolean includeEmpties) {
		CloseableIterator<T> iter=sortedIterator(region, false);
		return new WindowIterator<T>(iter,winSize,region,stepSize,includeEmpties);
	}
	
	@Override
	public int numOverlappers(Annotation region, boolean fullyContained) throws Exception {
		int counter=0;
		CloseableIterator<T> iter=sortedIterator(region, fullyContained);
		try {
			while(iter.hasNext()){
				iter.next();
				counter++;
			}
			iter.close();
			return counter;
		} catch(Exception e) {
			iter.close();
			throw(e);
		}
	}
	
	@Override
	public void writeToBAM(String fileName){
		writeToBAM(fileName, sortedIterator());
	}
	

	@Override
	public void writeToBAM(String fileName, Annotation region, boolean fullyContained){
		writeToBAM(fileName, sortedIterator(region, fullyContained));
	}
	
	private void writeToBAM(String fileName, CloseableIterator<T> iter){
		SAMFileHeader header=getReferenceCoordinateSpace().getBAMFileHeader();
		SAMFileWriter writer=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(header, false, new File(fileName));
			
		while(iter.hasNext()){
			T ann=iter.next();
			writer.addAlignment(ann.getSamRecord(header));
		}
		iter.close();
		writer.close();
		
	}
	


	/**
	 * This class requires that you have a sorted iterator of reads
	 * @author mguttman
	 *
	 * @param <T1>
	 */
	public class WindowIterator<T1 extends Annotation> implements CloseableIterator<PopulatedWindow<T1>>{

		IntervalTree<PopulatedWindow<T1>> windows;
		CloseableIterator<T1> iter;
		Iterator<PopulatedWindow<T1>> fullyFormedWindows;
		int windowLength;
		boolean hasNext;
		boolean assumeForward;
		Annotation region;
		int stepSize;
		boolean includeEmpties;
		
		PopulatedWindow<T1> nextWin;
		int nextPos;
		
		public WindowIterator(CloseableIterator<T1> iter, int windowLength, boolean assumeForward){
			this.iter=iter;
			this.windowLength=windowLength;
			this.windows=new IntervalTree<PopulatedWindow<T1>>();
			this.hasNext=false;
			this.assumeForward=assumeForward;
			this.stepSize=1;
			this.includeEmpties=false;
			
		}
		
		public WindowIterator(CloseableIterator<T1> iter, int windowLength)
		{
			this(iter,windowLength,true);
			this.region = null;
		}

		public WindowIterator(CloseableIterator<T1> iter, int windowLength,Annotation region) {
			this(iter,windowLength,true);
			this.region = region;
			this.nextPos = region.getReferenceStartPosition();
		}
		
		public WindowIterator(CloseableIterator<T1> iter, int windowLength,Annotation region, int stepSize) {
			this(iter,windowLength,region);
			this.stepSize = stepSize;
		}

		public WindowIterator(CloseableIterator<T1> iter, int windowLength,Annotation region, int stepSize, boolean includeEmpties)
		{
			this(iter,windowLength,region);
			this.stepSize = stepSize;
			this.includeEmpties = includeEmpties;
		}
		
		@Override
		public boolean hasNext() {
			if(fullyFormedWindows!=null && fullyFormedWindows.hasNext()){return true;}
			boolean hasEmpties = (includeEmpties && region!=null && nextPos < region.getReferenceEndPosition());
			//else if(iter.hasNext()){updateWindows(); return hasNext();}
			//else if(!iter.hasNext() && !windows.isEmpty())
			//{
			//	updateRemainingWindows();
			//	return hasNext();
			//}
			
			updateWindows();
			return (hasEmpties || fullyFormedWindows!=null && fullyFormedWindows.hasNext());
			
			//return false;
		}

		private void updateRemainingWindows(){
			fullyFormedWindows = windows.valueIterator();
			windows = new IntervalTree<PopulatedWindow<T1>>();
		}
		
		@Override
		public PopulatedWindow<T1> next() {
			PopulatedWindow<T1> rtrn;
			
			if(nextWin == null)
			{
				if(fullyFormedWindows!=null && fullyFormedWindows.hasNext())
					nextWin = fullyFormedWindows.next();
				else
				{
					nextPos++;
					return new ContiguousWindow<T1>(region.getName(), nextPos, nextPos+windowLength, Strand.BOTH);

				}
			}
			
			if(!includeEmpties || nextWin.getReferenceStartPosition()==nextPos)
			{
				rtrn = nextWin;
				nextWin = null;
			}
			else
				rtrn = new ContiguousWindow<T1>(region.getName(), nextPos, nextPos+windowLength, Strand.BOTH);
			
			nextPos++;
			return rtrn;
		}

		private void updateWindows(){
			while((fullyFormedWindows==null || !fullyFormedWindows.hasNext()) && iter.hasNext())
			{	T1 read=iter.next();
				//all windows with an end position before the start of this window
				fullyFormedWindows=removeFullyFormedWindows(read).iterator();
				addReadToWindows(read);
			}
			if((fullyFormedWindows==null || !fullyFormedWindows.hasNext()) && !windows.isEmpty())
			{
				updateRemainingWindows();
			}
		}

		private void addReadToWindows(T1 read){
			//Create all windows which overlap read blocks
			Iterator<SingleInterval> interval=read.getBlocks();
			while(interval.hasNext()){
				SingleInterval block=interval.next();
				int start = 0;
				if(region!=null)
					start=Math.max(roundUp(region.getReferenceStartPosition()-windowLength+1), roundUp(block.getReferenceStartPosition()-windowLength));
				else
					start = Math.max(0,roundUp(block.getReferenceStartPosition()-windowLength));
				int end = block.getReferenceEndPosition();
				if(region!=null)
					end = Math.min(region.getReferenceEndPosition(), block.getReferenceEndPosition());
				for(int i=start; i<end; i+=stepSize){
					PopulatedWindow<T1> window=windows.remove(i, i+windowLength);
					if(window==null){
						//make a window
						window=new ContiguousWindow<T1>(read.getReferenceName(), i, i+windowLength, Strand.BOTH);
					}
					window.addAnnotation(read);
					windows.put(window.getReferenceStartPosition(), window.getReferenceEndPosition(), window);
				}
			}
		}
		
		private int roundUp(int winStart)
		{
			return (winStart + stepSize -1) / stepSize *stepSize;
		}
		
		private Collection<PopulatedWindow<T1>> removeFullyFormedWindows(T1 read) {
			Iterator<PopulatedWindow<T1>> iter;
			if(assumeForward)
				{iter=windows.getNodesBeforeInterval(read.getReferenceStartPosition(), read.getReferenceStartPosition());}
			else
				{iter=windows.getNodesAfterInterval(read.getReferenceEndPosition(), read.getReferenceEndPosition());}
			ArrayList<PopulatedWindow<T1>> rtrn=new ArrayList<PopulatedWindow<T1>>();
			while(iter.hasNext()){
				PopulatedWindow<T1> w=iter.next();
				if((!w.getReferenceName().equalsIgnoreCase(read.getReferenceName())) || (w.getReferenceEndPosition()<read.getReferenceStartPosition())){
					rtrn.add(w);
					windows.remove(w.getReferenceStartPosition(), w.getReferenceEndPosition());
				}
			}
			this.hasNext=!rtrn.isEmpty();
			Comparator comparator = new WindowComparator();
			Collections.sort(rtrn,comparator);
			return rtrn;
		}

		
		@Override
		public void remove() {
			this.iter.remove();
		}

		@Override
		public void close() {
			this.iter.close();
		}


	}

	public class WindowComparator implements Comparator<PopulatedWindow<T>>{
		@Override
		public int compare(PopulatedWindow<T> win1, PopulatedWindow<T> win2)
		{
			return win1.getReferenceStartPosition()-win2.getReferenceEndPosition();
		}
	}
	
	@Override
	public CoordinateSpace getFeatureCoordinateSpace(){
		//Iterate through all records
		CloseableIterator<T> iter=sortedIterator();
		Map<String, Integer> sizes=new TreeMap<String, Integer>();
		while(iter.hasNext()){
			T annotation=iter.next();
			int size=annotation.size();
			String name=annotation.getName();
			sizes.put(name, size);
		}
		iter.close();
		return new CoordinateSpace(sizes);
	}
	
	@Override
	public int getNumAnnotations(){
		if(this.numAnnotations != 0)
			return this.numAnnotations;

		CloseableIterator<T> iter = this.sortedIterator();
		int count = 0;
		while(iter.hasNext())
		{
			iter.next();
			count++;
		}
		iter.close();
		this.numAnnotations = count;
		return this.numAnnotations;
	}
	
	
}
