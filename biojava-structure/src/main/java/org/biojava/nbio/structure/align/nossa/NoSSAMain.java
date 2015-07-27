package org.biojava.nbio.structure.align.nossa;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AlignmentTools;

/**
 * NoSSA stands for: Non-Sequential Structure Alignment.
 * <p>
 * It is a first implementation of a structure alignment algorithm that
 * uses geometrical hashing and some heuristics with a greedy approach
 * to align two protein structures.
 * <p>
 * It is a fast algorithm compared to the DP algorithms, with time complexity 
 * O(n*m), where n and m are the length of the two protein structures compared.
 * <p>
 * As a brief explanation, the algorithm calculates the translation vectors from
 * each Atom of the first structure to all other Atoms of the second structure,
 * as an n*m Matrix. After that, it determines which groups of Atom pairs share 
 * the same translation direction by clustering them given a threshold. That
 * allows to determine the equivalent positions (EQRs) between the two structures.
 * 
 * @author Aleix Lafita
 * @since 4.1.1
 *
 */
public class NoSSAMain implements StructureAlignment {

	/**
	 *  version history:
	 *  1.0 - Initial algorithm implementation
	 */
	public static final String version = "1.0";
	public static final String algorithmName = "jNoSSA";

	private NoSSAParameters params;
	private List<List<Vector3d>> matrix;
	private Map<Vector3d, List<Integer[]>> map;

	private Atom[] ca1;
	private Atom[] ca2;
	private AFPChain afpChain;

	public NoSSAMain(){
		params = new NoSSAParameters();
	}

	@Override
	public AFPChain align(Atom[] ca10, Atom[] ca20) throws StructureException {

		ca1 = ca10;
		ca2 = ca20;

		calculateVectorMatrix();
		generateVectorMap();

		//Select the representative vector 
		int maxNumber = 0;
		Vector3d repr = null;
		for (Vector3d v : map.keySet()){
			List<Integer[]> pairs = map.get(v);
			if (pairs.size() > maxNumber){
				repr = v;
				maxNumber = pairs.size();
			}
		}

		List<Integer[]> pairs = map.get(repr);
		
		Collections.sort(pairs, new Comparator<Integer[]>() {
			@Override
			public int compare(Integer[] o1, Integer[] o2) {
				return o1[0]-o2[0];
			}
		});
		
		List<List<Integer>> chain1 = new ArrayList<List<Integer>>();
		List<List<Integer>> chain2 = new ArrayList<List<Integer>>();
		List<Integer> block1 = new ArrayList<Integer>();
		List<Integer> block2 = new ArrayList<Integer>();
		int lastRes1 = 0;
		int lastRes2 = 0;
		for (int res=0; res<pairs.size(); res++){
			if (pairs.get(res)[0] == lastRes1) continue;
			else if (pairs.get(res)[1] < lastRes2){
				chain1.add(block1);
				chain2.add(block2);
				block1 = new ArrayList<Integer>();
				block2 = new ArrayList<Integer>();
				lastRes2 = pairs.get(res)[1];
			}
			block1.add(pairs.get(res)[0]);
			block2.add(pairs.get(res)[1]);
			lastRes2 = pairs.get(res)[1];
			lastRes1 = pairs.get(res)[0];
		}
		
		int[][][] optAln = new int[chain1.size()][][];
		for (int bk=0; bk<chain1.size(); bk++){
			optAln[bk] = new int[2][];
			optAln[bk][0] = new int[chain1.get(bk).size()];
			optAln[bk][1] = new int[chain1.get(bk).size()];
			for (int res=0; res<chain1.get(bk).size(); res++){
				optAln[bk][0][res] = chain1.get(bk).get(res);
				optAln[bk][1][res] = chain2.get(bk).get(res);
			}
		}
		
		afpChain = new AFPChain();
		AlignmentTools.replaceOptAln(optAln, afpChain, ca1, ca2);

		return afpChain;
	}

	@Override
	public AFPChain align(Atom[] ca1, Atom[] ca2, Object params)
			throws StructureException {

		setParameters((NoSSAParameters) params); 
		return align(ca1, ca2);
	}

	private void calculateVectorMatrix(){

		matrix = new ArrayList<List<Vector3d>>();
		for (int i=0; i<ca1.length; i++){
			matrix.add(new ArrayList<Vector3d>());
			for (int j=0; j<ca2.length; j++){
				Point3d a = new Point3d(ca1[i].getCoords());
				Point3d b = new Point3d(ca2[j].getCoords());
				Vector3d v = new Vector3d();
				v.sub(a, b);
				matrix.get(i).add(v);
			}
		}
	}

	private void generateVectorMap(){

		map = new HashMap<Vector3d, List<Integer[]>>();
		double epsilon = params.getEpsilon();

		for (int i=0; i<ca1.length; i++){
			for (int j=0; j<ca2.length; j++){

				Integer[] pair = {i,j};
				Vector3d w = (Vector3d) matrix.get(i).get(j).clone();
				w.normalize();
				
				double minDot = 0.0;
				Vector3d key = null;

				for (Vector3d v : map.keySet()){
					double dot = Math.abs(w.dot(v));
					if (dot>minDot){
						minDot = dot;
						key = v;
					}
				}
				
				if ((1-minDot) < epsilon){
					map.get(key).add(pair);
				}
				else {
					map.put(w, new ArrayList<Integer[]>());
					map.get(w).add(pair);
				}
			}
		}
	}

	@Override
	public ConfigStrucAligParams getParameters() {
		return params;
	}

	@Override
	public void setParameters(ConfigStrucAligParams parameters) {
		if (! (params instanceof NoSSAParameters )){
			throw new IllegalArgumentException(
					"provided parameter object is not of the correct type");
		}
		this.params = (NoSSAParameters) params;
	}

	@Override
	public String getAlgorithmName() {
		return algorithmName;
	}

	@Override
	public String getVersion() {
		return version;
	}

}
