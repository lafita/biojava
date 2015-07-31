package org.biojava.nbio.structure.align.nossa;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AlignmentTools;

/**
 * NoSSA stands for: Non-Sequential Structure Alignment.
 * <p>
 * This structure alignment algorithm uses geometrical properties and some 
 * heuristics to align two protein structures.
 * <p>
 * Its time complexity is O(n*m), where n and m are the length of the two 
 * structures compared.
 * <p>
 * The algorithm divides the structures into fragments and calculates the 
 * transformations from each AFP of the first structure to all other AFPs 
 * of the second structure, as an n*m Matrix. After that, it determines 
 * which groups of AFP pairs share the same rotation and translation direction 
 * by clustering them given a threshold. That allows to determine the equivalent 
 * positions (EQRs) between the two structures, which do not need to be 
 * sequential.
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
	public static final boolean debug = true;

	private NoSSAParameters params;
	
	private List<List<Matrix4d>> matrix;
	private Map<Matrix4d, List<Integer[]>> tmap;
	
	private List<List<Vector3d>> vectors1;
	private List<List<Vector3d>> vectors2;
	private Map<Vector3d, List<Integer[][]>> vmap;
	
	private Atom[] ca1;
	private Atom[] ca2;
	private List<Atom[]> afp1;
	private List<Atom[]> afp2;

	public NoSSAMain(){
		params = new NoSSAParameters();
	}

	@Override
	public AFPChain align(Atom[] ca10, Atom[] ca20) throws StructureException {

		ca1 = ca10;
		ca2 = ca20;

		afp1 = calculateFragments(ca1, true);
		afp2 = calculateFragments(ca2, true);

		calculateTransformsMatrix();
		generateTransformMap();

		//Select the representative transformation 
		int maxNumber = 0;
		Matrix4d repr = null;
		for (Matrix4d v : tmap.keySet()){
			List<Integer[]> pairs = tmap.get(v);
			if (pairs.size() > maxNumber){
				repr = v;
				maxNumber = pairs.size();
			}
		}

		List<Integer[]> pairs = tmap.get(repr);
		
		chainFragments(pairs);

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

		AFPChain afpChain = new AFPChain();
		AlignmentTools.replaceOptAln(optAln, afpChain, ca1, ca2);
		
		/*vectors1 = getVectors(ca1);
		vectors2 = getVectors(ca2);*/

		return null;
	}

	@Override
	public AFPChain align(Atom[] ca1, Atom[] ca2, Object params)
			throws StructureException {

		setParameters((NoSSAParameters) params); 
		return align(ca1, ca2);
	}

	private List<Atom[]> calculateFragments(Atom[] ca, boolean complete){

		List<Atom[]> afps = new ArrayList<Atom[]>();
		int afpSize = params.getAfpSize();
		
		if (complete){
			for (int i=0; i<ca.length-afpSize; i++){
				afps.add(new Atom[afpSize]);
			}
	
			for (int i=0; i<ca.length; i++){
				for (int f=0; f<afpSize; f++){
					if (i-f < 0) break;
					if (i-f > afps.size()-1) continue;
					afps.get(i-f)[f] = ca[i];
				}
			}
		} else {
			for (int i=0; i<ca.length/afpSize; i++){
				Atom[] afp = new Atom[afpSize];
				for (int f=0; f<afpSize; f++){
					afp[f] = ca[i*afpSize+f];
				}
				afps.add(afp);
			}
		}
		return afps;
	}

	private void calculateTransformsMatrix() throws StructureException{

		matrix = new ArrayList<List<Matrix4d>>();
		int nulls = 0;

		for (int i=0; i<afp1.size(); i++){
			matrix.add(new ArrayList<Matrix4d>());
			for (int j=0; j<afp2.size(); j++){
				Atom[] clone2 = StructureTools.cloneAtomArray(afp2.get(j));
				SVDSuperimposer svd = 
						new SVDSuperimposer(afp1.get(i), clone2);
				Matrix4d t = svd.getTransformation();
				Calc.transform(clone2, t);
				double rmsd = SVDSuperimposer.getRMS(afp1.get(i), clone2);
				if (rmsd > params.getMaxAfpRMSD()){
					matrix.get(i).add(null);
					nulls++;
				} else {
					matrix.get(i).add(t);
				}
			}
		}
		if (debug) System.out.println("Number of nulls: "+nulls+
				" out of "+afp1.size()*afp2.size()+" AFP comparisons.");
	}

	private void generateTransformMap(){

		tmap = new HashMap<Matrix4d, List<Integer[]>>();
		double epsilon = params.getEpsilon();

		for (int i=0; i<afp1.size(); i++){
			for (int j=0; j<afp2.size(); j++){

				Integer[] pair = {i,j};
				Matrix4d t = matrix.get(i).get(j);
				if (t == null) continue;

				boolean added = false;

				for (Matrix4d r : tmap.keySet()){
					double minDiff = Double.MAX_VALUE;
					for (Integer[] p : tmap.get(r)){
						Matrix4d pt = matrix.get(p[0]).get(p[1]);
						double diff = compareTransforms(pt, t);
						if (diff < minDiff) minDiff = diff;
					}
					if (minDiff < epsilon){
						tmap.get(r).add(pair);
						added = true;
					}
				}

				if (!added) {
					tmap.put(t, new ArrayList<Integer[]>());
					tmap.get(t).add(pair);
				}
			}
		}
	}
	
	private void chainFragments(List<Integer[]> pairs){
		
		Collections.sort(pairs, new Comparator<Integer[]>() {
			@Override
			public int compare(Integer[] o1, Integer[] o2) {
				return o1[0]-o2[0];
			}
		});
	}
	
	private void generateVectorMap(){
		
		//Will not work...
	}
	
	private static double compareVectors(Vector3d a, Vector3d b){
		
		Vector3d v = (Vector3d) a.clone();
		v.sub(b);
		return v.length();
	}

	private static double compareTransforms(Matrix4d a, Matrix4d b){

		//Error: translation, rotation axis, angle
		double error = 0;
		
		//Normalize the translation vector direction
		Vector3d v1 = new Vector3d(a.m03, a.m13, a.m23);
		Vector3d v2 = new Vector3d(b.m03, b.m13, b.m23);
		v1.normalize();
		v2.normalize();
		
		error += 1.0-v1.dot(v2);

		//Obtain the normalized axis vector
		AxisAngle4d axis1 = new AxisAngle4d();
		axis1.set(a);
		AxisAngle4d axis2 = new AxisAngle4d();
		axis2.set(b);
		
		Vector3d r1 = new Vector3d(axis1.x, axis1.y, axis1.z);
		Vector3d r2 = new Vector3d(axis2.x, axis2.y, axis2.z);
		r1.normalize();
		r2.normalize();
		
		error += 1.0-(Math.abs(r1.dot(r2)));
		
		error += Math.abs(axis1.angle - axis2.angle);

		return error;
	}
	
	private static List<List<Vector3d>> getVectors(Atom[] ca){
		
		List<List<Vector3d>> vectors = new ArrayList<List<Vector3d>>();
		
		for (int i=0; i<ca.length; i++){
			vectors.add(new ArrayList<Vector3d>());
			
			for (int j=0; j<ca.length; j++){
				Vector3d a = new Vector3d(ca[i].getCoords());
				Vector3d b = new Vector3d(ca[j].getCoords());
				a.sub(b);
				vectors.get(i).add(a);
			}
		}
		
		return vectors;
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
