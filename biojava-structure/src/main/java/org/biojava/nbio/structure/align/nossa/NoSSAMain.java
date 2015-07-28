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
import org.biojava.nbio.structure.AtomImpl;
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
	private Map<Matrix4d, List<Integer[]>> map;

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

		afp1 = calculateFragments(ca1);
		afp2 = calculateFragments(ca2);

		calculateTransformsMatrix();
		generateTransformMap();

		//Select the representative transformation 
		int maxNumber = 0;
		Matrix4d repr = null;
		for (Matrix4d v : map.keySet()){
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

		AFPChain afpChain = new AFPChain();
		AlignmentTools.replaceOptAln(optAln, afpChain, ca1, ca2);

		return afpChain;
	}

	@Override
	public AFPChain align(Atom[] ca1, Atom[] ca2, Object params)
			throws StructureException {

		setParameters((NoSSAParameters) params); 
		return align(ca1, ca2);
	}

	private List<Atom[]> calculateFragments(Atom[] ca){

		List<Atom[]> afps = new ArrayList<Atom[]>();
		int afpSize = params.getAfpSize();
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

		map = new HashMap<Matrix4d, List<Integer[]>>();
		double epsilon = params.getEpsilon();

		for (int i=0; i<afp1.size(); i++){
			for (int j=0; j<afp2.size(); j++){

				Integer[] pair = {i,j};
				Matrix4d t = matrix.get(i).get(j);
				if (t == null) continue;

				double minDiff = Double.MAX_VALUE;
				Matrix4d key = null;

				for (Matrix4d r : map.keySet()){
					double diff = compareTransforms(t, r);
					if (diff<minDiff){
						minDiff = diff;
						key = r;
					}
				}

				if (minDiff < epsilon){
					map.get(key).add(pair);
				}
				else {
					map.put(t, new ArrayList<Integer[]>());
					map.get(t).add(pair);
				}
			}
		}
	}

	private double compareTransforms(Matrix4d a, Matrix4d b){

		//Normalize the translation vector direction
		Vector3d v1 = new Vector3d(a.m03, a.m13, a.m23);
		Vector3d v2 = new Vector3d(b.m03, b.m13, b.m23);
		v1.normalize();
		v2.normalize();

		/*Matrix4d trans1 = (Matrix4d) a.clone();
		Matrix4d trans2 = (Matrix4d) b.clone();

		trans1.setColumn(3, new double[] {v1.x, v1.y, v1.z, 1.0});
		trans2.setColumn(3, new double[] {v2.x, v2.y, v2.z, 1.0});

		Atom at = new AtomImpl();
		at.setCoords(new double[] {0, 0, 0});

		Atom rotA = (Atom) at.clone();
		Calc.transform(rotA, trans1);
		Atom rotB = (Atom) at.clone();
		Calc.transform(rotB, trans2);

		//Compare transformations by Atom transformed difference
		return Calc.getDistanceFast(rotA, rotB);*/
		
		return 1.0-v1.dot(v2);
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
