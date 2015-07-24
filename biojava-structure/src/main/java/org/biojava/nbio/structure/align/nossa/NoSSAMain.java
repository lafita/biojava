package org.biojava.nbio.structure.align.nossa;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.nbio.structure.align.model.AFPChain;

/**
 * NoSSA stands for: Non-Sequential Structure Alignment.
 * <p>
 * It is a first implementation of a structure alignment algorithm that
 * uses geometrical hashing and some heuristics with a greedy approach
 * to align two protein structures.
 * <p>
 * It is a fast algorithm compared to the DP algorithms, with time complexity 
 * O(n*m), where n and m are the length of the two protein structures compared.
 * One of its advantages is that the core calculations can be distributed, 
 * using multiple threads or cores, thus reducing the constant factors of the
 * running time.
 * <p>
 * 
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

	@Override
	public AFPChain align(Atom[] ca1, Atom[] ca2) throws StructureException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public AFPChain align(Atom[] ca1, Atom[] ca2, Object params)
			throws StructureException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public ConfigStrucAligParams getParameters() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void setParameters(ConfigStrucAligParams parameters) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public String getAlgorithmName() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getVersion() {
		// TODO Auto-generated method stub
		return null;
	}

}
