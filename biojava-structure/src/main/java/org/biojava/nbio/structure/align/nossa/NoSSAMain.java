package org.biojava.nbio.structure.align.nossa;

import java.util.List;
import javax.vecmath.GMatrix;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.secstruc.SecStrucElement;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * NoSSA stands for: Non-Sequential Structure Alignment.
 * <p>
 * This structure alignment algorithm uses the secondary structure of
 * proteins, together with other geometrical properties, to align two 
 * protein structures.
 * <p>
 * Its time complexity is O(n*m), where n and m are the length of the two 
 * structures compared. However, because the structures are reduced into
 * coarse secondary structure elements (SSE), it is a factor faster than 
 * normal sequential CA-CA alignment algorithms.
 * <p>
 * The algorithm divides the structures into its SSE and calculates the 
 * transformations from each SSE of the first structure to all other 
 * equivalent SSE of the second structure, in a x*y Matrix. 
 * <p>
 * After that, it determines which groups of SSE pairs share the same 
 * transformation, so that they would be consistent in a global alignment.
 * This allows a non-sequential superposition of two structures, which can
 * be refined to obtain the equivalent positions (EQRs) with a Monte Carlo
 * sampling algorithm.
 * 
 * @author Aleix Lafita
 * @since 4.1.1
 *
 */
public class NoSSAMain implements StructureAlignment {

	/**
	 *  version history:
	 *  1.0 - Algorithm definition and implementation
	 */
	public static final String version = "1.0";
	public static final String algorithmName = "jNoSSA";
	
	private static final Logger logger = 
			LoggerFactory.getLogger(NoSSAMain.class);

	private NoSSAParameters params;
	
	private Atom[] ca1;
	private Atom[] ca2;
	
	private List<SecStrucElement> sse1;
	private List<SecStrucElement> sse2;
	
	private GMatrix simMatrix;

	public NoSSAMain(){
		params = new NoSSAParameters();
	}

	@Override
	public AFPChain align(Atom[] ca10, Atom[] ca20) throws StructureException {

		ca1 = ca10;
		ca2 = ca20;
		
		return null;
	}

	@Override
	public AFPChain align(Atom[] ca1, Atom[] ca2, Object params)
			throws StructureException {

		setParameters((ConfigStrucAligParams) params); 
		return align(ca1, ca2);
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
