/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on Jan 21, 2010
 *
 */
package demo;

import java.io.IOException;

import javax.vecmath.Matrix4d;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.nossa.NoSSAMain;
import org.biojava.nbio.structure.align.util.AtomCache;

/** 
 * Example of how to run a structure alignment using the NoSSA algorithm.
 * 
 * @author Aleix Lafita
 *
 */
public class DemoNoSSA {

	public static void main(String[] args) 
			throws IOException, StructureException {

		String name1 = "4hhb.A";
		String name2 = "4hhb.B";

		AtomCache cache = new AtomCache();

		StructureAlignment algorithm  = new NoSSAMain();

		Structure structure1 = cache.getStructure(name1);
		Structure structure2 = cache.getStructure(name2);

		Atom[] ca1 = StructureTools.getAtomCAArray(structure1);
		Atom[] ca2 = StructureTools.getAtomCAArray(structure2);
		
		/*Matrix4d t = new Matrix4d(new double[] {
			1, 0, 0,1,
			0, 0,-1,1,
			0, 1, 0,1,
			0, 0, 0,1});
		
		Calc.transform(structure2, t);*/

		AFPChain afpChain = algorithm.align(ca1,ca2);

		afpChain.setName1(name1);
		afpChain.setName2(name2);

		StructureAlignmentDisplay.display(afpChain, ca1, ca2);

	}
}
