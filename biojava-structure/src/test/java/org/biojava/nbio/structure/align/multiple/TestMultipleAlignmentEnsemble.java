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
 */
package org.biojava.nbio.structure.align.multiple;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentScorer;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Test methods in MultipleAlignmentEnsemle.
 *
 * @author Aleix Lafita
 * @since 4.2.0
 *
 */
public class TestMultipleAlignmentEnsemble {

	/**
	 * Test that all relevant information (except scores and cache) is correctly
	 * copied from the AFPChain to the generated MultipleAlignmentEnsemble
	 * object.
	 */
	@Test
	public static void testAFPconversion() {

		// Fill an AFPChain with the general information
		AFPChain afp = new AFPChain("algorithm");
		afp.setName1("name1");
		afp.setName2("name2");
		afp.setVersion("1.0");
		afp.setCalculationTime(System.currentTimeMillis());
		// Generate a optimal alignment with 3 blocks and 5 residues per block
		int[][][] optAln = new int[3][][];
		for (int b = 0; b < optAln.length; b++) {
			int[][] block = new int[2][];
			for (int c = 0; c < block.length; c++) {
				int[] residues = { b + 5, b + 6, b + 7, b + 8, b + 9 };
				block[c] = residues;
			}
			optAln[b] = block;
		}
		afp.setOptAln(optAln);
		afp.setBlockNum(optAln.length);
		// Set the transformation matrix to random numbers
		Matrix4d transform = new Matrix4d();
		transform.setRow(0, new double[] { 0.13, 1.5, 0.84, 0.44 });
		transform.setRow(1, new double[] { 1.3, 0.44, 2.3, 0.21 });
		transform.setRow(2, new double[] { 1.0, 1.2, 2.03, 0.89 });
		Matrix4d[] blockTransform = { transform, transform, transform };
		afp.setBlockTransformation(blockTransform);

		// Convert the AFPChain into a MultipleAlignment (without Atoms)
		MultipleAlignmentEnsemble ensemble = new MultipleAlignmentEnsembleImpl(
				afp, null, null, true);
		MultipleAlignment msa = ensemble.getMultipleAlignment(0);

		// Test for all the information
		assertEquals(afp.getName1(), ensemble.getStructureIdentifiers().get(0)
				.getIdentifier());
		assertEquals(afp.getName2(), ensemble.getStructureIdentifiers().get(1)
				.getIdentifier());
		assertEquals(afp.getAlgorithmName(), ensemble.getAlgorithmName());
		assertEquals(afp.getVersion(), ensemble.getVersion());
		assertTrue(ensemble.getCalculationTime().equals(
				afp.getCalculationTime()));
		assertEquals(afp.getBlockNum(), msa.getBlockSets().size());
		for (int b = 0; b < afp.getBlockNum(); b++) {
			assertEquals(afp.getBlockTransformation()[b], msa.getBlockSet(b)
					.getTransformations().get(1));
		}

		// Test for the scores
		assertEquals(msa.getScore(MultipleAlignmentScorer.CE_SCORE),
				(Double) afp.getAlignScore());
		assertEquals(msa.getScore(MultipleAlignmentScorer.AVGTM_SCORE),
				(Double) afp.getTMScore());
		assertEquals(msa.getScore(MultipleAlignmentScorer.RMSD),
				(Double) afp.getTotalRmsdOpt());

		// Test for the optimal alignment
		for (int b = 0; b < 3; b++) {
			for (int c = 0; c < 2; c++) {
				for (int res = 0; res < 5; res++) {
					Integer afpRes = afp.getOptAln()[b][c][res];
					assertEquals(afpRes, msa.getBlock(b).getAlignRes().get(c)
							.get(res));
				}
			}
		}
	}
}
