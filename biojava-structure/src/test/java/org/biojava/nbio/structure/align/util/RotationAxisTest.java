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
 * Created on Oct 7, 2013
 * Author: blivens
 *
 */

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
 * Created on Oct 7, 2013
 * Author: blivens
 *
 */
package org.biojava.nbio.structure.align.util;

import org.junit.Test;

import static org.junit.Assert.*;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

/**
 * Test the Rotation Axis.
 * 
 * @author Spencer Bliven
 */
public class RotationAxisTest {

	@Test
	public void testProjection() throws Exception {
		
		RotationAxis axis;
		
		// 180 around z
		Vector3d dir = new Vector3d(0,0,1);
		Point3d pos = new Point3d(0,0,0);
		
		axis = new RotationAxis(dir, pos, Math.PI);

		pos.set(1,2,3);
		Point3d projected = axis.getProjectedPoint(pos);
		assertTrue(new Point3d(0,0,3).epsilonEquals(projected, 1e-14));

		double dist = axis.getProjectedDistance(pos);
		assertEquals(Math.sqrt(5),dist,1e-14);


		// main diagonal through (1,1,0)
		dir.set(2,2,2);
		pos.set(1,1,0);
		axis = new RotationAxis(dir, pos, Math.PI);

		pos.set(1,1,0);
		projected = axis.getProjectedPoint(pos);
		assertTrue(new Point3d(1,1,0).epsilonEquals(projected, 1e-14));

		pos.set(0,0,-1);
		projected = axis.getProjectedPoint(pos);
		assertTrue(new Point3d(0,0,-1).epsilonEquals(projected, 1e-14));

		pos.set(-.5,-.5,0);
		projected = axis.getProjectedPoint(pos);
		assertTrue(new Point3d(0,0,-1).epsilonEquals(projected, 1e-14));

		dist = axis.getProjectedDistance(pos);
		assertEquals(Math.sqrt(3/2.),dist,1e-14);

		pos.set(0,0,0);
		projected = axis.getProjectedPoint(pos);
		assertTrue(new Point3d(1/3.,1/3.,-2/3.).epsilonEquals(projected, 1e-14));

	}

}
