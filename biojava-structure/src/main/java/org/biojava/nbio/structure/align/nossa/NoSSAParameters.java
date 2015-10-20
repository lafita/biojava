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
 * Created on Sep 15, 2009
 * Author: Andreas Prlic 
 *
 */
package org.biojava.nbio.structure.align.nossa;

import org.biojava.nbio.structure.align.ce.ConfigStrucAligParams;
import java.util.ArrayList;
import java.util.List;

/** 
 * Contains the parameters that can be sent to the NoSSA
 * algorithm.
 * 
 * @author Aleix Lafita
 *
 */
public class NoSSAParameters implements ConfigStrucAligParams {

	private int minSSEsize;

	public NoSSAParameters(){
		reset();
	}

	@Override
	public void reset(){
		minSSEsize = 3;
	}

	public int getMinSSEsize() {
		return minSSEsize;
	}

	public void setMinSSEsize(int minSSEsize) {
		this.minSSEsize = minSSEsize;
	}

	@Override
	public List<String> getUserConfigHelp() {
		List<String> params =new ArrayList<String>();
		String minSSEsize = "Minimum size of the SS elements.";
		params.add(minSSEsize);
		return params;
	}

	@Override
	public List<String> getUserConfigParameters() {
		List<String> params = new ArrayList<String>();
		params.add("MinSSEsize");
		return params;
	}

	@Override
	public List<String> getUserConfigParameterNames(){
		List<String> params = new ArrayList<String>();
		params.add("Min SSE Size");
		return params;
	}

	@Override
	@SuppressWarnings("rawtypes")
	public List<Class> getUserConfigTypes() {
		List<Class> params = new ArrayList<Class>();
		params.add(Integer.class);
		return params;
	}

}
