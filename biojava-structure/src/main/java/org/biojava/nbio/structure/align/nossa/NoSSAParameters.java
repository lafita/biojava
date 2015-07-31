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
 * Contains the parameters that can be sent to NoSSA.
 * 
 * @author Aleix Lafita
 *
 */
public class NoSSAParameters implements ConfigStrucAligParams {

	private double epsilon;
	private int afpSize;
	private double maxAfpRMSD;

	public NoSSAParameters(){
		reset();
	}

	@Override
	public String toString() {
		return "NoSSAParameters [epsilon=" + epsilon + "]";
	}

	@Override
	public void reset(){
		epsilon = 0.07;
		afpSize = 8;
		maxAfpRMSD = 2.0;
	}

	public double getEpsilon(){
		return epsilon;
	}
	
	public void setEpsilon(Double epsilon) {
		this.epsilon = epsilon;
	}

	public int getAfpSize() {
		return afpSize;
	}

	public void setAfpSize(Integer afpSize) {
		this.afpSize = afpSize;
	}

	public double getMaxAfpRMSD() {
		return maxAfpRMSD;
	}

	public void setMaxAfpRMSD(Double maxAfpRMSD) {
		this.maxAfpRMSD = maxAfpRMSD;
	}

	@Override
	public List<String> getUserConfigHelp() {
		List<String> params =new ArrayList<String>();
		String epsilon = "Maximum allowed transformation difference.";
		String afpSize = "Size of the structure fragments aligned.";
		String maxAfpRMSD = "Maximum allowed RMSD between the AFPs.";
		params.add(epsilon);
		params.add(afpSize);
		params.add(maxAfpRMSD);
		return params;
	}

	@Override
	public List<String> getUserConfigParameters() {
		List<String> params = new ArrayList<String>();
		params.add("Epsilon");
		params.add("AfpSize");
		params.add("MaxAfpRMSD");
		return params;
	}

	@Override
	public List<String> getUserConfigParameterNames(){
		List<String> params = new ArrayList<String>();
		params.add("Epsilon");
		params.add("AFP Size");
		params.add("Max AFP RMSD");
		return params;
	}

	@Override
	@SuppressWarnings("rawtypes")
	public List<Class> getUserConfigTypes() {
		List<Class> params = new ArrayList<Class>();
		params.add(Double.class);
		params.add(Double.class);
		params.add(Integer.class);
		return params;
	}

}
