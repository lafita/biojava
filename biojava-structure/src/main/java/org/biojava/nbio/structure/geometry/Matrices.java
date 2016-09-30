package org.biojava.nbio.structure.geometry;

import javax.vecmath.GMatrix;
import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;

import org.biojava.nbio.structure.jama.Matrix;

/**
 * Matrices contains static methods to operate and transform matrices used in 3D
 * geometry (transformation matrices, rotation matrices, distance matrices,
 * etc).
 * <p>
 * This class complements and extends the functionallity of vecmath and JAMA, as
 * well as providing methods to convert between them.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public class Matrices {

	/** Prevent instantiation */
	private Matrices() {
	}

	/**
	 * Convert a transformation matrix into a JAMA rotation matrix. Because the
	 * JAMA matrix is a pre-multiplication matrix and the Vecmath matrix is a
	 * post-multiplication one, the rotation matrix is transposed to ensure that
	 * the transformation they produce is the same.
	 *
	 * @param transform
	 *            Matrix4d with transposed rotation matrix
	 * @return rotation matrix as JAMA object
	 */
	@Deprecated
	public static Matrix getRotationJAMA(Matrix4d transform) {

		Matrix rot = new Matrix(3, 3);
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				rot.set(j, i, transform.getElement(i, j)); // transposed
			}
		}
		return rot;
	}

	/**
	 * Convert a transformation matrix into a rotation matrix.
	 *
	 * @param transform
	 *            Matrix4d
	 * @return rotation matrix
	 */
	public static Matrix3d getRotationMatrix(Matrix4d transform) {

		Matrix3d rot = new Matrix3d();
		transform.setRotationScale(rot);
		return rot;
	}

	/**
	 * Extract the translational vector of a transformation matrix.
	 *
	 * @param transform
	 *            Matrix4d
	 * @return Vector3d translation vector
	 */
	public static Vector3d getTranslationVector(Matrix4d transform) {
		Vector3d transl = new Vector3d();
		transform.get(transl);
		return transl;
	}

	/**
	 * Convert JAMA rotation and translation to a Vecmath transformation matrix.
	 * Because the JAMA matrix is a pre-multiplication matrix and the Vecmath
	 * matrix is a post-multiplication one, the rotation matrix is transposed to
	 * ensure that the transformation they produce is the same.
	 *
	 * @param rot
	 *            3x3 Rotation matrix
	 * @param trans
	 *            3x1 Translation matrix
	 * @return 4x4 transformation matrix
	 */
	@Deprecated
	public static Matrix4d getTransformation(Matrix rot, Matrix trans) {
		return new Matrix4d(new Matrix3d(rot.getColumnPackedCopy()),
				new Vector3d(trans.getColumnPackedCopy()), 1.0);
	}

	/**
	 * Convert a JAMA Matrix to a vecmath GMatrix.
	 * 
	 * @param jama
	 *            JAMA Matrix
	 * @return copy of the matrix as GMatrix
	 */
	public static GMatrix jamaToVecmath(Matrix jama) {
		return doubleToVecmath(jama.getArray());
	}

	/**
	 * Convert a double[][] matrix to a vecmath GMatrix.
	 * 
	 * @param matrix
	 *            double[][] matrix
	 * @return copy of the matrix as GMatrix
	 */
	public static GMatrix doubleToVecmath(double[][] matrix) {

		// Empty matrix case
		if (matrix.length == 0) {
			return new GMatrix(0, 0);
		}

		GMatrix mat = new GMatrix(matrix.length, matrix[0].length);
		for (int r = 0; r < matrix.length; r++)
			mat.setRow(r, matrix[r]);

		return mat;
	}

	/**
	 * Convert a vecmath GMatrix to a JAMA matrix.
	 * 
	 * @param gmatrix a vecmath GMatrix
	 * @return
	 */
	public static Matrix vecmathToJama(GMatrix gmatrix) {
		return new Matrix(vecmathToDouble(gmatrix));
	}

	/**
	 * Convert a vecmath GMatrix to a double[][] matrix.
	 * 
	 * @param gmatrix a vecmath GMatrix
	 * @return
	 */
	public static double[][] vecmathToDouble(GMatrix gmatrix) {
		
		double[][] matrix = new double[gmatrix.getNumRow()][];
		
		for (int r = 0; r < gmatrix.getNumRow(); r++) {
			matrix[r] = new double[gmatrix.getNumCol()];
			gmatrix.getColumn(r, matrix[r]);
		}
		
		return matrix;
	}
}
