/*
 * Copyright 2008 Rafael de Pelegrini Soares and Renan Pereira Gerber
 * 
 * This file is part of JCosmo.
 * 
 * JCosmo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * JCosmo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with JCosmo.  If not, see <http://www.gnu.org/licenses/>.
 */

package br.ufrgs.enq.jcosmo;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Locale;
import java.util.Scanner;

/**
 * Parses the contents of a MOL file and can return the atom connectivity
 * and bond type.
 * 
 * @author rafael
 *
 */
public class MolParser {
	
	int bondAtom1[];
	int bondAtom2[];
	int bondType[];
	int elementType[];
	
	/**
	 * Creates a new parser.
	 * @see #parseFile(String)
	 */
	public MolParser() {
	}
	
	/**
	 * Parse the given MOL file name.
	 * 
	 * @see #getBondAtom1()
	 * @see #getBondAtom2()
	 * @see #getBondType()
	 * 
	 * @throws Exception if there is a problem when reading the file
	 */
	public void parseFile(String filename) throws FileNotFoundException, Exception {

		Scanner input = new Scanner(new File(filename));
		input.useLocale(Locale.US);
		
		// first three lines are for comments, etc.
		input.nextLine();
		input.nextLine();
		input.nextLine();
		
		int natoms = input.nextInt();
		int nbonds = input.nextInt();
		input.nextLine();
		
		elementType = new int[natoms];
		for (int i = 0; i < natoms; i++) {
			input.next(); // x
			input.next(); // y
			input.next(); // z
			String atom = input.next();
			if(atom.equals("H"))
				elementType[i] = 1;
			if(atom.equals("He"))
				elementType[i] = 2;
			if(atom.equals("C"))
				elementType[i] = 6;
			if(atom.equals("N"))
				elementType[i] = 7;
			if(atom.equals("O"))
				elementType[i] = 8;
			if(atom.equals("F"))
				elementType[i] = 9;
			if(atom.equals("S"))
				elementType[i] = 16;
			if(atom.equals("Cl"))
				elementType[i] = 17;
			if(atom.equals("Br"))
				elementType[i] = 35;
			if(atom.equals("I"))
				elementType[i] = 53;
			// FIXME more atom types here
			input.nextLine();
		}
		
		bondAtom1 = new int[nbonds];
		bondAtom2 = new int[nbonds];
		bondType = new int[nbonds];
		
		for (int i = 0; i < nbonds; i++) {
			bondAtom1[i] = input.nextInt();
			bondAtom2[i] = input.nextInt();
			bondType[i] = input.nextInt();
			
			input.nextLine();
		}
		
	}

	/**
	 * @return the first bond atom atom list
	 */
	public int[] getBondAtom1() {
		return bondAtom1;
	}
	/**
	 * @return the second bond atom atom list
	 */
	public int[] getBondAtom2() {
		return bondAtom2;
	}

	/**
	 * @return the bond type list
	 */
	public int[] getBondType() {
		return bondType;
	}
	/**
	 * @return the element type list
	 */
	public int[] getElementType() {
		return elementType;
	}
}
