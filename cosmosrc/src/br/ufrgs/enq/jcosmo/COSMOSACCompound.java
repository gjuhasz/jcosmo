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


/**
 * This class is basically a structure which holds compound data.
 * The user should not instantiate objects of this class directly,
 * use {@link COSMOSACDataBase} as a factory instead.
 * 
 * @author rafael
 *
 */
public class COSMOSACCompound {
	
	public int ID;
	public String name, formula, CAS, family, preOptTool;
	public double Vcosmo, T, lnPvap;
	public double []charge;
	public double []area;
	
	public String toString(){
		return name + ' ' + formula + ' ' + CAS;
	}
}
