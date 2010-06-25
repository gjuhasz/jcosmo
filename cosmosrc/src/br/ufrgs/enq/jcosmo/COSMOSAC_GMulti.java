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
 * COSMO-SAC (GAMESS sigma profiles) activity model.
 *
 * @author Rafael de Pelegrini Soares
 * 
 */
public class COSMOSAC_GMulti extends COSMOPACMulti {
	public String toString(){
		return "COSMO-SAC(GAMESS)";
	}
	
	public COSMOSAC_GMulti(){
		type = SigmaProfileGenerator.FileType.GAMESS;
		extension = ".gout";
		
//		folder = "moltest/";
		folder = "gam6-31++G2d,p/";
//		folder = "gam6-31Gd/";
//		folder = "gamSTO3/";
		
		
		// nonHB
		// COST:0.2167219390934047 NP:68
		setBeta(2.050763884895556);
		setFpol(0.5513422645169235);
		
	
	}
}
