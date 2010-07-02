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
//		folder = "gam6-31++G2d,p/";
//		folder = "gamSTO3/";
		
		
		setAnorm(ANORM);
		setCHB(0);
		
		// nonHB
		// COST:0.21537997778086931 NP:68
		folder = "gam6-31Gd/";
		setBeta(2.139335747070768);
		setFpol(0.5204917405729581);
		// nonaqueous
		// COST:0.29569685222821074 NP:152
		setCHB(1, 1, 1700);
		setCHB(1, 2, 0);
		
		// nonHB
		// COST:0.3818185937017241 NP:68
		folder = "gam6-31Gd/";
		setBeta(1.0);
		setFpol(0.8965470231369211);
		// COST:0.29569685222821074 NP:152
		setCHB(1, 1, 1990);
		setCHB(1, 2, 0);
		
//		// aqueous
//		// 
//		setCHB(1, 4, 1700);
//		setCHB(3, 1, 1700);
//		setCHB(3, 2, 1700);
//		setCHB(3, 4, 1700);
		
	
	}
}
