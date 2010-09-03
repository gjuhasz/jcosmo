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

import br.com.vrtech.thermolab.Comp;
import br.com.vrtech.thermolab.Const;
import br.com.vrtech.thermolab.DataBase;
import br.com.vrtech.thermolab.Mixture;
import br.com.vrtech.thermolab.Phase;
import br.com.vrtech.thermolab.activity.UNIFAC;
import br.ufrgs.enq.jcosmo.COSMOPACMulti;
import br.ufrgs.enq.jcosmo.COSMOSACCompound;



/**
 * COSMO-SAC (MOPAC sigma profiles) activity model.
 *
 * @author Rafael de Pelegrini Soares
 * 
 */
public class COSMOPACMultiUNIFAC extends COSMOPACMulti {
	protected UNIFAC unifac;
	double lnGammaUnifac[];

	public String toString(){
		return "COSMO-SAC(MOPAC+UNIFAC)";
	}
	
	public COSMOPACMultiUNIFAC() {
		setIgnoreSG(true);
//		setIgnoreResidual(true);
		
		DataBase.getInstance("jar:(thermolab.jar)thermolabdbUNIFAC");
		unifac = new UNIFAC();
		unifac.setIgnoreResidual(true);
		
//		setBeta(1.6140672703121623);
//		setFpol(0.5932863313076369);
//		
//		setSigmaHB(0);
//		setCHB(0);
//		setCHB(1, 2, 4319);
//		setCHB(1, 3, 4000);
	}
	
	@Override
	public void setComponents(COSMOSACCompound[] comps) throws Exception {
		super.setComponents(comps);
		
		DataBase db = DataBase.getInstance();
		Mixture mix = new Mixture();
		
		for (int i = 0; i < comps.length; i++) {
			Comp c = db.getComp(comps[i].name);
			if(c == null)
				throw new Exception("Component not found:" + comps[i].name);
			mix.addComp(c);
		}
		unifac.setMixture(mix);
		lnGammaUnifac = new double[mix.getNComps()];
	}
	
	public void activityCoefficientLn(double[] lnGama, int start){
		for (int i = 0; i < lnGama.length; i++)
			lnGama[i] = 0;
		super.activityCoefficientLn(lnGama, start);
		
		unifac.set(T, Const.P1atm, z, Phase.CHANGED_ALL);
		unifac.activityCoefficientLn(lnGammaUnifac, start);
		
		for (int i = 0; i < lnGama.length; i++)
			lnGama[i] += lnGammaUnifac[i];
	}
}
