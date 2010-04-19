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

import java.util.HashMap;


/**
 * PCM-SAC (PCM/GAMAESS sigma profiles) activity model.
 *
 * @author Rafael de Pelegrini Soares
 * 
 */
public class PCMSAC extends COSMOSAC {
	public String toString(){
		return "PCM-SAC";
	}
	
	// Store already loaded compounds
	static HashMap<String, COSMOSACCompound> compList = new HashMap<String, COSMOSACCompound>();  

	public PCMSAC(int numberOfSegments) {
		super(numberOfSegments);

		
		// COST:0.5015689274373564, NP:521, all IDAC
		// $PCMCAV ALPHA(1)=1.2, RSOLV=1.2
		setBeta(1.4466214455291058);
		setCHB(31161.034262510475);
		setSigmaHB(0.0019266011890448748);
		setFpol(FPOL);
		setIgnoreSG(false);
		setAnorm(49.26756101791564);
		setVnorm(66.69);
		
		
		// delta HB only a small subset COST:0.4686047788338226 NP:106
		setBeta(2.315692083990375);
		setCHB(3109.542519003334);
		setSigmaHB(0.011500675324712325);
		setFpol(0.6461642756815558);
		setIgnoreSG(false);
		setCoord(0.0797380264399844);
		setAnorm(145.1060693551804);
		setVnorm(66.69);
	}

	public PCMSAC() {
		this(51);
	}
	
	protected void calculeDeltaW_HB(){
		for(int m=0; m<nsegments; ++m){

			for(int n=0; n<nsegments; ++n){
				int ACC = n, DON = m;
				if(charge[m]>=charge[n]){
					ACC = m;
					DON = n;
				}

				// Hydrogen Bond effect:
				double hb = 0.0;
				// New HB from Paul M. Mathias, Shiang-Tai Lin, Yuhua Song, Chau-Chyun Chen, Stanley I. Sandler
				// AIChE Annual Meeting Indianapolis, IN, 3-8 November 2002
//				double sigmaHB = 0.018;
				hb = 0;
				if(charge[ACC]>sigmaHB/3 && charge[DON]<sigmaHB/3 && Math.abs(charge[ACC] - charge[DON]) > sigmaHB){
//					hb = Math.max(0.0, Math.abs(charge[ACC] - charge[DON]) - sigmaHB);
					hb = Math.abs(charge[ACC] - charge[DON]);
					hb = -(hb*hb);
				}
				
				// Klamt, Fluid Phase Equilib. 2000
//				double cHBT_c = 1.5;
				double cHBT = 1; // Math.max(0, 1 + cHBT_c * (298.15/T - 1));

				deltaW_HB[m][n] = cHB*cHBT* hb;
			}
		}
	}


	public void setComponents(COSMOSACCompound comps[]) throws Exception {
		this.ncomps = comps.length;
		this.comps = comps;

		this.VCOSMO = new double[ncomps];
		this.area = new double[ncomps][];

		SigmaProfileGenerator s = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS_PCM,
				SigmaProfileGenerator.RAV,
				nsegments);
		for (int i = 0; i < comps.length; i++) {
			
			COSMOSACCompound c2 = compList.get(comps[i].name);
			if(c2!=null){
				comps[i].charge = c2.charge;
				this.VCOSMO[i] = comps[i].Vcosmo = c2.Vcosmo;
				this.area[i] = comps[i].area = c2.area;
				continue;
			}
			
			String name = comps[i].name.replace(' ','_');
			String extension = ".pcm.gout";
			
			try {
				s.parseFile("mopac/" + name + extension);												
			} catch (Exception e) {
				s.parseFile("mopac/" + name.replace('-','_') + extension);
			}
			
			comps[i].charge = s.getChargeDensity();
			this.VCOSMO[i] = comps[i].Vcosmo = s.getVolume();
			this.area[i] = comps[i].area = s.getSortedArea();
			
			compList.put(comps[i].name, comps[i]);
		}
		this.T = 300;

		z = new double[ncomps];
		for (int i = 0; i < ncomps; i++)
			z[i] = 1.0/ncomps;

		parametersChanged();
	}
}
