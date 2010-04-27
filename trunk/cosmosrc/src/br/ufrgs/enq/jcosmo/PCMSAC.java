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
public class PCMSAC extends COSMOPACM {
	public String toString(){
		return "PCM-SAC";
	}
	
	// Store already loaded compounds
	static HashMap<String, COSMOSACCompound> compList;  

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
		
		
		// delta HB only a small subset COST:0.39132272581177113 NP:163
		setBeta(1.5964708030478103);
		setCHB(1996.1201800255446);
		setSigmaHB(0.010166006038485578);
		setSigmaHB2(2.0322754257416418);
		setSigmaHB3(0.0);
		setFpol(0.6917);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(166.57659724842748);
		setVnorm(66.69);

		// only non-aqueous COST:0.4223700007983835 NP:227
		setBeta(2.1162691139604073);
		setCHB(959.0481495067659);
		setSigmaHB(0.016984297869683283);
		setSigmaHB2(2.6228348657791436);
		setSigmaHB3(0.0);
		setFpol(0.6917);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(134.64167457551045);
		setVnorm(66.69);
		
		
		// no C=O systems, COST:0.47509022374529086 NP:297
		setBeta(1.753295948035097);
		setCHB(1481.5019472726462);
		setSigmaHB(0.011150936306303545);
		setSigmaHB2(1.8584613097892642);
		setSigmaHB3(0.0);
		setFpol(0.6917);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(161.17891465949677);
		setVnorm(66.69);

		// Estimation of Fpol, Coord and Anorm for a very small subset:
//		idac/Alkane-Alkane.csv AARD:0.011508424352655757 NP:3
//		idac/Aromatic-Alkane.csv AARD:0.019907065297353543 NP:6
		setBeta(1);
		setCHB(0.0);
		setSigmaHB(0.0080);
		setSigmaHB2(2.0);
		setSigmaHB3(0.0);
		setFpol(0.7);
		setIgnoreSG(false);
		setCoord(7.5);
		setAnorm(80.0);
		setVnorm(66.69);
	}

	public PCMSAC() {
		this(51);
	}
	
	
	protected void calculeDeltaW(){
		double chargemn = 0;
		// polarizability corrections
		double fm = 1;
		double fn = 1;
		for(int m=0; m<nsegments; ++m){
			// initialize all SEGGAMMA (reused between calculations)
			SEGGAMMA[m] = 1.0;

			for(int n=0; n<nsegments; ++n){
				chargemn = charge[m]*fm+charge[n]*fn;
				deltaW[m][n] = (alphaPrime/2.0)*chargemn*chargemn;
			}
		}
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
				if(charge[ACC]>Math.abs(sigmaHB/sigmaHB2) && charge[DON]<-Math.abs(sigmaHB/sigmaHB2)
						&& Math.abs(charge[ACC] - charge[DON]) > Math.abs(sigmaHB)){
//					hb = Math.max(0.0, Math.abs(charge[ACC] - charge[DON]) - sigmaHB);
					hb = Math.abs(charge[ACC] - charge[DON]);
					hb = (hb*hb);
				}
				
				// Klamt, Fluid Phase Equilib. 2000
//				double cHBT_c = 1.5;
				double cHBT = 1; // Math.max(0, 1 + cHBT_c * (298.15/T - 1));

				deltaW_HB[m][n] = -cHB*cHBT*hb;
			}
		}
	}


	public void setComponents(COSMOSACCompound comps[]) throws Exception {
		this.ncomps = comps.length;
		this.comps = comps;

		this.VCOSMO = new double[ncomps];
		this.area = new double[ncomps][];
		
		if(compList==null)
			compList = new HashMap<String, COSMOSACCompound>();

		SigmaProfileGenerator s = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS_PCM,
				this.rav, nsegments);
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
			String folder = "moltest/";
			
			try {
				s.parseFile(folder + name + extension);												
			} catch (Exception e) {
				s.parseFile(folder + name.replace('-','_') + extension);
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
