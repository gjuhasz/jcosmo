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
	static HashMap<String, COSMOSACCompound> compList;  

	public PCMSAC(int numberOfSegments) {
		super(numberOfSegments);
		
		this.rav = RAV;

		// Estimation of Beta, Fpol, and Anorm: COST:0.08362529523602237, NP=109
//		idac/Alkane-Alkane.csv AARD:0.10012005066331824 NP:7
//		idac/Alkane-AlkylHalide.csv AARD:0.19950676462113945 NP:5
//		idac/Alkane-Ketone.csv AARD:0.1317831643908708 NP:21
//		idac/AlkylHalide-Alkane.csv AARD:0.10060078630976194 NP:32
//		idac/Aromatic-Alkane.csv AARD:0.02434364712056599 NP:39
//		idac/CycloAlkane-AlkylHalide.csv AARD:0.09614768531124265 NP:5
		setBeta(1.5313645197328083);
		setCHB(0.0);
		setSigmaHB(0.005949735966583021);
		setSigmaHB2(0.0);
		setSigmaHB3(1.0);
		setFpol(0.7323758674505356);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(361.45773761532666);
		setVnorm(66.69);
		
		setSigmaHB(0.0085);
		setCHB(80000);
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
				deltaW[m][n] = (fpol*alpha/2.0)*chargemn*chargemn;
			}
		}
	}
	
	protected void calculeDeltaW_HB2(){
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
				
//				if(charge[ACC]>0 && charge[DON]<0){
//					hb = Math.abs(charge[ACC] - charge[DON]);
//					hb = Math.pow(hb, sigmaHB3);
//				}
				
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
			String folder = "mopac/";
			
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
