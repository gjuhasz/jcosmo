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

import java.io.FileNotFoundException;
import java.util.HashMap;



/**
 * COSMO-SAC (MOPAC sigma profiles) activity model.
 *
 * @author Rafael de Pelegrini Soares
 * 
 */
public class COSMOPACMulti2 extends COSMOSACMulti {
	private String folder;
	private HashMap<String, COSMOSACCompound> compList;

	public String toString(){
		return "COSMO-SAC(MOPAC)";
	}
	
	public COSMOPACMulti2(){
		super(11, 5);
		
		// we use another averaging radius
		this.rav = COSMOSAC.RAV*1.1;
		
		folder = "mopac/";
//		folder = "moltest/";
//		folder = "mopAM1/";
//		folder = "mopAM1c/";
//		folder = "mopRM1/";
		
//		// article results
//		setResCorr(1);
//		setCHB(42700.7265672813);
//		setSigmaHB(0.0064);
//		setFpol(FPOL);
//		setIgnoreSG(false);
//		setAnorm(28.2);
//		setVnorm(66.69);
//		setAEff(7.5);
		
		// Estimating a base Fpol, Anorm and beta, later the other fpol's can be refined by groups
		// nonHB.csv (RAV*1.1), NP:74  COST:0.14125299308250233
		setBeta(1.7614059586590503);
		setCHB(0.0);
		setFpol(0.4487313504580587);
		setAnorm(68.84830662391695);
		
		// Now refining the fpol by using atom groups
		// idac/nonaqueous.csv AARD:0.30906048755539883 NP:166
		// idac/aqueous.csv AARD:0.5644266269440232 NP:80
		// Groups: Others, H-[X], O-H, [X]
		setFpol(0, 0, 0.5539202049925867);
		setFpol(0, 1, 0.3006546951718698);
		setFpol(0, 2, 0.44036655733523866);
		setFpol(0, 3, 0.9107567260382149);
		setFpol(1, 1, 0.2593918378709723);
		setFpol(1, 2, 0.39305951912580894);
		setFpol(1, 3, 0.522127832591752);
		setFpol(2, 2, 0.2829151539777929);
		setFpol(2, 3, 0.4273418078062349);
		setFpol(3, 3, 0.2348356781748074);
		
		// idac/nonaqueous.csv AARD:0.32142505824509987 NP:166
		// idac/aqueous.csv AARD:0.602308251715228 NP:82
		setBeta(1.442817285679101);
		setCHB(0.0);
		setAnorm(37.956541592840395);

		// refining using groups
		// idac/nonaqueous.csv AARD:0.17468664346591337 NP:166
		setFpol(0, 0, 0.5873804448893425);
		setFpol(0, 1, 0.4046809691780783);
		setFpol(0, 2, 0.45645447290586216);
		setFpol(0, 3, 0.8978885586648753);
		setFpol(0, 4, 0.040416710931485356);
		setFpol(1, 1, 0.019954367190134475);
		setFpol(1, 2, 0.7161438156521922);
		setFpol(1, 3, 2.796709799344939);
		setFpol(1, 4, 0.5517777396911522);
		setFpol(2, 2, 0.23687172097261877);
		setFpol(2, 3, 3.696452067701017);
		setFpol(2, 4, 3.8399582108072123);
		setFpol(3, 3, 0.006770465352492372);
		setFpol(3, 4, 0.01757967543675315);
		setFpol(4, 4, 0.2682355020982772);
		
		// idac/nonaqueous.csv AARD:0.2226227176807255 NP:166
		// idac/glycerol.csv AARD:1.1540519368661408 NP:7
		setFpol(0, 0, 0.6415731346889545);
		setFpol(0, 1, 0.20573703012675426);
		setFpol(0, 2, 0.48981671778781277);
		setFpol(0, 3, 0.9178808089795247);
		setFpol(0, 4, 0.3490764244862137);
		setFpol(1, 1, 0.03830351609984964);
		setFpol(1, 2, 0.8193485746870062);
		setFpol(1, 3, 0.5918798064940802);
		setFpol(1, 4, 0.17668904303788113);
		setFpol(2, 2, 0.5540332546953947);
		setFpol(2, 3, 1.172235878651001);
		setFpol(2, 4, 0.8130573774067755);
		setFpol(3, 3, 0.03501484559441628);
		setFpol(3, 4, 0.5961241288041279);
		setFpol(4, 4, 0.33574881273748186);		
		
		// Maybe we need more atom groups to get a better match
	}

	public void setComponents(COSMOSACCompound comps[]) throws Exception {
		this.comps = comps;
		this.ncomps = comps.length;

		this.VCOSMO = new double[ncomps];
		
		if(compList==null)
			compList = new HashMap<String, COSMOSACCompound>();

		SigmaProfileGenerator s = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.MOPAC,
				this.rav, nsegments);
		for (int i = 0; i < comps.length; i++) {
			
			String name = comps[i].name.replace(' ','_');
			String extension = ".cos";
			
			COSMOSACCompound c2 = compList.get(comps[i].name);
			if(c2!=null){
				comps[i].charge = c2.charge;
				this.VCOSMO[i] = comps[i].Vcosmo = c2.Vcosmo;
				comps[i].area = c2.area;
				comps[i].areaMulti = c2.areaMulti;
				continue;
			}
			
			try {
				s.parseFile(folder + name + extension);												
			} catch (FileNotFoundException e) {
				s.parseFile(folder + name.replace('-','_') + extension);
			}
			
			comps[i].charge = s.getChargeDensity();
			this.VCOSMO[i] = comps[i].Vcosmo = s.getVolume();
			
			// the multi-area
			double[] sigma1 = s.getAveragedChargeDensity();
			int[] atoms = s.getAtom();

			double[] area0 = s.getOriginalArea();
			double[] area1 = new double[area0.length];
			double[] area2 = new double[area0.length];
			double[] area3 = new double[area0.length];
			double[] area4 = new double[area0.length];

			comps[i].areaMulti = new double[ndescriptors][];
			
			MolParser molParser = new MolParser();
			molParser.parseFile(folder + name + ".mol");

			// lets filter the H-[N,O,...] atoms
			for (int m = 0; m < area0.length; m++) {
				int boundedType[] = {7, 8, 9, 17, 35, 53};
				if(molParser.matchType(atoms[m], 1, boundedType)){
					area1[m] += area0[m];
					area0[m] = 0;
				}
			}
			
			// lets filter the [N,O,...]-H atoms
			for (int m = 0; m < area0.length; m++) {
				int atomType[] = {7, 8, 9, 17, 35, 53};
//				int atomType[] = {8};
				if(molParser.matchType(atoms[m], atomType, 1)){
					area3[m] += area0[m];
					area0[m] = 0;
				}
			}
			
			// lets filter the C=C && C=-C atoms
			for (int m = 0; m < area0.length; m++) {
				if(molParser.matchBondType(atoms[m], 6, 2)){
					area4[m] += area0[m];
					area0[m] = 0;
				}
			}
			for (int m = 0; m < area0.length; m++) {
				if(molParser.matchBondType(atoms[m], 6, 3)){
					area4[m] += area0[m];
					area0[m] = 0;
				}
			}
			
			// lets filter the [N, O, ...] atoms
			for (int m = 0; m < area0.length; m++) {
				int atomType[] = {7, 8, 9, 17, 35, 53};
				if(molParser.matchType(atoms[m], atomType, 0)){
					area2[m] += area0[m];
					area0[m] = 0;
				}
			}
			
			s.simpleSorting(area0, sigma1);
			comps[i].areaMulti[0] = s.getSortedArea();
			s.simpleSorting(area1, sigma1);
			comps[i].areaMulti[1] = s.getSortedArea();
			s.simpleSorting(area2, sigma1);
			comps[i].areaMulti[2] = s.getSortedArea();
			s.simpleSorting(area3, sigma1);
			comps[i].areaMulti[3] = s.getSortedArea();
			s.simpleSorting(area4, sigma1);
			comps[i].areaMulti[4] = s.getSortedArea();
			
			compList.put(comps[i].name, comps[i]);
			
//			s.printProfile(System.out);
		}
		this.T = 300;

		z = new double[ncomps];
		for (int i = 0; i < ncomps; i++)
			z[i] = 1.0/ncomps;

		parametersChanged();
	}
}
