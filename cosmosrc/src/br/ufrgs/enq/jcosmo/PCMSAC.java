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
 * PCM-SAC (PCM/GAMAESS sigma profiles) activity model.
 *
 * @author Rafael de Pelegrini Soares
 * 
 */
public class PCMSAC extends COSMOSAC {
	public String toString(){
		return "PCM-SAC";
	}

	public PCMSAC(int numberOfSegments) {
		super(numberOfSegments);
		
		setBeta(1);
		setSigmaHB(0.02);
		
		// alcohol-water + water COST:0.4672116335560187
		setBeta(1.0);
		setCHB(18592.777394896664);
		setSigmaHB(0.012726877598985551);
		setFpol(0.35422711430612874);
		setIgnoreSG(false);
		setAnorm(24.470211102619594);
		setVnorm(66.69);
		
		// all IDAC COST:0.657490575262719
		setBeta(1.0);
		setCHB(45611.478357208674);
		setSigmaHB(0.0020111097044760917);
		setFpol(-0.2351557748422603);
		setIgnoreSG(false);
		setAnorm(30.550065238783525);
		setVnorm(66.69);
	}

	public PCMSAC() {
		this(51);
	}

	public void setComponents(COSMOSACCompound comps[]) throws Exception {
		this.ncomps = comps.length;

		this.VCOSMO = new double[ncomps];
		this.area = new double[ncomps][];

		SigmaProfileGenerator s = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS_PCM,
				SigmaProfileGenerator.RAV, nsegments);
		for (int i = 0; i < comps.length; i++) {
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
		}
		this.T = 300;

		z = new double[ncomps];
		for (int i = 0; i < ncomps; i++)
			z[i] = 1.0/ncomps;

		parametersChanged();
	}
}
