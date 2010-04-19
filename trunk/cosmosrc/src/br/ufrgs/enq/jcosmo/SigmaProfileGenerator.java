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
import java.io.PrintStream;
import java.util.Locale;
import java.util.Scanner;

/**
 * Create the sigma profile based on GAMESS LOG or MOPAC COS file.
 * 
 * <p>This class will extract the segments information from the file
 * to construct the sigma profile as needed by {@link COSMOSAC}.
 * 
 * <p>This code is based on the FORTRAN code Sigma-average available at
 * http://www.design.che.vt.edu/COSMO/. 
 * 
 * @author rafael
 *
 */
public class SigmaProfileGenerator {
	/** Default averaging radius */
	public static final double RAV = 0.81764200000000;
	
	private static final double CHARGE_LOWER = -0.025;
	protected static final double AU_ANGSTRON = 0.529177249;

	int[] atom,elem;
	double[] x, y, z, area, SIGMA, sigmaAveraged;
	double sortedArea[], CHGDEN[];
	
	double increment;
	double volume;
	double rav = RAV;
	
	public enum FileType {
		GAMESS,
		MOPAC,
		GAMESS_PCM,
	};
	
	FileType type;
	private int sigmaPoints;

	/**
	 * Creates the sigma profile given a COSMO/GAMESS output file with the default number of segments (51)
	 * and default averaging radius.
	 */
	public SigmaProfileGenerator(FileType type) throws Exception {
		this(type, RAV, 51);
	}
	
	/**
	 * Creates the sigma profile given a COSMO/GAMESS output and the number of segments,
	 * default averaging radius is used.
	 */
	public SigmaProfileGenerator(FileType type, int sigmaPoints) throws Exception {
		this(type, RAV, sigmaPoints);
	}
	
	/**
	 * Creates the sigma profile given a COSMO/GAMESS output file.
	 * 
	 * @param fileName the COSMO/GAMESS output file
	 * @param rav averaging radius
	 * @param sigmaPoints the number of elements on the resulting profile (usually 51)
	 */
	public SigmaProfileGenerator(FileType type, double rav, int sigmaPoints) {
		this.type = type;
		this.rav = rav;
		this.sigmaPoints = sigmaPoints;
	}


	/**
	 * Parses the given file to determine the sigma profile given a COSMO/GAMESS output file.
	 * 
	 * @param fileName the COSMO/GAMESS output file
	 * @throws Exception if there is a problem when reading the file
	 */
	public void parseFile(String fileName) throws Exception {
		// reset the volume
		volume = 0.0;
		switch (type) {
		case GAMESS:
			readSegmentChargesGAMESS(fileName);
			break;
		case GAMESS_PCM:
			readSegmentChargesGAMESS_PCM(fileName);
			break;
		case MOPAC:
			readSegmentChargesMOPAC(fileName);
			break;
		}
		
		increment = -(CHARGE_LOWER*2)/(double)(sigmaPoints-1);
		// SETTING CHGDEN MATRIX
		sortedArea = new double[sigmaPoints];
		CHGDEN = new double[sigmaPoints];

		for (int J = 0; J < sigmaPoints; J++) {
			sortedArea[J]=0.0;
			CHGDEN[J] = CHARGE_LOWER + increment*(double)J;
		}
		
		switch (type) {
		case GAMESS:
		case GAMESS_PCM:
//			simpleSorting(SIGMA);
			averageCharges();
//			averageCharges2();
			simpleSorting(sigmaAveraged);
			break;
		case MOPAC:
			averageCharges();
			simpleSorting(sigmaAveraged);
			break;
		}
	}

	/**
	 * @return the averaged and sorted area as a function of the charge density (sigma profile).
	 */
	public double[] getSortedArea(){
		return sortedArea;
	}

	/**
	 * @return the charge density
	 */
	public double[] getChargeDensity(){
		return CHGDEN;
	}

	void readSegmentChargesGAMESS(String filename) throws Exception{

		Scanner input = new Scanner(new File(filename));
		input.useLocale(Locale.US);

		// first lets try to find the "NUMBER OF SURFACE SEGMENTS IS"
		// volume information on "VOLUME INSIDE CONSTS:"
		int cosmoSegments = 0;
		while(input.hasNext()){
			if(volume==0 && input.next().equals("VOLUME") && input.next().equals("INSIDE") &&
					input.next().equals("CONSTS:")){
				volume = input.nextDouble();
			}
			if(input.next().equals("NUMBER") && input.next().equals("OF") &&
					input.next().equals("SURFACE") && input.next().equals("SEGMENTS") &&
					input.next().equals("IS")){
				cosmoSegments = input.nextInt();
				break;
			}
			input.nextLine();
		}

		if(cosmoSegments==0)
			throw new Exception("File " + filename + " does not have SEGMENTS information");

		x = new double[cosmoSegments];
		y = new double[cosmoSegments];
		z = new double[cosmoSegments];
		area = new double[cosmoSegments];
		SIGMA = new double[cosmoSegments];
		atom = new int[cosmoSegments];

		for (int i = 0; i < cosmoSegments; i++) {
			input.nextInt(); // segment id
			atom[i] = input.nextInt(); // atom
			x[i] = input.nextDouble()*AU_ANGSTRON;
			y[i] = input.nextDouble()*AU_ANGSTRON;
			z[i] = input.nextDouble()*AU_ANGSTRON;

			input.nextDouble(); // charge
			area[i] = input.nextDouble();
			SIGMA[i] = input.nextDouble();
			// input.nextLine(); // density (charge/area)
		}
		input.close();
	}

	void readSegmentChargesGAMESS_PCM(String filename) throws Exception{

		Scanner input = new Scanner(new File(filename));
		input.useLocale(Locale.US);

		// first lets try to find the "NUMBER OF SURFACE SEGMENTS IS"
		// volume information on "VOLUME INSIDE CONSTS:"
		int cosmoSegments = 0;
		while(input.hasNext()){
			if(cosmoSegments==0 && input.next().equals("TOTAL") && input.next().equals("NUMBER") &&
					input.next().equals("OF") && input.next().equals("TESSERAE=")){
				cosmoSegments = input.nextInt();
				
				// SURFACE AREA=    48.69984417(A**2)    CAVITY VOLUME=     30.55382014 (A**3)
				input.next();
				input.next();
				input.next();
				input.next();
				input.next();
				volume = input.nextDouble();
				
				input.nextLine();
				input.nextLine();
				input.nextLine();
				break;
			}
			input.nextLine();
		}

		if(cosmoSegments==0)
			throw new Exception("File " + filename + " does not have SEGMENTS information");

		x = new double[cosmoSegments];
		y = new double[cosmoSegments];
		z = new double[cosmoSegments];
		area = new double[cosmoSegments];
		SIGMA = new double[cosmoSegments];
		atom = new int[cosmoSegments];

		
		// first read the area information
		for (int i = 0; i < cosmoSegments; i++) {
			input.next(); // TESSERA
			input.next(); // SFERA
			area[i] = input.nextDouble()*AU_ANGSTRON*AU_ANGSTRON;
			x[i] = input.nextDouble()*AU_ANGSTRON;
			y[i] = input.nextDouble()*AU_ANGSTRON;
			z[i] = input.nextDouble()*AU_ANGSTRON;

			input.nextLine();
		}
		
		
		// then read the charge information
		while(input.hasNext()){
			if(input.next().equals("SFERA,TESSERA,X,Y,Z,CARICA")){
				break;
			}
		}
		for (int i = 0; i < cosmoSegments; i++) {
			input.next(); // TESSERA
			input.next(); // SFERA
			input.next(); // x
			input.next(); // y
			input.next(); // z
			SIGMA[i] = input.nextDouble();
			if(area[i]>0.0)
				SIGMA[i] /= area[i];

			input.nextLine();
		}
		
		input.close();
	}

	
	void readSegmentChargesMOPAC(String filename) throws Exception{

		Scanner input = new Scanner(new File(filename));
		input.useLocale(Locale.US);

		// first lets try to find the "NUMBER OF SURFACE SEGMENTS IS"
		int cosmoSegments = 0;
		while(input.hasNext()){
			if(volume==0 && input.next().equals("COSMO")
					 && input.next().equals("VOLUME")
					 && input.next().equals("=")
					){
				volume = input.nextDouble();
			}
			if(input.next().equals("SEGMENT")
					 && input.next().equals("DATA:")
					 && input.next().equals("NPS=")
					){
				cosmoSegments = input.nextInt();
				
				input.next();
				input.nextLine(); // skip the header file
				break;
			}
			input.nextLine();
		}

		if(cosmoSegments==0)
			throw new Exception("File " + filename + " does not have SEGMENTS information");

		atom = new int[cosmoSegments];
		elem = new int[cosmoSegments];
		x = new double[cosmoSegments];
		y = new double[cosmoSegments];
		z = new double[cosmoSegments];
		area = new double[cosmoSegments];
		SIGMA = new double[cosmoSegments];

		for (int i = 0; i < cosmoSegments; i++) {
			input.nextInt(); // segment id
			atom[i] = input.nextInt();
//			input.nextInt(); // elem
			elem[i] = input.nextInt(); // O=8, N=7, H=1
//			x[i] = input.nextDouble()*AU_ANGSTRON;
//			y[i] = input.nextDouble()*AU_ANGSTRON;
//			z[i] = input.nextDouble()*AU_ANGSTRON;
			x[i] = input.nextDouble();
			y[i] = input.nextDouble();
			z[i] = input.nextDouble();

			input.nextDouble(); // charge
			area[i] = input.nextDouble();
			SIGMA[i] = input.nextDouble();
			input.nextLine(); // skip potential
		}
		input.close();
	}

// from Mullins, Liu, Ghaderi, Fast, 2008
	void averageCharges() {
		sigmaAveraged = new double[x.length];

		double RAV2 = rav*rav;

		// BEGIN AVERAGING SURFACE CHARGES
		for(int J=0; J<x.length; ++J){
			double num = 0.0;
			double den = 0.0;

			for(int K=0; K<x.length; ++K){
				double RADK2 = area[K]/Math.PI;
				
				// only the same atom is averaged
//				if(atom!=null && atom[J]!=atom[K])
//					continue;
				
				double deltax = x[K]-x[J];
				double deltay = y[K]-y[J];
				double deltaz = z[K]-z[J];
				double DMN = Math.sqrt(deltax*deltax + deltay*deltay + deltaz*deltaz);
				
				double DMN2 = DMN*DMN;
				
				double exp = Math.exp(-DMN2/(RADK2+RAV2));
				
				double temp = ( (RADK2 * RAV2)/(RADK2 + RAV2) )*exp;
				num += SIGMA[K]*temp;
				den += temp;
//				num += SIGMA[K]*(RADK2 * REFF2)/(RADK2 + REFF2)*exp;
//				den += (RADK2*REFF2)/(RADK2+REFF2)*exp;
			}
			sigmaAveraged[J] = num/den;
		}
	}
	

// from Wang, Sandler, 2007 - using Aeff = 7.25
	void averageCharges2() {
		sigmaAveraged = new double[x.length];
		
		double Fdecay = 3.57;
		double Reff = rav;
		double Reff2 = Reff*Reff;

		// BEGIN AVERAGING SURFACE CHARGES
		for(int J=0; J<x.length; ++J){
			double num = 0.0;
			double den = 0.0;

			for(int K=0; K<x.length; ++K){
				double Rn2 = area[K]/Math.PI;
				
				double deltax = x[K]-x[J];
				double deltay = y[K]-y[J];
				double deltaz = z[K]-z[J];
				double Dmn = Math.sqrt(deltax*deltax + deltay*deltay + deltaz*deltaz);
				
				double Dmn2 = Dmn*Dmn;
				
				double exp = Math.exp(-Fdecay*(Dmn2/(Rn2+Reff2)));
				
				double temp = ( (Rn2 * Reff2)/(Rn2 + Reff2) )*exp;
				num += SIGMA[K]*temp;
				den += temp;
			}
			sigmaAveraged[J] = num/den;
		}
	}
	
	void simpleSorting(double sigmaAveraged[]){
		// SIGMA PROFILE SORTING TAKEN FROM LIN DISSERTATION**
		for (int J = 0; J < sigmaAveraged.length; J++) {
			int TMP = (int)((sigmaAveraged[J]-CHGDEN[0])/increment);
			
			if(TMP<0){
				sortedArea[0] += area[J];
				continue;
			}
			else if(TMP > (sortedArea.length-1)){
				sortedArea[sortedArea.length-1] += area[J];
				continue;
			}
			
			// Each point represents the center of an interval, so we distribute it
			// in the two adjacent points
			sortedArea[TMP]+= area[J]*(CHGDEN[TMP+1]-sigmaAveraged[J])/increment;
			sortedArea[TMP+1]+= area[J]*(sigmaAveraged[J]-CHGDEN[TMP])/increment;
		}
	}

	/**
	 * Prints the sigma profile to the given print stream
	 * @param out the print stream where to print to
	 */
	public void printProfile(PrintStream out){
		out.println("Sigma\tP(Sigma)");
		for (int i = 0; i < sortedArea.length; i++) {
			out.println(CHGDEN[i] + "\t" + sortedArea[i]);
		}
	}
	
	/**
	 * @return the cosmo cavity volume in A^3
	 */
	public double getVolume(){
		return volume;
	}
	
	public int[] getAtom() {
		return atom;
	}

	public int[] getElem() {
		return elem;
	}
}
