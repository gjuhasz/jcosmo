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

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

/**
 * Create the sigma profile based on COSMO/GAMESS LOG file.
 * 
 * <p>This class will extract the segments information of a COSMO/GAMESS LOG file
 * to construct the sigma profile as needed by {@link COSMOSAC}.
 * The log file can be passed directly to the constructor, it will look for the
 * segment information.
 * 
 * <p>This code is based on the FORTRAN code Sigma-average available at
 * http://www.design.che.vt.edu/COSMO/. 
 * 
 * @author rafael
 *
 */
public class SigmaAverageGAMESS {
	private static final double REFF = 0.81764200000000;
	private static final double CHARGE_LOWER = -0.025;
	private static final double AU_ANGSTRON = 0.529177249;

	double[] x, y, z, area, SIGMA, SIGMANEW;
	double SP[], CHGDEN[];
	
	double increment;

	/**
	 * Creates the sigma profile given a COSMO/GAMESS output file with the default number of segments (51).
	 */
	public SigmaAverageGAMESS(String fileName) throws Exception {
		this(fileName, 51);
	}

	/**
	 * Creates the sigma profile given a COSMO/GAMESS output file.
	 * 
	 * @param fileName the COSMO/GAMESS output file
	 * @param sigmaPoints the number of elements on the resulting profile (usually 51)
	 * @throws Exception if there is a problem when reading the file
	 */
	public SigmaAverageGAMESS(String fileName, int sigmaPoints) throws Exception {
		readSegmentCharges(fileName);
		
		increment = -(CHARGE_LOWER*2)/(double)(sigmaPoints-1);
		// SETTING CHGDEN MATRIX
		SP = new double[sigmaPoints];
		CHGDEN = new double[sigmaPoints];

		for (int J = 0; J < sigmaPoints; J++) {
			SP[J]=0.0;
			CHGDEN[J] = CHARGE_LOWER + increment*(double)J;
		}
		
		averageCharges(sigmaPoints);
//		simpleSorting();
		normalSorting(SIGMANEW);
//		normalSorting(SIGMA);
	}

	/**
	 * @return the averaged sigma profile.
	 */
	public double[] getSigmaProfile(){
		return SP;
	}

	/**
	 * The charge density from 
	 * @return
	 */
	public double[] getChargeDensity(){
		return CHGDEN;
	}

	void readSegmentCharges(String filename) throws Exception{

		Scanner input = new Scanner(new File(filename));
		input.useLocale(Locale.US);

		// first lets try to find the "NUMBER OF SURFACE SEGMENTS IS"
		int cosmoSegments = 0;
		while(input.hasNext()){
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

		for (int i = 0; i < cosmoSegments; i++) {
			input.nextInt(); // segment id
			input.nextInt(); // atom
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

	void averageCharges(int sigmaPoints) {
		SIGMANEW = new double[x.length];

		double REFF2 = REFF*REFF;

		// BEGIN AVERAGING SURFACE CHARGES
		for(int J=0; J<x.length; ++J){
			double num = 0.0;
			double den =0.0;

			for(int K=0; K<x.length; ++K){
				double RADK2 = area[K]/Math.PI;
				
				double deltax = x[K]-x[J];
				double deltay = y[K]-y[J];
				double deltaz = z[K]-z[J];
				double DMN = Math.sqrt(deltax*deltax + deltay*deltay + deltaz*deltaz);
				
				double DMN2 = DMN*DMN;
				
				double exp = Math.exp(-DMN2/(RADK2+REFF2));

				num += SIGMA[K]*(RADK2 * REFF2)/(RADK2 + REFF2)*exp;
				den += (RADK2*REFF2)/(RADK2+REFF2)*exp;
			}
			SIGMANEW[J] = num/den;
		}
	}
	
	void simpleSorting(){
		// SIGMA PROFILE SORTING TAKEN FROM LIN DISSERTATION**
		for (int J = 0; J < SIGMANEW.length; J++) {
			int TMP = (int)((SIGMANEW[J]-CHGDEN[0])/increment);
			SP[TMP]+= area[J]*(CHGDEN[TMP+1]-SIGMANEW[J])/increment;
			SP[TMP+1]+= area[J]*(SIGMANEW[J]-CHGDEN[TMP])/increment;
		}
	}

	void normalSorting(double [] SIGMANEW) throws MathException {

		double stdDev = increment;

//		stdDev = Math.sqrt(org.apache.commons.math.stat.StatUtils.variance(SIGMANEW));
//		stdDev /= 6;

		// BEGIN AVERAGING SURFACE CHARGES
		for(int J=0; J<SIGMANEW.length; ++J){
			double mean = SIGMANEW[J];
			NormalDistribution dist = new NormalDistributionImpl(mean, stdDev);
			double Ais = area[J];
			double cumProb = 0;
			for (int i = 0; i < SP.length; i++) {
				double prob = 0;
				if(Math.abs(CHGDEN[i]-mean)< 6*stdDev){
					prob = dist.cumulativeProbability(CHGDEN[i]+increment/2);
					SP[i] += (prob-cumProb)*Ais;
					cumProb = prob;
				}
			}
		}
	}

	/**
	 * Prints the sigma profile to the given print stream
	 * @param out the print stream where to print to
	 */
	public void printProfile(PrintStream out){
		out.println("Sigma\tP(Sigma)");
		for (int i = 0; i < SP.length; i++) {
			out.println(CHGDEN[i] + "\t" + SP[i]);
		}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String filename;

		try {
			int nseg = 51;
			if(args.length > 0){
				filename = args[0];
				if(args.length>1){
					try{
						nseg = Integer.parseInt(args[1]);
					}
					catch (NumberFormatException e) {
						System.out.println("Invalid number of segments");
						return;
					}
				}			
			}
			else{
				Scanner in = new Scanner(System.in);
				System.out.print("Enter with the GAMESS LOG FILE: ");
				filename = in.nextLine();
			}
			SigmaAverageGAMESS a = new SigmaAverageGAMESS(filename, nseg);
			a.printProfile(System.out);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
