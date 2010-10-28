package br.ufrgs.enq.jcosmo.test;

import java.io.File;

import org.apache.commons.math.stat.regression.SimpleRegression;

import br.ufrgs.enq.jcosmo.SigmaProfileGenerator;
import br.ufrgs.enq.jcosmo.SigmaProfileGenerator.FileType;

/**
 * Diagonal for averaging energy effects.
 * 
 * @author rafael
 *
 */
public class SigmaOrthogonal {

	public static void main(String[] args) {
		// Configuration
		String folder = "profiles/RM1_1.18/";
		String extension = ".cos";
		FileType type = SigmaProfileGenerator.FileType.MOPAC;
		double rav = 1.5;
		double rav2 = 2*rav;
		
//		String folder = "profiles/gamess/";
//		String extension = ".gout";
//		FileType type = SigmaProfileGenerator.FileType.GAMESS;
//		double rav = 1.5;
//		double rav2 = 1.5*rav;

		SimpleRegression line = new SimpleRegression();
		
		File folderObj = new File(folder);
		String names[] = folderObj.list();
		
		for(String name : names){
			if(!name.endsWith(extension))
				continue;
			
			name = name.substring(0, name.length() - extension.length());
			String fileName = folder + name + extension;
			
			System.out.println("Adding " + name + "...");

			SigmaProfileGenerator sigmaParser;

			double []sigma;
			double []sigma2;
			sigmaParser = new SigmaProfileGenerator(type);
			try {
				sigmaParser.parseFile(fileName, rav);
				sigma = sigmaParser.getAveragedChargeDensity();
				
				sigmaParser.parseFile(fileName, rav2);
				sigma2 = sigmaParser.getAveragedChargeDensity();
			} catch (Exception e) {
				e.printStackTrace();
				continue;
			}

			for (int i = 0; i < sigma.length; i++) {
				line.addData(sigma[i], sigma2[i]);
			}
		}

		System.out.println("Slope = " + line.getSlope() + " Intercept = " + line.getIntercept() + ", R2 = " + line.getR());
	}
}
