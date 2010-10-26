package br.ufrgs.enq.jcosmo2;

import java.util.ArrayList;
import java.util.List;

import com.csvreader.CsvReader;

/**
 * Class representing an set of infinite dilution activity coefficient (IDAC) experiments.
 * 
 * @author rafael e renan
 *
 */
public class IDACExperiments {
	
	String filename;

	double AARD;

	private String modelClass;
	private int NP;
	
	List<COSMOSAC2> models = new ArrayList<COSMOSAC2>();

	List<Double> gammaInfs = new ArrayList<Double>();
	List<Double> temperatures = new ArrayList<Double>();
	List<Boolean> valid = new ArrayList<Boolean>();

	private boolean printGammas = false;
	
	public void setPrintGammas(boolean printGammas) {
		this.printGammas = printGammas;
	}
	public IDACExperiments(String filename, String modelClass, boolean leastSquares) throws Exception {
		this.filename = filename;
		this.modelClass = modelClass;
		
		loadContents();
	}
	public IDACExperiments(String filename, String modelClass) throws Exception {
		this(filename, modelClass, false);
	}
	
	public List<COSMOSAC2> getModels() {
		return models;
	}
	public List<Boolean> getValid() {
		return valid;
	}
	public List<Double> getMeasures() {
		return gammaInfs;
	}

	
	public void loadContents() throws Exception{
		COSMOSAC2 model = null;
		
		CsvReader reader = new CsvReader(filename);
		
		// infinite dilution activity calculation
		double z[] = new double[2];
		z[0] = 0; z[1] = 1-z[0];
		
		// skip the headers
		reader.readHeaders();
		
		while(reader.readRecord()){
			model = (COSMOSAC2) Class.forName(modelClass).newInstance();
			
			String mixtureString = reader.get(0);
			double T, gammaInf;
			try{
				T = Double.parseDouble(reader.get(1));
				gammaInf = Double.parseDouble(reader.get(2));
			}
			catch (NumberFormatException e) {
				System.err.println(e.toString() + ", line " + (NP + 1));
				continue;
			}
			
			String[] compNames = mixtureString.split("/");
			if(compNames.length != 2) {
				System.err.println("Typo on mixture(" + mixtureString + ", line " + (NP + 1));
				continue;
			}
			
			boolean valid = true;
			Compound comps[] = new Compound[2];

			comps[0] = new Compound();
			comps[0].name = compNames[0].toUpperCase();
			comps[1] = new Compound();
			comps[1].name = compNames[1].toUpperCase();
			
			try{
				model.setComponents(comps);
			}
			catch (Exception e) {
				System.err.println(e.toString());
				valid = false;
			}
			if(valid){
				model.setTemperature(T);
				model.setComposition(z);
			}
			
			models.add(model);
			temperatures.add(T);
			gammaInfs.add(gammaInf);
			this.valid.add(valid);
		}
		reader.close();
	}
	
	public void calcDeviations() throws Exception{
		double z[] = new double[2];
		double lnGamma[] = new double[2];
		double T;
		// infinite dilution activity calculation
		z[0] = 0; z[1] = 1-z[0];
		
		AARD = 0;
		NP = 0;
		for (int i = 0; i < models.size(); i++) {
			boolean valid = this.valid.get(i);
			COSMOSAC2 model = models.get(i);
			T = temperatures.get(i);
			double gammaInf = gammaInfs.get(i);
			
			if(!valid){
//				if(printGammas){
//					System.out.println("ignored" + ',' + gammaInf + ',' + "-");
//				}
				continue;
			}
			
			model.setTemperature(T);
			model.setComposition(z);
			model.activityCoefficientLn(lnGamma, 0);
			
			double lngammaInf = Math.log(gammaInf);
			double rd = Math.abs(lngammaInf - lnGamma[0]);
			
			if(printGammas){
				System.out.println(
//						model.getComps()[0].name + '/' + model.getComps()[1].name + ',' +
//						gammaInf + "," + Math.exp(lnGamma[0]));
						Math.log(gammaInf) + "," + lnGamma[0]);
			}
			
			AARD += rd;
			++NP;
		}
		AARD /= NP;
		System.out.println(filename + " AARD:" + AARD + " NP:" + NP);
	}
	
	public double getAARD(){
		return AARD;
	}
	
	public int getNP() {
		return NP;
	}
}
