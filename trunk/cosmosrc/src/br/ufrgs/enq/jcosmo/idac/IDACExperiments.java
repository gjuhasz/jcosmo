package br.ufrgs.enq.jcosmo.idac;

import java.util.ArrayList;
import java.util.List;

import br.ufrgs.enq.jcosmo.COSMOPAC;
import br.ufrgs.enq.jcosmo.COSMOSAC;
import br.ufrgs.enq.jcosmo.COSMOSACCompound;
import br.ufrgs.enq.jcosmo.COSMOSACDataBase;

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
	
	List<COSMOSAC> models = new ArrayList<COSMOSAC>();

	List<Double> gammaInfs = new ArrayList<Double>();
	List<Double> temperatures = new ArrayList<Double>();

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
	
	public List<COSMOSAC> getModels() {
		return models;
	}
	public List<Double> getMeasures() {
		return gammaInfs;
	}

	
	public void loadContents() throws Exception{
		COSMOSAC model = null;
		
		CsvReader reader = new CsvReader(filename);
		
		// infinite dilution activity calculation
		double z[] = new double[2];
		z[0] = 0; z[1] = 1-z[0];
		
		// skip the headers
		reader.readHeaders();
		
		COSMOSACDataBase db = COSMOSACDataBase.getInstance();
		
		while(reader.readRecord()){
			model = (COSMOSAC) Class.forName(modelClass).newInstance();
			
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
			
			COSMOSACCompound comps[] = new COSMOSACCompound[2];
			if (modelClass == COSMOPAC.class.getName()){
				comps[0] = db.getComp("water");
				comps[0].name = compNames[0].toUpperCase();
				comps[1] = db.getComp("water");
				comps[1].name = compNames[1].toUpperCase();
			} 
			else {
				comps[0] = db.getComp(compNames[0].replace(' ', '-'));
				comps[1] = db.getComp(compNames[1].replace(' ', '-'));
				if(comps[0] == null){
					System.err.println("Component " + compNames[0] + " not found, it will be ignored.");
					continue;
				}
				if(comps[1] == null){
					System.err.println("Component " + compNames[1] + " not found, it will be ignored.");
					continue;
				}
			}
			
			try{
				model.setComponents(comps);
			}
			catch (Exception e) {
				System.err.println(e.toString());
				continue;
			}
			model.setTemperature(T);
			model.setComposition(z);
			
			models.add(model);
			temperatures.add(T);
			gammaInfs.add(gammaInf);
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
			COSMOSAC model = models.get(i);
			T = temperatures.get(i);
			double gammaInf = gammaInfs.get(i);
			
			model.setTemperature(T);
			model.setComposition(z);
			model.activityCoefficientLn(lnGamma, 0);
			
			double lngammaInf = Math.log(gammaInf);
			double rd = Math.abs(lngammaInf - lnGamma[0]);
			
			if(printGammas){
				System.out.println(model.getComps()[0].name + '/' + model.getComps()[1].name
						+ '\t' + gammaInf + '\t' + Math.exp(lnGamma[0]));
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
