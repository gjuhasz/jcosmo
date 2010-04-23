package br.ufrgs.enq.jcosmo.idac;

import br.ufrgs.enq.jcosmo.COSMOSAC;
import br.ufrgs.enq.jcosmo.COSMOSAC_G;
import br.ufrgs.enq.jcosmo.PCMSAC;

/**
 * Class representing an set of infinite dilution activity coefficient (IDAC) experiments.
 * 
 * @author rafael e renan
 *
 */
public class DiagonalCOSMOSAC {
	private static final long serialVersionUID = 1L;

	public static void main (String[] args) throws Exception{
//		String modelClass = COSMOSAC_G.class.getName();
//		String modelClass = COSMOSAC.class.getName();
//		String modelClass = COSMOPAC.class.getName();
		String modelClass = PCMSAC.class.getName();
		
		IDACDiagonal dig = new IDACDiagonal();
		
		dig.setTitle("IDAC for COSMOSAC model");
		
		dig.addIDACExperiments("idac/Alcohol-Water.csv", modelClass);
//		dig.addIDACExperiments("idac/Aldehyde-Water.csv", modelClass);
		dig.addIDACExperiments("idac/Alkane-Water.csv", modelClass);
		dig.addIDACExperiments("idac/Alkene-Water.csv", modelClass);
		dig.addIDACExperiments("idac/Alkyne-Water.csv", modelClass);
////		dig.addIDACExperiments("idac/AlkylHalide-Water.csv", modelClass);
		dig.addIDACExperiments("idac/Aromatic-Water.csv", modelClass);
////		dig.addIDACExperiments("idac/ArylHalide-Water.csv", modelClass);
//		dig.addIDACExperiments("idac/MultiringAromatics-Water.csv", modelClass);
//		dig.addIDACExperiments("idac/CarboxilicAcid-Water.csv", modelClass);
		dig.addIDACExperiments("idac/CycloAlkane-Water.csv", modelClass);
		dig.addIDACExperiments("idac/CycloAlkene-Water.csv", modelClass);
		dig.addIDACExperiments("idac/Ether-Water.csv", modelClass);
		dig.addIDACExperiments("idac/Ester-Water.csv", modelClass);
		dig.addIDACExperiments("idac/Ketone-Water.csv", modelClass);
//		dig.addIDACExperiments("idac/VinylHalide-Water.csv", modelClass);
		dig.addIDACExperiments("idac/Water.csv", modelClass);
//		
		dig.addIDACExperiments("idac/Alcohol-Alkane.csv", modelClass);
		dig.addIDACExperiments("idac/Alcohol-CycloAlkane.csv", modelClass);
//
		dig.addIDACExperiments("idac/Alkane-Alcohol.csv", modelClass);
//		dig.addIDACExperiments("idac/Alkane-AlkylHalide.csv", modelClass);
//		dig.addIDACExperiments("idac/Alkane-Amine.csv", modelClass);
//		dig.addIDACExperiments("idac/Alkane-CarboxilicAcid.csv", modelClass);
		dig.addIDACExperiments("idac/Alkane-Ketone.csv", modelClass);
		dig.addIDACExperiments("idac/Alkane-Phenol.csv", modelClass);
		
		dig.addIDACExperiments("idac/Alkene-Amine.csv", modelClass);
//		
//		dig.addIDACExperiments("idac/AlkylHalide-Alkane.csv", modelClass);
		dig.addIDACExperiments("idac/Amine-Alkane.csv", modelClass);
		dig.addIDACExperiments("idac/Aromatic-Alkane.csv", modelClass);
//		
		dig.addIDACExperiments("idac/CycloAlkane-Alcohol.csv", modelClass);
		dig.addIDACExperiments("idac/CycloAlkane-AlkylHalide.csv", modelClass);
		dig.addIDACExperiments("idac/CycloAlkane-Amine.csv", modelClass);
		dig.addIDACExperiments("idac/CycloAlkane-CarboxilicAcid.csv", modelClass);
		dig.addIDACExperiments("idac/CycloAlkane-Phenol.csv", modelClass);
//		
		dig.addIDACExperiments("idac/Ketone-Alcohol.csv", modelClass);
		dig.addIDACExperiments("idac/Ketone-Alkane.csv", modelClass);
//		
//		dig.addIDACExperiments("idac/Alkane-Alkane.csv", modelClass);
		
		// or just the families
//		dig.addIDACExperiments("idac/aqueous.csv", modelClass);
//		dig.addIDACExperiments("idac/nonaqueous.csv", modelClass);
		
		dig.showPlot();
	}
}