package br.ufrgs.enq.jcosmo.idac;

import br.ufrgs.enq.jcosmo.COSMOPACMulti;
import br.ufrgs.enq.jcosmo.COSMOPACMultiAtom;
import br.ufrgs.enq.jcosmo.COSMOSAC_GMulti;
import br.ufrgs.enq.jcosmo.COSMOSAC_GMultiAtom;
import br.ufrgs.enq.jcosmo.PCMSACMulti;

/**
 * Class representing an set of infinite dilution activity coefficient (IDAC) experiments.
 * 
 * @author rafael e renan
 *
 */
public class DiagonalCOSMOSACMultiFamilies {
	private static final long serialVersionUID = 1L;

	public static void main (String[] args) throws Exception{
//		String modelClass = COSMOPACMulti.class.getName();
		String modelClass = COSMOPACMultiAtom.class.getName();
//		String modelClass = PCMSACMulti.class.getName();
//		String modelClass = COSMOSAC_GMulti.class.getName();
//		String modelClass = COSMOSAC_GMultiAtom.class.getName();
		
		IDACDiagonalMulti dig = new IDACDiagonalMulti();
		
		dig.setTitle("IDAC for COSMOSAC model");
		
		// or just the families
//		dig.addIDACExperiments("idac/nonHB.csv", modelClass);
		dig.addIDACExperiments("idac/aqueous.csv", modelClass);
		dig.addIDACExperiments("idac/nonaqueous.csv", modelClass);
		
		dig.showPlot();
	}
}