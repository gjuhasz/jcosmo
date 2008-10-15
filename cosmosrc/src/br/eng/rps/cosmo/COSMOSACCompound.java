package br.eng.rps.cosmo;


/**
 * This class is basically a structure which holds compound data.
 * The user should not instantiate objects of this class directly,
 * use {@link COSMOSACDataBase} as a factory instead.
 * 
 * @author rafael
 *
 */
public class COSMOSACCompound {
	
	public int ID;
	public String name, formula, CAS, family, preOptTool;
	public double Vcosmo, T, lnPvap;
	public double []charge;
	public double []sigma;
	
	public String toString(){
		return name + ' ' + formula + ' ' + CAS;
	}
}
