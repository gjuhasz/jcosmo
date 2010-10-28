package br.ufrgs.enq.jcosmo2;

import br.ufrgs.enq.jcosmo.SigmaProfileGenerator;

public class COSMOSAC2_G extends COSMOSAC2 {
	
	public COSMOSAC2_G(){
		folder = "profiles/gamess/";
		extension = ".gout";
		type = SigmaProfileGenerator.FileType.GAMESS;
		
		rav = 1.5;
		rav2 = 1.5*rav;
		f_ortho = 0.7241079673777269;
		fcorr = -4.0;
//		skipAverage = true;
		
		setBeta(1.128202469186924);
		setFpol(0.8370064312795129);
	}

}
