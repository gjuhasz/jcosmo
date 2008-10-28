package br.ufrgs.enq.jcosmo.test;



public class Antoine {
	public static double CalcPsat(double A, double B, double C, double T)
	{
		double Psat;
		Psat = Math.exp(A - (B/(T + C)));
		return Psat;
	}	
}
