## Introduction ##

For each component in a binary mixture, one could make a vapor-liquid equilibrium diagram. Such a diagram would graph liquid/vapor mole fraction on a horizontal axis and pressure on a vertical axis, with constant temperature. In such VLE diagrams, liquid mole fractions for components 1 and 2 can be represented as x1 and x2 respectively, and vapor mole fractions of the corresponding components are commonly represented as y1 and y2. Using the modified Raoult's Law for not ideal liquids, through the Activity Coefficient (jcosmo) , it's possible to plot a vapor-liquid equilibrium diagrams.

### Raoult's Law ###

Raoult's law states: the vapor pressure of an ideal solution is dependent on the vapor pressure of each chemical component and the mole fraction of the component present in the solution. Changing the Raoult's Law using the activity coefficient to correct the liquid expression (right size) to a real liquid, enough to:

P.y<sub>i</sub> = x<sub>i</sub>.gamma<sub>i</sub>.P<sup>sat</sup><sub>i</sub> , (i = 1, 2, 3, ..., N)

### Antoine Equation ###

The Antoine equation is a vapor pressure equation and describes the relation of the saturated vapor pressure and the temperature for pure components. The equation uses three trial parameters A, B and C.

ln P<sup>sat</sup>/KPa = A - B / (T/°C + C)

## Test with jcosmo ##

We test two binary mixtures (Ethanol + Toluene and Methanol + Methil-Acetate) to check if the model represented the azeotrope. First was created compositions arrays of liquid (x1 and x2) ranging 0.05 and for each pair was calculated the activity coefficients and stored in another arrays (how to use jcosmo to calculate the activity coefficient
[here](http://code.google.com/p/jcosmo/wiki/ActivityCoefficient)).
```
int n = 21;
double[] x1 = new double [n];
double[] x2 = new double [n];
double[] gamma1 = new double [n];
double[] gamma2 = new double [n];
double[] z = new double[2];
double[] lnGamma = new double[2];
z[0] = 0.00;
int j = 0;
while(z[0] < 1.0001){
	z[1] = 1-z[0];
	x1[j] = z[0];
	x2[j] = z[1];
	cosmosac.setComposition(z);
	cosmosac.activityCoefficient(lnGamma);
	gamma1[j] = Math.exp(lnGamma[0]);
	gamma2[j] = Math.exp(lnGamma[1]);
	z[0] += 0.05;
	j++;
}
```

The saturated vapor pressures were obtained from the Antoine equation, so it's necessary to enter the experimental parameters of the equation (A, B and C) and the temperature you want to calculate. The parameters to Antoine equation were taken from _Introduction to chemist engineering thermodynamics 7th ed_. For example to Ethanol/Toluene mixture:
```
public class Antoine {
	public static double CalcPsat(double A, double B, double C, double T)
	{
		double Psat;
		Psat = Math.exp(A - (B/(T + C)));
		return Psat;
	}	
}
```
```
double T = 65; // °C
double[][] parAntoine = new double[3][3];
parAntoine[0][0] = 16.8958; // parameter A to Ethanol
parAntoine[0][1] = 3795.17; // parameter B to Ethanol
parAntoine[0][2] = 230.918; // parameter C to Ethanol
parAntoine[1][0] = 13.9320; // parameter A to Toluene
parAntoine[1][1] = 3056.96; // parameter B to Toluene
parAntoine[1][2] = 217.625; // parameter C to Toluene
Psat[0] = Antoine.CalcPsat(parAntoine[0][0], parAntoine[0][1], parAntoine[0][2], T);
Psat[1] = Antoine.CalcPsat(parAntoine[1][0], parAntoine[1][1], parAntoine[1][2], T);
```

Then you use the modified Raoult's Law to calculate the pressure for each composition and the vapor compositions.
```
double[] P = new double[n];
for (int k=0; k<n; k++)
	P[k] = x1[k]*gamma1[k]*Psat[0] + x2[k]*gamma2[k]*Psat[1];
```
```
double[] y1 = new double[n];
for (int k=0; k<21; k++)
	y1[k] = x1[k]*gamma1[k]*Psat[0]/P[k];
```

We plot the three curves (x1<sub>i</sub> vs. P<sub>i</sub>, y1<sub>i</sub> vs. P<sub>i</sub> and Raoult's Law) in the same diagram. The model used represented very well the Methanol/Methil-Acetate mixture, it's qualitatively correct, including the azeotrope. For the mixture Ethanol/Toluene it was not so good, although it has come close, not represented the azeotrope.


![http://i307.photobucket.com/albums/nn318/renanpgb/VLEdiagrams.jpg](http://i307.photobucket.com/albums/nn318/renanpgb/VLEdiagrams.jpg)