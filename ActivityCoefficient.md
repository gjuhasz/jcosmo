This page shows how to use JCosmoto calculate the activity coefficient of a binary mixture.

## Necessary procedure ##
1) Imports:
```
import br.ufrgs.enq.jcosmo.COSMOSAC;
import br.ufrgs.enq.jcosmo.COSMOSACCompound;
import br.ufrgs.enq.jcosmo.COSMOSACDataBase;
```
2) Creat variable COSMOSACDataBase:
```
COSMOSACDataBase db = COSMOSACDataBase.getInstance();
```
3) Select two componds of the database, the compound should be available in the JCosmo database and must be spelled correctly according to the database:
```
COSMOSACCompound c = db.getComp("compoundName");
```
4) Get of the database the "volume occupied by the molecular":
```
cavityVolume = c.Vcosmo;
```
5) Get of the database the sigma profiles:
```
sigma = c.sigma;
```
6) Creat variable COSMOSAC to calculate the activity coefficient:
```
COSMOSAC cosmosac = new COSMOSAC( );
```
7) Set parameters of calculation:
```
cosmosac.setParameters(cavityVolume, c1.charge, sigma);
```
8) Set temperature (Kelvin):
```
cosmosac.setTemperature(T);
```
9) Set composition (use array), the activity coefficient depends on the compositions:
```
double [] z = new double[2];
cosmosac.setComposition(z);
```
10) Use array to get ln(gamma):
```
double [] lnGamma = new double[2];
cosmosac.activityCoefficient(lnGamma);
```
**NOTE**: The JCosmo calculates the natural logarithm of the activity coefficient, to obtain the activity coefficient itself use the exponential.
```
gamma = Math.exp(lnGamma);
```

## Example Code ##

This is the simplest code:
```
COSMOSACDataBase db = COSMOSACDataBase.getInstance();
COSMOSACCompound c1 = db.getComp("water");
COSMOSACCompound c2 = db.getComp("ethanol");

double[] cavityVolume = new double[2];
cavityVolume[0] = c1.Vcosmo;
cavityVolume[1] = c2.Vcosmo;

double[][] sigma = new double[2][];
sigma[0] = c1.sigma;
sigma[1] = c2.sigma;

COSMOSAC cosmosac = new COSMOSAC();
cosmosac.setParameters(cavityVolume, c1.charge, sigma);
cosmosac.setTemperature(273.15); // T em Kelvin

double[] z = new double[2];
double[] lnGamma = new double[2];
z[0] = 0.50;
z[1] = 1 - z[0];

cosmosac.setComposition(z);
cosmosac.activityCoefficient(lnGamma);

double[] gamma = new double[2];
gamma[0] = Math.exp(lnGamma[0];
gamma[1] = Math.exp(lnGamma[1];
```