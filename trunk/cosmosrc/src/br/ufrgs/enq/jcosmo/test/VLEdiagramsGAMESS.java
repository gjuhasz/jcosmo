/*
 * Copyright 2008 Rafael de Pelegrini Soares and Renan Pereira Gerber
 * 
 * This file is part of JCosmo.
 * 
 * JCosmo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * JCosmo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with JCosmo.  If not, see <http://www.gnu.org/licenses/>.
 */

package br.ufrgs.enq.jcosmo.test;

import java.awt.BorderLayout;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.Range;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import br.ufrgs.enq.jcosmo.COSMOSAC;
import br.ufrgs.enq.jcosmo.COSMOSACCompound;
import br.ufrgs.enq.jcosmo.COSMOSACDataBase;
import br.ufrgs.enq.jcosmo.SigmaProfileGenerator;

public class VLEdiagramsGAMESS extends JFrame {
	private static final long serialVersionUID = 1L;

	int n = 21;

	public JPanel calcAcetoneWater() throws Exception{
		super.setTitle("P vs x1");
		double T = 298.15;
		setLayout(new BorderLayout());

		COSMOSACDataBase db = COSMOSACDataBase.getInstance();
		COSMOSACCompound c1 = db.getComp("ethyl-acetate");
		COSMOSACCompound c2 = db.getComp("water");

		double[] cavityVolume = new double[2];
		cavityVolume[0] = c1.Vcosmo;
		cavityVolume[1] = c2.Vcosmo;

		double [][] sigma = new double[2][];

		sigma[0] = c1.sigma;
		sigma[1] = c2.sigma;
		
		SigmaProfileGenerator c1Sigma = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS,
				c1.name.toLowerCase() + ".log");
		SigmaProfileGenerator c2Sigma = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS,
				c2.name.toLowerCase() + ".log");
		sigma[0] = c1Sigma.getSigmaProfile();
		sigma[1] = c2Sigma.getSigmaProfile();

		COSMOSAC cosmosac = new COSMOSAC();
		cosmosac.setParameters(cavityVolume, c1.charge, sigma);
		cosmosac.setIgnoreSG(true);
		cosmosac.setSigmaHB(COSMOSAC.SIGMAHB*1.2);

		cosmosac.setTemperature(T);

		double[] x1 = new double [n];
		double[] x2 = new double [n];
		double[] gamma1 = new double [n];
		double[] gamma2 = new double [n];
		double [] z = new double[2];
		double [] lnGamma = new double[2];
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

		double[] Psat = new double[2];
		Psat[0] = 13000;
		Psat[1] = 3000;
		double[] P = calcPx(x1, x2, gamma1, gamma2,Psat);
		double[] y1 = calcY(x1, gamma1, Psat, P);

		XYPlot plot1;
		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries liq = new XYSeries("liquid");
		XYSeries vap = new XYSeries("vapor");
		XYSeries raoult = new XYSeries("Raoult's Law");
		for(int i=0; i<n; i++){
			liq.add(x1[i], P[i]);
			vap.add(y1[i], P[i]);
		}
		raoult.add(0, Psat[1]);
		raoult.add(1, Psat[0]);
		dataset.addSeries(liq);
		dataset.addSeries(vap);
		dataset.addSeries(raoult);

		JFreeChart chart = ChartFactory.createXYLineChart(null, 
				"Mole Fraction: x1, y1", "P/KPa", null,
				PlotOrientation.VERTICAL, true, true, false);
		plot1 = (XYPlot) chart.getPlot();
		plot1.getDomainAxis().setRange(new Range(0.0, 1.0));

		plot1.setDataset(dataset);

//		XYSplineRenderer r = new XYSplineRenderer();
//		BasicStroke stroke = new BasicStroke(2f);
//		r.setStroke(stroke);
//		plot1.setRenderer(r);
//		r.setBaseShapesVisible(false);

		ChartPanel chartPanel = new ChartPanel(chart);
		JPanel jp1 = new JPanel(new BorderLayout());
		jp1.add(chartPanel);

		add(jp1, BorderLayout.CENTER);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setSize(400,500);

		return jp1;
	}
	
	public JPanel calcAcAcidOctanol() throws Exception{
		double T = 293.15;
		setLayout(new BorderLayout());

		COSMOSACDataBase db = COSMOSACDataBase.getInstance();
		COSMOSACCompound c1 = db.getComp("acetic-acid");
		COSMOSACCompound c2 = db.getComp("1-octanol");

		double[] cavityVolume = new double[2];
		cavityVolume[0] = c1.Vcosmo;
		cavityVolume[1] = c2.Vcosmo;

		double [][] sigma = new double[2][];

		sigma[0] = c1.sigma;
		sigma[1] = c2.sigma;
		
		SigmaProfileGenerator c1Sigma = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS,
				c1.name.toLowerCase() + ".log");
		SigmaProfileGenerator c2Sigma = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS,
				c2.name.toLowerCase() + ".log");
		sigma[0] = c1Sigma.getSigmaProfile();
		sigma[1] = c2Sigma.getSigmaProfile();

		COSMOSAC cosmosac = new COSMOSAC();
		cosmosac.setParameters(cavityVolume, c1.charge, sigma);
		cosmosac.setIgnoreSG(true);
		cosmosac.setSigmaHB(COSMOSAC.SIGMAHB*1.2);

		cosmosac.setTemperature(T);

		double[] x1 = new double [n];
		double[] x2 = new double [n];
		double[] gamma1 = new double [n];
		double[] gamma2 = new double [n];
		double [] z = new double[2];
		double [] lnGamma = new double[2];
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

		double[] Psat = new double[2];
		Psat[0] = 1600;
		Psat[1] = 10;
		double[] P = calcPx(x1, x2, gamma1, gamma2,Psat);
		double[] y1 = calcY(x1, gamma1, Psat, P);

		XYPlot plot1;
		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries liq = new XYSeries("liquid");
		XYSeries vap = new XYSeries("vapor");
		XYSeries raoult = new XYSeries("Raoult's Law");
		for(int i=0; i<n; i++){
			liq.add(x1[i], P[i]);
			vap.add(y1[i], P[i]);
		}
		raoult.add(0, Psat[1]);
		raoult.add(1, Psat[0]);
		dataset.addSeries(liq);
		dataset.addSeries(vap);
		dataset.addSeries(raoult);

		JFreeChart chart = ChartFactory.createXYLineChart(null, 
				"Mole Fraction: x1, y1", "P/KPa", null,
				PlotOrientation.VERTICAL, true, true, false);
		plot1 = (XYPlot) chart.getPlot();
		plot1.getDomainAxis().setRange(new Range(0.0, 1.0));

		plot1.setDataset(dataset);

		ChartPanel chartPanel = new ChartPanel(chart);
		JPanel jp1 = new JPanel(new BorderLayout());
		jp1.add(chartPanel);

		add(jp1, BorderLayout.CENTER);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setSize(400,500);

		return jp1;
	}
	
	public JPanel calcMethanolOctanol() throws Exception{
		double T = 293.15;
		setLayout(new BorderLayout());

		COSMOSACDataBase db = COSMOSACDataBase.getInstance();
		COSMOSACCompound c1 = db.getComp("methanol");
		COSMOSACCompound c2 = db.getComp("1-octanol");

		double[] cavityVolume = new double[2];
		cavityVolume[0] = c1.Vcosmo;
		cavityVolume[1] = c2.Vcosmo;

		double [][] sigma = new double[2][];

		sigma[0] = c1.sigma;
		sigma[1] = c2.sigma;
		
		SigmaProfileGenerator c1Sigma = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS,
				c1.name.toLowerCase() + ".log");
		SigmaProfileGenerator c2Sigma = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS,
				c2.name.toLowerCase() + ".log");
		sigma[0] = c1Sigma.getSigmaProfile();
		sigma[1] = c2Sigma.getSigmaProfile();

		COSMOSAC cosmosac = new COSMOSAC();
		cosmosac.setParameters(cavityVolume, c1.charge, sigma);
		cosmosac.setIgnoreSG(true);
		cosmosac.setSigmaHB(COSMOSAC.SIGMAHB*1.2);

		cosmosac.setTemperature(T);

		double[] x1 = new double [n];
		double[] x2 = new double [n];
		double[] gamma1 = new double [n];
		double[] gamma2 = new double [n];
		double [] z = new double[2];
		double [] lnGamma = new double[2];
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

		double[] Psat = new double[2];
		Psat[0] = 350;
		Psat[1] = 470;
		double[] P = calcPx(x1, x2, gamma1, gamma2,Psat);
		double[] y1 = calcY(x1, gamma1, Psat, P);

		XYPlot plot1;
		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries liq = new XYSeries("liquid");
		XYSeries vap = new XYSeries("vapor");
		XYSeries raoult = new XYSeries("Raoult's Law");
		for(int i=0; i<n; i++){
			liq.add(x1[i], P[i]);
			vap.add(y1[i], P[i]);
		}
		raoult.add(0, Psat[1]);
		raoult.add(1, Psat[0]);
		dataset.addSeries(liq);
		dataset.addSeries(vap);
		dataset.addSeries(raoult);

		JFreeChart chart = ChartFactory.createXYLineChart(null, 
				"Mole Fraction: x1, y1", "P/KPa", null,
				PlotOrientation.VERTICAL, true, true, false);
		plot1 = (XYPlot) chart.getPlot();
		plot1.getDomainAxis().setRange(new Range(0.0, 1.0));

		plot1.setDataset(dataset);

		ChartPanel chartPanel = new ChartPanel(chart);
		JPanel jp1 = new JPanel(new BorderLayout());
		jp1.add(chartPanel);

		add(jp1, BorderLayout.CENTER);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setSize(400,500);

		return jp1;
	}

	
	public JPanel calcEthanolWater() throws Exception{
		double T = 343.15;
		setLayout(new BorderLayout());

		COSMOSACDataBase db = COSMOSACDataBase.getInstance();
		COSMOSACCompound c1 = db.getComp("ethanol");
		COSMOSACCompound c2 = db.getComp("water");

		double[] cavityVolume = new double[2];
		cavityVolume[0] = c1.Vcosmo;
		cavityVolume[1] = c2.Vcosmo;

		double [][] sigma = new double[2][];

		sigma[0] = c1.sigma;
		sigma[1] = c2.sigma;
		
		SigmaProfileGenerator c1Sigma = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS,
				c1.name.toLowerCase() + ".log");
		SigmaProfileGenerator c2Sigma = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS,
				c2.name.toLowerCase() + ".log");
		sigma[0] = c1Sigma.getSigmaProfile();
		sigma[1] = c2Sigma.getSigmaProfile();

		COSMOSAC cosmosac = new COSMOSAC();
		cosmosac.setParameters(cavityVolume, c1.charge, sigma);
		cosmosac.setIgnoreSG(true);
		cosmosac.setSigmaHB(COSMOSAC.SIGMAHB*1.2);

		cosmosac.setTemperature(T);

		double[] x1 = new double [n];
		double[] x2 = new double [n];
		double[] gamma1 = new double [n];
		double[] gamma2 = new double [n];
		double [] z = new double[2];
		double [] lnGamma = new double[2];
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

		double[] Psat = new double[2];
		Psat[0] = 72388;
		Psat[1] = 32760;
		double[] P = calcPx(x1, x2, gamma1, gamma2,Psat);
		double[] y1 = calcY(x1, gamma1, Psat, P);

		XYPlot plot1;
		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries liq = new XYSeries("liquid");
		XYSeries vap = new XYSeries("vapor");
		XYSeries raoult = new XYSeries("Raoult's Law");
		for(int i=0; i<n; i++){
			liq.add(x1[i], P[i]);
			vap.add(y1[i], P[i]);
		}
		raoult.add(0, Psat[1]);
		raoult.add(1, Psat[0]);
		dataset.addSeries(liq);
		dataset.addSeries(vap);
		dataset.addSeries(raoult);

		JFreeChart chart = ChartFactory.createXYLineChart(null, 
				"Mole Fraction: x1, y1", "P/KPa", null,
				PlotOrientation.VERTICAL, true, true, false);
		plot1 = (XYPlot) chart.getPlot();
		plot1.getDomainAxis().setRange(new Range(0.0, 1.0));

		plot1.setDataset(dataset);

		ChartPanel chartPanel = new ChartPanel(chart);
		JPanel jp1 = new JPanel(new BorderLayout());
		jp1.add(chartPanel);

		add(jp1, BorderLayout.CENTER);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setSize(400,500);

		return jp1;
	}
	
	
	public void biEq() throws Exception{
		JTabbedPane tabbedPane = new JTabbedPane();
		
		
//		tabbedPane.addTab("Ehtyl-Acetate x Water, T=298.15",calcAcetoneWater());
		tabbedPane.addTab("Acetic Acid x Octanol, T=293.15",calcAcAcidOctanol());
		tabbedPane.addTab("Ethanol x Water, T=343.15",calcEthanolWater());

		
		add(tabbedPane, BorderLayout.CENTER);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setSize(400,500);
	}

	private double[] calcPx(double[] x1, double[] x2, double[] gamma1,
			double[] gamma2, double[] Psat){
		double[] P = new double[n];		
		for (int k=0; k<n; k++){
			P[k] = x1[k]*gamma1[k]*Psat[0] + x2[k]*gamma2[k]*Psat[1];}
		return P;
	}

	private double[] calcY(double[] x1, double[] gamma1, double[] Psat, double[] P){
		double[] y1 = new double[n];
		for (int k=0; k<y1.length; k++){
			y1[k] = x1[k]*gamma1[k]*Psat[0]/P[k];
		}
		return y1;
	}

	public static void main (String[] args) throws Exception	{
		VLEdiagramsGAMESS test = new VLEdiagramsGAMESS();
		test.biEq();
		test.setVisible(true);

	}
}
