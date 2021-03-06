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

import java.awt.BasicStroke;
import java.awt.BorderLayout;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYSplineRenderer;
import org.jfree.data.Range;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import br.ufrgs.enq.jcosmo.COSMOPAC;
import br.ufrgs.enq.jcosmo.COSMOSAC;
import br.ufrgs.enq.jcosmo.COSMOSACCompound;
import br.ufrgs.enq.jcosmo.COSMOSACDataBase;

public class VLEdiagrams extends JFrame {
	private static final long serialVersionUID = 1L;

	int n = 21;
	
	@SuppressWarnings("deprecation")
	public JPanel calcEthTol() throws Exception{
		super.setTitle("P vs x1");
		double T = 60;
		setLayout(new BorderLayout());

		COSMOSACDataBase db = COSMOSACDataBase.getInstance();
		
		COSMOSACCompound comps[] = new COSMOSACCompound[2];
		comps[0] = db.getComp("ethanol");
		comps[1] = db.getComp("toluene");		
		
		COSMOSAC cosmosac = new COSMOSAC();
		cosmosac.setComponents(comps);

		cosmosac.setTemperature(T+273.15);

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
		
		double[][] parAntoine = new double[3][3];
		parAntoine[0][0] = 16.8958;
		parAntoine[0][1] = 3795.17;
		parAntoine[0][2] = 230.918;
		parAntoine[1][0] = 13.9320;
		parAntoine[1][1] = 3056.96;
		parAntoine[1][2] = 217.625;

		double[] Psat = pSat(parAntoine, T);
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
		plot1.getRangeAxis().setRange(new Range(15.0, 50.0));
		
		plot1.setDataset(dataset);

		XYSplineRenderer r = new XYSplineRenderer();
		BasicStroke stroke = new BasicStroke(2f);
		r.setStroke(stroke);
		plot1.setRenderer(r);
		r.setBaseShapesVisible(false);
		
		ChartPanel chartPanel = new ChartPanel(chart);
		JPanel jp = new JPanel(new BorderLayout());
		jp.add(chartPanel);
						
		add(jp, BorderLayout.CENTER);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setSize(400,500);
		
		return jp;
	}
	
	@SuppressWarnings("deprecation")
	public JPanel calcMethanolMethilAcetate() throws Exception{
		super.setTitle("P vs x1");
		double T = 25;

		COSMOSACDataBase db = COSMOSACDataBase.getInstance();
		
		COSMOSACCompound comps[] = new COSMOSACCompound[2];
		comps[0] = db.getComp("methanol");
		comps[1] = db.getComp("methyl-acetate");		
		
		COSMOSAC cosmosac = new COSMOSAC();
		cosmosac.setComponents(comps);

		cosmosac.setTemperature(T+273.15);

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
		
		double[][] parAntoine = new double[3][3];
		parAntoine[0][0] = 16.5785;
		parAntoine[0][1] = 3638.27;
		parAntoine[0][2] = 239.500;
		parAntoine[1][0] = 14.2456;
		parAntoine[1][1] = 2662.78;
		parAntoine[1][2] = 219.690;

		double[] Psat = pSat(parAntoine, T);
		double[] P = calcPx(x1, x2, gamma1, gamma2,Psat);
		double[] y1 = calcY(x1, gamma1, Psat, P);
		
		XYPlot plot2;
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
		plot2 = (XYPlot) chart.getPlot();
		plot2.getDomainAxis().setRange(new Range(0.0, 1.0));
		plot2.getRangeAxis().setRange(new Range(15.0, 30.0));
		
		plot2.setDataset(dataset);

		XYSplineRenderer r = new XYSplineRenderer();
		BasicStroke stroke = new BasicStroke(2f);
		r.setStroke(stroke);
		plot2.setRenderer(r);
		r.setBaseShapesVisible(false);
		
		ChartPanel chartPanel = new ChartPanel(chart);
		JPanel jp = new JPanel(new BorderLayout());
		jp.add(chartPanel, BorderLayout.CENTER);
						
		add(jp);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setSize(400,500);
		
		return jp;
	}
	
/**
 * VLE diagrams using sigma profile by MOPAC
 * 
 * @return
 * @throws Exception 
 */	
	@SuppressWarnings("deprecation")
	public JPanel calcEthTolMOPAC() throws Exception{
		super.setTitle("P vs x1");
		double T = 60;
		setLayout(new BorderLayout());


		COSMOSACDataBase db = COSMOSACDataBase.getInstance();
		
		COSMOSACCompound comps[] = new COSMOSACCompound[2];
		comps[0] = db.getComp("ethanol");
		comps[1] = db.getComp("toluene");		
		
		COSMOSAC cosmosac = new COSMOPAC();
		cosmosac.setComponents(comps);

		cosmosac.setTemperature(T+273.15);

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
		
		double[][] parAntoine = new double[3][3];
		parAntoine[0][0] = 16.8958;
		parAntoine[0][1] = 3795.17;
		parAntoine[0][2] = 230.918;
		parAntoine[1][0] = 13.9320;
		parAntoine[1][1] = 3056.96;
		parAntoine[1][2] = 217.625;

		double[] Psat = pSat(parAntoine, T);
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
		};
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
		plot1.getRangeAxis().setRange(new Range(15.0, 50.0));
		
		plot1.setDataset(dataset);

		XYSplineRenderer r = new XYSplineRenderer();
		BasicStroke stroke = new BasicStroke(2f);
		r.setStroke(stroke);
		plot1.setRenderer(r);
		r.setBaseShapesVisible(false);
		
		ChartPanel chartPanel = new ChartPanel(chart);
		JPanel jp = new JPanel(new BorderLayout());
		jp.add(chartPanel);
						
		add(jp, BorderLayout.CENTER);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setSize(400,500);
		
		return jp;
	}
	
	@SuppressWarnings("deprecation")
	public JPanel calcMethanolMethilAcetateMOPAC() throws Exception{
		super.setTitle("P vs x1");
		double T = 25;
		
		COSMOSACDataBase db = COSMOSACDataBase.getInstance();

		COSMOSACCompound comps[] = new COSMOSACCompound[2];
		comps[0] = db.getComp("methanol");
		comps[1] = db.getComp("methyl-acetate");		
		
		COSMOSAC cosmosac = new COSMOPAC();
		cosmosac.setComponents(comps);

		cosmosac.setTemperature(T+273.15);

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
		
		double[][] parAntoine = new double[3][3];
		parAntoine[0][0] = 16.5785;
		parAntoine[0][1] = 3638.27;
		parAntoine[0][2] = 239.500;
		parAntoine[1][0] = 14.2456;
		parAntoine[1][1] = 2662.78;
		parAntoine[1][2] = 219.690;

		double[] Psat = pSat(parAntoine, T);
		double[] P = calcPx(x1, x2, gamma1, gamma2,Psat);
		double[] y1 = calcY(x1, gamma1, Psat, P);
		
		XYPlot plot2;
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
		plot2 = (XYPlot) chart.getPlot();
		plot2.getDomainAxis().setRange(new Range(0.0, 1.0));
		plot2.getRangeAxis().setRange(new Range(15.0, 30.0));
		
		plot2.setDataset(dataset);

		XYSplineRenderer r = new XYSplineRenderer();
		BasicStroke stroke = new BasicStroke(2f);
		r.setStroke(stroke);
		plot2.setRenderer(r);
		r.setBaseShapesVisible(false);
		
		ChartPanel chartPanel = new ChartPanel(chart);
		JPanel jp = new JPanel(new BorderLayout());
		jp.add(chartPanel, BorderLayout.CENTER);
						
		add(jp);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setSize(400,500);
		
		return jp;
	}
	
	public void biEq() throws Exception{
		JTabbedPane tabbedPane = new JTabbedPane();
		JPanel aba1 = calcEthTol();
		JPanel aba2 = calcMethanolMethilAcetate();
		JPanel aba3 = calcEthTolMOPAC();
		JPanel aba4 = calcMethanolMethilAcetateMOPAC();
		tabbedPane.addTab("Ethanol x Toluene, T=60°C",aba1);
		tabbedPane.addTab("Ethanol x Toluene, T=60°C, MOPAC",aba3);
		tabbedPane.addTab("Methanol x Methil-Acetate, T=25°C",aba2);		
		tabbedPane.addTab("Methanol x Methil-Acetate, T=25°C, MOPAC",aba4);
		add(tabbedPane, BorderLayout.CENTER);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setSize(400,500);
	}
	
	private double[] pSat(double[][] parAntoine, double T){
		double[] Psat = new double[2];
		
		Psat[0] = Antoine.CalcPsat(parAntoine[0][0], parAntoine[0][1],
									parAntoine[0][2], T);
		Psat[1] = Antoine.CalcPsat(parAntoine[1][0], parAntoine[1][1],
									parAntoine[1][2], T);
		return Psat;
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
		for (int k=0; k<21; k++){
			y1[k] = x1[k]*gamma1[k]*Psat[0]/P[k];
			}
		return y1;
	}
		
	public static void main (String[] args) throws Exception	{
		VLEdiagrams test = new VLEdiagrams();
		test.biEq();
		test.setVisible(true);
		
	}
}