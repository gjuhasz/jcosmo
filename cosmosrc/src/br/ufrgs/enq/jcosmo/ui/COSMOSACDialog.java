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

package br.ufrgs.enq.jcosmo.ui;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.sql.SQLException;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.ListSelectionModel;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.renderer.xy.XYStepRenderer;
import org.jfree.data.Range;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import br.ufrgs.enq.jcosmo.COSMOPAC;
import br.ufrgs.enq.jcosmo.COSMOSAC;
import br.ufrgs.enq.jcosmo.COSMOSACCompound;
import br.ufrgs.enq.jcosmo.COSMOSACDataBase;
import br.ufrgs.enq.jcosmo.COSMOSAC_G;
import br.ufrgs.enq.jcosmo.PCMSAC;

/**
 * Dialog for building charts of the activity coefficient using COSMO-SAC model.
 *  
 * @author Rafael de Pelegrini Soares
 */
public class COSMOSACDialog extends JFrame implements ActionListener {
	private static final long serialVersionUID = 1L;

	private static final String QUIT = "quit";
	private static final String ABOUT = "about";
	
	private JTextField temperature;
	private JTextField sigmaHB, sigmaHB2, sigmaHB3, chargeHB, beta, fpol, anorm;
	private JCheckBox ignoreSGButton;

	COSMOSACDataBase db;

	JList list;
	DefaultListModel listModel;

	JLabel gammaInf1Label;
	JLabel gammaInf2Label;
	JLabel lnGammaInf1Label;
	JLabel lnGammaInf2Label;

	XYPlot plot;
	XYPlot plotSegGamma;
	XYPlot sigmaProfilePlot;
	ChartPanel chartPanel;
	
	XYStepRenderer stepRenderer = new XYStepRenderer();

	JButton removeButton;
	
	boolean err = false;

	double cavityVolume[] = new double[2];
	double [] z = new double[2];
	double [] lnGamma = new double[2];
	
	JComboBox modelBox;

	public COSMOSACDialog() {
		super("JCOSMO Simple");
		setDefaultCloseOperation(EXIT_ON_CLOSE);
		setLayout(new BorderLayout());
		
		db = COSMOSACDataBase.getInstance();
		
		COSMOSAC models[] = new COSMOSAC[4];
		models[0] = new COSMOSAC();
		models[1] = new COSMOPAC();
		models[2] = new COSMOSAC_G();
		models[3] = new PCMSAC();
		modelBox = new JComboBox(models);
		modelBox.addActionListener(this);
		
		JPanel north = new JPanel(new GridLayout(0,2));
		add(north, BorderLayout.NORTH);
		JPanel northAba1 = new JPanel(new GridLayout(0,4));
		JPanel northAba2 = new JPanel(new GridLayout(0,2));

		//Where the GUI is created:
		JMenuBar menuBar;
		JMenu file, help;
		JMenuItem menuItem;

		//Create the menu bar.
		menuBar = new JMenuBar();

		// the file menu
		file = new JMenu("File");
		file.setMnemonic(KeyEvent.VK_F);
		menuBar.add(file);
		menuItem = new JMenuItem("Quit", KeyEvent.VK_Q);
		menuItem.setAccelerator(KeyStroke.getKeyStroke(
				KeyEvent.VK_F4, ActionEvent.ALT_MASK));
		menuItem.setActionCommand(QUIT);
		menuItem.addActionListener(this);
		file.add(menuItem);

		// the help menu
		help = new JMenu("Help");
		file.setMnemonic(KeyEvent.VK_H);
		menuBar.add(help);
		menuItem = new JMenuItem("About", KeyEvent.VK_A);
		menuItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_F1, 0));
		menuItem.setActionCommand(ABOUT);
		menuItem.addActionListener(this);
		help.add(menuItem);

		setJMenuBar(menuBar);

		listModel = new DefaultListModel();
		list = new JList(listModel);
		list.setBorder(BorderFactory.createTitledBorder("compounds"));
		list.setVisibleRowCount(2);
		list.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		JScrollPane listScrollPane = new JScrollPane(list);		

		JButton addButton = new JButton("Add");
		addButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				new AddCompoundDialog(COSMOSACDialog.this);
			}
		});

		removeButton = new JButton("Remove");
		removeButton.addActionListener(this);
		visibRemove(false);
		
		JButton calcButton = new JButton("Calculate");
		calcButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				rebuildChart();
				rebuildSigmaProfiles();
			}
		});

		JButton refreshButton = new JButton("Refresh");
		refreshButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				rebuildSigmaProfiles();
			}
		});

		ignoreSGButton = new JCheckBox("Ignore SG");
		ignoreSGButton.setToolTipText("Toogles the ignore flag for the Staverman-Guggenheim term");
		ignoreSGButton.addActionListener(this);

		JPanel but = new JPanel(new GridLayout(0,1));
		but.add(addButton, BorderLayout.EAST);
		but.add(removeButton, BorderLayout.EAST);
		but.add(modelBox);
		north.add(listScrollPane);
		north.add(but);
		
		northAba1.add(new JLabel("Temperature [K]"));
		northAba1.add(temperature = new JTextField(10));
		temperature.setText("298");

		northAba1.add(new JLabel("Sigma HB"));
		northAba1.add(sigmaHB = new JTextField(10));
		northAba1.add(new JLabel("Sigma HB2"));
		northAba1.add(sigmaHB2 = new JTextField(10));
		northAba1.add(new JLabel("Sigma HB3"));
		northAba1.add(sigmaHB3 = new JTextField(10));

		northAba1.add(new JLabel("Charge HB"));
		northAba1.add(chargeHB = new JTextField(10));
		northAba1.add(new JLabel("Beta"));
		northAba1.add(beta = new JTextField(10));
		northAba1.add(new JLabel("fpol"));
		northAba1.add(fpol = new JTextField(10));
		northAba1.add(new JLabel("Anorm"));
		northAba1.add(anorm = new JTextField(10));
		
		northAba1.add(ignoreSGButton);
		northAba1.add(calcButton);
		northAba2.add(new JLabel(""));
		northAba2.add(refreshButton);

		//		chart = new JLineChart();
		//		add(chart, BorderLayout.CENTER);
		//		chart.setTitle("Gamma Plot");
		//		chart.setSubtitle("");
		//		chart.setXAxisLabel("Mole Fraction, x_1");
		//		chart.setYAxisLabel("ln gamma, gE/RT");
		//		chart.setSource(getTitle());
		//		chart.setLegendPosition(LegendPosition.BOTTOM);
		//		chart.setShapesVisible(true);

		JFreeChart chart = ChartFactory.createXYLineChart(null, 
				"Mole Fraction, x_1", "ln gamma, gE/RT", null,
				PlotOrientation.VERTICAL, true, true, false);
		plot = (XYPlot) chart.getPlot();
		plot.getDomainAxis().setAutoRange(false);
		plot.getDomainAxis().setRange(new Range(0.0, 1.0));

		plot.setBackgroundPaint(Color.lightGray);
		plot.setDomainGridlinePaint(Color.white);
		plot.setRangeGridlinePaint(Color.white);

		XYLineAndShapeRenderer r = (XYLineAndShapeRenderer) plot.getRenderer();
		r.setUseFillPaint(true);
		r.setBaseFillPaint(Color.white);
		r.setBaseShapesVisible(true);

		JFreeChart sigmaProfileChart = ChartFactory.createXYLineChart(null, 
				"sigma", "P^x", null, PlotOrientation.VERTICAL, true, true, false);
		sigmaProfilePlot = sigmaProfileChart.getXYPlot();
		sigmaProfilePlot.getDomainAxis().setAutoRange(false);
		sigmaProfilePlot.getDomainAxis().setRange(new Range(-0.025, 0.025));

//		sigmaProfilePlot.setBackgroundPaint(Color.lightGray);
//		sigmaProfilePlot.setDomainGridlinePaint(Color.white);
//		sigmaProfilePlot.setRangeGridlinePaint(Color.white);

		JFreeChart chartSegGamma = ChartFactory.createXYLineChart(null, 
				"sigma", "Segment Gamma", null, PlotOrientation.VERTICAL, true, true, false);
		plotSegGamma = (XYPlot) chartSegGamma.getPlot();
		
		JPanel south = new JPanel();
		south.setLayout(new FlowLayout());
		south.add(new JLabel("<html>ln &gamma;<sup>&infin;</sup><sub>1</sub>:</html>"));
		south.add(lnGammaInf1Label = new JLabel());
		south.add(new JLabel("<html>ln &gamma;<sup>&infin;</sup><sub>2</sub>:</html>"));
		south.add(lnGammaInf2Label = new JLabel());
		south.add(Box.createHorizontalStrut(20));
		south.add(new JLabel("<html>&gamma;<sup>&infin;</sup><sub>1</sub>:</html>"));
		south.add(gammaInf1Label = new JLabel());
		south.add(new JLabel("<html>&gamma;<sup>&infin;</sup><sub>2</sub>:</html>"));
		south.add(gammaInf2Label = new JLabel());

		JPanel aba1 = new JPanel(new BorderLayout());
		aba1.add(northAba1, BorderLayout.NORTH);
		JPanel chartsPanel = new JPanel(new GridLayout(0,2));
		aba1.add(chartsPanel, BorderLayout.CENTER);
		chartsPanel.add(new ChartPanel(chart));
		chartsPanel.add(new ChartPanel(chartSegGamma));
		aba1.add(south, BorderLayout.SOUTH);

		JPanel aba2 = new JPanel(new BorderLayout());
		aba2.add(northAba2, BorderLayout.NORTH);
		aba2.add(chartPanel = new ChartPanel(sigmaProfileChart), BorderLayout.CENTER);

		JTabbedPane tabbedPane = new JTabbedPane();
		tabbedPane.addTab("gamma",aba1);
		tabbedPane.addTab("sigma",aba2);
		add(tabbedPane, BorderLayout.CENTER);

		//		cosmosac.setAEffPrime(6.596176570595075);
		//		cosmosac.setCoord(11.614599507917934);
		//		cosmosac.setVnorm(56.36966406129967);
		//		cosmosac.setAnorm(41.56058649432742);
		//		cosmosac.setCHB(65330.19484947528);
		//		cosmosac.setSigmaHB(0.008292411048046008);
		
		//Display the window.
		setSize(800, 600);
		setLocationRelativeTo(null);
		modelBox.setSelectedIndex(0);
		setVisible(true);

		// test for a mixture
//		addList("water");
//		addList("sec-butylamine");
//		addList("hydrogen-fluoride");
//		addList("ACETONE");
//		addList("METHANOL");
//		addList("ACETONE.opt");
//		addList("METHANOL.opt");
		addList("METHYL-ETHYL-KETONE");
		addList("ETHANOL");
		removeButton.setEnabled(true);
	}

	private void rebuildChart(){
		if (listModel.getSize()<2){
			JOptionPane.showMessageDialog(this, "Select 2 compounds.",
					"Error", JOptionPane.OK_OPTION);
			err = true;
			return;
		}
		if (listModel.getSize()>2){
			JOptionPane.showMessageDialog(this, "Select only 2 compounds.",
					"Error", JOptionPane.OK_OPTION);
			return;
		}		
		double T = Double.parseDouble(temperature.getText());
		if(T <= 0){
			JOptionPane.showMessageDialog(this, "Invalid Temperature.",
					"Error", JOptionPane.OK_OPTION);
			return;
		}

		COSMOSAC cosmosac = (COSMOSAC) modelBox.getSelectedItem();
		COSMOSACCompound comps[] = new COSMOSACCompound[2];
		try {
			comps[0] = db.getComp((String)listModel.getElementAt(0));
			comps[1] = db.getComp((String)listModel.getElementAt(1));
			cosmosac.setComponents(comps);
		} catch (Exception e1) {
			e1.printStackTrace();
			return;
		}
		if(comps[0] == null || comps[1] == null)
			return;

		cosmosac.setSigmaHB(Double.parseDouble(sigmaHB.getText()));
		cosmosac.setSigmaHB2(Double.parseDouble(sigmaHB2.getText()));
		cosmosac.setSigmaHB3(Double.parseDouble(sigmaHB3.getText()));
//		cosmosac.setSigmaHBUpper(Double.parseDouble(sigmaHBUpper.getText()));
		cosmosac.setCHB(Double.parseDouble(chargeHB.getText()));
		cosmosac.setBeta(Double.parseDouble(beta.getText()));
		cosmosac.setFpol(Double.parseDouble(fpol.getText()));
		cosmosac.setAnorm(Double.parseDouble(anorm.getText()));
		cosmosac.parametersChanged();

		setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));

//		cosmosac.setParameters(cavityVolume, c1.charge, sigma);

		cosmosac.setTemperature(T);

		// testing several compositions
		XYSeriesCollection dataset = new XYSeriesCollection();
		int n = 20;
		XYSeries lnGamma1 = new XYSeries(comps[0].name);
		XYSeries lnGamma2 = new XYSeries(comps[1].name);
		XYSeries ge_RT = new XYSeries("gE/RT");

		for(int i=0; i<=n; ++i){
			z[0] = (double)i/n; z[1] = 1-z[0];
			cosmosac.setComposition(z);
			cosmosac.activityCoefficient(lnGamma);

			lnGamma1.add(z[0], lnGamma[0]);
			lnGamma2.add(z[0], lnGamma[1]);
			ge_RT.add(z[0], z[0]*lnGamma[0] + z[1]*lnGamma[1]);

			if(z[0] == 0){
				lnGammaInf1Label.setText(String.format("%6.3g", lnGamma[0]));
				gammaInf1Label.setText(String.format("%6.3g", Math.exp(lnGamma[0])));
			}
			if(z[1] == 0){
				lnGammaInf2Label.setText(String.format("%6.3g", lnGamma[1]));
				gammaInf2Label.setText(String.format("%6.3g", Math.exp(lnGamma[1])));
			}
		}
		dataset.addSeries(lnGamma1);
		dataset.addSeries(lnGamma2);
		dataset.addSeries(ge_RT);

		plot.setDataset(dataset);

		
		// now the segment gamma
		dataset = new XYSeriesCollection();
		double [][]seggamma = cosmosac.getPureSegmentGamma(); 

		n = comps[0].charge.length;
		XYSeries g1= new XYSeries(comps[0].name);
		XYSeries g2= new XYSeries(comps[1].name);
		XYSeries g1s= new XYSeries(comps[0].name + " * sigma");
		XYSeries g2s= new XYSeries(comps[1].name + " * sigma");

		for(int j=0; j<n; ++j){
			g1.add(comps[0].charge[j], Math.log(seggamma[0][j]) );	
			g2.add(comps[1].charge[j], Math.log(seggamma[1][j]) );
			g1s.add(comps[0].charge[j], comps[0].area[j]*(Math.log(seggamma[1][j])-Math.log(seggamma[0][j])) );	
			g2s.add(comps[1].charge[j], comps[1].area[j]*(Math.log(seggamma[0][j])-Math.log(seggamma[1][j])) );
		}
		dataset.addSeries(g1);
		dataset.addSeries(g2);
		dataset.addSeries(g1s);
		dataset.addSeries(g2s);
		plotSegGamma.setDataset(dataset);
		
		// adjust the plot properties
		plotSegGamma.getDomainAxis().setAutoRange(false);
		plotSegGamma.getDomainAxis().setRange(new Range(-0.025, 0.025));
		XYLineAndShapeRenderer r = (XYLineAndShapeRenderer) plotSegGamma.getRenderer();
		r.setSeriesStroke(0, new BasicStroke(2.5f));
		r.setSeriesStroke(1, new BasicStroke(2.5f));
		BasicStroke dashed = new BasicStroke(
				2.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND,
				1.0f, new float[] {6.0f, 6.0f}, 0.0f
		);

		r.setSeriesStroke(2, dashed);
		r.setSeriesStroke(3, dashed);
		r.setSeriesPaint(0, Color.RED);
		r.setSeriesPaint(1, Color.BLUE);
		r.setSeriesPaint(2, Color.RED);
		r.setSeriesPaint(3, Color.BLUE);
//		plotSegGamma.setRenderer(stepRenderer);
//		plotSegGamma.setRenderer(3, stepRenderer);
		
		setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
	}

	private void rebuildSigmaProfiles(){
		if (listModel.getSize()==0){
			if (err == true){ 
				err = false;
				return;
			}
			JOptionPane.showMessageDialog(this, "Select compounds.",
					"Error", JOptionPane.OK_OPTION);
			return;
		}

		int num = listModel.getSize();
		COSMOSACCompound[] c = new COSMOSACCompound[num];
		double [][] area = new double[num][];
		double [][] charge = new double[num][];
		XYSeriesCollection dataset = new XYSeriesCollection();

		for (int i=0; i<num; i++){
			try {
				c[i] = db.getComp((String)listModel.getElementAt(i));			
			} catch (SQLException e1) {
				e1.printStackTrace();
				return;
			}
		}
		
		COSMOSAC cosmosac = (COSMOSAC) modelBox.getSelectedItem();
		try {
			cosmosac.setComponents(c);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		for (int i=0; i<num; i++){

			if(c[i] == null)
				return;

			setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));

			area[i] = c[i].area;			
			charge[i] = c[i].charge;

			int n = charge[0].length;
			XYSeries comp= new XYSeries(c[i].name);

			// charges represent the center of the segments
			comp.add(charge[i][0], area[i][0]);
			for(int j=1; j<n; ++j){
				comp.add(charge[i][j]-(charge[i][j]-charge[i][j-1])/2, area[i][j]);
			}

			dataset.addSeries(comp);
			sigmaProfilePlot.setRenderer(i, stepRenderer);
			sigmaProfilePlot.getRenderer().setSeriesStroke(i, new BasicStroke(2.5f));
		}
		sigmaProfilePlot.setDataset(dataset);

		setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
	}

	public void actionPerformed(ActionEvent e) {
		COSMOSAC cosmosac = (COSMOSAC) modelBox.getSelectedItem();
		if(modelBox == e.getSource()){
			
			chargeHB.setText(Double.toString(cosmosac.getCHB()));
			ignoreSGButton.setSelected(cosmosac.isIgnoreSG());
			sigmaHB.setText(Double.toString(cosmosac.getSigmaHB()));
			sigmaHB2.setText(Double.toString(cosmosac.getSigmaHB2()));
			sigmaHB3.setText(Double.toString(cosmosac.getSigmaHB3()));
			beta.setText(Double.toString(cosmosac.getBeta()));
			fpol.setText(Double.toString(cosmosac.getFpol()));
			anorm.setText(Double.toString(cosmosac.getAnorm()));
		}
		else if(ignoreSGButton == e.getSource()){
			cosmosac.setIgnoreSG(ignoreSGButton.isSelected());
		}
		else if(removeButton == e.getSource()){
			//This method can be called only if
			//there's a valid selection
			//so go ahead and remove whatever's selected.
			int index = list.getSelectedIndex();
			listModel.remove(index);

			int size = listModel.getSize();

			if (size == 0) { //Nobody's left, disable firing.
				visibRemove(false);

			} else { //Select an index.
				if (index == listModel.getSize()) {
					//removed item in last position
					index--;
				}

				list.setSelectedIndex(index);
				list.ensureIndexIsVisible(index);
			}
		}
		else if(e.getActionCommand() == QUIT){
			setVisible(false);
			db.fini();
			System.exit(0);
		}
		if(e.getActionCommand() == ABOUT){
			new AboutDialog(this);
		}
	}
	public void addList(String name) {
		listModel.addElement(name);
	}
	public void visibRemove(Boolean b) {
		removeButton.setEnabled(b);
		list.setSelectedIndex(0);
		list.ensureIndexIsVisible(0);
	}
	public boolean containsList(String name) {
		return listModel.contains(name);
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// create and show the GUI using a thread
		javax.swing.SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				//				try {
				//					UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
				//				} catch (Exception e) {
				//					e.printStackTrace();
				//				}
				new COSMOSACDialog();
			}
		});
	}
}