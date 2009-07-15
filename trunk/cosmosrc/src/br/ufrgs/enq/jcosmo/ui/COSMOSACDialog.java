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

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.renderer.xy.XYSplineRenderer;
import org.jfree.data.Range;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import br.ufrgs.enq.jcosmo.COSMOSAC;
import br.ufrgs.enq.jcosmo.COSMOSACCompound;
import br.ufrgs.enq.jcosmo.COSMOSACDataBase;

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

	COSMOSACDataBase db;

	JList list;
	DefaultListModel listModel;

	JLabel gammaInf1Label;
	JLabel gammaInf2Label;
	JLabel lnGammaInf1Label;
	JLabel lnGammaInf2Label;

	XYPlot plot;
	XYPlot plot2;
	ChartPanel chartPanel;

	JButton removeButton;
	
	boolean err = false;

	double cavityVolume[] = new double[2];
	double [] z = new double[2];
	double [] lnGamma = new double[2];
	COSMOSAC cosmosac = new COSMOSAC();

	@SuppressWarnings("deprecation")
	public COSMOSACDialog() {
		super("JCOSMO Simple");
		setDefaultCloseOperation(EXIT_ON_CLOSE);
		setLayout(new BorderLayout());

		db = COSMOSACDataBase.getInstance();

		JPanel north = new JPanel(new GridLayout(0,2));
		add(north, BorderLayout.NORTH);
		JPanel northAba1 = new JPanel(new GridLayout(0,2));
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
		JScrollPane listScrollPane = new JScrollPane(list);		

		JButton addButton = new JButton("Add");
		addButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				new AddClass(COSMOSACDialog.this);
			}
		});

		removeButton = new JButton("Remove");
		removeButton.addActionListener((ActionListener) new Remove());
		visibRemove(false);

		JButton calcButton = new JButton("Calculate");
		calcButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				rebuildChart();
				rebuildChart2();
			}
		});

		JButton refreshButton = new JButton("Refresh");
		refreshButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				rebuildChart2();
			}
		});

		JCheckBox ignoreSGButton = new JCheckBox("Ignore SG");
		ignoreSGButton.setToolTipText("Toogles the ignore flag for the Staverman-Guggenheim term");
		ignoreSGButton.setSelected(cosmosac.isIgnoreSG());
		ignoreSGButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				cosmosac.setIgnoreSG(((JCheckBox)e.getSource()).isSelected());
			}
		});

		JPanel but = new JPanel(new GridLayout(0,1));
		but.add(addButton, BorderLayout.EAST);
		but.add(removeButton, BorderLayout.EAST);
		north.add(listScrollPane);
		north.add(but);
		
		northAba1.add(new JLabel("Temperature [K]"));
		northAba1.add(temperature = new JTextField(10));
		temperature.setText("298");
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

		JFreeChart chart2 = ChartFactory.createXYLineChart(null, 
				"sigma", "P^x", null, PlotOrientation.VERTICAL, true, true, false);
		plot2 = (XYPlot) chart2.getPlot();
		plot2.getDomainAxis().setAutoRange(false);
		plot2.getDomainAxis().setRange(new Range(-0.025, 0.025));

		plot2.setBackgroundPaint(Color.lightGray);
		plot2.setDomainGridlinePaint(Color.white);
		plot2.setRangeGridlinePaint(Color.white);

		XYSplineRenderer r2 = new XYSplineRenderer();
		BasicStroke stroke = new BasicStroke(2.5f);
		r2.setStroke(stroke);

		plot2.setRenderer(r2);
//		XYLineAndShapeRenderer r2 = (XYLineAndShapeRenderer) plot2.getRenderer();
//		r2.setUseFillPaint(true);
//		r2.setBaseFillPaint(Color.white);
		r2.setBaseShapesVisible(false);

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
		aba1.add(chartPanel = new ChartPanel(chart), BorderLayout.CENTER);
		aba1.add(south, BorderLayout.SOUTH);

		JPanel aba2 = new JPanel(new BorderLayout());
		aba2.add(northAba2, BorderLayout.NORTH);
		aba2.add(chartPanel = new ChartPanel(chart2), BorderLayout.CENTER);

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
		setSize(500, 600);
		setLocationRelativeTo(null);
		setVisible(true);
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

		COSMOSACCompound c1, c2;
		try {
			c1 = db.getComp((String)listModel.getElementAt(0));
			c2 = db.getComp((String)listModel.getElementAt(1));
		} catch (SQLException e1) {
			e1.printStackTrace();
			return;
		}

		if(c1 == null || c2 == null)
			return;

		setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));

		cavityVolume[0] = c1.Vcosmo;
		cavityVolume[1] = c2.Vcosmo;

		double [][] sigma = new double[2][];
		sigma[0] = c1.sigma;
		sigma[1] = c2.sigma;

		cosmosac.setParameters(cavityVolume, c1.charge, sigma);

		cosmosac.setTemperature(T);

		// testing several compositions
		XYSeriesCollection dataset = new XYSeriesCollection();
		int n = 20;
		XYSeries lnGamma1 = new XYSeries(c1.name);
		XYSeries lnGamma2 = new XYSeries(c2.name);
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

		setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
	}

	private void rebuildChart2(){
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
		double [][] sigma = new double[num][];
		double [][] charge = new double[num][];
		XYSeriesCollection dataset = new XYSeriesCollection();

		for (int i=0; i<num; i++){
			try {
				c[i] = db.getComp((String)listModel.getElementAt(i));			
			} catch (SQLException e1) {
				e1.printStackTrace();
				return;
			}

			if(c[i] == null)
				return;

			setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));

			sigma[i] = c[i].sigma;			
			charge[i] = c[i].charge;

			int n = charge[0].length;
			XYSeries comp= new XYSeries(c[i].name);

			for(int j=0; j<n; ++j){
				comp.add(charge[i][j], sigma[i][j]);	
			}

			dataset.addSeries(comp);
			plot2.setDataset(dataset);
		}

		setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
	}

	class Remove implements ActionListener {
		public void actionPerformed(ActionEvent e) {
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
	}

	public void actionPerformed(ActionEvent e) {
		if(e.getActionCommand() == QUIT){
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