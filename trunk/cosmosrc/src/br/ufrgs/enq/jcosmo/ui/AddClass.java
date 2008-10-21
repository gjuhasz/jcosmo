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

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.sql.ResultSet;
import java.sql.SQLException;

import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.ListSelectionModel;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

import br.ufrgs.enq.jcosmo.COSMOSACDataBase;

/**
 * The About dialog.
 * 
 * @author renan, rafael
 *
 */
public class AddClass extends JDialog {
	private static final long serialVersionUID = 1L;

	JTextField componente;
	JList list;
	DefaultListModel listModel;
	JButton okButton;

	private COSMOSACDialog dlg;

	public AddClass(COSMOSACDialog parent) {
		super(parent, "Add compound...", true);
		this.dlg = parent;
		setLayout(new BorderLayout());

		JPanel north = new JPanel(new GridLayout(0,2));
		JPanel south = new JPanel(new FlowLayout(FlowLayout.RIGHT));

		listModel = new DefaultListModel();
		list = new JList(listModel);
		list.setVisibleRowCount(5);
		list.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		JScrollPane listScrollPane = new JScrollPane(list);
		
		okButton = new JButton("Add");
		okButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent evt) {
				add();				
			}	
		});
		okButton.setEnabled(false);		

		JButton searchButton = new JButton("Search");
		Search Search = new Search(searchButton);
		searchButton.setActionCommand("Search");
		searchButton.addActionListener(Search);

		JButton cancelButton = new JButton("Close");
		cancelButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent evt) {
				setVisible(false);
			}
		});
		
		// connect ESC to the cancelButton
		KeyStroke esc = KeyStroke.getKeyStroke(KeyEvent.VK_ESCAPE, 0);
		getRootPane().registerKeyboardAction(new ActionListener() {
			public void actionPerformed(ActionEvent evt) {
				setVisible(false);
			}
		}, esc, JComponent.WHEN_IN_FOCUSED_WINDOW); 

		componente = new JTextField(10);
		componente.addActionListener(Search);
		componente.getDocument().addDocumentListener((DocumentListener) Search);

		north.add(componente);
		north.add(searchButton);
		south.add(okButton);
		south.add(cancelButton);

		add(north, BorderLayout.NORTH);
		add(listScrollPane, BorderLayout.CENTER);
		add(south, BorderLayout.SOUTH);

		pack();
		setLocationRelativeTo(parent);
		setVisible(true);
	}

	class Search implements ActionListener, DocumentListener {
		public Search(JButton button) {

		}
		public void actionPerformed(ActionEvent e) {
			listModel.clear();

			String name = componente.getText().toUpperCase();

			String query = "SELECT name from COSMOINDEX where name like '%" + name + "%'";

			COSMOSACDataBase db = COSMOSACDataBase.getInstance();
			ResultSet res = null;
			try {
				res = db.executeQuery(query);

				// Add all results found
				while(res!=null && res.next()){
					listModel.addElement(res.getString(1));
				}
			} catch (SQLException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}			
			if (listModel.getSize()==0){
				JOptionPane.showMessageDialog(AddClass.this, "No results.",
						"Error", JOptionPane.OK_OPTION);
				okButton.setEnabled(false);
				return;
			}

			okButton.setEnabled(true);
			list.setSelectedIndex(0);
			list.ensureIndexIsVisible(0);
		}
		//Required by DocumentListener.
		public void insertUpdate(DocumentEvent e) {			
		}
		//Required by DocumentListener.
		public void removeUpdate(DocumentEvent e) {			
		}
		//Required by DocumentListener.
		public void changedUpdate(DocumentEvent e) {			
		}
	}

	public void add(){
		String name = listModel.getElementAt(list.getSelectedIndex()).toString();
		if (name.equals("") || dlg.containsList(name)) {
			Toolkit.getDefaultToolkit().beep();
			return;
		}
		dlg.addList(name);
		dlg.visibRemove(true);
//		setVisible(false);		
	}	
}
