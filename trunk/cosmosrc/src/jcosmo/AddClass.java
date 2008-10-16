package jcosmo;

import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.sql.ResultSet;
import java.sql.SQLException;

import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

import br.ufrgs.enq.jcosmo.COSMOSACDataBase;

/**
 * The About dialog.
 * 
 * @author rafael
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
		JPanel south = new JPanel(new GridLayout(0,2));

		listModel = new DefaultListModel();
		list = new JList(listModel);
		list.setVisibleRowCount(5);
		list.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		JScrollPane listScrollPane = new JScrollPane(list);

		JButton searchButton = new JButton("Search");
		Search Search = new Search(searchButton);
		searchButton.setActionCommand("Search");
		searchButton.addActionListener(Search);

		okButton = new JButton("Ok");
		okButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent evt) {
				ok();				
			}	
		});
		okButton.setEnabled(false);		

		JButton cancelButton = new JButton("Cancel");
		cancelButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent evt) {
				setVisible(false);
			}
		});	

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

				// The result should be exactly on row length
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

	public void ok(){
		String name = listModel.getElementAt(list.getSelectedIndex()).toString();
		if (name.equals("") || dlg.containsList(name)) {
			Toolkit.getDefaultToolkit().beep();
			return;
		}
		dlg.addList(name);
		dlg.visibRemove(true);
		setVisible(false);		
	}	
}
