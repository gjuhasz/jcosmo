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

package br.ufrgs.enq.jcosmo;


import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

/**
 * The compounds' data base for COSMO-SAC.
 * This database contains the data of the sigma profiles for a large set
 * of compounds.
 * 
 * TODO: there is an open-source package named GAMESS which has a
 * COSMO code (COSGMS) but accordingly to the author it cannot compute the required
 * sigma profiles.
 * 
 * @author Rafael de Pelegrini Soares
 *
 */
public class COSMOSACDataBase
{
	public static final String DEFAULT_DATABASE = "jar:(jcosmo.jar)cosmodb";

	private static final String DRIVER = "org.apache.derby.jdbc.EmbeddedDriver";

	/** The database connection */
	protected Connection conn = null;

	/** The statement */
	protected Statement statement;

	private ResultSet res;
	
	private static COSMOSACDataBase instance = null;
	
	private String connectionURL;
	
	
	/**
	 * @return the instance of the data base {@link #DEFAULT_DATABASE}.
	 */
	public static COSMOSACDataBase getInstance(){
		return getInstance(null);
	}
	
	/**
	 * @return the instance of the data base
	 * @param database the database location, if null {@link #DEFAULT_DATABASE} is used.
	 */
	public static COSMOSACDataBase getInstance(String database){
		if(instance == null){
			// Just in case we are running this compound alone
			instance = new COSMOSACDataBase();
			try {
				instance.init(database != null ? database : DEFAULT_DATABASE);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return instance;
	}

	/**
	 * Should be called when done with the database.
	 * This function disconnects from the database, the object should not
	 * be referenced after this function is called.
	 * 
	 * @throws SQLException
	 */
	public void fini() {
		/*
		 * In embedded mode, an application should shut down Derby.
		 * Shutdown throws the 08006 exception to confirm success.
		 */
		boolean gotSQLExc = false;
		try {
			DriverManager.getConnection(connectionURL + ";shutdown=true");
		} catch (SQLException se)  {
			// NOTE for entire Derby shutdown the expected code is "XJ015"
			if ( se.getSQLState().equals("08006") ) {		
				gotSQLExc = true;
			}
		}
		if (!gotSQLExc) {
			System.out.println("Database did not shut down normally");
		}
		instance = null;
	}


	private void init(String dbName) throws ClassNotFoundException, SQLException {
		connectionURL = "jdbc:derby:" + dbName;
		
//		Properties p = System.getProperties();
//		p.put("derby.storage.tempDirectory", "/tmp");
//		p.put("derby.stream.error.file", "/tmp/derby.log");

		/*
		 **  Load the Derby driver. 
		 **     When the embedded Driver is used this action start the Derby engine.
		 **  Catch an error and suggest a CLASSPATH problem
		 */
		Class.forName(DRIVER); 

		// Create (if needed) and connect to the database
		conn = DriverManager.getConnection(connectionURL);
		conn.setReadOnly(true);
		System.out.println("Connected to database " + dbName);

		statement = conn.createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
//		statement.setFetchSize(52);
		
		instance = this;
	}

	/**
	 * Execute the given query and returns the result set.
	 * @param query the query to be executed
	 * @return the result set for the given query
	 * @throws SQLException
	 */
	public ResultSet executeQuery(String query) throws SQLException{
		if(res!=null)
			res.close();
		return res = statement.executeQuery(query);
	}

	/**
	 * Returns the compound for the given name.
	 * If the compound was already loaded the previously loaded object
	 * is returned, otherwise the database is consulted to find the
	 * compound.
	 * 
	 * @param name the name of the compound
	 * @return the compound instance for the given name
	 * 
	 * @throws SQLException
	 */
	public COSMOSACCompound getComp(String name) throws SQLException{
		COSMOSACCompound c = new COSMOSACCompound();
		
		String query = "SELECT * from COSMOINDEX where name like '" + name.toUpperCase() + "'";
		
		executeQuery(query);
		
		// The result should be exactly on row length
		if(!res.next()){
			return null;
		}
		
		c.ID = res.getInt("ID");
		c.name = res.getString("Name");
		c.formula = res.getString("Formula");
		c.CAS = res.getString("CAS");
		c.family = res.getString("Family");
		c.preOptTool = res.getString("preOptTool");
		
		c.Vcosmo = res.getDouble("Vcosmo");
		c.T = res.getDouble("T");
		c.lnPvap = res.getDouble("lnPvap");
		
		
		// loading the sigma profile
		query = "SELECT COUNT(charge) from SIGMAPROFILE where CompID = " + c.ID;
		executeQuery(query);
		if(!res.next())
			return null;
		
		int chargelength = res.getInt(1);
		c.charge = new double[chargelength];
		c.area = new double[chargelength];

		query = "SELECT charge, sigma from SIGMAPROFILE where CompID = " + c.ID;
		executeQuery(query);

//		res.last();
//		int chargelength = res.getRow();
//		c.charge = new double[chargelength];
//		c.sigma = new double[chargelength];
//		res.beforeFirst();
		
		// The result should be exactly on row length
		for (int i = 0; i < chargelength; i++) {
			res.next();
			c.charge[i] = res.getDouble(1);
			c.area[i] = res.getDouble(2);
		}
		return c;
	}
}
