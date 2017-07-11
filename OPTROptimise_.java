/**
 * Plugin for fine reconstruction of 3D trajectory from MEVISLAB's Tubular tracking
 * Lionel Dupuy 2016
 * JHI - Plant Systems Modelling
 */
 
import ij.*;
import ij.measure.*;
import ij.plugin.*;
import ij.plugin.filter.*;
import ij.process.*;
import ij.process.ImageProcessor.*;
import ij.process.AutoThresholder.Method;
import ij.gui.*;
import ij.gui.PolygonRoi;

import ij.plugin.filter.GaussianBlur;
import ij.plugin.ContrastEnhancer;

import java.awt.*;
//import java.lang.Math.*;
import ij.plugin.frame.RoiManager;
import ij.plugin.filter.MaximumFinder;
import ij.process.ImageConverter;
import ij.process.EllipseFitter;
import ij.plugin.ImageCalculator;

import ij.gui.Roi;
import ij.gui.Roi.*;
import ij.gui.ShapeRoi;
import ij.plugin.frame.RoiManager;
import java.util.*;
import java.awt.font.*;
import java.io.*;

import java.awt.Color;

import java.util.Arrays;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.PowellOptimizer;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.FileUtils;
// Open dialog box
import javax.swing.JFileChooser;

//////////////////////////////////////////////////////////////////////////////
// Libraries for reading xml
import javax.xml.parsers.DocumentBuilder; 
import javax.xml.parsers.DocumentBuilderFactory;
//These classes are for the exceptions that can be thrown when the XML document is parsed:
import org.xml.sax.ErrorHandler;
import org.xml.sax.SAXException; 
import org.xml.sax.SAXParseException;
import org.xml.sax.helpers.*;

//These classes read the sample XML file and manage output:
import java.io.File;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;

//Finally, import the W3C definitions for a DOM, DOM exceptions, entities and nodes:
import org.w3c.dom.Document;
import org.w3c.dom.DocumentType;
import org.w3c.dom.Entity;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;


import Jama.Matrix; 

public class OPTROptimise_ implements PlugIn {
	@SuppressWarnings("unchecked") 
	/********************************************************************************************
	 * Parameters of the profile algorithm
	 * Change here for best fit to your image
	 */
	int w_profile = 200; 				// fwidth of the root profile in pixels
	int margin = 10;
	int l_profile = w_profile/4;		// length of root segments in pixels (after resampling)	
	int sign = 1;						// direction of the rotation
	
	int n = 100;							// number of projections used for calibration in Find_translation
	int m = 20;							// number of position along the root for calibration in Find_translation
	int dn = 40;						// step for resampling the initial trajectory
	double dl = 0.25;					// size of step in mm
	double totL = 0;					// total length of the root in mm
	double scale = 1.0/71.3;//66.8;				// size of one pixel in mm
	int ds = 20;						// slice increment for fitting
	
	/**
	 *********************************************************************************************
	 */
	
	// constants
	int n_slice;
	ImagePlus imp;
	File currentFile = null;
	// raw data from the tracing
	double[] x; // x coordinates of the tracing 
	double[] y; // y coordinate of the tracing
	double[] z; // z coordinate of the tracing
	int [][] CXik; //set of centerline point extracted projection data using local maxima. 
	double[] CXk; //current set of centerline point extracted projection data using local maxima.  a given point on the trajectory, n is number of slices
	
	double[] XXi;		// the corrected X coordinate of the trajectory
	double[] YYi;		// the corrected Y coordinate of the trajectory
	double[] ZZi;		// the corrected Z coordinate of the trajectory
	double[] Di;		// diameters of the root
	double[] A;		// angle distribution of the projections
	double[] ax_x;// centre of rotation of the object on the image
	int ex;		// translation in x due to the crop of the image
	int ey; 	// translation in y due to the crop of the image
	// Image data
	int[][][] IM;
	
/*	float[] X; // X coordinates of the tracing 
	float[] Y; // Y coordinate of the tracing
	float[] Z; // Z coordinate of the tracing
*/	
	float[] L; // distance from the tip (curvilinear abscisa)
	float[] NX; // normal vector coord x
	float[] NY; // normal vector coord y
	
	// float image
	float[][] I;
	float[][] PROFILES;

	
	ImageStack main_stack;
	int width = 0;
	int height = 0;
	public void run(String arg) 
	{
		if (IJ.versionLessThan("1.26i"))
			return;
		//IJ.run("32-bit");
		imp = IJ.getImage();
		
		// Image Data
		ImageProcessor IP = imp.getProcessor();
		width = IP.getWidth();
		height = IP.getHeight();

		
		Read_Profile_From_Text();
		calc_scale();
		Find_Translation();
		Refine_Profile();
		Export_Data();
		//Refresh_Roi();
	}


	/******************************************************************************************
	 * Read ROI data
	 */	
	private void Read_Profile_From_Text()
	{
		try
		{
			// settings
			//String filename = null;
			currentFile = null;
			boolean dtdValidate = false;
			boolean xsdValidate = false;
			String schemaSource = null;
		
			// Open File Dialog
			final JFileChooser fc = new JFileChooser();
			
			fc.setCurrentDirectory(new File("D://LIONEL//PROGRAMMING//CODE//CORE//IMAGEJ//OPT_OPTIMISE_TRACE//" ));// "E://Lionel//PROGRAMMING//IMAGEJ//OPT_OPTIMISE_TRACE//"
				//In response to a button click:
			int result = fc.showOpenDialog(null);
			if (result == JFileChooser.APPROVE_OPTION) {
				currentFile = fc.getSelectedFile();
				}
			
			// create xml objects
			DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
			dbf.setNamespaceAware(true);
			dbf.setValidating(dtdValidate || xsdValidate);
	
			DocumentBuilder db = dbf.newDocumentBuilder();
			Document doc = db.parse(currentFile);

			// Extract coordinates of the tracing
			NodeList itemlist = doc.getElementsByTagName("pos");
			x = new double[itemlist.getLength()];
			y = new double[itemlist.getLength()];
			z = new double[itemlist.getLength()];
			
			for (int temp = 0; temp < itemlist.getLength(); temp++) {
				Node nNode = itemlist.item(temp);
				if (nNode.getNodeType() == Node.ELEMENT_NODE) {
					String[] R = nNode.getTextContent().split(" ");
					x[temp] = (double)(Float.parseFloat(R[0]));
					y[temp] = (double)(Float.parseFloat(R[1]));
					z[temp] = (double)(Float.parseFloat(R[2]));
					//IJ.log("POS: x:" + x[temp] + " | y: " + y[temp] + " | z: "+ z[temp] );
				}
			}
			
		}
		catch (Exception e) {
			e.printStackTrace();
    	}
	}
	/******************************************************************************************
	//***********************************************************************************
	/ Determine optimal translation using pseudo inverse method (from least square)
	/ the result is ex ey, the translation between the center of the "traced image" and the center of the raw data
	/ the difference exist because the data is being croped and also because the reconstruction is an image 
	/ which dimension are a square inscribed in the circle induced by the rotation
	//*************************************************************************************
	 */
	private void calc_scale()
	{
		// determine total root length
		for (int i=1;i<x.length;i++)
		{
			totL = totL + scale*Math.sqrt( (x[i] - x[i-1])*(x[i] - x[i-1])  +  (y[i] - y[i-1])*(y[i] - y[i-1])  +  (z[i] - z[i-1])*(z[i] - z[i-1]));	
		}
		IJ.log("Root Length = " + totL);
		
		// determine dn
		dn = (int)(x.length / (totL / dl));
		IJ.log("ntot = " + x.length / dn);
	}	 
	private void Find_Translation()
	{
		IJ.log("step1");
		ImageStack STK = imp.getImageStack();


		int[] krota = new int[n];
		int[] x_img = new int[n];
		int[] x_trace = new int[n];
		double [][] CS = new double[n*m][2];
		double [][] Y = new double[n*m][2];

		for (int i=0;i<n;i++)
			{krota[i] = (int)(i*720./n+1.+0.5);}
		

		for (int i=0;i<n;i++)
			{
			imp.setSlice(krota[i]);
			ImageProcessor IP1 = STK.getProcessor(krota[i]);
			I = new float[width][height];
			I = IP1.getFloatArray();
				
			for (int j=0;j<m;j++)
				{
				int k = j*n+i;
				double X0 = x[(int)((x.length-1)*j/(m-1) )];
				double Y0 = y[(int)((x.length-1)*j/(m-1) )];
				double Z0 = z[(int)((x.length-1)*j/(m-1) )];
				double pi = 3.1519;
				double angle = sign*(2*pi*(krota[i] - 1)/360.)/2.;
				double c = Math.cos(angle);
				double s = Math.sin(angle);
				CS[k][0] = c;
				CS[k][1] = -s;
				
				int depth = (int)(Z0+0.5);
				
				if (depth < height & depth >=0)
				{
					int jmax1 = find_max(get_global_hprofile(margin, depth, I)) + margin; //find_max(get_local_hprofile(Xi[i], Yi[i])) + Xi[i] - (int)(w_profile/2+0.5);//
					
					PointRoi Proi1 = new PointRoi(jmax1,depth);
					imp.setRoi((Roi)Proi1);
	
					x_img[i] = jmax1;
	
					Y[k][0] = (double)jmax1 - c*(X0 - width/2.) + s*(Y0 - width/2.) - width/2.;
				}
				}
			
		}	
		
		Matrix YY = new Matrix(Y);
		Matrix CSCS = new Matrix(CS);
		Matrix res2 = CSCS.inverse().times(YY);
		Matrix res1 = new Matrix(new double[][] {{64},{64}});


		ex = (int)(res2.get(0,0)+0.5); // 64;//
		ey = (int)(res2.get(1,0)+0.5); // 64;//
		
		int n_slices = imp.getStackSize();
		

/*		for (int k=0;k<n_slices/ds;k++)
		{
			double pi = 3.1519;
			double angle = sign*(2.*pi* (double)k*ds)/((double)n_slices);
			int[] Xi = new int[(int)(x.length/dn)];
			int[] Yi = new int[(int)(x.length/dn)];
			
		}*/
		IJ.log("step2");
	}

	/******************************************************************************************
	 * Resample ROI data for refining the path
	 * Adjust the position of the root trace so that it is in the centre of the root
	 * adjsustment is similar to how ex and ey are determined except the adjustment is made on each point of the trajectory
	 */

	private void Refine_Profile()
	{
		IJ.log("step3");
		imp.deleteRoi();
		int n_slices = imp.getStackSize();
		ImageStack STK = imp.getImageStack();

		
		XXi = new double[(int)x.length/dn];		// the corrected X coordinate of the trajectory
		YYi = new double[(int)y.length/dn];		// the corrected Y coordinate of the trajectory
		ZZi = new double[(int)y.length/dn];		// the corrected Z coordinate of the trajectory
		Di =  new double[(int)y.length/dn];		// diamter of the root
		CXik = new int[(int)x.length/dn][n_slices];	// the centerline on each project obtained from the image
		CXk = new double[n_slices];							// same but only at one position along the trajectory
		A = new double[n_slices];							// angle distribution of the projections
		IM = new int[n_slices][XXi.length][w_profile];
		double [][] jmean = new double[(int)x.length/dn][1];
		double [][] Coeff = new double[(int)x.length/dn][2];	
		//double [][] Y0 = new double[(int)x.length/dn][1];
		// needed is more sophisticated fitting needed
		PROFILES = new float[x.length/dn][w_profile];
		
		//*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *

		
		// determines the position of the center line CXik[][]
		// determines the rotation  AB[], coefficient of the line defining the axis

		for (int i=0;i<XXi.length;i++) { Di[i] = 0;}	
			
		for (int k=0;k<n_slices;k++)
		{
			double pi = 3.1519;
			double angle = sign*(2.*pi* (double)k)/((double)n_slices);
			int[] Xi = new int[(int)x.length/dn];
			int[] Zi = new int[(int)z.length/dn];
			ImageProcessor IP1;
			IP1 = STK.getProcessor(k+1);
			if (k%(720/m) == 0)
			{
							
				I = new float[width][height];
				I = IP1.getFloatArray();
			}
			for (int i=0;i<XXi.length;i++)
			{
				if (k==0) { jmean[i][0] = 0; }
				double c = Math.cos(angle);
				double s = Math.sin(angle); 
				double xi = (x[i*dn] + ex - ((double)width)/2.);
				double yi = (y[i*dn] + ey - ((double)width)/2.);
				double xtmp = c*xi - s*yi;
				double ytmp = s*xi + c*yi;
				XXi[i] = (double)(x[i*dn] + ex);
				YYi[i] = (double)(y[i*dn] + ey);
				ZZi[i] = (double)(z[i*dn]);
				Xi[i] = (int)(xtmp + width/2. +0.5);
				Zi[i] = (int)(z[i*dn]+0.5);	

				// calculate diameter
				if (k%(720/m) == 0)
				{
					float[] N = {0,0};
					if (i==0)
					{
						N[0] = - (Zi[1] - Zi[0]);
						N[1] =   (Xi[1] - Xi[0]);
	
						
					}
					else
					{
						N[0] = - (Zi[i] - Zi[i-1]);
						N[1] =   (Xi[i] - Xi[i-1]);		
						//	
						double[] profile = get_local_Nprofile(Xi[i], Zi[i], N);	// uses I to get the profile	
						Di[i] = Di[i] + ((double)get_diameter_from_profile(profile))/((float)m);	
					}
				}
				 
				
				// calcualte coefficient for optimal position of the point
				//IM[k][i] = get_local_hprofile(I, Xi[i], Zi[i]);
				IP1.getRow(Xi[i] - (int)(w_profile/2.+0.5), Zi[i], IM[k][i], w_profile);
				int xmax = find_max(IM[k][i]) + Xi[i] - (int)(w_profile/2.+0.5);	//find_max(get_global_hprofile(margin, Zi[i], I)) + margin;//			

				CXik[i][k] = xmax ;
				jmean[i][0] += xmax/((double)n_slices);
				if (k == n_slices-1)
				{
					Coeff[i][0] = 1;
					Coeff[i][1] = z[i*dn];
				}
			}
		}

			
	IJ.log("step4");
	// fit rotation axis
	Matrix YY0 = new Matrix(jmean);
	Matrix C = new Matrix(Coeff);
	Matrix AB = C.inverse().times(YY0);
//	double a = AB.get(0,0);
//	double b = AB.get(1,0);
	ax_x = new double[2];
	ax_x[0] = AB.get(0,0);
	ax_x[1] = AB.get(1,0);
	IJ.log("step5");
	//*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
	// Global fit (optimisation) of the path
	for (int i=0;i<XXi.length;i++)
		{
		double [][] CS = new double[n_slices][2];
		double [][] Y = new double[n_slices][1];
		
		double axis_x = ax_x[0]+ax_x[1]*z[i*dn]; 
//		axis_x = width/2.;					// transformation brings no improvements
		for (int k=0;k<n_slices;k++)
			{
			double pi = 3.1519;
			double angle = sign*(2.*pi* (double)k)/((double)n_slices);
			double c = Math.cos(angle);
			double s = Math.sin(angle);
			CS[k][0] = c;
			CS[k][1] = -s;				
			Y[k][0] = CXik[i][k] + c*(axis_x) - s*(axis_x) - axis_x;//c*(width/2./Math.sqrt(2)) - s*(width/2./Math.sqrt(2)) - width/2.;
			
			// data for non linear optimisation
			CXk[k] = (double)(CXik[i][k]);
			A[k] = angle;
			}
			
		Matrix YY = new Matrix(Y);
		Matrix CSCS = new Matrix(CS);
		Matrix res = CSCS.inverse().times(YY);
		XXi[i] = res.get(0,0);
		YYi[i] = res.get(1,0);
			
		}
	IJ.log("step6");
	// quadratic optimisation
	// use previous results to initiate the non linear optimisation
	//double[] init_param = {ax_x[0], 0., 0., 0.,0.};
	double[] init_param = new double[XXi.length*2];	
	for (int i=0;i<XXi.length;i++)
	{
		init_param[i] = XXi[i]+20.;
		init_param[XXi.length+i] = YYi[i];
	}

	// create optimiser objects
	SimplexOptimizer optim1 = new SimplexOptimizer(1e-3, 1e-5);
	PowellOptimizer optim2 = new PowellOptimizer(1e-3, 1e-5);
	
	ObjFunction OFunc = new ObjFunction();
	
	PointValuePair optimum =
	optim1.optimize(
               new MaxEval(10000), 
               new ObjectiveFunction(OFunc), 
               GoalType.MINIMIZE, 
               new InitialGuess(init_param), //);//  
               new NelderMeadSimplex(2*XXi.length));

	
	
//    IJ.log(""+Arrays.toString(optimum.getPoint()));	
	double[] snake = optimum.getPoint();				// intrinsic parameters
	// (displacement only along the x axis)
	
	// distribute nodes of the tracine evenly
	IJ.log("step7");
	// Fine fit of the path 
	
	// Plot the reconstructed data
	int ds = 5;
	for (int k=0;k<n_slices/ds;k++)
		{
		double pi = 3.1519;
		double angle = sign*(2.*pi* (double)k*ds)/((double)n_slices);
		
		int[] Xi0 = new int[(int)x.length/dn];
		int[] Xi = new int[(int)x.length/dn];
		int[] Zi = new int[(int)y.length/dn];			// optimal position calculated from local maxima
		int[] PXi = new int[(int)x.length/dn];			// local maxima
		int[] SXi = new int[(int)x.length/dn];			// centerline adjsuted by snake computation
		for (int i=0;i<XXi.length;i++)			
			{
			double axis_x = ax_x[0]+ax_x[1]*z[i*dn];
			//axis_x = width/2.;					// transformation brings no improvements
			double c = Math.cos(angle);
			double s = Math.sin(angle); 
			double xi = XXi[i]- axis_x;//((double)width)/2./Math.sqrt(2);
			double yi = YYi[i]- axis_x;
			double xtmp = c*xi - s*yi;
			double ytmp = s*xi + c*yi;
			Xi[i] = (int)(xtmp + axis_x +0.5);
			Zi[i] = (int)(z[i*dn]+0.5);		
			
			
			double xi0 = x[i*dn]- width/2.;
			double yi0 = y[i*dn]- width/2.;	
			xtmp = c*xi0 - s*yi0;
			ytmp = s*xi0 + c*yi0;
			Xi0[i] = (int)(xtmp + width/2. +0.5);	
			
			// local maxima
			PXi[i] = CXik[i][k*ds];
			
			// snake coordinates
			xi = snake[i]- axis_x;//XXi[i] + snake[i]*0 				- w_profile/2. - axis_x;//((double)width)/2./Math.sqrt(2);
			yi = snake[XXi.length + i]- axis_x;//YYi[i] + snake[XXi.length + i]*0   - w_profile/2. - axis_x;			
		
			double xsnake = c*xi - s*yi;			
			SXi[i] = (int)(xsnake + axis_x +0.5);
			}
		imp.setSlice(k*ds+1);
		RoiManager manager = RoiManager.getInstance();
		if (manager == null)
			manager = new RoiManager();
		
/*		PolygonRoi Proi0 = new PolygonRoi(Xi0,Zi,Xi0.length,Roi.POLYLINE);
		Proi0.setStrokeColor(Color.red);
		imp.setRoi((Roi)Proi0);
		manager.addRoi(Proi0); 
			*/
				
/*		PolygonRoi Proi2 = new PolygonRoi(PXi,Zi,PXi.length,Roi.POLYLINE);
		Proi2.setStrokeColor(Color.blue);
		//Proi.setSize(large);
		imp.setRoi((Roi)Proi2);
		manager.addRoi(Proi2);
	*/	
		PolygonRoi Proi = new PolygonRoi(Xi,Zi,Xi.length,Roi.POLYLINE);
		//Proi.setStrokeColor(Color.red);
		//imp.setRoi((Roi)Proi);
		//manager.addRoi(Proi);
		
		
		PolygonRoi Proi3 = new PolygonRoi(SXi,Zi,PXi.length,Roi.POLYLINE);
		Proi3.setStrokeColor(Color.cyan);
		Proi.setLineWidth(1);
		imp.setRoi((Roi)Proi3);
		manager.addRoi(Proi3);
		
		}
		IJ.log("step8");
	}
	/******************************************************************************************
	 * Extract profile along a vertical line / exluding a margin
	 */	
	private double[] get_global_hprofile(int window, int depth, float[][] I)
	{

		double[] profile = new double[width-2*window];
		for (int j=window;j<width-window;j++)
		{
			double m = 0;
			for (int k=0;k<window;k++)
				{
				//IJ.log("PUTE: " + depth + " , " + height);
				m = m + (double)I[j-(int)(window/2.+0.5)+k][depth] / (double)window;
				//IJ.log("OK");
				}
			
			profile[j-window] = m;

		}
		int jmax = find_max(profile);

		return profile;
	}
	/******************************************************************************************
	 * Extract profile along the x axis perpendicular to the root
	 */
	private double[] get_local_hprofile(float[][] IMG, int xx, int yy)
		{
		double[] profile = new double[w_profile];
		//IJ.log("TEST"  + w_profile);
		for (int j=0;j<w_profile;j++)
			{
			int pos = xx+j-(int)(w_profile/2+0.5);
			if (pos>=0 && pos<width)
				{
					profile[j] = IMG[pos][yy];}
			else
				{profile[j] = 0;}
			}
		return profile;
		}	
	/******************************************************************************************
	 * Determine diameter from profile
	 */
	private int get_diameter_from_profile(double[] profile)
		{
		// determine the gradient
		double[] gradient = new double[profile.length-1];
		for (int j=0;j<profile.length-1;j++)
			{
			gradient[j] = Math.abs(profile[j+1] - profile[j]);
			}

		// find local maxima in the gradient	
		double[] local_max = new double[gradient.length];
		for (int j=1;j<gradient.length-1;j++)
			{
			if ((gradient[j+1] - gradient[j]) * (gradient[j] - gradient[j-1]) < 0 &  (gradient[j] - gradient[j-1]) > 0)
				{local_max[j] = gradient[j];}
			else
			{
				local_max[j] = 0;
			}
			}	

		// determine diameter
		int d1 = (int)(gradient.length/2. + 0.5);
		int d2 = (int)(gradient.length/2. + 0.5);
		for (int j=1;j<gradient.length;j++)
			{
			if (local_max[j] > 0 & j<= (int)(gradient.length/2. + 0.5))
				{
				d1 = j;
				}
			else if (local_max[j] > 0 & j>= (int)(gradient.length/2. + 0.5))
				{
				d2 = j;
				}
			}

		return d2 - d1;
		}	

		
	/******************************************************************************************
	 * Extract profile along the axis perpendicular to the root
	 */
	private double[] get_local_Nprofile(float xx, float yy, float[] N)
		{
		double[] profile = new double[w_profile];
		for (int j=0;j<w_profile;j++)
		{
			int xj = (int)xx+(int)((j-w_profile/2.)*N[0]+0.5);
			int yj = (int)yy+(int)((j-w_profile/2.)*N[1] +0.5);
			if (xj>=0 & yj >=0 & xj<I.length & yj < I[0].length)
				{
					profile[j] = I[xj][yj];
				}
			else
				{
					profile[j] = 0;
				}
		}
		return profile;
		}	
				
	/******************************************************************************************
	 * find_maxima
	 */		
	private int find_max(double[] profile)
	{
		int jmax = 0;
		double Imax = 0;
		for (int i=0;i<profile.length;i++)
		{
			if (profile[i] > Imax)
			{
				jmax = i;
				Imax = profile[i];
			}
		}
		return jmax;
	}
	private int find_max(int[] profile)
	{
		int jmax = 0;
		double Imax = 0;
		for (int i=0;i<profile.length;i++)
		{
			if (profile[i] > Imax)
			{
				jmax = i;
				Imax = profile[i];
			}
		}
		return jmax;
	}
	/******************************************************************************************
	 * Export Data
	 */
	private void Refresh_Roi()
	{
		imp.deleteRoi();
		int[] XXi = new int[x.length];
		int[] YYi = new int[y.length];
		for (int i=0;i<XXi.length;i++)
		{
			XXi[i] = (int)(x[i]+0.5);
			YYi[i] = (int)(y[i]+0.5);
		}
		PolygonRoi Proi = new PolygonRoi(XXi,YYi,XXi.length,Roi.POLYLINE);
		imp.setRoi((Roi)Proi);
	}
	
	/******************************************************************************************
	 * Export Data
	 */
	private void Export_Data()
	{

	try {
		String BaseFile = FilenameUtils.removeExtension(currentFile.getAbsolutePath());
		String CompFile = BaseFile + "_snake_scaled.txt";
		
		FileWriter outFile = new FileWriter(CompFile);
		PrintWriter out = new PrintWriter(outFile);

		// Write text to file
		for (int i=0;i<XXi.length;i++)
			{	
			//out.println("" + XXi[i] + "," + YYi[i]+ "\n");
			out.println("" + XXi[i]*scale + "," + YYi[i]*scale+ "," +ZZi[i]*scale + "," + Di[i]*scale);

			
			//get_local_Nprofile(float xx, float yy, float[] N);
			
			
			}
		out.close();

		} 
		catch (IOException e){
			e.printStackTrace();
		}		
	}


	/******************************************************************************************
	 * Get unit tangent vector
	 */	
	private float[] get_tangent(float x1, float y1, float x2, float y2)
		{
		float[] TT = new float[2];
		TT[0] = x2 - x1;
		TT[1] = y2 - y1;
		float l = (float)Math.sqrt(TT[0]*TT[0] + TT[1]*TT[1]);
		TT[0] = TT[0]/l;
		TT[1] = TT[1]/l;
	
		return TT;
	}

	/******************************************************************************************
	* Objective function for non linear optimisation
	*/	
	class ObjFunction implements MultivariateFunction

		{
			// need XXi,YYi,ZZi as global
			// need parameter for ax_x
			
			public double value(double[] var)
			{
				double E_im = 0;
				double E_int = 0;
				double alpha = 0.;
				double beta = 100000.;
				int n_nodes = var.length / 2;
				int n_slice = IM.length;
				int n_profile = IM[0][0].length;
				for (int i=0;i<n_nodes;i++)
				{
					E_im = 0;//Math.pow(var[i]-n_profile/2.,2) + Math.pow(var[n_nodes+i]-n_profile/2.,2);
					for (int j=0;j<n_slice;j++)
					{
						double co = Math.cos(A[j]);
						double si = Math.sin(A[j]);
						double cx = ax_x[0] + ZZi[i]*ax_x[1];
						double t_x = (var[i] - cx)*co - (var[n_nodes+i] - cx)*si + cx;
						double t_0x = (XXi[i] - cx)*co - (YYi[i] - cx)*si + cx;
						
						double dt = t_x-t_0x;
						double II = 0;
						if (dt>-n_profile/2.+0.5 & dt<n_profile/2.-0.5)
							{
							int I0 = IM[j][i][(int)(dt+n_profile/2.)];
							int I1 = IM[j][i][(int)(dt+n_profile/2.)+1];
							II = (double)I0 + dt * ((double)(I1-I0)); 
							//E_im -= I;//IM[j][i][(int)(dt+n_profile/2.+0.5)];
							}

						break;
						//E_im += Math.pow(CXik[i][j] - t_x, 2);

					}
					if (i>0 & i<n_nodes-1)  // elongation and bending energy
					{
						double v01 = Math.sqrt(Math.pow((XXi[i]-XXi[i-1]),2) + Math.pow((YYi[i]-YYi[i-1]),2));
						double v02 = Math.sqrt(Math.pow((XXi[i+1]-XXi[i]),2) + Math.pow((YYi[i+1]-YYi[i]),2));
						
						double v11 = Math.sqrt(Math.pow((var[i]-var[i-1]),2) + Math.pow((var[n_nodes+i]-var[n_nodes+i-1]),2));
						double v12 = Math.sqrt(Math.pow((var[i+1]-var[i]),2) + Math.pow((var[n_nodes+i+1]-var[n_nodes+i]),2));
						
						double cross = (var[i]-var[i-1])*(var[n_nodes+i+1]-var[n_nodes+i]);
						cross -= (var[n_nodes+i]-var[n_nodes+i-1])*(var[i+1]-var[i]); 
						cross /= (v11*v12);
						double curve = Math.abs(Math.asin(cross)/(v11+v12));
						E_int += beta*curve; //alpha*Math.pow((v01-v11)/v01,2) + 

					}
					if (i==n_nodes-1) // just the elongation energy at the last node
					{
						double v01 = Math.sqrt(Math.pow((XXi[i]-XXi[i-1]),2) + Math.pow((YYi[i]-YYi[i-1]),2));
						double v11 = Math.sqrt(Math.pow((var[i]-var[i-1]),2) + Math.pow((var[n_nodes+i]-var[n_nodes+i-1]),2));
	
						//E_int += alpha*Math.pow((v01-v11)/v01,2);
					}
				}
				//IJ.log("test: "  + E_im + " , " + E_int );
				return E_im ;//+ Math.abs(E_int*E_int*E_int);
			}
		}
		

}

