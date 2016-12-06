package fr.curie.ROMA;

import fr.curie.ROMA.GESignature;
import fr.curie.ROMA.GMTReader;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

import vdaoengine.analysis.*;
import vdaoengine.data.*;
import vdaoengine.data.io.*;
import vdaoengine.utils.Algorithms;
import vdaoengine.utils.Utils;
import vdaoengine.utils.VSimpleFunctions;
import vdaoengine.utils.VSimpleProcedures;


public class ModuleActivityAnalysis {

	// Constants
	
	static int STANDARD_GMT = 0;
	static int GMT_WITH_WEIGHTS = 1;
	static int WEIGHT_MATRIX = 2;
	public static int PCA_STANDARD = 0;
	public static int PCA_FIXED_CENTER = 1;
	public static int NON_LINEAR = 2;
	
	// Options
	
	static int typeOfModuleFile = GMT_WITH_WEIGHTS;
	//static int typeOfModuleFile = STANDARD_GMT;
	//static int typeOfPCAUsage = PCA_STANDARD;
	static int typeOfPCAUsage = PCA_FIXED_CENTER;
	
	static boolean robustPCAcalculation = true;
	static boolean robustPCAcalculationForSampling = false;
	
	static int minimalNumberOfGenesInModule = 10;
	static int maximalNumberOfGenesInModule = 1000;
	static int minimalNumberOfGenesInModuleFound = 8;	
	
	static String outputFolder = null;
	static String dataFile = null;
	static String sampleFile = null;
	static String moduleFile = null;
	static String activitySignsFile = null;
	
	static float mostContributingGenesZthreshold = 1f;
	static float diffSpotGenesZthreshold = 1f;
	
	static String fieldForAveraging = null;
	static String fieldValueForAveraging = null;
	
	static String fieldForDiffAnalysis = null;
	static String fieldValuesForDiffAnalysis = null;
	
	static boolean centerData = false;
	static boolean doubleCenterData = false;
	static int fillMissingValues = -1;
	
	 // Internal parameters
	 static float outlierThreshold = 2f;
	 static float correlationThreshold = 0.6f;
	 static String commonFactorForPartialCorrelation = null;
	 static float graphicalOutputThreshold = 0.05f;
	 static int outlierDimension = 3;
	 
	 // Sampling parameters
	 static int numberOfPermutations = 100;
	 static int samplingGeneSetSizes[] = {10,15,20,30,40,50,100,200};
	 static int numberOfGeneSetSizesToSample = 5;
	 
	 //static int 
	 
	 static boolean saveDecomposedFiles = false;
	 //static boolean produceTxtTables = false;

	

	// Working objects
	
	public static Vector<Metagene> signatures = new Vector<Metagene>();
	static VDataTable table = null;
	static VDataTable alldata = null;
	static VDataTable moduleTable = null;
	static VDataTable projectionTable = null;
	static VDataTable sampleTable = null;
	static HashMap moduleWeightsSum = null;
	
	static Vector<float[][]> As = new Vector<float[][]>();
	static Vector<float[][]> Ss = new Vector<float[][]>();
	static Vector<float[][]> globalProjections = new Vector<float[][]>();	
	static Vector<double[]> explainedVariances = new Vector<double[]>();
	static Vector<Vector<double[]>> explainedVariances_randomDistributions = null;
	
	
	static Vector<VDataTable> tables = new Vector<VDataTable>();
	
	public static String geneField = "GeneSymbol";
	

	public static void main(String[] args) {
		try{
			
			/*VDataTable vt = VDatReadWrite.LoadFromSimpleDatFile("C:/Datas/ROMA/data/mosaic/cliques/moduletable_simple.txt", true, "\t");
			VSimpleProcedures.findAllNumericalColumns(vt);
			VSimpleFunctions.makeCorrelationTableCorrectedForCommonFactor(vt, 0.6f, "Lessnick_EWS-FLI_All_signature", "C:/Datas/ROMA/data/mosaic/cliques/module_correlations_0.6_EWS_FLI1.txt");
			System.exit(0);*/

			//Compute metagene-based correlation graph
			/*File files[] = new File("C:/Datas/ROMA/data/single_cell/single_cell100/").listFiles();
			//Vector<VDataSet> vts = new Vector<VDataSet>();
			Vector<float[]> fs = new Vector<float[]>();
			Vector<String> names = new Vector<String>();
			int count=0;
			for(File f: files){
				if(f.getName().contains("_projs")){
					count++; 
					VDataTable vt = VDatReadWrite.LoadFromSimpleDatFile(f.getAbsolutePath(), true, "\t");
					float ff[] = new float[vt.rowCount];
					for(int k=0;k<vt.rowCount;k++)
						ff[k] = Float.parseFloat(vt.stringTable[k][0]);
					fs.add(ff);
					names.add(f.getName().substring(0, f.getName().length()-10));
					System.out.println(count+": "+names.get(names.size()-1));
				}
			}
			//System.out.println("Loaded");
			FileWriter fw = new FileWriter("C:/Datas/ROMA/data/single_cell/single_cell100/metagene_correlations0.8.txt");
			fw.write("SET1\tSET2\tCORR\tABSCORR\n");
			count = 0;
			for(int i=0;i<names.size();i++){
				
				//VDataTable vti = VDatReadWrite.LoadFromSimpleDatFile("C:/Datas/ROMA/data/single_cell/single_cell100/"+names.get(i)+"_projs.txt", true, "\t");
				//VSimpleProcedures.findAllNumericalColumns(vti);
				//VDataSet dsi = VSimpleProcedures.SimplyPreparedDatasetWithoutNormalization(vti, -1);
				
				for(int j=i+1;j<names.size();j++){
					
					//VDataTable vtj = VDatReadWrite.LoadFromSimpleDatFile("C:/Datas/ROMA/data/single_cell/single_cell100/"+names.get(j)+"_projs.txt", true, "\t");
					//VSimpleProcedures.findAllNumericalColumns(vtj);
					//VDataSet dsj = VSimpleProcedures.SimplyPreparedDatasetWithoutNormalization(vtj, -1);
				
				//float fi[] = new float[dsi.pointCount];
				//float fj[] = new float[dsj.pointCount];
				float fi[] = fs.get(i);
				float fj[] = fs.get(j);
				
				//for(int k=0;k<dsi.pointCount;k++) fi[k] = dsi.massif[k][0];
				//for(int k=0;k<dsj.pointCount;k++) fj[k] = dsj.massif[k][0];
				float corr = VSimpleFunctions.calcCorrelationCoeff(fi, fj);
				if(Math.abs(corr)>=0.8){
					fw.write(names.get(i)+"\t"+names.get(j)+"\t"+corr+"\t"+Math.abs(corr)+"\n");
				}
				}
				System.out.println(count+": "+names.get(i));
				count++;
			}
			fw.close();
			System.exit(0);*/
			
			// Reannotation
			//VDataTable vt1 = VDatReadWrite.LoadFromVDatFile("c:/datas/moduleactivities/data/bcpublic/bc4c.dat");
			/*VDataTable vt1 = VDatReadWrite.LoadFromVDatFile("c:/datas/moduleactivities/data/bcprivate/bc1na.dat");
			VDataTable annot = VDatReadWrite.LoadFromSimpleDatFile("c:/datas/moduleactivities/annot4.txt", true, "\t");
			VDataTable annot1 = VDatReadWrite.LoadFromSimpleDatFile("c:/datas/moduleactivities/annot_95av2.txt", true, "\t");
			annot.makePrimaryHash("Probeset");
			annot1.makePrimaryHash("Probeset");
			for(int i=0;i<vt1.rowCount;i++){
				String ps = vt1.stringTable[i][vt1.fieldNumByName("Probe")];
				String gs = vt1.stringTable[i][vt1.fieldNumByName("GeneSymbol")];
				String gt = vt1.stringTable[i][vt1.fieldNumByName("GeneTitle")];
				Vector<Integer> rows = annot.tableHashPrimary.get(ps);
				Vector<Integer> rows1 = annot1.tableHashPrimary.get(ps);
				if(((rows==null)||(rows.size()==0))&&((rows1==null)||(rows1.size()==0)))
					System.out.println(ps+" not found");
				else{
					String gs1 = null;
					String gt1 = null;
					if(rows!=null){
						gs1 = annot.stringTable[rows.get(0)][annot.fieldNumByName("GeneSymbol")];
						gt1 = annot.stringTable[rows.get(0)][annot.fieldNumByName("GeneTitle")];
					}else{
						gs1 = annot1.stringTable[rows1.get(0)][annot1.fieldNumByName("GeneSymbol")];
						gt1 = annot1.stringTable[rows1.get(0)][annot1.fieldNumByName("GeneTitle")];
					}
					vt1.stringTable[i][vt1.fieldNumByName("GeneTitle")] = gt1;					
					if(!gs.equals(gs1)){
						vt1.stringTable[i][vt1.fieldNumByName("GeneSymbol")] = gs1;
						System.out.println(gs+" -> "+gs1);
					}
				}
			}
			//VDatReadWrite.saveToVDatFile(vt1,"c:/datas/moduleactivities/data/bcpublic/bc4c_.dat");
			VDatReadWrite.saveToVDatFile(vt1,"c:/datas/moduleactivities/data/bcprivate/bc1na_.dat");
			System.exit(0);
			*/
			/*moduleFile = "c:/datas/ModuleActivities/modules_proteins_050209.txt";
			//moduleFile = "c:/datas/ModuleActivities/rbmodules.gmt";

			dataDatFileFolder = "c:/datas/ModuleActivities/data/bcpublic/";
			datFile = "bc4c.dat";
			fieldNumberForAveraging = 5;
			fieldValueForAveraging = "normal";
			fieldNumberForDiffAnalysis = 5;
			fieldValuesForDiffAnalysis = "noninvasive;invasive";
			
			/*dataDatFileFolder = "c:/datas/ModuleActivities/data/bcprivate/";
			datFile = "bc1na.dat";
			//fieldNumberForAveraging = 15;
			fieldValueForAveraging = null;
			fieldNumberForDiffAnalysis = 15;
			fieldValuesForDiffAnalysis = "noninvasive;invasive";*/
			//mostContributingGenesZthreshold = 0.8f;
			
			//activitySignsFile = dataDatFileFolder+"moduleSigns.txt";

			
			for(int i=0;i<args.length;i++){
				if(args[i].equals("-typeOfModuleFile"))
					typeOfModuleFile = Integer.parseInt(args[i+1]);
				if(args[i].equals("-typeOfPCAUsage"))
					typeOfPCAUsage = Integer.parseInt(args[i+1]);
				if(args[i].equals("-dataFile"))
					dataFile = args[i+1];
				if(args[i].equals("-outputFolder"))
					outputFolder = args[i+1];
				if(args[i].equals("-sampleFile"))
					sampleFile = args[i+1];
				if(args[i].equals("-moduleFile"))
					moduleFile = args[i+1];
				if(args[i].equals("-saveDecomposedFiles"))
					saveDecomposedFiles = Integer.parseInt(args[i+1])==1;
				if(args[i].equals("-activitySignsFile"))
					activitySignsFile = args[i+1];
				if(args[i].equals("-mostContributingGenesZthreshold"))
					mostContributingGenesZthreshold = Float.parseFloat(args[i+1]);
				if(args[i].equals("-diffSpotGenesZthreshold"))
					diffSpotGenesZthreshold = Float.parseFloat(args[i+1]);
				if(args[i].equals("-correlationThreshold"))
					correlationThreshold = Float.parseFloat(args[i+1]);
				if(args[i].equals("commonFactorForPartialCorrelation"))
					commonFactorForPartialCorrelation = args[i+1];
				if(args[i].equals("-graphicalOutputThreshold"))
					graphicalOutputThreshold = Float.parseFloat(args[i+1]);
				if(args[i].equals("-fieldForAveraging"))
					fieldForAveraging = args[i+1];
				if(args[i].equals("-fieldValueForAveraging"))
					fieldValueForAveraging = args[i+1];
				if(args[i].equals("-fieldForDiffAnalysis"))
					fieldForDiffAnalysis = args[i+1];
				if(args[i].equals("-fieldValuesForDiffAnalysis"))
					fieldValuesForDiffAnalysis = args[i+1];
				if(args[i].equals("-robustPCA"))
					robustPCAcalculation = Integer.parseInt(args[i+1])==1;
				if(args[i].equals("-robustPCASampling"))
					robustPCAcalculationForSampling = Integer.parseInt(args[i+1])==1;
				if(args[i].equals("-centerData"))
					centerData = Integer.parseInt(args[i+1])==1;
				if(args[i].equals("-doubleCenterData"))
					doubleCenterData = Integer.parseInt(args[i+1])==1;
				if(args[i].equals("-fillMissingValues"))
					fillMissingValues = Integer.parseInt(args[i+1]);
				if(args[i].equals("-minimalNumberOfGenesInModule"))
					minimalNumberOfGenesInModule = Integer.parseInt(args[i+1]);
				if(args[i].equals("-minimalNumberOfGenesInModuleFound"))
					minimalNumberOfGenesInModuleFound = Integer.parseInt(args[i+1]);
				if(args[i].equals("-maximalNumberOfGenesInModule"))
					maximalNumberOfGenesInModule = Integer.parseInt(args[i+1]);
			    if(args[i].equals("-outlierThreshold"))
			          outlierThreshold = (new Float(args[i+1])).floatValue();
			    if(args[i].equals("-outlierDimension"))
			          outlierDimension = (new Integer(args[i+1])).intValue();
			    if(args[i].equals("-numberOfPermutations"))
			    	numberOfPermutations = (new Integer(args[i+1])).intValue();
			    if(args[i].equals("-numberOfGeneSetSizesToSample"))
			    	numberOfGeneSetSizesToSample = (new Integer(args[i+1])).intValue();
			    //if(args[i].equals("-produceTxtTables"))
			    //	produceTxtTables = Integer.parseInt(args[i+1])==1;
			}
			

			if((args.length==0)||(dataFile==null)||(moduleFile==null)){
			
			System.out.println(":: OPTIONS");
			System.out.println("REQUIRED:");
			System.out.println(":: -dataFile : data in tab-delimited format (required)");
			System.out.println(":: -moduleFile : name of the gmt module file  (required)");
			System.out.println("Optional:");
			System.out.println(":: -outputFolder : folder name for keeping the resulting files (by default it will be the folder of the data file)");
			System.out.println(":: -sampleFile : description of samples in tab-delimited txt format");
			System.out.println(":: -saveDecomposedFiles : save dat files for each module in the output folder");
			System.out.println(":: -produceNumericalTables : in addition to dat files produce purely numerical tables (e.g., for analysis in Matlab)");
			System.out.println(":: -outlierThreshold : threshold for determining outliers");
			System.out.println(":: -typeOfModuleFile : 0 - for standard GMT, 1 - for GMT with weights (default)");
			System.out.println(":: -typeOfPCAUsage : 0 - for standard PCA, 1 - for PCA with fixed center (default)");
			System.out.println(":: -robustPCA : 0 - all points are used, 1 - leave one out -based removal of outliers");
			System.out.println(":: -robustPCASampling : 0 - all points are used in random sampling, 1 - leave one out -based removal of outliers in random sampling (slow down calculations)");
			System.out.println(":: -centerData : 0 - do not center, 1 - center each line");
			System.out.println(":: -doubleCenterData : 0 - do not center, 1 - center each line and each column");
			//System.out.println(":: -activitySignsFile : file with definition of the module signs (optional)");
			System.out.println(":: -mostContributingGenesZthreshold : (default 1) threshold (z-value) used to print out the names of the genes most contributing to the component");
			System.out.println(":: -diffSpotGenesZthreshold : (default 1) threshold (t-test) used to print out the names of the differentially expressed genes");
			System.out.println(":: -correlationThreshold: (default 0.6) threshold used to cut the edges of the module activities correlation graph");
			System.out.println(":: -commonFactorForPartialCorrelation: name of the gene set used for ");
			System.out.println(":: -graphicalOutputThreshold: (default 0.05) threshold on p-value used limit the number of projection files (_proj.txt) in the output");
			System.out.println(":: -fieldForAveraging : (optional) number of the field column used for computing the average module activities (ex: 5)");
			System.out.println(":: -fieldValueForAveraging : (optional) value of the field used for computing the reference value of module activity (ex: \"normal\")");
			System.out.println(":: -fieldForDiffAnalysis : (optional) number of the field column used for differential analysis (ex: 5)");
			System.out.println(":: -fieldValuesForDiffAnalysis : (optional) values of the field column used for differential analysis separated by % symbol (ex: \"invasive%noninvasive\")");
			System.out.println(":: -numberOfPermutations : (optional) number of samples for empirical p-values estimation");
			System.out.println(":: -numberOfGeneSetSizesToSample : (optional) number of random gene set sizes to test for empirical p-values estimation");
			System.out.println(":: -minimalNumberOfGenesInModule: minimal size of the gene set");
			System.out.println(":: -maximalNumberOfGenesInModule: maximal size of the gene set");
			System.out.println(":: -minimalNumberOfGenesInModuleFound: minimal number of genes in a gene set found in dataset");
			
			}else{

			
			VDatReadWrite.useQuotesEverywhere = false;
			loadData();
			
			System.out.println("============================================");
			System.out.println((new Date()).toString()+"\n");
			System.out.println("outputFolder= "+outputFolder);
			System.out.println("dataFile= "+dataFile);
			System.out.println("sampleFile= "+sampleFile);
			System.out.println("moduleFile= "+moduleFile);
			System.out.println("centerData= "+centerData);			
			System.out.println("doubleCenterData= "+doubleCenterData);
			System.out.println("fillMissingValues= "+fillMissingValues);
			System.out.println("activitySignsFile= "+activitySignsFile);
			System.out.println("mostContributingGenesZthreshold= "+mostContributingGenesZthreshold);
			System.out.println("diffSpotGenesZthreshold= "+diffSpotGenesZthreshold);
			System.out.println("correlationThreshold= "+correlationThreshold);
			System.out.println("graphicalOutputThreshold= "+graphicalOutputThreshold);
			System.out.println("typeOfPCAUsage= "+typeOfPCAUsage);
			System.out.println("robustPCACalculation= "+robustPCAcalculation);
			System.out.println("robustPCACalculationForSampling= "+robustPCAcalculationForSampling);
			System.out.println("outlierThreshold= "+outlierThreshold);
			System.out.println("typeOfModuleFile= "+typeOfModuleFile);			
			System.out.println("fieldForAveraging= "+fieldForAveraging);
			System.out.println("fieldValueForAveraging= "+fieldValueForAveraging);
			System.out.println("fieldForDiffAnalysis= "+fieldForDiffAnalysis);
			System.out.println("fieldValuesForDiffAnalysis= "+fieldValuesForDiffAnalysis);
			System.out.println("minimalNumberOfGenesInModule= "+minimalNumberOfGenesInModule);			
			System.out.println("maximalNumberOfGenesInModule= "+maximalNumberOfGenesInModule);
			System.out.println("minimalNumberOfGenesInModuleFound = "+minimalNumberOfGenesInModuleFound);						
			System.out.println("============================================");
			
			PCAMethod.verboseMode = false;
			PCAMethodFixedCenter.verboseMode = false;
			
			System.out.println();
			System.out.println("============================================");
			System.out.println("      Decomposing dataset into modules");
			System.out.println("============================================");
			decomposeDataSet();
			
			writeGeneModuleMatrix();			

			Date d1 = new Date();
			if(numberOfPermutations>0){
				System.out.println();
				System.out.println("============================================");
				System.out.println("       Calculating sampling");
				System.out.println("============================================");
				explainedVariances_randomDistributions = calcPermutedModuleActivities(table);
				savePermutedModuleActivities(explainedVariances_randomDistributions);
			}
			
			System.out.println(((new Date()).getTime()-d1.getTime())/1000f+" secs spent.");			
			
			// Central module activity computation function
			System.out.println();
			System.out.println("============================================");
			System.out.println("        Compute Module Activities");
			System.out.println("============================================");
			System.out.flush();
			Date d = new Date();
			ComputeModuleActivities();
			writeGeneProjectionMatrix();
			System.out.println(((new Date()).getTime()-d.getTime())/1000f+" secs spent.");
			
			
			System.out.println();
			System.out.println("============================================");
			System.out.println("        Produce graphical output");
			System.out.println("============================================");
			System.out.flush();
			d = new Date();
			ProduceGraphicalOutput();
			System.out.println(((new Date()).getTime()-d.getTime())/1000f+" secs spent.");
			
			
			System.out.println();
			System.out.println("============================================");
			System.out.println("        Top contributing genes determination");
			System.out.println("============================================");
			
			findMostContributingGenes();

			System.out.println();
			System.out.println("============================================");
			System.out.println("       Creating module activity table");
			System.out.println("============================================");
			
			moduleActivityTable();

			
			if(fieldForAveraging!=null){
			System.out.println();
			System.out.println("============================================");
			System.out.println("     Calculating average module activities");
			System.out.println("============================================");
			
			calcAverageModuleActivities(fieldForAveraging,fieldValueForAveraging);
			}
			
			if(fieldForDiffAnalysis!=null)if(fieldValuesForDiffAnalysis!=null){
			System.out.println();
			System.out.println("============================================");
			System.out.println("     Differential genes determination");
			System.out.println("============================================");
			
			StringTokenizer st = new StringTokenizer(fieldValuesForDiffAnalysis,"#");
			String lab1s = st.nextToken();
			String lab2s = st.nextToken();
			Vector lab1 = new Vector();
			Vector lab2 = new Vector();
			st = new StringTokenizer(lab1s,",");
			while(st.hasMoreTokens())
				lab1.add(st.nextToken());
			st = new StringTokenizer(lab2s,",");
			while(st.hasMoreTokens())
				lab2.add(st.nextToken());
			findDiffGenes(fieldForDiffAnalysis,lab1,lab2);			
			}
			}
			
			
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public static void loadData() throws Exception{
		//table = vdaoengine.data.io.VDatReadWrite.LoadFromVDatFile(dataDatFileFolder+dataFile);
		
		table = vdaoengine.data.io.VDatReadWrite.LoadFromSimpleDatFile(dataFile, true, "\t");
		TableUtils.findAllNumericalColumns(table);
		
		if(outputFolder==null){
			outputFolder = (new File(dataFile)).getParent()+"/out/";
		}
		
		if(!new File(outputFolder).exists()){
			if(!new File(outputFolder).mkdir())
				System.out.println("FATAL ERROR: can not create output folder "+outputFolder+" !");
		}
		
		if(!outputFolder.endsWith("\\"))if(!outputFolder.endsWith("/"))
			outputFolder = outputFolder+"/";
		
		
		if(fillMissingValues>0){
			System.out.println("Filling missing values using first "+fillMissingValues+" principal components...");
			table = TableUtils.fillMissingValues(table,fillMissingValues);
		}
		
		
		if(centerData){
			System.out.println("Centering data (zero mean for any gene)...");
			table = TableUtils.normalizeVDat(table, true, false);
			String _fn = dataFile;
			_fn = dataFile.substring(0, dataFile.length()-4);
			VDatReadWrite.saveToVDatFile(table, _fn+"_centered.dat");
		}
		if(doubleCenterData){
			System.out.println("Double centering data (zero mean for any gene and for any sample)...");
			VDataSet vd = VSimpleProcedures.SimplyPreparedDatasetWithoutNormalization(table, -1);
			float mas[][] = TableUtils.doubleCenterMatrix(vd.massif);
			for(int l=0;l<table.rowCount;l++){
				int k=0;
				for(int s=0;s<table.colCount;s++){
					table.stringTable[l][s] = ""+mas[l][k];
					if(table.fieldTypes[s]==table.NUMERICAL)
						k++;
				}
			}
			String _fn = dataFile;
			_fn = dataFile.substring(0, dataFile.length()-4);
			VDatReadWrite.saveToVDatFile(table, _fn+"_doublecentered.dat");
		}
		
		
		if(table.fieldNumByName(geneField)==-1)
			for(int i=0;i<table.colCount;i++){
				if(table.fieldTypes[i]==table.STRING)
					geneField = table.fieldNames[i];
			}
		table.makePrimaryHash(geneField);
		
		
		if(typeOfModuleFile==STANDARD_GMT){
		Vector<GESignature> sigs = GMTReader.readGMTDatabase(moduleFile,minimalNumberOfGenesInModule,maximalNumberOfGenesInModule);
		for(int i=0;i<sigs.size();i++){
			Metagene mg = new Metagene(sigs.get(i));
			mg.probeSets = mg.geneNames;			
			mg.initializeWeightsByOnes();
			signatures.add(mg);
		}
		}
		if(typeOfModuleFile==GMT_WITH_WEIGHTS){
			signatures = GMTReader.readGMTDatabaseWithWeights(moduleFile,minimalNumberOfGenesInModule,maximalNumberOfGenesInModule);
			for(int i=0;i<signatures.size();i++){
				Metagene mg = new Metagene(signatures.get(i));
				mg.probeSets = mg.geneNames;			
			}
		}
		if(activitySignsFile!=null){
			try{
				LineNumberReader lr = new LineNumberReader(new FileReader(activitySignsFile));
				String s = null;
				while((s=lr.readLine())!=null){
					StringTokenizer st = new StringTokenizer(s,"\t");
					String name = st.nextToken();
					String sign = st.nextToken();
					int k = GESignature.findSignatureByName(name, signatures);
					if(k!=-1){
						Metagene mg = signatures.get(k);
						mg.correctedSign = Float.parseFloat(sign);
					}
				}
			}catch(Exception e){
				e.printStackTrace();
			}
		}
		
		if(sampleFile!=null){
			sampleTable = VDatReadWrite.LoadFromSimpleDatFile(sampleFile, true, "\t");
			sampleTable.makePrimaryHash(sampleTable.fieldNames[0]);
		}
	}
	
	public static void decomposeDataSet(){
		Vector allgenes = new Vector();
		
		Vector<Metagene> signatureCompleted = new Vector<Metagene>();
		
		for(int i=0;i<signatures.size();i++){
			Metagene sig = (Metagene)signatures.get(i);
			System.out.println();
			System.out.println(sig.name);			
			VDataTable moddata = VSimpleProcedures.selectRowsFromList(table, sig.geneNames);
			if(moddata.rowCount>=minimalNumberOfGenesInModuleFound){

			if(robustPCAcalculation)	
				moddata = removeOutliers(moddata, sig, true);
			
			/*if(robustPCAcalculation){
				VDataSet ds = VSimpleProcedures.SimplyPreparedDataset(moddata,-1);
			    Vector<String> v = TableUtils.determinePC1Outliers(ds,geneField,outlierThreshold);
			    System.out.print(sig.name+" outliers:\t");
			    for(int j=0;j<v.size();j++){
			        	sig.geneNames.remove((String)v.get(j));
			        	System.out.print((String)v.get(j)+"\t");
			    }
			    System.out.println();
			    moddata = VSimpleProcedures.selectRowsFromList(table, sig.geneNames);
			
			}*/
			
			moddata = VSimpleProcedures.substituteRowsByStatistics(moddata, geneField, 4);
			tables.add(moddata);
			String fn = new String(sig.name);
			fn = Utils.replaceString(fn, ":", "_");
			fn = Utils.replaceString(fn, "/", "_");
			
			VDataSet vd = VSimpleProcedures.SimplyPreparedDataset(moddata, -1);
			assignWeights(sig,vd,moddata);
			if(moddata.fieldNumByName("WEIGHT")!=-1)
				moddata.fieldTypes[moddata.fieldNumByName("WEIGHT")] = moddata.NUMERICAL;
		
			if(saveDecomposedFiles){
			//VDatReadWrite.saveToVDatFile(moddata, outputFolder+fn+".dat");
			//if(produceDatTables){
				VDatReadWrite.saveToSimpleDatFile(moddata, outputFolder+fn+".txt");
			//}
			VDataTable transp = moddata.transposeTable(moddata.fieldNames[0]);
			if(sampleTable!=null)
				transp = VSimpleProcedures.MergeTables(transp, "NAME", sampleTable, sampleTable.fieldNames[0], "_");
			//VDatReadWrite.saveToVDatFile(transp, outputFolder+fn+"_T.dat");
			//if(produceDatTables){
				VDatReadWrite.saveToSimpleDatFile(transp, outputFolder+fn+"_T.txt");
			//}
			}
			
			if(moddata.fieldNumByName("WEIGHT")!=-1)
				moddata.fieldTypes[moddata.fieldNumByName("WEIGHT")] = moddata.STRING;
			
			Vector notFound = new Vector();
			for(int j=0;j<sig.geneNames.size();j++){
				String id = (String)sig.geneNames.get(j);
				if(allgenes.indexOf(id)<0)
					allgenes.add(id);
				if(table.tableHashPrimary.get(id)==null){
					notFound.add(id);
				}
			}
			if(notFound.size()>0){
				System.out.print("NOT FOUND:\t");
				for(int j=0;j<notFound.size();j++)
					System.out.print((String)notFound.get(j)+"\t");
				System.out.println();
			}
			signatureCompleted.add(sig);
			}else{
				System.out.println("Signature is not complete (less than "+minimalNumberOfGenesInModuleFound+" genes)");
			}
		}
		
		System.out.println("\n\n"+signatureCompleted.size()+" signatures are complete from "+signatures.size());
		signatures = signatureCompleted;
		
		System.out.println("\nSelecting genes from all modules...");
		alldata = VSimpleProcedures.selectRowsFromList(table, allgenes);
		alldata = VSimpleProcedures.substituteRowsByStatistics(alldata, geneField, 4);
		
		alldata.addNewColumn("MODULE", "", "", alldata.STRING, "");
		alldata.makePrimaryHash(geneField);
		
		for(int i=0;i<signatures.size();i++){
			Metagene sig = (Metagene)signatures.get(i);
			for(int j=0;j<sig.probeSets.size();j++){
				String gn = (String)sig.probeSets.get(j);
				Vector<Integer> rows = alldata.tableHashPrimary.get(gn);
				if(rows!=null)
				for(int k=0;k<rows.size();k++){
					alldata.stringTable[rows.get(k)][alldata.fieldNumByName("MODULE")] += sig.name+";";
				}
			}
		}
		
		VDatReadWrite.saveToVDatFile(alldata, outputFolder+"allgenes.dat");
		VDataTable transp = alldata.transposeTable(alldata.fieldNames[0]);
		if(sampleTable!=null)
			transp = VSimpleProcedures.MergeTables(transp, "NAME", sampleTable, sampleTable.fieldNames[0], "_");
		VDatReadWrite.saveToVDatFile(transp, outputFolder+"allgenes_T.dat");
	}
	
	public static void findMostContributingGenes(){
		// Create report on hot-spot genes
		VDataSet vdglobal = VSimpleProcedures.SimplyPreparedDataset(table, -1);
		moduleWeightsSum = new HashMap();
		
		for(int i=0;i<tables.size();i++){
			String mname = ((GESignature)signatures.get(i)).name;
			System.out.print("\nModule "+mname+"\t"+((GESignature)signatures.get(i)).geneNames.size()+" genes\t");
			VDataTable dat = (VDataTable)tables.get(i);
			
			VDataSet vdp = new VDataSet();
			vdp.coordCount = 1;
			vdp.pointCount = Ss.get(i).length;
			vdp.massif = new float[vdp.pointCount][vdp.coordCount];
			for(int k=0;k<vdp.pointCount;k++) vdp.massif[k][0] = (float)Ss.get(i)[k][0];
			
			vdp.calcStatistics();
			vdp.simpleStatistics.calcMedians();
			Vector genePlusNames = new Vector();
			Vector genePlusScores = new Vector();
			Vector geneMinusNames = new Vector();
			Vector geneMinusScores = new Vector();
			for(int j=0;j<vdp.pointCount;j++){
				float f[] = vdp.getVector(j);
				float median = vdp.simpleStatistics.getMedian(0);
				float stdev = vdp.simpleStatistics.getStdDev(0);
				if(typeOfPCAUsage == PCA_FIXED_CENTER){
					median = 0;
					stdev = vdp.simpleStatistics.getStdDev0(0);
				}
				float z = (f[0]-median)/stdev;
				if(z<=-mostContributingGenesZthreshold){
					geneMinusNames.add(dat.stringTable[j][dat.fieldNumByName(geneField)]);
					geneMinusScores.add(new Float(z));
					//geneMinusScores.add(new Float(f[0]));
				}
				if(z>=+mostContributingGenesZthreshold){
					genePlusNames.add(dat.stringTable[j][dat.fieldNumByName(geneField)]);
					genePlusScores.add(new Float(z));
					//genePlusScores.add(new Float(f[0]));
				}
			}
			// make output
			DecimalFormat fm = new DecimalFormat("#.#");
			float scores[] = new float[geneMinusScores.size()];
			for(int j=0;j<geneMinusScores.size();j++) scores[j] = ((Float)geneMinusScores.get(j)).floatValue();
			int ind_minus[] = Algorithms.SortMass(scores);
			scores = new float[genePlusScores.size()];
			for(int j=0;j<genePlusScores.size();j++) scores[j] = ((Float)genePlusScores.get(j)).floatValue();
			int ind_plus[] = Algorithms.SortMass(scores);
			for(int j=0;j<geneMinusScores.size();j++){
				float z = ((Float)geneMinusScores.get(ind_minus[j])).floatValue();
				System.out.print((geneMinusNames.get(ind_minus[j]))+"("+fm.format(z)+")"+",");
			}
			System.out.print("\t");
			for(int j=0;j<genePlusScores.size();j++){
				float z = ((Float)genePlusScores.get(ind_plus[j])).floatValue();
				System.out.print((genePlusNames.get(ind_plus[j]))+"("+fm.format(z)+")"+",");
			}
			System.out.println();
		}
	}
	
	public static void findDiffGenes(String fieldname, Vector labs1, Vector labs2){
		DecimalFormat fm = new DecimalFormat("#.##");		
		Vector allGeneNames = new Vector();
		Vector allGeneScores = new Vector();
		
		int field = sampleTable.fieldNumByName(fieldname);
		if(field==-1){
			System.out.println("ERROR: field "+fieldname+" not found in "+sampleFile);
		}else{
		
		for(int i=0;i<tables.size();i++){
			String mname = ((GESignature)signatures.get(i)).name;
			VDataTable dat = (VDataTable)tables.get(i);
			Vector genePlusNames = new Vector();
			Vector genePlusScores = new Vector();
			Vector geneMinusNames = new Vector();
			Vector geneMinusScores = new Vector();
			for(int j=0;j<dat.rowCount;j++){
				VStatistics stat1 = new VStatistics(1); 
				VStatistics stat2 = new VStatistics(1);
				for(int k=0;k<dat.colCount;k++)if(dat.fieldTypes[k]==dat.NUMERICAL){
					float f[] = new float[1];
					f[0] = Float.parseFloat(dat.stringTable[j][k]);
					//String label = dat.fieldInfo[k][field];
					String sample = dat.fieldNames[k];
					String label = null;
					if(sampleTable.tableHashPrimary.get(sample)!=null)
						label = sampleTable.stringTable[sampleTable.tableHashPrimary.get(sample).get(0)][field];
					//System.out.println(sample+"\t"+label);
					if(labs1.indexOf(label)>=0){
						stat1.addNewPoint(f);
						//System.out.println("Added "+f[0]);
					}
					if(labs2.indexOf(label)>=0){
						stat2.addNewPoint(f);
						//System.out.println("Added "+f[0]);
					}
				}
				stat1.calculate();
				stat2.calculate();
				float z = (float)((stat2.getMean(0)-stat1.getMean(0))/Math.sqrt(stat1.getStdDev(0)*stat1.getStdDev(0)+stat2.getStdDev(0)*stat2.getStdDev(0)));
				String gn = dat.stringTable[j][dat.fieldNumByName(geneField)];
				if(allGeneNames.indexOf(gn)<0){
					allGeneNames.add(gn);
					allGeneScores.add(new Float(z));
				}
				if(z<=-diffSpotGenesZthreshold){
					geneMinusNames.add(gn);
					geneMinusScores.add(new Float(z));
				}
				if(z>=+diffSpotGenesZthreshold){
					genePlusNames.add(gn);
					genePlusScores.add(new Float(z));
				}
			}
			// make output
			if((genePlusScores.size()>0)||(geneMinusScores.size()>0)){
			System.out.println("Module "+mname+"\t"+((GESignature)signatures.get(i)).geneNames.size()+" genes\t");
			float scores[] = new float[geneMinusScores.size()];
			for(int j=0;j<geneMinusScores.size();j++) scores[j] = ((Float)geneMinusScores.get(j)).floatValue();
			int ind_minus[] = Algorithms.SortMass(scores);
			scores = new float[genePlusScores.size()];
			for(int j=0;j<genePlusScores.size();j++) scores[j] = ((Float)genePlusScores.get(j)).floatValue();
			int ind_plus[] = Algorithms.SortMass(scores);
			for(int j=0;j<geneMinusScores.size();j++){
				float z = ((Float)geneMinusScores.get(ind_minus[j])).floatValue();
				System.out.print((geneMinusNames.get(ind_minus[j]))+"("+fm.format(z)+")"+"\t");
			}
			for(int j=0;j<genePlusScores.size();j++){
				float z = ((Float)genePlusScores.get(ind_plus[j])).floatValue();
				System.out.print((genePlusNames.get(ind_plus[j]))+"("+fm.format(z)+")"+"\t");
			}
			System.out.println();
			}

		}
		float absscores[] = new float[allGeneNames.size()];
		for(int i=0;i<allGeneScores.size();i++){
			float z = ((Float)allGeneScores.get(i)).floatValue();
			//absscores[i] = Math.abs(z);
			absscores[i] = z;
		}
		int inds[] = vdaoengine.utils.Algorithms.SortMass(absscores);
		try{
			FileWriter fw = new FileWriter(outputFolder+"diffgenes.xls");
			for(int i=inds.length-1;i>=0;i--){
				fw.write(((String)allGeneNames.get(inds[i]))+"\t"+fm.format(((Float)allGeneScores.get(inds[i])).floatValue())+"\n");
			}
			fw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
		}
	}
	
	public static Vector<Vector<double[]>> calcPermutedModuleActivities(VDataTable global_table){
		
		int minsize = Integer.MAX_VALUE;
		int maxsize = -1;
		for(Metagene m: signatures){
			if(m.geneNames.size()>maxsize) maxsize = m.geneNames.size();
			if(m.geneNames.size()<minsize) minsize = m.geneNames.size();
		}
		float dsize = ((float)Math.log10(maxsize)-(float)Math.log10(minsize))/(float)numberOfGeneSetSizesToSample;
		samplingGeneSetSizes = new int[numberOfGeneSetSizesToSample+1];
		for(int i=0;i<numberOfGeneSetSizesToSample+1;i++)
			samplingGeneSetSizes[i] = (int)(minsize*Math.pow(10,i*dsize));
		
		
		Random r = new Random();
		Vector<Vector<double[]>> randomExplainedVariances = new Vector<Vector<double[]>>();
		VDataSet vdglobal = VSimpleProcedures.SimplyPreparedDataset(global_table, -1);
		Vector<float[][]> rAs = new Vector<float[][]>();
		Vector<float[][]> rSs = new Vector<float[][]>();
		for(int i=0;i<samplingGeneSetSizes.length;i++){
			int size = samplingGeneSetSizes[i];
			System.out.println("Random set size = "+size);
			Vector<Metagene> random_sets = new Vector<Metagene>();
			Vector<VDataTable> random_tables = new Vector<VDataTable>();
			for(int j=0;j<numberOfPermutations;j++){
				if(j==((int)((float)j/10f))*10)
					System.out.print((j)+"\t");
				Metagene mg = new Metagene();
				mg.name = "RND"+size+"_"+(j+1);
				mg.geneNames = new Vector<String>();
				int k=0;
				while(k<size){
					int row = r.nextInt(global_table.rowCount);
					String geneName = global_table.stringTable[row][global_table.fieldNumByName(geneField)];
					if(!mg.geneNames.contains(geneName)){
						mg.geneNames.add(geneName);
						k++;
					}
				}
				mg.probeSets = mg.geneNames;
				mg.initializeWeightsByOnes();
				random_sets.add(mg);
				VDataTable rtab = VSimpleProcedures.selectRowsFromList(global_table, mg.geneNames, false);
				
				// removing outliers
				if(robustPCAcalculationForSampling)
					rtab = removeOutliers(rtab,mg, false);
				
				random_tables.add(rtab);
			}
			System.out.println();
			Vector<double[]> rev = calcModuleActivities(vdglobal, random_tables, random_sets, rAs, rSs, null, true,false);
			Vector<double[]> rev1 = new Vector<double[]>();
			for(int k=0;k<rev.size();k++){
				double l1l2[] = rev.get(k);
				double l1l2l12[] = new double[3];
				l1l2l12[0] = l1l2[0];
				l1l2l12[1] = l1l2[1];
				l1l2l12[2] = l1l2[0]/l1l2[1];
				rev1.add(l1l2l12);
			}
			randomExplainedVariances.add(rev1);
		}
		return randomExplainedVariances;
	}
	
	public static void savePermutedModuleActivities(Vector<Vector<double[]>> randomExplainedVariances) throws Exception{
		FileWriter fw = new FileWriter(outputFolder+"sampling.txt");
		for(int i=0;i<randomExplainedVariances.size();i++){
			fw.write("SIZE "+samplingGeneSetSizes[i]+"\n");
			fw.write("================================\n");
			fw.write("L1\tL2\tL1/L2\n");
			for(int j=0;j<randomExplainedVariances.get(i).size();j++){
				double l1 = randomExplainedVariances.get(i).get(j)[0];
				double l2 = randomExplainedVariances.get(i).get(j)[1];
				fw.write(l1+"\t"+l2+"\t"+(l1/l2)+"\n");
			}
			fw.write("\n\n");
		}
		fw.close();
	}
	
	/*
	 * The values argument contains values[0]=l1 values[1]=l2 values[2]=l1/l2
	 */
	public static float[] calcPValue4ExplainedVariance(Vector<Vector<double[]>> randomExplainedVariances, int setSize, double values[]){
		float pvalues[] = null;
		if(randomExplainedVariances!=null)if(randomExplainedVariances.size()>0){
		pvalues = new float[randomExplainedVariances.get(0).get(0).length];
		// first, determine the closest set size at logarithmic scale
		int k = -1;
		double mindist = Double.MAX_VALUE;
		for(int i=0;i<samplingGeneSetSizes.length;i++){
			double dist = Math.abs(Math.log(samplingGeneSetSizes[i])-Math.log(setSize));
			if(dist<mindist){
				k = i;
				mindist = dist;
			}
		}
		// Then compute the p-values
		Vector<double[]> distribution = randomExplainedVariances.get(k);
		pvalues = calcPValues(distribution, values);
		}else{
			pvalues = new float[values.length];
		}
		return pvalues;
	}
	
	public static float[] calcPValues(Vector<double[]> distribution, double values[]){
		float pvals[] = new float[distribution.get(0).length];
		for(int i=0;i<distribution.size();i++){
			double vector[] = distribution.get(i);
			for(int j=0;j<vector.length;j++){
				if(vector[j]>values[j])
					pvals[j] = pvals[j]+1;
			}
		}
		for(int i=0;i<pvals.length;i++)
			pvals[i]/=distribution.size();
		return pvals;
	}
	
	/*
	 * Central procedure for estimating module activities
	 */
	public static void ComputeModuleActivities(){
		VDataSet vdglobal = VSimpleProcedures.SimplyPreparedDataset(table, -1);
		explainedVariances = calcModuleActivities(vdglobal, tables, signatures, As, Ss, globalProjections, true, true);
	}
	
	/*
	 * computing module activities
	 */
	public static Vector<double[]> calcModuleActivities(VDataSet vdglobal, Vector<VDataTable> tableList, Vector<Metagene> metagenes, Vector<float [][]> As, Vector<float[][]> Ss, Vector<float[][]> globalSs, boolean dealWithSign, boolean verbose){
		Vector<double[]> acts = new Vector<double[]>();
		
		for(int i=0;i<tableList.size();i++){
			VDataTable dat = (VDataTable)tableList.get(i);
			if(verbose)
				System.out.println(""+(i+1)+"/"+tableList.size()+": "+metagenes.get(i).name+" "+dat.rowCount);
			Metagene mg = metagenes.get(i);
			
			//System.out.println("TP53 weight = "+mg.weightSpecified.get(mg.geneNames.indexOf("TP53")));
			
			//if(mg.name.equals("IC1_down"))
			//	System.out.println("IC1_down size = "+mg.geneNames.size()+", points in dataset = "+dat.rowCount);
			VDataSet vd = VSimpleProcedures.SimplyPreparedDataset(dat, -1);
			VDataSet vdnn = VSimpleProcedures.SimplyPreparedDataset(dat, -1,false,false);
			
			//System.out.println(mg.name);
			
			switch(typeOfPCAUsage){
			case 0/*PCA_STANDARD*/:
				assignWeights(mg,vd,dat);				
				PCAMethod pca = new PCAMethod();
				pca.setDataSet(vd);
				pca.calcBasis(2);
				vd.calcStatistics();				
				double explained_var[] = pca.calcDispersionsRelative(vd.simpleStatistics.totalDispersion*vd.simpleStatistics.totalDispersion);
				acts.add(explained_var);
				// correct the sign
				if(dealWithSign){
				for(int k=0;k<vd.coordCount;k++)
					pca.linBasis.basis[0][k]*=mg.correctedSign;
				VDataSet vdp = pca.getProjectedDataset();
				float weightSum = calcSpecifiedWeightSum(dat,vdp,mg);
				if(weightSum<0){
					for(int k=0;k<vd.coordCount;k++) pca.linBasis.basis[0][k]*=-1;
					vdp = pca.getProjectedDataset();
				}
				float A[][] = new float[vd.coordCount][2];
				float S[][] = new float[vd.pointCount][2];
				for(int k=0;k<vd.coordCount;k++) A[k][0]=(float)pca.linBasis.basis[0][k];
				for(int k=0;k<vd.coordCount;k++) A[k][1]=(float)pca.linBasis.basis[1][k];
				for(int k=0;k<vd.pointCount;k++) S[k][0]=vdp.massif[k][0];
				for(int k=0;k<vd.pointCount;k++) S[k][1]=vdp.massif[k][1];
				As.add(A); Ss.add(S);
				if(globalSs!=null){
					VDataSet vdgp = pca.getProjectedDataset(vdglobal);
					globalSs.add(vdgp.massif);
				}
				}
			break;
			case 1/*PCA_FIXED_CENTER*/:
				assignWeights(mg,vd,dat);
	             for(int kk=0;kk<vdnn.pointCount;kk++)
	                 vdglobal.processIntoSpace(vdnn.getVector(kk));
				PCAMethodFixedCenter pcaf = new PCAMethodFixedCenter();
				pcaf.setDataSet(vdnn);
				pcaf.calcBasis(2);
				vdnn.calcStatistics();
                double df[] = pcaf.calcDispersionsRelative(vdnn.simpleStatistics.totalSDVtozero*vdnn.simpleStatistics.totalSDVtozero);
                acts.add(df);
                
                VDataSet vdp_test = pcaf.getProjectedDataset();
                //float sumw = calcAllWeightSum(dat,vdp_test,mg);
                //System.out.println("Sum of weights = "+sumw);
                //System.out.println("["+pcaf.linBasis.basis[0][0]+"\t"+pcaf.linBasis.basis[0][1]+"\t"+pcaf.linBasis.basis[0][2]+"...]");
                
    			// correct the sign
    			if(dealWithSign){
    			for(int k=0;k<vd.coordCount;k++)
    				pcaf.linBasis.basis[0][k]*=mg.correctedSign;
    			VDataSet vdp = pcaf.getProjectedDataset();
    			float weightSum = calcSpecifiedWeightSum(dat,vdp,mg);
    			if(weightSum<0){
    				for(int k=0;k<vd.coordCount;k++) pcaf.linBasis.basis[0][k]*=-1;
    				vdp = pcaf.getProjectedDataset();
    			}
				float A[][] = new float[vd.coordCount][2];
				float S[][] = new float[vd.pointCount][2];
				for(int k=0;k<vd.coordCount;k++) A[k][0]=(float)pcaf.linBasis.basis[0][k];
				for(int k=0;k<vd.coordCount;k++) A[k][1]=(float)pcaf.linBasis.basis[1][k];
				for(int k=0;k<vd.pointCount;k++) S[k][0]=vdp.massif[k][0];
				for(int k=0;k<vd.pointCount;k++) S[k][1]=vdp.massif[k][1];
				As.add(A); Ss.add(S);
				if(globalSs!=null){
					VDataSet vdgp = pcaf.getProjectedDataset(vdglobal);
					globalSs.add(vdgp.massif);
				}
    			}
			break;
			}
			
		}
		return acts;
	}
	
	public static void moduleActivityTable() throws Exception{
		moduleTable = new VDataTable();
		moduleTable.copyHeader(table);
		moduleTable.rowCount = tables.size();
		moduleTable.stringTable = new String[moduleTable.rowCount][moduleTable.colCount];

		for(int i=0;i<tables.size();i++){
			
			Metagene mg = (Metagene)signatures.get(i);
			String mname = mg.name;
			moduleTable.stringTable[i][moduleTable.fieldNumByName(geneField)] = mname;
			VDataTable dat = (VDataTable)tables.get(i);
			VDataSet vd = VSimpleProcedures.SimplyPreparedDataset(dat, -1);
			int k=0;
			for(int j=0;j<dat.colCount;j++)if(vd.selector.isColumnSelected(j)){
				moduleTable.stringTable[i][j] = ""+As.get(i)[k][0];
				k++;
			}
		}
		
		
		if(explainedVariances_randomDistributions==null){
		System.out.println("MODULE\tL1\tL1/L2\tNUMBER_OF_GENES");
		for(int i=0;i<tables.size();i++){
			Metagene mg = (Metagene)signatures.get(i);
			float l1 = (float)explainedVariances.get(i)[0];
			float l2 = (float)explainedVariances.get(i)[1];
			System.out.println(mg.name+"\t"+l1+"\t"+(l1/l2)+"\t"+mg.geneNames.size());
		}
		
		FileWriter fw = new FileWriter(outputFolder+"module_scores.xls");
		fw.write("MODULE\tL1\tL1/L2\tNUMBER_OF_GENES\n");
		for(int i=0;i<tables.size();i++){
			Metagene mg = (Metagene)signatures.get(i);
			float l1 = (float)explainedVariances.get(i)[0];
			float l2 = (float)explainedVariances.get(i)[1];
			fw.write(mg.name+"\t"+l1+"\t"+(l1/l2)+"\t"+mg.geneNames.size()+"\n");
		}
		fw.close();
		}
		if(explainedVariances_randomDistributions!=null){
			FileWriter fw = new FileWriter(outputFolder+"module_scores.xls");
			System.out.println("MODULE\tL1\tL1_pv\tL1/L2\tL1/L2_pv\tNUMBER_OF_GENES");
			fw.write("MODULE\tL1\tL1_pv\tL1/L2\tL1/L2_pv\tNUMBER_OF_GENES\n");
			for(int i=0;i<tables.size();i++){
				Metagene mg = (Metagene)signatures.get(i);
				double values[] = new double[3];
				values[0] = explainedVariances.get(i)[0];
				values[1] = explainedVariances.get(i)[1];
				values[2] = values[0]/values[1];
				float pvalues[] = calcPValue4ExplainedVariance(explainedVariances_randomDistributions,mg.geneNames.size(),values);
				System.out.println(mg.name+"\t"+values[0]+"\t"+pvalues[0]+"\t"+(values[0]/values[1])+"\t"+pvalues[2]+"\t"+mg.geneNames.size());
				fw.write(mg.name+"\t"+values[0]+"\t"+pvalues[0]+"\t"+(values[0]/values[1])+"\t"+pvalues[2]+"\t"+mg.geneNames.size()+"\n");
			}
			fw.close();
		}
		
		VDatReadWrite.saveToVDatFile(moduleTable, outputFolder+"moduletable_simple.dat");
		VDatReadWrite.saveToSimpleDatFile(moduleTable, outputFolder+"moduletable_simple.txt", true);
		
		VDataTable tablescores = VDatReadWrite.LoadFromSimpleDatFile(outputFolder+"module_scores.xls", true, "\t");
		VDataTable moduleTableWithScores = VSimpleProcedures.MergeTables(moduleTable, moduleTable.fieldNames[0], tablescores, tablescores.fieldNames[0], "0");
		VDatReadWrite.saveToVDatFile(moduleTableWithScores, outputFolder+"moduletable_withscores.dat");
		VDatReadWrite.saveToSimpleDatFile(moduleTableWithScores, outputFolder+"moduletable_withscores.txt", true);
		
		VDataTable vt = VDatReadWrite.LoadFromSimpleDatFile(outputFolder+"moduletable_simple.txt", true, "\t");
		VSimpleProcedures.findAllNumericalColumns(vt);
		VSimpleFunctions.makeCorrelationTable(vt,correlationThreshold,outputFolder+"module_correlations"+correlationThreshold+".txt");


		
		VDataTable moduleTableT = moduleTable.transposeTable(geneField);
		if(sampleTable!=null)
			moduleTableT = VSimpleProcedures.MergeTables(moduleTableT, "NAME", sampleTable, sampleTable.fieldNames[0], "_");
		VDatReadWrite.saveToVDatFile(moduleTableT, outputFolder+"moduletable_simple_T.dat");
		
	}
	
	public static void calcAverageModuleActivities(String field, String normalClassLabel) throws Exception{
		FileWriter fw = new FileWriter(outputFolder+"moduleActivities.xls");
		FileWriter fwz = new FileWriter(outputFolder+"moduleActivitiesZscores.xls");
		Vector classes = new Vector();
		int fieldN = sampleTable.fieldNumByName(field);
		if(fieldN==-1){
			System.out.println("ERROR: field "+field+" is not found in "+sampleFile);
		}else{
		for(int i=0;i<moduleTable.fieldNames.length;i++){
			String sample = moduleTable.fieldNames[i];
			if(sampleTable.tableHashPrimary.containsKey(sample)){
			String cl = sampleTable.stringTable[sampleTable.tableHashPrimary.get(sample).get(0)][fieldN];
			if(cl!=null)
			if(classes.indexOf(cl)<0){
				classes.add(cl);
			}
			}
		}
		classes.add("_");
		int normalClass = classes.indexOf(normalClassLabel);		
		System.out.println("Normal class = "+normalClass);
		Vector stats = new Vector();
		VStatistics glstat = new VStatistics(moduleTable.rowCount);
		fw.write("NAME\t");
		fwz.write("NAME\t");
		for(int i=0;i<classes.size();i++){
			stats.add(new VStatistics(moduleTable.rowCount));
			fw.write((String)classes.get(i)+"\t");
		} 
		for(int i=0;i<classes.size();i++){
			fw.write((String)classes.get(i)+"_s\t");
		}		
		fw.write("\n"); fwz.write("\n");
		VDataTable transp = moduleTable.transposeTable(geneField);
		transp = VSimpleProcedures.MergeTables(transp, "NAME", sampleTable, sampleTable.fieldNames[0], "_");
		VDatReadWrite.saveToVDatFile(transp, outputFolder+"temp.dat");
		VDataSet transpd = VSimpleProcedures.SimplyPreparedDatasetWithoutNormalization(transp, -1);
		for(int i=0;i<transpd.pointCount;i++){
			float f[] = transpd.getVector(i);
			glstat.addNewPoint(f);
			String cl = transp.stringTable[i][transp.fieldNumByName(field)];
			if(cl!=null){
				int k = classes.indexOf(cl);
				if(k==-1)
					System.out.println("Lavel "+cl+" not found");
				((VStatistics)stats.get(k)).addNewPoint(f);
			}
		}
		glstat.calculate(); glstat.calcMedians();
		for(int i=0;i<classes.size();i++){
			((VStatistics)stats.get(i)).calculate();
			((VStatistics)stats.get(i)).calcMedians();
		}
		for(int i=0;i<moduleTable.rowCount;i++){
			fw.write(moduleTable.stringTable[i][moduleTable.fieldNumByName(geneField)]+"\t");
			fwz.write(moduleTable.stringTable[i][moduleTable.fieldNumByName(geneField)]+"\t");			
			float min = glstat.getMin(i);
			float max = glstat.getMax(i);
			float glmedian = glstat.getMedian(i);
			float glmean = glstat.getMean(i);
			float glstdev = glstat.getStdDev(i);
			//System.out.println(min+"\t"+max);
			float stdev = 0f;
			for(int k=0;k<classes.size();k++){
				float sd = ((VStatistics)stats.get(k)).getStdDev(i);
				stdev += sd*sd;
			}
			stdev = (float)Math.sqrt(stdev);
			if(normalClass>=0){
				glmean = ((VStatistics)stats.get(normalClass)).getMean(i);
				stdev = ((VStatistics)stats.get(normalClass)).getStdDev(i);
				if(stdev==0) stdev = 1;
				glmedian = ((VStatistics)stats.get(normalClass)).getMedian(i);
			}
			for(int k=0;k<classes.size();k++){
				//float mean = ((VStatistics)stats.get(k)).getMedian(i);
				float mean = ((VStatistics)stats.get(k)).getMean(i);
				//System.out.println((String)classes.get(k)+":"+mean+"\t"+((VStatistics)stats.get(k)).pointsNumber);				
				//float val = (mean-min)/(max-min);
				//float val = (float)(mean-glmedian)/glstdev;
				//float val = mean-glmedian;
				//float val = (float)(mean-glmean)/stdev;
				float val = (float)(mean-glmean);
				//float val = mean;
				fw.write(""+val+"\t");
			}
			for(int k=0;k<classes.size();k++){
				float sv = ((VStatistics)stats.get(k)).getStdDev(i);
				fw.write(""+sv+"\t");
			}
			
			if(stats.size()>1)
				fwz.write(""+(((VStatistics)stats.get(0)).getMean(i)-((VStatistics)stats.get(1)).getMean(i))/stdev+"\t");
			
			fw.write("\n");
			fwz.write("\n");
		}
		fw.close();
		fwz.close();
		}
	}
	
	public static void writeGeneModuleMatrix(){
		Vector allGeneNames = new Vector();
		for(int i=0;i<signatures.size();i++){
			Metagene mg = signatures.get(i);
			for(int j=0;j<mg.geneNames.size();j++)
				if(allGeneNames.indexOf((String)mg.geneNames.get(j))<0)
					allGeneNames.add((String)mg.geneNames.get(j));
		}
		Collections.sort(allGeneNames);
		try{
			FileWriter fw = new FileWriter(outputFolder+"GeneModuleMatrix.xls");
			fw.write("GENESYMBOL\t");
			for(int i=0;i<signatures.size();i++){
				Metagene mg = signatures.get(i);
				fw.write(mg.name+"\t");
			}
			fw.write("\n");
			
			for(int i=0;i<allGeneNames.size();i++){
				String gname = (String)allGeneNames.get(i);
				fw.write(gname+"\t");
				for(int j=0;j<signatures.size();j++){
					Metagene mg = signatures.get(j);
					int k = mg.geneNames.indexOf(gname); 
					if(k>=0){
						fw.write(((Float)mg.weights.get(k)).floatValue()+"\t");
					}else{
						fw.write("0\t");
					}
				}
				fw.write("\n");	
			}
			
			
			fw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public static void writeGeneProjectionMatrix(){
		Vector allGeneNames = new Vector();
		for(int i=0;i<tables.size();i++){
			VDataTable mg = tables.get(i);
			mg.makePrimaryHash(mg.fieldNames[0]);
			for(int j=0;j<mg.rowCount;j++)
				if(allGeneNames.indexOf(mg.stringTable[j][0])<0)
					allGeneNames.add(mg.stringTable[j][0]);
		}
		Collections.sort(allGeneNames);
		try{
			FileWriter fw = new FileWriter(outputFolder+"GeneProjectionMatrix.xls");
			fw.write("GENESYMBOL\t");
			for(int i=0;i<signatures.size();i++){
				Metagene mg = signatures.get(i);
				fw.write(mg.name+"\t");
			}
			fw.write("\n");
			
			for(int i=0;i<allGeneNames.size();i++){
				String gname = (String)allGeneNames.get(i);
				fw.write(gname+"\t");
				for(int j=0;j<tables.size();j++){
					VDataTable mg = tables.get(j);
					if(mg.tableHashPrimary.get(gname)==null){
						fw.write("N/A\t");
					}else{
					int k = mg.tableHashPrimary.get(gname).get(0);
					//if(mg.geneNames.size()!=Ss.get(j).length)
					//	System.out.println("ERROR: size ("+mg.geneNames.size()+") of "+mg.name+" does not correspond to the length of projection vector ("+Ss.get(j).length+")");
						fw.write(""+Ss.get(j)[k][0]+"\t");
					}
				}
				fw.write("\n");	
			}
			
			
			fw.close();
			
			VDataTable vt = VDatReadWrite.LoadFromSimpleDatFile(outputFolder+"GeneProjectionMatrix.xls", true, "\t");
			VSimpleProcedures.findAllNumericalColumns(vt);
			VDatReadWrite.saveToVDatFile(vt, outputFolder+"GeneProjectionMatrix.dat");
			
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	
	public static void assignWeights(Metagene mg, VDataSet vd, VDataTable dat){
		for(int k=0;k<mg.weights.size();k++){
			float w = ((Float)mg.weights.get(k)).floatValue();
			if(Math.abs(w-1f)>1e-6f){
				vd.weighted = true;
				//System.out.println("Weight "+mg.geneNames.get(k)+"="+w);
				break;
			}
		}
		
		if(vd.weighted){
			
			if(dat.fieldNumByName("WEIGHT")==-1)
				dat.addNewColumn("WEIGHT", "", "", dat.STRING, "1");
			
			float wSum = 0f;
			vd.weights = new float[vd.pointCount];
			for(int k=0;k<vd.pointCount;k++){
				String gname = dat.stringTable[k][dat.fieldNumByName(geneField)];
				int kk = mg.geneNames.indexOf(gname);
				float w = ((Float)mg.weights.get(kk)).floatValue();
				if(dat.fieldNumByName("WEIGHT")!=-1)
					dat.stringTable[k][dat.fieldNumByName("WEIGHT")] = ""+w;				
				if(Math.abs(w-1f)>1e-6f){
					vd.weights[k] = Math.abs(w);
				}else
					vd.weights[k] = 1;
				wSum+=w;
			}
			vd.weightSum = wSum;
		}
	}
	
	public static float calcSpecifiedWeightSum(VDataTable vt, VDataSet vdp, Metagene sig){
		float res = 0f;
		//System.out.println(sig.name);
		for(int i=0;i<vt.rowCount;i++){
			int k = sig.geneNames.indexOf(vt.stringTable[i][vt.fieldNumByName(geneField)]);
			if((Boolean)sig.weightSpecified.get(k)){
				res+=(Float)sig.weights.get(k)*vdp.massif[i][0];
				//System.out.println(vt.stringTable[i][vt.fieldNumByName("GeneSymbol")]+"\tProj="+vdp.massif[i][0]+"Weight=\t"+(Float)sig.weights.get(k)+"\tRes="+res);
			}
		}
		return res;
	}
	
	public static float calcAllWeightSum(VDataTable vt, VDataSet vdp, Metagene sig){
		float res = 0f;
		//System.out.println(sig.name);
		for(int i=0;i<vt.rowCount;i++){
			int k = sig.geneNames.indexOf(vt.stringTable[i][vt.fieldNumByName(geneField)]);
				res+=(Float)sig.weights.get(k)*vdp.massif[i][0];
				//System.out.println(vt.stringTable[i][vt.fieldNumByName("GeneSymbol")]+"\tProj="+vdp.massif[i][0]+"Weight=\t"+(Float)sig.weights.get(k)+"\tRes="+res);
		}
		return res;
	}
	
	
	  public static Vector<String> determinePC1Outliers(VDataSet vd, String probeID, double thresh){
		    Vector<String> vr = new Vector<String>();
		    return vr;
		  }
	
	  public static Vector determineOutliers(VDataSet vd, String probeID, double thresh, int comp){
		    Vector vr = new Vector();
		    double dt[][] = new double[vd.pointCount][vd.pointCount];

		    PCAMethod PCA = new PCAMethod();
		    PCA.setDataSet(vd);
		    PCA.calcBasis(comp);
		    VDataSet vdp = PCA.getProjectedDataset();

		    VStatistics sgn = new VStatistics(1);
		    float d[] = new float[1];
		    for(int i=0;i<vd.pointCount;i++)
		      for(int j=i+1;j<vd.pointCount;j++){
		         float xi[] = vdp.getVector(i);
		         float xj[] = vdp.getVector(j);
		         dt[i][j] = 0;
		         for(int k=0;k<xi.length;k++) dt[i][j]+=(xi[k]-xj[k])*(xi[k]-xj[k]);
		         dt[i][j] = Math.sqrt(dt[i][j]);
		         d[0] = (float)dt[i][j];
		         sgn.addNewPoint(d);
		      }
		    sgn.calculate();

		    int pid = vd.table.fieldNumByName(probeID);
		    for(int i=0;i<vd.pointCount;i++){
		      VStatistics vst = new VStatistics(1);
		      for(int j=0;j<vd.pointCount;j++) if(i!=j){
		        d[0] = (float)(dt[i][j]+dt[j][i]);
		        vst.addNewPoint(d);
		      }
		      vst.calculate();
		      float z = Math.abs((vst.getMean(0)-sgn.getMean(0))/sgn.getStdDev(0));
		      float x[] = vdp.getVector(i);
		      float s = 0f;
		      for(int k=1;k<x.length;k++) s+=Math.abs(x[k]);
		      //z = z*s/Math.abs(x[0]);
		      if(z>thresh){
		        vr.add(vd.table.stringTable[i][pid]);
		        //System.out.println(vd.table.stringTable[i][pid]+"\t"+z);
		      }
		    }
		    return vr;
		  }

	  public static VDataTable removeOutliers(VDataTable rtab, Metagene mg, boolean verbose){
			VDataSet ds = VSimpleProcedures.SimplyPreparedDataset(rtab,-1);
		    //Vector<String> v = determineOutliers(ds,geneField,outlierThreshold,outlierDimension);
			Vector<String> v = leaveOneOutPCAOutliers(ds,geneField,outlierThreshold);
		    if(v.size()>0)if(verbose)
		    	System.out.print("outliers: ");
		    for(int jj=0;jj<v.size();jj++){
		    		if(verbose)
		    			System.out.print(v.get(jj)+"\t");
		    		int k = mg.geneNames.indexOf((String)v.get(jj));
		        	mg.geneNames.remove((String)v.get(jj));
		        	mg.weights.remove(k);
		        	mg.weightSpecified.remove(k);
		    }
		    if(v.size()>0)if(verbose)
		    	System.out.println();
		    VDataTable rtab1 = VSimpleProcedures.selectRowsFromList(table, mg.geneNames,false);
		    return rtab1;
	  }
	  
	  public static Vector<String> leaveOneOutPCAOutliers(VDataSet ds, String idfield, float threshold){
		  Vector<String> outliers = new Vector<String>();
		  ds.weighted = true;
		  ds.weights = new float[ds.pointCount];

		  VStatistics stat = new VStatistics(1);
		  float f[] = new float[1];
		  float vars[] = new float[ds.pointCount];
		  
		  for(int k=0;k<ds.pointCount;k++){
			  for(int i=0;i<ds.pointCount;i++) ds.weights[i] = 1f;
			  ds.weights[k] = 0f;
				PCAMethod pca = new PCAMethod();
				pca.setDataSet(ds);
				pca.calcBasis(1);
				ds.calcStatistics();				
				double explained_var[] = pca.calcDispersions();
				f[0] = (float)explained_var[0];
				vars[k] = f[0];
				stat.addNewPoint(f);
		  }
		  stat.calculate();
		  for(int k=0;k<ds.pointCount;k++){
			  float z = Math.abs((vars[k]-stat.getMean(0))/stat.getStdDev(0));
			  if(z>threshold)
				  outliers.add(ds.table.stringTable[k][ds.table.fieldNumByName(idfield)]);
		  }
		  
		  return outliers;
	  }
	  
	  public static void ProduceGraphicalOutput() throws Exception{
		  
		  FileWriter fw = new FileWriter(outputFolder+"allnames.txt");
		  for(int i=0;i<table.rowCount;i++)
			  fw.write(table.stringTable[i][0]+"\n");
		  fw.close();
		  
		  for(int i=0;i<tables.size();i++){
			  String name = signatures.get(i).name;
			  tables.get(i).makePrimaryHash(tables.get(i).fieldNames[0]);
			  
				Metagene mg = (Metagene)signatures.get(i);
				double values[] = new double[3];
				values[0] = explainedVariances.get(i)[0];
				values[1] = explainedVariances.get(i)[1];
				values[2] = values[0]/values[1];
				
				float pvalues[] = new float[2];
				if(explainedVariances_randomDistributions!=null)
					pvalues = calcPValue4ExplainedVariance(explainedVariances_randomDistributions,mg.geneNames.size(),values);
			  
			  if(saveDecomposedFiles)if((pvalues[0]<=graphicalOutputThreshold)||(pvalues[2]<=graphicalOutputThreshold)){
				  
			  int density[] = analyseGlobalProjectionDensity(globalProjections.get(i),10,0.15f,(int)(0.1f*table.rowCount));
				  
			  FileWriter fw1 = new FileWriter(outputFolder+name+"_projs.txt");
			  FileWriter fwst = new FileWriter(outputFolder+name+"_stats.txt");
			  
			  for(int k=0;k<table.rowCount;k++){
				  fw1.write(globalProjections.get(i)[k][0]+"\t"+globalProjections.get(i)[k][1]+"\t");
				  String gname = table.stringTable[k][0];
				  if(tables.get(i).tableHashPrimary.get(gname)!=null){
					  int ind = tables.get(i).tableHashPrimary.get(gname).get(0);
					  fw1.write("1\t");
					  //fw1.write("\t"+Ss.get(i)[ind][0]+"\t"+Ss.get(i)[ind][1]+"\t"+gname);
				  }else{
					  fw1.write("0\t");
				  }
				  fw1.write(""+density[k]);
				  fw1.write("\n");
			  }
			  fw1.close();
			  fwst.close();}
		  }
		  
	  }
	  
	  public static int[] analyseGlobalProjectionDensity(float points2D[][], int numberOfLevels, float windowWidth, int reducedSampleSize){
		  int levels[] = new int[points2D.length];
		  double dens[] = new double[points2D.length];
		  VStatistics vst = new VStatistics(2);
		  for(int i=0;i<points2D.length;i++)
			  vst.addNewPoint(points2D[i]);
		  vst.calculate();
		  float stx = vst.getStdDev(0);
		  float sty = vst.getStdDev(1);
		  
		  double minval = Double.MAX_VALUE;
		  double maxval = Double.MIN_VALUE;
		  
		  int inds[] = new int[reducedSampleSize];
		  Random r = new Random();
		  HashSet<Integer> set = new HashSet<Integer>(); 
		  for(int i=0;i<reducedSampleSize;i++){
			  int k = r.nextInt(points2D.length);
			  if(!set.contains(k))
				  set.add(k);
		  }
		  
		  for(int i=0;i<points2D.length;i++){
			  float xi = points2D[i][0];
			  float yi = points2D[i][1];
			  double dist = xi*xi/stx/stx+yi/sty*yi/sty;
			  dens[i] = -Math.sqrt(dist);
			  if(dens[i]>maxval) maxval = dens[i];
			  if(dens[i]<minval) minval = dens[i];
		  }
		  
		  /*for(int i=0;i<points2D.length;i++){
			  float xi = points2D[i][0];
			  float yi = points2D[i][1];
			  
			  //for(int j=0;j<points2D.length;j++)if(i!=j){
			  for(int j=0;j<inds.length;j++)if(inds[j]!=i){
				  float xj = points2D[inds[j]][0];
				  float yj = points2D[inds[j]][1];
				  double dx = (xi-xj)*(xi-xj)/stx/stx;
				  double dy = (yi-yj)*(yi-yj)/sty/sty;
				  dens[i]+=Math.exp(-dx-dy);
			  }
			  dens[i] = Math.log(dens[i]);
			  if(dens[i]>maxval) maxval = dens[i];
			  if(dens[i]<minval) minval = dens[i];
		  }*/
		  
		  for(int i=0;i<dens.length;i++){
			  levels[i] = (int)((dens[i]-minval)/(maxval-minval)*numberOfLevels);
		  }
		  
		  return levels;
	  }

}
