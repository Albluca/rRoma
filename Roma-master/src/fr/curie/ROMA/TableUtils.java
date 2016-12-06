package fr.curie.ROMA;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

import vdaoengine.data.*;
import vdaoengine.analysis.*;
import vdaoengine.utils.*;
import vdaoengine.data.io.*;


public class TableUtils {

  public static void main(String[] args) {
	  
	  try{
		  
			VDataTable vtt = VDatReadWrite.LoadFromVDatFile("C:/Datas/NotchP53/data/tcga/filtered_vida_parsed_crc_220_microarray.dat");
			VDataTable vt1 = TableUtils.fillMissingValues(vtt, 10);
			VDatReadWrite.saveToVDatFile(vt1, "C:/Datas/NotchP53/data/tcga/filtered_restored.dat");
			System.exit(-1);
		  
		  
	  VDataTable vt = vdaoengine.data.io.VDatReadWrite.LoadFromSimpleDatFile("C:/Datas/BreastCancer/IVOIRE/normDataRMA.txt",true,"\t");
	  VDataTable annot = VDatReadWrite.LoadFromSimpleDatFile("c:/datas/moduleactivities/annot4.txt", true, "\t");
      for(int i=1;i<vt.colCount;i++)
    	  vt.fieldTypes[i] = vt.NUMERICAL;
      vt = TableUtils.normalizeVDat(vt, true, false);
      VDataTable mt = vdaoengine.utils.VSimpleProcedures.MergeTables(vt,"Probe",annot,"Probeset","0");
      VDatReadWrite.saveToSimpleDatFile(mt, "c:/datas/BreastCancer/IVOIRE/normDataRMAn.txt");
      
		  
	  /*VDataTable vt = vdaoengine.data.io.VDatReadWrite.LoadFromVDatFile("C:/Datas/EWING/Kinetics/kinetic_data/geneviewer/tsc1n.dat");
	  //vt = TableUtils.normalizeVDat(vt, true, false);
	  //vt = TableUtils.filterByVariation(vt, 1000, false);
	  Vector ids = Utils.loadStringListFromFile("C:/Datas/EWING/Kinetics/kinetic_data/618best");
	  vt = VSimpleProcedures.selectRowsFromList(vt, ids, "CHIP");
	  vdaoengine.data.io.VDatReadWrite.saveToVDatFile(vt, "C:/Datas/EWING/Kinetics/kinetic_data/geneviewer/tsc1n_618sel.dat");*/
	  }catch(Exception e){
		  e.printStackTrace();
	  }
	  
  }

  public static void printFieldInfoSummary(VDataTable vt, int fieldClass, int fieldSubclass){
    Vector cl = new Vector();
    Vector subcl = new Vector();
    for(int i=0;i<vt.colCount;i++){
      String scl = vt.fieldInfo[i][fieldClass];
      String sscl = vt.fieldInfo[i][fieldSubclass];
      if((scl!=null)&&(!scl.equals("")))
        if(cl.indexOf(scl)<0)
          cl.add(scl);
      if((sscl!=null)&&(!sscl.equals("")))
        if(subcl.indexOf(sscl)<0)
          subcl.add(sscl);
    }
    float table[][] =  new float[cl.size()][subcl.size()];
    float tablea[][] =  new float[cl.size()][subcl.size()];
    float totalsubcl[] =  new float[subcl.size()];
    for(int i=0;i<vt.colCount;i++){
      String scl = vt.fieldInfo[i][fieldClass];
      String sscl = vt.fieldInfo[i][fieldSubclass];
      if((scl!=null)&&(!scl.equals("")))
        if((sscl!=null)&&(!sscl.equals(""))){
             int i1 = cl.indexOf(scl);
             int j1 = subcl.indexOf(sscl);
             table[i1][j1]+=1f;
             tablea[i1][j1]+=1f;
             totalsubcl[j1]+=1f;
        }
    }
    System.out.print("Class");
    for(int j=0;j<subcl.size();j++)
      System.out.print("\t"+(String)subcl.elementAt(j)+"\t"+(String)subcl.elementAt(j));
    System.out.println();
    for(int i=0;i<cl.size();i++){
      System.out.print((String)cl.elementAt(i));
      float sum = 0f;
      for(int j=0;j<subcl.size();j++)
        sum+=table[i][j];
      for(int j=0;j<subcl.size();j++){
        table[i][j]/=sum;
        System.out.print("\t"+tablea[i][j]+"\t"+table[i][j]);
      }
    System.out.println();
    }
    System.out.print("Total");
    float sum = 0f;
    for(int j=0;j<subcl.size();j++)
      sum+=totalsubcl[j];
    for(int j=0;j<subcl.size();j++){
      System.out.print("\t"+totalsubcl[j]);
      totalsubcl[j]/=sum;
      System.out.print("\t"+totalsubcl[j]);
    }
    System.out.println();
  }

  public static void printClassSpecificSignaturePrediction(VDataTable vt, String fieldClass, String fieldPrediction, String fieldAnswer, String rightAnswer){
    Vector cl = new Vector();
    float thresh = 0.1f;
    int idfc = vt.fieldNumByName(fieldClass);
    int idfp = vt.fieldNumByName(fieldPrediction);
    int idfa = vt.fieldNumByName(fieldAnswer);
    for(int i=0;i<vt.rowCount;i++){
      String scl = vt.stringTable[i][idfc];
      if((scl!=null)&&(!scl.equals("")))
        if(cl.indexOf(scl)<0)
          cl.add(scl);
    }
    float TP[] = new float[cl.size()];
    float TN[] = new float[cl.size()];
    float FP[] = new float[cl.size()];
    float FN[] = new float[cl.size()];
    float TPt = 0f, TNt = 0f, FPt = 0f, FNt = 0f;
    float TPr = 0f, TNr = 0f, FPr = 0f, FNr = 0f;
    float ratio = 0f; float total = 0f;

    for(int i=0;i<vt.rowCount;i++){
     String scl = vt.stringTable[i][idfc];
     float pred = Float.parseFloat(vt.stringTable[i][idfp]);
     boolean predicted = pred>thresh;
     boolean answered = vt.stringTable[i][idfa].equals(rightAnswer);
     if((predicted)&&(answered)) { TP[cl.indexOf(scl)]+=1f; TPt+=1f; }
     if((predicted)&&(!answered)) { FP[cl.indexOf(scl)]+=1f; FPt+=1f; }
     if((!predicted)&&(answered)) { FN[cl.indexOf(scl)]+=1f; FNt+=1f; }
     if((!predicted)&&(!answered)) { TN[cl.indexOf(scl)]+=1f; TNt+=1f; }
     if(answered) ratio+=1f;
     total+=1f;
    }
    ratio/=total;

    System.out.println("Class\tSn\tSp\tPPV\tNPV\tAcc\tTP\tFP\tTN\tFN");
    for(int i=0;i<cl.size();i++){
      System.out.print((String)cl.elementAt(i));
      System.out.print("\t"+(TP[i]/(TP[i]+FN[i]))+"\t"+(TN[i]/(FP[i]+TN[i])));
      System.out.print("\t"+(TP[i]/(TP[i]+FP[i]))+"\t"+(TN[i]/(FN[i]+TN[i])));
      System.out.print("\t"+((TP[i]+TN[i])/(TP[i]+FP[i]+TN[i]+FN[i])));
      System.out.println("\t"+TP[i]+"\t"+FP[i]+"\t"+TN[i]+"\t"+FN[i]);
    }
    System.out.print("Total");
    System.out.print("\t"+(TPt/(TPt+FNt))+"\t"+(TNt/(FPt+TNt)));
    System.out.print("\t"+(TPt/(TPt+FPt))+"\t"+(TNt/(FNt+TNt)));
    System.out.print("\t"+((TPt+TNt)/(TPt+FPt+TNt+FNt)));
    System.out.println("\t"+TPt+"\t"+FPt+"\t"+TNt+"\t"+FNt);

    Random r = new Random();
    for(int k=0;k<1000;k++){
      for(int i=0;i<vt.rowCount;i++){
        boolean predicted = r.nextFloat()>(1-ratio);
        boolean answered = vt.stringTable[i][idfa].equals(rightAnswer);
        if((predicted)&&(answered)) { TPr+=1f; }
        if((predicted)&&(!answered)) { FPr+=1f; }
        if((!predicted)&&(answered)) { FNr+=1f; }
        if((!predicted)&&(!answered)) { TNr+=1f; }
      }
    }
    System.out.print("Random_guess");
    System.out.print("\t"+(TPr/(TPr+FNr))+"\t"+(TNr/(FPr+TNr)));
    System.out.print("\t"+(TPr/(TPr+FPr))+"\t"+(TNr/(FNr+TNr)));
    System.out.print("\t"+((TPr+TNr)/(TPr+FPr+TNr+FNr)));
    System.out.println("\t"+TPr+"\t"+FPr+"\t"+TNr+"\t"+FNr);

  }

  public static VDataTable filterByVariation(VDataTable vt, int numOfGenes, boolean doScaling){
    float var[] = new float[vt.rowCount];
    VDataTable vtr = new VDataTable();
    vtr.copyHeader(vt);
    vtr.rowCount = numOfGenes;
    vtr.colCount = vt.colCount;
    vtr.stringTable = new String[numOfGenes][vt.colCount];
    VDataSet vd = null;
    if(doScaling)
      vd = vdaoengine.utils.VSimpleProcedures.SimplyPreparedDataset(vt,-1);
    else
      vd = vdaoengine.utils.VSimpleProcedures.SimplyPreparedDatasetWithoutNormalization(vt,-1);
    for(int i=0;i<vt.rowCount;i++){
      float v[] = vd.getVector(i);
      var[i] = (float)norm(v);
    }
    int ord[] = Sort(var);
    for(int i=0;i<numOfGenes;i++){
      int k = ord[i];
      vtr.stringTable[i] = vt.stringTable[k];
    }
    return vtr;
  }
  
  public static VDataTable filterByAverageValue(VDataTable vt, float threshold){
	  
	  	VDataSet vd = VSimpleProcedures.SimplyPreparedDatasetWithoutNormalization(vt.transposeTable(vt.fieldNames[0]), -1);
	  	vd.calcStatistics();
	  	vd.simpleStatistics.calcMedians();
	  	int numOfGenes = 0;
	  	for(int i=0;i<vd.simpleStatistics.medians.length;i++){
	  		if(vd.simpleStatistics.medians[i]>threshold)
	  			numOfGenes++;
	  	}
	  	System.out.println("Number of genes = "+numOfGenes);
	  
	    float var[] = new float[vt.rowCount];
	    VDataTable vtr = new VDataTable();
	    vtr.copyHeader(vt);
	    vtr.rowCount = numOfGenes;
	    vtr.colCount = vt.colCount;
	    vtr.stringTable = new String[numOfGenes][vt.colCount];

	    int k=0;
	    for(int i=0;i<vt.rowCount;i++){
	    	if(vd.simpleStatistics.medians[i]>threshold){
	    	  vtr.stringTable[k] = vt.stringTable[i];
	    	  k++;
	      }
	    }
	    
	    return vtr;
	  }
  
  public static VDataTable PCAtable(VDataTable vt, boolean doScaling){
	  return PCAtable(vt, doScaling, 3);
  }

  public static VDataTable PCAtable(VDataTable vt, boolean doScaling, int numOfComponents){
    VDataTable vtr = new VDataTable();
    int numOfStringFields = 0;
    for(int i=0;i<vt.colCount;i++)
      if(vt.fieldTypes[i]==vt.STRING)
        numOfStringFields++;
    vtr.fieldClasses = new String[numOfStringFields+numOfComponents];
    vtr.fieldDescriptions = new String[numOfStringFields+numOfComponents];
    vtr.fieldNames = new String[numOfStringFields+numOfComponents];
    vtr.fieldTypes = new int[numOfStringFields+numOfComponents];
    if(vt.fieldInfo!=null)
      vtr.fieldInfo = new String[numOfStringFields+numOfComponents][vt.fieldInfo[0].length];
    int k=0;
    for(int i=0;i<vt.colCount;i++)
      if(vt.fieldTypes[i]==vt.STRING){
        vtr.fieldClasses[k] = vt.fieldClasses[i];
        vtr.fieldNames[k] = vt.fieldNames[i];
        vtr.fieldDescriptions[k] = vt.fieldDescriptions[i];
        vtr.fieldTypes[k] = vt.fieldTypes[i];
        if(vt.fieldInfo!=null)
          vtr.fieldInfo[k] = vt.fieldInfo[i];
        k++;
      }
    for(int i=0;i<numOfComponents;i++){
      vtr.fieldNames[i+numOfStringFields] = "PC"+(i+1);
      vtr.fieldTypes[i+numOfStringFields] = vt.NUMERICAL;
    }
    vtr.rowCount = vt.rowCount;
    vtr.colCount = numOfStringFields+numOfComponents;
    vtr.stringTable = new String[vtr.rowCount][vtr.colCount];
    for(int i=0;i<vtr.rowCount;i++){
      k = 0;
      for(int j=0;j<vt.colCount;j++){
        if(vt.fieldTypes[j]==vt.STRING){
            vtr.stringTable[i][k] = vt.stringTable[i][j];
            k++;
        }
      }
    }

    VDataSet vd = null;
      if(doScaling)
        vd = vdaoengine.utils.VSimpleProcedures.SimplyPreparedDataset(vt,-1);
      else
        vd = vdaoengine.utils.VSimpleProcedures.SimplyPreparedDatasetWithoutNormalization(vt,-1);
    PCAMethod pca = new PCAMethod();
    pca.setDataSet(vd);
    pca.calcBasis(numOfComponents);
    VDataSet vdp = pca.getProjectedDataset();
    for(int i=0;i<vtr.rowCount;i++)
    for(int j=0;j<numOfComponents;j++){
      vtr.stringTable[i][j+numOfStringFields] = "" + vdp.massif[i][j];
    }
    return vtr;
  }

  public static double norm( float[] data )
      {
       double d = 0;
       for(int i=0;i<data.length;i++)
           d+=data[i]*data[i];
       return Math.sqrt(d);
      }

      public static double mean( float[] data )
          {
           double d = 0;
           for(int i=0;i<data.length;i++)
               d+=data[i];
           return d/data.length;
          }

          public static double stddev( float[] data )
              {
               double d = 0;
               float mn = (float)mean(data);
               for(int i=0;i<data.length;i++)
                   d+=(data[i]-mn)*(data[i]-mn);
               d/=data.length;
               return Math.sqrt(d);
              }


   public static int[] Sort(float cais[]){
    int res[]=new int[cais.length];
    for (int i = 0; i < res.length; i++) res[i]=i;

    int i,j,k,inc,n=cais.length;
    float v;

    inc=1;
    do {
            inc *= 3;
            inc++;
    } while (inc <= n);

    do {
            inc /= 3;
            for (i=inc+1;i<=n;i++) {
                    v=cais[res[i-1]];
                    j=i;
                    k=res[i-1];
                    while (cais[res[j-inc-1]]<v) {
                            //cais[j]=cais[j-inc];
                            res[j-1]=res[j-inc-1];
                            j -= inc;
                            if (j <= inc) break;
                    }
                    //cais[j]=v;
                    res[j-1]=k;
            }
    } while (inc > 0);

    return res;
    }


    public static VDataTable substituteRowsByStatistics(VDataTable vt, String classField, int statType){ // 0 - mean, 1 - stdv, 2 - min, 3 - max, 4 - median
      VDataTable vres = new VDataTable();
      Vector classes = new Vector();
      for(int i=0;i<vt.rowCount;i++){
        String cl = vt.stringTable[i][vt.fieldNumByName(classField)];
        if(classes.indexOf(cl)<0)
          classes.add(cl);
      }

      VDataSet vd = vdaoengine.utils.VSimpleProcedures.SimplyPreparedDatasetWithoutNormalization(vt,-1);
      Vector stats = new Vector();
      for(int i=0;i<classes.size();i++){
        VStatistics vst = new VStatistics(vd.coordCount);
        stats.add(vst);
      }
      for(int i=0;i<vd.pointCount;i++){
        float f[] = vd.getVector(i);
        String cl = vt.stringTable[i][vt.fieldNumByName(classField)];
        int k = classes.indexOf(cl);
        VStatistics vs = (VStatistics)stats.elementAt(k);
        vs.addNewPoint(f);
      }
      for(int i=0;i<stats.size();i++){
        VStatistics vs = (VStatistics)stats.elementAt(i);
        vs.calculate();
      }

      vres.copyHeader(vt);
      vres.rowCount = classes.size();
      vres.stringTable = new String[vres.rowCount][vres.colCount];
      String pref = "";
      if(statType==0) pref = "MEAN";
      if(statType==1) pref = "STV";
      if(statType==2) pref = "MIN";
      if(statType==3) pref = "MAX";
      for(int i=0;i<vres.rowCount;i++){
        vres.stringTable[i][vres.fieldNumByName(classField)] = (String)classes.elementAt(i)+"_"+pref;
        VStatistics vs = (VStatistics)stats.elementAt(i);
        for(int j=0;j<vd.selector.selectedColumns.length;j++){
          float num = 0f;
          if(statType==0) num = vs.getMean(j);
          if(statType==1) num = vs.getStdDev(j);
          if(statType==2) num = vs.getMin(j);
          if(statType==3) num = vs.getMax(j);
          vres.stringTable[i][vd.selector.selectedColumns[j]] = ""+num;
        }
      }

      return vres;
    }
    
    

    public static VDataTable normalizeVDat(VDataTable vt, boolean centr, boolean norm) throws Exception{
    	return normalizeVDat(vt, centr, norm, false, false, false);
    }
    
    public static VDataTable normalizeVDat(VDataTable vt, boolean centr, boolean norm, boolean useMedian, boolean addColumnWithMeanValue, boolean addColumnStdDev) throws Exception{

      boolean num[] = new boolean[vt.colCount];
      int numfields = 0;
      for(int i=0;i<vt.colCount;i++){
        if(vt.fieldTypes[i]==vt.NUMERICAL) {numfields++; num[i]=true; }
        else num[i]=false;
      }

      VDataTable res = new VDataTable();
      res.copyHeader(vt);
      res.rowCount = vt.rowCount;
      res.stringTable = new String[res.rowCount][res.colCount];
      if(addColumnWithMeanValue){
    	  res.addNewColumn(useMedian?"MEDIAN":"MEAN", "", "", res.NUMERICAL, "0");
      }
      if(addColumnStdDev){
    	  res.addNewColumn("STDEV", "", "", res.NUMERICAL, "0");
      }

      
      VStatistics vs = new VStatistics(vt.rowCount);
      for(int j=0;j<vt.colCount;j++){
        float f[] = new float[vt.rowCount];
        if(num[j]){
          for(int i=0;i<vt.rowCount;i++){
             f[i] = Float.parseFloat(vt.stringTable[i][j]);
          }
          vs.addNewPoint(f);
        }
      }
      vs.calculate();
      if(useMedian)
    	  vs.calcMedians();
      

      for(int i=0;i<vt.rowCount;i++){
        float f[] = new float[numfields];
        int k=0;
        for(int j=0;j<vt.colCount;j++){
          if(num[j]){
          float x = Float.parseFloat(vt.stringTable[i][j]);
          //f[j-1] = (x-vs.getMean(i))/vs.getStdDev(i);
          if(centr){
        	  if(useMedian)
        		  f[k] = (x-vs.getMedian(i));
        	  else
        		  f[k] = (x-vs.getMean(i));
          }
          else
            f[k] = x;
          if(norm)
            f[k] = (x-vs.getMean(i))/vs.getStdDev(i);
          k++;
          }
        }
        k = 0;
        for(int j=0;j<vt.colCount;j++){
          if(num[j]) { res.stringTable[i][j] = ""+f[k];  k++; }
          else res.stringTable[i][j] = vt.stringTable[i][j];
        }
        if(addColumnWithMeanValue){
      	  res.stringTable[i][res.fieldNumByName(useMedian?"MEDIAN":"MEAN")] = ""+(useMedian?vs.getMedian(i):vs.getMean(i));
        }
        if(addColumnStdDev){
      	  res.stringTable[i][res.fieldNumByName("STDEV")] = ""+vs.getStdDev(i);
        }
        
      }
      return res;
    }

    // strNames - selected fields to test for associations 
    public static void printAssosiationTable(VDataTable vt, String fout, Vector<String> strNames, float thresh, int minNumberOfDistinctValuesInNumericals, boolean transposed){
      Vector<String> fieldNames = new Vector<String>();
      Vector<String[]> fieldClasses = new Vector<String[]>();
      Vector<String> valNames = new Vector<String>();
      Vector<float[]> vals = new Vector<float[]>();
      
      int numOfDistinctValues[] = TableUtils.countNumberOfDistinctValuesInColumns(vt);
      for(int i=0;i<vt.colCount;i++)
    	  System.out.println(vt.fieldNames[i]+"\t"+numOfDistinctValues[i]);
      
      // Make list of all found string columns or those numericals which contains small number of labels (categorical)
      for(int i=0;i<vt.colCount;i++){
    	  if((vt.fieldTypes[i]==vt.STRING)||(numOfDistinctValues[i]<minNumberOfDistinctValuesInNumericals)){
          fieldNames.add((String)vt.fieldNames[i]);
          String cl[] = new String[vt.rowCount];
          for(int j=0;j<vt.rowCount;j++)
            cl[j] = vt.stringTable[j][i];
          fieldClasses.add(cl);
        }
      }

      // Make list of all found numerical values
      for(int i=0;i<vt.colCount;i++){
        if((vt.fieldTypes[i]==vt.NUMERICAL)&&(numOfDistinctValues[i]>=minNumberOfDistinctValuesInNumericals)){
          valNames.add(vt.fieldNames[i]);
          float f[] = new float[vt.rowCount];
          for(int j=0;j<vt.rowCount;j++){
        	if(vt.stringTable[j][i].equals("_")||vt.stringTable[j][i].equals("@")||vt.stringTable[j][i].equals("")||vt.stringTable[j][i].equals("NA")||vt.stringTable[j][i].equals("\"\""))
        		f[j] = Float.NaN;
        	else
        		f[j] = Float.parseFloat(vt.stringTable[j][i]);
          }
          vals.add(f);
        }
      }

      float f[][] = findAssosiations(strNames, fieldNames, fieldClasses, valNames, vals);
      if(!transposed)
        printAssociations(f,strNames, fieldNames,valNames,fout,thresh);
      else
        printAssociationsT(f,strNames, fieldNames,valNames,fout,thresh);

    }

    public static float[][] findAssosiationsInDataTable(VDataTable vt, Vector strNames){
      Vector fieldNames = new Vector();
      Vector fieldClasses = new Vector();
      Vector valNames = new Vector();
      Vector vals = new Vector();
      Vector<String> fieldsToTest = strNames;
      
      for(int i=0;i<strNames.size();i++){
        int id = vt.fieldNumByName((String)strNames.elementAt(i));
        if(id!=-1){
          fieldNames.add((String)strNames.elementAt(i));
          String cl[] = new String[vt.rowCount];
          for(int j=0;j<vt.rowCount;j++)
            cl[j] = vt.stringTable[j][id];
          fieldClasses.add(cl);
        }
      }

      for(int i=0;i<vt.colCount;i++){
        if(vt.fieldTypes[i]==vt.NUMERICAL){
          valNames.add(vt.fieldNames[i]);
          float f[] = new float[vt.rowCount];
          for(int j=0;j<vt.rowCount;j++)
            f[j] = Float.parseFloat(vt.stringTable[j][i]);
          vals.add(f);
        }
      }

      return findAssosiations(fieldsToTest, fieldNames, fieldClasses, valNames, vals);
    }

    public static float[][] findAssosiations(Vector<String> fieldsToTest, Vector<String> fieldNames, Vector<String[]> fieldClasses, Vector<String> valNames, Vector<float[]> vals){
    	
        float res[][] = new float[fieldsToTest.size()][fieldNames.size()+valNames.size()];
        
        // For now we assume that all fieldsToTest are numerical, extend with Fisher test in the future
        for(int kk=0;kk<fieldsToTest.size();kk++){

        float val[] = (float[])vals.get(valNames.indexOf(fieldsToTest.get(kk)));
        
        // Associations with categorical values
        for(int i=0;i<fieldClasses.size();i++){
          String cl[] = (String[])fieldClasses.elementAt(i);
          Vector<String> lbls = new Vector<String>();
          for(int j=0;j<cl.length;j++){
            String lb = cl[j];
            if((!lb.equals("_"))&&(!lb.equals("NA"))&&(!lb.equals(""))&&(!lb.equals("\"\""))){
              if(lbls.indexOf(lb)<0) lbls.add(lb);
            }
          }
          
          int countt = (int)(0.5f*(lbls.size()-1)*lbls.size());
          //System.out.println(fieldNames.get(i)+"\t"+lbls.size());
          if(countt<50){
          float tvalues[] = new float[countt];
          int k = 0;
          double maxVal = 0f;
          String compMax1 = "";
          String compMax2 = "";
          for(int k1=0;k1<lbls.size();k1++)
            for(int k2=k1+1;k2<lbls.size();k2++){
               String lb1 = (String)lbls.elementAt(k1);
               String lb2 = (String)lbls.elementAt(k2);
                 Vector set1 = new Vector();
                 Vector set2 = new Vector();
                 for(int jj=0;jj<cl.length;jj++){
                   if(cl[jj].equals(lb1)) if(!Float.isNaN(val[jj])) set1.add(new Float(val[jj]));
                   if(cl[jj].equals(lb2)) if(!Float.isNaN(val[jj])) set2.add(new Float(val[jj]));
                 }
               double tvalue = calcTTest(set1,set2);
               if(tvalue>maxVal){
            	   maxVal = tvalue;
               }
               tvalues[k++] = (float)Math.abs(tvalue);
            }
          float tval = (float)max(tvalues);
          res[kk][i] = tval;
          
          // too many comparisons, does not make sense
          }else{
        	  res[kk][i] = Float.NaN; 
          }
        	  
          //descriptions.add(compMax1+"_vs_"+compMax2);
        }
        
        //Associations with numerical values (by Spearmann correlation)
        for(int i=0;i<valNames.size();i++){
        	float val1[] = (float[])vals.get(i);
        	float corr = VSimpleFunctions.calcSpearmanCorrelationCoeffMissingValues(val, val1);
        	res[kk][fieldClasses.size()+i] = corr;
        }
        
        }
        
        return res;
    }

    public static void printAssociations(float vals[][], Vector<String> fieldsToTest, Vector fieldNames, Vector valNames, String fout, float thresh){
      try{
        FileWriter fw = new FileWriter(fout);
        fw.write("FIELD");
        for(int i=0;i<fieldNames.size();i++){
          fw.write((String)fieldNames.elementAt(i)+"\t");
        }
        for(int i=0;i<valNames.size();i++)
            fw.write((String)valNames.elementAt(i)+"\t");
        fw.write("\n");
        
        for(int i=0;i<fieldsToTest.size();i++){
        	fw.write(fieldsToTest.get(i)+"\t");
            for(int j=0;j<fieldNames.size();j++){
           	 float f = vals[i][j];
        	 DecimalFormat df = new DecimalFormat("#.##");
        	 String sf = df.format(f); 
           	 if(f>=thresh)
           		 fw.write(sf+"\t");
           	 else
           		 fw.write("_\t");
            }
             for(int j=0;j<valNames.size();j++){
            	 float f = vals[i][j+fieldNames.size()];
            	 DecimalFormat df = new DecimalFormat("#.##");
            	 String sf = df.format(f); 
            	 fw.write(sf+"\t");
             }
        fw.write("\n");
        }
        fw.close();
      }catch(Exception e){
        e.printStackTrace();
      }
    }

    public static void printAssociationsT(float vals[][], Vector<String> fieldsToTest, Vector fieldNames, Vector valNames, String fout, float thresh){
      try{
        FileWriter fw = new FileWriter(fout);
        
        fw.write("VAL\t");
        for(int i=0;i<fieldsToTest.size();i++)
          fw.write((String)fieldsToTest.elementAt(i)+"\t");
        fw.write("\n");
        
        for(int i=0;i<fieldNames.size();i++){
        	fw.write(fieldNames.get(i)+"\t");
        	for(int j=0;j<fieldsToTest.size();j++){
        		float f = vals[j][i];
           	 	DecimalFormat df = new DecimalFormat("#.##");
           	 	String sf = df.format(f); 
                if(f>=thresh)
                    fw.write(sf+"\t");
                  else
                    fw.write("_\t");
        	}
            fw.write("\n");
        }
        
        fw.write("\n");
        
        for(int i=0;i<valNames.size();i++){
        	fw.write((String)valNames.elementAt(i)+"\t");
        	for(int j=0;j<fieldsToTest.size();j++){
        		float f = vals[j][i+fieldNames.size()];
       	 		DecimalFormat df = new DecimalFormat("#.##");
       	 		String sf = df.format(f);
       	 		fw.write(sf+"\t");
        	}
        	fw.write("\n");
        }
             
        fw.close();
      }catch(Exception e){
        e.printStackTrace();
      }
    }


    public static double calcTTest(Vector set1, Vector set2){
      double r = 0;
      VStatistics stat1 = new VStatistics(1);
      VStatistics stat2 = new VStatistics(1);
      float d[] = new float[1];
      for(int i=0;i<set1.size();i++){
        d[0] = ((Float)set1.elementAt(i)).floatValue();
        stat1.addNewPoint(d);
      }
      for(int i=0;i<set2.size();i++){
        d[0] = ((Float)set2.elementAt(i)).floatValue();
        stat2.addNewPoint(d);
      }
      stat1.calculate();
      stat2.calculate();
      r = (stat1.getMean(0)-stat2.getMean(0))/Math.sqrt(stat1.getStdDev(0)*stat1.getStdDev(0)/set1.size()+stat2.getStdDev(0)*stat2.getStdDev(0)/set2.size());
      return r;
    }

    public static float max(float[] t) {
      float maximum = 0;
        if(t.length!=0)
        {
        maximum = t[0];   // start with the first value
        for (int i=1; i<t.length; i++) {
            if (t[i] > maximum) {
                maximum = t[i];   // new maximum
            }
        }
        }
        return maximum;
    }

    public static void reformatForEisenCluster(VDataTable vt, String labelField, String fn, boolean conserveOrder) throws Exception{
      FileWriter fw = new FileWriter(fn);
      fw.write("YORF\tNAME\tGWEIGHT\t");
      if(conserveOrder)
        fw.write("GORDER\t");
      int k = vt.fieldNumByName(labelField);
      int lastcolumn = -1;
      for(int i=0;i<vt.colCount;i++)
        if(vt.fieldTypes[i]==vt.NUMERICAL)
          lastcolumn = i;
      for(int i=0;i<vt.colCount;i++)
        if(vt.fieldTypes[i]==vt.NUMERICAL)
          if(i!=lastcolumn)
            fw.write(vt.fieldNames[i]+"\t");
          else{
            fw.write(vt.fieldNames[i]);
          }
      fw.write("\n");
      fw.write("EWEIGHT\t\t\t");
      for(int i=0;i<vt.colCount;i++)
        if(vt.fieldTypes[i]==vt.NUMERICAL)
          fw.write("1\t");
      fw.write("\n");
      for(int i=0;i<vt.rowCount;i++){
        fw.write(vt.stringTable[i][k]+"\t"+vt.stringTable[i][k]+"\t1\t");
        if(conserveOrder)
          fw.write(""+(i+1)+"\t");
        for(int j=0;j<vt.colCount;j++)
          if(vt.fieldTypes[j]==vt.NUMERICAL)
            if(j!=lastcolumn)
              fw.write(vt.stringTable[i][j]+"\t");
            else
              fw.write(vt.stringTable[i][j]);
        fw.write("\n");
      }
      fw.close();
    }
    
    public static VDataTable fillMissingValues(VDataTable vt, int numberOfComponentsToUse) throws Exception{

        boolean num[] = new boolean[vt.colCount];
        int numfields = 0;
        for(int i=0;i<vt.colCount;i++){
          if(vt.fieldTypes[i]==vt.NUMERICAL) {numfields++; num[i]=true; }
          else num[i]=false;
        }

        VDataTable res = new VDataTable();
        res.copyHeader(vt);
        res.rowCount = vt.rowCount;
        res.stringTable = new String[res.rowCount][res.colCount];

        VDataSet ds = VSimpleProcedures.SimplyPreparedDataset(vt, -1); 
        
        PCAMethod pca = new PCAMethod();
        pca.setDataSet(ds);
        pca.calcBasis(numberOfComponentsToUse);
        VDataSet ds1 = pca.getProjectedDataset();
        
        for(int i=0;i<vt.rowCount;i++){
          float f[] = new float[numfields];
          float initialVector[] = ds.getVector(i);
          float projections[] = ds1.getVector(i);
          float projectedVector[] = pca.projectFromInToOut(projections);
          ds.unprocessFromSpace(projectedVector);
          
          int k=0;
          for(int j=0;j<vt.colCount;j++){
            if(num[j]){
            	// process the lines
            	if(Float.isNaN(initialVector[k]))
            		f[k] = projectedVector[k];
            	else
            		f[k] = Float.parseFloat(vt.stringTable[i][j]);
            	k++;
            }
          }
          k = 0;
          for(int j=0;j<vt.colCount;j++){
            if(vt.fieldTypes[j]==vt.NUMERICAL) { res.stringTable[i][j] = ""+f[k];  k++; }
            else res.stringTable[i][j] = vt.stringTable[i][j];
          }
        }
        return res;
      }

	public static void AnnotateWithSubsetOfGenes(VDataTable vt, String fileName, String fieldName){
		Vector<String> list = Utils.loadStringListFromFile(fileName);
		vt.addNewColumn(fieldName, "", "", vt.NUMERICAL, "0");
		vt.addNewColumn(fieldName+"_NAME", "", "", vt.STRING, "");		
		for(int i=0;i<vt.rowCount;i++){
			String gene = vt.stringTable[i][vt.fieldNumByName("GENE")];
			StringTokenizer st = new StringTokenizer(gene,";");
			Vector<String> gene_list = new Vector<String>();
			while(st.hasMoreTokens()) gene_list.add(st.nextToken());
			Vector<String> names = new Vector<String>(); 
			for(int j=0;j<list.size();j++){
				if(gene_list.contains(list.get(j))){
					names.add(list.get(j));
				}
			}
			String name = ""; for(int j=0;j<names.size();j++) name+=names.get(j)+";"; if(name.length()>0) name = name.substring(0, name.length()-1);			
			vt.stringTable[i][vt.fieldNumByName(fieldName)] = ""+names.size();
			vt.stringTable[i][vt.fieldNumByName(fieldName+"_NAME")] = name;					
		}
	}
	
	public static void AnnotateWithGMTfile(VDataTable vt, String fileName, String fieldName){
		Vector<GESignature> sets = GMTReader.readGMTDatabase(fileName);
		vt.addNewColumn(fieldName, "", "", vt.STRING, "");		
		vt.addNewColumn(fieldName+"_NUM", "", "", vt.NUMERICAL, "0");				
		for(int i=0;i<vt.rowCount;i++){
			String gene = vt.stringTable[i][vt.fieldNumByName("GENE")];
			Vector<String> names = new Vector<String>(); 
			for(int j=0;j<sets.size();j++){
				StringTokenizer st = new StringTokenizer(gene,";");
				while(st.hasMoreTokens()){
				String g = st.nextToken();
				if(sets.get(j).geneNames.contains(g)){
					if(!names.contains(g))
						names.add(sets.get(j).name);
				}
				}
			}
			String name = ""; for(int j=0;j<names.size();j++) name+=names.get(j)+";"; if(name.length()>0) name = name.substring(0, name.length()-1);
			vt.stringTable[i][vt.fieldNumByName(fieldName)] = name;
			vt.stringTable[i][vt.fieldNumByName(fieldName+"_NUM")] = ""+names.size();								
		}
	}
	
	public static void AnnotateWithMMetagene(VDataTable vt, String fileName, String fieldName){
		Vector<String> lines = Utils.loadStringListFromFile(fileName);
		Vector<String> gnames = new Vector<String>();
		Vector<Float> gvalues = new Vector<Float>();
		for(String s: lines){
			StringTokenizer st = new StringTokenizer(s,"\t");
			gnames.add(st.nextToken());
			gvalues.add(Float.parseFloat(st.nextToken()));
		}
		vt.addNewColumn(fieldName, "", "", vt.STRING, "");		
		vt.addNewColumn(fieldName+"_VALUE", "", "", vt.NUMERICAL, "0");
		for(int i=0;i<vt.rowCount;i++){
			String gene = vt.stringTable[i][vt.fieldNumByName("GENE")];
			Vector<String> names = new Vector<String>();
			Vector<Float> values = new Vector<Float>();
			for(int j=0;j<gnames.size();j++){
				StringTokenizer st = new StringTokenizer(gene,";");
				while(st.hasMoreTokens()){
				String g = st.nextToken();
				if(gnames.get(j).equals(g)){
					if(!names.contains(g)){
						names.add(gnames.get(j));
						values.add(gvalues.get(j));
					}
				}
				}
			}
			String name = ""; float value = 0f; 
			for(int j=0;j<names.size();j++){ 
				name+=names.get(j)+";";
				if(Math.abs(values.get(j))>value)
					value = values.get(j);
			}
			if(name.length()>0) name = name.substring(0, name.length()-1);
			vt.stringTable[i][vt.fieldNumByName(fieldName)] = name;
			vt.stringTable[i][vt.fieldNumByName(fieldName+"_VALUE")] = ""+value;								
		}		
	}
	
	public static void makeTableCorrelationGraph(String table1, String prefix1, String table2, String prefix2, String folder, float correlationThreshold) throws Exception{
		VDataTable vt1 = VDatReadWrite.LoadFromVDatFile(table1);
		vt1.makePrimaryHash(vt1.fieldNames[0]);
		VDataTable vt2 = VDatReadWrite.LoadFromVDatFile(table2);
		vt2.makePrimaryHash(vt2.fieldNames[0]);
		// make common list of objects
		HashSet<String> names_set = new HashSet<String>();
		Vector<String> names = new Vector<String>();
		for(int i=0;i<vt1.rowCount;i++){
			String name1 = vt1.stringTable[i][0];
			if(vt2.tableHashPrimary.get(name1)!=null){
				if(!names_set.contains(name1))
					names_set.add(name1);
			}
		}
		for(String s: names_set)
			names.add(s);
		Collections.sort(names);
		System.out.println(names.size()+" common objects are found.");
		
		FileWriter fw  = new FileWriter(folder+prefix1+"_"+prefix2+".txt");
		fw.write("FIELD1\tFIELD2\tCORRELATION\tCORRELATION_ABS\n");
		
		for(int i=0;i<vt1.colCount;i++)if(vt1.fieldTypes[i]==vt1.NUMERICAL){
			String fni = vt1.fieldNames[i];
			for(int j=0;j<vt2.colCount;j++)if(vt2.fieldTypes[j]==vt2.NUMERICAL){
				String fnj = vt2.fieldNames[j];
				float xi[] = new float[names.size()];
				float xj[] = new float[names.size()];
				for(int k=0;k<names.size();k++){
					xi[k] = Float.parseFloat(vt1.stringTable[vt1.tableHashPrimary.get(names.get(k)).get(0)][i]);
					xj[k] = Float.parseFloat(vt2.stringTable[vt2.tableHashPrimary.get(names.get(k)).get(0)][j]);
				}
				float corr = VSimpleFunctions.calcCorrelationCoeff(xi, xj);
				float abscorr = Math.abs(VSimpleFunctions.calcCorrelationCoeff(xi, xj));
				if(abscorr>=correlationThreshold){
					fw.write(fni+"_"+prefix1+"\t"+fnj+"_"+prefix2+"\t"+corr+"\t"+Math.abs(corr)+"\n");
				}
			}
		}
		
		fw.close();
	}
	
	public static void findAllNumericalColumns(VDataTable vt){
		for(int i=0;i<vt.colCount;i++){
			boolean isNumerical = true;
			for(int j=0;j<vt.rowCount;j++){
				String s = vt.stringTable[j][i].trim();
				if(s.equals("\"\"")||s.equals("NA")||s.equals("")||s.equals("_")){
					s = "@";
					vt.stringTable[j][i] = s;
				}else{
				try{
					Float f = Float.parseFloat(s);
				}catch(Exception e){
					//e.printStackTrace();
					isNumerical = false;
				}}
			}
			if(isNumerical)
				vt.fieldTypes[i] = vt.NUMERICAL;
		}
	}
	
	public static int[] countNumberOfDistinctValuesInColumns(VDataTable vt){
		int num[] = new int[vt.colCount];
		for(int i=0;i<vt.colCount;i++){
			HashSet<String> set = new HashSet<String>();
			for(int j=0;j<vt.rowCount;j++){
				String s = vt.stringTable[j][i].trim();
				if((!s.equals("\"\""))&&(!s.equals("NA"))&&(!s.equals(""))&&(!s.equals("@")&&(!s.equals("_")))){
					if(!set.contains(s))
						set.add(s);
				}
			}
			num[i] = set.size();
		}
		return num;
	}
	
	public static float[][] doubleCenterMatrix(float matrix[][]){
		float res[][] = matrix.clone();
		
		int n = matrix.length;
		int m = matrix[0].length;
		float sumPerRow[] = new float[n];
		int countInRow[] = new int[n];
		float sumPerColumn[] = new float[m];
		int countInColumn[] = new int[m];
		int countTotal = 0;
		float totalSum = 0f;
		for(int i=0;i<n;i++){
			for(int j=0;j<m;j++)if(!Float.isNaN(matrix[i][j])){
				totalSum+=matrix[i][j];
				sumPerRow[i]+=matrix[i][j];
				countInRow[i]++;
				countTotal++;
			}
		}
		for(int j=0;j<m;j++){
			for(int i=0;i<n;i++)if(!Float.isNaN(matrix[i][j])){
				sumPerColumn[j]+=matrix[i][j];
				countInColumn[j]++;
			}
		}
		for(int i=0;i<n;i++){
			for(int j=0;j<m;j++){
				res[i][j] = matrix[i][j]-1/(float)countInRow[i]*sumPerRow[i]-1/(float)countInColumn[j]*sumPerColumn[j]+1/countTotal*totalSum;
			}
		}
		
		return res;
	}



}