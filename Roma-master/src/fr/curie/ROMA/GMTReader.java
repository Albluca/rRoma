package fr.curie.ROMA;

import java.io.*;
import java.util.*;
import vdaoengine.data.*;
import vdaoengine.utils.*;
import vdaoengine.analysis.*;

public class GMTReader {
  public static void main(String[] args) {

	MergeUPDNSignatures("C:/Datas/ROMA/data/mosaic/oncosig.gmt",true); System.exit(0);
	  
	
    //Vector db = readGMTDatabase("c:/datas/msigdb/c2.symbols.gmt");
    Vector db = readGMTDatabase("c:/datas/msigdb/c2.symbols.gmt");
    //writeNameSpaceGMT(db,"c:/datas/msigdb/c2.symbols_u133ab.gmt","c:/datas/msigdb/annots/HG_U133AB_annot.txt");
    //VDataTable vt = vdaoengine.data.io.VDatReadWrite.LoadFromVDatFile("c:/datas/breastcancer/wang/wangn_a.dat");
    //writeNameSpaceGMTVariationFiltered(db,"c:/datas/msigdb/c2.symbols_u133ab_vf_wang.gmt","c:/datas/msigdb/annots/HG_U133AB_annot.txt",vt,"CHIP");

    //VDataTable vt = vdaoengine.data.io.VDatReadWrite.LoadFromVDatFile("c:/datas/breastcancer/onlya/an_a.dat");
    //writeNameSpaceGMTVariationFiltered(db,"c:/datas/msigdb/c2.symbols_u133ab_vf_curie.gmt","c:/datas/msigdb/annots/HG_U133AB_annot.txt",vt,"CHIP");

    //VDataTable vt = vdaoengine.data.io.VDatReadWrite.LoadFromVDatFile("c:/datas/breastcancer/160606/wang_basal/wang_basal_na.dat");
    //writeNameSpaceGMTVariationFiltered(db,"c:/datas/msigdb/c2.symbols_u133ab_vf_wangbasal.gmt","c:/datas/msigdb/annots/HG_U133AB_annot.txt",vt,"CHIP");
    //GMTReader.writeGCTandCLS(vt,"c:/datas/breastcancer/160606/wang_basal/wang_basal_na",0,"CHIP","GeneSymbol");

    //VDataTable vt = vdaoengine.data.io.VDatReadWrite.LoadFromVDatFile("c:/datas/breastcancer/160606/curie_basal/curie_basal_na.dat");
    //writeNameSpaceGMTVariationFiltered(db,"c:/datas/msigdb/c2.symbols_u133ab_vf_curiebasal.gmt","c:/datas/msigdb/annots/HG_U133AB_annot.txt",vt,"CHIP");
    //GMTReader.writeGCTandCLS(vt,"c:/datas/breastcancer/160606/curie_basal/curie_basal_na",0,"CHIP","GeneSymbol");

    String prefix = "c:/datas/breastcancer/160606/wang_basal/wang_basal_na";
    VDataTable vt = vdaoengine.data.io.VDatReadWrite.LoadFromVDatFile(prefix+".dat");
    VDataTable vtt = vt.transposeTable("CHIP");
    vdaoengine.data.io.VDatReadWrite.saveToVDatFile(vtt,prefix+"_transposed.dat");
    VDataTable vtav = TableUtils.substituteRowsByStatistics(vtt,"D1",0);
    VDataTable vtavt = vtav.transposeTable("D1");
    vdaoengine.data.io.VDatReadWrite.saveToVDatFile(vtavt,prefix+"_averaged.dat");

    VDataTable vtsv = TableUtils.substituteRowsByStatistics(vtt,"D1",1);
    VDataTable vtsvt = vtsv.transposeTable("D1");
    vdaoengine.data.io.VDatReadWrite.saveToVDatFile(vtsvt,prefix+"_stdv.dat");



    //MetageneLoader ml = new MetageneLoader();
    //ml.loadSignatures("c:/datas/breastcancer/signatures");
    //saveSignaturesTOGMTFormat(ml.listOfGESignatures,"c:/datas/msigdb/signatures.gmt");

  }
  
  public static Vector<GESignature> readGMTDatabase(String fin){
	  return readGMTDatabase(fin, -1);
  }
  
  public static Vector<GESignature> readGMTDatabase(String fin, int minimalNumberOfGenesInModule){
	  return readGMTDatabase(fin, minimalNumberOfGenesInModule, Integer.MAX_VALUE);
  }

  public static Vector<GESignature> readGMTDatabase(String fin, int minimalNumberOfGenesInModule, int maximalNumberOfGenesInModule){
    Vector res = new Vector();
    try{
      LineNumberReader lr = new LineNumberReader(new FileReader(fin));
      String s = null;
      while((s=lr.readLine())!=null)if(!s.trim().equals("")){
        StringTokenizer st = new StringTokenizer(s,"\t");
        GESignature gs = new GESignature();
        gs.name = st.nextToken();
        gs.description = st.nextToken();
        String sname = null;
        while(st.hasMoreTokens()){
          String name = st.nextToken();
          sname = name;
          float weight = 1f;
          if(name.indexOf("[")>=0){
          	sname = name.substring(0,name.indexOf("["));
          	String sweight = name.substring(name.indexOf("[")+1,name.length()-1);
          	weight = Float.parseFloat(sweight);
          }
          if(!gs.geneNames.contains(sname))
        	  gs.geneNames.add(sname);
        }
        if((minimalNumberOfGenesInModule<0)||((gs.geneNames.size()>=minimalNumberOfGenesInModule)&&(gs.geneNames.size()<=maximalNumberOfGenesInModule)))
        	res.add(gs);
      }
    }catch(Exception e){
      e.printStackTrace();
    }
    return res;
  }
  
  public static Vector readGMTDatabaseWithWeights(String fin){
	  return readGMTDatabaseWithWeights(fin,-1);
  }
  
  public static Vector readGMTDatabaseWithWeights(String fin, int minimalNumberOfGenesInModule){
	  return readGMTDatabaseWithWeights(fin,minimalNumberOfGenesInModule, Integer.MAX_VALUE);
  }
  
  public static Vector readGMTDatabaseWithWeights(String fin, int minimalNumberOfGenesInModule, int maximalNumberOfGenesInModule){
	    Vector res = new Vector();
	    try{
	      LineNumberReader lr = new LineNumberReader(new FileReader(fin));
	      String s = null;
	      while((s=lr.readLine())!=null){
	        StringTokenizer st = new StringTokenizer(s,"\t");
	        Metagene gs = new Metagene();
	        gs.weights = new Vector();
	        gs.weightSpecified = new Vector(); 
	        gs.name = st.nextToken();
	        gs.description = st.nextToken();
	        while(st.hasMoreTokens()){
	          String name = st.nextToken();
	          String sname = name;
	          float weight = 1f;
	          if(name.indexOf("[")>=0){
	          	sname = name.substring(0,name.indexOf("["));
	          	String sweight = name.substring(name.indexOf("[")+1,name.length()-1);
	          	weight = Float.parseFloat(sweight);
	          	gs.weightSpecified.add(new Boolean(true));
	          }else
	        	gs.weightSpecified.add(new Boolean(false));
	          if(!gs.geneNames.contains(sname)){
	        	  gs.geneNames.add(sname);
	        	  gs.weights.add(new Float(weight));
	          }
	        }
	        if((minimalNumberOfGenesInModule<0)||((gs.geneNames.size()>=minimalNumberOfGenesInModule))&&((gs.geneNames.size()<=maximalNumberOfGenesInModule)))
	        	res.add(gs);
	      }
	    }catch(Exception e){
	      e.printStackTrace();
	    }
	    return res;
	  }
  
  

  /*public static void writeNameSpaceGMT(Vector db, String fout, String annotf){
     HashMap ann = MetageneLoader.hashAnnotation(annotf,false);
     Vector notFound = new Vector();
     try{
       FileWriter fw = new FileWriter(fout);
       for(int i=0;i<db.size();i++){
         GESignature gs = (GESignature)db.elementAt(i);
         convertToNameSpace(gs,ann,notFound);
         fw.write(gs.name+"\t"+gs.description);
         for(int j=0;j<gs.probeSets.size();j++)
           fw.write("\t"+(String)gs.probeSets.elementAt(j));
         fw.write("\r\n");
       }
       fw.close();
     }catch(Exception e){
       e.printStackTrace();
     }
  }

  public static void writeNameSpaceGMTVariationFiltered(Vector db, String fout, String annotf, VDataTable vt, String chipID){
     HashMap ann = MetageneLoader.hashAnnotation(annotf,false);
     VTableQuery vtq = new VTableQuery();
     vtq.table = vt;
     vtq.hashTable(chipID);
     VDataSet vd = vdaoengine.utils.VSimpleProcedures.SimplyPreparedDataset(vt,-1);
     VStatistics vst = new VStatistics(vd.pointCount);
     for(int i=0;i<vd.coordCount;i++){
       float x[] = new float[vd.pointCount];
       for(int j=0;j<vd.pointCount;j++)
         x[j] = vd.massif[j][i];
       vst.addNewPoint(x);
     }
     vst.calculate();


     Vector notFound = new Vector();
     try{
       FileWriter fw = new FileWriter(fout);
       for(int i=0;i<db.size();i++){
         GESignature gs = (GESignature)db.elementAt(i);
         convertToNameSpaceVariationFiltered(gs,ann,notFound,vtq,vst,chipID);
         fw.write(gs.name+"\t"+gs.description);
         for(int j=0;j<gs.probeSets.size();j++)
           fw.write("\t"+(String)gs.probeSets.elementAt(j));
         fw.write("\r\n");
       }
       fw.close();
     }catch(Exception e){
       e.printStackTrace();
     }
  }
	*/

  public static void convertToNameSpace(GESignature gs, HashMap ann, Vector notFound){
    Vector gnames = new Vector();
    Vector gps = new Vector();
    for(int i=0;i<gs.geneNames.size();i++){
      String gn = (String)gs.geneNames.elementAt(i);
      GESignature gsig = (GESignature)ann.get(gn);
      if(gsig!=null)
      for(int j=0;j<gsig.probeSets.size();j++){
        String ps = (String)gsig.probeSets.elementAt(j);
        if(gps.indexOf(ps)<0){
          gnames.add(gn);
          gps.add(ps);
        }
      }
      if(gsig==null){
        if(notFound.indexOf(gn)<0){
          notFound.add(gn);
          System.out.println(""+(notFound.size())+"\t"+gn+"\t was not found");
        }
      }
    }
    gs.geneNames = gnames;
    gs.probeSets = gps;
  }

  public static void convertToNameSpaceVariationFiltered(GESignature gs, HashMap ann, Vector notFound, VTableQuery vtq, VStatistics vst, String chipID){
    Vector gnames = new Vector();
    Vector gps = new Vector();
    for(int i=0;i<gs.geneNames.size();i++){
      String gn = (String)gs.geneNames.elementAt(i);
      GESignature gsig = (GESignature)ann.get(gn);
      if(gsig!=null){
      selectProbeSetWithHighestVariation(gsig,vtq,vst,chipID);
      for(int j=0;j<gsig.probeSets.size();j++){
        String ps = (String)gsig.probeSets.elementAt(j);
        if(gps.indexOf(ps)<0){
          gnames.add(gn);
          gps.add(ps);
        }
      }
      }
      if(gsig==null){
        if(notFound.indexOf(gn)<0){
          notFound.add(gn);
          System.out.println(""+(notFound.size())+"\t"+gn+"\t was not found");
        }
      }
    }
    gs.geneNames = gnames;
    gs.probeSets = gps;
  }


  public static void selectProbeSetWithHighestVariation(GESignature gs, VTableQuery vtq, VStatistics vst, String chipID){
    float maxx = -1; int maxi = -1;
    for(int i=0;i<gs.probeSets.size();i++){
      String ps = (String)gs.probeSets.elementAt(i);
      int k = vtq.queryHash(ps);
      if(ps.startsWith("AFFX"))
        k = -1;
      if(k!=-1){
        float x = vst.getStdDev(k);
        if(x>maxx){ maxx = x; maxi = i;}
      }else{
        System.out.println(ps+"\tnot found");
      }
    }
    if(maxi!=-1){
      Vector newList = new Vector();
      newList.add(gs.probeSets.elementAt(maxi));
      gs.probeSets = newList;
    }
  }

  public static void saveSignaturesTOGMTFormat(Vector sig, String fout){
    try{
      FileWriter fw = new FileWriter(fout);
      for(int i=0;i<sig.size();i++){
        GESignature gs = (GESignature)sig.elementAt(i);
        fw.write(gs.name+"\t"+"signature");
        for(int j=0;j<gs.probeSets.size();j++)
            fw.write("\t"+(String)gs.probeSets.elementAt(j));
        fw.write("\r\n");
      }
      fw.close();
    }catch(Exception e){
      e.printStackTrace();
    }
  }

  public static void writeGCTandCLS(VDataTable vt, String fn, int descrNum, String nameField, String descriptionField){
    try{
      Vector v = new Vector();
      for(int i=0;i<vt.colCount;i++){
        if(vt.fieldTypes[i]==vt.NUMERICAL)
          v.add(new Integer(i));
      }

      FileWriter fwt = new FileWriter(fn+".gct");
      fwt.write("#1.2\r\n");
      fwt.write(vt.rowCount+"\t"+v.size()+"\r\n");
      fwt.write("Name\tDescription");
      for(int i=0;i<v.size();i++)
        fwt.write("\t"+vt.fieldNames[((Integer)v.elementAt(i)).intValue()]);
      fwt.write("\r\n");
      for(int i=0;i<vt.rowCount;i++){
        fwt.write(vt.stringTable[i][vt.fieldNumByName(nameField)]+"\t"+vt.stringTable[i][vt.fieldNumByName(descriptionField)]);
        for(int j=0;j<v.size();j++){
          fwt.write("\t"+vt.stringTable[i][((Integer)v.elementAt(j)).intValue()]);
        }
      fwt.write("\r\n");
      }
      fwt.close();

      FileWriter fwc = new FileWriter(fn+".cls");
      Vector vc = new Vector();
      Vector vcn = new Vector();
      for(int i=0;i<vt.colCount;i++){
        String s = vt.fieldInfo[i][descrNum];
        if((s!=null)&&(!s.equals("")))
          if(vc.indexOf(s)<0){
            vc.add(s);
            vcn.add(new Integer(1));
          }else{
            vcn.setElementAt(new Integer(((Integer)vcn.elementAt(vc.indexOf(s))).intValue()),vc.indexOf(s));
          }
      }
      fwc.write(""+v.size()+" "+vc.size()+" 1\r\n");
      fwc.write("#");
      for(int i=0;i<vc.size();i++)
        fwc.write(" "+(String)vc.elementAt(i));
      fwc.write("\r\n");
      for(int i=0;i<v.size();i++){
        String s = vt.fieldInfo[((Integer)v.elementAt(i)).intValue()][descrNum];
        fwc.write(""+vc.indexOf(s)+" ");
      }
      fwc.write("\r\n");
      fwc.close();

    }catch(Exception e){
      e.printStackTrace();
    }
  }
  
	  public static void MergeUPDNSignatures(String fn, boolean putSigns) {

	      GMTReader gmt = new GMTReader();
	      Vector sigs = gmt.readGMTDatabase(fn);

	      Vector newSig = new Vector();
	      String suffixdn = "";
	      String suffixup = "";
	      if(putSigns) suffixdn+="[-1]";
	      if(putSigns) suffixup+="[1]";


	      for(int i=0;i<sigs.size();i++){
	        String fndn = ((GESignature)sigs.get(i)).name;
	        if(fndn.endsWith("_DN")){
	          fndn = fndn.substring(0,fndn.length()-3);
	          for(int j=0;j<sigs.size();j++){
	            String fnup = ((GESignature)sigs.get(j)).name;
	            if(fnup.equals(fndn+"_UP")){
	              System.out.println("Found "+fndn);
	              GESignature gdn = (GESignature)sigs.get(i);
	              GESignature gup = (GESignature)sigs.get(j);
	              GESignature nsig = new GESignature();
	              for(int k=0;k<gdn.geneNames.size();k++)
	                nsig.geneNames.add(gdn.geneNames.get(k)+suffixdn);
	              for(int k=0;k<gup.geneNames.size();k++)
	                nsig.geneNames.add(gup.geneNames.get(k)+suffixup);
	              nsig.name = fndn;
	              nsig.probeSets = nsig.geneNames;
	              newSig.add(nsig);
	            }
	          }
	        }
	      }

	      String fn1 = fn.substring(0, fn.length()-4)+"_mergedDNUP.gmt";
	      GMTReader.saveSignaturesTOGMTFormat(newSig,fn1);
	  }

  

}