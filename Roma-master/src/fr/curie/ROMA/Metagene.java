package fr.curie.ROMA;

import java.util.*;
import java.io.*;

public class Metagene extends GESignature {

  public Vector weights = null;
  public Vector weightSpecified = null;  
  public Vector sampleNames = null;
  public Vector samplePattern = null;
  
  public float correctedSign = 1f;

  public Metagene() {
  }

  public Metagene(GESignature gs) {
    for(int i=0;i<gs.probeSets.size();i++){
      probeSets.add((String)gs.probeSets.elementAt(i));
    }
    for(int i=0;i<gs.geneNames.size();i++){
    	geneNames.add((String)gs.geneNames.elementAt(i));
    }
    name = gs.name;
  }

  public void initializeWeightsByOnes(){
    weights = new Vector();
    weightSpecified = new Vector();
    for(int i=0;i<probeSets.size();i++){
      Float f = new Float(1f);
      weights.add(f);
      weightSpecified.add(new Boolean(false));
    }
  }

  public void initializeWeightsByZeros(){
    weights = new Vector();
    weightSpecified = new Vector();    
    for(int i=0;i<probeSets.size();i++){
      Float f = new Float(0f);
      weights.add(f);
      weightSpecified.add(new Boolean(false));
    }
  }

  public void saveToFile(String fn){
    try{
      FileWriter fw = new FileWriter(fn);
      for(int i=0;i<probeSets.size();i++){
        fw.write(probeSets.elementAt(i)+"\t"+((Float)weights.elementAt(i)).floatValue()+"\r\n");
      }
      fw.close();
    }catch(Exception e){
      e.printStackTrace();
    }
  }


}