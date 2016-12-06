package fr.curie.ROMA;

import java.util.*;

public class GESignature {

  public String name = "unknown";
  public String description = "";
  public Vector<String> probeSets = new Vector();
  public Vector<String> geneNames = new Vector();
  public Vector locusPosition = new Vector();
  public Vector annotation = new Vector();

  public static int findSignatureByName(String name, Vector sigs){
	  int r = -1;
	  for(int i=0;i<sigs.size();i++){
		  if(((GESignature)sigs.get(i)).name.equals(name))
			  r = i;
	  }
	  return r;
  }
  
  public static Vector<GESignature> getAllSignaturesContainingGene(String id, Vector<GESignature> sigs){
	  Vector<GESignature> res = new Vector<GESignature>();
	  for(int i=0;i<sigs.size();i++)
		  if(sigs.get(i).geneNames.contains(id))
			  res.add(sigs.get(i));
	  return res;
  }
}