package compression_algorithms.sequitur;

import java.util.ArrayList;
import java.util.Vector;

public class rule{


  // Guard symbol to mark beginning
  // and end of rule.

  public guard theGuard;

  // Counter keeps track of how many
  // times the rule is used in the
  // grammar.

  public int count;

  // The total number of rules.

  public static int numRules = 0;

  // The rule's number.
  // Used for identification of
  // non-terminals.

  public int number;

  // Index used for printing.

  public int index;
  
  // Number of total characters.
  
  public int length;

  rule(){
    number = numRules;
    numRules++;
    theGuard = new guard(this);
    count = 0;
    index = 0;
    length = 0; 
  }

  public symbol first(){
    return theGuard.n;
  }
  
  public symbol last(){
    return theGuard.p;
  }

  public void getRules(ArrayList<String> s, boolean compute){
	Vector<rule> rules = new Vector<rule>(numRules);
    rule currentRule;
    rule referedTo;
    symbol sym;
    int index;
    int processedRules = 0;
    StringBuffer text = new StringBuffer();

    rules.addElement(this);
    currentRule = (rule)rules.elementAt(processedRules);
    text.append(' ');
    for (sym=currentRule.first();(!sym.isGuard());sym=sym.n){
    	if (sym.isNonTerminal()){
    		referedTo = ((nonTerminal)sym).r;
    		if ((rules.size() > referedTo.index) &&
    				((rule)rules.elementAt(referedTo.index) ==
    				referedTo)){
    			index = referedTo.index;
    		}else{
    			index = rules.size();
    			referedTo.index = index;
    			rules.addElement(referedTo);
    		}
    		text.append(index);
    	}else{
    		if (sym.value == ' '){
    			text.append('_');
    		}else{
    			if (sym.value == '\n'){
    				text.append("\\n");
    			}else if(sym.value == -1) {
    				text.deleteCharAt(0);
    				s.add(new String(text));
    				text.setLength(0);
    			} else {
    				text.append((char)sym.value);
    			}
    		}
    	}
    	text.append(' ');
    }
	text.setLength(0);
	if(!compute) return;
	Sequitur.firstRule.getRulesVerbose();
	Sequitur.similarity = new double[Sequitur.rules.size()][Sequitur.rules.size()];
	for(double[] d : Sequitur.similarity) {
		for(int i = 0; i < d.length; i++) {
			d[i] = -1;
		}
	}
	Sequitur.firstRule.calculateLength();
  }
  
  public String getRule0(){
		Vector<rule> rules = new Vector<rule>(numRules);
	    rule currentRule;
	    rule referedTo;
	    symbol sym;
	    int index;
	    int processedRules = 0;
	    StringBuffer text = new StringBuffer();

	    rules.addElement(this);
	    currentRule = (rule)rules.elementAt(processedRules);
	    text.append(' ');
	    for (sym=currentRule.first();(!sym.isGuard());sym=sym.n){
	    	if (sym.isNonTerminal()){
	    		referedTo = ((nonTerminal)sym).r;
	    		if ((rules.size() > referedTo.index) &&
	    				((rule)rules.elementAt(referedTo.index) ==
	    				referedTo)){
	    			index = referedTo.index;
	    		}else{
	    			index = rules.size();
	    			referedTo.index = index;
	    			rules.addElement(referedTo);
	    		}
	    		text.append(index);
	    	}else{
	    		if (sym.value == ' '){
	    			text.append('_');
	    		}else{
	    			if (sym.value == '\n'){
	    				text.append("\\n");
	    			}else if(sym.value == -1) {
	    				text.deleteCharAt(0);
	    				return text.toString();
	    			} else {
	    				text.append((char)sym.value);
	    			}
	    		}
	    	}
	    	text.append(' ');
	    }
	    return text.toString();
  }
  
  public String getRulesVerbose(){
	    Vector<rule> rules = new Vector<rule>(numRules);
	    rule currentRule;
	    rule referedTo;
	    symbol sym;
	    int index;
	    int processedRules = 0;
	    StringBuffer text = new StringBuffer();

	    rules.addElement(this);
	    Sequitur.rules.add(this);
	    while (processedRules < rules.size()){
	      currentRule = (rule)rules.elementAt(processedRules);
	      if(processedRules > 0) {
	    	  text.append("R");
		      text.append(processedRules + " " + currentRule.length);
		      text.append(" -> ");
	      }
	      for (sym=currentRule.first();(!sym.isGuard());sym=sym.n){
	        if (sym.isNonTerminal()){
	          referedTo = ((nonTerminal)sym).r;
	          if ((rules.size() > referedTo.index) &&
	              ((rule)rules.elementAt(referedTo.index) ==
	               referedTo)){
	            index = referedTo.index;
	          }else{
	            index = rules.size();
	            referedTo.index = index;
	            rules.addElement(referedTo);
	    	    Sequitur.rules.add(referedTo);
	          }
	          if(processedRules > 0) {
		          text.append('R');
		          text.append(index);
	          }
	        }else if(processedRules > 0) {
	          if (sym.value == ' '){
	            text.append('_');
	          }else{
	            if (sym.value == '\n'){
	              text.append("\\n");
	            }else
	              text.append((char)sym.value);
	          }
	        }
	        if(processedRules > 0) {
	        	text.append(' ');
	        }
	      }
	      if(processedRules > 0) {
	    	  text.append('\n');
	      }
	      processedRules++;
	    }
	    return new String(text);
	  }
  
  /**
   * 
   * @param rule1
   * @param rule2
   * @return Percentage of rule1 in rule2
   */
  public double similarity(String rule1, String rule2) {
	  int r1, r2;
	  rule currentRule;
	  symbol sym;
	  int sum = 0, counter = 0;

	  r1 = Integer.parseInt(rule1);
	  r2 = Integer.parseInt(rule2);

	  if(Sequitur.similarity[r1][r2] != -1) {
		  return Sequitur.similarity[r1][r2];
	  } else if(r1 == r2) {
		  Sequitur.similarity[r1][r2] = 1;
		  return 1;
	  } else {
		  currentRule = Sequitur.rules.get(r2);
		  for (sym=currentRule.first();(!sym.isGuard());sym=sym.n){
			  if(sym.isNonTerminal()) {
				  sum += similarity(rule1, ((nonTerminal)sym).r.index + "") * ((nonTerminal)sym).r.length;
				  counter += ((nonTerminal)sym).r.length;
			  } else {
				  counter++;
			  }
		  }
		  
		  Sequitur.similarity[r1][r2] = sum/(double)counter;
		  return sum/(double)counter;
	  }		  
  }
  
  public void calculateLength() {
	  symbol sym;
	  for(sym = first(); (!sym.isGuard());sym=sym.n) {
		  if(sym.isNonTerminal()) {
			  if(((nonTerminal)sym).r.length == 0) ((nonTerminal)sym).r.calculateLength();
			  length += ((nonTerminal)sym).r.length;
		  } else {
			  length++;
		  }
	  }
  }
  
  @Override
	public String toString() {
	  return "R" + index;
	}

}

