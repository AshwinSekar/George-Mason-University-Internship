package compression_algorithms.sequitur;

public class nonTerminal extends symbol implements Cloneable{

	  rule r;

	  nonTerminal(rule theRule){
	    r = theRule;
	    r.count++;
	    value = numTerminals+r.number;
	    p = null;
	    n = null;
	  }

	  /**
	   * Extra cloning method necessary so that
	   * count in the corresponding rule is
	   * increased.
	   */

	  protected Object clone(){

	    nonTerminal sym = new nonTerminal(r);

	    sym.p = p;
	    sym.n = n;
	    return sym;
	  }

	  public void cleanUp(){
	    join(p,n);
	    deleteDigram();
	    r.count--;
	  }

	  public boolean isNonTerminal(){
	    return true;
	  }

	  /**
	   * This symbol is the last reference to
	   * its rule. The contents of the rule
	   * are substituted in its place.
	   */

	  public void expand(){
	    join(p,r.first());
	    join(r.last(),n);

	    // Necessary so that garbage collector
	    // can delete rule and guard.

	    r.theGuard.r = null;
	    r.theGuard = null;
	  }
	}
