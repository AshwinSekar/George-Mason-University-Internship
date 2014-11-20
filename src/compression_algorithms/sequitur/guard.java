package compression_algorithms.sequitur;

public class guard extends symbol{

	  rule r;
	  boolean isMiddle;

	  guard(rule theRule){	
	    r = theRule;
	    value = -1;
	    p = this;
	    n = this;
	    isMiddle = false;
	  }
	  
	  public guard(rule theRule, boolean t){
		    r = theRule;
		    value = -1;
		    p = this;
		    n = this;
		    isMiddle = t;
		  }


	  public void cleanUp(){
	    join(p,n);
	  }

	  public boolean isGuard(){
	    return !isMiddle;
	  }

	  public void deleteDigram(){
	    
	    // Do nothing
	  }
	  
	  public boolean check(){
	    return false;
	  }
	}