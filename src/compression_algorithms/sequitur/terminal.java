package compression_algorithms.sequitur;

public class terminal extends symbol implements Cloneable{

	  terminal(int theValue){
	    value = theValue;
	    p = null;
	    n = null;
	  }
	  
	  public void cleanUp(){
	    join(p,n);
	    deleteDigram();
	  }
	}