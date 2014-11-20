package compression_algorithms.sequitur;

import java.util.ArrayList;
import java.util.Hashtable;

public class Sequitur {
	public static rule firstRule;
	public static double[][] similarity;
	public static ArrayList<rule> rules;
	
	public static void compress(String s) {
		for(int i = 0; i < s.length(); i++){
			firstRule.last().insertAfter(new terminal(s.charAt(i)));
			firstRule.last().p.check();
		}

		firstRule.last().insertAfter(new guard(firstRule, true));
	};
	
	public static void init() {
		rules = new ArrayList<rule>();
		rule.numRules = 0;
		firstRule = new rule();
		symbol.theDigrams = new Hashtable<symbol, symbol>();
	}
	
	public static void main(String[] args) {
		firstRule = new rule();
		rules = new ArrayList<rule>();
		compress("he pushed her from here to there");
		//compress("The bank manager greeted the banker");
		System.out.println(firstRule.getRule0());
		System.out.println(firstRule.getRulesVerbose());
	}
}
