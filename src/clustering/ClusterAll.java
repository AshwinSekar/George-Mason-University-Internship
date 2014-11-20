package clustering;

import static compression_algorithms.sequitur.Sequitur.compress;
import static compression_algorithms.sequitur.Sequitur.init;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;

import net.jpountz.lz4.LZ4;

import compression_algorithms.sequitur.Sequitur;


public class ClusterAll {
	static BufferedReader f =  null;
	static PrintWriter o = null;

	static ArrayList<String> sequences;
	static ArrayList<String> origSequences;
	static ArrayList<LZ4> sequencesLz4;
	static ArrayList<LZ4> origSequencesLz4;
	static ArrayList<String> sequencesSequitur;
	static ArrayList<String> origSequencesSequitur;
	static int size;

	static final int MATCHFACTOR = 2;
	static final int MISMATCHFACTOR = -1;
	static final int GAPPENALTY = -1;
	static final int GAPEXTENSION = -1;
	static final int[][] DNASUBSITUTIONMATRIX = {{0,2,1,2},{2,0,2,1},{1,2,0,2},{2,1,2,0}};

	static String fileName;
	static HashMap<String, Double> LSA;
	
	static int max, maxI,maxJ, h;
	static int[][] score = new int[61][61];
	static int[][] pointers = new int[61][61];

	static StringBuilder a = new StringBuilder();
	static StringBuilder b = new StringBuilder();

	public static void main(String[] args) {
		LSA = new HashMap<String, Double>();
		String[] asdf = new String[2];
		asdf[1] = "LSHDIV_DataFiles";
		fileName = args[0];
		try {
			System.out.println("Loading Sequences");
			BufferedReader a = new BufferedReader(new FileReader("PreComputed/" + fileName));
			String s;
			double d;
			while((s = a.readLine()) != null) {
				d = Double.parseDouble(a.readLine());
				LSA.put(s, d);
			}
			a.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		File s = new File(fileName);
		if(!s.exists()) s.mkdir();
		asdf[0] = fileName;
		ClusterCDMLz4.main(asdf);
		ClusterCDMSequitur.main(asdf);
		ClusterLz4.main(asdf);

		try {
			System.out.println("Reprinting LSA");
			PrintWriter o = new PrintWriter(new File("PreComputed/" + fileName));
			for(String b : LSA.keySet()) {
				o.println(b);
				o.println(LSA.get(b));
			}
			o.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

/**
 * Populates the sequences ArrayLis	t with the first n sequences from the dataset
 * 
 * @param n - Number of sequences
 */
public static void populateSequences(int n, boolean wsim) {
	sequences = new ArrayList<String>();
	String sequence;
	String curLine = " ";

	try {
		f.readLine(); // Skips past first header
		if(n != -1) {
			for(int i = 0; i < n ; i++) {
				sequence = ""; 
				while((curLine = f.readLine()).charAt(0) != '>') {
					curLine = curLine.replaceAll("n", ""); // Filter out the n's in the DNA sequence
					sequence += curLine;
				}
				sequences.add(sequence);
			}
		} else {
			while(true) {
				sequence = ""; 
				curLine = f.readLine();
				if(curLine == null) break;
				while((curLine).charAt(0) != '>') {
					curLine = curLine.replaceAll("n", ""); // Filter out the n's in the DNA sequence
					sequence += curLine;
					curLine = f.readLine();
					if(curLine == null) break;
				}
				sequences.add(sequence);
			}
		}
	} catch (IOException e) {
		e.printStackTrace();
	}
	size = sequences.size();
	origSequences = new ArrayList<String>(sequences);

	init();
	for(int i = 0; i < sequences.size(); i++) {
		compress(sequences.get(i));
	}
	sequencesSequitur = new ArrayList<String>();
	Sequitur.firstRule.getRules(sequencesSequitur, false);

	origSequencesSequitur =  new ArrayList<String>(sequencesSequitur);

	sequencesLz4 = new ArrayList<LZ4>();
	for(int i = 0; i < sequences.size(); i++) {
		try {
			sequencesLz4.add(LZ4.compress(sequences.get(i)));
		} catch (UnsupportedEncodingException e) {
			e.printStackTrace();
		}

	}
	origSequencesLz4 =  new ArrayList<LZ4>(sequencesLz4);

}

/**
 * Generates the statistics from the clustering in clusters
 * @param clusters
 * @return [Singletons,Doubletons, W. Sim, Chao1, Shannon, Sace]
 */
public static double[] generateStatistics(ArrayList<ArrayList<String>> clusters, boolean wsim) {
	double[] stat = new double[6];
	int[] freq =  new int[size+1];
	double total = 0;
	for(ArrayList<String> c : clusters) {
		if(c.size() == 1) stat[0]++;
		if(c.size() == 2) stat[1]++;
		freq[c.size()]++;
		total += (c.size()*(c.size()-1))/2;
	}

	if(wsim) {
		stat[2] = 0;
		double sim,den;
		System.out.print("       ");
		for(ArrayList<String> c : clusters) {
			sim = 0;
			den = 0;
			for(int i = 0; i < c.size(); i++) {
				for(int j = i+1; j < c.size(); j++) {
					if(LSA.containsKey(c.get(i) + "_" + c.get(j)))
						sim += LSA.get(c.get(i) + "_" + c.get(j));
					else if(LSA.containsKey(c.get(j) + "_" + c.get(i)))
						sim += LSA.get(c.get(j) + "_" + c.get(i));
					else {
						LSA.put(c.get(i) + "_" + c.get(j), localSeqAlignmentSimilarity(c.get(i),c.get(j)));
						sim += LSA.get(c.get(i) + "_" + c.get(j));
					}	
					den++;
					System.out.printf("\r %.2f", den/total);
				}
			}
			if(c.size() != 1)
				sim = sim/den;
			else
				sim = 1.0;
			stat[2] += (c.size()/((double) size)) * sim;
		}

		stat[2] = stat[2]*100;
		System.out.println();

	}

	int obs = clusters.size();
	stat[3] = obs + ((stat[0]*(stat[0]-1.0))/(2.0*(stat[1]+1.0)));


	for(int i = 0 ; i < clusters.size(); i++) {
		stat[4] += ( clusters.get(i).size()/ (double) size) * Math.log(clusters.get(i).size()/ (double) size);
	}
	stat[4] = -stat[4];
	int threshold = 10;
	int sRare = 0, sAbund = 0;

	for(ArrayList<String> s : clusters) {
		if(s.size() > threshold) sAbund++;
		else sRare++;
		freq[s.size()]++;
	}

	double nRare = 0;
	double a = 0;
	for(int i = 1 ; i <= threshold; i++) {
		nRare += i * freq[i];
		a += i*(i-1)*freq[i];
	}

	double cACE = 1 - freq[1]/nRare;

	double lamdaACE2 = Math.max((sRare/cACE) * (a/(nRare*(nRare -1)))-1, 0);
	double lamdaACE = Math.max(lamdaACE2*(1 + (nRare*(1 - cACE)*a)/(nRare*(nRare-cACE))),0);

	if(Math.sqrt(lamdaACE2) < .8) {
		stat[5] = sAbund + (sRare/cACE) + (freq[1]/cACE)*lamdaACE2;
	} else {
		stat[5] = sAbund + (sRare/cACE) + (freq[1]/cACE)*lamdaACE;
	}
	sequences = new ArrayList<String>(origSequences);
	sequencesLz4 = new ArrayList<LZ4>(origSequencesLz4);
	sequencesSequitur = new ArrayList<String>(origSequencesSequitur);
	System.out.println();
	System.out.println();
	return stat;
}

public static double localSeqAlignmentSimilarity(String s1, String s2) {
	//System.out.println(s1);
	//System.out.println(s2);
	//System.out.println();
	s1 = " " + s1;
	s2 = " " + s2;
	max = 0;
	h = 0;
	maxI = 0;
	maxJ = 0;
	//int[][] score = new int[s1.length()][s2.length()];
	//int[][] pointers = new int[s1.length()][s2.length()];
	//  Pointer convention. 1 = left, 2 = up, 3 = diagonal left up 
	for(int i = 1; i < s1.length(); i++) {
		pointers[i][0] = 2;
	}
	for(int i = 1; i < s2.length(); i++) {
		pointers[0][i] = 1;
	}
	boolean inGap = false;
	for(int i = 1; i < s1.length(); i++) {
		for(int j = 1; j < s2.length();  j++) {
			h = -99;
			if(score[i-1][j-1] + match(s1.charAt(i),s2.charAt(j)) > h) {
				h = score[i-1][j-1] + match(s1.charAt(i),s2.charAt(j));
				pointers[i][j] = 3;
				inGap = false;
			} 
			if(!inGap) {
				if(score[i-1][j] + GAPPENALTY > h) {
					h = score[i-1][j] + GAPPENALTY;
					pointers[i][j] = 2;
					inGap = true;
				} 
				if(score[i][j-1] + GAPPENALTY > h) {
					h = score[i][j-1] + GAPPENALTY;
					pointers[i][j] = 1;
					inGap = true;
				}
			} else {
				if(score[i-1][j] + GAPEXTENSION > h) {
					h = score[i-1][j] + GAPEXTENSION;
					pointers[i][j] = 2;
					inGap = true;
				} 
				if(score[i][j-1] + GAPEXTENSION > h) {
					h = score[i][j-1] + GAPEXTENSION;
					pointers[i][j] = 1;
					inGap = true;
				}
			}
			if(0 > h)
				h = 0;

			score[i][j] = h;
			if(h >= max) {
				max = h;
				maxI = i;
				maxJ = j;
			}
		}
	}

	/*
		for(int i = 0 ; i < s1.length(); i++) {
			for(int j = 0 ; j < s2.length(); j++) {
				System.out.print(score[i][j] + " ");
			}
			System.out.print("     ");
			for(int j = 0 ; j < s2.length(); j++) {
				System.out.print(pointers[i][j] + " ");
			}
			System.out.println();
		}
	 */


	double matches = 0;
	String o1 = "",  o2 = "";
	while(!(maxI == 0 && maxJ == 0)) {
		if(pointers[maxI][maxJ] == 3) {
			o1 += s1.charAt(maxI);
			o2 += s2.charAt(maxJ);
			maxI--;
			maxJ--;
		} else if(pointers[maxI][maxJ] == 2) {
			o1 += s1.charAt(maxI);
			o2 += "_";
			maxI--;
		} else if(pointers[maxI][maxJ] == 1) {
			o1 += "_";
			o2 += s2.charAt(maxJ);
			maxJ--;
		}
	}
	a.append(o1);
	b.append(o2);
	o1 = a.reverse().toString();
	o2 = b.reverse().toString();
	a.setLength(0);
	b.setLength(0);
	//System.out.println(o1);
	//System.out.println(o2);
	for(int i = 0; i < Math.min(o1.length(), o2.length()); i++) {
		if(o1.charAt(i) == o2.charAt(i)) matches++;
	}
	return matches/Math.min(o1.length(), o2.length());
}

public static int match(char c1, char c2) {
	if(c1 == c2) return MATCHFACTOR;
	else return MISMATCHFACTOR;
	//return 2 - DNASUBSITUTIONMATRIX["ACGT".indexOf(c1)]["ACGT".indexOf(c1)];
}
}

