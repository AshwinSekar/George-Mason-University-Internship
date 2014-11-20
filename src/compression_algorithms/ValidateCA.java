package compression_algorithms;

import static compression_algorithms.RLE.RunLengthEncoding.encode;
import static compression_algorithms.sequitur.Sequitur.compress;
import static compression_algorithms.sequitur.Sequitur.init;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import net.jpountz.lz4.LZ4;

import compression_algorithms.sequitur.Sequitur;

public class ValidateCA {
	static BufferedReader f;
	static PrintWriter o;
	
	static ArrayList<String> sequences;
	static ArrayList<String> sequencesSequitur;
	static ArrayList<LZ4> sequenceslz4;
	static ArrayList<String> sequencesRLE;
	
	static ArrayList<ArrayList<LocalSequenceAlignment>> rankings;
	
	static ArrayList<ArrayList<JaccardCoefficient>> rankingsSequitur;
	static int[] matchesSequitur;
	
	static ArrayList<ArrayList<EditDistance>> rankingslz4;
	static int[] matcheslz4;
	
	static ArrayList<ArrayList<EditDistance>> rankingsRLE;
	static int[] matchesRLE;
	
	static ArrayList<ArrayList<CDM>> rankingsCDMSequitur;
	static int[] matchesCDMSequitur;
	
	static ArrayList<ArrayList<CDM>> rankingsCDMLz4;
	static int[] matchesCDMLz4;
	
	static ArrayList<ArrayList<CDM>> rankingsCDMRLE;
	static int[] matchesCDMRLE;
	
	static final int MATCHFACTOR = 2;
	static final int MISMATCHFACTOR = -1;
	static final int GAPPENALTY = -1;
	static final int GAPEXTENSION = -1;
	static final int[][] DNASUBSITUTIONMATRIX = {{0,2,1,2},{2,0,2,1},{1,2,0,2},{2,1,2,0}};

	/**
	 * @param args	
	 */
	public static void main(String[] args) {
		f = null;
		o = null;
		try {
			f = new BufferedReader(new FileReader("LSHDIV_DataFiles/60_Equal_53R.fasta"));
			o = new PrintWriter(new BufferedWriter(new FileWriter("Compression_Results.csv")));
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		System.out.println("");
		populateSequences(-1);
		System.out.println("Loaded sequences");
		compressSequitur();
		System.out.println("Sequitur finished");
		compresslz4();
		System.out.println("Lz4 finished");
		compressRLE();
		System.out.println("RLE finished");
		

		
		
		outputResults();
				
		try {
			f.close();
			o.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}	
	
	/**
	 * Populates the sequences ArrayList with the first n sequences from the dataset
	 * 
	 * @param n - Number of sequences
	 */
	public static void populateSequences(int n) {
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
	}
	
	/**
	 *  Runs all the sequences through the Sequitur algorithm and populates
	 */
	public static void compressSequitur() {
		init();
		for(int i = 0; i < sequences.size(); i++) {
			compress(sequences.get(i));
		}
		sequencesSequitur = new ArrayList<String>();
		Sequitur.firstRule.getRules(sequencesSequitur, true);
	}
	
	/**
	 *  Runs all the sequences through the lz4 algorithm and populates
	 */
	public static void compresslz4() {
		sequenceslz4 = new ArrayList<LZ4>();
		
		for(int i = 0; i < sequences.size(); i++) {
			try {
				sequenceslz4.add(LZ4.compress(sequences.get(i)));
			} catch (UnsupportedEncodingException e) {
				e.printStackTrace();
			}

		}

	}

	/**
	 *  Runs all the sequences through the RLE algorithm and populates
	 */
	public static void compressRLE() {
		sequencesRLE = new ArrayList<String>();
		for(int i = 0; i < sequences.size(); i++) {
			sequencesRLE.add(encode(sequences.get(i)));
		}

	}

	/**
	 * Computes the edit distances between all the sequences in S
	 */
	public static void computeDistances() {
		HashMap<String,Double> LSA = new HashMap<String, Double>();
		try {
			System.out.println("Loading Precomputed");
			BufferedReader a = new BufferedReader(new FileReader("PreComputed/" + "60_Equal_53R.fasta"));
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
		
		System.out.println("Computing distances and rankings");
		System.out.println();
		System.out.println("LSA");
		rankings = new ArrayList<ArrayList<LocalSequenceAlignment>>();
		LocalSequenceAlignment temp;
		for(int i = 0; i < sequences.size(); i++) {
			System.out.print("\r" + (i+1) + "/" + sequences.size()  + "         ");
			rankings.add(new ArrayList<LocalSequenceAlignment>());
			for(int z = 0; z < 5; z++) rankings.get(i).add(new LocalSequenceAlignment((double) 0,0));
			for(int j = 0; j < i; j++) {
				if(LSA.containsKey(sequences.get(i) + "_" + sequences.get(j)))
					temp = new LocalSequenceAlignment(LSA.get(sequences.get(i) + "_" + sequences.get(j)),j);
				else if(LSA.containsKey(sequences.get(j) + "_" + sequences.get(i)))
					temp = new LocalSequenceAlignment(LSA.get(sequences.get(j) + "_" + sequences.get(i)),j);
				else {
					temp = new LocalSequenceAlignment((double) localSeqAlignment(sequences.get(i),sequences.get(j)),j);
				}
				temp = new LocalSequenceAlignment((double) localSeqAlignment(sequences.get(i), sequences.get(j)),j);
				if(temp.compareTo(rankings.get(i).get(4)) < 0) {
					rankings.get(i).set(4, new LocalSequenceAlignment(temp));
					Collections.sort(rankings.get(i));
				}
			}
			for(int j = i+1; j < sequences.size(); j++) {
				if(LSA.containsKey(sequences.get(i) + "_" + sequences.get(j)))
					temp = new LocalSequenceAlignment(LSA.get(sequences.get(i) + "_" + sequences.get(j)),j);
				else if(LSA.containsKey(sequences.get(j) + "_" + sequences.get(i)))
					temp = new LocalSequenceAlignment(LSA.get(sequences.get(j) + "_" + sequences.get(i)),j);
				else {
					temp = new LocalSequenceAlignment((double) localSeqAlignment(sequences.get(i),sequences.get(j)),j);
				}
				temp = new LocalSequenceAlignment((double) localSeqAlignment(sequences.get(i), sequences.get(j)),j);
				if(temp.compareTo(rankings.get(i).get(4)) < 0) {
					rankings.get(i).set(4, new LocalSequenceAlignment(temp));
					Collections.sort(rankings.get(i));
				}
			}
		}
		System.out.println("\nFinished LSA for S");
		
		System.out.println();
		System.out.println("Sequitur");
		JaccardCoefficient tempJ;
		rankingsSequitur = new ArrayList<ArrayList<JaccardCoefficient>>();
		for(int i = 0; i < sequencesSequitur.size(); i++) {
			System.out.print("\r" + (i+1) + "/" + sequences.size()  + "         ");
			rankingsSequitur.add(new ArrayList<JaccardCoefficient>());
			for(int z = 0; z < 5; z++) rankingsSequitur.get(i).add(new JaccardCoefficient(0,0));

			for(int j = 0; j < i; j++) {
				tempJ = new JaccardCoefficient(jaccard(sequencesSequitur.get(i), sequencesSequitur.get(j)),j);
				if(tempJ.compareTo(rankingsSequitur.get(i).get(4)) < 0) {
					rankingsSequitur.get(i).set(4, new JaccardCoefficient(tempJ));
					Collections.sort(rankingsSequitur.get(i));
				}
			}
			for(int j = i+1; j < sequencesSequitur.size(); j++) {
				tempJ = new JaccardCoefficient(jaccard(sequencesSequitur.get(i), sequencesSequitur.get(j)),j);
				if(tempJ.compareTo(rankingsSequitur.get(i).get(4)) < 0) {
					rankingsSequitur.get(i).set(4, new JaccardCoefficient(tempJ));
					Collections.sort(rankingsSequitur.get(i));
				}			}
		}
		System.out.println("\nFinished Jaccard for Sequitur");
		
		System.out.println();
		System.out.println("LZ4");
		rankingslz4 = new ArrayList<ArrayList<EditDistance>>();
		EditDistance tempD;
		for(int i = 0; i < sequenceslz4.size(); i++) {
			System.out.print("\r" + (i+1) + "/" + sequences.size()  + "         ");
			rankingslz4.add(new ArrayList<EditDistance>());
			for(int z = 0; z < 5; z++) rankingslz4.get(i).add(new EditDistance(Integer.MAX_VALUE,0));

			for(int j = 0; j < i; j++) {
				tempD = new EditDistance(editDistance(sequenceslz4.get(i), sequenceslz4.get(j)),j);
				if(tempD.compareTo(rankingslz4.get(i).get(4)) < 0) {
					rankingslz4.get(i).set(4, new EditDistance(tempD));
					Collections.sort(rankingslz4.get(i));
				}
			}
			for(int j = i+1; j < sequenceslz4.size(); j++) {
				tempD = new EditDistance(editDistance(sequenceslz4.get(i), sequenceslz4.get(j)),j);
				if(tempD.compareTo(rankingslz4.get(i).get(4)) < 0) {
					rankingslz4.get(i).set(4, new EditDistance(tempD));
					Collections.sort(rankingslz4.get(i));
				}
			}
		}
		System.out.println("\nFinished Edit distance for LZ4");
		
		System.out.println();
		System.out.println("RLE");
		rankingsRLE = new ArrayList<ArrayList<EditDistance>>();
		for(int i = 0; i < sequencesRLE.size(); i++) {
			System.out.print("\r" + (i+1) + "/" + sequences.size()  + "         ");
			rankingsRLE.add(new ArrayList<EditDistance>());
			for(int z = 0; z < 5; z++) rankingsRLE.get(i).add(new EditDistance(Integer.MAX_VALUE,0));

			for(int j = 0; j < i; j++) {
				tempD = new EditDistance(editDistance(sequencesRLE.get(i), sequencesRLE.get(j)),j);
				if(tempD.compareTo(rankingsRLE.get(i).get(4)) < 0) {
					rankingsRLE.get(i).set(4, new EditDistance(tempD));
					Collections.sort(rankingsRLE.get(i));
				}
			}
			for(int j = i+1; j < sequencesRLE.size(); j++) {
				tempD = new EditDistance(editDistance(sequencesRLE.get(i), sequencesRLE.get(j)),j);
				if(tempD.compareTo(rankingsRLE.get(i).get(4)) < 0) {
					rankingsRLE.get(i).set(4, new EditDistance(tempD));
					Collections.sort(rankingsRLE.get(i));
				}
			}
		}
		
		System.out.println("\nFinished Edit distance for RLE");
		
		init();
		String s1, s2, s1s2, s2s1;
		double concat, indiv;
		rankingsCDMSequitur = new ArrayList<ArrayList<CDM>>();
		System.out.println();
		System.out.println("CDM Sequitur");
		CDM tempC;
		for(int i = 0; i < sequences.size(); i++) {
			System.out.print("\r" + (i+1) + "/" + sequences.size()  + "         ");
			rankingsCDMSequitur.add(new ArrayList<CDM>());
			for(int z = 0; z < 5; z++) rankingsCDMSequitur.get(i).add(new CDM(1,0));

			for(int j = 0; j < i; j++) {
				tempC = new CDM(cdmSequitur(sequences.get(i), sequences.get(j)),j);
				if(tempC.compareTo(rankingsCDMSequitur.get(i).get(4)) < 0) {
					rankingsCDMSequitur.get(i).set(4, new CDM(tempC));
					Collections.sort(rankingsCDMSequitur.get(i));
				}
			}
			for(int j = i+1; j < sequences.size(); j++) {
				tempC = new CDM(cdmSequitur(sequences.get(i), sequences.get(j)),j);
				if(tempC.compareTo(rankingsCDMSequitur.get(i).get(4)) < 0) {
					rankingsCDMSequitur.get(i).set(4, new CDM(tempC));
					Collections.sort(rankingsCDMSequitur.get(i));
				}
			}
		}
		System.out.println("\nFinished CMD for Sequitur");
		System.out.println();
		System.out.println("CDM Lz4");
		rankingsCDMLz4 = new ArrayList<ArrayList<CDM>>();
		for(int i = 0; i < sequences.size(); i++) {
			System.out.print("\r" + (i+1) + "/" + sequences.size()  + "         ");
			rankingsCDMLz4.add(new ArrayList<CDM>());
			for(int z = 0; z < 5; z++) rankingsCDMLz4.get(i).add(new CDM(1,0));

			for(int j = 0; j < i; j++) {
				tempC = new CDM(cdmLZ4(sequences.get(i), sequences.get(j)),j);
				if(tempC.compareTo(rankingsCDMLz4.get(i).get(4)) < 0) {
					rankingsCDMLz4.get(i).set(4, new CDM(tempC));
					Collections.sort(rankingsCDMLz4.get(i));
				}
			}
			for(int j = i+1; j < sequences.size(); j++) {
				tempC = new CDM(cdmLZ4(sequences.get(i), sequences.get(j)),j);
				if(tempC.compareTo(rankingsCDMLz4.get(i).get(4)) < 0) {
					rankingsCDMLz4.get(i).set(4, new CDM(tempC));
					Collections.sort(rankingsCDMLz4.get(i));
				}
			}
		}

		System.out.println("\nFinished CMD for LZ4");
		System.out.println();
		System.out.println("CDM RLE");
		rankingsCDMRLE = new ArrayList<ArrayList<CDM>>();
		for(int i = 0; i < sequences.size(); i++) {
			System.out.print("\r" + (i+1) + "/" + sequences.size()  + "         ");
			rankingsCDMRLE.add(new ArrayList<CDM>());
			for(int z = 0; z < 5; z++) rankingsCDMRLE.get(i).add(new CDM(Integer.MAX_VALUE,0));

			for(int j = 0; j < i; j++) {
				s1 = sequencesRLE.get(i);
				s2 = sequencesRLE.get(j);
				s1s2 = encode(sequences.get(i) + sequences.get(j));
				s2s1 = encode(sequences.get(j) + sequences.get(i));
				concat = (s1s2.length() + s2s1.length())/2.0;
				indiv = s1.length() + s2.length();
				
				tempC = new CDM(concat/indiv,j);
				if(tempC.compareTo(rankingsCDMRLE.get(i).get(4)) < 0) {
					rankingsCDMRLE.get(i).set(4, new CDM(tempC));
					Collections.sort(rankingsCDMRLE.get(i));
				}
			}
			for(int j = i+1; j < sequencesRLE.size(); j++) {
				s1 = sequencesRLE.get(i);
				s2 = sequencesRLE.get(j);
				s1s2 = encode(sequences.get(i) + sequences.get(j));
				s2s1 = encode(sequences.get(j) + sequences.get(i));
				concat = (s1s2.length() + s2s1.length())/2.0;
				indiv = s1.length() + s2.length();
				
				tempC = new CDM(concat/indiv,j);
				if(tempC.compareTo(rankingsCDMRLE.get(i).get(4)) < 0) {
					rankingsCDMRLE.get(i).set(4, new CDM(tempC));
					Collections.sort(rankingsCDMRLE.get(i));
				}
			}
		}
		System.out.println("\nFinished CMD for RLE");
		System.out.println();
		
		
		System.out.println("Matches Sequitur");
		matchesSequitur = new int[3];
		boolean[] matches = new boolean[3];
		int totalMatches;
		for(int i =0; i < sequences.size(); i++) {
			System.out.print("\r" + (i+1) + "/" + sequences.size()  + "         ");
			matches[0] = rankings.get(i).get(0).sequenceNumber == rankingsSequitur.get(i).get(0).sequenceNumber;
			matches[1] = rankings.get(i).get(2).sequenceNumber == rankingsSequitur.get(i).get(2).sequenceNumber;
			matches[2] = rankings.get(i).get(4).sequenceNumber == rankingsSequitur.get(i).get(4).sequenceNumber;
			totalMatches = 0;
			for(boolean b : matches) if(b) totalMatches++;
			if(totalMatches != 0)
				matchesSequitur[totalMatches-1]++;
		}
		System.out.println("\n");
		System.out.println("Matches Lz4");
		matcheslz4 = new int[3];
		for(int i =0; i < sequences.size(); i++) {
			System.out.print("\r" + (i+1) + "/" + sequences.size()  + "         ");
			matches[0] = rankings.get(i).get(0).sequenceNumber == rankingslz4.get(i).get(0).sequenceNumber;
			matches[1] = rankings.get(i).get(2).sequenceNumber == rankingslz4.get(i).get(2).sequenceNumber;
			matches[2] = rankings.get(i).get(4).sequenceNumber == rankingslz4.get(i).get(4).sequenceNumber;
			totalMatches = 0;
			for(boolean b : matches) if(b) totalMatches++;
			if(totalMatches != 0)
				matcheslz4[totalMatches-1]++;
		}
		System.out.println("\n");
		System.out.println("Matches RLE");
		matchesRLE = new int[3];
		for(int i =0; i < sequences.size(); i++) {
			System.out.print("\r" + (i+1) + "/" + sequences.size()  + "         ");
			matches[0] = rankings.get(i).get(0).sequenceNumber == rankingsRLE.get(i).get(0).sequenceNumber;
			matches[1] = rankings.get(i).get(2).sequenceNumber == rankingsRLE.get(i).get(2).sequenceNumber;
			matches[2] = rankings.get(i).get(4).sequenceNumber == rankingsRLE.get(i).get(4).sequenceNumber;
			totalMatches = 0;
			for(boolean b : matches) if(b) totalMatches++;
			if(totalMatches != 0)
				matchesRLE[totalMatches-1]++;
		}
		System.out.println("\n");
		System.out.println("Matches CDM Sequitur");
		matchesCDMSequitur = new int[3];
		for(int i =0; i < sequences.size(); i++) {
			System.out.print("\r" + (i+1) + "/" + sequences.size()  + "         ");
			matches[0] = rankings.get(i).get(0).sequenceNumber == rankingsCDMSequitur.get(i).get(0).sequenceNumber;
			matches[1] = rankings.get(i).get(2).sequenceNumber == rankingsCDMSequitur.get(i).get(2).sequenceNumber;
			matches[2] = rankings.get(i).get(4).sequenceNumber == rankingsCDMSequitur.get(i).get(4).sequenceNumber;
			totalMatches = 0;
			for(boolean b : matches) if(b) totalMatches++;
			if(totalMatches != 0)
				matchesCDMSequitur[totalMatches-1]++;
		}
		System.out.println("\n");
		System.out.println("Matches CDM LZ4");
		matchesCDMLz4 = new int[3];
		for(int i =0; i < sequences.size(); i++) {
			System.out.print("\r" + (i+1) + "/" + sequences.size()  + "         ");
			matches[0] = rankings.get(i).get(0).sequenceNumber == rankingsCDMLz4.get(i).get(0).sequenceNumber;
			matches[1] = rankings.get(i).get(2).sequenceNumber == rankingsCDMLz4.get(i).get(2).sequenceNumber;
			matches[2] = rankings.get(i).get(4).sequenceNumber == rankingsCDMLz4.get(i).get(4).sequenceNumber;
			totalMatches = 0;
			for(boolean b : matches) if(b) totalMatches++;
			if(totalMatches != 0)
				matchesCDMLz4[totalMatches-1]++;
		}
		System.out.println("\n");
		System.out.println("Matches CDM RLE");
		matchesCDMRLE = new int[3];
		for(int i =0; i < sequences.size(); i++) {
			System.out.print("\r" + (i+1) + "/" + sequences.size()  + "         ");
			matches[0] = rankings.get(i).get(0).sequenceNumber == rankingsCDMRLE.get(i).get(0).sequenceNumber;
			matches[1] = rankings.get(i).get(2).sequenceNumber == rankingsCDMRLE.get(i).get(2).sequenceNumber;
			matches[2] = rankings.get(i).get(4).sequenceNumber == rankingsCDMRLE.get(i).get(4).sequenceNumber;
			totalMatches = 0;
			for(boolean b : matches) if(b) totalMatches++;
			if(totalMatches != 0)
				matchesCDMRLE[totalMatches-1]++;
		}
		System.out.println("\nFinished distances");
	}	

	/**
	 * Outputs the results to an excel compatible doc
	 */
	public static void outputResults() {
		computeDistances();
		System.out.println("\nOutputting to excel");
		o.println(" , , ,Set S, , , ,Sequitur, , , , , ,Lz4, , , , , ,RLE, , , , , ,CDM Sequitur, , , , , ,CDM Lz4, , , , , ,CDM RLE");
		o.println("Sequence #, ,Match 1, Match 3, Match 5, ,Match 1, Match 3, Match 5, , , ,Match 1, Match 3, Match 5, , , ,Match 1, Match 3, Match 5, , , ,Match 1, Match 3, Match 5, , , ,Match 1, Match 3, Match 5, , , ,Match 1, Match 3, Match 5");
		for(int i = 0; i < 10; i++) {
			o.print((i+1) + ", ");
			for(int j = 0; j < 5; j+=2)
				o.print("," + (rankings.get(i).get(j).sequenceNumber+1));
			o.print(", ");
			
			for(int j = 0; j < 5; j+=2)
				o.print("," + (rankingsSequitur.get(i).get(j).sequenceNumber+1));
			if(i < 3) {
				o.print("," + (i+1) + "/3 match," + matchesSequitur[i] + ", ");
			} else {
				o.print(", , , ");
			}
			
			for(int j = 0; j < 5; j+=2)
				o.print("," + (rankingslz4.get(i).get(j).sequenceNumber+1));
			if(i < 3) {
				o.print("," + (i+1) + "/3 match," + matcheslz4[i] + ", ");
			} else {
				o.print(", , , ");
			}
			
			for(int j = 0; j < 5; j+=2)
				o.print("," + (rankingsRLE.get(i).get(j).sequenceNumber+1));
			if(i < 3) {
				o.print("," + (i+1) + "/3 match," + matchesRLE[i] + ", ");
			} else {
				o.print(", , , ");
			}
			
			for(int j = 0; j < 5; j+=2)
				o.print("," + (rankingsCDMSequitur.get(i).get(j).sequenceNumber+1));
			if(i < 3) {
				o.print("," + (i+1) + "/3 match," + matchesCDMSequitur[i] + ", ");
			} else {
				o.print(", , , ");
			}
			
			for(int j = 0; j < 5; j+=2)
				o.print("," + (rankingsCDMLz4.get(i).get(j).sequenceNumber+1));
			if(i < 3) {
				o.print("," + (i+1) + "/3 match," + matchesCDMLz4[i] + ", ");
			} else {
				o.print(", , , ");
			}
			
			for(int j = 0; j < 5; j+=2)
				o.print("," + (rankingsCDMRLE.get(i).get(j).sequenceNumber+1));
			if(i < 3) {
				o.print("," + (i+1) + "/3 match," + matchesCDMRLE[i] + ", ");
			} else {
				o.print(", , , ");
			}
			
			o.println();
		}
		o.println();
		
		/*o.println("Sequences:");
		for(int i = 0 ; i < sequences.size(); i++) {
			o.println("#" + (i+1) + "," + sequences.get(i).length() + "," + sequences.get(i));
			o.println("," + sequencesSequitur.get(i).length() + "," + sequencesSequitur.get(i));
			o.println("," + sequenceslz4.get(i).length + "," + sequenceslz4.get(i));
			o.println("," + sequencesRLE.get(i).length() + "," + sequencesRLE.get(i));
			o.println();
		}*/
		
		/*o.println("Local sequence alignment in set S:");
		
		for(int i = 1; i <= sequences.size(); i++) {
			o.print("," + i);
		}
		o.println();

		for(int i = 0; i < sequences.size(); i++) {
			o.print(i+1);
			for(int j = 0; j < sequences.size(); j++) {
				o.print("," + (distances.get(i).get(j).distance == Integer.MIN_VALUE ? "" : distances.get(i).get(j).distance) );
			}
			o.println();
		}
		o.println();		

		o.println("Jaccard Coefficients in Sequitur:");
		for(int i = 1; i <= sequences.size(); i++) {
			o.print("," + i);
		}
		o.println();
		
		for(int i = 0; i < sequences.size(); i++) {
			o.print(i+1);
			for(int j = 0; j < sequences.size(); j++) {
				o.print("," + (jaccardSequitur.get(i).get(j).jaccard == Integer.MIN_VALUE ? "" : jaccardSequitur.get(i).get(j).jaccard) );
			}
			o.println();
		}
		o.println();
		
		o.println("Edit distances in Lz4:");
		for(int i = 1; i <= sequenceslz4.size(); i++) {
			o.print("," + i);
		}
		o.println();
		
		for(int i = 0; i < sequenceslz4.size(); i++) {
			o.print(i+1);
			for(int j = 0; j < sequenceslz4.size(); j++) {
				o.print("," + (distanceslz4.get(i).get(j).distance == Integer.MAX_VALUE ? "" : distanceslz4.get(i).get(j).distance) );
			}
			o.println();
		}
		o.println();
		
		o.println("Edit distances in RLE:");
		for(int i = 1; i <= sequencesRLE.size(); i++) {
			o.print("," + i);
		}
		o.println();
		
		for(int i = 0; i < sequencesRLE.size(); i++) {
			o.print(i+1);
			for(int j = 0; j < sequencesRLE.size(); j++) {
				o.print("," + (distancesRLE.get(i).get(j).distance == Integer.MAX_VALUE ? "" : distancesRLE.get(i).get(j).distance) );
			}
			o.println();
		}
		o.println();
		
		o.println("CDM in Sequitur:");
		for(int i = 1; i <= sequences.size(); i++) {
			o.print("," + i);
		}
		o.println();
		
		for(int i = 0; i < sequences.size(); i++) {
			o.print(i+1);
			for(int j = 0; j < sequences.size(); j++) {
				o.print("," + (CDMSequitur.get(i).get(j).CDM == Integer.MAX_VALUE ? "" : CDMSequitur.get(i).get(j).CDM) );
			}
			o.println();
		}
		o.println();
		
		o.println("CDM in Lz4:");
		for(int i = 1; i <= sequences.size(); i++) {
			o.print("," + i);
		}
		o.println();
		
		for(int i = 0; i < sequences.size(); i++) {
			o.print(i+1);
			for(int j = 0; j < sequences.size(); j++) {
				o.print("," + (CDMLz4.get(i).get(j).CDM == Integer.MAX_VALUE ? "" : CDMLz4.get(i).get(j).CDM) );
			}
			o.println();
		}
		o.println();
		
		o.println("CDM in RLE:");
		for(int i = 1; i <= sequences.size(); i++) {
			o.print("," + i);
		}
		o.println();
		
		for(int i = 0; i < sequences.size(); i++) {
			o.print(i+1);
			for(int j = 0; j < sequences.size(); j++) {
				o.print("," + (CDMRLE.get(i).get(j).CDM == Integer.MAX_VALUE ? "" : CDMRLE.get(i).get(j).CDM) );
			}
			o.println();
		}
		
		o.println();*/
		System.out.println("Done!");
	}
	
	/**
	 * Computes the edit distance between two strings
	 * Given two strings s1 and s2, the edit distance between s1 and s2 is the minimum number of operations required to convert string s1 to s2. 
	 * 
	 * The following operations are allowed:
	 *  - Replacing one character of string by another character.
	 *  - Deleting a character from string
	 *  - Adding a character to string
	 * 
	 * @param s1 - String 1
	 * @param s2 - String 2
	 * @return - The edit distance between the two strings.
	 */
	public static int editDistance(String s1, String s2) {
	        int m = s1.length();
	        int n = s2.length();
	        int distances[][] = new int[m+1][n+1];
	        
	        for (int i = 0; i <= m; i++) {
	                distances[i][0] = i;
	        }
	        for (int j = 0; j <= n; j++) {
	                distances[0][j] = j;
	        }
	 
	        for (int i = 1; i <= m; i++) {
	                for (int j = 1; j <= n; j++) {
	                        if (s1.charAt(i-1) == s2.charAt(j-1))
	                        	distances[i][j] = distances[i-1][j-1];
	                        else 
	                        	distances[i][j] = 1 + Math.min(Math.min(distances[i][j-1],distances[i-1][j]),distances[i-1][j-1]);
	                }
	        }
	 
	        return distances[m][n];
	}

	/**
	 * Computes the edit distance between two byte arrays
	 * Given two strings s1 and s2, the edit distance between s1 and s2 is the minimum number of operations required to convert byte array s1 to s2. 
	 * 
	 * The following operations are allowed:
	 *  - Replacing one character of string by another character.
	 *  - Deleting a character from string
	 *  - Adding a character to string
	 * 
	 * @param s1 - LZ4 1 which contains a byte array
	 * @param s2 - LZ4 2 which contains a byte array
	 * @return - The edit distance between the two byte arrays.
	 */
	public static int editDistance(LZ4 s1, LZ4 s2) {
	        int m = s1.length;
	        int n = s2.length;
	        int distances[][] = new int[m+1][n+1];
	        
	        for (int i = 0; i <= m; i++) {
	                distances[i][0] = i;
	        }
	        for (int j = 0; j <= n; j++) {
	                distances[0][j] = j;
	        }
	 
	        for (int i = 1; i <= m; i++) {
	                for (int j = 1; j <= n; j++) {
	                        if (s1.data[i-1] == s2.data[j-1])
	                        	distances[i][j] = distances[i-1][j-1];
	                        else 
	                        	distances[i][j] = 1 + Math.min(Math.min(distances[i][j-1],distances[i-1][j]),distances[i-1][j-1]);
	                }
	        }
	 
	        return distances[m][n];
	}
	
	public static double jaccard(String s1, String s2) {
		HashSet<String> seq1 = new HashSet<String>(), seq2 = new HashSet<String>(), intersection, union;
		
		for(String s : s1.split(" ")) {
			if(Character.isDigit(s.charAt(0)))
				seq1.add(s);
		}
		for(String s : s2.split(" ")) {
			if(Character.isDigit(s.charAt(0)))
				seq2.add(s);
		}
		
		intersection = new HashSet<String>();
		for(String s : seq1) {
			for(String t : seq2) {
				if(Math.max(Sequitur.firstRule.similarity(s,t),Sequitur.firstRule.similarity(t, s)) >= .4) {
					intersection.add(s);
					break;
				}
			}
		}
		
		union = new HashSet<String>(seq1);
		union.addAll(seq2);
		
		return ((double)intersection.size())/union.size();
	}

	/**
	 * Computes the Local Sequence Alignment coefficient between two dna sequences
	 *  
	 * @param s1 - String 1
	 * @param s2 - String 2
	 * @return - The Local Sequence Alignment between the two dna sequences.
	 */
	public static int localSeqAlignment(String s1, String s2) {
		//System.out.println(s1);
		//System.out.println(s2);
		//System.out.println();
		s1 = " " + s1;
		s2 = " " + s2;
		int max = 0, h = 0 ;//maxI = 0, maxJ = 0;
		//String o1 = "", o2 = "";
		int[][] score = new int[s1.length()][s2.length()];
		int[][] pointers = new int[s1.length()][s2.length()];
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
					//maxI = i;
					//maxJ = j;
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
		/*
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
		
		System.out.println(new StringBuilder(o1).reverse().toString());
		System.out.println(new StringBuilder(o2).reverse().toString());
		System.out.println("Score: " + max);*/
		return max;
	}
	
	public static int match(char c1, char c2) {
		if(c1 == c2) return MATCHFACTOR;
		else return MISMATCHFACTOR;
		//return 2 - DNASUBSITUTIONMATRIX["ACGT".indexOf(c1)]["ACGT".indexOf(c1)];
	}
	
	public static double cdmLZ4(String s1, String s2) {
		int b1 = 0 , b2 = 0 , b1b2 = 0 , b2b1 = 0;
		try {
			b1 = LZ4.compress(s1).length;
			b2 = LZ4.compress(s2).length;
			b1b2 = LZ4.compress(s1 + s2).length;
			b2b1 = LZ4.compress(s1 + s2).length;
		} catch (UnsupportedEncodingException e) {
			e.printStackTrace();
		}
		double concat = (b1b2 + b2b1)/2.0;
		double indiv = b1 + b2;
		return concat/indiv;
	}
	
	public static double cdmSequitur(String a, String b) {
		init();
		compress(a);
		String s1 = Sequitur.firstRule.getRule0();
		init();
		compress(b);
		String s2 = Sequitur.firstRule.getRule0();
		init();
		
		ArrayList<String> cdm = new ArrayList<String>();
		compress(a);
		compress(b);
		Sequitur.firstRule.getRules(cdm, false);
		String s1s2 = cdm.get(0) + cdm.get(1);
		init();
		compress(b);
		compress(a);

		Sequitur.firstRule.getRules(cdm, false);
		String s2s1 = cdm.get(0) + cdm.get(1);
		double concat = (s1s2.length() + s2s1.length())/2.0;
		double indiv = s1.length() + s2.length();
		return(concat/indiv);
	}
}
