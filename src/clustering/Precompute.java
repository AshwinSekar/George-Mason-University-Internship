package clustering;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

public class Precompute {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String file = "60_Equal_112R.fasta";

		ArrayList<String> sequences = new ArrayList<String>();
		String sequence;
		String curLine = " ";

		try {
			BufferedReader f = new BufferedReader(new FileReader("LSHDIV_DataFiles/" + file));
			f.readLine(); // Skips past first header
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
			PrintWriter o = new PrintWriter(new File("PreComputed/" + file));
			long total = sequences.size();
			total = total * (total-1);
			total = total/ 2;
			long counter = 0;
			System.out.print("           ");
			for(int i = 0; i < sequences.size(); i++) {
				for(int j = i +1; j < sequences.size(); j++) {
					o.println(sequences.get(i) + "_" + sequences.get(j));
					o.println(ClusterAll.localSeqAlignmentSimilarity(sequences.get(i),sequences.get(j)));
					counter++;
					System.out.print("\r" + counter + "/" + total + " " + (double)counter/total + "%         ");
				}
			}
			System.out.println();
			f.close();
			o.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

}
