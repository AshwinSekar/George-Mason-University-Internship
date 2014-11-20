package clustering;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Random;

import compression_algorithms.ValidateCA;

public class ClusterSequitur {
	static ArrayList<ArrayList<String>> clusters;
	static double avgSize = 0;
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		try {
			ClusterAll.o = new PrintWriter(new BufferedWriter(new FileWriter("Cluster_Sequitur.csv")));
			ClusterAll.o.println("File Name,Threshold,# Clu, Run Time (s), Singletons,Doubletons, W. Sim, Chao1, Shannon, Ace, avg Cluster size");
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		
		//File folder = new File("LSHDIV_DataFiles");
		//File[] listOfFiles = folder.listFiles(); 
		long timeBefore, timeAfter;

		//for (int i = 0; i < listOfFiles.length; i++) {
		//	if (listOfFiles[i].isFile()) {
		//		String files = listOfFiles[i].getName();
				String files = "returned_seqs_3.fa";
				System.out.println(files);
				ClusterAll.fileName = files;
				try {
					//ClusterAll.f = new BufferedReader(new FileReader("LSHDIV_DataFiles/" + files));
					ClusterAll.f = new BufferedReader(new FileReader("Simulated Data Sets/returned_seqs_3.fa"));
				} catch (IOException e) {
					e.printStackTrace();
				}
				ClusterAll.populateSequences(-1,false);
				for(double threshold = .6; threshold <= .8; threshold += .05) {
					System.out.println("Threshold "  + threshold);
					System.out.println("Clustering");
					System.out.print("           ");
					timeBefore = System.currentTimeMillis();
					try {
						cluster(threshold);
					} catch (UnsupportedEncodingException e) {
						e.printStackTrace();
					}
					timeAfter = System.currentTimeMillis();
					System.out.println((timeAfter - timeBefore)/100.0);
					System.out.println("Outputting");
					outputResults(files,timeAfter - timeBefore,threshold);
				}
		//	}	
		//}
		
		try {
			ClusterAll.f.close();
			ClusterAll.o.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
	
	/**
	 * Clusters the dna sequences in the ArrayList sequences
	 * @throws UnsupportedEncodingException 
	 */
	public static void cluster(double threshold) throws UnsupportedEncodingException {
		Random r = new Random();
		String representative, representativeSeq;
		clusters = new ArrayList<ArrayList<String>>();
		int index;
		while(!ClusterAll.sequences.isEmpty()) {
			index = r.nextInt(ClusterAll.sequences.size());
			representative = ClusterAll.sequences.remove(index);
			representativeSeq = ClusterAll.sequencesSequitur.remove(index);
			clusters.add(new ArrayList<String>());
			clusters.get(clusters.size() - 1).add(representative);
			for(int i = 0; i < ClusterAll.sequences.size(); i++) {
				System.out.println(ValidateCA.jaccard(representative, ClusterAll.sequences.get(i)));
				if(ValidateCA.jaccard(representativeSeq, ClusterAll.sequencesSequitur.get(i)) < threshold) {
					clusters.get(clusters.size() - 1).add(ClusterAll.sequences.remove(i));
					ClusterAll.sequencesSequitur.remove(i);
					i--;
					System.out.print("\r" + (double)ClusterAll.sequences.size()/ClusterAll.size + "%         ");
				}
			}
			avgSize += clusters.get(clusters.size() - 1).size();
		}
		avgSize = avgSize / clusters.size()	;
		System.out.println();
	}
	
	/**
	 * Outputs results of clustering
	 */
	public static void outputResults(String fileName,long time, double threshold) {
		ClusterAll.o.print(fileName +"," +threshold + "," + clusters.size() + "," + (time)/100.0 + ",");
		for(double d : ClusterAll.generateStatistics(clusters,false)) {
			ClusterAll.o.print(d + ",");
		}
		ClusterAll.o.println(avgSize);
		/*for(int i = 0; i < clusters.size(); i++) {
			ClusterAll.o.println("Cluster " + i + " size: " + clusters.get(i).size());
			for(int j = 0; j < clusters.get(i).size();j++) {
				ClusterAll.o.println(clusters.get(i).get(j));
			}
			ClusterAll.o.println();
		}*/
	}


}
