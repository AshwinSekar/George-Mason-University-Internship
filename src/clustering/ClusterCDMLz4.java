package clustering;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Random;

import compression_algorithms.ValidateCA;

public class ClusterCDMLz4 {
	static ArrayList<ArrayList<String>> clusters;
	static double avgSize = 0;
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String files = args[0];
		String dir = args[1];
		try {
			ClusterAll.o = new PrintWriter(new BufferedWriter(new FileWriter(files + "/Cluster_CDMLz4.csv")));
			ClusterAll.o.println("File Name,Threshold,# Clu, Run Time (s), Singletons,Doubletons, W. Sim, Chao1, Shannon, Ace, avg Cluster size");
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		long timeBefore, timeAfter;


		System.out.println(files);
		ClusterAll.fileName = files;
		try {
			ClusterAll.f = new BufferedReader(new FileReader(dir + "/" + files));
		} catch (IOException e) {
			e.printStackTrace();
		}
		ClusterAll.populateSequences(-1,true);
		for(double threshold = .6; threshold < 1.01; threshold += .02) {
			System.out.println("Threshold "  + threshold);
			System.out.println("Clustering");
			timeBefore = System.currentTimeMillis();
			cluster(threshold);
			timeAfter = System.currentTimeMillis();
			System.out.println((timeAfter - timeBefore)/1000.0 + " secs");
			System.out.println("Outputting");
			outputResults(files,timeAfter - timeBefore,threshold);
		}

		try {
			ClusterAll.f.close();
			ClusterAll.o.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	/**
	 * Clusters the dna sequences in the ArrayList sequences
	 */
	public static void cluster(double threshold) {
		Random r = new Random();
		String representative;
		clusters = new ArrayList<ArrayList<String>>();
		while(!ClusterAll.sequences.isEmpty()) {
			representative = ClusterAll.sequences.remove(r.nextInt(ClusterAll.sequences.size()));
			clusters.add(new ArrayList<String>());
			clusters.get(clusters.size() - 1).add(representative);
			for(int i = 0; i < ClusterAll.sequences.size(); i++) {
				if(ValidateCA.cdmLZ4(representative, ClusterAll.sequences.get(i)) < threshold) {
					clusters.get(clusters.size() - 1).add(ClusterAll.sequences.remove(i));
					i--;
					//System.out.print("\r" + (double)ClusterAll.sequences.size()/ClusterAll.size + "%         ");
				}
			}
			avgSize += clusters.get(clusters.size() - 1).size();
		}
		avgSize = avgSize / clusters.size()	;
	}

	/**
	 * Outputs results of clustering
	 */
	public static void outputResults(String fileName,long time, double threshold) {
		ClusterAll.o.print(fileName +"," +threshold + "," + clusters.size() + "," + (time)/1000.0 + ",");
		for(double d : ClusterAll.generateStatistics(clusters,true)) {
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
