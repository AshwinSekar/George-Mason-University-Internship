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

import net.jpountz.lz4.LZ4;

import compression_algorithms.ValidateCA;

public class ClusterLz4 {
	static ArrayList<ArrayList<String>> clusters;
	static double avgSize = 0;
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String files = args[0];
		String dir = args[1];
		try {
			ClusterAll.o = new PrintWriter(new BufferedWriter(new FileWriter(files + "/Cluster_Lz4.csv")));
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
		for(double threshold = .1; threshold < .501; threshold += .02) {
			System.out.println("Threshold "  + threshold);
			System.out.println("Clustering");
			timeBefore = System.currentTimeMillis();
			try {
				cluster(threshold);
			} catch (UnsupportedEncodingException e) {
				e.printStackTrace();
			}
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
	 * @throws UnsupportedEncodingException 
	 */
	public static void cluster(double threshold) throws UnsupportedEncodingException {
		Random r = new Random();
		String representative;
		LZ4 representativeLz4;
		clusters = new ArrayList<ArrayList<String>>();
		int index;
		while(!ClusterAll.sequences.isEmpty()) {
			index = r.nextInt(ClusterAll.sequences.size());
			representative = ClusterAll.sequences.remove(index);
			representativeLz4 = ClusterAll.sequencesLz4.remove(index);
			clusters.add(new ArrayList<String>());
			clusters.get(clusters.size() - 1).add(representative);
			for(int i = 0; i < ClusterAll.sequences.size(); i++) {
				if(ValidateCA.editDistance(representativeLz4, ClusterAll.sequencesLz4.get(i)) < Math.max(representative.length(), ClusterAll.sequences.get(i).length())*(threshold)) {
					clusters.get(clusters.size() - 1).add(ClusterAll.sequences.remove(i));
					ClusterAll.sequencesLz4.remove(i);
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
