package compression_algorithms;

public class JaccardCoefficient implements Comparable<JaccardCoefficient>{
	double jaccard;
	int sequenceNumber;

	public JaccardCoefficient(double jaccard, int sequenceNumber) {
		this.jaccard = jaccard;
		this.sequenceNumber = sequenceNumber;
	}
	
	public JaccardCoefficient(JaccardCoefficient copy) {
		this.jaccard = copy.jaccard;
		this.sequenceNumber = copy.sequenceNumber;
	}

	@Override
	public int compareTo(JaccardCoefficient arg0) {
		if(jaccard - arg0.jaccard == 0) return 0;
		return (jaccard - arg0.jaccard > 0)?-1:1;
	}

}
