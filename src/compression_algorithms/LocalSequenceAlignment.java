package compression_algorithms;

public class LocalSequenceAlignment implements Comparable<LocalSequenceAlignment>{
	int distance;
	int sequenceNumber;

	public LocalSequenceAlignment(Double double1, int sequenceNumber) {
		this.distance = double1.intValue();
		this.sequenceNumber = sequenceNumber;
	}
	
	public LocalSequenceAlignment(LocalSequenceAlignment copy) {
		this.distance = copy.distance;
		this.sequenceNumber = copy.sequenceNumber;
	}

	@Override
	public int compareTo(LocalSequenceAlignment arg0) {
		return (distance > arg0.distance)?-1:1;
	}
}
