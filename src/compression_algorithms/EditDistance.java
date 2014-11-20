package compression_algorithms;

public class EditDistance implements Comparable<EditDistance> {
	int distance;
	int sequenceNumber;

	public EditDistance(int distance, int sequenceNumber) {
		this.distance = distance;
		this.sequenceNumber = sequenceNumber;
	}
	
	public EditDistance(EditDistance copy) {
		this.distance = copy.distance;
		this.sequenceNumber = copy.sequenceNumber;
	}

	@Override
	public int compareTo(EditDistance arg0) {
		return distance - arg0.distance;
	}

}
