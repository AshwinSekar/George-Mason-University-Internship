package compression_algorithms;

public class CDM implements Comparable<CDM>{
	double CDM;
	int sequenceNumber;

	public CDM(double CDM, int sequenceNumber) {
		this.CDM = CDM;
		this.sequenceNumber = sequenceNumber;
	}
	
	public CDM(CDM copy) {
		this.CDM = new Double(copy.CDM);
		this.sequenceNumber = copy.sequenceNumber;
	}


	@Override
	public int compareTo(CDM arg0) {
		if(CDM - arg0.CDM == 0) return 0;
		return (CDM - arg0.CDM > 0)?1:-1;
	}
	
	@Override
	public String toString() {
		return sequenceNumber + ": " + CDM;
	}
}
