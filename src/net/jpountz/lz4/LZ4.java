package net.jpountz.lz4;

import java.io.UnsupportedEncodingException;

public class LZ4 {
	static LZ4Factory factory = LZ4Factory.fastestInstance();	
	static LZ4Compressor compressor = factory.fastCompressor();
	public int length;
	public byte[] data;
	
	public LZ4(int length, byte[] data) {
		super();
		this.length = length;
		this.data = data;
	}
	
	@Override
	public String toString() {
		String ans = "";
		for(int i = 0 ; i < length; i++) {
			ans += data[i] + " ";
		}
		return ans;
	}
	
	public static LZ4 compress(String s) throws UnsupportedEncodingException {
		int maxCompressedLength;
		int length;
		byte[] data;
		
		maxCompressedLength = compressor.maxCompressedLength(s.length());
		data = new byte[maxCompressedLength];
		length = compressor.compress(s.getBytes("UTF-8"), 0, s.length(), data, 0);
		return new LZ4(length,data);
	}
	
	public static void main(String[] args) throws UnsupportedEncodingException {
		System.out.println(compress("CAGACTTGACATTGAAGGAACCTGCTAGAAATAGCGGGGTGCTCTTCGGAGCCCATGAAA"));
	}
}
