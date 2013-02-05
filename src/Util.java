import java.io.BufferedWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.RNATools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;


public class Util {
	
	/**
	 * denominator of RCA
	 * 
	 * @param posFreqs	positional frequencies
	 * @return			positional frequency product
	 */
	public static double[] calcPositionalCount_64(double[][] posFreqs){
		// fill positional codon counts
		double[] codonPosFreq64 = new double[64];

		int codonNum;

		// letters at each position (i j k --> 0123 = acgt)
		// pos 0
		for(int i = 0; i < 4; i++)
			// pos 1
			for(int j = 0; j < 4; j++)
				// pos 2
				for(int k = 0; k < 4; k++)
				{
					// calculate which codon
					codonNum = i*16 + j*4 + k;
					// product of positional frequencies of the 3 bases
					codonPosFreq64[codonNum] = posFreqs[0][i] * posFreqs[1][j] * posFreqs[2][k];
				}		
		// table of positional frequency product
		return codonPosFreq64;
	}
	/**
	 * 
	 * @param codonCounts
	 * @return the positional counts of the given sequence's codon counts
	 */
	public static int[][] calcPositionalCount(int[] codonCounts){
		// fill positional codon counts
		int[][] posCodonCounts = new int[3][4];

		int codonNum;
		int codonCount;
		// letters at each position (i j k --> 0123 = acgt)
		// pos 0
		for(int i = 0; i < 4; i++)
			// pos 1
			for(int j = 0; j < 4; j++)
				// pos 2
				for(int k = 0; k < 4; k++)
				{
					codonNum = i*16 + j*4 + k;
					codonCount = codonCounts[codonNum];
					
					// add to positional counts
					posCodonCounts[0][i] += codonCount;
					posCodonCounts[1][j] += codonCount;
					posCodonCounts[2][k] += codonCount;
				}
		return posCodonCounts;
	}


	/**
	 * 
	 * @param seq Sequence
	 * @return codon counts
	 * @throws Exception from TripletNumber()
	 */
	public static int[] calcTripletCounts(String seq) throws Exception
	{
		int[] freq = new int[64];

		// chop off extra letters
		int seqsize = seq.length();
		while(seqsize % 3 != 0)
			seqsize --;
		
		// for each triplet
		for(int i = 0; i < seqsize; i += 3)
		{
			String sTr = seq.substring(i, i + 3);
			int n = TripletNumber(sTr.toString());
			
			if(n >= 0)
				freq[n]+=1;
		}

		return freq;

	}

	
	/**
	 * author: Carbone
	 * 
	 * @param st 3 character string of bases
	 * @return codon number [0-63]
	 * @throws Exception from BaseNumber()
	 */
	public static int TripletNumber(String st) throws Exception
	{
		// put st's chars into char array cA
		char[] cA = new char[3];
		st.getChars(0,3,cA,0);
		
		// convert chars into base number
		return BaseNumber(new String(cA,0,1))*16+BaseNumber(new String(cA,1,1))*4+BaseNumber(new String(cA,2,1));
	}

	/**
	 * convert letter to base number
	 * 
	 * @param	cn	String
	 * @return		[0-3] int corresponding to base
	 */
	public static int BaseNumber(String cn) throws Exception
	{
		// normal bases (acgt)
		if (cn.equals("a")) return 0;
		else if (cn.equals("c")) return 1;
		else if (cn.equals("g")) return 2;
		else if (cn.equals("t")) return 3;

		// other letters
		else if (cn.equals("n")) return -100;
		else if (cn.equals("k")) return -100;
		else if (cn.equals("s")) return -100;
		else if (cn.equals("r")) return -100;
		else if (cn.equals("y")) return -100;
		else if (cn.equals("m")) return -100;
		else if (cn.equals("w")) return -100;
		
		// unrecognized symbol
		else
			throw new Exception("error in Util.BaseNumber: " + cn);
	}
//--------------------------------
	/**
	 * sum the total number of Gs and Cs in the seq
	 */
	public static int GCsum(int[][] posCounts) throws Exception
	{
		int sum = 0;
		// for each of the 3 positions
		for(int i = 0; i < 3; i ++)
		{
			sum += posCounts[i][BaseNumber("c")];
			sum += posCounts[i][BaseNumber("g")];
		}
		return sum;
	}

//--------------------------------
	/**
	 * writes size 64 array of double data in an 8x8 table to the passed BufferedWriter
	 * labeled with codon names (aaa etc)
	 * 
	 * @param writeConsole BufferedWriter to write data to
	 * @param data			data to write
	 * @throws IOException from BufferedWriter functions
	 */
	public static void print8x8String(BufferedWriter writeConsole, double[] data) throws IOException
	{
		// format
		DecimalFormat df = new DecimalFormat("0.000");
		// for each of the 64 elements
		for(int i = 0; i < 64; i ++)
		{
			// write the data element
			writeConsole.write(Util.TripletName(i) + ":" + df.format(data[i]) + "\t");
			// after every 8, add a newline character
			if(i%8 == 7)
				writeConsole.newLine();
		}
	}
	/**
	 * prints a 64 array of double data in an 8x8 table to the console
	 * labeled with codon names (aaa etc)
	 * 
	 * @param printMe	size 64 data to print
	 */
	public static void print8x8(double[] printMe)
	{
		// for each of the 64 elements
		for(int i = 0; i < 64; i ++)
		{
			// print the formatted element (with triplet label)
			System.out.format("%s:%.6f\t", Util.TripletName(i), printMe[i]);
			// after every 8 elements, print a newline
			if(i%8 == 7)
				System.out.println();
		}
	}
	/**
	 * prints a 64 array of integer data in an 8x8 table to the console
	 * labeled with codon names (aaa etc)
	 * 
	 * @param printMe	size 64 data to print
	 */
	public static void print8x8(int[] printMe)
	{
		// for each of the 64 elements
		for(int i = 0; i < 64; i ++)
		{
			// print the formatted element (with triplet label)
			System.out.format("%s:%d\t\t", Util.TripletName(i), printMe[i]);
			// after every 8 elements, print a newline
			if(i%8 == 7)
				System.out.println();
		}
	}
	
	/**
	 * prints a 3x4 array of double data into the console
	 * with base and position labels
	 * 
	 * @param printMe data to print
	 */
	public static void print3x4(double[][] printMe)
	{
		// for each of the 3 positions
		for(int i = 0; i < 3; i++)
		{
			// label the position
			System.out.print("Position " + i + ":\t");
			// for each of the bases
			for(int j = 0; j < 4; j++)
				// print the base data
				System.out.format("%.3f\t", printMe[i][j]);
			
			// new line
			System.out.println();
		}
	}
	
//--------------------------------
	/**
	 * Vector of Vectors of synonymous codons (int format) made
	 * 
	 * @return table of synonymous codons
	 * @throws IllegalAlphabetException from codonToAA()
	 * @throws IllegalSymbolException from codonToAA()
	 */
	public static Vector<Vector<Integer>> calcSynonCodonTable() throws IllegalAlphabetException, IllegalSymbolException
	{
		Vector<Vector<Integer>> aaSynonVV = new Vector<Vector<Integer>>();
		
		// make hashmap
		Map<String, Vector<Integer>> aaSynonMap = new HashMap<String, Vector<Integer>>();
		String aa;
		for(int i = 0; i < 64; i++)
		{
			aa = codonToAA(new Integer(i));
			
			// if aa DNA in aaSynonMap, make a new vector. else, grab the old
			Vector<Integer> codons;
			if(!aaSynonMap.containsKey(aa))
				codons = new Vector<Integer>();
			else
				codons = aaSynonMap.get(aa);
			codons.add(i);
			
			// put aa and associated vector of synonymous codons in the map
			aaSynonMap.put(aa, codons);
		}
		
		// convert map to vector of vectors of synonymous codons
		for(String eachKey: aaSynonMap.keySet())
			aaSynonVV.add(aaSynonMap.get(eachKey));
		
		return aaSynonVV;
	}
	
//--------------------------------
	/**
	 * generate double[] of max synonymous codon frequency
	 * 
	 * Vector containing Vectors of synonymous codons (ex. 0, 10, 26)
	 * loop through the families of synonymous codons
	 * 		find the max codon usage within the family
	 * 		place the max in maxCodonCount[] location for all synonymous codons
	 * return maxCodonCount[]
	 * 
	 * @param allCodonFreq	double[64] of codon usage
	 * @return				double[64] of maximum codon usage between synonyms
	 * @throws IllegalAlphabetException
	 * @throws IllegalSymbolException
	 */
	public static double[] maxCodonFreq(double[] allCodonFreq) throws IllegalAlphabetException, IllegalSymbolException
	{
		Vector<Vector<Integer>> aaSynonVV = calcSynonCodonTable();

		double[] maxCodonFreq = new double[64];
		
		// for each set of synonymous codons
		for(Vector<Integer> synonymousVector: aaSynonVV)
		{
			// find max codon usage
			double maxValue = 0;
			for(Integer codonLoc: synonymousVector)
			{
				double codonFreq = allCodonFreq[codonLoc];
				
				// new max found
				if(codonFreq > maxValue)
					// record max-thus-far of codon usage
					maxValue 	= codonFreq;
			}
			
			// put max count in each synonymous codon max count
			for(Integer codonLoc: synonymousVector)
				maxCodonFreq[codonLoc] = maxValue;
		}
		// return max codon frequency
		return maxCodonFreq;
	}

//--------------------------------

	
	
	//--------------------------------
	/**
	 * generate double[] of variance among synonymous codon frequencies
	 * 
	 * Vector containing Vectors of synonymous codons (ex. 0, 10, 26)
	 * loop through the families of synonymous codons
	 * 		find the max codon usage within the family
	 * 		place the max in maxCodonCount[] location for all synonymous codons
	 * return maxCodonCount[]
	 * 
	 * @param allCodonFreq	double[64] of codon usage
	 * @return				double[64] of maximum codon usage between synonyms
	 * @throws IllegalAlphabetException
	 * @throws IllegalSymbolException
	 */
	public static double[] varianceCodonFreq(double[] allCodonFreq) throws IllegalAlphabetException, IllegalSymbolException
	{
		Vector<Vector<Integer>> aaSynonVV = calcSynonCodonTable();

		double[] varianceCodonFreq = new double[64];
		
		// for each set of synonymous codons
		for(Vector<Integer> synonymousVector: aaSynonVV)
		{
			// find mean and squared mean (for variance) of codon usage on each set of syn codons
			double meansq = 0;
			double mean = 0;
			int counter=1;
			for(Integer codonLoc: synonymousVector)
			{
				double codonFreq = allCodonFreq[codonLoc];
				
				//compute mean
				mean=mean+codonFreq;
				//compute mean squared
				meansq=meansq+Math.pow(codonFreq, 2);
				counter=counter+1;
			}
			double variance=0;
			//assign variance
			variance=( (meansq/(double)counter) - Math.pow( (mean/(double)counter), 2) );
			
			// put variance in each synonymous codon variance count
			for(Integer codonLoc: synonymousVector)
				varianceCodonFreq[codonLoc] = variance;
		}
		// return variance vector
		return varianceCodonFreq;
	}

	//--------------------------------
	/**
	 * generate double[] of variance normalization factor for synonymous codons
	 * 
	 * Vector containing Vectors of synonymous codons (ex. 0, 10, 26)
	 * loop through the families of synonymous codons
	 * 		find the max codon usage within the family
	 * 		place the max in maxCodonCount[] location for all synonymous codons
	 * return maxCodonCount[]
	 * 
	 * @param allCodonFreq	double[64] of codon usage
	 * @return				double[64] of maximum codon usage between synonyms
	 * @throws IllegalAlphabetException
	 * @throws IllegalSymbolException
	 */
	public static double[] correct_variance(double[] uncorr_variance) throws IllegalAlphabetException, IllegalSymbolException
	{
		Vector<Vector<Integer>> aaSynonVV = calcSynonCodonTable();

		double[] corr_variance = new double[64];	//corrected variance
		
		// for each set of synonymous codons
		for(Vector<Integer> synonymousVector: aaSynonVV)
		{
			int counter=0;
			// count # of synonymous codons 
			for(Integer codonLoc: synonymousVector)
			{
				counter=counter+1;
			}
			
			//assign correction
			//(x-min)/(max-min) --- because min=0 always, it becomes x/max
			//times 0.5 and plus 0.5 set the normalized value on the range 0.5-1.0
			for(Integer codonLoc: synonymousVector)
			{
				switch (counter)
				{
				case 1: 
					corr_variance[codonLoc]=uncorr_variance[codonLoc]; 
					break;
				case 2: //divide by 0.5 (equals to multiply by 2)
					corr_variance[codonLoc]=uncorr_variance[codonLoc]*(double)2*0.5+0.5;
					break;
				case 3: //divide by 0.333 (multiply by 3)
					corr_variance[codonLoc]=uncorr_variance[codonLoc]*(double)3*0.5+0.5;
					break;
				case 4: //divide by 0.333 (multiply by 3)
					corr_variance[codonLoc]=uncorr_variance[codonLoc]*(double)3*0.5+0.5	;
					break;
				case 6: //divide by 0.3 (multiply by 3.33333)
					corr_variance[codonLoc]=uncorr_variance[codonLoc]*((double) 10 / (double) 3) *0.5+0.5;
					break;
				}
			}
		}
		
		// return corrected variance
		return corr_variance;
	}
	/**
	 * returns amino acid name of given codon
	 * 
	 * @param codon codon in integer format
	 * @return String amino acid name
	 * @throws IllegalAlphabetException
	 * @throws IllegalSymbolException
	 */
	public static String codonToAA(int codon) throws IllegalAlphabetException, IllegalSymbolException
	{
		// convert int --> DNA --> RNA --> translate to AA --> name of AA
		return RNATools.translate(DNATools.toRNA(DNATools.createDNA(Util.TripletName(codon)))).symbolAt(1).getName();
	}

	/**
	 * author: Carbone
	 * 
	 * convert [0-63] int name of codon to string
	 * @param n integer codon
	 * @return leter format of codon
	 */
	public static String TripletName(int n){
		// decompress integer --> letter information (int format)
		int i1 = (int)((float)n/16);
		int i2 = (int)((float)(n-i1*16)/4);
		int i3 = n-i1*16-i2*4;
		
		// add bases to return string		
		// return String name of codon
		return BaseName(i1) + BaseName(i2) + BaseName(i3);
	}
	
	/**
	 * 
	 * convert 0-3 integer name of nucleotide into string
	 * (ex: 0 --> "a")
	 * 
	 * @param	n	integer name of nucleotide
	 * @return		String name of nucleotide
	 */
	public static String BaseName(int n){
		String s = new String("");
		
		// integer to string
		switch(n){
		case 0: s = new String("a"); break;
		case 1: s = new String("c"); break;
		case 2: s = new String("g"); break;
		case 3: s = new String("t"); break;
		}
		
		return s;
	}
	

	
}
