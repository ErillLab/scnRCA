
import java.util.Arrays;
import java.util.Vector;

public class CDS implements Comparable<CDS>
{
	private double			score;
	private boolean[] 		roundrobinRefset;
	
	// final 
	private final int[]		codonCounts;
	private final int[][]	posCounts;
	private final int		codonTotal;
	private final double	gc3content;
	private final double	gccontent;
	// base sequence; cdsName, cdsLocation, cdsDirection, cdsSequence, cdsFunction, cdsProtein
	private final String 	sequence;
	private Vector<String> 	extraInfo 	= new Vector<String>(7);
	private boolean 		isRibosomal = false;
	private int				location	= 0;
	
	/**
	 * constructor
	 * 
	 * @param seq				String sequences of bases
	 * @throws Exception		due to Util functions
	 */
	public CDS(String seq) throws Exception
	{
		this.sequence = seq;
		
		// fill codonCounts
		codonCounts = Util.calcTripletCounts(seq);
		
		int codonsum = 0;
		// sum codon total
		for(int i = 0; i < 64; i++)
			codonsum += codonCounts[i];
		codonTotal = codonsum;
		
		// fill positional codon counts (for RCA)
		posCounts = Util.calcPositionalCount(codonCounts);
		
		gc3content = (double)(posCounts[2][Util.BaseNumber("c")] + posCounts[2][Util.BaseNumber("g")]) / (double)codonTotal;
		gccontent  = (double) Util.GCsum(posCounts) / (double)codonTotal;

		score = 1;
		
	}// end constructor
	/**
	 * 
	 * @param roundrobinRef initial reference set participation
	 */
	public void setRoundRobin(boolean[] roundrobinRefset)
	{
		this.roundrobinRefset = roundrobinRefset;
	}


//---- getters
	
	/**
	 * @return	64 array of codon counts (shallow copy)
	 */
	public int[]	getCodonCounts()
	{
		return codonCounts;
	}
	
	/**
	 * @return	3x4 array positional base counts (shallow copy)
	 */
	public int[][]	getPositionalCounts()
	{
		return posCounts;
	}

	/**
	 * @return GC3 content of the CDS
	 */
	public double	getGC3content()
	{
		return gc3content;
	}
	/**
	 * @return GC content of the CDS
	 */
	public double	getGCcontent()
	{
		return gccontent;
	}
	
	/**
	 * @return	score (CAI/gCAI/RCA/gRCA)
	 */	
	public double	getScore()
	{
		return score;
	}

	/**
	 * @return codon total count
	 */
	public int		getCodonTotal()
	{
		return codonTotal;
	}
	/**
	 * @return shallow copy vector of CDS attributes (name, product, etc)
	 */
	public Vector<String> getInfo()
	{
		return extraInfo;
	}

	/**
	 * @return String of base sequence
	 */
	public String 	getSequence()
	{
		return sequence;
	}
	
	/**
	 * @return whether the CDS is in the reference set specified
	 */
	public boolean	getRefSetContributes(int roundrobinStep)
	{
		return this.roundrobinRefset[roundrobinStep];
	}
	/**
	 * @return whether the CDS is ribosomal
	 */
	public boolean	isRibosomal()
	{
		return isRibosomal;
	}
	public int		location()
	{
		return location;
	}
	
	/**
	 * 
	 * @return sum of final round robin reference sets the CDS is in
	 */
	public int		getTotalRefsetContributes()
	{
		int sum = 0;
		for(int i = 0; i < roundrobinRefset.length; i++)
			if(roundrobinRefset[i])
				sum++;
		return sum;
	}
	
	public String	getRefsetContributes()
	{
		String str = "[";
		for(int i = 0; i < roundrobinRefset.length; i++)
			str += (roundrobinRefset[i] ? 1 : 0) + " ";
		str += "]";
		return str;
	}

//----- setters
	/**
	 * @param	score	score (CAI/gCAI/RCA/gRCA/etc)
	 */
	public void setScore(double score)
	{
		this.score = score;
	}
	/**
	 * calculates and sets the score of the CDS based on the log of w values
	 * @param wValues_Log log values of the w values
	 */
	public void setScore(double[] wValues_Log)
	{
		int uncounted = 0;
		// sum scores of codons
		double temp_score = 0;
		
		// calculate scores based on (NOT stop, only-child codons)
		for(int i = 0; i < 64; i++)
		{
			// take out counts for stop codons and only-child codons
			switch(i){
			// stop (terminator) -- TAA, TAG, TGA
			case 48:
				uncounted += codonCounts[i];
				break;
			case 50:
				uncounted += codonCounts[i];
				break;
			case 56:
				uncounted += codonCounts[i];
				break;
			case 14:
				// met -- ATG
				uncounted += codonCounts[i];
				break;
			case 58:
				// trp -- TGG
				uncounted += codonCounts[i];
				break;
			default:
				temp_score += wValues_Log[i] * (double)codonCounts[i] ;
				break;
			}
		}

		temp_score = temp_score / (this.codonTotal - uncounted);
		temp_score = Math.exp(temp_score);
		
		this.score = temp_score;

	}

	/**
	 * mutator - adds attributes to CDS object
	 * 
	 * 
	 * @param LocusTag
	 * @param Name
	 * @param GeneSynonym
	 * @param Location
	 * @param Direction			+/- strand
	 * @param Function			function of protein product
	 * @param Protein
	 * @param Product
	 * @param Note
	 */
	public void addInfo(String LocusTag, String Name, String GeneSynonym, String Location, String Direction, String Function, String Protein, String Product, String Note)
	{
		// add passed attributes to extraInfo (vector of attributes)
		extraInfo.add(LocusTag);
		extraInfo.add(GeneSynonym);
		extraInfo.add(Name);
		extraInfo.add(Location);
		extraInfo.add(Direction);
		extraInfo.add(Function);
		extraInfo.add(Protein);
		extraInfo.add(Product);
		extraInfo.add(Note);
		
		// if the function mentions ribosomal, set the cds to be ribosomal
		if(Function.contains("Ribosom") || Product.contains("Ribosom") || Note.contains("Ribosom") || Function.contains("ribosom") || Product.contains("ribosom")|| Note.contains("ribosom"))
			isRibosomal = true;
	}
	/**
	 * 
	 * @param roundrobinStep	round robin step number
	 * @param shouldContribute	whether the CDS is in the reference set for that round robin step
	 */
	public void setRefSetContributes(int roundrobinStep, boolean shouldContribute)
	{
		this.roundrobinRefset[roundrobinStep] = shouldContribute;
	}

	public void	setLocation(int location)
	{
		this.location = location;
	}
	
//---- comparable
	/**
	 * tests equality of additional information / attributes. if none exist, return comparison of scores
	 * 
	 * @param	other	other CDS object to compare attributes / scores of
	 * @return	boolean does it equal?
	 */
	public boolean equals(CDS other)
	{
		// if no additional info, return comparison of scores
		if(extraInfo.isEmpty()) 
			return Arrays.equals(codonCounts, other.codonCounts);
		// else, compare the additional info!
		else
			return extraInfo.equals(other.getInfo());
	}

	/**
	 * compare scores of CDS and passed CDS
	 * 
	 * @param	other	CDS to compare to
	 * @return			comparison (-1 less than, 0 equal to, 1 more than)
	 */
	public int compareTo(CDS other)
	{
		return Double.compare(this.score, other.getScore());
	}

}
