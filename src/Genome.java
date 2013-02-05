

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.Vector;

import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;


public class Genome {
	// Vector of CDSs comprising the genome
	private Vector<CDS> genome = new Vector<CDS>();
	//
	private int			genomeSize;
	private int			global_codonCountTotal;

	// used each round robin round
	private int					refsetRRIndex;
	// used each iteration
	private BufferedWriter		iterationWriter;
	private Vector<double[]> 	wValStore;

	
	//---- STATISTICS of Reference Set
	private int			refsetCDStotal;
	private	int			codonCountTotal;
	// reference set codon count / positional counts / total codon counts
	private int[]		codonCounts			= new int[64];
	private int[][]		codonPosCounts		= new int[3][4];
	// frequency of codon in reference set
	private double[]	codonFreq			= new double[64];
	private double[]	codonFreqMAX		= new double[64];
	private double[][]	codonPosFreq		= new double[3][4];
	// added by Erill?
	private double[]    rcaMAX				= new double[64];
	private double[]	RCAvar				= new double[64];
	private double[]	corrRCAvar			= new double[64];

	//---- DERIVED FIGURES
	// RCC - number of CDSs in reference set which contain that codon
	private	int[]		codonRefSetCDScount	= new int[64];
	private double[]	rcc					= new double[64];
	// W (preliminary and final - g, etc applied)
	private double[]	prelimW				= new double[64];
	private double[]	wValues				= new double[64];
	
	// CRITERIONS	
	private double		scoreCriterion_ribosomal;
	private double 		strengthCriterion_Freq;
	private double 		contentCriterion;
	private double		gc3RefsetAvg;
	private double		gcRefsetAvg;
	private double		aaCriterion;
	// used for the strength criterion
	private int[]		global_codonCounts	= new int[64];
	private double[]	global_codonFreq	= new double[64];
	

	//----- PARAMETERS
	private char		whichIndex;
	private boolean		useRCC;
	//
	private String		folderOut	= "";
	
	
	/**
	 * constructor
	 * 
	 * @param genome 		Vector of CDSs comprising the genome
	 * @param whichIndex 	which index to use
	 * @param useRCC		whether to use rare codon correction
	 * @param organism		name of the organism
	 * @throws Exception	from FileWriter and calcStats()
	 */
	public Genome(Vector<CDS> genome, char whichIndex, boolean useRCC, String organism, int rrSize) throws Exception
	{
		this.whichIndex	= whichIndex;
		this.useRCC 		= useRCC;
		this.refsetRRIndex = 0;

		this.genome 	= genome;
		this.genomeSize = genome.size();

		
		// create output folder, iteration file writer
		(new File(folderOut)).mkdir();
		this.folderOut += "/" + organism + "_/";
		
		for(int i = 0; i < 64; i++)
			wValues[i] = 1;
		
		calcStats();
		
	}


	public double aaCrit() throws IllegalAlphabetException, IllegalSymbolException
	{
		calcCodonCounts();
		
		return aaCriterion();
	}
	public void printAACrit() throws IllegalAlphabetException, IllegalSymbolException
	{
		System.out.println("AA Crit:\t" + aaCriterion());
	}
	
	
// ---------------- CALCULATE REFERENCE SET STATS
	/**
	 * controls calculation of statistics for the reference set
	 * 
	 * calculate codon frequency, positional frequency, rare codon correction values, # cds contained.
	 * Additional Parameters Applied Here (ex g)
	 */
	private void	calcStats() throws Exception
	{
		// count codons --> codonCounts, codonCountTotal, refsetContains (for g), positional codon counts
		calcCodonCounts();
		calcCodonFrequency();
		calcPositionalFrequency();
		// max frequency
		codonFreqMAX = Util.maxCodonFreq(codonFreq);

		// build prelim W
		if(whichIndex == 'c')
			calcW_CAI();
		else if(whichIndex == 'r' /*PLUS SCN modifier*/)
			//calcW_nRCA_VAR();
			calcW_nRCA();
			//calcW_RCA();
		else
			throw new Exception("REFERENCE SET: passed whichIndex is buggy");
		
		// build final W
		calcRCC();
		// update w values
		if(useRCC)
			for(int i = 0; i < 64; i++)
				wValues[i] = rcc[i] * prelimW[i];
		else
			for(int i = 0; i < 64; i++)
				wValues[i] = prelimW[i];
		
//		printW();
	}
	/**
	 * calculate reference set codon counts, total number of codons in reference set,
	 * total number of CDSs in reference set which contain each codon, positional counts
	 */
	private void	calcCodonCounts()
	{
		// reinitialize values
		refsetCDStotal 		= 0;
		codonCountTotal 	= 0;
		global_codonCountTotal	= 0;
		for(int i = 0; i < 64; i++)
		{
			codonCounts[i] = 0;
			codonRefSetCDScount[i] = 0;
			global_codonCounts[i] = 0;
		}
		
		// for each CDS in the reference set
		int[] eachcdsCount;
		int cdsCodonCount;
		
		for(CDS eachCDS: this.genome)
		{
			eachcdsCount 	= eachCDS.getCodonCounts();
			global_codonCountTotal += eachCDS.getCodonTotal();
			
			//---- COUNT REFERENCE SET CODONS
			if(eachCDS.getRefSetContributes(refsetRRIndex))
			{
				// add to refset cds total count
				refsetCDStotal ++;
				// add each CDS's codon total to the ref set codon total
				codonCountTotal += eachCDS.getCodonTotal();
				
				// for each codon, if it's contained in the CDS
					// add the count to the reference set codon counts, and one count to the reference set cds counts
				for(int i = 0; i < 64; i++)
				{
					cdsCodonCount 	 = eachcdsCount[i];
					if(cdsCodonCount != 0)
					{
						codonCounts[i]	+= cdsCodonCount;
						codonRefSetCDScount[i] ++;
					} //end if
				}// end for
			}// end if
			
			// ---- COUNT GENOME CODONS
			for(int i = 0; i < 64; i++)
			{
				cdsCodonCount 		   = eachcdsCount[i];
				global_codonCounts[i] += cdsCodonCount;
			}// end for
		}// end for
		
		
		
		codonPosCounts = Util.calcPositionalCount(codonCounts);
	}
	/**
	 * calculate codon frequencies in the reference set
	 * (codon counts in reference set / total codon count in reference set)
	 * 
	 * @note uses Laplace's Rule of Succession
	 */
	private void 	calcCodonFrequency()
	{
		// for each codon, codon frequency = codon counts in reference set / codon count total
		for(int i = 0; i < 64; i ++)
		{
			codonFreq[i] 		= ((double) codonCounts[i] 			+ (2d / 59d)) / ((double) codonCountTotal 		 + 2d);
			global_codonFreq[i] = ((double) global_codonCounts[i] 	+ (2d / 59d)) / ((double) global_codonCountTotal + 2d);
		}
	}
	/**
	 * calculate positional frequency in the reference set
	 * (positional counts in ref set / total codon count in reference set)
	 * 
	 * @note uses Laplace's rule of succession
	 */
	private void 	calcPositionalFrequency()
	{
		// for each position
		for(int i = 0; i < 3; i ++)
			// for each base
			for(int j = 0; j < 4; j ++)
				// calculate positional frequency (using Laplace's rule of succession)
				codonPosFreq[i][j] = ((double)codonPosCounts[i][j] + (2d/4d)) / ((double)codonCountTotal + 2d);

//		System.out.println("POSITIONAL FREQ");
//		Util.print3x4(codonPosFreq);
	}
	
	/**
	 * calculate rare codon correction value
	 * (as seen in carbone's method = number of CDSs in reference set that contain the codon / reference set size in CDSs)
	 */
	private void 	calcRCC()
	{
		for(int i = 0; i < 64; i++)
			rcc[i] = ((double)codonRefSetCDScount[i] + 1d) / ((double)refsetCDStotal + 2d);
	}

	/**
	 * calculate w value (codon frequency in reference set / max codon frequency)
	 */
	private void 	calcW_CAI()
	{
		// for each of 64 codons
		for(int i = 0; i < 64; i++)
			// calculate the w values
			prelimW[i] = codonFreq[i] / codonFreqMAX[i];
	}
	/**
	 * calculate w value (codon frequency in reference set / max codon frequency (position based))
	 */
	private void 	calcW_RCA()
	{
		// convert 3x4 positional frequency matrix --> 64x1
			// fx1 * fy2 * fz3
		double[] codonPosFreq64 = Util.calcPositionalCount_64(codonPosFreq);
		
		for(int i = 0; i < 64; i++)
			prelimW[i] = codonFreq[i] / codonPosFreq64[i];
	}
	
	
	/**
	 * calculate w value (codon frequency in reference set / max codon frequency (position based))
	 * @throws IllegalSymbolException 
	 * @throws IllegalAlphabetException 
	 */
	private void 	calcW_nRCA_VAR() throws IllegalAlphabetException, IllegalSymbolException
	{
		// convert 3x4 positional frequency matrix --> 64x1
			// fx1 * fy2 * fz3
		double[] codonPosFreq64 = Util.calcPositionalCount_64(codonPosFreq);

		
		for(int i = 0; i < 64; i++)
			prelimW[i] = codonFreq[i] / codonPosFreq64[i];
		
		// normalize RCA by amino acid (like in CAI)
		
		//get maximum RCAxyz codon for every amino acid (and assign it to each codon)
		//using the maxCodonFreq method used by CAI, but giving it the prelim RCAxyz values
		//instead of the codcon frequencies the method originally used
		rcaMAX = Util.maxCodonFreq(prelimW);					//get max RCA for AA

		///*** TAKEN OUT BY MINDY TO TRY SOMETHING
		
		//normalize RCA
		// for each of 64 codons, divide the RCAxyz value by the max RCAxyz value obtained
		for(int i = 0; i < 64; i++)
		{
			// calculate the w values
			//prelimW[i] = prelimW[i] / rcaMAX[i];
			prelimW[i] = prelimW[i] / rcaMAX[i] ;
		}
		
		//VARnRCA computation
		//compute variances and normalize them (for VARnRCA)
		RCAvar = Util.varianceCodonFreq(prelimW);		//get variance for AA
		corrRCAvar=Util.correct_variance(RCAvar);		//correct variance for AA
		
		// for each of 64 codons
		for(int i = 0; i < 64; i++)
		{			
			//weight the codons according to their variance (VARnRCA)
			//but adjusted to max possible  variance in each codon class
			prelimW[i] = prelimW[i]* corrRCAvar[i];
		}

	}
	
	/**
	 * calculate w value
	 * 		RCA = (codon frequency in reference set / max codon frequency (position based)
	 * 		nRCA = RCA / max synonymous RCA
	 * @throws IllegalAlphabetException
	 * @throws IllegalSymbolException
	 */
	private void	calcW_nRCA() throws IllegalAlphabetException, IllegalSymbolException
	{
		// convert 3x4 positional frequency matrix --> 64x1
		// fx1 * fy2 * fz3
		double[] codonPosFreq64 = Util.calcPositionalCount_64(codonPosFreq);
	
		for(int i = 0; i < 64; i++)
			prelimW[i] = codonFreq[i] / codonPosFreq64[i];
		
		rcaMAX = Util.maxCodonFreq(prelimW);
		
		for(int j = 0; j < 64; j++)
			prelimW[j] = prelimW[j] / rcaMAX[j];
			
	}
	
	/// set w values (porting over a new starting reference set
	public void		setW(double[] wValues)
	{
		this.wValues = wValues;
	}
	
	
	
	

// ---------------- SCORE CDSs
	/**
	 * set CDS scores based on the w values
	 */
	public void 	scoreCDS() throws Exception
	{
		// take the logged form of all the w values
		double[] wValues_Log = new double[64];
		for(int i = 0; i < 64; i++)
			wValues_Log[i] = Math.log(wValues[i]);
//			wValues_Log[i] = codonFreq[i] * Math.log(wValues[i]);
		
		// pass the logged forms to each CDS to calculate their score
		for(CDS eachCDS: genome)
			eachCDS.setScore(wValues_Log);
	}
	
// ---------------- RECALCULATE REFERENCE SET	
	/**
	 * sorts genome by index, resets the reference set to top #cutoffNumCDSs CDSs
	 * write in iterationWriter
	 * 
	 * @param cutoffNumCDSs number of CDSs to be in reference set
	 * @throws IOException	from iteration writer
	 */
	public void		recutRefSet(int cutoffNumCDSs) throws IOException
	{
		// select reference set
		Collections.sort(genome, Collections.reverseOrder());
		
		// for each CDS, reset whether it is in the reference set
		CDS cds;
		for(int i = 0; i < genome.size(); i++)
		{
			cds = genome.get(i);
			// the first set of CDSs are 
			if(i < cutoffNumCDSs)
				cds.setRefSetContributes(refsetRRIndex, true);
			else
				cds.setRefSetContributes(refsetRRIndex, false);

		}
		
		// update iteration writer
		iterationWriter.write(cutoffNumCDSs + " CDSs");
		iterationWriter.newLine();
	}
	/**
	 * calculate the statistics for the reference set (counts, etc)
	 * record w values (to determine whether reference set is stable/oscillating)
	 * write in iteration writer
	 * 
	 * @param iterationCount
	 * @throws Exception from calcRefSetStats and iteration writer
	 */
	public void 	recalcRefSet(int iterationCount) throws Exception
	{			
		double[] oldW = new double[64];
		for(int i = 0; i < 64; i++)
			oldW[i] = wValues[i];
		wValStore.add(oldW);
		
		// calculate stats
		calcStats();
//		Util.print8x8(wValues);
		/*
		System.out.println("Nonseedy");
		Util.print8x8(codonCounts);
		System.out.println();
		*/
		
		// print stats
		iterationWriter.write("Iteration#: " + iterationCount);
		iterationWriter.newLine();
		printf_Stats();
	}
	public void 	recalcRefset_simple() throws Exception
	{
		calcStats();
		
		/*
		System.out.println("SEEDY:");
		Util.print8x8(codonCounts);
//		Util.print8x8(wValues);
		System.out.println();
		*/
	}

		
// ----------------- calculate + PRINT REFERENCE SET STATS (to console and files)
	/**
	 * plus refset GC + GC3 content
	 * content criterion	:
	 * 		sum : for each CDS in genome, difference of score from average * difference of GC3 content from average
	 * 		/ (# cds in genome - 1) * standard deviation of scores * standard deviation of GC3 content
	 * 
	 * ribosomal criterion	= 
	 * 		sum : for each ribosomal CDS in genome, 
	 * 			difference of log(score) from genome (log) average / standard deviation of genome CDS (log) scores
	 * 		/ # ribosomal CDSs
	 * strength criterion	=
	 * 		sum: for each reference set CDS in genome,
	 * 			difference of log(score) from genome (log)average / standard deviation of genome CDS (log)scores
	 * 		/ # reference set CDSs
	 * 
	 * technically this code could be a little faster... I could forget calculating genome log score avg, 
	 * and just make it log score sum... but it is more readable this way
	 */
	private void	calcCriterion() throws Exception
	{
		
		// ------------ calculate arithmetic mean for genome CDS scores and GC3 content
		double	genomeScoreAvg		= 0.0d;
		double	genomeLogScoreAvg	= 0.0d;
		double	gc3contentAvg		= 0.0d;
		for(CDS eachCDS: genome)
		{
			genomeScoreAvg		+= eachCDS.getScore();
			genomeLogScoreAvg	+= Math.log(eachCDS.getScore());
			gc3contentAvg		+= eachCDS.getGC3content();
		}
		genomeScoreAvg		/= (double)genomeSize;
		genomeLogScoreAvg	/= (double)genomeSize;
		gc3contentAvg		/= (double)genomeSize;
		
		// ------------- calculate standard deviation of scores and GC3 content
		double	genomeScoreStDEV	= 0.0d;
		double	genomeLogScoreStDEV	= 0.0d;
		double	gc3contentStDEV		= 0.0d;
		
		for(CDS eachCDS: genome)
		{
			genomeScoreStDEV	+= Math.pow(genomeScoreAvg		- eachCDS.getScore(),			2);
			genomeLogScoreStDEV	+= Math.pow(genomeLogScoreAvg	- Math.log(eachCDS.getScore()),	2);
			gc3contentStDEV		+= Math.pow(gc3contentAvg		- eachCDS.getGC3content(),		2);			
		}
		genomeScoreStDEV 	= Math.sqrt(genomeScoreStDEV	/ (double) (genomeSize -1));
		genomeLogScoreStDEV = Math.sqrt(genomeLogScoreStDEV	/ (double) (genomeSize -1));
		gc3contentStDEV		= Math.sqrt(gc3contentStDEV		/ (double) (genomeSize -1));
		
		// ------------- calculate Ribosomal Criterion, Content Criterion
		int		ribosomeCDSsize 	= 0;
		int		refsetSize 			= 0;
		scoreCriterion_ribosomal 	= 0.0d;
		strengthCriterion_Freq 		= 0.0d;
		contentCriterion 			= 0.0d;
		gc3RefsetAvg				= 0.0d;
		gcRefsetAvg					= 0.0d;
		for(CDS eachCDS: genome)
		{
			if(eachCDS.isRibosomal())
			{
				scoreCriterion_ribosomal 	+= (Math.log(eachCDS.getScore()) - genomeLogScoreAvg);
				ribosomeCDSsize ++;
			}
			if(eachCDS.getRefSetContributes(refsetRRIndex))
			{
				strengthCriterion_Freq 		+= (Math.log(eachCDS.getScore()) - genomeLogScoreAvg);
				refsetSize ++;
				gc3RefsetAvg += eachCDS.getGC3content();
				gc3RefsetAvg += eachCDS.getGCcontent();
			}
			
			contentCriterion	+= (eachCDS.getScore() - genomeScoreAvg) * (eachCDS.getGC3content() - gc3contentAvg);
		}
		scoreCriterion_ribosomal	/= genomeLogScoreStDEV;
		scoreCriterion_ribosomal	/= (double) ribosomeCDSsize;
		strengthCriterion_Freq		/= genomeLogScoreStDEV;
		strengthCriterion_Freq 		/= (double) refsetSize;
		contentCriterion			/= ((double)(genomeSize - 1) * genomeScoreStDEV * gc3contentStDEV);
		gc3RefsetAvg				/= (double) refsetSize ++;
		aaCriterion					= aaCriterion();
		
	}
	/**
	 * aa criterion - euclidean distance of aa usage RR vs Genome
	 * @throws IllegalSymbolException 
	 * @throws IllegalAlphabetException 
	 */
	private double aaCriterion() throws IllegalAlphabetException, IllegalSymbolException
	{
		double distance = 0;

		double refsetAA;
		double genomeAA;
		
		Vector<Vector<Integer>> aaTable = Util.calcSynonCodonTable();
		for(Vector<Integer> aa: aaTable)
		{
			refsetAA = 0;
			genomeAA = 0;
			
			for(Integer codon: aa)
			{
//				refsetAA += codonCounts[codon];
//				genomeAA += global_codonCounts[codon];
				refsetAA += codonFreq[codon];
				genomeAA += global_codonFreq[codon];
			}
			
			distance += Math.pow((double)(refsetAA - genomeAA), 2);
		}
		
		distance = Math.sqrt(distance);
		
		return distance;
	}
	
// ----------------- PRINT REF SET STATS
	// to console
	/**
	 * print reference set codon counts to console
	 */
	public void 	printCodonCounts()
	{
		System.out.println("Codon Counts:");
		Util.print8x8(codonCounts);
	}
	/**
	 * print reference set positional frequencies to console
	 */
	public void		printPositionalFreq()
	{
		System.out.println("Positional Frequences:");
		Util.print3x4(codonPosFreq);
	}
	/**
	 * print reference set codon frequencies to console
	 */
	public void 	printCodonFreq()
	{
		System.out.println("Codon Frequencies:");
		Util.print8x8(codonFreq);
	}
	/**
	 * print reference set max synonymous codon frequencies
	 */
	public void		printCodonFreqMAX()
	{
		System.out.println("Max Synonymous Codon Frequencies:");
		Util.print8x8(codonFreqMAX);
	}
	/**
	 * print rare codon correction values
	 */
	public void		printRCC()
	{
		System.out.println("Rare Codon Correction: ");
		Util.print8x8(rcc);
	}
	/**
	 * print w values (without g)
	 */
	public void		printPrelimW()
	{
		System.out.println("Preliminary W Values:");
		Util.print8x8(prelimW);
	}
	/**
	 * print iteration's w values (may include g)
	 */
	public void 	printW()
	{
		System.out.println("W Values:");
		Util.print8x8(wValues);
	}
	/**
	 * print CDS scores, separating reference set scores from non reference set
	 */
	public void		printCDSrefsetSCORES()
	{
		// refset scores header
		System.out.println("Reference Set scores:");

		// for each CDS in the genome
		boolean inRefSet = true;
		for(CDS eachCDS: genome)
		{
			// genome scores header
			if(inRefSet && !eachCDS.getRefSetContributes(refsetRRIndex))
			{
				inRefSet = false;
				System.out.println("\nrest of genome scores:");
			}
			
			// print scores
			System.out.format("%.3f ", eachCDS.getScore());
		}
		System.out.println();
	}
	/**
	 * print reference set CDS information, 1 line each
	 */
	public void		printCDSrefsetINFO()
	{
		
		System.out.println("Reference Set CDS information");
		
		// for each cds in the reference set
		for(CDS eachCDS: genome)
			if(eachCDS.getRefSetContributes(refsetRRIndex))
			{
				// print each piece of information, separated by tab
				for(String eachInfoPiece: eachCDS.getInfo())
					System.out.print(eachInfoPiece+ "\t");
				// new line
				System.out.println();
			}
		System.out.println();
	}
	/**
	 * print criterions (ribosomal, strength, content) to console
	 */
	public void		printCriterion()
	{
		// print ribosomal score criterion
		System.out.println("Ribosomal Score Criterion:\t" + scoreCriterion_ribosomal);
		// print strength criterion
		System.out.println("Strength Criterion Freq:\t" + strengthCriterion_Freq);
		// print content criterion
		System.out.print("Content Criterion:\t\t" + contentCriterion);
		// label content criterion with interpretation (based on Carbone interpretation
		if(contentCriterion > .7)
			System.out.println("\t-GC3 BIAS");
		else if(contentCriterion < -.7)
			System.out.println("\t-AT3 BIAS");
		else
			System.out.println("\t-unremarkable");	
		// print amino acid bias criterion
		System.out.println("AA Bias Criterion:\t\t" + aaCriterion);
	}
	
	/**
	 * print for each CDS how many reference sets they are in
	 */
	public void		printRefset()
	{
		System.out.println();
		for(CDS cds: genome)
				System.out.print("[" + cds.getTotalRefsetContributes() + "] ");
	}
	
	// to files
	/**
	 * for each iteration, print all statistics to a BufferedWriter (iterationWriter)
	 */
	public void		printf_Stats() throws IOException
	{
		// codon frequencies
		iterationWriter.write("Codon Frequencies:");
		iterationWriter.newLine();
		Util.print8x8String(iterationWriter, codonFreq);
		
		// max synonymous codon frequencies
		iterationWriter.write("Max Synonymous Codon Frequencies:");
		iterationWriter.newLine();
		Util.print8x8String(iterationWriter, codonFreqMAX);
		
		// codon frequencies
		iterationWriter.write("Codon Frequencies:");
		iterationWriter.newLine();
		Util.print8x8String(iterationWriter, codonFreq);

		// preliminary w values
		iterationWriter.write("Preliminary W Values:");
		iterationWriter.newLine();
		Util.print8x8String(iterationWriter, prelimW);

		// g values
		iterationWriter.write("G Values:");
		iterationWriter.newLine();
		Util.print8x8String(iterationWriter, rcc);

		// iteration's w values
		iterationWriter.write("W Values:");
		iterationWriter.newLine();
		Util.print8x8String(iterationWriter, wValues);
		
		// new line separating iterations
		iterationWriter.newLine();
	}
	
	/**
	 * controller for printing data to files after final iteration
	 * 
	 * @throws Exception
	 */
	public void		printf_FINAL() throws Exception
	{
		// print entire genome, refset
		// abbreviated and in fasta form
		printf_GenomeRefset();
		
		// print strength
		calcCriterion();
		printf_Criterion();	
	}
	/**
	 * print ultimate reference set and genome to refset_info.txt and genome_info.txt
	 * and produce fasta file as genome_fasta.fas and refset_fasta.fas
	 * 
	 * @throws IOException from FileWriter
	 */
	public void 	printf_GenomeRefset() throws IOException
	{
		// open BufferedWriters
		BufferedWriter genomeWriter = new BufferedWriter(new FileWriter(folderOut + "/genome_info.txt"));
		BufferedWriter refsetWriter = new BufferedWriter(new FileWriter(folderOut + "/refset_info.txt"));
		//
		BufferedWriter genomeWriter_FASTA = new BufferedWriter(new FileWriter(folderOut + "/genome_fasta.fas"));
		BufferedWriter refsetWriter_FASTA = new BufferedWriter(new FileWriter(folderOut + "/refset_fasta.fas"));

		// sort
		Collections.sort(genome, Collections.reverseOrder());
		
		// for each CDS in the gnome
		String cdsInfoString = "";
		for(CDS eachCDS: genome)
		{
			
			// build string with name, score, other information
			cdsInfoString = "";
			Vector<String> cdsInfo = eachCDS.getInfo();
			for(String eachInfoBit: cdsInfo)
				if(!eachInfoBit.isEmpty())
					cdsInfoString += eachInfoBit + "\t";

			// write info to genome writer
			genomeWriter.write(eachCDS.getScore() + "\t");
			genomeWriter.write(cdsInfoString);
			genomeWriter.newLine();
			
			// add additional information to verbose genome writer
			genomeWriter_FASTA.write("> " + cdsInfoString);
			genomeWriter_FASTA.newLine();
			genomeWriter_FASTA.write(eachCDS.getSequence());
			genomeWriter_FASTA.newLine();
			
			// if the CDS is in the reference set, add information to refsetWriter as well
			if(eachCDS.getRefSetContributes(refsetRRIndex))
			{
				refsetWriter.write(cdsInfoString);
				refsetWriter.newLine();
				
				refsetWriter_FASTA.write("> " + cdsInfoString);
				refsetWriter_FASTA.newLine();
				refsetWriter_FASTA.write(eachCDS.getSequence());
				refsetWriter_FASTA.newLine();
			}
		}
		
		// close BufferedWriters
		genomeWriter.close();
		refsetWriter.close();
		//
		genomeWriter_FASTA.close();
		refsetWriter_FASTA.close();
	}
	/**
	 * print to a file an appended genome_info.txt file
	 * each CDS line will also contain the number of round robin reference sets
	 * the CDS is in
	 * 
	 * @param roundrobinRRN total # round robin runs
	 * @throws IOException	from the buffered writer
	 */
	public void 	printf_RRFinal(int roundrobinRRN) throws IOException
	{
		BufferedWriter genomeWriter = new BufferedWriter(new FileWriter(folderOut + "/genome_info_RR.txt"));

		String cdsInfoString = "";
		for(CDS eachCDS: genome)
		{
			// build string with name, score, other information
			cdsInfoString = "";
			Vector<String> cdsInfo = eachCDS.getInfo();
			for(String eachInfoBit: cdsInfo)
				if(!eachInfoBit.isEmpty())
					cdsInfoString += eachInfoBit + "\t";

			// write info to genome writer
//			genomeWriter.write(eachCDS.getScore() + "\t");
			genomeWriter.write(eachCDS.getRefsetContributes() + "\t");
			genomeWriter.write("RRTotal:" + eachCDS.getTotalRefsetContributes() + "/" + roundrobinRRN + "\t");
			genomeWriter.write(cdsInfoString);
			genomeWriter.newLine();
		}
		genomeWriter.close();
	}
	public void		printf_Refset() throws IOException
	{
		BufferedWriter refset_init_Writer = new BufferedWriter(new FileWriter(folderOut + "/refset_info_" + refsetRRIndex + "_.txt"));
		
		String cdsInfoString;
		for(CDS eachCDS: genome)
			if(eachCDS.getRefSetContributes(refsetRRIndex))
			{
				// build info string
				cdsInfoString = "";
				Vector<String> cdsInfo = eachCDS.getInfo();
				for(String eachInfoBit: cdsInfo)
					if(!eachInfoBit.isEmpty())
						cdsInfoString += eachInfoBit + "\t";
				
				// write
				refset_init_Writer.write(cdsInfoString);
				refset_init_Writer.newLine();
			}
				
		refset_init_Writer.close();
	}
	
	/**
	 * print criterions (ribosomal, strength, and content) to ribo_strength_criterion.txt
	 * @throws Exception from the BufferedWriter
	 */
	public void		printf_Criterion() throws Exception
	{
		// write to file
		BufferedWriter criterionWriter = new BufferedWriter(new FileWriter(folderOut + "/ribo_strength_criterion.txt"));
		criterionWriter.write("Ribosomal Score Criterion: " + scoreCriterion_ribosomal);
		criterionWriter.newLine();
		criterionWriter.write("Strength Criterion Freq: " + strengthCriterion_Freq);
		criterionWriter.newLine();
		criterionWriter.write("Content Criterion: " + contentCriterion);

		// close criterion writer
		criterionWriter.close();
	}

	/**
	 * Open the buffered writer to take notes on each iteration
	 * @throws IOException from BufferedWriter
	 */
	public void		openIterationWriter() 	throws IOException
	{
		iterationWriter = new BufferedWriter(new FileWriter(this.folderOut + "/iterationVerbose.txt"));		
	}
	/**
	 * Closes the buffered writer which takes notes on each iteration
	 * @throws IOException from BufferedWriter
	 */
	public void		closeIterationWriter() 	throws IOException
	{
		iterationWriter.close();
	}

	// ----------------- accessors and mutations // setters and getters
	/**
	 * starts new round robin (resets the round robin round #)
	 * resets w value storage
	 */
	public void		setrefsetRRIndex_OscillationMeasure(int refsetRRIndex) throws Exception
	{
		this.refsetRRIndex = refsetRRIndex;
		wValStore = new Vector<double[]>();
	}
	/**
	 * 
	 * @param folderOut path to folder
	 */
	public void		setOutputFolder(String folderOut)
	{
		this.folderOut = folderOut;
	}
	/**
	 * @return size of genome
	 */
	public int 		getGenomeSize()
	{
		return genomeSize;
	}
	/**
	 * 
	 * @return frequency of codons in reference set (shallow)
	 */
	public double[] getCodonFreq()
	{
		return codonFreq;
	}
	/**
	 * 
	 * @return w values (shallow)
	 */
	public double[] getW()
	{
		return wValues;
	}
	/**
	 * compares the current w values with the last set of w values
	 * if the w values match, then the reference set must have stabilized
	 * 
	 * @return whether reference set has stabilized
	 */
	public boolean	isStable()
	{
		return Arrays.equals(wValStore.lastElement(), wValues);
	}
	/**
	 * compares the current w values with all other sets of w values
	 * if the w values match, then the reference set must be oscillating
	 * 
	 * return -1 if not oscillating
	 * otherwise return the number of iterations back the 
	 * 
	 * @return size of oscillations of reference set (-1 if not oscillating)
	 */
	public int		isOscillating()
	{
		for(double[] pastWval: wValStore)
			if(Arrays.equals(pastWval, wValues))
				return wValStore.indexOf(pastWval);
		return -1;
	}
	/**
	 * set all contributing
	 */
	public void		setAllContributing(int roundnum)
	{
		for(CDS cds: genome)
			cds.setRefSetContributes(roundnum, true);
	}
}
