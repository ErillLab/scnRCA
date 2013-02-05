/*
 * FILL PARAMS MAY BE ODD? CHECK iteration branching
 * 
 * 
 * made out of OrCAI 4
 * 
 * Round robin starting point change
 * 
 * take out stop-start-met-trp
 * laplace pseudocounts
 * 
 * 
 * 
 */


//-i c -g true -d 2.0 -p 0.5 -f C:\Documents and Settings\Mindy\My Documents\workspace\OrCAI_3\genome_source\Mycoplasma_gallisepticum_str_R_low_.gb -m -1

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Scanner;
import java.util.Vector;

import org.biojava.bio.Annotation;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojavax.bio.seq.RichFeature;

	
@SuppressWarnings("deprecation")
public class OrCAI_GUI_main 
{
	// user set
	private static char 	index;
	private static boolean	rcc;
	private static double 	divisor;
	private static double 	finalRefsetPercent;
	private static String 	genomePath;
	private static String 	fileOrg;
	private static int		maxIterations;
	//
	private static boolean	roundrobin;
	private static int		roundrobinRRN;
	private static int		roundrobinISS;
	private static boolean	roundrobinRandom;
	//
	private static boolean	selectRefset;
	private static String 	refsetPath;
	
	// derived
	private static String	outputFolderName;
	private static Genome	genome;
	private static int		genomeSize;
	
	private static double[] aaCrits;

	/**
	 * main method - program driver
	 * 
	 * @param args command line
	 */
	public static void main(String[] args) throws InterruptedException, IllegalAlphabetException, IllegalSymbolException
	{
		try
		{
			// initialize variables
			index 				= 'c';
			rcc 				= true;
			divisor 			= 2d;
			finalRefsetPercent 	= 1;
			genomePath 			= "";
			fileOrg 			= "";
			maxIterations 		= -1;
			roundrobin 			= false;
			roundrobinRRN		= 10;
			roundrobinISS		= 1;
			roundrobinRandom	= false;
			
			
			// gui to get params (if not already given in cmd line)
			OrCmd_GUI cmdGUI = new OrCmd_GUI();
			
			do{
				//------------------------------- SET UP
				String[] cmds;
				if(args.length == 0)
				{
					System.out.println("Awaiting Orders from Graphical User Interface");
					System.out.println("Please press \"Compute\" or \"Exit\"");
					cmdGUI.setVisible();
					while(!cmdGUI.readyToRead())
						Thread.sleep(100);
					System.out.println("---\n" + cmdGUI.cmdline() + "\n");
					cmds = cmdGUI.cmdline().split(" ");
				}else
					cmds = args;
				
				// fill run parameters
				fillParams(cmds);

				genome = new Genome(parsefile(genomePath), index, rcc, fileOrg, roundrobinRRN);
//				genomeSize = genome.getGenomeSize();
				aaCrits = new double[roundrobinRRN];
				
				//-------------------------------- ROUND ROBIN
				for(int rrStep = 0; rrStep < roundrobinRRN; rrStep++)
				{
//					genome.printRefset();

					genome.setrefsetRRIndex_OscillationMeasure(rrStep);
					outputFolderName = generateFolderName(rrStep);
					generateFolder(outputFolderName);
					genome.setOutputFolder(outputFolderName);
					genome.openIterationWriter();
					
					runRound(rrStep);

					genome.closeIterationWriter();
				}
				// print round robin results
				outputFolderName = generateFolderName(-99999);
				outputFolderName = outputFolderName.replace("-99999", "FINAL");
				System.out.println(outputFolderName);
				generateFolder(outputFolderName);
				genome.setOutputFolder(outputFolderName);
				genome.printf_RRFinal(roundrobinRRN);
				
				// print amino acid criterion
				/*
				System.out.println("AA Crits");
				for(int i = 0; i < roundrobinRRN; i++)
				{
					System.out.println(aaCrits[i]);
				}
				*/
				
				//------------------------------- FINAL REFSET
//				genome.printRefset();
				
				System.out.println();				
			// keep going if the gui is invoked
			}while(args.length == 0);
			
			System.out.println("\n" + "no errors!");
			
		}catch(Exception e)
		{
			System.out.println("EXCEPTION!: ---> " + e.getMessage());
			e.printStackTrace();
		}
		System.out.println("Done");
	}// end main
	

	
	// BIG BOSS
		private static void runRound(int rrStep) throws Exception
		{
			/**
			genome.setrefsetRRIndex_OscillationMeasure(rrStep);
			genome.printf_Refset();
			**/
			genome.printAACrit();
			aaCrits[rrStep] = genome.aaCrit();
			
			// ------------------- RUN ITERATIONS
			//-------------------- ITERATIONS
			System.out.print("\n------------------IN THE ITERATIONS");
			if(roundrobin)
				System.out.print(" :: ROUND " + rrStep);
			System.out.println();
			
			// calculate initial reference set size
			int refsetSize;
			double refsetPercent;
			
			refsetSize = genomeSize;
			refsetPercent = 100d;
			if(roundrobin || selectRefset)
			{
				refsetPercent = finalRefsetPercent;
				refsetSize = (int) (genomeSize * (refsetPercent / 100));
			}

			int 	iterationCount 		= 1;
			boolean continueIterating 	= false;
			String 	whyStop 			= "";
			
			boolean firstRun = true;
			do
			{
				//---------------------- BUILD REF SET and score the CDSs
				//-- select refset calculate and port in w values
				if(selectRefset && firstRun)
				{
					// calculate w (from selected reference set)
					Genome newRefset = new Genome(parsefile(refsetPath), index, rcc, "", 1);
					newRefset.setAllContributing(0);
					newRefset.recalcRefset_simple();
					// set w
					genome.setW(newRefset.getW());
					firstRun = false;				
				}
				//-- cut the reference set(avoid recutting the reference set in the first run in the round robin)
				//-- to make new w values
				else
				{
					if((roundrobin && ! firstRun) || !roundrobin)
						genome.recutRefSet(refsetSize);
					else
						firstRun = false;
					
					genome.recalcRefSet(iterationCount);
				}
				
				//----- score the CDSs
				genome.scoreCDS();
				genome.printW();
				
				//---------------------- CONTINUE ITERATING? determination
				continueIterating = false;
				whyStop = "";
				// if max iterations is specified and max iterations has been reached, do not continue
				if(maxIterations != -1 && iterationCount >= maxIterations)
					whyStop = "Iteration Maxed!";
				// CHECK STABILITY
				else if(iterationCount > 2 && genome.isStable())
					whyStop = "Stabilized";
				else if(iterationCount > 2 && genome.isOscillating() != -1)
					whyStop = "Oscillating with " + (iterationCount - genome.isOscillating());
				else if(refsetSize == 1)
					whyStop = "About to Run Out of CDSs";				
				//---------------------- NEXT ITERATION PREP (as long as continuing iteration is good)
				else
				{
					continueIterating = true;
					// step number of iterations
					iterationCount ++;
					// calculate # of CDSs in next reference set if using the carbone method
					// recalc percent
					if(! roundrobin && ! selectRefset){
						refsetPercent /= divisor;
						if(refsetPercent < finalRefsetPercent)
							refsetPercent = finalRefsetPercent;
						
						refsetSize = (int) (genomeSize * ((double)refsetPercent / (double)100));
					}
				}
				
				if(!whyStop.isEmpty())
					continueIterating = false;
				
			}while(continueIterating);
			
			System.out.println("\n----------------- " + whyStop + " at " + (iterationCount - 1));			
			printFinalStatsToFile(whyStop, iterationCount);
				// print to file
			genome.printf_FINAL();
				// to console - must be after printing to file
			genome.printCriterion();
		}
	
	
	// ----------------- print stats -------------------------------------------------------------------
		private static void 	printFinalStatsToFile(String whyStop, int iterationCount) throws IOException
		{
			BufferedWriter runtimeNotesWriter = new BufferedWriter(new FileWriter(new File(outputFolderName + "/" + "runtime_log.txt")));

			String userStats[] = userStats().split("\n");
			for(String eachStat: userStats)
			{
				runtimeNotesWriter.write(eachStat);
				runtimeNotesWriter.newLine();
			}

			runtimeNotesWriter.write("Why Stop: " + whyStop);
			
			runtimeNotesWriter.close();
		}
		private static void 	eachIterPrint()
		{
//			genome.printCodonCounts();
//			refset.printPositionalFreq();
//			genome.printCodonFreq();
//			refset.printCodonFreqMAX();
//			genome.printPrelimW();
//			genome.printG();
//			genome.printW();
//			genome.printCDSscores();
		}
		private static String 	userStats()
		{
			String userStats = "";
			userStats += "File Org:\t\t\t" + fileOrg + "\n";
			userStats += "Index:\t\t\t\t" + index + "\n";
			userStats += "g factor use:\t\t\t" + rcc + "\n";
			userStats += "Division Factor:\t\t" + divisor + "\n";
			userStats += "Final Refset Percent-Size:\t" + finalRefsetPercent + "\n";
			userStats += "File Name:\t\t\t" + genomePath + "\n";
			userStats += "Max Iterations:\t\t\t" + maxIterations + "\n";
			userStats += "Round robin:\t\t\t" + roundrobin + "\n";
			return userStats;
		}

	// ----------------- general tools -------------------------------------------------------------------
		private static void		generateFolder(String path)
		{
			File outputDir = new File(path);
			if(!outputDir.exists())
				outputDir.mkdirs();
		}
		private static String 	generateFolderName(int rrStep) throws Exception
	{
		String foldername = fileOrg + "-";
		
		foldername += finalRefsetPercent + "_";
		
		if(rcc)
			foldername += "rcc_";
		
		switch(index)
		{
		case 'c': foldername += "CAI"; break;
		case 'r': foldername += "nRCA"; break;
		default:
			throw new Exception("bad whichW in generate folder name");
		}
		
		if(roundrobin)
			foldername += "/round_" + rrStep;


		
		return foldername;
	}
		private static char 	detectFileType(String filePath) throws Exception
		{
			char fileType;
			
			filePath = filePath.trim();			
			String[] fileNameSplit = filePath.split("\\.");
			String extension = fileNameSplit[fileNameSplit.length - 1];
			if(extension == "fas" || extension == "fasta" || extension == "fma")
				fileType = 'f';
			else if(extension == "gb" || extension == "gbk" || extension == "gbwithparts")
				fileType = 'g';
			else
			{
				Scanner readFile = new Scanner(new FileReader(filePath));
				char firstChar = readFile.nextLine().charAt(0);
				if(firstChar == '>')
					fileType = 'f';
				else if(firstChar == 'L')
					fileType = 'g';
				else
					throw new Exception("What Kind of File is this?! Only FASTA/Genbank allowed.");
				readFile.close();
			}
			
			return fileType;
		}

		
	// ----------------- PARSING for CDS ------------------------------------------------------------------
		/**
		 * parsing controls
		 * 
		 * @param 	fileformat 	format that the file is in (f-fasta, g-genbank)
		 * @param	filename	name of file to pull CDS information from
		 * @return				vector of CDS objects
		 */
		private static Vector<CDS>	parsefile(String filePath) throws Exception
		{

			Vector<CDS> cdsVector = new Vector<CDS>();
			
			// detect file type
			filePath = filePath.trim();
			char fileformat = detectFileType(filePath);
			
			// parse files in their respective formats
			switch(fileformat)
			{
			case 'f': cdsVector = parsefileFASTA(filePath); 	break;
			case 'g': cdsVector = parsefileGENBANK_2(filePath);	break;
			default: throw new Exception("Invalid File Format to Parse: " + fileformat);
			}
			
			System.out.println("PARSED in parsefile: " + cdsVector.size() + "# CDSs");
			if(cdsVector.size() == 0)
				throw new Exception("PARSED 0 CDSs in parsefile from " + filePath);
			genomeSize = cdsVector.size();
			
			// fill the round robin ness
			if(roundrobinRandom)
				fillRoundRobin_Random(cdsVector);
			else
				fillRoundRobin(cdsVector);
			
			
			return cdsVector;
		}
		
		/**
		 * parses genome files (FASTA) --> genome vector of all CDS
		 * 
		 * FIX FASTA DIRECTIONAL
		 * 
		 * 
		 * @param	filename	name of file
		 * @param	fileFormat	format of file ('f' - FASTA, 'g' - GenBank)
		 * 
		 * @throws	IOException
		 * @throws	NoSuchElementException
		 * @throws	BioException
		 */
		private static Vector<CDS>	parsefileFASTA(String filePath) throws Exception
		{
			Vector<CDS> cdsFASTA = new Vector<CDS>(30000);
			
			// open file
			System.out.println("Opening File " + filePath);
			

			// initialize cds sequence values
			String cdsSequence = "";
			String cdsName 		= "";
			String cdsLocation	= "";
			String cdsDirection	= "";
			String cdsFunction	= "";
			String cdsProtein	= "";
			String cdsProduct	= "";
			String cdsNote		= "";
			
//			boolean[] roundrobinRefset = new boolean[roundrobinRRN];
			boolean firstseq = true;
			
			// run through each line in the FASTA file
			String line = "";
			Scanner fastaFile = new Scanner(new FileReader(filePath));
			while(fastaFile.hasNext())
			{
				line = fastaFile.nextLine();
				
				// header
				if(line.charAt(0) == '>')
				{
					if(!firstseq)
					{
						// make new CDS
						CDS cds = new CDS(cdsSequence.toLowerCase());
						// add information
						cds.addInfo("", cdsName, "", cdsLocation, cdsDirection, cdsFunction, cdsProtein, cdsProduct, cdsNote);
						// add the CDS object to genome vector of CDSs
						cdsFASTA.add(cds);
					}
					firstseq = false;
					
					cdsNote = line;
					cdsSequence = "";
				}
				// comment
				else if(line.charAt(0) == ';')
				{}
				// sequence
				else
				{
					cdsSequence += line.trim();
				}
			}//end loop
			
			// add last sequence
			// make new CDS
			CDS cds = new CDS(cdsSequence.toLowerCase());
			// add information
			cds.addInfo("", cdsName, "", cdsLocation, cdsDirection, cdsFunction, cdsProtein, cdsProduct, cdsNote);
			// add the CDS object to genome vector of CDSs
			cdsFASTA.add(cds);
			
			
			return cdsFASTA;
		}	
		/**
		 * 
		 * @param filename
		 * @return vector of CDSs comprising the GenBank file
		 * @throws Exception from FileReader and various BioJava functions
		 */
		private static Vector<CDS> 	parsefileGENBANK_2(String filePath) throws Exception
		{
			Vector<CDS> cdsGENBANK = new Vector<CDS>(30000);

			// open file
			System.out.println("Opening File \"" + filePath + "\"");

			int totalGenes = 0;
			Map<String, String> geneNoteMap = new HashMap<String, String>();


			
			BufferedReader 		geneFile1 = new BufferedReader(new FileReader(filePath));
			SequenceIterator	iterG1 = SeqIOTools.readGenbank(geneFile1);
			
			Feature ft1;
			Annotation ftAnnotation1;
			Sequence seq1;
			
			String geneName = "";
			String geneNote = "";
			while(iterG1.hasNext())
			{
				seq1 = iterG1.nextSequence();
				for(Iterator<Feature> i = seq1.features(); i.hasNext();)
				{
					ft1 = i.next();
					if(ft1.getType().equals("gene"))
					{
						
						ftAnnotation1 = ft1.getAnnotation();
						if(ftAnnotation1.containsProperty("note") && ftAnnotation1.containsProperty("gene"))
						{
			              geneName = ftAnnotation1.getProperty("gene").toString();
			              geneNote = ftAnnotation1.getProperty("note").toString();
			              
//			              System.out.println(geneName+"\t: "+geneNote);
			              geneNoteMap.put(geneName, geneNote);
						}
					}
					
					else if(ft1.getType().equals("CDS"))
						totalGenes ++;
				}
			}
			geneFile1.close();
			
			
			BufferedReader 		geneFile = new BufferedReader(new FileReader(filePath));
			//GenBank
			SequenceIterator iterG = SeqIOTools.readGenbank(geneFile);
			
			
//			boolean[] roundrobinRefset = new boolean[roundrobinNumRounds];
			
			Feature ft;
			String cdsLocusTag;
			String cdsName;
			String cdsGeneSynonym;
			int cdsLocationInt;
			String cdsLocation;
			String cdsDirection;
			String cdsSequence;
			String cdsFunction;
			String cdsProtein;
			String cdsProduct;
			String cdsNote;
			CDS cds;
			
			while(iterG.hasNext())
			{

				// each seq
				Sequence seq = iterG.nextSequence();

				// iterate through features
				for(Iterator<Feature> seqFeat = seq.features(); seqFeat.hasNext();)
				{						
					ft = seqFeat.next();

					// for all protein coding sequence
					if(ft.getType().equals("CDS"))
					{
						Annotation ftAnnotation = ft.getAnnotation();
						
						// acquire cds's name, location, direction, sequence, function, protein product
						cdsLocusTag		= "";
						cdsName 		= "";
						cdsGeneSynonym	= "";
						cdsLocationInt 	= ft.getLocation().getMin();
						cdsLocation 	= ft.getLocation().toString();
						try{
							cdsDirection	= Character.toString(((RichFeature)ft).getStrand().getToken());
						}catch(Exception e)
						{
							cdsDirection = "";
						}
						cdsSequence 	= ft.getSymbols().seqString();
						cdsFunction		= "";
						cdsProtein		= "";
						cdsProduct		= "";
						cdsNote			= "";

						if(ftAnnotation.containsProperty("locus_tag"))
							cdsLocusTag		= ft.getAnnotation().getProperty("locus_tag").toString();
						if(ftAnnotation.containsProperty("gene"))
							cdsName 		= ft.getAnnotation().getProperty("gene").toString();
						if(ftAnnotation.containsProperty("gene_synonym"))
							cdsGeneSynonym		= ft.getAnnotation().getProperty("gene_synonym").toString();
						if(ftAnnotation.containsProperty("function"))
							cdsFunction = ft.getAnnotation().getProperty("function").toString();
						if(ftAnnotation.containsProperty("protein_id"))
							cdsProtein		= ft.getAnnotation().getProperty("protein_id").toString();
						if(ftAnnotation.containsProperty("product"))
							cdsProduct		= ft.getAnnotation().getProperty("product").toString();
						if(geneNoteMap.containsKey(cdsName))
							cdsNote = geneNoteMap.get(cdsName);


						
						 // print attributes
						/*
						 		System.out.println("Locus Tag: "	+ cds:LocusTag);
								System.out.println("Name: " 		+ cdsName);
								System.out.println("Location: " 	+ cdsLocation);
								System.out.println("Direction: "	+ cdsDirection);
								System.out.println("Function: " 	+ cdsFunction);
								System.out.println("Protein: " 		+ cdsProtein);
								System.out.println("Sequence: " 	+ cdsSequence);
								System.out.println("Product: "		+ cdsProduct);
						*/

						// create round robin contribution
						/*
						if(roundrobin)
						{
							for(int i = 0; i < roundrobinNumRounds; i++)
							{
								if(i == roundrobinStep)
									roundrobinRefset[i] = true;
								else
									roundrobinRefset[i] = false;
							}
							roundrobinStep ++;
							roundrobinStep %= roundrobinNumRounds;
							
						}
						else
							roundrobinRefset[0] = true;
						*/
						
						// create CDS object (with deep copy of the refset indicator)
						cds = new CDS(cdsSequence);
//						cds = new CDS(cdsSequence, Arrays.copyOf(roundrobinRefset, roundrobinNumRounds));
						// add extra information to CDS object
						cds.addInfo(cdsLocusTag, cdsName, cdsGeneSynonym, cdsLocation, cdsDirection, cdsFunction, cdsProtein, cdsProduct, cdsNote);
						cds.setLocation(cdsLocationInt);
						// add CDS object to genome vector of CDSs
						cdsGENBANK.add(cds);
						

					}// end if
				}// end for seq features

			}// end sequence iterator

			geneFile.close();
			
			Collections.sort(cdsGENBANK, new CDSComparatorLocation());

			return cdsGENBANK;
		}
		
		private static void 		fillRoundRobin(Vector<CDS> genomeBlank) throws Exception
		{
			int[] rrTable = fillRRHelper();
			
			int tableSize = rrTable.length;		
			int initSet;
			
			boolean[] rrOk;
			int i = 0;
			for(CDS cds: genomeBlank)
			{
				initSet = rrTable[i % tableSize];
				rrOk = new boolean[roundrobinRRN];
				for(int j = 0; j < roundrobinRRN; j ++)
				{
					if(j == initSet)
						rrOk[j] = true;
					else
						rrOk[j] = false;
				}
				cds.setRoundRobin(rrOk);
				i++;
			}
			
			System.out.println("Ello poppet NON RANDOM");
			
			// printing purposes
			/**
			for(int i = 0; i < genomeSize; i++)
			{
				for(int j = 0; j < rrSetNum; j++)
					System.out.print((genomeRR[i][j]?1:0) + " ");
				System.out.println();
			}
			**/
			
		}
		
		private static void			fillRoundRobin_Random(Vector<CDS> genomeBlank)
		{			
			System.out.println("Ello poppet RANDOM");

			
			double randomD;
			
			boolean[] rrOk;
			for(CDS cds: genomeBlank)
			{
				rrOk = new boolean[roundrobinRRN];
				for(int i = 0; i < roundrobinRRN; i++)
					rrOk[i] = false;
				cds.setRoundRobin(rrOk);
			}
			
			
			int whichCDS;
			CDS cds;
			// for all the round robin rounds rounds
			for(int i = 0; i < roundrobinRRN; i++)
			{
				// for every CDS in the starting reference set, set it!
				for(int j = 0; j < (roundrobinISS * genomeSize / 100); j++)
				{
					// generate random #
					randomD = Math.random();
					whichCDS = (int) Math.floor(randomD * genomeSize);
//					System.out.println("i j: " + i + " " + j + " " + whichCDS);
					
					// check if it's already set
					cds = genomeBlank.get(whichCDS);

//					System.out.println("CDS: " + whichCDS + " " + cds.getRefSetContributes(i));
					if(!cds.getRefSetContributes(i))
					{
						// set it
						cds.setRefSetContributes(i, true);
//						System.out.println("VICTORY IS MINE " + whichCDS);
					}
					// else don't up the counter
					else
					{
						j --;
//						System.out.println("Alreadydone");
					}
				}
			}
			
			/*
			boolean setRR;
			for(CDS cds2: genomeBlank)
			{
				System.out.println(cds2.getRefsetContributes());
			}
			*/
			
			// exit. just for now.
//			System.out.println("Exit for funzies");
//			System.exit(0);
		}
		
		private static int[] 		fillRRHelper() throws Exception
		{
			// print params
			System.out.println("Round Robin Number of Sets: " + roundrobinRRN);
			System.out.println("Reference Set Size (% of genome): " + roundrobinISS);
			
			// calculate chunk and step size
			int chunkSize	= (int) (100 / roundrobinISS);
			int stepSize	= chunkSize / roundrobinRRN;

			// print chunk and step size
			System.out.println("Chunk Size: " + chunkSize);
			System.out.println("Step Size: " + stepSize);
			
			if(stepSize < 1) throw new Exception ("chunkRR Zero Attack");
			
			int[] rr = new int[chunkSize];
			
			for(int i = 0; i < chunkSize; i ++)
				rr[i] = -1;
			for(int i = 0; i < roundrobinRRN; i++)
				rr[i * stepSize] = i;
			
			/*
			System.out.println("Table:");
			for(int i = 0; i < chunkSize; i ++)
			{
				if(i % stepSize == 0) System.out.println();
				System.out.print(rr[i] + " ");
			}
			System.out.println();
			*/
			
			return rr;
		}


	// ------- FILL PARAMS
		private static void 		fillParams(String[] cmds) throws Exception
		{
			// loop through commands fill params
			char iParam = '0';
			String iCMD = "";
			boolean onF = false;
			boolean onS = false;
			genomePath = "";
			refsetPath = "";
			for(int i = 0; i < cmds.length; i++)
			{
				if(cmds[i].charAt(0) == '-' && !cmds[i].equals("-1"))
				{
					iParam	= cmds[i].charAt(1);
					iCMD 	= cmds[i+1];

					onF = false;
					onS = false;
					
					switch(iParam)
					{
					case 'i': 	index 				= iCMD.charAt(0);			break;
					case 'g': 	rcc 				= iCMD.equals("true");		break;
					case 'd': 	divisor 			= Double.parseDouble(iCMD); break;
					case 'p': 	finalRefsetPercent 	= Double.parseDouble(iCMD);	break;
					case 'f': 	onF 				= true;						break;
					case 'm':	maxIterations		= Integer.parseInt(iCMD);	break;
					case 'r':	roundrobin 			= true;
								roundrobinRRN		= Integer.parseInt(iCMD);
								roundrobinISS		= Integer.parseInt(cmds[i+2]);
								roundrobinRandom	= (boolean) Boolean.parseBoolean(cmds[i+3]);
																				break;
					case 's':	selectRefset		= true;
								onS					= true;
																				break;
					}// end switch
				}
				else if(onF)
					genomePath += " " + cmds[i];
				else if(onS)
					refsetPath += " " + cmds[i];

			}// end for
			
			if(!roundrobin)
			{
				roundrobinRRN = 1;
				roundrobinISS = 1;
			}
			
			genomePath = genomePath.trim();
			refsetPath = refsetPath.trim();

							
			// detect file organism
			fileOrg = ((new File(genomePath)).getName().split("\\."))[0];
		}

}
