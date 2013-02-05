/**
 * rr corrected, stage 1
 * 
 * class ExtensionFileFilter lacks documentation
 * otherwise done
 * 
 * GUI for scRCA
 * interface for user to determine parameters to run program
 * by generating a command line to pass to the main program to parse
 * 
 * command line can only be accessed once per "Compute" button press
 */

import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.*;
import javax.swing.filechooser.FileFilter;


public class OrCmd_GUI{
	// main frame
	private final JFrame frame = new JFrame();
	// buttons
	private JButton computeButton	= new JButton("Compute");
	private JButton exitButton		= new JButton("Exit");
	
	// output
	private String	cmdLine			= "";
	private boolean ready			= false;
	// parameters
	private char	initialRefset		= 'w';
	private int		roundrobinRNN		= 10;
	private int		roundrobinISS		= 1;
	private boolean roundrobinRandom	= false;
	private char	index				= 'c';
	private boolean rcc					= true;
	private double	divisor				= 2d;
	private double	finalRefsetPercent 	= 1;
	private String	filePath			= "";
	private String	dirPath				= "";
	private String 	refsetPath			= "";
	private int		maxIterCount		= -1;
	
	/**
	 * constructor
	 * builds GUI
	 */
	public OrCmd_GUI()
	{
		// frame settings
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setLayout(new GridLayout(4,0));
		frame.setPreferredSize(new Dimension(500,600));
		frame.setTitle("scnRCA");
		
		// panels in frame
		JPanel paneInOut	= new JPanel(new GridLayout(2,0));
		JPanel paneIndex	= new JPanel(new GridLayout(0,2));
		JPanel paneRefsetStart = new JPanel(new GridLayout(3,0));
		JPanel paneButtons	= new JPanel(new GridLayout(2,0));
		
		//------------------------ SELECT FILE, FOLDER
		// file
		JPanel paneInFile	= new JPanel(new GridLayout(0,2));	
		paneInFile.setBorder(BorderFactory.createTitledBorder("Select Genome File"));
		final JTextField filepathField = new JTextField();
		JButton fileButton = new JButton("Browse");
		fileButton.addActionListener(new ActionListener()
		{
			public void actionPerformed(ActionEvent e)
			{
				// filter for FASTA and GenBank files only
				JFileChooser fileChooser = new JFileChooser();
				FileFilter filefilter = new ExtensionFileFilter("FASTA AND GENBANK ONLY", new String[]{"gb","gbk","fas","fasta","fma"});
				fileChooser.setFileFilter(filefilter);
				fileChooser.setCurrentDirectory(new java.io.File("."));

				// set text field to chosen file
				int returnVal = fileChooser.showOpenDialog(null);
		        if (returnVal == JFileChooser.APPROVE_OPTION) {
		            File file = fileChooser.getSelectedFile();
		            filepathField.setText(file.getPath());
		        }
			}
		});
		// put file selecting elements into the file panel
		paneInFile.add(filepathField);
		paneInFile.add(fileButton);
		
		//-- folder
		JPanel paneOutDir	= new JPanel(new GridLayout(0,2));
		paneOutDir.setBorder(BorderFactory.createTitledBorder("Select Output Folder"));
		// text field
		final JTextField dirpathField = new JTextField();
		dirpathField.setText((new java.io.File(".")).getAbsolutePath());
		// button
		JButton dirButton = new JButton("Browse");
		dirButton.addActionListener(new ActionListener()
		{
			public void actionPerformed(ActionEvent e)
			{
				JFileChooser dirChooser = new JFileChooser();
				dirChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				
				dirChooser.setCurrentDirectory(new java.io.File("."));
		        int returnVal = dirChooser.showOpenDialog(null);

		        if (returnVal == JFileChooser.APPROVE_OPTION) {
		            File file = dirChooser.getSelectedFile();
		            dirpathField.setText(file.getPath());
		        }
			}
		});
		//-- put file selecting elements into the file panel
		paneOutDir.add(dirpathField);
		paneOutDir.add(dirButton);
		
		// in and out together
		paneInOut.add(paneInFile);
		paneInOut.add(paneOutDir);
		
		//------------------------ Index
		JPanel paneIndexLeft = new JPanel(new GridLayout(2,0));
		paneIndexLeft.setBorder(BorderFactory.createTitledBorder("Basic"));
		final JRadioButton rcaRadio = new JRadioButton("nRCA");
		final JRadioButton caiRadio = new JRadioButton("CAI");
		rcaRadio.setSelected(true);
		
		ButtonGroup group = new ButtonGroup();
		group.add(rcaRadio);
		group.add(caiRadio);
		
		final JCheckBox rccCheck = new JCheckBox("Rare Codon Correction");
		rccCheck.setSelected(true);
		
		JPanel paneIndexRadio = new JPanel(new GridLayout(0,2));
		paneIndexRadio.setBorder(BorderFactory.createTitledBorder("Index"));
		paneIndexRadio.add(rcaRadio);	
		paneIndexRadio.add(caiRadio);
		
		paneIndexLeft.add(paneIndexRadio);
		paneIndexLeft.add(rccCheck);
		
		///
		JPanel paneIndexTweaks	= new JPanel(new GridLayout(3,0));
		paneIndexTweaks.setBorder(BorderFactory.createTitledBorder("Tweaks"));
		
		// tweak RIGHT PANE
//		JPanel paneTweakRight = new JPanel(new GridLayout(0,3));
		// divisor
		final JTextField divisorField = new JTextField();
		divisorField.setText("2");
		divisorField.setBorder(BorderFactory.createTitledBorder("Divisor"));
		// final set size
		final JTextField refsetPercent = new JTextField();
		refsetPercent.setText("1");
		refsetPercent.setBorder(BorderFactory.createTitledBorder("Final Set Size (%)"));
		// max iterations
		final JTextField maxIterations = new JTextField();
		maxIterations.setText("-1");
		maxIterations.setBorder(BorderFactory.createTitledBorder("Max Iterations"));
		//
		paneIndexTweaks.add(divisorField);
		paneIndexTweaks.add(refsetPercent);
		paneIndexTweaks.add(maxIterations);
		
		
		paneIndex.add(paneIndexLeft);
		paneIndex.add(paneIndexTweaks);
		
		//------------------------ STARTING REFERENCE SET
		// tweak LEFT PANE -- starting reference set pattern
		paneRefsetStart.setBorder(BorderFactory.createTitledBorder("Starting Reference Set(s)"));
		//
		final JRadioButton 	wholegenomeRadio 	= new JRadioButton("Whole Genome");
		final JRadioButton 	roundrobinRadio 	= new JRadioButton("Round Robin");
		final JRadioButton	selectrefsetRadio	= new JRadioButton("Select Reference Set");
		final JTextField 	roundrobinRRNField 	= new JTextField("10");
		roundrobinRRNField.setBorder(BorderFactory.createTitledBorder("RNN"));
		final JTextField	roundrobinISSField 	= new JTextField("1");
		roundrobinISSField.setBorder(BorderFactory.createTitledBorder("ISS"));
		final JCheckBox		roundrobinRANDOMcheck = new JCheckBox("Random");
		wholegenomeRadio.setSelected(true);
		final JTextField	selectRefsetField	= new JTextField();
		// button
		JButton selectrefsetButton = new JButton("Browse");
		selectrefsetButton.addActionListener(new ActionListener()
		{
			public void actionPerformed(ActionEvent e)
			{
				// filter for FASTA and GenBank files only
				JFileChooser fileChooser = new JFileChooser();
				FileFilter filefilter = new ExtensionFileFilter("FASTA AND GENBANK ONLY", new String[]{"gb","gbk","fas","fasta","fma"});
				fileChooser.setFileFilter(filefilter);
				fileChooser.setCurrentDirectory(new java.io.File("."));

				// set text field to chosen file
				int returnVal = fileChooser.showOpenDialog(null);
		        if (returnVal == JFileChooser.APPROVE_OPTION) {
		            File file = fileChooser.getSelectedFile();
		            selectRefsetField.setText(file.getPath());
		        }
			}
		});
		//
		JPanel paneRadio = new JPanel(new GridLayout(2,0));
		paneRadio.add(wholegenomeRadio);
		paneRadio.add(roundrobinRadio);
		JPanel paneFieldRR = new JPanel(new GridLayout(0,3));
		paneFieldRR.add(roundrobinRRNField);
		paneFieldRR.add(roundrobinISSField);
		paneFieldRR.add(roundrobinRANDOMcheck);
		JPanel paneSelectRefset = new JPanel(new GridLayout(2,0));
		paneSelectRefset.add(selectrefsetRadio);
		JPanel paneSelectRefsetField = new JPanel(new GridLayout(0,2));
		paneSelectRefsetField.add(selectRefsetField);
		paneSelectRefsetField.add(selectrefsetButton);
		paneSelectRefset.add(paneSelectRefsetField);
		//
		ButtonGroup refsetGroup = new ButtonGroup();
		refsetGroup.add(wholegenomeRadio);
		refsetGroup.add(roundrobinRadio);
		refsetGroup.add(selectrefsetRadio);
		//
		paneRefsetStart.add(paneRadio);
		paneRefsetStart.add(paneFieldRR);
		paneRefsetStart.add(paneSelectRefset);
		
		//------------------------ TWEAKS

		
//		paneTweaks.add(paneRefsetStart);
//		paneTweaks.add(paneTweakRight);
		
		//------------------------ ACTION
		paneButtons.setBorder(BorderFactory.createTitledBorder("Action!"));
		computeButton.addActionListener(new ActionListener()
		{
			public void actionPerformed(ActionEvent e)
			{
				// filepath
				filePath 	= filepathField.getText();
				dirPath		= dirpathField.getText();
				refsetPath	= selectRefsetField.getText();
				
				// refset initial
				if(wholegenomeRadio.isSelected())
					initialRefset = 'w';
				else if(roundrobinRadio.isSelected())
					initialRefset = 'r';
				else if(selectrefsetRadio.isSelected())
					initialRefset = 's';
				
				roundrobinRNN = Integer.parseInt(roundrobinRRNField.getText());
				roundrobinISS = Integer.parseInt(roundrobinISSField.getText());
				roundrobinRandom = roundrobinRANDOMcheck.isSelected();
				
				// index
				if(caiRadio.isSelected())
					index = 'c';
				else if(rcaRadio.isSelected())
					index = 'r';
				
				// rcc
				rcc = rccCheck.isSelected();
				
				// divison factor
				divisor = Double.parseDouble(divisorField.getText());

				// reference set percent
				finalRefsetPercent = Double.parseDouble(refsetPercent.getText());
				
				// max iterations
				maxIterCount = Integer.parseInt(maxIterations.getText());
				
				if(filePath.isEmpty())
					JOptionPane.showMessageDialog(null, "File Input Path is Empty");
				else if(dirPath.isEmpty())
					JOptionPane.showMessageDialog(null, "Output Directory Path is Empty");
				else if(divisor < 1)
					JOptionPane.showMessageDialog(null, "Divisor cannot be less than 1");
				else if(finalRefsetPercent > 100 || finalRefsetPercent < 0)
					JOptionPane.showMessageDialog(null, "Minimum Reference Set Size (Percentage) must be from 0-100");
				else if(maxIterCount != -1 && maxIterCount < 0)
					JOptionPane.showMessageDialog(null, "Max Iterations must be noted as -1 or greater than 0");
				else if(initialRefset == 's' && refsetPath.isEmpty())
					JOptionPane.showMessageDialog(null, "Reference Set Path is Empty.");
				else
				{
					cmdLine  = "-i " + index;
					cmdLine += " -g " + rcc;
					cmdLine += " -d " + divisor;
					cmdLine += " -p " + finalRefsetPercent;
					cmdLine += " -f " + filePath;
					cmdLine += " -m " + maxIterCount;
					if(initialRefset == 'r')
						cmdLine += " -r " + roundrobinRNN + " " + roundrobinISS + " " + roundrobinRandom;
					else if (initialRefset == 's')
						cmdLine += " -s " + refsetPath;
					
					ready = true;
				}
			}
		});
		paneButtons.add(computeButton);
		exitButton.addActionListener(new ActionListener()
		{
			public void actionPerformed(ActionEvent e)
			{
				System.exit(0);
			}
		});
		paneButtons.add(exitButton);
		
		// ----------------- finally
		frame.add(paneInOut);
		frame.add(paneIndex);
		frame.add(paneRefsetStart);
//		frame.add(paneTweaks);
		frame.add(paneButtons);
		
		frame.pack();
	}
	
	/**
	 * accessor for command line
	 * @return generated command line
	 */
	public String cmdline()
	{
		// reset whether command line has finished preparation
		ready = false;
		// return generated command line
		return cmdLine;
	}
	/**
	 * 
	 * @return whether the command line has been generated
	 */
	public boolean readyToRead()
	{
		// whether the command line is ready to be accessed
		return ready;
	}
	/**
	 * make the GUI visible
	 */
	public void setVisible()
	{
		// set the frame of the GUI to visible
		frame.setVisible(true);
	}
}


/**
 * 
 * @author www.java2s.com
 * http://www.java2s.com/Code/JavaAPI/javax.swing/JFileChoosersetFileFilterFileFilterfilter.htm
 *
 */
class ExtensionFileFilter extends FileFilter {
	  String description;
	  String extensions[];

	  public ExtensionFileFilter(String description, String extension) {
	    this(description, new String[] { extension });
	  }

	  public ExtensionFileFilter(String description, String extensions[]) {
	    if (description == null) {
	      this.description = extensions[0];
	    } else {
	      this.description = description;
	    }
	    this.extensions = (String[]) extensions.clone();
	    toLower(this.extensions);
	  }

	  private void toLower(String array[]) {
	    for (int i = 0, n = array.length; i < n; i++) {
	      array[i] = array[i].toLowerCase();
	    }
	  }

	  public String getDescription() {
	    return description;
	  }

	  public boolean accept(File file) {
	    if (file.isDirectory()) {
	      return true;
	    } else {
	      String path = file.getAbsolutePath().toLowerCase();
	      for (int i = 0, n = extensions.length; i < n; i++) {
	        String extension = extensions[i];
	        if ((path.endsWith(extension) && (path.charAt(path.length() - extension.length() - 1)) == '.')) {
	          return true;
	        }
	      }
	    }
	    return false;
	  }
	}
