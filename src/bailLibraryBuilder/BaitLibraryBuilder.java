package bailLibraryBuilder;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import commandLineInterface.CLI;
import commandLineInterface.CLIParseException;
import support.BioSequence;
import support.BlastHSP;
import support.BlastHit;
import support.BlastIteration;
import support.BlastReader;
import support.BlastXMLReader;
import support.FastaReader;


/**
 * 
 * @author steuernb
 *
 */
public class BaitLibraryBuilder {
	File selfblast;
	File fastaFile;
	
	Hashtable<String,StringBuilder> sequences;
	Vector<String> prefixOrder;
	
	HashSet<String> sequencesUsed;
	
	BufferedWriter out;
	
	
	int baitsWritten;
	
	
	
	public BaitLibraryBuilder(File selfblast, File fastaFile, File outputFile)throws IOException {
		super();
		this.selfblast = selfblast;
		this.fastaFile = fastaFile;
		
		
		sequences = new Hashtable<String,StringBuilder>();
		sequencesUsed = new HashSet<String>();
		
		
		this.out = new BufferedWriter(new FileWriter(outputFile));

		baitsWritten = 0;
	}
	
	
	

	
	
	
	
	public void execute(double identity, int minHSPLength, int baitLength,int overlap, double fractionOfBaitLengthToIgnore)throws IOException{
		this.readSequences();
		
		int counta = 0;
		
		
		BlastReader reader;
		reader = new BlastXMLReader(this.selfblast);
			
		
		
		
		 
		for(BlastIteration it = reader.readIteration(); it != null; it = reader.readIteration()){
		
			
			writeMyBaits(it, identity, minHSPLength, baitLength,overlap, fractionOfBaitLengthToIgnore);
			counta ++;
			if(counta%100==0){
				System.out.print(".");
				if(counta%10000==0){
					System.out.println();
				}
			}
			
		}
		System.out.println();
		reader.close();
		this.close();
	}
	
	
	
	public void execute(Vector<String> myorder, Hashtable<String,Double> minIdentities, int minHSPLength, int baitLength,int overlap, double fractionOfBaitLengthToIgnore)throws IOException{
		this.readSequences();
		
		int counta = 0;
		
		
		
		
		for(Enumeration<String> myenum = myorder.elements(); myenum.hasMoreElements();){
			String currentPrefix = myenum.nextElement();
			System.out.println(currentPrefix);
			
			BlastReader reader;
			reader = new BlastXMLReader(this.selfblast);
			
			
			
			for(BlastIteration it = reader.readIteration(); it != null; it = reader.readIteration()){
				if( !it.getQueryID().startsWith(currentPrefix)){
					continue;
				}
				writeMyBaits(it, minIdentities, minHSPLength, baitLength,overlap, fractionOfBaitLengthToIgnore);
				counta ++;
				if(counta%100==0){
					System.out.print(".");
					if(counta%10000==0){
						System.out.println();
					}
				}
				
			}
			System.out.println();
			reader.close();
			
			
			
			
		}
		
		
		this.close();
	}
	
	
	
	public void writeMaskedSources(File outputFile)throws IOException{
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		
		for(Enumeration< String> myenum= sequences.keys(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			BioSequence seq = new BioSequence(key, sequences.get(key).toString());
			out.write(seq.getFastaString(100));
		}
		
		out.close();
	}
	
	
	private void readSequences()throws IOException{
		
		
		FastaReader fastaReader = new FastaReader(fastaFile);
		for (BioSequence seq = fastaReader.readEntry(); seq != null; seq = fastaReader.readEntry()) {
			sequences.put(seq.getIdentifier(), new StringBuilder(seq.getSequence().toUpperCase()));
			
			
		}
		fastaReader.close();
	}
	
	
	
	
	
	
	private void maskSequence(BlastIteration it , double minIdentity, int minHSPLength, int baitStart, int baitEnd){
		for(Enumeration<BlastHit > myenum= it.getHits().elements(); myenum.hasMoreElements(); ){
			
			
			BlastHit hit = myenum.nextElement();
			if(hit.getQueryID().equalsIgnoreCase(hit.getHitID())){
				continue;
			}
			
			this.maskSequence(hit, minIdentity, minHSPLength, baitStart, baitEnd);
		}
	}
	
	
	private void maskSequence(BlastHit hit, double minIdentity, int minHSPLength, int baitStart, int baitEnd){
		if( this.sequencesUsed.contains(hit.getHitID())){
			return;
		}
		for(Enumeration<BlastHSP> myenum = hit.getHSPs().elements(); myenum.hasMoreElements();){
			BlastHSP hsp = myenum.nextElement();
			if( hsp.getIdentityPercentage()< minIdentity){
				continue;
			}
			
			
			
			// get overlap
			//find out if the sequence overlaps with query_hsp over at least minHSPLength bases
			int queryStart = Math.max(hsp.getQueryStart(), baitStart);
			int queryEnd   = Math.min(hsp.getQueryEnd(), baitEnd);
			
			if( queryEnd-queryStart < minHSPLength){
				continue;
			}
			
			int maskStart = Math.max(hsp.getHitStart(), hsp.getHitStart() + (queryStart - hsp.getQueryStart()));
			int maskEnd = Math.min(hsp.getHitEnd(),  hsp.getHitEnd()+ baitEnd - hsp.getQueryEnd() );
			
			if(hsp.getHitStrand() ==-1){  //if hit strand is -1 hitstart > hitend. so swap around
				 maskStart = Math.max(hsp.getHitEnd(), hsp.getHitEnd() + (queryStart - hsp.getQueryStart()));
				 maskEnd = Math.min(hsp.getHitStart(), hsp.getHitEnd() + baitEnd - hsp.getQueryEnd());
				
			}
			
			StringBuilder builder = sequences.get(hsp.getHitName());
			for( int i = maskStart; i < maskEnd; i++ ){
				try{
					builder.setCharAt(i,  Character.toLowerCase(builder.charAt(i)));
				}catch (StringIndexOutOfBoundsException e){
					e.printStackTrace();
					System.out.println(hsp.getHitName() + "\t" +  builder.toString() ) ;
				}
				
			}
			sequences.put(hsp.getHitName(), builder);
		}
		
		
	}
	
	
	
	private int getNextBaitStartPosition(int start, StringBuilder sb, int baitLength, double fractionOfBaitLengthToIgnore){
		
	
		
		Pattern p = Pattern.compile("[ATGC]+");
		Matcher m = p.matcher(sb.toString());
		
		if(m.find(start)){
			
			if(m.end()-m.start()>=baitLength){
				return m.start();
			}
			
			
			
			int mystart = m.start();
			int myend = m.end();
			
			if(myend == sb.length() && myend-mystart < baitLength * fractionOfBaitLengthToIgnore){
				return sb.length()+1;
				
			}
			
			boolean startCanBeElongated = true;
			boolean endCanBeElongated = true;
			boolean workAtStart = true;
			while(myend-mystart <baitLength && (startCanBeElongated || endCanBeElongated)){
				if(workAtStart){
					if(startCanBeElongated && mystart >0 && Pattern.matches("[ATGCatgc]+", ""+sb.charAt(mystart-1))){
						
							mystart = mystart-1;
							
					}else{
							
						startCanBeElongated=false;
					}
					workAtStart = false;
				}else{
					if(endCanBeElongated && myend < sb.length()-1 && Pattern.matches("[ATGCatgc]+", ""+sb.charAt(myend+1))){
						
						myend = myend + 1;
					}else{
						endCanBeElongated = false;
					}
						
					
					workAtStart = true;
				}
				
				
			}
			
			
			
			if(myend - mystart >= baitLength){
				int upperCaseBases = 0;
				for( int i = mystart; i <mystart+baitLength; i++ ){
					if( Character.isUpperCase(sb.charAt(i) )){
						upperCaseBases++;
					}
				}
				if( upperCaseBases >baitLength * fractionOfBaitLengthToIgnore){
					return mystart;
				}
			}
				
			return getNextBaitStartPosition(m.end(), sb, baitLength, fractionOfBaitLengthToIgnore);
			
			
		}
		
		
		
		
		
		return sb.length()+1;
		
	}
	
	
	
	private void writeMyBaits(BlastIteration it, double minIdentity, int minHSPLength, int baitLength,int overlap, double fractionOfBaitLengthToIgnore)throws IOException{
		
		String id = it.getQueryID();
		
		StringBuilder sb = sequences.get(id);
		
		int start = (int) Math.round((Math.random() * baitLength));
		start = getNextBaitStartPosition(start, sb, baitLength, fractionOfBaitLengthToIgnore);
		while(start + baitLength < sb.length()){
			out.write(new BioSequence(id+"_"+start, sb.substring(start, start+baitLength)).getFastaString());
			this.baitsWritten++;
			maskSequence(it, minIdentity, minHSPLength, start, start+baitLength);
			start = start + baitLength-overlap;
			start = getNextBaitStartPosition(start, sb, baitLength, fractionOfBaitLengthToIgnore);
			
		}
		
		this.sequencesUsed.add(id);
		
	}
	
	private void writeMyBaits(BlastIteration it, Hashtable<String,Double> minIdentities, int minHSPLength, int baitLength,int overlap, double fractionOfBaitLengthToIgnore)throws IOException{
		
		String id = it.getQueryID();
		
		double minIdentity = 100.0;
		boolean found  = false;
		for(Enumeration<String> myenum = minIdentities.keys(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			if( id.startsWith(key)){
				minIdentity = minIdentities.get(key);
				found = true;
				break;
			}
		}
		if(!found){
			System.err.println("Warning: did not find prefix from "+id+" in list of min identities. Using 100% instead");
		}
		
		
		StringBuilder sb = sequences.get(id);
		
		int start = (int) Math.round((Math.random() * baitLength));
		
		start = getNextBaitStartPosition(start, sb, baitLength, fractionOfBaitLengthToIgnore);
		
		while(start + baitLength < sb.length()){
			out.write(new BioSequence(id+"_"+start, sb.substring(start, start+baitLength)).getFastaString());
			this.baitsWritten++;
			maskSequence(it, minIdentity, minHSPLength, start, start+baitLength);
			int lastStart = start;
			start = start + baitLength-overlap;
			
			
			start = getNextBaitStartPosition(start, sb, baitLength, fractionOfBaitLengthToIgnore);
			
			while( start  <= lastStart){
				lastStart++;
				start = getNextBaitStartPosition(lastStart, sb, baitLength, fractionOfBaitLengthToIgnore);
			}
			
			}
		
		
		this.sequencesUsed.add(id);
	}
	
	
	
	
	
	
	
	private void close()throws IOException{
		this.out.close();
	}
	
	public int getNumBaitsWritte(){
		return this.baitsWritten;
	}
	
	
	
	
	
	public static void checkReverseComplementaryBaits(File blast, File sequenceFile, File outputFile, int percentquerycoverage)throws IOException{
		Hashtable<String,Integer> h = new Hashtable<String,Integer>();
		BlastReader reader = new BlastXMLReader(blast);
		
		HashSet<String> wasfound = new HashSet<String>();
		HashSet<String> remove = new HashSet<String>();
		
		for(BlastIteration it = reader.readIteration(); it != null; it = reader.readIteration()){
			for(Enumeration<BlastHit> myenum = it.getHits().elements(); myenum.hasMoreElements();){
				BlastHSP hsp= myenum.nextElement().getHSPs().get(0);
				if(hsp.getHitStrand() == -1){
					
					String s = hsp.getQueryName() + "\t" + hsp.getHitName();
					if( wasfound.contains(s)){
						continue;
					}
					wasfound.add(hsp.getHitName()+ "\t" +hsp.getQueryName());
					System.out.println(hsp.getQueryName() + "\t" + hsp.getHitName() + "\t" + hsp.getPercentQueryCoverage());
					if(hsp.getPercentQueryCoverage() > percentquerycoverage){
						remove.add(hsp.getHitName());
					}
					
					int num = 0;
					if(h.containsKey(hsp.getQueryName())){
						num = h.get(hsp.getQueryName());
					}
					num++;
					h.put(hsp.getQueryName(), num);
					
				}
			}
		}
		for(Enumeration<String> myenum = h.keys();myenum.hasMoreElements();){
			String key = myenum.nextElement();
			int i = h.get(key);
			if(i>1){
				System.out.println(key + "\t" + i);
			}
		}
		System.out.println("\n\nREMOVED:\n");
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		FastaReader fastaReader = new FastaReader(sequenceFile);
		for (BioSequence seq = fastaReader.readEntry(); seq != null; seq = fastaReader.readEntry()) {
			if(remove.contains(seq.getIdentifier())){
				System.out.println(seq.getFastaString());
			}else{
				out.write(seq.getFastaString());
			}
			
		}
		fastaReader.close();
		out.close();
		
	}
	
	
	
	
	
	
	
	public static void mainTest(String[] args){
		File dir = new File("/Users/steuernb/Documents/projects/shortjobs/2015_07_29_alpha_gliadin");
		try {
			int baitLength = 120;
			int overlap = 0;
			double fractionOfBaitLenghToIgnrore = 1.5;
			double identity = 80.0;
			int minHSPLength = 100;
			
			File inputFile = new File(dir,"all_alpha_gliadins.fasta");
			File inputBlast = new File(dir,"alpha_gliadins_selfblast.parsed");
					
			File outputFile = new File(dir, "test4.fasta");
			BaitLibraryBuilder builder = new BaitLibraryBuilder(inputBlast, inputFile, outputFile);
			builder.execute(identity, minHSPLength, baitLength, overlap, fractionOfBaitLenghToIgnrore);
			builder.writeMaskedSources(new File(dir ,"test4.masked.fasta"));
			builder.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
	}
	
	
	
	

	public static void main(String[] args){
		
		try {
		
			
			System.out.println("java -jar BaitLibraryBuilder.jar -i inputFile -o outputFile, -b blastFile");
			
			
			
			
			
			CLI cli = new CLI();
			
			cli.parseOptions(args);
			
			
			
			
			if( !cli.hasOption("i") || !cli.hasOption("o") || !cli.hasOption("b")) {
				throw new CLIParseException("parameters -i -o and -b are mandatory");
			}
				
			File inputFile = new File(cli.getArg("i"));
			File outputFile = new File(cli.getArg("o"));
			File blastFile = new File(cli.getArg("b"));
			double identity = 95.0;
			
			int overlap = 0;
			int baitlength = 120;
			int hsplength = 100;
			double minFractionOfBaitLengthToIgnore = 1.5;
			
			
			if (cli.hasOption("v")) {
				overlap = Integer.parseInt(cli.getArg("v"));
			}	
			
			if (cli.hasOption("a")) {
				baitlength = Integer.parseInt(cli.getArg("a"));
			}
			
			if (cli.hasOption("s")) {
				hsplength = Integer.parseInt(cli.getArg("s"));
			}
			
			if (cli.hasOption("f")) {
				minFractionOfBaitLengthToIgnore = Double.parseDouble(cli.getArg("f"));
			}
			
				
				
			Hashtable<String,Double> identities = null;
			Vector<String> order = null;
			if( cli.hasOption("d")){
				File identityFile = new File(cli.getArg("d"));
				if( identityFile.exists()){
					identities = new Hashtable<String,Double>();
					order = new Vector<String>();
					BufferedReader in = new BufferedReader(new FileReader(identityFile));
					for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
						identities.put(inputline.split("\t")[0], Double.parseDouble(inputline.split("\t")[1]));
						order.add(inputline.split("\t")[0]);
					}
					in.close();
				}else{
					try {
						identity = Double.parseDouble(cli.getArg("d"));
					} catch (NumberFormatException e) {
						throw new CLIParseException("The argument for parameter d is neither a float number nor a file");
					}
				}
				
			}
			BaitLibraryBuilder builder = new BaitLibraryBuilder(blastFile, inputFile, outputFile);
			
			if(identities != null){
				builder.execute(order, identities, hsplength, baitlength, overlap, minFractionOfBaitLengthToIgnore);
			}else{
				builder.execute(identity, hsplength, baitlength, overlap, minFractionOfBaitLengthToIgnore);
			}
			
			if( cli.hasOption("m")) {
				File maskedOutput = new File(cli.getArg("m"));
				builder.writeMaskedSources(maskedOutput);
			}
			
		
			builder.close();
		}catch(Exception e) {
			e.printStackTrace();
		}
		
		
		
		
	
	}
	
	
	

	
	
	
	
	
	
	
	
	
	
	
	
	

	

	
	
	
	
	
	
	
	
}
