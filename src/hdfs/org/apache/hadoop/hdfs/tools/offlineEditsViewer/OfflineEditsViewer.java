/**
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.apache.hadoop.hdfs.tools.offlineEditsViewer;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Map;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.conf.Configured;
import org.apache.hadoop.util.Tool;
import org.apache.hadoop.util.ToolRunner;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

/**
 * This class implements an offline edits viewer, tool that
 * can be used to view edit logs.
 */
public class OfflineEditsViewer extends Configured implements Tool {

  private EditsLoader  editsLoader;
  private final static String defaultProcessor = "xml";

  /**
   * Set editsLoader
   *
   * @param editsLoader EditsLoader
   */
  private void setEditsLoader(EditsLoader editsLoader) {
    this.editsLoader = editsLoader;
  }

  /**
   * Process EditLog file.
   *
   * @param visitor use this visitor to process the file
   */
  public void go(EditsVisitor visitor) throws IOException  {
    setEditsLoader(EditsLoader.LoaderFactory.getLoader(visitor));
    editsLoader.loadEdits();
  }

  /**
   * Print help.
   */  
  private void printHelp() {
    String summary =
      "Usage: bin/hdfs oev [OPTIONS] -i INPUT_FILE -o OUTPUT_FILE\n" +
      "Offline edits viewer\n" +
      "Parse a Hadoop edits log file INPUT_FILE and save results\n" +
      "in OUTPUT_FILE.\n" +
      "Required command line arguments:\n" +
      "-i,--inputFile <arg>   edits file to process, xml (case\n" +
      "                       insensitive) extension means XML format,\n" +
      "                       any other filename means binary format\n" +
      "-o,--outputFile <arg>  Name of output file. If the specified\n" +
      "                       file exists, it will be overwritten,\n" +
      "                       format of the file is determined\n" +
      "                       by -p option\n" +
      "\n" + 
      "Optional command line arguments:\n" +
      "-p,--processor <arg>   Select which type of processor to apply\n" +
      "                       against image file, currently supported\n" +
      "                       processors are: binary (native binary format\n" +
      "                       that Hadoop uses), xml (default, XML\n" +
      "                       format), stats (prints statistics about\n" +
      "                       edits file)\n" +
      "-h,--help              Display usage information and exit\n" +
      "-v,--verbose           More verbose output, prints the input and\n" +
      "                       output filenames, for processors that write\n" +
      "                       to a file, also output to screen. On large\n" +
      "                       image files this will dramatically increase\n" +
      "                       processing time (default is false).\n";


    System.out.println(summary);
    System.out.println();
    ToolRunner.printGenericCommandUsage(System.out);
  }

  /**
   * Build command-line options and descriptions
   *
   * @return command line options
   */
  public static Options buildOptions() {
    Options options = new Options();

    // Build in/output file arguments, which are required, but there is no 
    // addOption method that can specify this
    OptionBuilder.isRequired();
    OptionBuilder.hasArgs();
    OptionBuilder.withLongOpt("outputFilename");
    options.addOption(OptionBuilder.create("o"));
    
    OptionBuilder.isRequired();
    OptionBuilder.hasArgs();
    OptionBuilder.withLongOpt("inputFilename");
    options.addOption(OptionBuilder.create("i"));
    
    options.addOption("p", "processor", true, "");
    options.addOption("v", "verbose", false, "");
    options.addOption("h", "help", false, "");

    return options;
  }

  /**
   * Main entry point for ToolRunner (see ToolRunner docs)
   *
   * @param argv The parameters passed to this program.
   * @return 0 on success, non zero on error.
   */
  @Override
  public int run(String[] argv) throws Exception {
    int exitCode = 0;

    Options options = buildOptions();
    if(argv.length == 0) {
      printHelp();
      return -1;
    }

    CommandLineParser parser = new PosixParser();
    CommandLine cmd;

    try {
      cmd = parser.parse(options, argv);
    } catch (ParseException e) {
      System.out.println(
        "Error parsing command-line options: " + e.getMessage());
      printHelp();
      return -1;
    }

    if(cmd.hasOption("h")) { // print help and exit
      printHelp();
      return -1;
    }

    boolean printToScreen    = false;
    String inputFilenameArg  = cmd.getOptionValue("i");
    String outputFilenameArg = cmd.getOptionValue("o");
    String processor         = cmd.getOptionValue("p");
    if(processor == null) { processor = defaultProcessor; }

    if(cmd.hasOption("v")) { // print output to screen too
      printToScreen = true;
      System.out.println("input  [" + inputFilenameArg  + "]");
      System.out.println("output [" + outputFilenameArg + "]");
    }

    try {
      go(EditsVisitorFactory.getEditsVisitor(
        outputFilenameArg,
        processor,
        TokenizerFactory.getTokenizer(inputFilenameArg),
        printToScreen));
    } catch (EOFException e) {
      System.err.println("Input file ended unexpectedly. Exiting");
      exitCode = 255;
    } catch(IOException e) {
      System.err.println("Encountered exception. Exiting: " + e.getMessage());
      exitCode = 1;
    }

    return exitCode;
  }

  /**
   * main() runs the offline edits viewer using ToolRunner
   *
   * @param argv Command line parameters.
   */
  public static void main(String[] argv) throws Exception {
    int res = ToolRunner.run(new OfflineEditsViewer(), argv);
    System.exit(res);
  }
}
