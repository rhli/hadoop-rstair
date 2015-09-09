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

package org.apache.hadoop.raid;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.Set;
import java.util.HashSet;
import java.util.List;

import java.io.IOException;
import java.lang.IllegalArgumentException;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;

import org.apache.hadoop.util.NativeStair;

/**
 * This StairCode class implement a transposed stair code
 * @author Mingqiang Li
 */
public class StairCode extends ErasureCode {
  public static final Log LOG = LogFactory.getLog(StairCode.class);
  
  /**
   * Total number of data symbols in each stripe
   */
  private int stripeSize;
  
  /**
   * Total number of parity symbols in each stripe
   */
  private int paritySize;

  /**
   * Number of local groups
   */
  private int numLocalGroup;
  
  /**
   * Size of each local group
   */
  private int sizeLocalGroup;  

  /**
   * Number of parity groups
   */
  private int numParityGroup;
  
  /**
   * Number of local parity symbols in each local group
   */
  private int localParitySize;
  
  /**
   * Global parity vector
   */
  private int[] globalParityVector;

  public StairCode(int stripeSize, int paritySize, String jsonStr) {
    init(stripeSize, paritySize, jsonStr);
  }
  
  @Deprecated
  public StairCode(int stripeSize, int paritySize) {
    init(stripeSize, paritySize);
  }

  public StairCode() {
  }
  
  /**
   * Parse the input value of globalParityVector from globalParityVectorStr.
   * @param globalParityVectorStr
   *          The string that stores globalParityVector.
   * @return An array that stores globalParityVector. 
   */
  private int[] parseglobalParityVector(String globalParityVectorStr) 
  throws IllegalArgumentException {	
	  // The expected format of globalParityVectorStr is as follows:
	  // 	"(e_0,e_1,...,e_{numLocalGroup-1})".
	  
	  // filter out blank spaces in globalParityVectorStr
	  String globalParityVectorStrNeat = globalParityVectorStr.replace(" ", "");
	  
	  // change globalParityVectorStr into "e_0,e_1,...,e_{numLocalGroup-1}"
	  if ((!globalParityVectorStrNeat.startsWith("(")) || 
			  (!globalParityVectorStrNeat.endsWith(")"))) {
		  throw new IllegalArgumentException("Invalid globalParityVector<" + 
			  globalParityVectorStr + ">. Pls check if it is put into (xxx)!");
	  }	
    /* fixed by RH Dec 12th 2014 begins 
     * remove "(" and ")"
     */
	  globalParityVectorStrNeat = globalParityVectorStrNeat.replace("(", "");
	  globalParityVectorStrNeat = globalParityVectorStrNeat.replace(")", "");
		//globalParityVectorStrNeat = globalParityVectorStrNeat.substring(1, 
		//	  globalParityVectorStrNeat.length() - 1);
    /* fixed by RH Dec 12th 2014 ends */
	  
	  // check the correction of the length and each e_i (i=0,1,...,numLocalGroup-1) 
	  // in globalParityVectorStr 
	  int globalParityVectorLen = 0;
    /* fixed by RH Jan 15th 2014 begins 
     * consider one element array */
    if(!globalParityVectorStrNeat.contains(",")){
		  int element = Integer.valueOf(globalParityVectorStrNeat).intValue();
		  if ((element <= 0) || (element > (this.sizeLocalGroup - this.localParitySize))) {
		    throw new IllegalArgumentException("Invalid globalParityVector<" + 
		      globalParityVectorStr + ">. Pls check the value range of its elements!");
		  }
		  // check the length
		  globalParityVectorLen++;
		  if (globalParityVectorLen > this.numLocalGroup) {
		    throw new IllegalArgumentException("Invalid globalParityVector<" + 
		      globalParityVectorStr + ">. Pls check the number of its elements!");
		  }
    } else {
	    for (String str : globalParityVectorStrNeat.split(",")) {
		    // check e_i
		    int element = Integer.valueOf(str).intValue();
		    if ((element <= 0) || (element > (this.sizeLocalGroup - this.localParitySize))) {
		  	  throw new IllegalArgumentException("Invalid globalParityVector<" + 
		        globalParityVectorStr + ">. Pls check the value range of its elements!");
		    }
		    
		    // check the length
		    globalParityVectorLen++;
		    if (globalParityVectorLen > this.numLocalGroup) {
		  	  throw new IllegalArgumentException("Invalid globalParityVector<" + 
		        globalParityVectorStr + ">. Pls check the number of its elements!");
		    }
	    }
    }
    /* fixed by RH Jan 15th 2014 ends */ 
	  
	  // store globalParityVector into an array 	  
    int[] globalParityVectorTmp = new int[globalParityVectorLen];
    int i = 0;
	  for (String str : globalParityVectorStrNeat.split(",")) {
		  globalParityVectorTmp[i] = Integer.valueOf(str).intValue();
      //LOG.info("globalParityVectorTmp[" + i + "]:" + globalParityVectorTmp[i]);
		  i++;
	  }
	  
	  return globalParityVectorTmp;  
  }
  
  public void init(int stripeSize, int paritySize, String jsonStr) {	  
	  // get additional parameters for STAIR codes
    try {
    	JSONObject json = new JSONObject(jsonStr); 
      this.numLocalGroup = json.getInt("num_local_group"); 
      this.sizeLocalGroup = json.getInt("size_local_group"); 
      this.numParityGroup = json.getInt("num_parity_group"); 
      this.localParitySize = json.getInt("local_parity_length"); 
      this.globalParityVector = parseglobalParityVector(json.getString("global_parity_vector"));
    } catch (JSONException ex) {
      ex.printStackTrace();
    } catch (IllegalArgumentException ex) {
      System.out.println("Just caught an exception..." + ex.getMessage());
    }
    
    // do initialization
    //LOG.info("pre-initialize()");
    init(stripeSize, paritySize);
    //LOG.info("post-initialize()");
    
    // log the initialization info
    LOG.info(" Initialized " + StairCode.class +
             " stripeLength:" + stripeSize +
             " parityLength:" + paritySize +
             " numLocalGroup:" + this.numLocalGroup + 
             " sizeLocalGroup:" + this.sizeLocalGroup + 
             " numParityGroup:" + this.numParityGroup + 
             " localParitySize:" + this.localParitySize + 
             " globalParityVector:" + Arrays.toString(this.globalParityVector));
  }

  @Override
  public void init(int stripeSize, int paritySize) {
    this.stripeSize = stripeSize;
    this.paritySize = paritySize;
    
    int numGlobalParity = 0;
    for (int i = 0; i < this.globalParityVector.length; i++) {
    	numGlobalParity += this.globalParityVector[i];
    }
    
    if ((this.stripeSize + this.paritySize != this.numLocalGroup * this.sizeLocalGroup) ||
    		(this.paritySize != numGlobalParity + this.sizeLocalGroup * this.numParityGroup  + 
    		this.localParitySize * (this.numLocalGroup - this.numParityGroup))) {
    	LOG.error("Invalid stripeSize<" + stripeSize + 
    			"> & paritySize<" + paritySize + ">. Pls double-check them!");
    } else {
    	if (NativeStair.isAvailable()) { 
        if (!NativeStair.nativeInit(this.stripeSize, this.paritySize, this.numLocalGroup, 
              this.sizeLocalGroup, this.numParityGroup, this.localParitySize, 
              this.globalParityVector, this.globalParityVector.length)) { 
          LOG.info("Fail to initialize stair code"); 
        } 
      } else { 
        LOG.info("Cannot use the native C implementation of stair code");
   	 }
    }
  }

  @Override
  public int stripeSize() {
    return this.stripeSize;
  }

  @Override
  public int paritySize() {
    return this.paritySize;
  }

  @Override
  public int symbolSize() {
    return NativeStair.nativeSymbolSize();
  }

  @Deprecated
  public void encode(int[] message, int[] parity) {
  }
  
  /**
   * encodeBulk - encode data in bulk
   * @param inputs - input data blocks
   * @param outputs - output parity blocks
   */
  @Override
  public void encodeBulk(byte[][] inputs, byte[][] outputs) throws IOException {
    encodeBulk(inputs, outputs, true);
  }
  
  @Override
  public void encodeBulk(byte[][] inputs, byte[][] outputs, boolean useNative)
  throws IOException {
    assert(inputs.length == this.stripeSize);
    assert(outputs.length == this.paritySize);
    
    if (NativeStair.isAvailable() && useNative) {
      if (!NativeStair.encodeBulk(inputs, outputs)) {
	      throw new IOException("Fail to encode data in bulk");
      }
    } else {
      throw new IOException("Stair code is using native C implementation");
    }
  }
   
  /**
   * locationsToReadForDecode - figure out which locations need to be read 
   * 			to decode erased locations
   * @param erasedLocations - erased locations
   * @return - locations to read
   */
  @Override
  public List<Integer> locationsToReadForDecode(List<Integer> erasedLocations)
      throws TooManyErasedLocations {
    LOG.info("Erased locations: "+erasedLocations.toString());
    int bitmapLen = this.stripeSize + this.paritySize;
    int[] erasedLocationsBitmap = new int[bitmapLen];
    int[] locationsToReadBitmap = new int[bitmapLen];
    for (int i = 0; i < bitmapLen; i++) {
    	erasedLocationsBitmap[i] = 0;
    	locationsToReadBitmap[i] = 0;
    }
    
    for (int i = 0; i < erasedLocations.size(); i++) {
    	erasedLocationsBitmap[erasedLocations.get(i)] = 1;
    }
    
    //if (NativeStair.nativeLocationsToReadForDecode(erasedLocationsBitmap, locationsToReadBitmap)) {
    locationsToReadBitmap=
      NativeStair.nativeLocationsToReadForDecode(erasedLocationsBitmap,locationsToReadBitmap);
    if (locationsToReadBitmap==null) {
    	String locationsStr = ""; 
      for (Integer erasedLocation : erasedLocations) { 
        locationsStr += " " + erasedLocation; 
      }
    	throw new TooManyErasedLocations("Locations: " + locationsStr);
    }

    List<Integer> locationsToRead = new ArrayList<Integer>();
    for (int i = 0; i < bitmapLen; i++) {
      //LOG.info("locationsToReadBitmap[" + i + "]: "+ locationsToReadBitmap[i]);
    	if (locationsToReadBitmap[i] == 1) {
    		locationsToRead.add(i);
    	}
    }
    LOG.info("locationsToRead: "+ locationsToRead.toString());
    
    return locationsToRead;
  }

  /*
   * Performs Reed Solomon decoding, assuming that all positions not included in
   * erasedLocations are available.
   */
  @Deprecated
  public void decode(int[] data, int[] erasedLocations, int[] erasedValues) {
  }

  /**
   * decodeBulk - decode erased blocks in bulk
   * @param readBufs - a stripe containing erased blocks
   * @param writeBufs - decoded blocks
   * @param erasedLocations - location of erased blocks
   */
  @Override
  public void decodeBulk(byte[][] readBufs, byte[][] writeBufs, 
  	int[] erasedLocations)  throws IOException {
    decodeBulk(readBufs, writeBufs, erasedLocations, 0, readBufs[0].length);
  }

  @Override
  public void decodeBulk(byte[][] readBufs, byte[][] writeBufs, 
      int[] erasedLocations, int dataStart, int dataLen) throws IOException {
    //LOG.info("before decodeBulk" + System.currentTimeMillis());
    decodeBulk(readBufs, writeBufs, erasedLocations, dataStart, dataLen, true); 
    //LOG.info("after decodeBulk" + System.currentTimeMillis());
  }
  
  public void decodeBulk(byte[][] readBufs, byte[][] writeBufs, int[] erasedLocations,
      int dataStart, int dataLen, boolean useNative) throws IOException{
    //LOG.info("decodeBulk checkpoint1 " + System.currentTimeMillis());
    assert (readBufs.length == this.stripeSize + this.paritySize);
    assert (writeBufs.length >= erasedLocations.length);
    //LOG.info("decodeBulk checkpoint2 " + System.currentTimeMillis());

    if (erasedLocations.length == 0) {
      return;
    }
    
    //LOG.info("decodeBulk checkpoint3 " + System.currentTimeMillis());
    if (NativeStair.isAvailable() && useNative) {
      //LOG.info("decodeBulk checkpoint3.1 " + System.currentTimeMillis());
      if (!NativeStair.decodeBulk(readBufs, writeBufs, erasedLocations, dataStart, dataLen)) {
	      throw new IOException("Fail to decode data in bulk");
      }
    } else {
      throw new IOException("Stair code is using native C implementation");
    }
    //LOG.info("decodeBulk checkpoint4 " + System.currentTimeMillis());
  }
  
  @Override
  public void decodeOneBlock(byte[][] readBufs, byte[] decodeVec, int dataLen, 
      int[] erasedLocations, int decodeLocation, int decodePos, int decodeLen,
      boolean useNative) throws IOException {
    int numErasedLocation = erasedLocations.length; 
    if (numErasedLocation == 0) {
      return;
    }
    int pos = -1;
    for (int i = 0; i < numErasedLocation; i++) {
      if (erasedLocations[i] == decodeLocation) {
        pos = i;  
        break;
      }
    }
    if (pos == -1) {
      LOG.error("Location " + decodeLocation + " is not in the erasedLocation");
      return; 
    }
    byte[][] tmpBufs = new byte[numErasedLocation][];
    for (int i = 0; i < numErasedLocation; i++) {
      tmpBufs[i] = new byte[dataLen];
    }
    decodeBulk(readBufs, tmpBufs, erasedLocations, decodePos, decodeLen, useNative);
    System.arraycopy(tmpBufs[pos], decodePos, decodeVec, decodePos, decodeLen);
  }
}


 
