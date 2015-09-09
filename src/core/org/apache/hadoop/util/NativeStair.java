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
package org.apache.hadoop.util;

import java.nio.ByteBuffer;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/**
 * This NativeStair class provides JNI-based native Stair extensions
 * @author Mingqiang Li
 */
public class NativeStair {
  public static final Log LOG = LogFactory.getLog(NativeStair.class);

  /**
   * Return true if the JNI-based native Stair extensions are available.
   */
  public static boolean isAvailable() {
    return NativeCodeLoader.isNativeCodeLoaded();
  }

  /* Added by RH Dec 13th 2014, begins */
  //static {
  //  // Try to load native hadoop library and set fallback flag appropriately
  //  LOG.debug("Trying to load the custom-built native-hadoop library...");
  //  try {
  //    System.load("/usr/local/lib/libgf_complete.so");
  //    LOG.info("Loaded the native-gfcomplete library");
  //  } catch (Throwable t) {
  //    // Ignore failure to load
  //    LOG.debug("Failed to load native-gfcomplete with error: " + t);
  //  }
  //}
  /* Added by RH Dec 13th 2014, ends */
  
  public static native boolean nativeInit(int stripeSize, int paritySize, int numLocalGroup, 
		  int sizeLocalGroup, int numParityGroup, int localParitySize, 
		  int[] globalParityVector, int globalParityVectorLen);

  public static native void nativeCleanup();

  public static native int nativeSymbolSize();
  
  public static boolean encodeBulk(byte[][] inputs, byte[][] outputs) {
    ByteBuffer[] inputBuffers = new ByteBuffer[inputs.length];
    ByteBuffer[] outputBuffers = new ByteBuffer[outputs.length];
    int bufferLen = inputs[0].length;
	
    for (int i = 0; i < inputs.length; i++) {
      inputBuffers[i] = directify(inputs[i], 0, bufferLen);
    }
    for (int i = 0; i < outputs.length; i++) {
      outputBuffers[i] = ByteBuffer.allocateDirect(bufferLen);
    }
	
    //LOG.info("inputLength:" + inputs.length +
    //    " outputLength:" + outputs.length +
    //    " bufferLength:" + bufferLen);
    if (!nativeEncodeBulk(inputBuffers, inputs.length, outputBuffers, outputs.length, bufferLen)) { 
      LOG.info("nativeEncodeBulk returns false");
      return false;
    }
    //LOG.info("nativeEncodeBulk returns true");
	
    for (int i = 0; i < outputs.length; i++) {
      outputBuffers[i].get(outputs[i]);
    }

    return true;    
  }
  
  private static native boolean nativeEncodeBulk(ByteBuffer[] inputBuffers, int numInputBuffers,
      ByteBuffer[] outputBuffers, int numOutputBuffers, int dataLen);
  
  /* Fixed by RH Dec 15th 2014 begins */
  //public static native boolean nativeLocationsToReadForDecode
  //	(int[] erasedLocationsBitmap, int[] locationsToReadBitmap);
  public static native int[] nativeLocationsToReadForDecode
  	(int[] erasedLocationsBitmap, int[] locationsToReadBitmap);
  /* Fixed by RH Dec 15th 2014 ends */
  
  public static boolean decodeBulk(byte[][] readBufs, byte[][] writeBufs, 
      int[] erasedLocation, int dataStart, int dataLen) {
    //LOG.info("decodeBulk(): checkpoint1 " + System.currentTimeMillis());
    ByteBuffer[] inputBuffers = new ByteBuffer[readBufs.length];
    ByteBuffer[] outputBuffers = new ByteBuffer[writeBufs.length];
	
    //LOG.info("decodeBulk(): checkpoint1.1 " + System.currentTimeMillis());
    //LOG.info("decodeBulk(): checkpoint3 " + System.currentTimeMillis());
    /* Added by RH Dec 16th 2014 begins
     * We use -1 to mark placeholder */
    int erasedLen=0;
    int erasedLocationsBitmap[] = new int[20];
    int locationsToReadBitmap[] = new int[20];
    for (int i = 0; i < erasedLocation.length; i++) {
      //LOG.info("erasedLocation[" + i + "]: " + erasedLocation[i]);
      if(erasedLocation[i]!=-1) {
        //LOG.info("erasedLocationsBitmap[" + erasedLocation[i] + "]=1");
        erasedLocationsBitmap[erasedLocation[i]] = 1;
        //inputBuffers[i] = ByteBuffer.allocateDirect(dataLen);
        erasedLen++;
      }else{
        break;
      }
    }
    //LOG.info("decodeBulk(): checkpoint2 " + System.currentTimeMillis());
    //for (int i = 0; i < erasedLocation.length; i++) {
    //for (int i = 0; i < writeBufs.length; i++) {
    for (int i = 0; i < erasedLen; i++) {
      outputBuffers[i] = ByteBuffer.allocateDirect(dataLen);
    }
    /* TODO:!! This is hardcoded!! De-hard code after deadline */
    locationsToReadBitmap = 
      nativeLocationsToReadForDecode(erasedLocationsBitmap,locationsToReadBitmap);
    for (int i = 0; i < readBufs.length; i++) {
      //LOG.info("locationsToReadBitmap[" + i + "]=" + locationsToReadBitmap[i]);
      if(locationsToReadBitmap[i] == 1 || erasedLocationsBitmap[i] == 1) {
        //LOG.info("directify()" + i);
        inputBuffers[i] = directify(readBufs[i], dataStart, dataLen);
      } else {
        /**
         * modified by RH Jul 2nd, 2015 begins
         */
        //inputBuffers[i] = ByteBuffer.allocateDirect(dataLen);
        inputBuffers[i] = null;
        /**
         * modified by RH Jul 2nd, 2015 ends
         */
        //inputBuffers[i] = directify(readBufs[i], dataStart, dataLen);
      }
    }
    //LOG.info("decodeBulk(): checkpoint4 " + System.currentTimeMillis());
    /* Added by RH Dec 16th 2014 ends */
	
    /* Fixed by RH Dec 16th 2014 begins */
    //if (!nativeDecodeBulk(inputBuffers, readBufs.length, outputBuffers, writeBufs.length, 
    //      erasedLocation, erasedLocation.length, dataLen)) { 
    //  return false;
    //}
    //LOG.info("decodeBulk(): before nativeDecodeBulk()");
    if (!nativeDecodeBulk(inputBuffers, readBufs.length, outputBuffers, writeBufs.length, 
          erasedLocation, erasedLen, dataLen)) { 
      return false;
    }
    //LOG.info("decodeBulk(): checkpoint5 " + System.currentTimeMillis());
    /* Fixed by RH Dec 16th 2014 ends */
	
    //for (int i = 0; i < writeBufs.length; i++) {
    //for (int i = 0; i < writeBufs.length; i++) {
    for (int i = 0; i < erasedLen; i++) {
      //LOG.info(i);
      outputBuffers[i].get(writeBufs[i], dataStart, dataLen);
    }
    //LOG.info("decodeBulk(): checkpoint6 " + System.currentTimeMillis());

    return true;
  }
  
  private static native boolean nativeDecodeBulk(ByteBuffer[] inputBuffers, int numInputBuffers,
      ByteBuffer[] outputBuffers, int numOutputBuffers, int[] erasedLocation, int erasedLocationCount, 
      int dataLen);
  
  private static ByteBuffer directify(byte[] readBufs, int dataStart, int dataLen) {
    //LOG.info("directify starts: " + System.nanoTime());
    ByteBuffer newBuf = null;
    newBuf = ByteBuffer.allocateDirect(dataLen);
    newBuf.position(0);
    newBuf.mark();
    newBuf.put(readBufs, dataStart, dataLen);
    newBuf.reset();
    newBuf.limit(dataLen);
    //LOG.info("directify ends: " + System.nanoTime());
    return newBuf;
  }
  
  /* Added by RH Jan 22nd 2015 begins */
  //private static ByteBuffer directify2(byte[] readBufs, int dataStart, int dataLen) {
  //  //LOG.info("directify starts: " + System.nanoTime());
  //  ByteBuffer newBuf = null;
  //  newBuf = ByteBuffer.allocateDirect(dataLen);
  //  newBuf.position(0);
  //  newBuf.mark();
  //  newBuf.put(readBufs, dataStart, dataLen);
  //  newBuf.reset();
  //  newBuf.limit(dataLen);
  //  //LOG.info("directify ends: " + System.nanoTime());
  //  return newBuf;
  //}
  /* Added by RH Jan 22nd 2015 ends */
}
