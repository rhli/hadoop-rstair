/**
 * @file NativeStair.h
 * @brief Native C implementation of stair codes
 * @author Mingqiang Li, mingqiangli.cn@gmail.com
 * @version 0.1
 * @date 2014-12-04
 * @note revised on 2015-01-09
 */
#ifndef NativeStair_H_INCLUDED
#define NativeStair_H_INCLUDED

#include <arpa/inet.h>
#include <time.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

/* Added by RH Dec 24th 2014 begins */
#include <pthread.h>
/* Added by RH Dec 24th 2014 ends */

#include "org_apache_hadoop.h"
#include "org_apache_hadoop_util_NativeStair.h"

#include "stair.h"

/**
 * @brief Java_org_apache_hadoop_util_NativeStair_nativeInit - initialize all related variables
 *
 * @param env - JNI interface pointer
 * @param clazz - Java class object
 * @param stripeSize - number of data symbols in a stripe
 * @param paritySize - number of parity symbols in a stripe
 * @param numLocalGroup - number of local groups
 * @param sizeLocalGroup - size of each local group
 * @param numParityGroup - number of parity groups
 * @param localParitySize - number of local parity symbols in each local group
 * @param globalParityVector - global parity vector
 * @param globalParityVectorLen - length of global parity vector
 *
 * @return - a boolean value to indicate if the operation is successful
 */
JNIEXPORT jboolean JNICALL Java_org_apache_hadoop_util_NativeStair_nativeInit
(JNIEnv *env, jclass clazz, jint stripeSize, jint paritySize, jint numLocalGroup, jint sizeLocalGroup, 
 jint numParityGroup, jint localParitySize, jintArray globalParityVector, jint globalParityVectorLen);

/**
 * @brief Java_org_apache_hadoop_util_NativeStair_nativeCleanup - clean up all allocated space
 *
 * @param env - JNI interface pointer
 * @param clazz - Java class object
 */
JNIEXPORT void JNICALL Java_org_apache_hadoop_util_NativeStair_nativeCleanup
(JNIEnv *env, jclass clazz);

/**
 * @brief Java_org_apache_hadoop_util_NativeStair_nativeSymbolSize - get symbol size
 *
 * @param env - JNI interface pointer
 * @param clazz - Java class object
 *
 * @return - symbol size
 */
JNIEXPORT int JNICALL Java_org_apache_hadoop_util_NativeStair_nativeSymbolSize
(JNIEnv *env, jclass clazz);

/**
 * @brief Java_org_apache_hadoop_util_NativeStair_nativeEncodeBulk - encode data into parity in bulk
 *
 * @param env - JNI interface pointer
 * @param clazz - Java class object
 * @param inputBuffers - input buffers that store data 
 * @param numInputBuffers - number of input buffers
 * @param outputBuffers - output buffers that store parity 
 * @param numOutputBuffers - number of output buffers
 * @param dataLength - length of each data/parity block
 *
 * @return - a boolean value to indicate if the operation is successful
 */
JNIEXPORT jboolean JNICALL Java_org_apache_hadoop_util_NativeStair_nativeEncodeBulk
(JNIEnv *env, jclass clazz, jobjectArray inputBuffers, jint numInputBuffers, 
 jobjectArray outputBuffers, jint numOutputBuffers, jint dataLength);

/**
 * @brief Java_org_apache_hadoop_util_NativeStair_nativeLocationsToReadForDecode - determine locations 
 *                 to be read for decoding
 *
 * @param env - JNI interface pointer
 * @param clazz - Java class object
 * @param erasedLocationsBitmap - bitmap of erased locations
 * @param locationsToReadBitmap - bitmap of locations to be read
 *
 * @return - a boolean value to indicate if the operation is successful
 *
 * fixed by RH Dec 15th 2014 
 * We need to return the array instead of modifying directly
 */
//JNIEXPORT jboolean JNICALL Java_org_apache_hadoop_util_NativeStair_nativeLocationsToReadForDecode
//(JNIEnv *env, jclass clazz, jintArray erasedLocationsBitmap, jintArray locationsToReadBitmap);
JNIEXPORT jintArray JNICALL Java_org_apache_hadoop_util_NativeStair_nativeLocationsToReadForDecode
(JNIEnv *env, jclass clazz, jintArray erasedLocationsBitmap, jintArray locationsToReadBitmap);
/* fixed by RH Dec 15th 2014 ends */

/**
 * @brief Java_org_apache_hadoop_util_NativeStair_nativeDecodeBulk - decode erased blocks in bulk
 *
 * @param env - JNI interface pointer
 * @param clazz - Java class object
 * @param inputBuffers - input buffers that store the whole stripe 
 * @param numInputBuffers - number of input buffers
 * @param outputBuffers - output buffers that store the decoded blocks 
 * @param numOutputBuffers - number of output buffers
 * @param erasedLocationsArray - an array storing erased locations
 * @param numErasedLocations - number of erased locations
 * @param dataLength - length of each data/parity block
 *
 * @return - a boolean value to indicate if the operation is successful
 */
JNIEXPORT jboolean JNICALL Java_org_apache_hadoop_util_NativeStair_nativeDecodeBulk
(JNIEnv *env, jclass clazz, jobjectArray inputBuffers, jint numInputBuffers, jobjectArray outputBuffers,
 jint numOutputBuffers, jbooleanArray erasedLocationsArray, jint numErasedLocations, jint dataLength);

#endif
