/**
 * @file NativeStair.c
 * @brief Native C implementation of stair codes
 * @author Mingqiang Li, mingqiangli.cn@gmail.com
 * @version 0.1
 * @date 2014-12-04
 * @note revised on 2015-01-09
 */
#include "NativeStair.h"

int n_stair;
int r_stair;
int m_stair;
int l_stair;
int m_partial_stair;
int *error_vector_stair = NULL;

bool *parity_layout_stair = NULL;

/*gf_t object for accelerating GF calculation*/
int word_size_in_gf;
gf_t gf_obj;

int *row_coding_matrix_stair = NULL;
int *column_coding_matrix_stair = NULL;

int best_encoding_method;

/* Added by RH Dec 22nd 2014 begins 
 * We add a mapping from the coding implementation and the HDFS inplementation */
int num_of_parity_symbols;
int* stripeDataMap;
int* stripeParityMap;
int* parityStripeMap;
int* dataStripeMap;
/* Added by RH Dec 22nd 2014 ends */

/* Added by RH Dec 24th 2014 begins
 * initialized flag */
bool initialized=false;
pthread_mutex_t* initMutex=NULL;
/* Added by RH Dec 24th 2014 ends */

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
 jint numParityGroup, jint localParitySize, jintArray globalParityVector, jint globalParityVectorLen) {
	//fprintf(stderr,"nativeInit(): checkpoint 1\n");
  /* Added by RH Jan 17th 2015 begins */
  //if(initMutex==NULL){
  //  pthread_mutex_init(initMutex,NULL);
  //}
  //pthread_mutex_lock(initMutex);
  //if(initialized){
  //  pthread_mutex_unlock(initMutex);
  //  return JNI_TRUE;
  //}
  /* Added by RH Jan 17th 2015 ends */
	int *global_parity_vector = (int *)(*env)->GetIntArrayElements(env, globalParityVector, 0);
	int i, j;

	/*set and check configuration parameters*/
	n_stair = numLocalGroup;
	if (n_stair <= 0) {
		fprintf(stderr, "Error: bad n_stair (%d), which should be > 0!\n", n_stair);
		return JNI_FALSE;
	}
	r_stair = sizeLocalGroup;
	if (r_stair <= 0) {
		fprintf(stderr, "Error: bad r_stair (%d), which should be > 0!\n", r_stair);
		return JNI_FALSE;
	}
	m_stair = numParityGroup;
	if ((m_stair < 0) || (m_stair >= n_stair)) {
		fprintf(stderr, "Error: bad m_stair (%d), which should be in [0, %d)!\n", m_stair, n_stair);
		return JNI_FALSE;
	}	
	l_stair = localParitySize;
	if ((l_stair < 0) || (l_stair >= r_stair)) {
		fprintf(stderr, "Error: bad l_stair (%d), which should be in [0, %d)!\n", l_stair, r_stair);
		return JNI_FALSE;
	}
	m_partial_stair = globalParityVectorLen;	
	if ((m_partial_stair <= 0) || (m_partial_stair >= n_stair - m_stair)) {
		fprintf(stderr, "Error: bad m_partial_stair (%d), which should be in (0, %d)!\n", 
				m_partial_stair, n_stair - m_stair);
		return JNI_FALSE;
	}	
	error_vector_stair = TypeAlloc(int, m_partial_stair);
	for (i = 0; i < m_partial_stair; i++) {	
		error_vector_stair[i] = global_parity_vector[i];
	}
	for (i = 0; i < m_partial_stair; i++) {	
		if ((error_vector_stair[i] <= 0) || (error_vector_stair[i] >= r_stair - l_stair)) {
			fprintf(stderr, "Error: bad %d-th element (%d) in the error vector, which should be in (0, %d)!\n", 
					i, error_vector_stair[i], r_stair - l_stair);
			return JNI_FALSE;
		}
		if (i > 0) {
			if (error_vector_stair[i] < error_vector_stair[i - 1]) {
				fprintf(stderr, "Error: the %d-th element (%d) in the error vector is not in \
						monotonically increasing order!\n", i, error_vector_stair[i]);
				return JNI_FALSE;
			}
		}
	}

	/* Added by RH Dec 22nd begins */
	num_of_parity_symbols = 0;
	num_of_parity_symbols += ((r_stair - l_stair) * m_stair);
	for (i = 0; i < m_partial_stair; i++) {
		num_of_parity_symbols += error_vector_stair[i];
	}
	num_of_parity_symbols += (n_stair * l_stair);
	/* Added by RH Dec 22nd ends */

	InitParityLayout(n_stair, r_stair, m_stair, l_stair, 
			m_partial_stair, error_vector_stair, &parity_layout_stair,
			&parityStripeMap,&stripeParityMap,&dataStripeMap,&stripeDataMap);

	if(!InitGF(n_stair, r_stair, m_stair, l_stair, m_partial_stair, 
				error_vector_stair, &word_size_in_gf, &gf_obj)) {
		return JNI_FALSE;
	}
	/*return JNI_TRUE;*/
	//fprintf(stderr,"nativeInit(): checkpoint 4\n");

	InitCodingMatrices(gf_obj, n_stair, r_stair, m_stair, l_stair, m_partial_stair,
			error_vector_stair, &row_coding_matrix_stair, &column_coding_matrix_stair);
	/*fprintf("nativeInit(): row_coding_matrix_stair %x\n",row_coding_matrix_stair);*/
	//if(row_coding_matrix_stair==NULL){
	//	fprintf(stderr,"NativeStair.c: row_coding_matrix_stair not initialized!!\n");
	//}

	best_encoding_method = ChooseEncodingMethod(n_stair, r_stair, m_stair, 
			l_stair, m_partial_stair, error_vector_stair, parity_layout_stair);

  /* Commented by RH Jan 21th 2015 begins */
	//fprintf(stderr, "\nSucceed in initializing a transposed stair code with the 
  //        following parameters:-)\n");
	//fprintf(stderr, "Parameters:\n");
	//fprintf(stderr, "      n_stair: %d\n", n_stair);			
	//fprintf(stderr, "      r_stair: %d\n", r_stair);							
	//fprintf(stderr, "      m_stair: %d\n", m_stair);								
	//fprintf(stderr, "      l_stair: %d\n", l_stair);					
	//fprintf(stderr, "      error_vector_stair: (");
	//for (i = 0; i < m_partial_stair - 1; i++) {		
	//	fprintf(stderr, "%d, ", error_vector_stair[i]);	
	//}					
	//fprintf(stderr, "%d)\n", error_vector_stair[m_partial_stair - 1]);	
	//fprintf(stderr, "      word_size_in_gf: %d\n", word_size_in_gf);					
	//fprintf(stderr, "      row_coding_matrix_stair: (see below)\n");				
	//for (i = 0; i < m_stair + m_partial_stair; i++) {			
	//	fprintf(stderr, "         | ");			
	//	for (j = 0; j < n_stair - m_stair; j++) {				
	//		fprintf(stderr, "%3d ", row_coding_matrix_stair[(n_stair - m_stair) * i + j]);			
	//	}			
	//	fprintf(stderr, "|\n");		
	//}							
	//fprintf(stderr, "      column_coding_matrix_stair: (see below)\n");				
	//for (i = 0; i < l_stair + error_vector_stair[m_partial_stair - 1]; i++) {			
	//	fprintf(stderr, "         | ");			
	//	for (j = 0; j < r_stair - l_stair; j++) {				
	//		fprintf(stderr, "%3d ", column_coding_matrix_stair[(r_stair - l_stair) * i + j]);			
	//	}			
	//	fprintf(stderr, "|\n");
	//}	
	//if (best_encoding_method == STD_METHOD) {
	//	fprintf(stderr, "      best_encoding_method: STD_METHOD\n");	
	//}

	//if (best_encoding_method == UP_METHOD) {
	//	fprintf(stderr, "      best_encoding_method: UP_METHOD\n");	
	//}

	//if (best_encoding_method == DOWN_METHOD) {
	//	fprintf(stderr, "      best_encoding_method: DOWN_METHOD\n");	
	//}			
	//fprintf(stderr, "      decoding method: DOWN_METHOD\n");		
	//fprintf(stderr, "\n");	
  /* Commented by RH Jan 21th 2015 ends */

  /* Added by RH Jan 17th 2015 begins */
  //initialized=true;
  //pthread_mutex_unlock(initMutex);
  /* Added by RH Jan 17th 2015 ends */

	return JNI_TRUE;
}

/**
 * @brief Java_org_apache_hadoop_util_NativeStair_nativeCleanup - clean up all allocated space
 *
 * @param env - JNI interface pointer
 * @param clazz - Java class object
 */
JNIEXPORT void JNICALL Java_org_apache_hadoop_util_NativeStair_nativeCleanup
(JNIEnv *env, jclass clazz) {
	free(error_vector_stair);
  error_vector_stair=NULL;
	free(parity_layout_stair);
  parity_layout_stair=NULL;
  free(stripeDataMap);
  stripeDataMap=NULL;
  free(stripeParityMap);
  stripeDataMap=NULL;
  free(parityStripeMap);
  stripeDataMap=NULL;
  free(dataStripeMap);
  stripeDataMap=NULL;

	/*free the gf_t object*/		
	gf_free(&gf_obj, 1);

	free(row_coding_matrix_stair);
  row_coding_matrix_stair = NULL;
	free(column_coding_matrix_stair);
  column_coding_matrix_stair = NULL;

	fprintf(stderr, "\nWe have cleaned up all space allocated for the stair code!\n");
	fprintf(stderr, "\n");	
}

/**
 * @brief Java_org_apache_hadoop_util_NativeStair_nativeSymbolSize - get symbol size
 *
 * @param env - JNI interface pointer
 * @param clazz - Java class object
 *
 * @return - symbol size
 */
JNIEXPORT int JNICALL Java_org_apache_hadoop_util_NativeStair_nativeSymbolSize
(JNIEnv *env, jclass clazz) {
	return word_size_in_gf;
}

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
 jobjectArray outputBuffers, jint numOutputBuffers, jint dataLength) {
	char* data[numInputBuffers];
	char* coding[numOutputBuffers];
	int i;

	for (i = 0; i < numInputBuffers; i++) {
		jobject j_inputBuffer = (*env)->GetObjectArrayElement(env, inputBuffers, i);
		data[i] = (char*)(*env)->GetDirectBufferAddress(env, j_inputBuffer);
	}

	for (i = 0; i < numOutputBuffers; i++) {
		jobject j_outputBuffer = (*env)->GetObjectArrayElement(env, outputBuffers, i);
		coding[i] = (char*)(*env)->GetDirectBufferAddress(env, j_outputBuffer);
		//memset(coding[i], 0, dataLength);
	}

	if ((dataLength <= 0) || (dataLength % (word_size_in_gf/8) != 0)) {
		fprintf(stderr, "Error: bad dataLength (%d), which should be a positive multiple of %d!\n", 
				dataLength, word_size_in_gf/8);
		return JNI_FALSE;
	}

	if(!EncodeBulk(gf_obj, n_stair, r_stair, m_stair, l_stair, m_partial_stair, 
				error_vector_stair, parity_layout_stair, best_encoding_method, 
				row_coding_matrix_stair, column_coding_matrix_stair, 
				data, coding, dataLength)) {
		return JNI_FALSE;
	}

	return JNI_TRUE;
}

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
 */
JNIEXPORT jintArray JNICALL Java_org_apache_hadoop_util_NativeStair_nativeLocationsToReadForDecode
(JNIEnv *env, jclass clazz, jintArray erasedLocationsBitmap, jintArray locationsToReadBitmap) {
	//JNIEXPORT jboolean JNICALL Java_org_apache_hadoop_util_NativeStair_nativeLocationsToReadForDecode
	//(JNIEnv *env, jclass clazz, jintArray erasedLocationsBitmap, jintArray locationsToReadBitmap) {
	jint *erased_layout = (*env)->GetIntArrayElements(env, erasedLocationsBitmap, 0);
	jint *read_layout = (*env)->GetIntArrayElements(env, locationsToReadBitmap, 0);
	/* added by RH begins */
	jint locBitmapSize = (*env)->GetArrayLength(env, locationsToReadBitmap);
	jintArray retVal = (*env)->NewIntArray(env,locBitmapSize);
	/* added by RH ends */

	/* added by RH Dec 22nd begins */
	int i;
	bool* stripeErasedLayout=(bool*)calloc(sizeof(bool),locBitmapSize);
	int* stripeReadLayout=(int*)calloc(sizeof(int),locBitmapSize);
	for(i=0;i<locBitmapSize;i++){
		if(erased_layout[i]==1){
			if(i<num_of_parity_symbols){
				stripeErasedLayout[parityStripeMap[i]]=true;
				//printf("stripeErasedLayout[%d]=true\n",parityStripeMap[i]);
			}else{
				stripeErasedLayout[dataStripeMap[i-num_of_parity_symbols]]=true;
				//printf("stripeErasedLayout[%d]=true\n",dataStripeMap[i-num_of_parity_symbols]);
			}
		}
	}
	/* added by RH Dec 22nd ends */

	if (!FindLocationsToReadForDecode(n_stair, r_stair, m_stair, l_stair, 
				/* fixed by RH Dec 15th 2014 begins 
				 * updated by RH Dec 22nd 2014 passed mapped erased layout begins */
				//		m_partial_stair, error_vector_stair, erased_layout, read_layout)) {
				//return JNI_FALSE;
		m_partial_stair, error_vector_stair, stripeErasedLayout, stripeReadLayout)) {
			return NULL;
			/* updated by RH Dec 22nd 2014 passed mapped erased layout ends 
			 * fixed by RH Dec 15th 2014 ends */
		}

	/* fixed by RH Dec 22nd 2014 begins */
	//return JNI_TRUE;
	for(i=0;i<locBitmapSize;i++){
		if(stripeReadLayout[i]==1){
			//printf("stripeReadLayout[%d]=1\n",i);
			if(stripeDataMap[i]!=-1){
				//printf("stripeDataMap[%d]=%d\n",i,stripeDataMap[i]);
				read_layout[stripeDataMap[i]+num_of_parity_symbols]=1;
			}else{
				//printf("stripeParityMap[%d]=%d\n",i,stripeParityMap[i]);
				read_layout[stripeParityMap[i]]=1;
			}
		}
	}
	(*env)->SetIntArrayRegion(env,retVal,0,locBitmapSize,read_layout);
	free(stripeReadLayout);
	free(stripeErasedLayout);
	return retVal;
	/* fixed by RH Dec 22nd 2014 ends */
}

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
 jint numOutputBuffers, jbooleanArray erasedLocationsArray, jint numErasedLocations, jint dataLength) {
  //struct timeval timer;
  //gettimeofday(&timer,NULL);
  //printf("native decodeBulk checkpoint1 %lu\n",timer.tv_sec*1000000+timer.tv_usec);
	char* stripe[numInputBuffers];
	char* outputBufs[numOutputBuffers];
	bool erased_layout[numInputBuffers];
	int i;

	if ((dataLength <= 0) || (dataLength % (word_size_in_gf/8) != 0)) {
		fprintf(stderr, "Error: bad dataLength (%d), which should be a positive multiple of %d!\n", 
				dataLength, word_size_in_gf/8);
		return JNI_FALSE;
	}

	/* fixed by RH Dec 22nd 2014 begins 
	 * We pass mapped bufs */
	for (i = 0; i < numInputBuffers; i++) {
		jobject j_inputBuffer = (*env)->GetObjectArrayElement(env, inputBuffers, i);
    /* hard code begins */
    if(j_inputBuffer==NULL){
      //puts("NULL buffer!");
      continue;
    }
    /* hard code ends */
		if(i<num_of_parity_symbols){
			stripe[parityStripeMap[i]] = (char*)(*env)->GetDirectBufferAddress(env, j_inputBuffer);
		}else{
			stripe[dataStripeMap[i-num_of_parity_symbols]] = (char*)(*env)->GetDirectBufferAddress(env, j_inputBuffer);
		}
	}
	/* fixed by RH Dec 22nd 2014 ends */
  //gettimeofday(&timer,NULL);
  //printf("native decodeBulk checkpoint2 %lu\n",timer.tv_sec*1000000+timer.tv_usec);
  //printf("native decodeBulk checkpoint2\n");

	for (i = 0; i < numInputBuffers; i++) {
		erased_layout[i] = false;
	}
	int *erasedLocations = (int *)(*env)->GetIntArrayElements(env, erasedLocationsArray, 0);
	for (i = 0; i < numErasedLocations; i++) {
		jobject j_outputBuffer = (*env)->GetObjectArrayElement(env, outputBuffers, i);
		outputBufs[i] = (char*)(*env)->GetDirectBufferAddress(env, j_outputBuffer);
		//memset(outputBufs[i], 0, dataLength);
		//fprintf(stderr,"erasedLocation[%d]=%d\n",i,erasedLocations[i]);
		erased_layout[erasedLocations[i]] = true;
	}
	/* added by RH Dec 22nd begins */
	bool* stripeErasedLayout=(bool*)calloc(sizeof(bool),n_stair*r_stair);
	//int* stripe=(int*)calloc(sizeof(int),n_stair*r_stair);
	for(i=0;i<n_stair*r_stair;i++){
		if(erased_layout[i]==true){
			if(i<num_of_parity_symbols){
				stripeErasedLayout[parityStripeMap[i]]=true;
				//printf("stripeErasedLayout[%d]=true\n",parityStripeMap[i]);
			}else{
				stripeErasedLayout[dataStripeMap[i-num_of_parity_symbols]]=true;
				//printf("stripeErasedLayout[%d]=true\n",dataStripeMap[i-num_of_parity_symbols]);
			}
		}
	}
	/* added by RH Dec 22nd ends */
  //gettimeofday(&timer,NULL);
  //printf("native decodeBulk checkpoint3 %lu\n",timer.tv_sec*1000000+timer.tv_usec);

	/* added by RH Dec 24th begins 
	 * Generate a re-mapped output bufs */
	int* stripeIdxs=(int*)calloc(sizeof(int),numErasedLocations);
	char** stripeOutputs=(char**)calloc(sizeof(char*),numErasedLocations);
	for(i=0;i<numErasedLocations;i++){
		if(erasedLocations[i]<num_of_parity_symbols){
			stripeIdxs[i]=parityStripeMap[erasedLocations[i]];
		}else{
			stripeIdxs[i]=dataStripeMap[erasedLocations[i]-num_of_parity_symbols];
		}
	}
	int inputIdx;
	for(inputIdx=0;inputIdx<numErasedLocations;inputIdx++){
		int minimal=n_stair*r_stair;
		int minIdx=-1;
		for(i=0;i<numErasedLocations;i++){
			//get minimal
			if(stripeIdxs[i]<minimal){
				minimal=stripeIdxs[i];
				minIdx=i;
			}
		}
		if(minIdx!=-1){
			stripeOutputs[inputIdx]=outputBufs[minIdx];
			stripeIdxs[minIdx]=n_stair*r_stair;
		}
	}
  //gettimeofday(&timer,NULL);
  //printf("native decodeBulk checkpoint4 %lu\n",timer.tv_sec*1000000+timer.tv_usec);
  //printf("native decodeBulk checkpoint4 %lu\n",timer.tv_sec*1000000+timer.tv_usec);
	/* added by RH Dec 24th ends */

	//puts("nativeDecodeBulk() before DecodeBulk()");
	//if(row_coding_matrix_stair==NULL){
	//  fprintf(stderr,"row_coding_matrix_stair not initiated!!!\n");
	//}
	if (!DecodeBulk(gf_obj, n_stair, r_stair, m_stair, l_stair, m_partial_stair, 
				error_vector_stair, parity_layout_stair, row_coding_matrix_stair, 
				/* fixed by RH Dec 22nd 2014 by passing the mapped bufs begins
				 * updated by RH Dec 24th 2014 re-mapping output bufs begins */
				//column_coding_matrix_stair, stripe, dataLength, erased_layout, outputBufs)) {
				//column_coding_matrix_stair, stripe, dataLength, stripeErasedLayout, outputBufs)) {
		column_coding_matrix_stair, stripe, dataLength, stripeErasedLayout, stripeOutputs)) {
			/* updated by RH Dec 24th 2014 re-mapping output bufs ends
			 * fixed by RH Dec 22nd 2014 by passing the mapped bufs */
			free(stripeErasedLayout);
			free(stripeIdxs);
			free(stripeOutputs);
			return JNI_FALSE;
		}
  //gettimeofday(&timer,NULL);
  //printf("native decodeBulk checkpoint5 %lu\n",timer.tv_sec*1000000+timer.tv_usec);

	//puts("nativeDecodeBulk() after DecodeBulk()");
	free(stripeErasedLayout);
	free(stripeIdxs);
	free(stripeOutputs);

	return JNI_TRUE;
}
