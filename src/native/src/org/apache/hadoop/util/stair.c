/**
 * @file stair.c
 * @brief C implementation of stair codes
 * @author Mingqiang Li, mingqiangli.cn@gmail.com
 * @version 0.1
 * @date 2014-12-04
 * @note revised on 2015-01-09
 */
#include "stair.h"

/**
 * @brief InitParityLayout - initialize parity layout
 *
 * @param n_stair - parameter "n" in the stair code
 * @param r_stair - parameter "r" in the stair code
 * @param m_stair - parameter "m" in the stair code
 * @param l_stair - parameter "l" in the stair code
 * @param m_partial_stair - parameter "m_partial" in the stair code
 * @param error_vector_stair - parameter "error_vector" in the stair code
 * @param parity_layout_stair - parameter "parity_layout" in the stair code
 */
void InitParityLayout(int n_stair, int r_stair, int m_stair, int l_stair,
		int m_partial_stair, int *error_vector_stair, bool **parity_layout_stair,
		int** parityStripeMap,int** stripeParityMap,int** dataStripeMap,int** stripeDataMap) {
	int array_size;
	int i, j;

	array_size = n_stair * r_stair;

	/*initialize the parity layout*/
	(*parity_layout_stair) = TypeAlloc(bool, array_size);
	for (i = 0; i < array_size; i++) {
		(*parity_layout_stair)[i] = 0;
	}
	for (i = 0; i < m_stair; i++) {
		j = n_stair - m_stair + i;
		while (j < array_size) {
			(*parity_layout_stair)[j] = 1;
			j += n_stair;
		}
	}
	for (i = 0; i < m_partial_stair; i++) {
		for (j = 0; j < error_vector_stair[i]; j++) {
			(*parity_layout_stair)[n_stair * (r_stair - l_stair - 1 - j) + 
				(n_stair - m_stair - m_partial_stair + i)] = 1;
		}
	}
	for (i = 0; i < l_stair; i++) {
		for (j = 0; j < n_stair - m_stair; j++) {
			(*parity_layout_stair)[n_stair * (r_stair - l_stair + i) + j] = 1;
		}
	}

	/* Added by RH Dec 20th begins 
	 * We create maps mapping the index of stripe to parity/data and vice versa */
	(*parityStripeMap)=(int*)calloc(sizeof(int),array_size);
	(*stripeParityMap)=(int*)calloc(sizeof(int),array_size);
	(*dataStripeMap)=(int*)calloc(sizeof(int),array_size);
	(*stripeDataMap)=(int*)calloc(sizeof(int),array_size);
	int parityIdx=0;
	int dataIdx=0;
	for(i=0;i<array_size;i++){
		if((*parity_layout_stair)[i]){
			//printf("stripeIdx: %d parityIdx: %d\n",i,parityIdx);
			(*stripeParityMap)[i]=parityIdx;
			(*parityStripeMap)[parityIdx]=i;
			(*stripeDataMap)[i]=-1;
			parityIdx++;
		}else{
			//printf("stripeIdx: %d dataIdx: %d\n",i,dataIdx);
			(*stripeDataMap)[i]=dataIdx;
			(*dataStripeMap)[dataIdx]=i;
			(*stripeParityMap)[i]=-1;
			dataIdx++;
		}
	}
	/* Added by RH Dec 20th ends */
}

/**
 * @brief InitGF - initialize the gf_t object
 *
 * @param n_stair - parameter "n" in the stair code
 * @param r_stair - parameter "r" in the stair code
 * @param m_stair - parameter "m" in the stair code
 * @param l_stair - parameter "l" in the stair code
 * @param m_partial_stair - parameter "m_partial" in the stair code
 * @param error_vector_stair - parameter "error_vector" in the stair code
 * @param word_size_pointer - pointer to word_size_in_gf 
 * @param gf_pointer - pointer to gf_obj
 *
 * @return - a boolean value to indicate if the operation is successful
 */
bool InitGF(int n_stair, int r_stair, int m_stair, int l_stair, int m_partial_stair, 
		int *error_vector_stair, int *word_size_pointer, gf_t *gf_pointer) {
	/*initialize the GF-Complete instance gf_obj*/
	(*word_size_pointer) = 8; /*such a GF word size setting is enough for most practical configurations*/

	if (((n_stair + m_partial_stair) > 
				(1<<(*word_size_pointer))) || 
			((r_stair + error_vector_stair[m_partial_stair - 1]) > 
			 (1<<(*word_size_pointer)))) {
		fprintf(stderr, "Error: too large input parameters ---the internal GF word size setting (%d) is not enough for \
				the input parameters!\n", (*word_size_pointer));
		return 0;
	}

	if (!gf_init_easy(gf_pointer, (*word_size_pointer))) {
		fprintf(stderr, "Error: bad gf specification!\n");
		return 0;
	}

	return 1;
}

/**
 * @brief InitCodingMatrices - initialize coding matrices
 *
 * @param gf_obj - gf_t object
 * @param n_stair - parameter "n" in the stair code
 * @param r_stair - parameter "r" in the stair code
 * @param m_stair - parameter "m" in the stair code
 * @param l_stair - parameter "l" in the stair code
 * @param m_partial_stair - parameter "m_partial" in the stair code
 * @param error_vector_stair - parameter "error_vector" in the stair code
 * @param row_coding_matrix_stair - parameter "row_coding_matrix" in the stair code
 * @param column_coding_matrix_stair - parameter "column_coding_matrix" in the stair code
 */
void InitCodingMatrices(gf_t gf_obj, int n_stair, int r_stair, int m_stair, int l_stair, int m_partial_stair, 
		int *error_vector_stair, int **row_coding_matrix_stair, int **column_coding_matrix_stair) {
	int xor_sum;
	int i, j;

	/*initialize basic coding matrices*/
	//fprintf(stderr,"InitCodingMatrices() checkpoint1\n");
	(*row_coding_matrix_stair) = TypeAlloc(int, (n_stair - m_stair) * (m_stair + m_partial_stair));
	for (i = 0; i < m_stair + m_partial_stair; i++) {
		for (j = 0; j < n_stair - m_stair; j++) {
			xor_sum = i ^ ((m_stair + m_partial_stair) + j);
			(*row_coding_matrix_stair)[(n_stair - m_stair) * i + j] = gf_obj.divide.w32(&gf_obj, 1, xor_sum);
		}
	}
	//fprintf(stderr,"InitCodingMatrices() checkpoint2\n");

	(*column_coding_matrix_stair) = 
		TypeAlloc(int, (r_stair - l_stair) * (l_stair + error_vector_stair[m_partial_stair - 1]));
	for (i = 0; i < l_stair + error_vector_stair[m_partial_stair - 1]; i++) {
		for (j = 0; j < r_stair - l_stair; j++) {
			xor_sum = i ^ ((l_stair + error_vector_stair[m_partial_stair - 1]) + j);
			(*column_coding_matrix_stair)[(r_stair - l_stair) * i + j] = gf_obj.divide.w32(&gf_obj, 1, xor_sum);
		}
	}	
	//fprintf(stderr,"InitCodingMatrices() checkpoint3\n");
}

/**
 * @brief ChooseEncodingMethod - choose the encoding method that can cause the smallest number of Mult_XORs
 *
 * @param n_stair - parameter "n" in the stair code
 * @param r_stair - parameter "r" in the stair code
 * @param m_stair - parameter "m" in the stair code
 * @param l_stair - parameter "l" in the stair code
 * @param m_partial_stair - parameter "m_partial" in the stair code
 * @param error_vector_stair - parameter "error_vector" in the stair code
 * @param parity_layout_stair - parameter "parity_layout" in the stair code
 *
 * @return - the best encoding method among STD_METHOD, UP_METHOD, and DOWN_METHOD
 */
int ChooseEncodingMethod(int n_stair, int r_stair, int m_stair, int l_stair, 
		int m_partial_stair, int *error_vector_stair, bool *parity_layout_stair) {
	int array_size, num_of_global_parity_symbols;
	int mult_xor_std, mult_xor_up, mult_xor_down;
	int row_max, column_max, *riser, *tread, step_limit, step_cnt;
	int i, j, k;

	array_size = n_stair * r_stair;
	num_of_global_parity_symbols = 0;
	for (i = 0; i < m_partial_stair; i++) {
		num_of_global_parity_symbols += error_vector_stair[i];
	}
	//fprintf(stderr,"ChooseEncodingMethod() checkpoint 1\n");

	/*1. calculate the number of Mult_XORs involved by standard encoding*/

	/*(a) a parity symbol in Row i_0 and Column j_0 depends only on the data symbols d_{i,j}'s where i <= i_0 and j <= j_0*/

	mult_xor_std = 0;	
	for (i = 0; i < array_size; i++) {
		if (parity_layout_stair[i]) {
			row_max = i/n_stair;
			column_max = i%n_stair;

			for (j = 0; j <= row_max; j++) {
				for (k = 0; k <= column_max; k++) {
					if (!parity_layout_stair[n_stair * j + k]) mult_xor_std++;
				}
			}
		}
	}
	//fprintf(stderr,"ChooseEncodingMethod() checkpoint 2\n");

	/*(b) each parity symbol is unrelated to any data symbol in any other column (row) spanned by the same tread (riser)*/
	if (n_stair <= r_stair) {
		step_limit = n_stair;
	}
	else {
		step_limit = r_stair;
	}
	riser = TypeAlloc(int, step_limit);
	tread = TypeAlloc(int, step_limit);
	for (i = 0; i < step_limit; i++) {
		riser[i] = -1;
		tread[i] = -1;
	}

	i = 0;
	if (l_stair > 0) {
		riser[i] = l_stair;
		i++;
	}
	j = 0;
	riser[i] = error_vector_stair[j];
	i++;
	for (j = 1; j < m_partial_stair; j++) {
		if (error_vector_stair[j] > error_vector_stair[j - 1]) {
			riser[i] = error_vector_stair[j] - error_vector_stair[j - 1];
			i++;
		}
	}
	if (m_stair > 0) {
		riser[i] = r_stair - error_vector_stair[m_partial_stair - 1];
	}

	i = 0;
	if (l_stair > 0) {
		tread[i] = n_stair - m_stair - m_partial_stair;
		i++;
	}
	k = 1;
	for (j = 1; j < m_partial_stair; j++) {
		if (error_vector_stair[j] > error_vector_stair[j - 1]) {
			tread[i] = k;
			i++;
			k = 1;
		}
		else {
			k++;
		}			
	}
	tread[i] = k;
	i++;	
	if (m_stair > 0) {
		tread[i] = m_stair;
	}

	step_cnt = 0;
	while ((step_cnt < step_limit) && (riser[step_cnt] != -1)) step_cnt++;

	k = 0;
	for (i = 0; i < step_cnt; i++) {
		k += riser[i];		
		mult_xor_std -= (r_stair - k) * k * ((tread[i] - 1) * tread[i] / 2);		
	}

	k = 0;
	for (i = step_cnt - 1; i >= 0 ; i--) {
		k += tread[i];		
		mult_xor_std -= (n_stair - k) * k * ((riser[i] - 1) * riser[i] / 2);		
	}	

	free(riser);
	free(tread);
	//fprintf(stderr,"ChooseEncodingMethod() checkpoint 3\n");

	/*2. calculate the number of Mult_XORs involved by upstairs encoding*/

	mult_xor_up = (n_stair - m_stair) * 
		(m_stair * r_stair + num_of_global_parity_symbols) + 
		(r_stair - l_stair) * 
		((n_stair - m_stair) * 
		 (l_stair + error_vector_stair[m_partial_stair - 1]));
	//fprintf(stderr,"ChooseEncodingMethod() checkpoint 4\n");

	/*3. calculate the number of Mult_XORs involved by downstairs encoding*/

	mult_xor_down = (n_stair - m_stair) * 
		((r_stair - l_stair) * 
		 (m_stair + m_partial_stair)) + 
		(r_stair - l_stair) * 
		(n_stair * l_stair + num_of_global_parity_symbols);

#if DEBUG_INFO
	printf("The number of GF Mult_XORs in standard encoding: %d\n", mult_xor_std);
	printf("The number of GF Mult_XORs in upstairs encoding: %d\n", mult_xor_up);
	printf("The number of GF Mult_XORs in downstairs encoding: %d\n", mult_xor_down);
#endif

	if ((mult_xor_up <= mult_xor_std) && (mult_xor_up <= mult_xor_down)) {
#if DEBUG_INFO
		printf("Thus, we choose upstairs encoding\n\n");
#endif

		return UP_METHOD;
	}
	else if ((mult_xor_down <= mult_xor_std) && (mult_xor_down <= mult_xor_up)) {
#if DEBUG_INFO
		printf("Thus, we choose downstairs encoding\n\n");
#endif

		return DOWN_METHOD;
	}
	else {
#if DEBUG_INFO
		printf("Thus, we choose standard encoding\n\n");
#endif

		return STD_METHOD;
	}
}

/**
 * @brief InvertMatrix - invert a square matrix in GF 
 *
 * @param gf_obj - gf_t object
 * @param square_matrix - the square matrix to be inverted 
 *                 (note: this function will eventually change square_matrix into an identiy matrix)
 * @param inverse_matrix - the matrix that stores the result of the inverse matrix
 * @param matrix_order - the order of the square matrix
 *
 * @return - a boolean value to indicate if the operation is successful
 */
bool InvertMatrix(gf_t gf_obj, int *square_matrix, int *inverse_matrix, int matrix_order) {
	int matrix_size;
	int index_row_start1, index_row_start2;
	int tmp, mult_factor;
	int i, j, k, l;

	matrix_size = matrix_order * matrix_order;

	/*first store an identity matrix in inverse_matrix*/
	i = 0;
	while (i < matrix_size) {
		if (i / matrix_order == i % matrix_order) {
			inverse_matrix[i] = 1;
		}
		else {
			inverse_matrix[i] = 0;
		}

		i++;
	}

	/*convert square_matrix into an upper triangular matrix*/
	for (i = 0; i < matrix_order; i++) {
		index_row_start1 = matrix_order * i;

		/*if the i-th element in the i-th row is zero, we need to swap the i-th row with 
		  a row (below it) whose i-th element is non-zero*/
		if (square_matrix[index_row_start1 + i] == 0) { 
			j = i + 1;
			while ((j < matrix_order) && (square_matrix[matrix_order * j + i] == 0)) {
				j++;
			}
			/*if we cannot find such a row below the i-th row, we can judge that square_matrix is noninvertible*/
			if (j == matrix_order) { 
				return 0;
			}

			/*swap the i-th row with the j-th row for both square_matrix and inverse_matrix*/
			index_row_start2 = matrix_order * j;

			for (k = 0; k < matrix_order; k++) {
				tmp = square_matrix[index_row_start1 + k];
				square_matrix[index_row_start1 + k] = square_matrix[index_row_start2 + k];
				square_matrix[index_row_start2 + k] = tmp;

				/*do the same for inverse_matrix*/
				tmp = inverse_matrix[index_row_start1 + k];
				inverse_matrix[index_row_start1 + k] = inverse_matrix[index_row_start2 + k];
				inverse_matrix[index_row_start2 + k] = tmp;
			}
		}

		tmp = square_matrix[index_row_start1 + i];
		/*if the i-th element in the i-th row is not equal to 1, 
		  divide each element in this row by the i-th element*/
		if (tmp != 1) { 
			mult_factor = gf_obj.divide.w32(&gf_obj, 1, tmp);

			for (j = 0; j < matrix_order; j++) { 
				square_matrix[index_row_start1 + j] = gf_obj.multiply.w32(&gf_obj, 
						square_matrix[index_row_start1 + j], mult_factor);

				/*do the same for inverse_matrix*/
				inverse_matrix[index_row_start1 + j] = gf_obj.multiply.w32(&gf_obj, 
						inverse_matrix[index_row_start1 + j], mult_factor);
			}
		}

		/*multiply the i-th row with a factor and add it to each row below it 
		  such that the i-th element in each row becomes zero*/
		for (j = i + 1; j < matrix_order; j++) {
			index_row_start2 = matrix_order * j;
			k = index_row_start2 + i;

			if (square_matrix[k] != 0) { 				
				if (square_matrix[k] == 1) {
					for (l = 0; l < matrix_order; l++) {
						square_matrix[index_row_start2 + l] ^= square_matrix[index_row_start1 + l];

						/*do the same for inverse_matrix*/
						inverse_matrix[index_row_start2 + l] ^= inverse_matrix[index_row_start1 + l];
					}					
				}
				else {
					mult_factor = square_matrix[k];

					for (l = 0; l < matrix_order; l++) {
						square_matrix[index_row_start2 + l] ^= gf_obj.multiply.w32(&gf_obj, 
								square_matrix[index_row_start1 + l], mult_factor);

						/*do the same for inverse_matrix*/
						inverse_matrix[index_row_start2 + l] ^= gf_obj.multiply.w32(&gf_obj, 
								inverse_matrix[index_row_start1 + l], mult_factor);
					}					
				}
			}
		}
	}	

	/*based on the upper triangular matrix, make square_matrix become an identity matrix. 
	  then, inverse_matrix become the final inverse matrix*/
	for (i = matrix_order - 1; i >= 0; i--) {
		index_row_start1 = matrix_order * i;

		for (j = 0; j < i; j++) {
			index_row_start2 = matrix_order * j;
			k = index_row_start2 + i;

			if (square_matrix[k] != 0) { 	
				if (square_matrix[k] == 1) {
					for (l = 0; l < matrix_order; l++) {
						/*square_matrix[index_row_start2 + l] ^= square_matrix[index_row_start1 + l];*/

						/*do the same for inverse_matrix*/
						inverse_matrix[index_row_start2 + l] ^= inverse_matrix[index_row_start1 + l];
					}					
				}
				else {
					mult_factor = square_matrix[k];

					for (l = 0; l < matrix_order; l++) {
						/*square_matrix[index_row_start2 + l] ^= gf_obj.multiply.w32(&gf_obj, 
						  square_matrix[index_row_start1 + l], mult_factor);*/

						/*do the same for inverse_matrix*/
						inverse_matrix[index_row_start2 + l] ^= gf_obj.multiply.w32(&gf_obj, 
								inverse_matrix[index_row_start1 + l], mult_factor);
					}					
				}

				/*we simply zero this element since square_matrix will eventually become an identity matrix*/
				square_matrix[k] = 0; 
			}
		}
	}

	return 1;
}

/**
 * @brief StandardizePCM - standardize the parity check matrix, such that
 *                all columns corresponding to parity symbols can together form an indentity matrix  
 *
 * @param gf_obj - gf_t object
 * @param n_stair - parameter "n" in the stair code
 * @param r_stair - parameter "r" in the stair code
 * @param m_stair - parameter "m" in the stair code
 * @param l_stair - parameter "l" in the stair code
 * @param m_partial_stair - parameter "m_partial" in the stair code
 * @param error_vector_stair - parameter "error_vector" in the stair code
 * @param parity_layout_stair - parameter "parity_layout" in the stair code
 * @param original_parity_check_matrix - the original parity check matrix
 * @param std_parity_check_matrix - the standardized parity check matrix
 *
 * @return - a boolean value to indicate if the operation is successful
 */
bool StandardizePCM(gf_t gf_obj, int n_stair, int r_stair, int m_stair, int l_stair, 
		int m_partial_stair, int *error_vector_stair, 
		bool *parity_layout_stair, int *original_parity_check_matrix, 
		int *std_parity_check_matrix) {
	int array_size, num_of_parity_symbols;
	int *encoder, *inverse;
	int i, j, k, h;

	array_size = n_stair * r_stair;
	num_of_parity_symbols = 0;
	num_of_parity_symbols += ((r_stair - l_stair) * m_stair);
	for (i = 0; i < m_partial_stair; i++) {
		num_of_parity_symbols += error_vector_stair[i];
	}
	num_of_parity_symbols += (n_stair * l_stair);

	encoder = TypeAlloc(int, num_of_parity_symbols * num_of_parity_symbols);
	inverse = TypeAlloc(int, num_of_parity_symbols * num_of_parity_symbols);

	j = 0;
	k = 0;
	for (i = 0; i < array_size; i++) {
		if (parity_layout_stair[i] == 1) {
			if((i % n_stair < n_stair - m_stair) && 
					(i / n_stair < r_stair - l_stair)) { /*an inside global parity symbol*/
				for (h = 0; h < num_of_parity_symbols; h++) {
					encoder[num_of_parity_symbols * h + 
						m_stair * (r_stair - l_stair) + 
						n_stair * l_stair + j] = 
						original_parity_check_matrix[array_size * h + i];
				}
				j++;
			}
			else { /*a row parity symbol*/
				for (h = 0; h < num_of_parity_symbols; h++) {
					encoder[num_of_parity_symbols * h + k] = 
						original_parity_check_matrix[array_size * h + i];
				}
				k++;
			}		
		}
	}

	if (!InvertMatrix(gf_obj, encoder, inverse, num_of_parity_symbols)) {
		fprintf(stderr, "Error: encounter a noninvertible matrix during the standardization of \
				the parity check matrix!\n");

		free(encoder);
		free(inverse);

		return 0;
	}

	for (i = 0; i < num_of_parity_symbols; i++) {
		for (j = 0; j < array_size; j++) {
			std_parity_check_matrix[array_size * i + j] = 0;
			for (k = 0; k < num_of_parity_symbols; k++) {
				std_parity_check_matrix[array_size * i + j] ^= gf_obj.multiply.w32(&gf_obj, 
						inverse[num_of_parity_symbols * i + k], 
						original_parity_check_matrix[array_size * k + j]);
			}				
		}
	}	

	free(encoder);
	free(inverse);

	return 1;
}

/**
 * @brief DoStdEncoding - calculate parity symbols from data symbols in a stripe using standard encoding method 
 *
 * @param gf_obj - gf_t object
 * @param n_stair - parameter "n" in the stair code
 * @param r_stair - parameter "r" in the stair code
 * @param m_stair - parameter "m" in the stair code
 * @param l_stair - parameter "l" in the stair code
 * @param m_partial_stair - parameter "m_partial" in the stair code
 * @param error_vector_stair - parameter "error_vector" in the stair code
 * @param parity_layout_stair - parameter "parity_layout" in the stair code
 * @param std_parity_check_matrix - the standardized parity check matrix
 * @param stripe - the stripe to be encoded
 * @param block_size - size of each data/parity block
 *
 * @return - a boolean value to indicate if the operation is successful
 */
bool DoStdEncoding(gf_t gf_obj, int n_stair, int r_stair, int m_stair, int l_stair, 
		int m_partial_stair, int *error_vector_stair, 
		bool *parity_layout_stair, int *std_parity_check_matrix, 
		char **stripe, int block_size) {
	int array_size, num_of_parity_symbols;
	int *pindex;
	int coef;
	int i, j;

	array_size = n_stair * r_stair;
	num_of_parity_symbols = 0;
	num_of_parity_symbols += ((r_stair - l_stair) * m_stair);
	for (i = 0; i < m_partial_stair; i++) {
		num_of_parity_symbols += error_vector_stair[i];
	}
	num_of_parity_symbols += (n_stair * l_stair);

	/*use pindex to record the position of parity symbol in each row*/
	pindex = TypeAlloc(int, num_of_parity_symbols);
	for (i = 0; i < num_of_parity_symbols; i++) {
		/*in each row of the standard parity check matrix of a stair code, there is 
		  an interesting phenomenon we can leverage: all the elements on the right of 
		  the parity symbol are zeros, because of the parity relations in STAIR codes*/
		j = array_size - 1;
		while((j > 0) && (std_parity_check_matrix[array_size * i + j] == 0)) {
			j--;
		}

		if (std_parity_check_matrix[array_size * i + j] != 1) { /*this element should be 1*/
			fprintf(stderr, "Error: the input parity check matrix is not a standardized one!\n");

			free(pindex);

			return 0;
		}

		pindex[i] = j;		
	}

	/*set all parity symbols in the stripe to be zeros before encoding*/
	for (i = 0; i < array_size; i++) {
		if (parity_layout_stair[i]) {
			memset(stripe[i], 0, block_size);
		}
	}

	/*calculate each parity symbol from data symbols*/
	for (i = 0; i < num_of_parity_symbols; i++) {
		for (j = 0; j < pindex[i]; j++) {
			coef = std_parity_check_matrix[array_size * i + j];
			if (coef != 0) {
				gf_obj.multiply_region.w32(&gf_obj, stripe[j], stripe[pindex[i]], coef, block_size, 1);
			}
		}		
	}	

	free(pindex);

	return 1;
}

/**
 * @brief StdEncodeStripe - encode a stripe using standard encoding method 
 *
 * @param gf_obj - gf_t object
 * @param n_stair - parameter "n" in the stair code
 * @param r_stair - parameter "r" in the stair code
 * @param m_stair - parameter "m" in the stair code
 * @param l_stair - parameter "l" in the stair code
 * @param m_partial_stair - parameter "m_partial" in the stair code
 * @param error_vector_stair - parameter "error_vector" in the stair code
 * @param parity_layout_stair - parameter "parity_layout" in the stair code
 * @param row_coding_matrix_stair - parameter "row_coding_matrix" in the stair code
 * @param column_coding_matrix_stair - parameter "column_coding_matrix" in the stair code
 * @param stripe - the stripe to be encoded
 * @param block_size - size of each data/parity block
 *
 * @return - a boolean value to indicate if the operation is successful
 */
bool StdEncodeStripe(gf_t gf_obj, int n_stair, int r_stair, int m_stair, int l_stair, int m_partial_stair, 
		int *error_vector_stair, bool *parity_layout_stair, int *row_coding_matrix_stair, 
		int *column_coding_matrix_stair, char **stripe, int block_size) {
	int array_size, num_of_parity_symbols;
	int *original_parity_check_matrix, *std_parity_check_matrix;
	int i, j, k, row_id;
#if DEBUG_INFO
	struct timeval tv;
	double t_start, t_end, t_transform, t_do_encoding;
#endif

	array_size = n_stair * r_stair;
	num_of_parity_symbols = 0;
	num_of_parity_symbols += ((r_stair - l_stair) * m_stair);
	for (i = 0; i < m_partial_stair; i++) {
		num_of_parity_symbols += error_vector_stair[i];
	}
	num_of_parity_symbols += (n_stair * l_stair);

	/*1. allocate two parity check matrices*/
	original_parity_check_matrix = TypeAlloc(int, array_size * num_of_parity_symbols);	
	std_parity_check_matrix = TypeAlloc(int, array_size * num_of_parity_symbols);

	/*2. create the original parity check matrix*/	
	for (i = 0; i < array_size * num_of_parity_symbols; i++) {
		original_parity_check_matrix[i] = 0;
	}
	/*(a) the first m * (r - l) rows corresponding to horizontal parity symbols*/
	for (i = 0; i < r_stair - l_stair; i++) {
		for (j = 0; j < m_stair; j++) {
			row_id = m_stair * i + j;
			for (k = 0; k < n_stair - m_stair; k++) {
				original_parity_check_matrix[array_size * row_id + n_stair * i + k] = 
					row_coding_matrix_stair[(n_stair - m_stair) * j + k];
			}
			original_parity_check_matrix[array_size * row_id + n_stair * i + 
				(n_stair - m_stair) + j] = 1;
		}
	}
	/*(b) the next n * l rows corresponding to vertical parity symbols*/
	for (i = 0; i < l_stair; i++) {
		for (j = 0; j < n_stair; j++) {
			row_id = m_stair * (r_stair - l_stair) + n_stair * i + j;
			for (k = 0; k < r_stair - l_stair; k++) {
				original_parity_check_matrix[array_size * row_id + n_stair * k + j] = 
					column_coding_matrix_stair[(r_stair - l_stair) * i + k];
			}
			original_parity_check_matrix[array_size * row_id + 
				n_stair * (r_stair - l_stair + i) + j] = 1;
		}
	}
	/*(c) the final s rows corresponding to global parity symbols */
	row_id = m_stair * (r_stair - l_stair) + n_stair * l_stair;
	for (i = 0; i < m_partial_stair; i++) {
		for (j = 0; j < error_vector_stair[i]; j++) {
			for (k = 0; k < array_size; k++) {
				if ((k % n_stair < n_stair - m_stair) && 
						(k / r_stair < r_stair - l_stair)) {
					original_parity_check_matrix[array_size * row_id + k] = gf_obj.multiply.w32(&gf_obj, 
							row_coding_matrix_stair[
							(n_stair - m_stair) * (m_stair + i) + 
							k % n_stair], 
							column_coding_matrix_stair[
							(r_stair - l_stair) * 
							(l_stair + error_vector_stair[m_partial_stair - 1] - 1 - j) + 
							k / n_stair]);
				}
			}
			row_id++;
		}
	}		

#if DEBUG_INFO
	gettimeofday(&tv, NULL);
	t_start = (double) tv.tv_sec + (double) tv.tv_usec * 1e-6;
#endif

	/*3. transform the original parity check matrix into a standard one*/
	if (!StandardizePCM(gf_obj, n_stair, r_stair, m_stair, l_stair, m_partial_stair, 
				error_vector_stair, parity_layout_stair, 
				original_parity_check_matrix, std_parity_check_matrix)) {
		fprintf(stderr, "Error: fail to transform the original parity check matrix into a standard one!\n");

		free(original_parity_check_matrix);
		free(std_parity_check_matrix);

		return 0;
	}

#if DEBUG_INFO
	gettimeofday(&tv, NULL);
	t_end = (double) tv.tv_sec + (double) tv.tv_usec * 1e-6;
	t_transform = t_end - t_start;
#endif

#if DEBUG_INFO
	gettimeofday(&tv, NULL);
	t_start = (double) tv.tv_sec + (double) tv.tv_usec * 1e-6;
#endif

	/*4. use the standard parity check matrix to do encoding*/
	if (!DoStdEncoding(gf_obj, n_stair, r_stair, m_stair, l_stair, m_partial_stair, 
				error_vector_stair, parity_layout_stair, 
				std_parity_check_matrix, stripe, block_size)) {
		fprintf(stderr, "Error: fail to encode parity symbols from data symbols in a standard way!\n");

		free(original_parity_check_matrix);
		free(std_parity_check_matrix);

		return 0;
	}

#if DEBUG_INFO
	gettimeofday(&tv, NULL);
	t_end = (double) tv.tv_sec + (double) tv.tv_usec * 1e-6;
	t_do_encoding = t_end - t_start;
#endif

#if DEBUG_INFO
	printf("Total data: %d * %d = %d bytes\n", block_size, array_size - num_of_parity_symbols, 
			block_size * (array_size - num_of_parity_symbols));	
	printf("Time: %.6lf seconds for standardizing the parity check matrix\n", t_transform); 
	printf("      %.6lf seconds for doing encoding \n", t_do_encoding); 
	printf("Speed: %.6lf MB/s (including standardizing)\n", 
			block_size * (array_size - num_of_parity_symbols) / (t_transform + t_do_encoding) / (1024 *1024));
	printf("       %.6lf MB/s (excluding standardizing)\n\n", 
			block_size * (array_size - num_of_parity_symbols) / t_do_encoding / (1024 *1024));
#endif

	free(original_parity_check_matrix);
	free(std_parity_check_matrix);

	return 1;
}

/**
 * @brief CheckMergeableSubops - check if two sub-operations can be merged 
 *
 * @param op1 - the 1st sub-operation
 * @param op2 - the 2nd sub-operation
 *
 * @return - a boolean value to indicate if the two sub-operations can be merged
 */
inline bool CheckMergeableSubops(int *op1, int *op2) {
	int i = 0;

	while ((op1[i] != -1) && (op2[i] != -1) && (op1[i] == op2[i])) i++;

	if ((op1[i] != -1) || (op2[i] != -1)) return 0;

	return 1;	
}

/**
 * @brief DoUpstairsScheduling - do a scheduling on sub-operations for upstairs decoding 
 *
 * @param n_stair - parameter "n" in the stair code
 * @param r_stair - parameter "r" in the stair code
 * @param m_stair - parameter "m" in the stair code
 * @param l_stair - parameter "l" in the stair code
 * @param m_partial_stair - parameter "m_partial" in the stair code
 * @param error_vector_stair - parameter "error_vector" in the stair code
 * @param erased_layout - the layout of erased blocks in a stripe
 * @param subop_list - a list that records sub-operations, each of which has the following format:
 *			             {
 *                                 direction(row/column: 0/1), 
 *				        row/column index, 
 *				        a list of column/row indices for the input symbols, 
 *				        -1 (separator), 
 *				        a list of column/row indices for the output symbols, 
 *				        -1 (separator)
 *                                }
 *
 * @return - a boolean value to indicate if the operation is successful
 *
 * fixed by RH Dec 17th 2014 change erased_layout to bool type
 */
bool DoUpstairsScheduling(int n_stair, int r_stair, int m_stair, int l_stair, int m_partial_stair, 
		int *error_vector_stair, bool *erased_layout, int **subop_list) {
	int array_size, cnt_row_canonical_stripe, cnt_column_canonical_stripe;
	int **bitmap_canonical_stripe, *available_cnt_row, *available_cnt_column;
	int erased_cnt, *erased_cnt_row, *erased_cnt_column;   
	int outside_global_cnt_row, outside_global_index_start_row; 
	int *column_sorted_by_erased_cnt, first_column_recovered_by_global; 	
	int virtual_row_cnt_to_be_decoded, virtual_row_index; 
	int iteration_continue; /*if 1, continue to do the iteration; otherwize stop*/
	int subop_order; /*record the execution order of the sub-operation*/	
	int input_index, output_index;
	int tmp, i, j, k, l;

	array_size = n_stair * r_stair;		
	cnt_row_canonical_stripe = r_stair + error_vector_stair[m_partial_stair - 1];
	cnt_column_canonical_stripe = n_stair + m_partial_stair;

	/*# of erased symbols in each row*/
	erased_cnt_row = TypeAlloc(int, r_stair); 
	/*# of erased symbols in each column*/
	erased_cnt_column = TypeAlloc(int, n_stair); 
	/*# of available symbols in each row of the canonical stripe*/
	available_cnt_row = TypeAlloc(int, cnt_row_canonical_stripe); 
	/*# of available symbols in each column of the canonical stripe*/
	available_cnt_column = TypeAlloc(int, cnt_column_canonical_stripe); 

	column_sorted_by_erased_cnt = TypeAlloc(int, n_stair);

	for (i = 0; i < r_stair; i++) erased_cnt_row[i] = 0;
	for (i = 0; i < n_stair; i++) erased_cnt_column[i] = 0; 
	for (i = 0; i < cnt_row_canonical_stripe; i++) available_cnt_row[i] = 0;
	for (i = 0; i < cnt_column_canonical_stripe; i++) available_cnt_column[i] = 0; 

	bitmap_canonical_stripe = TypeAlloc(int *, cnt_row_canonical_stripe);
	for (i = 0; i < cnt_row_canonical_stripe; i++) {
		bitmap_canonical_stripe[i] = TypeAlloc(int, cnt_column_canonical_stripe);
		for (j = 0; j < cnt_column_canonical_stripe; j++) bitmap_canonical_stripe[i][j] = 0;			
	}

	/*at first, all (data and row parity) symbols are available*/
	for (i = 0; i < r_stair; i++) {
		for (j = 0; j < n_stair; j++) {
			bitmap_canonical_stripe[i][j] = 1;
			available_cnt_row[i]++;
			available_cnt_column[j]++;
		}
	}

	/*also, outside global parity symbols are available*/
	for (i = 0; i < m_partial_stair; i++) {
		for (j = 0; j < error_vector_stair[i]; j++) {
			bitmap_canonical_stripe[r_stair + j][n_stair + i] = 1;
			available_cnt_row[r_stair + j]++;
			available_cnt_column[n_stair + i]++;
		}
	}

	/*some symbols are erased according to erased_layout[]*/
	erased_cnt = 0;
	for (i = 0; i < array_size; i++) {
		if (erased_layout[i]) {
			bitmap_canonical_stripe[i / n_stair][i % n_stair] = 0;
			erased_cnt++;
			erased_cnt_row[i / n_stair]++;
			erased_cnt_column[i % n_stair]++;
			available_cnt_row[i / n_stair]--;
			available_cnt_column[i % n_stair]--;
		}
	}

	iteration_continue = 1; /*the first iteration is always performed*/
	subop_order = 0; /*first find the first sub-operation*/

	/*Step 1: Try to reconstruct lost symbols using vertical parity symbols*/
	for (i = 0; i < n_stair; i++) {
		/*if this column can be recovered using vertical parity symbols, add a sub-operation*/
		if ((erased_cnt_column[i] > 0) && (available_cnt_column[i] >= r_stair - l_stair)) { 						
			subop_list[subop_order][0] = 1;
			subop_list[subop_order][1] = i;

			input_index = 2;
			output_index = 2 + (r_stair - l_stair) + 1;
			for (j = 0; j < r_stair; j++) {
				if (bitmap_canonical_stripe[j][i] == 1) { /*we meet an available one*/						
					if (input_index - 2 < r_stair - l_stair) {
						subop_list[subop_order][input_index] = j;

						input_index++;
					}
				}
				else { /*we meet an erased one*/
					subop_list[subop_order][output_index] = j;
					bitmap_canonical_stripe[j][i] = 1; 
					erased_cnt--;
					erased_cnt_row[j]--;
					erased_cnt_column[i]--;
					available_cnt_row[j]++;
					available_cnt_column[i]++;

					output_index++;
				}
			}

			subop_order++;
		}
	}

	/*Step 2: Try to reconstruct lost symbols using global and vertical parity symbols column by column*/
	while ((erased_cnt > 0) && (iteration_continue == 1)) {
		/*turn off the iteration flag first. later, if we can recover a symbol in column direction, 
		  we will turn on this flag again*/
		iteration_continue = 0; 

		/*sort columns by their numbers of erased symbols in a non-increasing order and then 
		  store the indices into column_sorted_by_erased_cnt[]*/
		for (i = 0; i < n_stair; i++) column_sorted_by_erased_cnt[i] = n_stair - 1 - i;
		for (i = 0; i < n_stair; i++) {
			k = i;
			for (j = k + 1; j < n_stair; j++) {
				if (erased_cnt_column[column_sorted_by_erased_cnt[j]] > 
						erased_cnt_column[column_sorted_by_erased_cnt[k]]) {
					k = j;
				}
			}

			if (k != i) {
				tmp = column_sorted_by_erased_cnt[i];
				column_sorted_by_erased_cnt[i] = column_sorted_by_erased_cnt[k];
				column_sorted_by_erased_cnt[k] = tmp;
			}
		}

		i = m_stair; /*the m columns with the most erased symbols are left for reconstruction using row parity symbols*/
		/*but, what we can recover with global parity symbols cannot exceed error_vector[m_partial-1] + l*/
		while ((i < n_stair) && 
				(erased_cnt_column[column_sorted_by_erased_cnt[i]] > 
				 error_vector_stair[m_partial_stair - 1] + l_stair)) {			
			i++; 
		}
		first_column_recovered_by_global = i;
		virtual_row_cnt_to_be_decoded = erased_cnt_column[column_sorted_by_erased_cnt[i]] - l_stair;

		for (i = 0; i < virtual_row_cnt_to_be_decoded; i++) {
			virtual_row_index = r_stair + error_vector_stair[m_partial_stair - 1] - i -1;

			/*for each virtual row, we first count the number of outside global parity symbols*/
			outside_global_cnt_row = 0;
			for (j = 0; j < m_partial_stair; j++) {
				if (error_vector_stair[j] > i) outside_global_cnt_row++;
			}
			outside_global_index_start_row = m_partial_stair - outside_global_cnt_row;

			/*before reconstruct the virtual row, 
			  we need to generate [(n-m) - outside_global_cnt_row] virtual symbols from good chunks*/
			k = 0;
			for (j = 0; (j < n_stair) && 
					(k < (n_stair - m_stair) - outside_global_cnt_row); j++) {
				if (erased_cnt_column[j] == 0){ /*we meet a good chunk*/
					subop_list[subop_order][0] = 1;
					subop_list[subop_order][1] = j;

					input_index = 2;
					for (input_index = 2; input_index < 2 + r_stair - l_stair; input_index++) {
						subop_list[subop_order][input_index] = input_index - 2;
					}

					output_index = 2 + r_stair - l_stair + 1;
					subop_list[subop_order][output_index] = virtual_row_index;
					bitmap_canonical_stripe[virtual_row_index][j] = 1; 
					available_cnt_row[virtual_row_index]++;
					available_cnt_column[j]++;

					subop_order++;
					k++;
				}
			}
			if (k < (n_stair - m_stair) - outside_global_cnt_row) { 
				/*we cannot generate enough virtual symbols*/
				fprintf(stderr, "Error: encounter an unrecoverable failure scenario in the %d-th virtual row\n", i);	

				for (i = 0; i < cnt_row_canonical_stripe; i++) free(bitmap_canonical_stripe[i]);
				free(bitmap_canonical_stripe);	

				free(erased_cnt_row);	
				free(erased_cnt_column);	
				free(available_cnt_row);	
				free(available_cnt_column);	

				free(column_sorted_by_erased_cnt);

				return 0;
			}

			/*now, we begin to reconstruct the virtual row*/
			subop_list[subop_order][0] = 0;
			subop_list[subop_order][1] = virtual_row_index;
			input_index = 2;
			for (j = 0; (j < n_stair) && 
					(input_index - 2 < (n_stair - m_stair) - outside_global_cnt_row); j++) {
				if (bitmap_canonical_stripe[virtual_row_index][j] == 1){ 
					/*we have generated a virtual symbol from this chunk just now. so, 
					  we can use the virtual symbol here directly*/
					subop_list[subop_order][input_index] = j;

					input_index++;
				}
			}
			for (j = outside_global_index_start_row; j < m_partial_stair; j++) { 
				/*also use outside global parity symbols*/
				subop_list[subop_order][input_index] = n_stair + j;

				input_index++;
			}

			/*reconstruct all virtual symbols in this row is a bad choice. so, 
			  reconstruct only what may help our reconstruction*/
			output_index = 2 + (n_stair - m_stair) + 1;
			j = first_column_recovered_by_global;
			while ((j < n_stair) && 
					/*this column shouldn't have been recovered in the previous loop of i*/
					(erased_cnt_column[column_sorted_by_erased_cnt[j]] > i) && 
					/*the corresponding virtual symbol is unavailable*/
					(bitmap_canonical_stripe[virtual_row_index][column_sorted_by_erased_cnt[j]] == 0)) {
				subop_list[subop_order][output_index] = column_sorted_by_erased_cnt[j];
				bitmap_canonical_stripe[virtual_row_index][column_sorted_by_erased_cnt[j]] = 1;
				available_cnt_row[virtual_row_index]++;
				available_cnt_column[column_sorted_by_erased_cnt[j]]++;

				output_index++;
				j++;
			}
			subop_order++;

			/*check recoverable columns. recover it if find one*/
			for (j = 0; j < n_stair; j++) {
				/*some symbols are erased, but the available symbols are more than the threshold of r-l*/
				if ((erased_cnt_column[j] > 0) && (available_cnt_column[j] >= r_stair - l_stair)) { 
					iteration_continue = 1;

					subop_list[subop_order][0] = 1;
					subop_list[subop_order][1] = j;

					input_index = 2;
					output_index = 2 + r_stair - l_stair + 1;
					for (k = 0; k < r_stair; k++) {
						if (bitmap_canonical_stripe[k][j] == 1) { /*available symbols in the chunk*/
							subop_list[subop_order][input_index] = k;

							input_index++;
						}
						else {/*recover it*/
							subop_list[subop_order][output_index] = k;
							bitmap_canonical_stripe[k][j] = 1; 
							erased_cnt--;
							erased_cnt_row[k]--;
							erased_cnt_column[j]--;
							available_cnt_row[k]++;
							available_cnt_column[j]++;

							output_index++;
						}
					}

					k = r_stair + error_vector_stair[m_partial_stair - 1] - 1;
					while ((k >= virtual_row_index) && (input_index - 2 < r_stair - l_stair)) {
						subop_list[subop_order][input_index] = k; /*available virtual symbols in the same column*/

						input_index++;
						k--;
					}

					subop_order++;
				}
			}
		}
	}

	/*Step 3: Try to reconstruct lost symbols using horizontal parity symbols*/
	for (i = 0; i < r_stair; i++) {
		/*if this row can be recovered using horizontal parity symbols, add a sub-operation*/
		if ((erased_cnt_row[i] > 0) && (available_cnt_row[i] >= n_stair - m_stair)) { 						
			subop_list[subop_order][0] = 0;
			subop_list[subop_order][1] = i;

			input_index = 2;
			output_index = 2 + (n_stair - m_stair) + 1;
			for (j = 0; j < n_stair; j++) {
				if (bitmap_canonical_stripe[i][j] == 1) { /*we meet an available one*/						
					if (input_index - 2 < n_stair - m_stair) {
						subop_list[subop_order][input_index] = j;

						input_index++;
					}
				}
				else { /*we meet an erased one*/
					subop_list[subop_order][output_index] = j;
					bitmap_canonical_stripe[i][j] = 1; 
					erased_cnt--;
					erased_cnt_row[i]--;
					erased_cnt_column[j]--;
					available_cnt_row[i]++;
					available_cnt_column[j]++;

					output_index++;
				}
			}

			subop_order++;
		}
	}

	if (erased_cnt > 0) {
		fprintf(stderr, "Error: %d erased symbols cannot be recovered!\n", erased_cnt);

		for (i = 0; i < cnt_row_canonical_stripe; i++) free(bitmap_canonical_stripe[i]);
		free(bitmap_canonical_stripe);	

		free(erased_cnt_row);	
		free(erased_cnt_column);	
		free(available_cnt_row);	
		free(available_cnt_column);

		free(column_sorted_by_erased_cnt);	

		return 0;
	}

	/*we'd better merge the sub-operations that use the same set of symbols as the input*/	
	for (i = subop_order - 1; i > 0; i--){
		j = i - 1;
		while ((j >= 0) && (!CheckMergeableSubops(subop_list[i], subop_list[j]))) j--;

		if (j >= 0) { /*find a mergeable sub-operation*/
			k = 0;
			while (subop_list[j][k] != -1) {
				subop_list[i][k] = -1;  /*erase this sub-operation*/

				k++;
			}

			k++;
			l = k;
			while (subop_list[j][k] != -1) k++;

			while (subop_list[i][l] != -1) {
				subop_list[j][k] = subop_list[i][l];
				subop_list[i][l] = -1;  /*erase this sub-operation*/

				l++;
				k++;
			}

			subop_order--;
		}		
	}
	/*move the remaining sub-operations to the front of subop_list[] 
	  since there may be some empty entries in the middle of the list*/
	j = 0;
	for (i = 0; i < subop_order; i++) {
		while (subop_list[j][0] == -1) j++;

		if (j > i) {
			k = 0;
			while (subop_list[j][k] != -1) {
				subop_list[i][k] = subop_list[j][k];
				subop_list[j][k] = -1;  /*erase this sub-operation*/

				k++;
			}

			k++;
			while (subop_list[j][k] != -1) {
				subop_list[i][k] = subop_list[j][k];
				subop_list[j][k] = -1;  /*erase this sub-operation*/

				k++;
			}			
		}

		j++;
	}

	for (i = 0; i < cnt_row_canonical_stripe; i++) free(bitmap_canonical_stripe[i]);
	free(bitmap_canonical_stripe);	

	free(erased_cnt_row);	
	free(erased_cnt_column);	
	free(available_cnt_row);	
	free(available_cnt_column);	

	free(column_sorted_by_erased_cnt);

	return 1;
}

/**
 * @brief DoDownstairsScheduling - do a scheduling on sub-operations for downstairs decoding 
 *
 * @param n_stair - parameter "n" in the stair code
 * @param r_stair - parameter "r" in the stair code
 * @param m_stair - parameter "m" in the stair code
 * @param l_stair - parameter "l" in the stair code
 * @param m_partial_stair - parameter "m_partial" in the stair code
 * @param error_vector_stair - parameter "error_vector" in the stair code
 * @param erased_layout - the layout of erased blocks in a stripe
 * @param subop_list - a list that records sub-operations, each of which has the following format:
 *			             {
 *                                 direction(row/column: 0/1), 
 *				        row/column index, 
 *				        a list of column/row indices for the input symbols, 
 *				        -1 (separator), 
 *				        a list of column/row indices for the output symbols, 
 *				        -1 (separator)
 *                                }
 *
 * @return - a boolean value to indicate if the operation is successful
 *
 * fixed by RH Dec 17th 2014 change erased_layout to bool type
 */
bool DoDownstairsScheduling(int n_stair, int r_stair, int m_stair, int l_stair, int m_partial_stair, 
		int *error_vector_stair, bool *erased_layout, int **subop_list) {	
	int *t_error_vector_stair;
	int array_size, cnt_row_canonical_stripe, cnt_column_canonical_stripe;
	int **bitmap_canonical_stripe, *available_cnt_row, *available_cnt_column;
	int erased_cnt, *erased_cnt_row, *erased_cnt_column;  
	int outside_global_cnt_column, outside_global_index_start_column; 
	int *row_sorted_by_erased_cnt, first_row_recovered_by_global; 	
	int virtual_column_cnt_to_be_decoded, virtual_column_index; 
	int iteration_continue; /*if 1, continue to do the iteration; otherwize stop*/
	int subop_order; /*record the execution order of the sub-operation*/	
	int input_index, output_index;
	int tmp, i, j, k, l;

	/*generate the transposed error vector*/
	t_error_vector_stair = TypeAlloc(int, error_vector_stair[m_partial_stair - 1]); 
	for (i = 0; i < error_vector_stair[m_partial_stair - 1]; i++) {
		j = 0;
		while ((j < m_partial_stair) && 
				(error_vector_stair[m_partial_stair - 1 - j] >= 
				 (error_vector_stair[m_partial_stair - 1] - i))) {
			j++;
		}
		t_error_vector_stair[i] = j;
	}

	array_size = n_stair * r_stair;		
	cnt_row_canonical_stripe = r_stair + error_vector_stair[m_partial_stair - 1];
	cnt_column_canonical_stripe = n_stair + m_partial_stair;

	/*# of erased symbols in each row*/
	erased_cnt_row = TypeAlloc(int, r_stair); 
	/*# of erased symbols in each column*/
	erased_cnt_column = TypeAlloc(int, n_stair); 
	/*# of available symbols in each row of the canonical stripe*/
	available_cnt_row = TypeAlloc(int, cnt_row_canonical_stripe); 
	/*# of available symbols in each column of the canonical stripe*/
	available_cnt_column = TypeAlloc(int, cnt_column_canonical_stripe); 

	row_sorted_by_erased_cnt = TypeAlloc(int, r_stair);

	for (i = 0; i < r_stair; i++) erased_cnt_row[i] = 0;
	for (i = 0; i < n_stair; i++) erased_cnt_column[i] = 0; 
	for (i = 0; i < cnt_row_canonical_stripe; i++) available_cnt_row[i] = 0;
	for (i = 0; i < cnt_column_canonical_stripe; i++) available_cnt_column[i] = 0; 

	bitmap_canonical_stripe = TypeAlloc(int *, cnt_row_canonical_stripe);
	for (i = 0; i < cnt_row_canonical_stripe; i++) {
		bitmap_canonical_stripe[i] = TypeAlloc(int, cnt_column_canonical_stripe);
		for (j = 0; j < cnt_column_canonical_stripe; j++) bitmap_canonical_stripe[i][j] = 0;			
	}

	/*at first, all (data and row parity) symbols are available*/
	for (i = 0; i < r_stair; i++) {
		for (j = 0; j < n_stair; j++) {
			bitmap_canonical_stripe[i][j] = 1;
			available_cnt_row[i]++;
			available_cnt_column[j]++;
		}
	}

	/*also,outside global parity symbols are available*/
	for (i = 0; i < m_partial_stair; i++) {
		for (j = 0; j < error_vector_stair[i]; j++) {
			bitmap_canonical_stripe[r_stair + j][n_stair + i] = 1;
			available_cnt_row[r_stair + j]++;
			available_cnt_column[n_stair + i]++;
		}
	}

	/*some symbols are erased according to erased_layout[]*/
	erased_cnt = 0;
	for (i = 0; i < array_size; i++) {
		if (erased_layout[i]) {
			bitmap_canonical_stripe[i / n_stair][i % n_stair] = 0;
			erased_cnt++;
			erased_cnt_row[i / n_stair]++;
			erased_cnt_column[i % n_stair]++;
			available_cnt_row[i / n_stair]--;
			available_cnt_column[i % n_stair]--;
		}
	}

	iteration_continue = 1; /*the first iteration is always performed*/
	subop_order = 0; /*first find the first sub-operation*/

	/*Step 1: Try to reconstruct lost symbols using horizontal parity symbols*/
	for (i = 0; i < r_stair; i++) {
		/*if this row can be recovered using horizontal parity symbols, add a sub-operation*/
		if ((erased_cnt_row[i] > 0) && (available_cnt_row[i] >= n_stair - m_stair)) { 						
			subop_list[subop_order][0] = 0;
			subop_list[subop_order][1] = i;

			input_index = 2;
			output_index = 2 + (n_stair - m_stair) + 1;
			for (j = 0; j < n_stair; j++) {
				if (bitmap_canonical_stripe[i][j] == 1) { /*we meet an available one*/						
					if (input_index - 2 < n_stair - m_stair) {
						subop_list[subop_order][input_index] = j;

						input_index++;
					}
				}
				else { /*we meet an erased one*/
					subop_list[subop_order][output_index] = j;
					bitmap_canonical_stripe[i][j] = 1; 
					erased_cnt--;
					erased_cnt_row[i]--;
					erased_cnt_column[j]--;
					available_cnt_row[i]++;
					available_cnt_column[j]++;

					output_index++;
				}
			}

			subop_order++;
		}
	}

	/*Step 2: Try to reconstruct lost symbols using global and horizontal parity symbols row by row*/
	while ((erased_cnt > 0) && (iteration_continue == 1)) {
		/*turn off the iteration flag first. later, if we can recover a symbol in row direction, 
		  we will turn on this flag again*/
		iteration_continue = 0; 

		/*sort rows by their numbers of erased symbols in a non-increasing order and then 
		  store the indices into row_sorted_by_erased_cnt[]*/
		for (i = 0; i < r_stair; i++) row_sorted_by_erased_cnt[i] = r_stair - 1 - i;
		for (i = 0; i < r_stair; i++) {
			k = i;
			for (j = k + 1; j < r_stair; j++) {
				if (erased_cnt_row[row_sorted_by_erased_cnt[j]] > 
						erased_cnt_row[row_sorted_by_erased_cnt[k]]) {
					k = j;
				}
			}

			if (k != i) {
				tmp = row_sorted_by_erased_cnt[i];
				row_sorted_by_erased_cnt[i] = row_sorted_by_erased_cnt[k];
				row_sorted_by_erased_cnt[k] = tmp;
			}
		}

		i = l_stair; /*the l rows with the most erased symbols are left for reconstruction using vertical parity symbols*/
		/*but, what we can recover with global parity symbols cannot exceed m_partial + m*/
		while ((i < r_stair) && 
				(erased_cnt_row[row_sorted_by_erased_cnt[i]] > m_partial_stair+ m_stair)) {			
			i++; 
		}
		first_row_recovered_by_global = i;
		virtual_column_cnt_to_be_decoded = erased_cnt_row[row_sorted_by_erased_cnt[i]] - m_stair;

		for (i = 0; i < virtual_column_cnt_to_be_decoded; i++) {
			virtual_column_index = n_stair + m_partial_stair - i -1;

			/*for each virtual column, we first count the number of outside global parity symbols*/
			outside_global_cnt_column = 0;
			for (j = 0; j < error_vector_stair[m_partial_stair - 1]; j++) {
				if (t_error_vector_stair[j] > i) outside_global_cnt_column++;
			}
			outside_global_index_start_column = 
				error_vector_stair[m_partial_stair - 1] - outside_global_cnt_column;
      fprintf(stderr,"outside_global_cnt_column=error_vector_stair[%d]=%d\n",
          m_partial_stair-1,error_vector_stair[m_partial_stair - 1] - outside_global_cnt_column);

			/*before reconstruct the virtual column, 
			  we need to generate [(r-l) - outside_global_cnt_column] virtual symbols from good rows*/
			k = 0;
			for (j = 0; (j < r_stair) && 
					(k < (r_stair - l_stair) - outside_global_cnt_column); j++) {
        fprintf(stderr,"i:%d k:%d erased_cnt_row[%d]:%d\n",i,k,j,erased_cnt_row[i]);
				if (erased_cnt_row[j] == 0){ /*we meet a good row*/
          fprintf(stderr," i:%d k:%d erased_cnt_row[%d]:%d\n",i,k,j,erased_cnt_row[i]);
					subop_list[subop_order][0] = 0;
					subop_list[subop_order][1] = j;

					input_index = 2;
					for (input_index = 2; input_index < 2 + n_stair - m_stair; input_index++) {
						subop_list[subop_order][input_index] = input_index - 2;
					}

					output_index = 2 + n_stair - m_stair + 1;
					subop_list[subop_order][output_index] = virtual_column_index;
					bitmap_canonical_stripe[j][virtual_column_index] = 1; 
					available_cnt_row[j]++;
					available_cnt_column[virtual_column_index]++;

					subop_order++;
					k++;
          fprintf(stderr,"  i:%d k:%d erased_cnt_row[%d]:%d\n",i,k,j,erased_cnt_row[i]);
          fprintf(stderr,"\n");
				}
			}
			if (k < (r_stair - l_stair) - outside_global_cnt_column) { 
				/*we cannot generate enough virtual symbols*/
        fprintf(stderr,"outside_global_cnt_column:%d k:%d\n",outside_global_cnt_column,k);
        for(i=0;i<error_vector_stair[m_partial_stair-1];i++){
          fprintf(stderr,"t_error_vector_stair[%d]:%d\n",i,t_error_vector_stair[i]);
        }
				fprintf(stderr, "Error: encounter an unrecoverable failure scenario in the %d-th virtual column\n", i);	

				free(t_error_vector_stair);

				for (i = 0; i < cnt_row_canonical_stripe; i++) free(bitmap_canonical_stripe[i]);
				free(bitmap_canonical_stripe);	

				free(erased_cnt_row);	
				free(erased_cnt_column);	
				free(available_cnt_row);	
				free(available_cnt_column);	

				free(row_sorted_by_erased_cnt);

				return 0;
			}

			/*now, we begin to reconstruct the virtual column*/
			subop_list[subop_order][0] = 1;
			subop_list[subop_order][1] = virtual_column_index;
			input_index = 2;
			for (j = 0; (j < r_stair) && 
					(input_index - 2 < (r_stair - l_stair) - outside_global_cnt_column); j++) {
				if (bitmap_canonical_stripe[j][virtual_column_index] == 1){ 
					/*we have generated a virtual symbol from this row just now. so, 
					  we can use the virtual symbol here directly*/
					subop_list[subop_order][input_index] = j;

					input_index++;
				}
			}
			for (j = outside_global_index_start_column; j < error_vector_stair[m_partial_stair - 1]; j++) { 
				/*also use outside global parity symbols*/
				subop_list[subop_order][input_index] = r_stair + j;

				input_index++;
			}

			/*reconstruct all virtual symbols in this column is a bad choice. so, 
			  reconstruct only what may help our reconstruction*/
			output_index = 2 + (r_stair - l_stair) + 1;
			j = first_row_recovered_by_global;
			while ((j < r_stair) && 
					/*this row shouldn't have been recovered in the previous loop of i*/
					(erased_cnt_row[row_sorted_by_erased_cnt[j]] > i) && 
					/*the corresponding virtual symbol is unavailable*/
					(bitmap_canonical_stripe[row_sorted_by_erased_cnt[j]][virtual_column_index] == 0)) {
				subop_list[subop_order][output_index] = row_sorted_by_erased_cnt[j];
				bitmap_canonical_stripe[row_sorted_by_erased_cnt[j]][virtual_column_index] = 1;
				available_cnt_row[row_sorted_by_erased_cnt[j]]++;
				available_cnt_column[virtual_column_index]++;

				output_index++;
				j++;
			}
			subop_order++;

			/*check recoverable rows. recover it if find one*/
			for (j = 0; j < r_stair; j++) {
				/*some symbols are erased, but the available symbols are more than the threshold of n-m*/
        fprintf(stderr,"* recoverable: j:%d erased_cnt_row[%d]=%d,available_cnt_column[%d]=%d\n",j);
				if ((erased_cnt_row[j] > 0) && (available_cnt_row[j] >= n_stair - m_stair)) { 
          fprintf(stderr,"** recoverable: j:%d\n",j);
					iteration_continue = 1;

					subop_list[subop_order][0] = 0;
					subop_list[subop_order][1] = j;

					input_index = 2;
					output_index = 2 + n_stair - m_stair + 1;
					for (k = 0; k < n_stair; k++) {
						if (bitmap_canonical_stripe[j][k] == 1) { /*available symbols in the row*/
							subop_list[subop_order][input_index] = k;

							input_index++;
						}
						else {/*recover it*/
							subop_list[subop_order][output_index] = k;
							bitmap_canonical_stripe[j][k] = 1; 
							erased_cnt--;
							erased_cnt_row[j]--;
							erased_cnt_column[k]--;
							available_cnt_row[j]++;
							available_cnt_column[k]++;

							output_index++;
						}
					}

					k = n_stair + m_partial_stair - 1;
					while ((k >= virtual_column_index) && (input_index - 2 < n_stair - m_stair)) {
						subop_list[subop_order][input_index] = k; /*available virtual symbols in the same row*/

						input_index++;
						k--;
					}

					subop_order++;
				}
			}
		}
	}

	/*Step 3: Try to reconstruct lost symbols using vertical parity symbols*/
	for (i = 0; i < n_stair; i++) {
		/*if this column can be recovered using vertical parity symbols, add a sub-operation*/
		if ((erased_cnt_column[i] > 0) && (available_cnt_column[i] >= r_stair - l_stair)) { 						
			subop_list[subop_order][0] = 1;
			subop_list[subop_order][1] = i;

			input_index = 2;
			output_index = 2 + (r_stair - l_stair) + 1;
			for (j = 0; j < r_stair; j++) {
				if (bitmap_canonical_stripe[j][i] == 1) { /*we meet an available one*/						
					if (input_index - 2 < r_stair - l_stair) {
						subop_list[subop_order][input_index] = j;

						input_index++;
					}
				}
				else { /*we meet an erased one*/
					subop_list[subop_order][output_index] = j;
					bitmap_canonical_stripe[j][i] = 1; 
					erased_cnt--;
					erased_cnt_row[j]--;
					erased_cnt_column[i]--;
					available_cnt_row[j]++;
					available_cnt_column[i]++;

					output_index++;
				}
			}

			subop_order++;
		}
	}

	if (erased_cnt > 0) {
		fprintf(stderr, "Error: %d erased symbols cannot be recovered!\n", erased_cnt);

		free(t_error_vector_stair);

		for (i = 0; i < cnt_row_canonical_stripe; i++) free(bitmap_canonical_stripe[i]);
		free(bitmap_canonical_stripe);	

		free(erased_cnt_row);	
		free(erased_cnt_column);	
		free(available_cnt_row);	
		free(available_cnt_column);

		free(row_sorted_by_erased_cnt);

		return 0;
	}

	/*we'd better merge the sub-operations that use the same set of symbols as the input*/	
	for (i = subop_order - 1; i > 0; i--){
		j = i - 1;
		while ((j >= 0) && (!CheckMergeableSubops(subop_list[i], subop_list[j]))) j--;

		if (j >= 0) { /*find a mergeable sub-operation*/
			k = 0;
			while (subop_list[j][k] != -1) {
				subop_list[i][k] = -1;  /*erase this sub-operation*/

				k++;
			}

			k++;
			l = k;
			while (subop_list[j][k] != -1) k++;

			while (subop_list[i][l] != -1) {
				subop_list[j][k] = subop_list[i][l];
				subop_list[i][l] = -1;  /*erase this sub-operation*/

				l++;
				k++;
			}

			subop_order--;
		}		
	}
	/*move the remaining sub-operations to the front of subop_list[] 
	  since there may be some empty entries in the middle of the list*/
	j = 0;
	for (i = 0; i < subop_order; i++) {
		while (subop_list[j][0] == -1) j++;

		if (j > i) {
			k = 0;
			while (subop_list[j][k] != -1) {
				subop_list[i][k] = subop_list[j][k];
				subop_list[j][k] = -1;  /*erase this sub-operation*/

				k++;
			}

			k++;
			while (subop_list[j][k] != -1) {
				subop_list[i][k] = subop_list[j][k];
				subop_list[j][k] = -1;  /*erase this sub-operation*/

				k++;
			}			
		}

		j++;
	}	

	free(t_error_vector_stair);

	for (i = 0; i < cnt_row_canonical_stripe; i++) free(bitmap_canonical_stripe[i]);
	free(bitmap_canonical_stripe);	

	free(erased_cnt_row);	
	free(erased_cnt_column);	
	free(available_cnt_row);	
	free(available_cnt_column);	

	free(row_sorted_by_erased_cnt);

	return 1;
}

/**
 * @brief DoDecodingSubops - do sub-operations involved in iterative decoding 
 *
 * @param gf_obj - gf_t object
 * @param n_stair - parameter "n" in the stair code
 * @param r_stair - parameter "r" in the stair code
 * @param m_stair - parameter "m" in the stair code
 * @param l_stair - parameter "l" in the stair code
 * @param m_partial_stair - parameter "m_partial" in the stair code
 * @param error_vector_stair - parameter "error_vector" in the stair code
 * @param row_coding_matrix_stair - parameter "row_coding_matrix" in the stair code
 * @param column_coding_matrix_stair - parameter "column_coding_matrix" in the stair code
 * @param canonical_stripe - the canonical stripe to be encoded
 * @param block_size - size of each data/parity block
 * @param erased_layout - the layout of erased blocks in a stripe
 * @param subop_list - a list that records sub-operations, each of which has the following format:
 *			             {
 *                                 direction(row/column: 0/1), 
 *				        row/column index, 
 *				        a list of column/row indices for the input symbols, 
 *				        -1 (separator), 
 *				        a list of column/row indices for the output symbols, 
 *				        -1 (separator)
 *                                }
 *
 * @return - a boolean value to indicate if the operation is successful
 *
 * fixed by RH Dec 17th 2014 change erased_layout to bool type
 */
bool DoDecodingSubops(gf_t gf_obj, int n_stair, int r_stair, int m_stair, int l_stair, 
		int m_partial_stair, int *error_vector_stair, 
		int *row_coding_matrix_stair, int *column_coding_matrix_stair, 
		char **canonical_stripe, int block_size,  bool *erased_layout, int **subop_list) {
	int array_size, cnt_row_canonical_stripe, cnt_column_canonical_stripe;
	int *row_std_parity_check_matrix, *column_std_parity_check_matrix;
	int *eindex, *sindex, syndrome_cnt, syndrome_cnt_max, cnt_symbols_to_be_decoded;
	int *encoder, *inverse;
	/*in std_decoding_matrix, all columns corresponding to erased symbols can together form an indentity matrix*/
	int *std_decoding_matrix, std_decoding_matrix_size_max, coef; 
	int *referred_symbol;
	int subop_order, index_in_a_subop;
	int i, j, k;

	array_size = n_stair * r_stair;		
	cnt_row_canonical_stripe = r_stair + error_vector_stair[m_partial_stair - 1];
	cnt_column_canonical_stripe = n_stair + m_partial_stair;

	//puts("DoDecodingSubOp() checkpoint1");
	/*generate standard parity check matrices for the systematic MDS codes in both row and column directions*/
	row_std_parity_check_matrix = TypeAlloc(int, 
			cnt_column_canonical_stripe * (m_stair + m_partial_stair));
	for (i = 0; i < m_stair + m_partial_stair; i++) {
		for (j = 0; j < n_stair - m_stair; j++) {
			//fprintf(stderr,"DoDecodingSubOp() %d %d\n",i,j);
			row_std_parity_check_matrix[cnt_column_canonical_stripe * i + j] = 
				row_coding_matrix_stair[(n_stair - m_stair) * i + j];
		}
		for (j = 0; j < m_stair + m_partial_stair; j++) {
			//fprintf(stderr,"DoDecodingSubOp() %d %d\n",i,j);
			if (i == j) {
				row_std_parity_check_matrix[cnt_column_canonical_stripe * i + 
					(n_stair - m_stair) + j] = 1;
			}
			else{
				row_std_parity_check_matrix[cnt_column_canonical_stripe * i + 
					(n_stair - m_stair) + j] = 0;
			}
		}
	}
	//puts("DoDecodingSubOp() checkpoint1.5");
	column_std_parity_check_matrix = TypeAlloc(int, 
			cnt_row_canonical_stripe * (l_stair + error_vector_stair[m_partial_stair - 1]));
	for (i = 0; i < l_stair + error_vector_stair[m_partial_stair - 1]; i++) {
		for (j = 0; j < r_stair - l_stair; j++) {
			column_std_parity_check_matrix[cnt_row_canonical_stripe * i + j] = 
				column_coding_matrix_stair[(r_stair - l_stair) * i + j];
		}
		for (j = 0; j < l_stair + error_vector_stair[m_partial_stair - 1]; j++) {
			if (i == j) {
				column_std_parity_check_matrix[cnt_row_canonical_stripe * i + 
					(r_stair - l_stair) + j] = 1;
			}
			else{
				column_std_parity_check_matrix[cnt_row_canonical_stripe * i + 
					(r_stair - l_stair) + j] = 0;
			}
		}
	}
	//puts("DoDecodingSubOp() checkpoint2");

	/*allocate some arrays that will be used during decoding*/
	syndrome_cnt_max = 
		(m_stair + m_partial_stair >= l_stair + error_vector_stair[m_partial_stair - 1]? 
		 m_stair + m_partial_stair : l_stair + error_vector_stair[m_partial_stair - 1]);

	eindex = TypeAlloc(int, syndrome_cnt_max);
	sindex = TypeAlloc(int, syndrome_cnt_max);

	encoder = TypeAlloc(int, syndrome_cnt_max * syndrome_cnt_max);
	inverse = TypeAlloc(int, syndrome_cnt_max * syndrome_cnt_max);

	std_decoding_matrix_size_max = 
		(cnt_column_canonical_stripe * (m_stair + m_partial_stair) >= 
		 cnt_row_canonical_stripe * (l_stair + error_vector_stair[m_partial_stair - 1])? 
		 cnt_column_canonical_stripe * (m_stair + m_partial_stair) : 
		 cnt_row_canonical_stripe * (l_stair + error_vector_stair[m_partial_stair - 1]));
	std_decoding_matrix = TypeAlloc(int, std_decoding_matrix_size_max);

	referred_symbol = TypeAlloc(int, (n_stair - m_stair > r_stair - l_stair? 
				n_stair - m_stair : r_stair - l_stair));
	//puts("DoDecodingSubOp() checkpoint3");

	/*set all erased symbols in the canonical stripe to be zeros before decoding*/
	//fprintf(stderr,"array size %d\n",array_size);
	for (i = 0; i < array_size; i++) {
		if (erased_layout[i]) {
      //puts("erased_layout");
			//fprintf(stderr,"array size: %d failed: %d\n",array_size,i);
			memset(canonical_stripe[cnt_column_canonical_stripe * (i / n_stair) + 
					i % n_stair], 0, block_size);
		}
	}	
	//puts("DoDecodingSubOp() checkpoint3.1");
	/*also set all virtual symbols in the canonical stripe to be zeros*/
	for (i = 0; i < r_stair; i++) {
		for (j = n_stair; j < cnt_column_canonical_stripe; j++) {
			//fprintf(stderr,"canonical_stripe[%d] set zero\n",cnt_column_canonical_stripe * i + j);
			memset(canonical_stripe[cnt_column_canonical_stripe * i + j], 0, block_size);
		}
	}
	//puts("DoDecodingSubOp() checkpoint3.2");
	for (i = r_stair; i < cnt_row_canonical_stripe; i++) {
		for (j = 0; j < cnt_column_canonical_stripe; j++) {
			//fprintf(stderr," canonical_stripe[%d] set zero\n",cnt_column_canonical_stripe * i + j);
			memset(canonical_stripe[cnt_column_canonical_stripe * i + j], 0, block_size);
		}
	}
	//puts("DoDecodingSubOp() checkpoint4");

	/*now, we begin to do each sub-operation*/
	subop_order = 0;
	while (subop_list[subop_order][0] != -1){
		if(subop_list[subop_order][0] == 0) { 
			/*this is a decoding in row direction*/
			//puts("DoDecodingSubOp() checkpoint4.1");
			for (i = 0; i < m_stair + m_partial_stair; i++) sindex[i] = -1;
			for (i = 0; i < n_stair - m_stair; i++) referred_symbol[i] = 0;

			k = 0;
			syndrome_cnt = 0;

			index_in_a_subop = 2 + (n_stair - m_stair) + 1;
			while (subop_list[subop_order][index_in_a_subop] != -1) { /*output part*/
				eindex[k] = subop_list[subop_order][index_in_a_subop]; /*really erased symbol*/

				if (subop_list[subop_order][index_in_a_subop] >= n_stair - m_stair) {
					/*pay attention to the index of sindex here*/
					sindex[k] = subop_list[subop_order][index_in_a_subop] - 
						(n_stair - m_stair); 
					syndrome_cnt++;
				}
				else {
					referred_symbol[subop_list[subop_order][index_in_a_subop]] = 1;
				}

				k++;
				index_in_a_subop++;
			}	
			cnt_symbols_to_be_decoded = k;

			index_in_a_subop = 2;
			i = 0;
			while (subop_list[subop_order][index_in_a_subop] != -1) { /*input part*/
				if (subop_list[subop_order][index_in_a_subop] >= n_stair - m_stair) {
					/*some may have been occupied*/
					while ((i < syndrome_cnt_max) && (sindex[i] != -1)) i++; 
					sindex[i] = subop_list[subop_order][index_in_a_subop] - 
						(n_stair - m_stair);
					i++;
					syndrome_cnt++;
				}
				else {
					referred_symbol[subop_list[subop_order][index_in_a_subop]] = 1;
				}

				index_in_a_subop++;
			}

			/*clearly, cnt_symbols_to_be_decoded <= syndrome_cnt*/
			if (cnt_symbols_to_be_decoded < syndrome_cnt) {
				/*add some virtually erased symbols to ensure the following encoder is a square matrix*/	
				k = cnt_symbols_to_be_decoded;
				for (j = 0; (j < n_stair - m_stair) && (k < syndrome_cnt); j++){ 
					if (referred_symbol[j] == 0){
						eindex[k] = j; /*virtually erased symbol, which will not be really reconstructed*/
						k++;
					}
				}
			}

			for (i = 0; i < syndrome_cnt; i++) {
				for (j = 0; j < syndrome_cnt; j++) {
					encoder[syndrome_cnt * i + j] = 
						row_std_parity_check_matrix[cnt_column_canonical_stripe * sindex[i] + eindex[j]];
				}
			}

			/*generate the inverse of encoder*/
			if (InvertMatrix(gf_obj, encoder, inverse, syndrome_cnt) == 0) {
				fprintf(stderr, "Error: a square submatrix of the row parity check matrix cannot be inverted!\n");

				free(row_std_parity_check_matrix);
				free(column_std_parity_check_matrix);

				free(referred_symbol);

				free(std_decoding_matrix);

				if (syndrome_cnt_max > 0) {
					free(encoder);
					free(inverse);

					free(eindex);
					free(sindex);
				}

				return 0;
			}

			/*calculate std_decoding_matrix, in which all columns corresponding to (both really and virtually) 
			  erased symbols can together form an indentity matrix*/
			for (i = 0; i <syndrome_cnt; i++) {
				for (j = 0; j <cnt_column_canonical_stripe; j++) {
					std_decoding_matrix[cnt_column_canonical_stripe * i + j] = 0;
					for (k = 0; k <syndrome_cnt; k++) {
						std_decoding_matrix[cnt_column_canonical_stripe * i + j] ^= 
							gf_obj.multiply.w32(&gf_obj, inverse[syndrome_cnt * i + k], 
									row_std_parity_check_matrix[cnt_column_canonical_stripe * sindex[k] + j]);
					}				
				}
			}

			/*reconstruct really erased symbols based on std_decoding_matrix*/
			for (i = 0; i < cnt_symbols_to_be_decoded; i++) {
				for (j = 0; j < cnt_column_canonical_stripe; j++) {
					if (eindex[i] != j) { /*not itself*/
						coef = std_decoding_matrix[cnt_column_canonical_stripe * i + j];
						if (coef != 0) {
							gf_obj.multiply_region.w32(&gf_obj, 
									canonical_stripe[
									cnt_column_canonical_stripe * subop_list[subop_order][1] + j
									], 
									canonical_stripe[
									cnt_column_canonical_stripe * subop_list[subop_order][1] + eindex[i]
									], 
									coef, block_size, 1);
						}
					}
				}		
			}				
		} 
		else { 
			/*this is a decoding in column direction*/
			//puts("DoDecodingSubOp() checkpoint4.2");
			for (i = 0; i < l_stair + error_vector_stair[m_partial_stair - 1]; i++) sindex[i] = -1;
			for (i = 0; i < r_stair - l_stair; i++) referred_symbol[i] = 0;

			k = 0;
			syndrome_cnt = 0;

			index_in_a_subop = 2 + (r_stair - l_stair) + 1;
			while (subop_list[subop_order][index_in_a_subop] != -1) { /*output part*/
				eindex[k] = subop_list[subop_order][index_in_a_subop]; /*really erased symbol*/

				if (subop_list[subop_order][index_in_a_subop] >= r_stair - l_stair) {
					/*pay attention to the index of sindex here*/
					sindex[k] = subop_list[subop_order][index_in_a_subop] - 
						(r_stair - l_stair); 
					syndrome_cnt++;
				}
				else {
					referred_symbol[subop_list[subop_order][index_in_a_subop]] = 1;
				}

				k++;
				index_in_a_subop++;
			}			
			cnt_symbols_to_be_decoded = k;

			index_in_a_subop = 2;
			i = 0;
			while (subop_list[subop_order][index_in_a_subop] != -1) { /*input part*/
				if (subop_list[subop_order][index_in_a_subop] >= r_stair - l_stair) {
					/*some may have been occupied*/
					while ((i < syndrome_cnt_max) && (sindex[i] != -1)) i++; 
					sindex[i] = subop_list[subop_order][index_in_a_subop] - 
						(r_stair - l_stair);
					i++;
					syndrome_cnt++;
				}
				else {
					referred_symbol[subop_list[subop_order][index_in_a_subop]] = 1;
				}

				index_in_a_subop++;
			}

			/*clearly, cnt_symbols_to_be_decoded <= syndrome_cnt*/
			if (cnt_symbols_to_be_decoded < syndrome_cnt) {
				/*add some virtually erased symbols to ensure the following encoder is a square matrix*/			
				k = cnt_symbols_to_be_decoded;
				for (j = 0; (j < r_stair - l_stair) && (k < syndrome_cnt); j++){ 
					if (referred_symbol[j] == 0){
						eindex[k] = j; /*virtually erased symbol, which will not be really reconstructed*/
						k++;
					}
				}
			}

			for (i = 0; i < syndrome_cnt; i++) {
				for (j = 0; j < syndrome_cnt; j++) {
					encoder[syndrome_cnt * i + j] = 
						column_std_parity_check_matrix[cnt_row_canonical_stripe * sindex[i] + eindex[j]];
				}
			}

			/*generate the inverse of encoder*/
			if (InvertMatrix(gf_obj, encoder, inverse, syndrome_cnt) == 0) {
				fprintf(stderr, "Error: a square submatrix of the column parity check matrix cannot be inverted!\n");

				free(row_std_parity_check_matrix);
				free(column_std_parity_check_matrix);

				free(referred_symbol);

				free(std_decoding_matrix);

				if (syndrome_cnt_max > 0) {
					free(encoder);
					free(inverse);

					free(eindex);
					free(sindex);
				}

				return 0;
			}

			/*calculate std_decoding_matrix, in which all columns corresponding to (both really and virtually) 
			  erased symbols can together form an indentity matrix*/
			for (i = 0; i < syndrome_cnt; i++) {
				for (j = 0; j < cnt_row_canonical_stripe; j++) {
					std_decoding_matrix[cnt_row_canonical_stripe * i + j] = 0;
					for (k = 0; k < syndrome_cnt; k++) {
						std_decoding_matrix[cnt_row_canonical_stripe * i + j] ^= 
							gf_obj.multiply.w32(&gf_obj, inverse[syndrome_cnt * i + k], 
									column_std_parity_check_matrix[cnt_row_canonical_stripe * sindex[k] + j]);
					}				
				}
			}

			/*reconstruct really erased symbols based on std_decoding_matrix*/
			for (i = 0; i < cnt_symbols_to_be_decoded; i++) {
				for (j = 0; j < cnt_row_canonical_stripe; j++) {
					if (eindex[i] != j) { /*not itself*/
						coef = std_decoding_matrix[cnt_row_canonical_stripe * i + j];
						if (coef != 0) {
							gf_obj.multiply_region.w32(&gf_obj, 
									canonical_stripe[
									cnt_column_canonical_stripe * j + subop_list[subop_order][1]
									], 
									canonical_stripe[
									cnt_column_canonical_stripe * eindex[i] + subop_list[subop_order][1]
									], 
									coef, block_size, 1);
						}
					}
				}		
			}				
		}

		subop_order++;
	}
	//puts("DoDecodingSubOp() checkpoint5");

	free(row_std_parity_check_matrix);
	free(column_std_parity_check_matrix);

	free(referred_symbol);

	free(std_decoding_matrix);

	if (syndrome_cnt_max > 0) {
		free(encoder);
		free(inverse);

		free(eindex);
		free(sindex);
	}

	return 1;
}

/**
 * @brief IterativeDecodeStripe - decode a stripe using an iterative (either upstairs or downstairs) method  
 *
 * @param gf_obj - gf_t object
 * @param n_stair - parameter "n" in the stair code
 * @param r_stair - parameter "r" in the stair code
 * @param m_stair - parameter "m" in the stair code
 * @param l_stair - parameter "l" in the stair code
 * @param m_partial_stair - parameter "m_partial" in the stair code
 * @param error_vector_stair - parameter "error_vector" in the stair code
 * @param row_coding_matrix_stair - parameter "row_coding_matrix" in the stair code
 * @param column_coding_matrix_stair - parameter "column_coding_matrix" in the stair code
 * @param stripe - the stripe to be encoded
 * @param block_size - size of each data/parity block
 * @param erased_layout - the layout of erased blocks in a stripe
 * @param method - a specified method, either upstairs (UP_METHOD) or downstairs (DOWN_METHOD)
 *
 * @return - a boolean value to indicate if the operation is successful
 * 
 * fixed by RH Dec 16th 2014 
 */
bool IterativeDecodeStripe(gf_t gf_obj, int n_stair, int r_stair, int m_stair, int l_stair, 
		int m_partial_stair, int *error_vector_stair, 
		int *row_coding_matrix_stair, int *column_coding_matrix_stair, 
		char **stripe, int block_size, bool *erased_layout, int method) {
	/* fixed by RH Dec 16th 2014 ends */
	int cnt_row_canonical_stripe, cnt_column_canonical_stripe;
	int array_size, canonical_array_size, num_of_parity_symbols,  num_of_virtual_symbols;
	char **canonical_stripe, **virtual_symbol_buffers;
	int row_index, column_index, canonical_stripe_index, virtual_symbol_index;
	int **subop_list, subop_cnt_max, subop_len_max;
	int i, j;
#if DEBUG_INFO
	struct timeval tv;
	double t_start, t_end, t_scheduling, t_do_subops;
#endif

	if ((method != UP_METHOD) && (method != DOWN_METHOD)) {
		fprintf(stderr, "Error: wrong method argument for IterativeDecodeStripe(): \
				neither UP_METHOD nor DOWN_METHOD!\n");		
			return 0;
	}

	array_size = n_stair * r_stair;	
	cnt_row_canonical_stripe = r_stair + error_vector_stair[m_partial_stair - 1];
	cnt_column_canonical_stripe = n_stair + m_partial_stair;
	canonical_array_size = cnt_column_canonical_stripe * cnt_row_canonical_stripe;
	num_of_parity_symbols = 0;
	num_of_parity_symbols += ((r_stair - l_stair) * m_stair);
	for (i = 0; i < m_partial_stair; i++) {
		num_of_parity_symbols += error_vector_stair[i];
	}
	num_of_parity_symbols += (n_stair * l_stair);
	num_of_virtual_symbols = canonical_array_size - array_size;

	canonical_stripe = TypeAlloc(char *, canonical_array_size);
	virtual_symbol_buffers = TypeAlloc(char *, num_of_virtual_symbols);
	for (virtual_symbol_index = 0; virtual_symbol_index < num_of_virtual_symbols; virtual_symbol_index++) {
		virtual_symbol_buffers[virtual_symbol_index] = TypeAlloc(char, block_size);
	}
	virtual_symbol_index = 0;
  //puts("IterativeDecodeStripe() checkpoint1");
	for (canonical_stripe_index = 0; canonical_stripe_index < canonical_array_size; canonical_stripe_index++) {
		/* fixed by RH Dec 17th 2014 begins 
		 * fatal bug... */
		/*row_index = canonical_array_size / cnt_column_canonical_stripe;*/
		/*column_index = canonical_array_size % cnt_column_canonical_stripe;*/
		row_index = canonical_stripe_index / cnt_column_canonical_stripe;
		column_index = canonical_stripe_index % cnt_column_canonical_stripe;
		/* fixed by RH Dec 17th 2014 ends */
		if ((row_index < r_stair) && (column_index < n_stair)) { /*a real block*/
			canonical_stripe[canonical_stripe_index] = stripe[n_stair * row_index + column_index];
		}
		else { /*a virtual block*/
			canonical_stripe[canonical_stripe_index] = virtual_symbol_buffers[virtual_symbol_index];
			virtual_symbol_index++;
		}
	}

	/*allocate space for recording sub-operations*/
	subop_cnt_max = 
		((n_stair - m_stair) * error_vector_stair[m_partial_stair - 1] >= 
		 (r_stair - l_stair) * m_partial_stair? 
		 (n_stair - m_stair) * error_vector_stair[m_partial_stair - 1] : 
		 (r_stair - l_stair) * m_partial_stair) + 
		num_of_parity_symbols + 1; 
	subop_len_max = 2 + 
		(cnt_column_canonical_stripe >= cnt_row_canonical_stripe? 
		 cnt_column_canonical_stripe : cnt_row_canonical_stripe) + 
		2; 
	subop_list = TypeAlloc(int *, subop_cnt_max); 
	for (i = 0; i < subop_cnt_max; i++) {
		subop_list[i] = TypeAlloc(int, subop_len_max);
		for (j = 0; j < subop_len_max; j++) subop_list[i][j] = -1;
	}

#if DEBUG_INFO
	gettimeofday(&tv, NULL);
	t_start = (double) tv.tv_sec + (double) tv.tv_usec * 1e-6;	
#endif

  //puts("IterativeDecodeStripe() checkpoint2");
	/*first make a scheduling on sub-operations*/
	if (method == UP_METHOD) {
		if (!DoUpstairsScheduling(n_stair, r_stair, m_stair, l_stair, m_partial_stair, 
					error_vector_stair, erased_layout, subop_list)) {
			fprintf(stderr, "Error: fail to recover all erased symbols!\n");

			free(canonical_stripe);

			for (virtual_symbol_index = 0; virtual_symbol_index < num_of_virtual_symbols; 
					virtual_symbol_index++) {
				free(virtual_symbol_buffers[virtual_symbol_index]);
			}
			free(virtual_symbol_buffers);

			for (i = 0; i < subop_cnt_max; i++) free(subop_list[i]);
			free(subop_list);

			return 0;
		}
	} else{
		if (!DoDownstairsScheduling(n_stair, r_stair, m_stair, l_stair, m_partial_stair, 
					error_vector_stair, erased_layout, subop_list)) {
			fprintf(stderr, "Error: fail to recover all erased symbols!\n");

			free(canonical_stripe);

			for (virtual_symbol_index = 0; virtual_symbol_index < num_of_virtual_symbols; 
					virtual_symbol_index++) {
				free(virtual_symbol_buffers[virtual_symbol_index]);
			}
			free(virtual_symbol_buffers);

			for (i = 0; i < subop_cnt_max; i++) free(subop_list[i]);
			free(subop_list);

			return 0;
		}
	}

#if DEBUG_INFO
	gettimeofday(&tv, NULL);
	t_end = (double) tv.tv_sec + (double) tv.tv_usec * 1e-6;
	t_scheduling = t_end - t_start;
#endif


#if DEBUG_INFO
	gettimeofday(&tv, NULL);
	t_start = (double) tv.tv_sec + (double) tv.tv_usec * 1e-6;
#endif	

	/*now, we begin to do each sub-operation*/
  //puts("IterativeDecodeStripe() checkpoint3");
	if (!DoDecodingSubops(gf_obj, n_stair, r_stair, m_stair, l_stair, m_partial_stair, 
				error_vector_stair, row_coding_matrix_stair, column_coding_matrix_stair, 
				canonical_stripe, block_size,  erased_layout, subop_list)) {
		fprintf(stderr, "Error: fail to execute all decoding subops!\n");

		free(canonical_stripe);

		for (virtual_symbol_index = 0; virtual_symbol_index < num_of_virtual_symbols; 
				virtual_symbol_index++) {
			free(virtual_symbol_buffers[virtual_symbol_index]);
		}
		free(virtual_symbol_buffers);

		for (i = 0; i < subop_cnt_max; i++) free(subop_list[i]);
		free(subop_list);

		return 0;
	}
  //puts("IterativeDecodeStripe() checkpoint4");

#if DEBUG_INFO
	gettimeofday(&tv, NULL);
	t_end = (double) tv.tv_sec + (double) tv.tv_usec * 1e-6;
	t_do_subops = t_end - t_start;
#endif

#if DEBUG_INFO
	printf("Total data: %d * %d = %d bytes\n", block_size, array_size-num_of_parity_symbols, block_size * (array_size - num_of_parity_symbols));	
	printf("Time: %.6lf seconds for scheduling sub-operations\n", t_scheduling); 
	printf("      %.6lf seconds for doing all sub-operations\n", t_do_subops); 
	printf("Speed: %.6lf MB/s (including scheduling)\n", block_size * (array_size - num_of_parity_symbols) / (t_scheduling + t_do_subops) / (1024 * 1024));
	printf("       %.6lf MB/s (excluding scheduling)\n\n", block_size * (array_size - num_of_parity_symbols) / t_do_subops / (1024 * 1024));
#endif

	free(canonical_stripe);

	for (virtual_symbol_index = 0; virtual_symbol_index < num_of_virtual_symbols; 
			virtual_symbol_index++) {
		free(virtual_symbol_buffers[virtual_symbol_index]);
	}
	free(virtual_symbol_buffers);

	for (i = 0; i < subop_cnt_max; i++) free(subop_list[i]);
	free(subop_list);

	return 1;
}

/**
 * @brief EncodeBulk - encode data into parity in bulk
 *
 * @param gf_obj - gf_t object
 * @param n_stair - parameter "n" in the stair code
 * @param r_stair - parameter "r" in the stair code
 * @param m_stair - parameter "m" in the stair code
 * @param l_stair - parameter "l" in the stair code
 * @param m_partial_stair - parameter "m_partial" in the stair code
 * @param error_vector_stair - parameter "error_vector" in the stair code
 * @param parity_layout_stair - parameter "parity_layout" in the stair code
 * @param best_encoding_method - the chosen best encoding method
 * @param row_coding_matrix_stair - parameter "row_coding_matrix" in the stair code
 * @param column_coding_matrix_stair - parameter "column_coding_matrix" in the stair code
 * @param data - data blocks to be encoded
 * @param coding - parity blocks to be generated
 * @param block_size - size of each data/parity block
 *
 * @return - a boolean value to indicate if the operation is successful
 */
bool EncodeBulk(gf_t gf_obj, int n_stair, int r_stair, int m_stair, int l_stair, int m_partial_stair, 
		int *error_vector_stair, bool *parity_layout_stair, int best_encoding_method, 
		int *row_coding_matrix_stair, int *column_coding_matrix_stair, 
		char **data, char **coding, int block_size) {
	int array_size;
	char **stripe;
	int data_index, coding_index, stripe_index;

	array_size = n_stair * r_stair;
	stripe = TypeAlloc(char *, array_size);

	data_index = 0;
	coding_index = 0;
	stripe_index = 0;
	for (stripe_index = 0; stripe_index < array_size; stripe_index++) {
		if(parity_layout_stair[stripe_index] == 0) { /*a data block*/
			stripe[stripe_index] = data[data_index];
			data_index++;
		}
		else { /*a coding block*/
			stripe[stripe_index] = coding[coding_index];
			coding_index++;
		}
	}

#if DEBUG_INFO
	PrintStripe(n_stair, r_stair, stripe, block_size);
#endif

	if (best_encoding_method == STD_METHOD) {
		if (!StdEncodeStripe(gf_obj, n_stair, r_stair, m_stair, l_stair,
					m_partial_stair, error_vector_stair, parity_layout_stair, 
					row_coding_matrix_stair, column_coding_matrix_stair,
					stripe, block_size)) {
			fprintf(stderr, "Error: fail to encode the stripe using the standard method!\n");

			free(stripe);

			return 0;
		}
	}

	if (best_encoding_method == UP_METHOD) {
		if (!IterativeDecodeStripe(gf_obj, n_stair, r_stair, m_stair, l_stair,
					m_partial_stair, error_vector_stair,
					row_coding_matrix_stair, column_coding_matrix_stair,
					stripe, block_size, parity_layout_stair, UP_METHOD)) {
			fprintf(stderr, "Error: fail to encode the stripe using the upstairs method!\n");

			free(stripe);

			return 0;
		}
	}

	if (best_encoding_method == DOWN_METHOD) {
		if (!IterativeDecodeStripe(gf_obj, n_stair, r_stair, m_stair, l_stair,
					m_partial_stair, error_vector_stair, 
					row_coding_matrix_stair, column_coding_matrix_stair,
					stripe, block_size, parity_layout_stair,  DOWN_METHOD)) {
			fprintf(stderr, "Error: fail to encode the stripe using the downstairs method!\n");

			free(stripe);

			return 0;
		}
	}

#if DEBUG_INFO
	PrintStripe(n_stair, r_stair, stripe, block_size);
#endif

	free(stripe);

	return 1;		
}

/**
 * @brief FindLocationsToReadForDecode - determine the locations to be read for decoding
 *
 * @param n_stair - parameter "n" in the stair code
 * @param r_stair - parameter "r" in the stair code
 * @param m_stair - parameter "m" in the stair code
 * @param l_stair - parameter "l" in the stair code
 * @param m_partial_stair - parameter "m_partial" in the stair code
 * @param error_vector_stair - parameter "error_vector" in the stair code
 * @param erased_layout - the layout of erased blocks in a stripe
 * @param read_layout - the layout of blocks to be read
 *
 * @return - a boolean value to indicate if the operation is successful
 *
 * fixed by RH Dec 17th 2014 change erased_layout to bool type
 */
bool FindLocationsToReadForDecode(int n_stair, int r_stair, int m_stair, int l_stair, 
		int m_partial_stair, int *error_vector_stair, bool *erased_layout, int *read_layout) {
	int array_size, cnt_column_canonical_stripe, cnt_row_canonical_stripe, num_of_parity_symbols;
	int **subop_list, subop_cnt_max, subop_len_max;
	int subop_order;
	int i, j;

	array_size = n_stair * r_stair;
	cnt_column_canonical_stripe = n_stair + m_partial_stair;
	cnt_row_canonical_stripe = r_stair + error_vector_stair[m_partial_stair - 1];
	num_of_parity_symbols = 0;
	num_of_parity_symbols += ((r_stair - l_stair) * m_stair);
	for (i = 0; i < m_partial_stair; i++) {
		num_of_parity_symbols += error_vector_stair[i];
	}
	num_of_parity_symbols += (n_stair * l_stair);

	/* Added by RH for debug Dec 17th 2014 begins */
	//for(i=0;i<array_size;i++){
	//  if(erased_layout[i]){
	//    printf("erased_layout[%d]=true\n",i);
	//  }else{
	//    printf("erased_layout[%d]=false\n",i);
	//  }
	//}
	/* Added by RH for debug Dec 17th 2014 ends */

	/*allocate space for recording sub-operations*/
	subop_cnt_max = 
		((n_stair - m_stair) * error_vector_stair[m_partial_stair - 1] >= 
		 (r_stair - l_stair) * m_partial_stair? 
		 (n_stair - m_stair) * error_vector_stair[m_partial_stair - 1] : 
		 (r_stair - l_stair) * m_partial_stair) + 
		num_of_parity_symbols + 1; 
	subop_len_max = 2 + 
		(cnt_column_canonical_stripe >= cnt_row_canonical_stripe? 
		 cnt_column_canonical_stripe : cnt_row_canonical_stripe) + 
		2; 
	subop_list = TypeAlloc(int *, subop_cnt_max); 
	for (i = 0; i < subop_cnt_max; i++) {
		subop_list[i] = TypeAlloc(int, subop_len_max);
		for (j = 0; j < subop_len_max; j++) subop_list[i][j] = -1;
	}

	/*first make a scheduling on sub-operations*/
	if (!DoUpstairsScheduling(n_stair, r_stair, m_stair, l_stair, m_partial_stair, 
				error_vector_stair, erased_layout, subop_list)) {
		fprintf(stderr, "Error: fail to recover all erased symbols!\n");

		for (i = 0; i < subop_cnt_max; i++) free(subop_list[i]);
		free(subop_list);

		return 0;
	}

	for (i = 0; i < array_size; i++) {
		read_layout[i] = 0;
	}

	/*now, we begin to find locations to be read for decoding*/
	subop_order = 0;
	while (subop_list[subop_order][0] != -1){
		if(subop_list[subop_order][0] == 0) { 
			/*this is a decoding in row direction*/
			i = 2;
			while (subop_list[subop_order][i] != -1) {
				if ((subop_list[subop_order][1] < r_stair) && 
						(subop_list[subop_order][i] < n_stair)) {
					read_layout[n_stair * subop_list[subop_order][1] + 
						subop_list[subop_order][i]] = 1;
					/* Added by RH for debug begins */
					//fprintf(stderr,"n_stair: %d subop_list[%d][1]: %d subop_list[%d][%d]: %d\n",
					//		subop_order,subop_list[subop_order][1],
					//		n_stair,subop_order,i, subop_list[subop_order][i]); 
					//fprintf(stderr,"read_layout[%d]=1\n",n_stair * subop_list[subop_order][1] + 
					//		subop_list[subop_order][i]);
					/* Added by RH for debug ends */
				}

				i++;
			}
		} else { 
			/*this is a decoding in column direction*/
			i = 2;
			while (subop_list[subop_order][i] != -1) {
				if ((subop_list[subop_order][1] < n_stair) && 
						(subop_list[subop_order][i] < r_stair)) {
					read_layout[n_stair * subop_list[subop_order][i] + 
						subop_list[subop_order][1]] = 1;
					/* Added by RH for debug begins */
					//fprintf(stderr,"read_layout[%d]=1\n",n_stair * subop_list[subop_order][i] + 
					//		subop_list[subop_order][1]);
					/* Added by RH for debug ends */
				}

				i++;
			}
		}

		subop_order++;
	}

	for (i = 0; i < subop_cnt_max; i++) free(subop_list[i]);
	free(subop_list);

	return 1;
}

/**
 * @brief DecodeBulk - decode erased blocks in bulk 
 *
 * @param gf_obj - gf_t object
 * @param n_stair - parameter "n" in the stair code
 * @param r_stair - parameter "r" in the stair code
 * @param m_stair - parameter "m" in the stair code
 * @param l_stair - parameter "l" in the stair code
 * @param m_partial_stair - parameter "m_partial" in the stair code
 * @param error_vector_stair - parameter "error_vector" in the stair code
 * @param parity_layout_stair - parameter "parity_layout" in the stair code
 * @param row_coding_matrix_stair - parameter "row_coding_matrix" in the stair code
 * @param column_coding_matrix_stair - parameter "column_coding_matrix" in the stair code
 * @param stripe - the stripe to be decoded
 * @param block_size - size of each data/parity block
 * @param erased_layout - the layout of erased blocks in a stripe
 * @param output_buffers - output buffers for storing the decoded blocks
 *
 * @return - a boolean value to indicate if the operation is successful
 *
 * fixed by RH Dec 17th 2014 change erased_layout to bool type
 */
bool DecodeBulk(gf_t gf_obj, int n_stair, int r_stair, int m_stair, int l_stair, int m_partial_stair, 
		int *error_vector_stair, bool *parity_layout_stair, int *row_coding_matrix_stair, 
		int *column_coding_matrix_stair, char **stripe, int block_size, 
		bool *erased_layout, char **output_buffers) {
	int array_size, output_index;
	int i;

#if DEBUG_INFO
	PrintStripe(n_stair, r_stair, stripe, block_size);
#endif

	//puts("DecodeBulk() checkpoint1");
	//fprintf(stderr,"DecodeBulk(): checkpoint1 %d\n",row_coding_matrix_stair[0]);
	if (!IterativeDecodeStripe(gf_obj, n_stair, r_stair, m_stair, l_stair, m_partial_stair, 
				error_vector_stair, row_coding_matrix_stair, column_coding_matrix_stair,	
				stripe, block_size, erased_layout,  UP_METHOD)) {
		fprintf(stderr, "Error: fail to decode the stripe using the downstairs method!\n");

		/* fixed by RH Dec 17th 2014 begins */
		//free(stripe);
		/* fixed by RH Dec 17th 2014 ends */

		return 0;
	}
	//puts("DecodeBulk(): checkpoint2");

	array_size = n_stair * r_stair;
	output_index = 0;
	for (i = 0; i < array_size; i++) {
		if (erased_layout[i]) {
			memcpy(output_buffers[output_index], stripe[i], block_size);
			output_index++;
		}
	}	
	//fprintf(stderr,"DecodeBulk(): checkpoint3");

#if DEBUG_INFO
	PrintStripe(n_stair, r_stair, stripe, block_size);
#endif

	/* fixed by RH Dec 17th 2014 begins */
	//free(stripe);
	/* fixed by RH Dec 17th 2014 ends */

	return 1;		
}

/**
 * @brief PrintStripe - print all data and parity blocks in a stripe
 *
 * @param n_stair - parameter "n" in the stair code
 * @param r_stair - parameter "r" in the stair code
 * @param stripe - the stripe to be printed
 * @param block_size - size of each data/parity block
 *
 * fixed by RH Dec 18th 2014 to print unsigned char
 */
void PrintStripe(int n_stair, int r_stair, char **stripe, int block_size) {
	int i, j;

	for (i = 0; i < r_stair; i++) {
		for (j = 0; j < n_stair; j++) {
			if (block_size == 1) {
				printf("%02x ", (unsigned char)stripe[n_stair * i + j][0]);	
			}
			else if (block_size == 2) {
				printf("%02x%02x ", (unsigned char)stripe[n_stair * i + j][0], 
						(unsigned char)stripe[n_stair * i + j][1]);
			}
			else {
				printf("%02x...%02x ", (unsigned char)stripe[n_stair * i + j][0], 
						(unsigned char)stripe[n_stair * i + j][block_size - 1]);	
			}
		}				
		printf("\n");		
	}		
	printf("\n");
}

