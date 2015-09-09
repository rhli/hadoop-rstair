/**
 * @file stair.h
 * @brief C implementation of stair codes
 * @author Mingqiang Li, mingqiangli.cn@gmail.com
 * @version 0.1
 * @date 2014-12-04
 * @note revised on 2015-01-09
 */
#ifndef _STAIR_H_
#define _STAIR_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

/*use GF-Complete library*/
/* fixed by RH Dec 13th, 2014 begins */
//#include "gf_complete.h"
#include "gf_complete.h"
/* fixed by RH Dec 13th, 2014 ends */

#define TypeAlloc(type, num) (type *) malloc(sizeof(type)*(num))

/*macro for standard encoding method*/
#define STD_METHOD 0 
/*macro for upstairs encoding method*/
#define UP_METHOD 1
/*macro for downstairs encoding method*/
#define DOWN_METHOD 2

#define ROW_FIRST_MODE 0
#define COLUMN_FIRST_MODE 1

#define DEBUG_INFO 0
#define WORST_ERASURE 1

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
void InitParityLayout(int n_stair, int r_stair, int m_stair, int l_stair, int m_partial_stair, 
		int *error_vector_stair, bool **parity_layout_stair,
		int** parityStripeMap,int** stripeParityMap,int** dataStripeMap,int** stripeDataMap) ;

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
		int *error_vector_stair, int *word_size_pointer, gf_t *gf_pointer);

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
		int *error_vector_stair, int **row_coding_matrix_stair, int **column_coding_matrix_stair);

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
int ChooseEncodingMethod(int n_stair, int r_stair, int m_stair, int l_stair, int m_partial_stair, 
		int *error_vector_stair, bool *parity_layout_stair);

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
bool InvertMatrix(gf_t gf_obj, int *square_matrix, int *inverse_matrix, int matrix_order);

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
bool StandardizePCM(gf_t gf_obj, int n_stair, int r_stair, int m_stair, int l_stair, int m_partial_stair, 
		int *error_vector_stair, bool *parity_layout_stair, int *original_parity_check_matrix, 
		int *std_parity_check_matrix);

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
bool DoStdEncoding(gf_t gf_obj, int n_stair, int r_stair, int m_stair, int l_stair, int m_partial_stair, 
		int *error_vector_stair, bool *parity_layout_stair, int *std_parity_check_matrix, 
		char **stripe, int block_size);

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
		int *column_coding_matrix_stair, char **stripe, int block_size);

/**
 * @brief CheckMergeableSubops - check if two sub-operations can be merged 
 *
 * @param op1 - the 1st sub-operation
 * @param op2 - the 2nd sub-operation
 *
 * @return - a boolean value to indicate if the two sub-operations can be merged
 */
bool CheckMergeableSubops(int *op1, int *op2);

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
		int *error_vector_stair, bool *erased_layout, int **subop_list);

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
		int *error_vector_stair, bool *erased_layout, int **subop_list);

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
bool DoDecodingSubops(gf_t gf_obj, int n_stair, int r_stair, int m_stair, int l_stair, int m_partial_stair, 
		int *error_vector_stair, int *row_coding_matrix_stair, int *column_coding_matrix_stair, 
		char **canonical_stripe, int block_size,  bool *erased_layout, int **subop_list);

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
bool IterativeDecodeStripe(gf_t gf_obj, int n_stair, int r_stair, int m_stair, int l_stair, int m_partial_stair, 
		int *error_vector_stair, int *row_coding_matrix_stair, int *column_coding_matrix_stair, 
		char **stripe, int block_size, bool *erased_layout, int method);
/* fixed by RH Dec 16th 2014 ends */

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
		char **data, char **coding, int block_size);

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
		int m_partial_stair, int *error_vector_stair, bool *erased_layout, int *read_layout);

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
		bool *erased_layout, char **output_buffers);

/**
 * @brief PrintStripe - print all data and parity blocks in a stripe
 *
 * @param n_stair - parameter "n" in the stair code
 * @param r_stair - parameter "r" in the stair code
 * @param stripe - the stripe to be printed
 * @param block_size - size of each data/parity block
 */
void PrintStripe(int n_stair, int r_stair, char **stripe, int block_size);

#endif
